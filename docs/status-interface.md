# 运行状态接口约定（schema_version 1）

> 适用：DBReducer（本仓库）。pAnno 第 2 步托管 pFind 时沿用同一套规则，只把 `tool`、日志前缀换成 `pAnno` / `PANNO`。
> 读取方：pAnno（判断是否继续流程）、GUI（显示失败阶段、原因与处理建议）。
> 状态：**已确认（2026-10-06）**，退出码采用 10–16 的新编号。

一次运行通过三条通道对外报告，各司其职：

| 通道 | 用途 | 谁读 |
|---|---|---|
| **退出码** | 一个整数说清"成功 / 哪一类失败" | pAnno（`system()` 的返回值，第 2 步起改为 `CreateProcess`）|
| **状态文件** `dbr_status.json` | 完整的结构化信息：阶段、进度、失败进程、参数、建议 | GUI、pAnno |
| **日志行** `[DBR:事件] key=value …` | 实时事件流，GUI 边跑边解析 | GUI（实时日志面板）|

---

## 1. 退出码

| 码 | `result` 名称 | 含义 | 状态文件 |
|---|---|---|---|
| 0 | `ok` | 成功，`DBReducer.fasta` 已写出且非空 | 有 |
| 1 | `param_error` | cfg 参数错误、输入文件不存在（保留现有 `exit(1)` 的语义）| 有 |
| 10 | `precheck_failed` | 启动前预检判定无法运行（**预留**，本期 D8 只告警）| 有 |
| 11 | `pfind_failed` | pFind 主进程 `Searcher.exe` 返回非 0 | 有 |
| 12 | `worker_crashed` | pFind 子进程（`TagSearcher.exe` / `SearchEngine.exe` 等）异常退出 | 有 |
| 13 | `stalled` | 停滞超时（`stall_timeout_min` > 0 时才可能出现）：输出目录无文件更新**且** pFind 进程树 CPU 几乎不增长，判据见 2.8 | 有 |
| 14 | `output_missing` | pFind **没有写出**结果文件（`pFind.protein` / `pFind.spectra` 不存在或 0 字节）。2026-10-07 起，结果文件只有表头、或过滤后 0 条蛋白**不再算失败**，改为告警 `no_hits` | 有 |
| 15 | `internal_error` | DBReducer 自身捕获到异常（如 `std::bad_alloc`），在退出前写好状态 | 有 |
| 16 | `cancelled` | 被用户 / GUI 主动中止（**预留**）| 尽量有 |
| 其他 | —— | **DBReducer 自身没来得及收尾就死了**，见第 4 节 | 可能停在 `running` |

**为什么跳过 2–9**：MinGW 运行时的 `abort()`、未捕获的 C++ 异常都以 **3** 退出（2026-10-06 实测），
系统级崩溃是 `0xC0000005` 这类大数。我们自己的失败码放在 10 以后，读取方看到 3 或大数就能确定是"程序自身崩溃"，
不会和"pFind 失败"混淆。

---

## 2. 状态文件 `dbr_status.json`

### 2.1 位置与生命周期

- 写在 `SearchRes` 目录下（与 `DBReducer.fasta` 同目录）。
- **启动时写一次**（`result=running`），**每次切换阶段、每次心跳都更新**，**结束时写最终结果**。
- 每次都整份重写，不追加。

### 2.2 写入方式（原子化）

1. 写到 `dbr_status.json.tmp`；
2. 用 `MoveFileExW(tmp, json, MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)` 覆盖。
   不能用 C 的 `rename()`：Windows 上目标文件已存在时它直接失败；
3. 覆盖失败（读取方正好打开着文件）时间隔 200 ms 重试 3 次，仍失败就跳过本次。
   **写状态文件失败永远不能导致运行失败。**

读取方：打开、读完、立即关闭；文件不存在或解析失败就当作"暂时没有"，下一轮再读。

### 2.3 编码

- **UTF-8，无 BOM**。
- 来自 cfg 或文件系统的字符串（路径等）是系统 ANSI 编码，写入前要从 `CP_ACP` 转成 UTF-8，
  并对 `\`、`"`、控制字符做 JSON 转义。

### 2.4 字段

| 字段 | 类型 | 何时有 | 说明 |
|---|---|---|---|
| `schema_version` | int | 总是 | 本约定为 1 |
| `tool` | string | 总是 | `"DBReducer"` |
| `tool_build` | string | 总是 | 编译日期 + git 短哈希，如 `"2026-10-08 6e45449"`，用来确认是哪一版产出的 |
| `pid` | int | 总是 | DBReducer 自己的进程号 |
| `result` | string | 总是 | 运行中为 `running`；结束后取第 1 节的名称 |
| `exit_code` | int / null | 结束后 | 与进程退出码一致；运行中为 `null` |
| `stage` | string | 总是 | `SearchRes` 的目录名：`DNARed` / `RNARed` / `Exter` |
| `option` | int | 总是 | 命令行第 2 个参数（0/1/2/3）|
| `phase` | string | 总是 | **当前（或失败时）所处的步骤**，见 2.5 |
| `start_time` / `updated_at` / `end_time` | string | 前两个总是有；`end_time` 结束后才有 | 本地时间，ISO 8601 带时区，如 `2026-10-08T14:02:11+08:00` |
| `elapsed_sec` | int | 总是 | |
| `heartbeat_sec` | int | 总是 | 心跳间隔，来自 cfg 的 `heartbeat_sec`（默认 300）。读取方据此判断"多久没更新算异常" |
| `progress` | object | `phase=search` 起 | `qry_done`、`qry_total`、`newest_file`、`newest_mtime`、`newest_age_sec`。`qry_total` 是谱图文件数（数不出来为 -1）；`qry_done` 是输出目录里的 `*.qry.res` 个数——`resume=0` 时只数本次搜库开始后写出的，`resume=1` 时连上次留下的一起数 |
| `failed_process` | object | 失败码 11/12 | `name`（exe 文件名，保留原大小写，如 `"TagSearcher.exe"`）、`pid`、`exit_code`（`GetExitCodeProcess` 取到的值，8 位十六进制字符串，如 `"0x00000003"`）。注意：内存分配失败导致的 `abort()` 退出码是 **3**；Windows 事件日志里看到的 `0x40000015`（`STATUS_FATAL_APP_EXIT`）是**异常码**，不是进程退出码，两者不要混用 |
| `proteins` | int | 成功后 | `DBReducer.fasta` 的蛋白条数 |
| `params` | object | 总是 | 实际生效的参数：`multip` `thread` `maxprolen` `maxspec` `search_mode` `resume` `stall_timeout_min` `threshold` `ProteinDatabase` `db_bytes` |
| `precheck` | object | 预检之后 | 搜库预检（仅 option 0）：`est_mem_gb`、`warned`、`workers`（= min(multip, 谱图文件数)）、`spectra_batch`（= min(maxspec, 单个文件谱图数)）、`per_worker_gb`；BuildDict 估算（**仅 option 3**，超没超都记；D9 之后其余 option 流式读库、不估算也不写这两项）：`builddict_est_gb`、`builddict_warned`；以及 `commit_avail_gb` |
| `warnings` | array | 有告警时 | 每条 `{ "what": "...", "info": "key=value ..." }`，与 `[DBR:WARN]` 日志行一一对应（`info` 即日志行 `what` 之后的部分），最多保留最近 50 条 |
| `detail` | string | 失败时 | 一句话说明发生了什么（中文，供人看）|
| `suggestion_code` | string | 失败或告警时 | 建议的固定代码，见 2.6 |
| `suggestion` | string | 同上 | `suggestion_code` 对应的中文文案，GUI 可直接显示 |

### 2.5 `phase` 取值

按顺序推进，**失败时停在出错的那一步**，这就是"在哪一步终止"：

| `phase` | 含义 | 对应代码 |
|---|---|---|
| `init` | 读 cfg、建输出目录 | `main.cpp` 开头 |
| `precheck` | 资源预检 | D8 |
| `search` | pFind 搜库（`Searcher.exe` 及其子进程在跑）| `main.cpp:73`、`:160` |
| `verify_search` | 检查 pFind 产物 | D6 |
| `rebuild` | 读原始库、读 `.protein` / `.spectra`、挑蛋白 | `Rebuild.cpp`、`Rebuild-spec.cpp` |
| `write` | 写出并校验 `DBReducer.fasta` | |
| `done` | 成功结束。**失败时不会出现 `done`**，`phase` 保持在出错的那一步 | |

`option` 为 1/2/3 时不搜库，直接从 `init` 跳到 `rebuild`。

### 2.6 `suggestion_code`

DBReducer 和 GUI 共用这张表。日志行里只出现代码，中文文案只写在状态文件里（原因见 3.2）。

| 代码 | 文案 |
|---|---|
| `reduce_multip` | 预估内存超出可提交余额，建议降低 multip 或将库分片 |
| `reduce_maxprolen` | 单个进程的蛋白分批过大，建议降低 maxprolen |
| `reduce_maxspec` | 单个进程的谱图批次过大，建议降低 maxspec |
| `free_memory` | DBReducer 读取原始库时内存不足，请关闭其他占用内存的程序或换用内存更大的机器 |
| `shard_db` | 翻译库过大，缩库最后读取原始库时可能内存不足，建议将库分片（D9 之后只有 option 3 会给出）|
| `check_pfind_license` | pFind 没有产出任何结果就结束了，很可能是 pFind 已过期（按系统日期判断），请更新 pFind 或检查系统日期 |
| `abort_reduce_multip` | 进程以退出码 3 中止（abort），可能是内存不足，也可能是程序内部错误：可先降低 multip 重跑；若仍失败，请保留日志与状态文件联系开发者 |
| `check_input` | 请检查 cfg 中的路径与输入文件是否存在 |
| `check_pfind` | 请检查 pFind 安装是否完整（thirdparty\pFind3\bin）|
| `rerun_resume` | 可设置 resume=1 后重跑，已完成的文件会被跳过 |
| `see_pfind_log` | 请查看 SearchRes 下 pFind 的日志定位原因 |
| `report_bug` | 程序内部错误，请保留日志与 dbr_status.json 联系开发者 |

`check_pfind_license` 的判断：`Searcher.exe` 返回 0，但结果文件不存在、输出目录里没有任何 `*.qry.res`、且搜库耗时不到 5 分钟。

`abort_reduce_multip`：pFind 子进程以 3 退出时给出。MinGW 的 `abort()` 与未捕获异常都是 3，内存分配失败只是可能原因之一，所以文案不下定论。

### 2.7 示例

```json
{
  "schema_version": 1,
  "tool": "DBReducer",
  "tool_build": "2026-10-08 6e45449",
  "pid": 18244,
  "result": "worker_crashed",
  "exit_code": 12,
  "stage": "DNARed",
  "option": 0,
  "phase": "search",
  "start_time": "2026-10-08T14:02:11+08:00",
  "updated_at": "2026-10-08T17:27:40+08:00",
  "end_time": "2026-10-08T17:27:40+08:00",
  "elapsed_sec": 12329,
  "heartbeat_sec": 300,
  "progress": { "qry_done": 5, "qry_total": 36, "newest_file": "test.F6.L1.9.ot", "newest_mtime": "2026-10-08T17:20:03+08:00" },
  "failed_process": { "name": "TagSearcher.exe", "pid": 23552, "exit_code": "0x00000003" },
  "params": { "multip": 6, "thread": 2, "maxprolen": 600000000, "maxspec": 10000, "search_mode": -1,
              "resume": 0, "stall_timeout_min": 0, "threshold": 0.1,
              "ProteinDatabase": "D:\\data\\trans_DNA_dedup.fasta", "db_bytes": 20830000000 },
  "precheck": { "est_mem_gb": 174.0, "warned": true, "workers": 1, "spectra_batch": 70000, "per_worker_gb": 174.0,
                "builddict_est_gb": 87.3, "builddict_warned": false, "commit_avail_gb": 118.5 },
  "detail": "pFind 子进程 TagSearcher.exe 异常退出（退出码 0x00000003），可能是内存不足，也可能是程序内部错误，已终止整个 pFind 进程树",
  "suggestion_code": "abort_reduce_multip",
  "suggestion": "进程以退出码 3 中止（abort），可能是内存不足，也可能是程序内部错误：可先降低 multip 重跑；若仍失败，请保留日志与状态文件联系开发者"
}
```

### 2.8 停滞判据（13）

`stall_timeout_min` > 0 时，以下两条**同时**持续 `stall_timeout_min` 分钟才判停滞、终止 pFind 进程树：

1. `SearchRes` 目录（含一层子目录）没有任何文件更新（看最新的 mtime）；
2. pFind 整个进程树的累计 CPU 时间（Job 的 `TotalUserTime + TotalKernelTime`，含已退出的进程）几乎不增长：按窗口结算，窗口内增量 < 单核的 5% 记为空闲。窗口长度取阈值的 1/4，限制在 1–10 秒。

只看第 1 条会误判：`ReRanker` 等阶段会满负荷计算很久而不写文件（测试方 Q11a：71 秒无文件更新）。阈值不取 0：主进程轮询等待 worker 时也有零星 CPU 消耗，worker 挂起时仍应判停滞。
不能把进程放进 Job 时（`job_unavailable` 告警），取不到 CPU 时间，退回只看第 1 条。
已知局限：死循环空转的进程 CPU 一直在涨，不会被判停滞。

---

## 3. 日志行

### 3.1 格式

```
[DBR:<事件>] key=value key=value …
```

- 每个事件独占一行，打印后立即 `fflush(stdout)`。
- **每个前缀行前面先输出一个换行**（与该行在同一次写出里），保证它总在行首、可以按行匹配。原因：pFind 刷新进度只输出 `\r` 不换行，不加的话前缀行会接在进度后面（测试方 R_Q11a：第二轮 10 条心跳有 5 条粘在 `[pFind] <Rerank>: 30%` 之后）。代价是输出里多出空行，读取方忽略空行即可。
- 读取方切行时建议把 `\r` 也当行分隔符（pFind 的进度行以 `\r` 结尾）。
- **只用 ASCII**。值里不含空格；确实需要空格时用双引号括起，内部的 `"` 写成 `\"`。
- 时长一律给秒数（`elapsed_sec=`），需要人读的形式由 GUI 自己换算。

### 3.2 为什么日志行只用 ASCII

GUI 按 UTF-8 解码子进程输出（`gui/app/core/runner.py:115`），而 DBReducer 打印的路径是系统 ANSI 编码
（中文 Windows 上是 GBK），pFind 的输出也是 ANSI。中文路径或中文文案一旦进入日志行就会变成乱码。
所以日志行里只放代码和数字，中文只出现在 UTF-8 的状态文件里。

### 3.3 事件

| 事件 | 何时 | 必带的 key |
|---|---|---|
| `START` | 启动、读完 cfg 后 | `stage` `option` `pid` `db_bytes` `multip` `thread` `maxprolen` `maxspec` `resume` `stall_timeout_min` `heartbeat_sec` |
| `PHASE` | 每次切换 `phase` | `phase` `elapsed_sec` |
| `WARN` | 告警，不影响运行 | `what`，其余 key 视 `what` 而定，见下表 |
| `CLEANUP` | 搜库前删掉了残留的 `<id>.complete`（只在删了文件时打印）| `removed_complete` |
| `RESUME` | `resume=1` 时清理残缺文件之后 | `skip`（已完成、将被跳过的文件数）`removed_incomplete`，有删除时加 `file`（逗号分隔）|
| `HEARTBEAT` | `search` 阶段每 `heartbeat_sec` 一次 | `elapsed_sec` `qry_done` `qry_total` `newest_age_sec`，有文件时加 `newest`；取得到进程树 CPU 时间时加 `job_cpu_sec`（累计）`cpu_idle_sec`（CPU 几乎不增长已持续的秒数）|
| `ERROR` | 判定失败、退出前 | `result` `exit_code` `phase`，失败进程相关时加 `process` `pid` `proc_exit` |
| `DONE` | 最后一行，不论成败 | `result` `exit_code` `elapsed_sec`，成功时加 `proteins` |

`WARN` 的 `what` 取值：

| `what` | 含义 | 其余 key |
|---|---|---|
| `memory` | 预检估算的内存超出可提交余额 | `est_mem_gb` `commit_avail_gb` `per_worker_gb` `multip` `workers` `spectra_batch` `suggestion_code`；单个 worker 装得下时另有 `suggest_multip`（≥1），装不下时不给 |
| `precheck_skipped` | 没能完成预检（读不了库等）| — |
| `child_exit` | 名单外的 pFind 子进程非 0 退出（不终止运行）| `process` `pid` `proc_exit` |
| `job_unavailable` | 没能把 pFind 放进 Job Object，子进程崩溃检测失效 | — |
| `spectra_bad_lines` | `.spectra` 里有无法解析的行，已跳过 | `count` `in` |
| `no_hits` | pFind 正常跑完，但没有任何可保留的蛋白（结果文件只有表头，或过滤后 0 条）。**不算失败**：以 0 退出，`proteins=0`，写出空的 `DBReducer.fasta` | `reason`（`header_only` / `filtered_out`）|
| `builddict_memory` | 预检（**仅 option 3**）：缩库最后整库读取原始库（`BuildDict`）的内存估算超出可提交余额。系数为**估算、未实测**：库文件大小 × 4.5。其余 option 自 D9 起流式读库，不再出现 | `est_gb` `commit_avail_gb` `suggestion_code=shard_db` |

示例：

```
[DBR:START] stage=DNARed option=0 pid=18244 db_bytes=20830000000 multip=1 thread=4 maxprolen=50000000 maxspec=70000 resume=1 stall_timeout_min=0
[DBR:PHASE] phase=precheck elapsed_sec=1
[DBR:WARN] what=memory est_mem_gb=174.0 commit_avail_gb=118.5 suggestion_code=reduce_multip
[DBR:RESUME] skip=5 removed_incomplete=1
[DBR:PHASE] phase=search elapsed_sec=3
[DBR:HEARTBEAT] elapsed_sec=12300 qry_done=5 qry_total=36 newest_age_sec=420 job_cpu_sec=98400 cpu_idle_sec=0
[DBR:ERROR] result=worker_crashed exit_code=12 phase=search process=TagSearcher.exe pid=23552 proc_exit=0x00000003
[DBR:DONE] result=worker_crashed exit_code=12 elapsed_sec=12329
```

---

## 4. 读取方的判定规则

pAnno 和 GUI 按下表解读，**以退出码为准，状态文件补充细节**：

| 情况 | 判定 | 显示 / 处理 |
|---|---|---|
| 退出码 0 且 `result=ok` | 成功 | 继续流程 |
| 退出码在第 1 节表内（1、10–16）| 按 `result` 归类 | 显示 `phase`、`detail`、`suggestion`；pAnno 停止流程 |
| 退出码 3，或 ≥ `0xC0000000`，或不在表内 | **DBReducer 自身崩溃** | 状态文件若停在 `running`，其中的 `phase` 就是崩溃时所处的步骤；建议 `report_bug` |
| 进程还在，但 `updated_at` 已超过 `3 × heartbeat_sec` 没更新 | **疑似挂死** | 只提示，不自动终止（是否终止由 `stall_timeout_min` 决定）|
| 状态文件根本不存在 | **没能正常启动** | 多半是路径或 exe 本身的问题；建议 `check_input` / `check_pfind` |

---

## 5. 对实现的要求（D5）

- **所有退出路径**都要先写最终状态再退出，包括 `utils` 里 `ReadFilebyStream` 等函数中现有的 `exit(1)`
  （可改为抛异常，由 `main` 统一收尾）。
- `main` 最外层 `try/catch`：`std::bad_alloc` 和其他 `std::exception` 都归为 `internal_error`（15），
  `detail` 写入 `what()`。
- 建目录失败目前是 `return 0`（`main.cpp:55`、`:146`），改为 `param_error`（1）。
- 进程托管模块（D2）要能把失败子进程的名字、PID、退出码交给状态文件。
- 判定为"pFind 子进程崩溃"（12）的名单：`TagSearcher.exe`、`SearchEngine.exe`、`PSSearcher.exe`、`ReRanker.exe`、`Group.exe`。
  依据：pFind 2021 版源码中这些程序正常结束都返回 0，非 0 只出现在出错分支；且 SearchEngine 串行调用
  TagSearcher → ReRanker → PSSearcher → ReRanker，不检查它们的返回值，任何一个崩了都会静默产出残缺结果。
  名单外的子进程（cmd.exe、Summary.exe、ACGen.exe 等）非 0 退出只打 `WARN what=child_exit`。

## 6. 本期不定的内容

- pAnno 自己的状态文件（第 2 步）：沿用本约定，增加"流水线阶段"一层，届时再补字段。
- `cancelled`（16）的触发方式：GUI 如何通知 DBReducer 中止，留到第 3 步。
- `precheck_failed`（10）：等内存系数实测标定后再启用。
