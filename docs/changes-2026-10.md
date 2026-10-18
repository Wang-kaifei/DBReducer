# 2026-10 改动说明（李治紧急更新：D0–D9）

基线：HEAD `6e45449`。未提交，由用户确认后手动提交。
接口约定见 [status-interface.md](status-interface.md)。

## 改了什么

| 编号 | 内容 | 文件 |
|---|---|---|
| D0 | CMake 不再链接写死的他人目录 `D:/DBReducer/...`，改为把本仓库 `utils` 作为子工程编译 | `DBReducer/CMakeLists.txt` |
| D1 | 新 cfg key `resume`（默认 0）。为 1 时调用 `Searcher.exe <cfg> 1`，pFind 跳过已有 `qry.res` 的谱图文件；启动前删掉写了一半的 `qry.res`（空文件或末字节不是换行）| `main.cpp` |
| D2 | `system("Searcher.exe ...")` 换成进程托管：Job Object 收纳整棵进程树（`KILL_ON_JOB_CLOSE`，DBReducer 退出 / 被杀时 pFind 全部随之终止）；`TagSearcher` / `SearchEngine` / `PSSearcher` / `ReRanker` / `Group` 任一非 0 退出即终止整棵树、退出码 12 | `ProcessGuard.hpp/.cpp`（可复用，pAnno 第 2 步直接拷）、`main.cpp` |
| D3 | 心跳：新 cfg key `heartbeat_sec`（默认 300），打印 `[DBR:HEARTBEAT]` 并刷新状态文件 | `main.cpp` |
| D4 | 停滞超时：新 cfg key `stall_timeout_min`（默认 0 = 关，可写小数），输出目录无文件更新且 pFind 进程树 CPU 几乎不增长、两者都持续 N 分钟即终止，退出码 13（第 2 步测试后由"只看文件"改来，见下文） | `main.cpp` |
| D5 | 退出码 / `dbr_status.json` / `[DBR:事件]` 日志行。所有退出路径都先写最终状态：`utils/IOUtil.h` 的 `exit(1)` 改为抛 `FileError`；`main` 外层统一 `try/catch`，`bad_alloc` 等归为 15 | `Status.hpp/.cpp`、`main.cpp`、`utils/IOUtil.h` |
| D6 | 产物校验：pFind 结束后 `pFind.protein`（或 `.spectra`）必须存在且非空；`DBReducer.fasta` 必须至少 1 条。否则退出码 14 | `main.cpp` |
| D7 | 读 `.spectra` 时检查列数、`stof` 加异常保护。区分"无鉴定谱图"（6 列，pFind 正常输出，静默跳过）与真正的残缺行（计数后 `WARN`）| `Rebuild-spec.cpp` |
| D8 | 启动 pFind 前资源预检，只告警不改参数。实际 worker 数取 min(multip, 谱图文件数)；谱图侧取 min(maxspec, 单个文件谱图数)，后者在最大的谱图文件里均匀取 8 段采样外推（yeast −7%、梨 +14%）| `Precheck.hpp/.cpp`、`main.cpp` |

顺带：

- 搜库前总是删掉残留的 `<id>.complete`：`multip>1` 时 Searcher 靠它判断 worker 是否结束，上次中途被杀留下的会让本次不等对应 worker 就进入蛋白推断。
- `Params` 里原先没有初始值的 `multip`/`thread`/`maxspec`/`maxprolen`/`msmsnum` 等给了默认值；`len_lower`/`len_upper` 写反的默认值改正（pAnno 总会写这些键，实际运行不受影响）。
- 建目录失败原先 `return 0`，现为退出码 1；pFind 目录缺 `modification.ini` 原先实际以 1 退出，现在给出明确原因（仍为 1）。
- 结束时的耗时打印原先把 `duration` 对象直接传给 `%f`（未定义行为），已改正。
- 新增 `build.ps1`：不依赖 make，直接用 msys2 g++ 编译。
- 状态文件新增 `warnings` 数组，所有 `[DBR:WARN]` 同步记入。
- `verify_search` 阶段除了"存在且非空"，还要求表头之后至少有一行数据（`.protein` 2 行表头、`.spectra` 1 行），只有表头时在此报 14。

## 第 2 步追加（2026-10-07，配合 pAnno 流程打通）

1. **0 命中不再算失败**：pFind 结果文件只有表头，或过滤后 0 条蛋白 → 以 0 退出、`proteins=0`、`WARN what=no_hits reason=header_only|filtered_out`，写出空的 `DBReducer.fasta`；只有表头时直接跳过 `BuildDict`。结果文件**不存在**仍为 14。option 3 在 0 命中时不再提取 top2。
2. **D8b**：所有 option 在 init 之后估算 `BuildDict` 内存（库大小 × 4.5，按数据结构推算、**未实测**），超出可提交余额时 `WARN what=builddict_memory`，建议码 `shard_db`。
3. **pFind 疑似过期**：本次搜过库、结果文件缺失、输出目录无任何 `*.qry.res`、搜库不到 5 分钟 → 建议码 `check_pfind_license`。
4. `build.ps1`：`-OutDir` 支持绝对路径。

## 第一轮测试后的修正（2026-10-07）

测试方（user-questions 会话）在 yeast + 2026 版 pFind 上的结果：正常运行无误判（各子进程退出码均为 0）、死锁场景 1 秒内退出码 12、续跑与残缺 qry.res、停滞、残缺 spectra、结果缺失、接入 pAnno 均通过。据此修正：

1. 正常结束时 `progress` 没刷新（只在心跳时更新）→ pFind 结束后再扫一次输出目录。
2. 预检严重高估（yeast、multip=64、maxspec=500 万时估出 18 TB）→ 按谱图文件数封顶 worker 数、按单个文件谱图数封顶谱图侧；`suggest_multip` 不再出现 0。
3. （可选项，已做）`warnings` 数组；只有表头的结果文件提前在 `verify_search` 报错。

## 第 2 步测试后的修正（2026-10-07）

1. **停滞误判（测试方 Q11a）**：ReRanker 满负荷 71 秒不写文件被判停滞。改为"输出目录无文件更新**且** pFind 进程树 CPU 几乎不增长"同时持续 N 分钟才算（判据见 status-interface.md 2.8）。`ProcessGuard` 每 tick 附带 Job 的累计 CPU 时间（`TickInfo`），新增 `CpuActivity` 做窗口结算；`HEARTBEAT` 行加 `job_cpu_sec` `cpu_idle_sec`。冒充程序验证：满负荷 12 s / 45 s（阈值 3 s / 30 s）均正常完成；worker 挂起、主进程每 50 ms 轮询，3 s / 30 s 准时判 13；原有 ok/crash/stall/fail/killdbr 不变，均无残留进程。
2. **退出码 3 不下定论**：pFind 子进程以 3 退出，原文案"多为内存分配失败后 abort"+`reduce_multip`，改为"可能是内存不足，也可能是程序内部错误"+新建议码 `abort_reduce_multip`（降低 multip 重跑；仍失败则保留日志联系开发者）。
3. **BuildDict 估算总是记录**：状态文件 `precheck` 新增 `builddict_est_gb` `builddict_warned`，所有 option 都写（原先只在超出时有 WARN）。
4. 修正：`SuggestionText` 漏了 `shard_db`、`check_pfind_license`，状态文件里 `suggestion` 文案为空，已补。
5. **`[DBR:*]` 行粘在 pFind 进度后面**（测试方 R_Q11a）：pFind 刷进度只输出 `\r`，前缀行接在同一行，按行匹配会漏。现在每个前缀行前先输出一个换行（同一次 printf）。冒充 pFind 以 `\r` 刷进度验证：旧版 14 行有 6 行不在行首，新版 0 行。

## D9 / I-06：缩库最后一步流式读库（2026-10-16）

原先缩库最后一步 `BuildDict` 把**整个原始翻译库**读进两张哈希表（名字→序列、名字→完整表头），再按 pFind 结果取出要写的蛋白。
梨 DNA 库 2.7 GB 时这一步提交量涨到十几 GB；李治那种二十多 GB 的库在这一步就可能内存不足。

现在先确定要写出哪些蛋白，再边扫 fasta 边只保留这些蛋白（`Rebuild.cpp::LoadSelected`，解析规则与 `BuildDict` 逐行相同，同样 first-wins：
完整表头取该名字第一次出现的表头，序列取第一个非空序列）。

- `.protein` 路径（option 2、option 0 且 `output_pro_info=1`，即 pAnno 的 RNA/DNA 缩库）：先 `ReadProteinName` 再 `LoadSelected`。
  输出按名字排序（`std::set`）写出，与库的顺序无关 → **逐字节不变**。
- `.spectra` 路径（option 1、option 0 且 `output_pro_info=0`，即 pAnno 的外部库缩库）：先 `GroupPros`（与原始库无关），再只读入非 subset 蛋白。
  按 pros 的顺序写出 → 与改前相同（该路径改前就因 OpenMP 顺序不固定，见下）。
- option 3：后面的 top2 要扫全库，**仍整库读入**，D8b 估算与告警保留。
- 其余 option 不再做 D8b 估算：日志改为 `Rebuild reads the database by streaming ...`，状态文件 `precheck.builddict_*` 只在 option 3 出现；
  `fasta size: N` 换成 `fasta scanned: X entries, selected Y proteins, loaded Z`。

验证（`d:\npc\test\panno\d9_stream\`，新 769b2771 / 旧 c5cb8551，同一份 pFind 结果）：

| 用例 | 结果 | 提交峰值 新 / 旧 |
|---|---|---|
| yeast RNA / DNA option 2 | 逐字节相同，且与 Q0_full 原目录相同 | 4 MB / 13 MB、4 MB / 41 MB |
| yeast RNA / DNA option 1 | 记录集合相同（703 / 697 条）；新旧各连跑 5 次，**旧版 5 次 5 个 md5**——顺序本来就不固定 | — |
| yeast DNA option 3 | DBReducer.fasta 记录集合相同，top2.fasta 记录集合相同，top2.pep 逐字节相同 | 两版相同（都整库读）|
| 边界用例（首个表头前的游离行、同名两次、同名首次空序列、`\|` 截名、制表符表头、多行序列、空行、CRLF、末行无换行、subset 蛋白）| option 2 逐字节相同；option 1 记录集合相同、subset 蛋白不写出 | — |
| 梨 RNA option 2（库 53 MB）| 逐字节相同，且与 Pyrus_released2 原目录相同 | 6 MB / 181 MB |
| 梨 DNA option 2（库 2.7 GB，2167 万条，选 2849 条）| 逐字节相同，且与 Pyrus_released2 原目录相同；耗时 18.5 s / 205 s | **6 MB / 12.22 GB**（物理内存峰值 9 MB / 11.83 GB）|
| 梨 DNA option 1（.spectra 路径，选 8844 条）| 记录集合相同；耗时 35 s / 146 s | 0.50 GB / 12.71 GB（新版这 0.5 GB 在 yeast option 1 上同样有，与库大小无关）|

提交峰值取自操作系统记录的进程峰值（`PeakPagefileUsage`，即任务管理器"提交大小"的峰值）。梨 DNA 旧版实测 12.22 GB ÷ 库 2.53 GB ≈ **4.8 倍**，
D8b 原先估的 4.5 倍略低（option 3 仍用这一估算，系数未改）。

**最坏情况：名单 = 全库**（`make_worst.py` 给库里每个名字造一行 pFind.protein）。名单最多等于库的条目数（名单里的名字必须在库里，否则新旧都报错）。
这时新版内存里留下的就是旧版的全部内容，理论上与旧版持平。第一版（769b2771）把名单另复制成 `unordered_set` 查找，实测比旧版高约 20%
（梨 RNA 26 万条：0.282 / 0.233 GB；yeast DNA 6.7 万条：0.060 / 0.051 GB）。改为直接在原名单上查（`.protein` 路径查 `ReadProteinName` 的 `std::set`；
`.spectra` 路径的查找表只存指向 pros 里名字的 `string_view`）后（**eed8eba0**）：

| 用例 | 结果 | 提交峰值 新 / 旧 | 耗时 新 / 旧 |
|---|---|---|---|
| 最坏：梨 RNA 全库 26 万条 | 逐字节相同 | 0.232 / 0.233 GB | 3.5 / 4.0 s |
| 最坏：yeast DNA 全库 6.7 万条 | 逐字节相同 | 0.051 / 0.051 GB | 1.1 / 1.5 s |
| 梨 DNA option 2 | 逐字节相同（e7d29162）| 5 MB | 18.0–18.5 s（769b2771 为 17.0 s：std::set 查找略慢，可忽略）|
| 梨 DNA option 1、yeast option 1/3、边界用例 | 记录集合相同 | — | — |

即最坏情况与旧版持平，不会变成负优化。（另：复制 20 GB 文件的同时跑梨 DNA option 2 曾用 317 s，复制结束后重测为 18 s，是磁盘争用，与程序无关。）
**交付版为 eed8eba0**（`dist\d9b\`），已放入 overlay 与 `d9_stream\pAnnobin`；769b2771（`dist\d9\`）仅作记录。

## 已做的验证（2026-10-06，本机）

- **基线**：未修改源码 + 本机 g++ 12.2.0 + 发布版 `libUtils.dll`，编出的 `DBReducer.exe` 与发布版只差 6 字节（PE 时间戳与校验和）。
- **CMake（D0）**：`cmake -G "MinGW Makefiles"` 能独立编出，产物行为与 `build.ps1` 一致。
- **缩库结果不变**：梨 RNARed 数据，option 2 输出与基线、与原始运行三者 MD5 相同；option 1 输出与基线条目集合相同（该路径原本就因 OpenMP 而顺序不固定，基线自己两次运行 MD5 也不同）。
- **D7**：`.spectra` 末尾加空行与半行、阈值放到读到文件尾 → 基线 `bad_alloc` 崩溃（退出码 3），修改版正常完成并报告 1 条残缺行。
- **进程托管（用冒充的 Searcher/TagSearcher）**：正常完成、worker 以 3 退出、停滞、主进程返回非 0、续跑、预检告警、强杀 DBReducer 七个场景，退出码与状态文件均符合约定，均无残留进程。**真实 pFind 上的验收由测试方完成。**
