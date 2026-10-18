/**
 * @file Status.hpp
 * @brief 运行状态对外报告：退出码、状态文件 dbr_status.json、[DBR:事件] 日志行。
 *        约定见仓库根目录 docs/status-interface.md（schema_version 1）。
 *
 * 字符串约定：本模块内部一律保存 UTF-8。来自 cfg / 文件系统的路径是系统 ANSI 编码，
 * 传进来之前用 AnsiToUtf8 转换；源码里的中文字面量本身就是 UTF-8。
 */
#ifndef DBR_STATUS_HPP
#define DBR_STATUS_HPP

#include <stdexcept>
#include <string>

namespace dbr {

enum class Result {
    running, ok, param_error, precheck_failed, pfind_failed,
    worker_crashed, stalled, output_missing, internal_error, cancelled
};
const char *ResultName(Result r);
int ResultExitCode(Result r);

/// 失败的子进程（状态文件的 failed_process）
struct ProcInfo {
    bool valid = false;
    std::string name;
    unsigned long pid = 0;
    unsigned long exit_code = 0;
};

/// 由 main 统一收尾的失败。detail 必须是 UTF-8。
class DbrError : public std::runtime_error {
public:
    DbrError(Result r, const std::string &detail_utf8, const std::string &suggestion_code = "",
             const ProcInfo &proc = ProcInfo())
        : std::runtime_error(detail_utf8), result(r), suggestion_code(suggestion_code), proc(proc) {}
    Result result;
    std::string suggestion_code;
    ProcInfo proc;
};

/// search 阶段的进度（状态文件的 progress）
struct Progress {
    int qry_done = -1;
    int qry_total = -1;
    std::string newest_file;     // UTF-8
    std::string newest_mtime;    // ISO 8601
    long long newest_age_sec = -1;
};

std::string AnsiToUtf8(const std::string &ansi);
std::string NowIso();
std::string FileTimeIso(unsigned long long filetime);
std::string Hex8(unsigned long v);   // "0x0000000C"

/// 打印一行 "[DBR:<event>] <kv>" 并立即 flush。kv 只允许 ASCII。
void Log(const char *event, const std::string &kv);
/// 告警：打印 "[DBR:WARN] what=<what> <kv>"，同时记入状态文件的 warnings 数组并刷新状态文件。
void Warn(const std::string &what, const std::string &kv = "");

/// 启动时调用一次。search_res 为 ANSI 路径（状态文件写在它下面）。
void StatusInit(const std::string &search_res_ansi, int option, int heartbeat_sec);
/// 记录实际生效的参数（状态文件 params）。value_json 为已经是 JSON 的值。
void StatusParamRaw(const std::string &key, const std::string &value_json);
void StatusParamNum(const std::string &key, double v);
void StatusParamStr(const std::string &key, const std::string &value_utf8);
/// 预检结果（状态文件 precheck）
void StatusPrecheck(double est_mem_gb, double commit_avail_gb, bool warned, int workers, long long spectra_batch,
                    double per_worker_gb);
/// BuildDict 内存估算（状态文件 precheck.builddict_*），仅 option 3（D9 之后其余 option 流式读库）
void StatusBuildDict(double est_gb, double commit_avail_gb, bool warned);
/// 切换阶段：打印 [DBR:PHASE] 并写状态文件
void StatusPhase(const char *phase);
const char *CurrentPhase();
/// 心跳：更新进度并写状态文件（日志行由调用方打印）
void StatusHeartbeat(const Progress &p);
/// 成功时写出的蛋白条数（[DBR:DONE] 的 proteins）
void StatusSetProteins(long long n);
/// 不改变结果、只附加告警建议（状态文件 suggestion_code / suggestion）
void StatusSuggest(const std::string &code);
long long ElapsedSec();

/// 收尾：写最终状态、打印 [DBR:ERROR]（失败时）与 [DBR:DONE]，返回进程退出码。
int StatusFinish(Result r, const std::string &detail_utf8 = "", const std::string &suggestion_code = "",
                 const ProcInfo *proc = nullptr);

}  // namespace dbr

#endif
