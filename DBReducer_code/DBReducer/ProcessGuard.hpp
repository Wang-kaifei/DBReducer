/**
 * @file ProcessGuard.hpp
 * @brief 启动一个外部程序并托管它的整棵进程树（Windows Job Object）。
 *
 *  - 子进程及其派生的所有进程（包括经 cmd / start /b 拉起的）都被收进同一个 Job；
 *  - Job 设了 KILL_ON_JOB_CLOSE：托管方自己退出（哪怕被强杀）时，整棵树随之终止，不留孤儿；
 *  - 名单内的子进程以非 0 退出 → 立即终止整棵树并报告是谁、退出码多少；
 *  - 每 tick 回调一次（附带整棵树的累计 CPU 时间），调用方在里面做心跳 / 停滞判断，返回 true 即终止整棵树。
 *
 * 本模块不依赖 DBReducer 的其他代码，pAnno 托管 pFind 两轮搜库时可直接复制使用。
 */
#ifndef PROCESS_GUARD_HPP
#define PROCESS_GUARD_HPP

#include <functional>
#include <string>
#include <vector>

namespace procguard {

struct ChildExit {
    std::string name;           // exe 文件名，保留原大小写，如 "TagSearcher.exe"
    unsigned long pid = 0;
    unsigned long exit_code = 0;
};

enum class Outcome {
    Finished,       ///< 主进程自己结束了，退出码见 main_exit_code
    ChildFailed,    ///< 名单内的子进程以非 0 退出，已终止整棵树，详情见 failed
    Aborted,        ///< on_tick 要求终止，已终止整棵树
    LaunchFailed    ///< CreateProcess 失败，错误码见 launch_error
};

struct Result {
    Outcome outcome = Outcome::Finished;
    unsigned long main_pid = 0;
    unsigned long main_exit_code = 0;
    ChildExit failed;
    unsigned long launch_error = 0;
    bool job_ok = true;         ///< false 表示没能把进程放进 Job，只能监视主进程（子进程崩溃检测失效）
};

/// 每次 tick 回调时提供的信息
struct TickInfo {
    bool cpu_valid = false;              ///< 能否取到 Job 的 CPU 时间（进程没能放进 Job 时为 false）
    unsigned long long job_cpu_100ns = 0;///< Job 内所有进程（含已退出的）累计 CPU 时间（用户态+内核态，100 ns 为单位）
};

struct Options {
    std::string exe_path;                     ///< 可执行文件完整路径（ANSI）
    std::string args;                         ///< 参数串，不含 exe 本身；需要时自己加引号
    std::string work_dir;                     ///< 工作目录，空则继承
    std::vector<std::string> fatal_children;  ///< 小写 exe 名（比较时不分大小写）；这些子进程非 0 退出即判定失败
    unsigned tick_ms = 1000;                  ///< 轮询与回调间隔
    std::function<bool(const TickInfo &)> on_tick;  ///< 返回 true → 终止整棵树，Outcome::Aborted
    std::function<void(const ChildExit &)> on_other_child_failed;  ///< 名单外的子进程非 0 退出（只通知）
};

Result Run(const Options &opt);

/// 由 Job 的 CPU 时间判断进程树是否在"干活"（供停滞判断用）。
/// 每 window_ms 结算一次：窗口内 CPU 增量 ≥ 单核的 busy_frac 即视为活跃。
/// 阈值不取 0：主进程轮询等待 worker 时也有零星 CPU 消耗，不能因此判为活跃。
class CpuActivity {
public:
    explicit CpuActivity(unsigned long long now_ms, unsigned window_ms = 10000, double busy_frac = 0.05)
        : window_ms_(window_ms), busy_frac_(busy_frac), last_active_ms_(now_ms) {}
    void Feed(const TickInfo &ti, unsigned long long now_ms);  ///< 每次 on_tick 调用
    bool valid() const { return have_base_; }                  ///< 取到过 CPU 时间
    double cpu_sec() const { return cpu_100ns_ / 1e7; }        ///< 最近一次读到的累计 CPU 秒数
    unsigned long long idle_ms(unsigned long long now_ms) const {  ///< 距最近一个活跃窗口结束
        return now_ms > last_active_ms_ ? now_ms - last_active_ms_ : 0;
    }

private:
    unsigned window_ms_;
    double busy_frac_;
    unsigned long long last_active_ms_;
    bool have_base_ = false;
    unsigned long long base_cpu_ = 0, base_ms_ = 0, cpu_100ns_ = 0;
};

}  // namespace procguard

#endif
