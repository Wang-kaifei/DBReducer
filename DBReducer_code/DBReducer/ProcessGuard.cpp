/**
 * @file ProcessGuard.cpp
 * @brief 见 ProcessGuard.hpp。
 *
 * 发现子进程用两条路：Job 的完成端口消息（JOB_OBJECT_MSG_NEW_PROCESS，及时）+ 每 tick 列一次
 * Job 内进程（补漏，完成端口消息不保证送达）。发现后立刻 OpenProcess 持有句柄，
 * 这样子进程退出后仍能用 GetExitCodeProcess 取到退出码。
 */
#ifdef _WIN32
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <sys/ptrace.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "ProcessGuard.hpp"

#include <algorithm>
#include <map>

namespace procguard {

namespace {

std::string BaseName(const std::string &path) {
    size_t pos = path.find_last_of("\\/");
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

std::string Lower(std::string n) {
    std::transform(n.begin(), n.end(), n.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return n;
}

}  // namespace

#ifdef _WIN32
// ================================================================ Windows：Job Object

namespace {

struct Tracked {
    std::string name;
    HANDLE handle = nullptr;
};

/// 复制一份可继承的标准句柄给子进程。复制失败就原样传（与 system() 的行为一致）。
HANDLE InheritableStd(DWORD which, bool *dup) {
    *dup = false;
    HANDLE h = GetStdHandle(which);
    if (h == nullptr || h == INVALID_HANDLE_VALUE) return h;
    HANDLE d = nullptr;
    if (DuplicateHandle(GetCurrentProcess(), h, GetCurrentProcess(), &d, 0, TRUE, DUPLICATE_SAME_ACCESS)) {
        *dup = true;
        return d;
    }
    return h;
}

}  // namespace

Result Run(const Options &opt) {
    Result res;

    // ---- Job + 完成端口 ----
    HANDLE job = CreateJobObjectA(nullptr, nullptr);
    HANDLE port = nullptr;
    if (job) {
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION li;
        ZeroMemory(&li, sizeof(li));
        li.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
        SetInformationJobObject(job, JobObjectExtendedLimitInformation, &li, sizeof(li));
        port = CreateIoCompletionPort(INVALID_HANDLE_VALUE, nullptr, 0, 1);
        if (port) {
            JOBOBJECT_ASSOCIATE_COMPLETION_PORT acp;
            acp.CompletionKey = job;
            acp.CompletionPort = port;
            if (!SetInformationJobObject(job, JobObjectAssociateCompletionPortInformation, &acp, sizeof(acp))) {
                CloseHandle(port);
                port = nullptr;
            }
        }
    }

    // ---- 启动（挂起），放进 Job 后再恢复，避免它在入 Job 前就派生子进程 ----
    STARTUPINFOA si;
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    bool d0, d1, d2;
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = InheritableStd(STD_INPUT_HANDLE, &d0);
    si.hStdOutput = InheritableStd(STD_OUTPUT_HANDLE, &d1);
    si.hStdError = InheritableStd(STD_ERROR_HANDLE, &d2);

    std::string cmd = "\"" + opt.exe_path + "\"";
    if (!opt.args.empty()) cmd += " " + opt.args;
    std::vector<char> buf(cmd.begin(), cmd.end());
    buf.push_back('\0');

    PROCESS_INFORMATION pi;
    ZeroMemory(&pi, sizeof(pi));
    BOOL ok = CreateProcessA(opt.exe_path.c_str(), buf.data(), nullptr, nullptr, TRUE, CREATE_SUSPENDED, nullptr,
                             opt.work_dir.empty() ? nullptr : opt.work_dir.c_str(), &si, &pi);
    DWORD create_err = ok ? 0 : GetLastError();
    if (d0) CloseHandle(si.hStdInput);
    if (d1) CloseHandle(si.hStdOutput);
    if (d2) CloseHandle(si.hStdError);
    if (!ok) {
        res.outcome = Outcome::LaunchFailed;
        res.launch_error = create_err;
        if (port) CloseHandle(port);
        if (job) CloseHandle(job);
        return res;
    }
    res.main_pid = pi.dwProcessId;
    if (!job || !AssignProcessToJobObject(job, pi.hProcess)) res.job_ok = false;
    ResumeThread(pi.hThread);
    CloseHandle(pi.hThread);

    std::map<DWORD, Tracked> tracked;
    auto track = [&](DWORD pid) {
        if (pid == 0 || pid == pi.dwProcessId || tracked.count(pid)) return;
        HANDLE h = OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION | SYNCHRONIZE, FALSE, pid);
        if (!h) return;   // 已经退出且句柄全关了，取不到退出码
        char path[MAX_PATH * 2];
        DWORD n = sizeof(path);
        std::string name = QueryFullProcessImageNameA(h, 0, path, &n) ? BaseName(std::string(path, n)) : "";
        tracked[pid] = Tracked{name, h};
    };
    auto list_job = [&]() {
        if (!job || !res.job_ok) return;
        // 先查个数，再按个数分配
        JOBOBJECT_BASIC_PROCESS_ID_LIST probe;
        ZeroMemory(&probe, sizeof(probe));
        QueryInformationJobObject(job, JobObjectBasicProcessIdList, &probe, sizeof(probe), nullptr);
        DWORD cap = probe.NumberOfAssignedProcesses + 16;
        std::vector<char> mem(sizeof(JOBOBJECT_BASIC_PROCESS_ID_LIST) + cap * sizeof(ULONG_PTR));
        auto *lst = reinterpret_cast<JOBOBJECT_BASIC_PROCESS_ID_LIST *>(mem.data());
        if (!QueryInformationJobObject(job, JobObjectBasicProcessIdList, lst, (DWORD)mem.size(), nullptr)) return;
        for (DWORD i = 0; i < lst->NumberOfProcessIdsInList; ++i) track((DWORD)lst->ProcessIdList[i]);
    };
    auto is_fatal = [&](const std::string &name) {
        return std::find(opt.fatal_children.begin(), opt.fatal_children.end(), Lower(name)) != opt.fatal_children.end();
    };
    auto kill_tree = [&]() {
        if (job && res.job_ok) TerminateJobObject(job, 1);
        else TerminateProcess(pi.hProcess, 1);
    };

    ULONGLONG last_tick = GetTickCount64();
    bool done = false;
    while (!done) {
        // 1) 等完成端口消息（最多 tick_ms），顺手登记新进程
        if (port) {
            DWORD msg = 0;
            ULONG_PTR key = 0;
            LPOVERLAPPED ov = nullptr;
            DWORD wait = opt.tick_ms;
            while (GetQueuedCompletionStatus(port, &msg, &key, &ov, wait)) {
                wait = 0;
                if (msg == JOB_OBJECT_MSG_NEW_PROCESS) track((DWORD)(ULONG_PTR)ov);
            }
        } else {
            WaitForSingleObject(pi.hProcess, opt.tick_ms);
        }
        // 2) 补漏
        list_job();

        // 3) 子进程：先于主进程检查。主进程常常是因为子进程失败才结束的
        for (auto it = tracked.begin(); it != tracked.end();) {
            if (WaitForSingleObject(it->second.handle, 0) != WAIT_OBJECT_0) {
                ++it;
                continue;
            }
            DWORD code = 0;
            GetExitCodeProcess(it->second.handle, &code);
            ChildExit ce{it->second.name, it->first, code};
            CloseHandle(it->second.handle);
            it = tracked.erase(it);
            if (code == 0) continue;
            if (is_fatal(ce.name)) {
                kill_tree();
                res.outcome = Outcome::ChildFailed;
                res.failed = ce;
                done = true;
                break;
            }
            if (opt.on_other_child_failed) opt.on_other_child_failed(ce);
        }
        if (done) break;

        // 4) 主进程
        if (WaitForSingleObject(pi.hProcess, 0) == WAIT_OBJECT_0) {
            DWORD code = 0;
            GetExitCodeProcess(pi.hProcess, &code);
            res.main_exit_code = code;
            res.outcome = Outcome::Finished;
            kill_tree();   // 主进程走了还残留的子进程一并清掉
            break;
        }

        // 5) 回调
        ULONGLONG now = GetTickCount64();
        if (now - last_tick >= opt.tick_ms) {
            last_tick = now;
            TickInfo ti;
            if (job && res.job_ok) {
                JOBOBJECT_BASIC_ACCOUNTING_INFORMATION acc;
                if (QueryInformationJobObject(job, JobObjectBasicAccountingInformation, &acc, sizeof(acc), nullptr)) {
                    ti.cpu_valid = true;
                    ti.job_cpu_100ns = (unsigned long long)acc.TotalUserTime.QuadPart +
                                       (unsigned long long)acc.TotalKernelTime.QuadPart;
                }
            }
            if (opt.on_tick && opt.on_tick(ti)) {
                kill_tree();
                res.outcome = Outcome::Aborted;
                break;
            }
        }
    }

    for (auto &t : tracked) CloseHandle(t.second.handle);
    CloseHandle(pi.hProcess);
    if (port) CloseHandle(port);
    if (job) CloseHandle(job);
    return res;
}

bool IsAbortExit(unsigned long code) { return code == 3; }   // MinGW 的 abort() 以 3 退出

#else
// ================================================================ Linux：ptrace（2026-11）
//
// 子进程先 PTRACE_TRACEME 再停在 SIGSTOP，托管方设好选项后放行，然后 exec。
// 只跟踪进程的派生 / exec / 退出三类事件，不跟踪系统调用，开销只在这几个时刻。
// 子孙进程由内核自动挂上跟踪（PTRACE_O_TRACEFORK / VFORK；glibc 的 system() 走 posix_spawn，即 vfork），
// 于是能拿到整棵树里每个进程的退出码。PTRACE_O_EXITKILL：托管方退出时，被跟踪的进程一律被 SIGKILL。
// 被跟踪进程的事件要及时 waitpid 处理，否则它停着不动（它的父进程也就等不到它退出），所以主循环每 10 ms 收一次。

namespace {

struct Tracked {
    std::string name;
    unsigned long long ticks = 0;   ///< 最近一次读到的 CPU 时间（时钟滴答）
    bool started = false;           ///< 已经放行过初始的 SIGSTOP
};

std::string ProcExeName(pid_t pid) {
    char path[64], buf[4096];
    snprintf(path, sizeof(path), "/proc/%d/exe", (int)pid);
    ssize_t n = readlink(path, buf, sizeof(buf) - 1);
    if (n <= 0) return "";
    std::string p(buf, (size_t)n);
    const std::string del = " (deleted)";
    if (p.size() > del.size() && p.compare(p.size() - del.size(), del.size(), del) == 0) p.resize(p.size() - del.size());
    return BaseName(p);
}

/// /proc/<pid>/stat 里的 utime + stime（时钟滴答，含该进程全部线程）。取不到返回 false
bool ProcTicks(pid_t pid, unsigned long long *ticks) {
    char path[64];
    snprintf(path, sizeof(path), "/proc/%d/stat", (int)pid);
    FILE *fp = fopen(path, "r");
    if (!fp) return false;
    char buf[2048];
    size_t n = fread(buf, 1, sizeof(buf) - 1, fp);
    fclose(fp);
    buf[n] = '\0';
    const char *p = strrchr(buf, ')');   // 进程名可能含空格和括号，从最后一个 ')' 之后数
    if (!p) return false;
    unsigned long long ut = 0, st = 0;
    // ) 之后：state ppid pgrp session tty_nr tpgid flags minflt cminflt majflt cmajflt utime stime
    if (sscanf(p + 1, " %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %llu %llu", &ut, &st) != 2) return false;
    *ticks = ut + st;
    return true;
}

/// 把 "\"a b\" 1 1" 这样的参数串拆成 {"a b", "1", "1"}：按空白分隔，双引号内的空白不拆，引号本身去掉。
/// 只处理调用方自己拼的参数串（路径加引号 + 数字），不支持转义。
std::vector<std::string> SplitArgs(const std::string &s) {
    std::vector<std::string> out;
    std::string cur;
    bool in_q = false, has = false;
    for (char c : s) {
        if (c == '"') { in_q = !in_q; has = true; continue; }
        if (!in_q && (c == ' ' || c == '\t')) {
            if (has) out.push_back(cur);
            cur.clear();
            has = false;
            continue;
        }
        cur += c;
        has = true;
    }
    if (has) out.push_back(cur);
    return out;
}

unsigned long long TickMs() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (unsigned long long)ts.tv_sec * 1000ULL + (unsigned long long)ts.tv_nsec / 1000000ULL;
}

unsigned long ExitCodeOf(int status) {
    if (WIFEXITED(status)) return (unsigned long)WEXITSTATUS(status);
    if (WIFSIGNALED(status)) return 128UL + (unsigned long)WTERMSIG(status);
    return 0;
}

}  // namespace

bool IsAbortExit(unsigned long code) { return code == 128UL + SIGABRT; }

Result Run(const Options &opt) {
    Result res;
    const bool want_tree = !opt.fatal_children.empty() || (bool)opt.on_other_child_failed;

    std::vector<std::string> args = SplitArgs(opt.args);
    std::vector<char *> argv;
    argv.push_back(const_cast<char *>(opt.exe_path.c_str()));
    for (auto &a : args) argv.push_back(&a[0]);
    argv.push_back(nullptr);

    int errpipe[2];   // exec 失败时子进程把 errno 写回来；成功时 O_CLOEXEC 让它在 exec 时自动关闭
    if (pipe2(errpipe, O_CLOEXEC) != 0) {
        res.outcome = Outcome::LaunchFailed;
        res.launch_error = (unsigned long)errno;
        return res;
    }
    fflush(stdout);
    fflush(stderr);
    const pid_t pid = fork();
    if (pid < 0) {
        res.outcome = Outcome::LaunchFailed;
        res.launch_error = (unsigned long)errno;
        close(errpipe[0]);
        close(errpipe[1]);
        return res;
    }
    if (pid == 0) {   // ---- 子进程：只调异步信号安全的函数
        close(errpipe[0]);
        int e = 0;
        if (!opt.work_dir.empty() && chdir(opt.work_dir.c_str()) != 0) {
            e = errno;
            if (write(errpipe[1], &e, sizeof(e)) < 0) {}
            _exit(127);
        }
        ptrace(PTRACE_TRACEME, 0, nullptr, nullptr);   // 被禁止也照常运行，托管方会发现没跟踪上
        raise(SIGSTOP);                                // 等托管方设好选项
        execv(opt.exe_path.c_str(), argv.data());
        e = errno;
        if (write(errpipe[1], &e, sizeof(e)) < 0) {}
        _exit(127);
    }
    close(errpipe[1]);
    res.main_pid = (unsigned long)pid;

    // ---- 等它停在 SIGSTOP，设跟踪选项
    int st = 0;
    while (waitpid(pid, &st, __WALL | WUNTRACED) < 0 && errno == EINTR) {}
    bool traced = false;
    if (WIFSTOPPED(st)) {
        long o = PTRACE_O_TRACEEXEC | PTRACE_O_TRACEEXIT | PTRACE_O_EXITKILL;
        if (want_tree) o |= PTRACE_O_TRACEFORK | PTRACE_O_TRACEVFORK;
        traced = ptrace(PTRACE_SETOPTIONS, pid, nullptr, (void *)o) == 0;
        if (traced) ptrace(PTRACE_CONT, pid, nullptr, nullptr);
        else kill(pid, SIGCONT);   // 没跟踪上：它是被 raise 停住的普通进程
    }
    res.job_ok = traced;

    // ---- exec 成功与否
    int exec_err = 0;
    ssize_t got;
    while ((got = read(errpipe[0], &exec_err, sizeof(exec_err))) < 0 && errno == EINTR) {}
    close(errpipe[0]);
    if (got == (ssize_t)sizeof(exec_err) || !WIFSTOPPED(st)) {
        res.outcome = Outcome::LaunchFailed;
        res.launch_error = got == (ssize_t)sizeof(exec_err) ? (unsigned long)exec_err : 0;
        kill(pid, SIGKILL);
        while (true) {   // 收尸（被跟踪时中间还会有 exit 事件停一下）
            if (waitpid(pid, &st, __WALL) < 0) { if (errno == EINTR) continue; break; }
            if (WIFEXITED(st) || WIFSIGNALED(st)) break;
            if (traced) ptrace(PTRACE_CONT, pid, nullptr, nullptr);
        }
        return res;
    }

    std::map<pid_t, Tracked> tracked;              // 主进程以外、还活着的被跟踪进程
    Tracked main_t;
    main_t.name = BaseName(opt.exe_path);
    main_t.started = true;
    unsigned long long done_ticks = 0;             // 已退出进程的 CPU 时间
    bool main_alive = true;
    const long clk = sysconf(_SC_CLK_TCK) > 0 ? sysconf(_SC_CLK_TCK) : 100;

    auto is_fatal = [&](const std::string &name) {
        return std::find(opt.fatal_children.begin(), opt.fatal_children.end(), Lower(name)) != opt.fatal_children.end();
    };
    auto kill_tree = [&]() {
        for (auto &t : tracked) kill(t.first, SIGKILL);
        if (main_alive) kill(pid, SIGKILL);
    };
    auto info_of = [&](pid_t w) -> Tracked & {
        if (w == pid) return main_t;
        auto it = tracked.find(w);
        if (it == tracked.end()) it = tracked.emplace(w, Tracked()).first;   // 事件比父进程的 fork 通知先到
        return it->second;
    };

    // 处理一个 waitpid 结果。子进程失败记进 child_fail（多个只记第一个），主进程退出记进 main_status
    bool have_main_exit = false, have_child_fail = false;
    int main_status = 0;
    ChildExit child_fail;
    auto handle = [&](pid_t w, int s) {
        if (WIFSTOPPED(s)) {
            const int sig = WSTOPSIG(s), ev = (s >> 16) & 0xff;
            Tracked &t = info_of(w);
            if (sig == SIGTRAP && ev != 0) {
                unsigned long msg = 0;
                ptrace(PTRACE_GETEVENTMSG, w, nullptr, &msg);
                if (ev == PTRACE_EVENT_FORK || ev == PTRACE_EVENT_VFORK) {
                    Tracked &c = info_of((pid_t)msg);
                    if (c.name.empty()) c.name = t.name;   // exec 之前沿用父进程的名字
                } else if (ev == PTRACE_EVENT_EXEC) {
                    std::string n = ProcExeName(w);
                    if (!n.empty()) t.name = n;
                } else if (ev == PTRACE_EVENT_EXIT) {
                    ProcTicks(w, &t.ticks);   // 退出前最后一刻的 CPU 时间
                }
                ptrace(PTRACE_CONT, w, nullptr, nullptr);
            } else if (sig == SIGSTOP && !t.started) {
                t.started = true;                          // 新挂上的子进程的初始停顿，不转发
                if (t.name.empty()) t.name = ProcExeName(w);
                ptrace(PTRACE_CONT, w, nullptr, nullptr);
            } else if (sig == SIGSTOP || sig == SIGTSTP || sig == SIGTTIN || sig == SIGTTOU) {
                // 暂停类信号不转发：转发后进程进入 group-stop，又报一次同样的停顿，会来回打转。批处理程序不需要暂停
                ptrace(PTRACE_CONT, w, nullptr, nullptr);
            } else {
                ptrace(PTRACE_CONT, w, nullptr, (void *)(long)sig);   // 普通信号原样转交
            }
            return;
        }
        if (!WIFEXITED(s) && !WIFSIGNALED(s)) return;
        if (w == pid) {
            have_main_exit = true;
            main_status = s;
            main_alive = false;
            done_ticks += main_t.ticks;
            return;
        }
        auto it = tracked.find(w);
        std::string name = it != tracked.end() ? it->second.name : std::string();
        if (it != tracked.end()) {
            done_ticks += it->second.ticks;
            tracked.erase(it);
        }
        const unsigned long code = ExitCodeOf(s);
        if (code == 0) return;
        ChildExit ce{name, (unsigned long)w, code};
        if (is_fatal(name)) {
            if (!have_child_fail) { have_child_fail = true; child_fail = ce; }
        } else if (opt.on_other_child_failed) {
            opt.on_other_child_failed(ce);
        }
    };
    auto drain = [&]() {
        while (true) {
            int s = 0;
            const pid_t w = waitpid(-1, &s, __WALL | WNOHANG);
            if (w > 0) { handle(w, s); continue; }
            if (w < 0 && errno == EINTR) continue;
            break;
        }
    };
    // 终止整棵树并收完尸（最多 5 秒），免得被跟踪的进程一直挂在僵尸状态
    auto kill_and_reap = [&]() {
        kill_tree();
        const unsigned long long t0 = TickMs();
        while ((main_alive || !tracked.empty()) && TickMs() - t0 < 5000) {
            drain();
            kill_tree();
            usleep(10000);
        }
    };

    unsigned long long last_tick = TickMs();
    while (true) {
        drain();
        // 子进程先于主进程判断：主进程常常是因为子进程失败才结束的
        if (have_child_fail) {
            res.outcome = Outcome::ChildFailed;
            res.failed = child_fail;
            kill_and_reap();
            break;
        }
        if (have_main_exit) {
            res.outcome = Outcome::Finished;
            res.main_exit_code = ExitCodeOf(main_status);
            kill_and_reap();   // 主进程走了还残留的子进程一并清掉
            break;
        }
        const unsigned long long now = TickMs();
        if (now - last_tick >= opt.tick_ms) {
            last_tick = now;
            TickInfo ti;
            if (traced) {
                unsigned long long ticks = done_ticks;
                if (main_alive && ProcTicks(pid, &main_t.ticks)) ticks += main_t.ticks;
                for (auto &t : tracked)
                    if (ProcTicks(t.first, &t.second.ticks)) ticks += t.second.ticks;
                ti.cpu_valid = true;
                ti.job_cpu_100ns = ticks * (10000000ULL / (unsigned long long)clk);
            }
            if (opt.on_tick && opt.on_tick(ti)) {
                res.outcome = Outcome::Aborted;
                kill_and_reap();
                break;
            }
        }
        usleep(10000);
    }
    if (!have_main_exit) {   // 被我们杀掉的主进程：退出码按实际收到的
        res.main_exit_code = ExitCodeOf(main_status);
    }
    return res;
}
#endif

void CpuActivity::Feed(const TickInfo &ti, unsigned long long now_ms) {
    if (!ti.cpu_valid) return;
    cpu_100ns_ = ti.job_cpu_100ns;
    if (!have_base_) {
        have_base_ = true;
        base_cpu_ = cpu_100ns_;
        base_ms_ = now_ms;
        return;
    }
    if (now_ms - base_ms_ < window_ms_) return;
    const double cpu_ms = (cpu_100ns_ - base_cpu_) / 1e4;
    if (cpu_ms >= busy_frac_ * (double)(now_ms - base_ms_)) last_active_ms_ = now_ms;
    base_cpu_ = cpu_100ns_;
    base_ms_ = now_ms;
}

}  // namespace procguard
