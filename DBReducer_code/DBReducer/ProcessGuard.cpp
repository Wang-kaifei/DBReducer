/**
 * @file ProcessGuard.cpp
 * @brief 见 ProcessGuard.hpp。
 *
 * 发现子进程用两条路：Job 的完成端口消息（JOB_OBJECT_MSG_NEW_PROCESS，及时）+ 每 tick 列一次
 * Job 内进程（补漏，完成端口消息不保证送达）。发现后立刻 OpenProcess 持有句柄，
 * 这样子进程退出后仍能用 GetExitCodeProcess 取到退出码。
 */
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include "ProcessGuard.hpp"

#include <algorithm>
#include <map>

namespace procguard {

namespace {

struct Tracked {
    std::string name;
    HANDLE handle = nullptr;
};

std::string BaseName(const std::string &path) {
    size_t pos = path.find_last_of("\\/");
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

std::string Lower(std::string n) {
    std::transform(n.begin(), n.end(), n.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return n;
}

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
