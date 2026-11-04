/**
 * @file Status.cpp
 * @brief 见 Status.hpp。
 */
#ifdef _WIN32
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0601
#endif
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <ctime>
#endif
#include "../utils/Platform.h"

#include "Status.hpp"

#include <cstdio>
#include <utility>
#include <vector>

// 版本信息：build.ps1 生成 dbr_build_info.h；CMake 用 target_compile_definitions 直接定义 DBR_BUILD
#if !defined(DBR_BUILD) && defined(__has_include)
#if __has_include("dbr_build_info.h")
#include "dbr_build_info.h"
#endif
#endif
#ifndef DBR_BUILD
#define DBR_BUILD "unknown"
#endif

namespace dbr {

namespace {

struct State {
    bool inited = false;
    std::string status_path;      // ANSI，文件 API 用
    std::string stage;            // UTF-8
    int option = 0;
    int heartbeat_sec = 300;
    std::string phase = "init";
    Result result = Result::running;
    int exit_code = -1;           // -1 表示 null
    std::string start_time;
    std::string end_time;
    unsigned long long start_tick = 0;
    std::vector<std::pair<std::string, std::string>> params;   // key -> 已是 JSON 的值
    bool has_progress = false;
    Progress progress;
    bool has_precheck = false;
    double est_mem_gb = 0, commit_avail_gb = 0, per_worker_gb = 0;
    bool warned = false;
    int workers = 0;
    long long spectra_batch = 0;
    bool has_builddict = false, builddict_warned = false;
    double builddict_est_gb = 0;
    std::vector<std::pair<std::string, std::string>> warnings;   // what -> kv（ASCII）
    ProcInfo failed;
    std::string detail;
    std::string suggestion_code;
    long long proteins = -1;
};
State g;

/// suggestion_code -> 中文文案（与 docs/status-interface.md 2.6 节同步）
const char *SuggestionText(const std::string &code) {
    if (code == "reduce_multip") return "预估内存超出可提交余额，建议降低 multip 或将库分片";
    if (code == "reduce_maxprolen") return "单个进程的蛋白分批过大，建议降低 maxprolen";
    if (code == "reduce_maxspec") return "单个进程的谱图批次过大，建议降低 maxspec";
    if (code == "check_input") return "请检查 cfg 中的路径与输入文件是否存在";
    if (code == "check_pfind") return "请检查 pFind 安装是否完整（thirdparty" PSEP "pFind3" PSEP "bin）";
    if (code == "rerun_resume") return "可设置 resume=1 后重跑，已完成的文件会被跳过";
    if (code == "see_pfind_log") return "请查看 SearchRes 下 pFind 的日志定位原因";
    if (code == "free_memory") return "DBReducer 读取原始库时内存不足，请关闭其他占用内存的程序或换用内存更大的机器";
    if (code == "report_bug") return "程序内部错误，请保留日志与 dbr_status.json 联系开发者";
    if (code == "shard_db") return "翻译库过大，缩库最后读取原始库时可能内存不足，建议将库分片";
    if (code == "check_pfind_license")
        return "pFind 没有产出任何结果就结束了，很可能是 pFind 已过期（按系统日期判断），请更新 pFind 或检查系统日期";
    if (code == "abort_reduce_multip")
        return "进程以退出码 3 中止（abort），可能是内存不足，也可能是程序内部错误：可先降低 multip 重跑；若仍失败，请保留日志与状态文件联系开发者";
    return "";
}

std::string JsonStr(const std::string &utf8) {
    std::string o = "\"";
    for (unsigned char c : utf8) {
        switch (c) {
            case '"': o += "\\\""; break;
            case '\\': o += "\\\\"; break;
            case '\n': o += "\\n"; break;
            case '\r': o += "\\r"; break;
            case '\t': o += "\\t"; break;
            default:
                if (c < 0x20) {
                    char b[8];
                    snprintf(b, sizeof(b), "\\u%04x", c);
                    o += b;
                } else {
                    o += (char)c;
                }
        }
    }
    return o + "\"";
}

std::string Num(double v) {
    char b[64];
    if (v == (double)(long long)v && v > -1e15 && v < 1e15) snprintf(b, sizeof(b), "%lld", (long long)v);
    else snprintf(b, sizeof(b), "%.4g", v);
    return b;
}

std::string BaseName(const std::string &path) {
    std::string p = path;
    while (!p.empty() && (p.back() == '\\' || p.back() == '/')) p.pop_back();
    size_t pos = p.find_last_of("\\/");
    return pos == std::string::npos ? p : p.substr(pos + 1);
}

std::string BuildJson() {
    std::string j = "{\n";
    auto field = [&](const char *k, const std::string &v, bool last = false) {
        j += "  \"";
        j += k;
        j += "\": ";
        j += v;
        j += last ? "\n" : ",\n";
    };
    field("schema_version", "1");
    field("tool", JsonStr("DBReducer"));
    field("tool_build", JsonStr(DBR_BUILD));
    field("pid", Num(plat::Pid()));
    field("result", JsonStr(ResultName(g.result)));
    field("exit_code", g.exit_code < 0 ? "null" : Num(g.exit_code));
    field("stage", JsonStr(g.stage));
    field("option", Num(g.option));
    field("phase", JsonStr(g.phase));
    field("start_time", JsonStr(g.start_time));
    field("updated_at", JsonStr(NowIso()));
    if (!g.end_time.empty()) field("end_time", JsonStr(g.end_time));
    field("elapsed_sec", Num((double)ElapsedSec()));
    field("heartbeat_sec", Num(g.heartbeat_sec));
    if (g.has_progress) {
        std::string p = "{ \"qry_done\": " + Num(g.progress.qry_done) +
                        ", \"qry_total\": " + Num(g.progress.qry_total) +
                        ", \"newest_file\": " + JsonStr(g.progress.newest_file) +
                        ", \"newest_mtime\": " + JsonStr(g.progress.newest_mtime) +
                        ", \"newest_age_sec\": " + Num((double)g.progress.newest_age_sec) + " }";
        field("progress", p);
    }
    if (g.failed.valid) {
        field("failed_process", "{ \"name\": " + JsonStr(g.failed.name) + ", \"pid\": " + Num(g.failed.pid) +
                                    ", \"exit_code\": " + JsonStr(Hex8(g.failed.exit_code)) + " }");
    }
    std::string pr = "{";
    for (size_t i = 0; i < g.params.size(); ++i) {
        pr += (i ? ", " : " ") + JsonStr(g.params[i].first) + ": " + g.params[i].second;
    }
    pr += g.params.empty() ? "}" : " }";
    field("params", pr);
    if (g.has_precheck || g.has_builddict) {
        // 搜库预检只在 option 0 时做；BuildDict 估算只在 option 3 时做（D9 之后其余 option 流式读库）
        std::string pc = "{";
        char b[256];
        if (g.has_precheck) {
            snprintf(b, sizeof(b),
                     " \"est_mem_gb\": %.1f, \"warned\": %s, \"workers\": %d, \"spectra_batch\": %lld, "
                     "\"per_worker_gb\": %.2f,",
                     g.est_mem_gb, g.warned ? "true" : "false", g.workers, g.spectra_batch, g.per_worker_gb);
            pc += b;
        }
        if (g.has_builddict) {
            snprintf(b, sizeof(b), " \"builddict_est_gb\": %.1f, \"builddict_warned\": %s,", g.builddict_est_gb,
                     g.builddict_warned ? "true" : "false");
            pc += b;
        }
        snprintf(b, sizeof(b), " \"commit_avail_gb\": %.1f }", g.commit_avail_gb);
        pc += b;
        field("precheck", pc);
    }
    if (!g.warnings.empty()) {
        std::string w = "[";
        for (size_t i = 0; i < g.warnings.size(); ++i) {
            w += (i ? ", " : " ") + std::string("{ \"what\": ") + JsonStr(g.warnings[i].first) +
                 ", \"info\": " + JsonStr(g.warnings[i].second) + " }";
        }
        w += " ]";
        field("warnings", w);
    }
    if (!g.detail.empty()) field("detail", JsonStr(g.detail));
    if (!g.suggestion_code.empty()) {
        field("suggestion_code", JsonStr(g.suggestion_code));
        field("suggestion", JsonStr(SuggestionText(g.suggestion_code)));
    }
    if (g.proteins >= 0) field("proteins", Num((double)g.proteins));
    // 去掉最后一个逗号
    if (j.size() >= 2 && j[j.size() - 2] == ',') j.erase(j.size() - 2, 1);
    j += "}\n";
    return j;
}

/// 原子写：先写 .tmp，再 MoveFileEx 覆盖。写不成功只放弃本次，绝不影响运行。
void WriteStatus() {
    if (!g.inited || g.status_path.empty()) return;
    const std::string json = BuildJson();
    const std::string tmp = g.status_path + ".tmp";
    FILE *fp = fopen(tmp.c_str(), "wb");
    if (!fp) return;   // 目录还不存在等情况
    fwrite(json.data(), 1, json.size(), fp);
    fclose(fp);
    for (int i = 0; i < 3; ++i) {
        if (plat::ReplaceFile(tmp, g.status_path))
            return;
        plat::SleepMs(200);   // 读取方正好打开着文件（Windows）
    }
    plat::RemoveFile(tmp);
}

}  // namespace

const char *ResultName(Result r) {
    switch (r) {
        case Result::running: return "running";
        case Result::ok: return "ok";
        case Result::param_error: return "param_error";
        case Result::precheck_failed: return "precheck_failed";
        case Result::pfind_failed: return "pfind_failed";
        case Result::worker_crashed: return "worker_crashed";
        case Result::stalled: return "stalled";
        case Result::output_missing: return "output_missing";
        case Result::internal_error: return "internal_error";
        case Result::cancelled: return "cancelled";
    }
    return "internal_error";
}

int ResultExitCode(Result r) {
    switch (r) {
        case Result::ok: return 0;
        case Result::param_error: return 1;
        case Result::precheck_failed: return 10;
        case Result::pfind_failed: return 11;
        case Result::worker_crashed: return 12;
        case Result::stalled: return 13;
        case Result::output_missing: return 14;
        case Result::internal_error: return 15;
        case Result::cancelled: return 16;
        case Result::running: return 15;
    }
    return 15;
}

std::string AnsiToUtf8(const std::string &ansi) {
#ifndef _WIN32
    return ansi;   // Linux 上文件名、cfg 都按 UTF-8，不用转
#else
    if (ansi.empty()) return ansi;
    int wn = MultiByteToWideChar(CP_ACP, 0, ansi.data(), (int)ansi.size(), nullptr, 0);
    if (wn <= 0) return ansi;
    std::wstring w(wn, L'\0');
    MultiByteToWideChar(CP_ACP, 0, ansi.data(), (int)ansi.size(), &w[0], wn);
    int un = WideCharToMultiByte(CP_UTF8, 0, w.data(), wn, nullptr, 0, nullptr, nullptr);
    if (un <= 0) return ansi;
    std::string u(un, '\0');
    WideCharToMultiByte(CP_UTF8, 0, w.data(), wn, &u[0], un, nullptr, nullptr);
    return u;
#endif
}

namespace {
#ifdef _WIN32
std::string FormatLocal(const SYSTEMTIME &st) {
    TIME_ZONE_INFORMATION tz;
    DWORD r = GetTimeZoneInformation(&tz);
    long bias = tz.Bias;
    if (r == TIME_ZONE_ID_DAYLIGHT) bias += tz.DaylightBias;
    else if (r == TIME_ZONE_ID_STANDARD) bias += tz.StandardBias;
    long off = -bias;   // 分钟，东八区为 +480
    char sign = off >= 0 ? '+' : '-';
    if (off < 0) off = -off;
    char b[40];
    snprintf(b, sizeof(b), "%04u-%02u-%02uT%02u:%02u:%02u%c%02ld:%02ld", st.wYear, st.wMonth, st.wDay, st.wHour,
             st.wMinute, st.wSecond, sign, off / 60, off % 60);
    return b;
}
#else
std::string FormatLocal(time_t t) {
    struct tm lt;
    if (!localtime_r(&t, &lt)) return "";
    long off = lt.tm_gmtoff / 60;   // 分钟，东八区为 +480
    char sign = off >= 0 ? '+' : '-';
    if (off < 0) off = -off;
    char b[40];
    snprintf(b, sizeof(b), "%04d-%02d-%02dT%02d:%02d:%02d%c%02ld:%02ld", lt.tm_year + 1900, lt.tm_mon + 1, lt.tm_mday,
             lt.tm_hour, lt.tm_min, lt.tm_sec, sign, off / 60, off % 60);
    return b;
}
#endif
}  // namespace

std::string NowIso() {
#ifndef _WIN32
    return FormatLocal(time(nullptr));
#else
    SYSTEMTIME st;
    GetLocalTime(&st);
    return FormatLocal(st);
#endif
}

std::string FileTimeIso(unsigned long long filetime) {
    if (filetime == 0) return "";
#ifndef _WIN32
    if (filetime < plat::kUnixEpochAsFileTime) return "";
    return FormatLocal((time_t)((filetime - plat::kUnixEpochAsFileTime) / 10000000ULL));
#else
    FILETIME ft, lft;
    ft.dwLowDateTime = (DWORD)(filetime & 0xFFFFFFFFULL);
    ft.dwHighDateTime = (DWORD)(filetime >> 32);
    SYSTEMTIME st;
    if (!FileTimeToLocalFileTime(&ft, &lft) || !FileTimeToSystemTime(&lft, &st)) return "";
    return FormatLocal(st);
#endif
}

std::string Hex8(unsigned long v) {
    char b[16];
    snprintf(b, sizeof(b), "0x%08lX", v);
    return b;
}

void Log(const char *event, const std::string &kv) {
    fflush(stdout);   // 先把其他输出冲出去，保证顺序
    // 前导换行：pFind 刷新进度只输出 \r 不换行，不加的话本行会接在进度后面，按行匹配就漏了。
    // 与本行放在同一次 printf / flush 里，避免子进程的输出插进两者之间
    printf("\n[DBR:%s] %s\n", event, kv.c_str());
    fflush(stdout);
}

void StatusInit(const std::string &search_res_ansi, int option, int heartbeat_sec) {
    g.inited = true;
    std::string dir = search_res_ansi;
    while (!dir.empty() && (dir.back() == '\\' || dir.back() == '/')) dir.pop_back();
    g.status_path = dir.empty() ? "" : dir + PSEP "dbr_status.json";
    g.stage = AnsiToUtf8(BaseName(dir));
    g.option = option;
    g.heartbeat_sec = heartbeat_sec > 0 ? heartbeat_sec : 300;
    g.start_time = NowIso();
    g.start_tick = plat::TickMs();
}

void StatusParamRaw(const std::string &key, const std::string &value_json) {
    for (auto &p : g.params) {
        if (p.first == key) {
            p.second = value_json;
            return;
        }
    }
    g.params.emplace_back(key, value_json);
}
void StatusParamNum(const std::string &key, double v) { StatusParamRaw(key, Num(v)); }
void StatusParamStr(const std::string &key, const std::string &value_utf8) { StatusParamRaw(key, JsonStr(value_utf8)); }

void StatusPrecheck(double est_mem_gb, double commit_avail_gb, bool warned, int workers, long long spectra_batch,
                    double per_worker_gb) {
    g.has_precheck = true;
    g.est_mem_gb = est_mem_gb;
    g.commit_avail_gb = commit_avail_gb;
    g.warned = warned;
    g.workers = workers;
    g.spectra_batch = spectra_batch;
    g.per_worker_gb = per_worker_gb;
}

void StatusBuildDict(double est_gb, double commit_avail_gb, bool warned) {
    g.has_builddict = true;
    g.builddict_est_gb = est_gb;
    g.builddict_warned = warned;
    if (!g.has_precheck) g.commit_avail_gb = commit_avail_gb;
}

void Warn(const std::string &what, const std::string &kv) {
    Log("WARN", "what=" + what + (kv.empty() ? "" : " " + kv));
    // 同一类告警反复出现（如多个子进程非 0 退出）时只留最近 50 条，防止状态文件无限变大
    if (g.warnings.size() >= 50) g.warnings.erase(g.warnings.begin());
    g.warnings.emplace_back(what, kv);
    WriteStatus();
}

void StatusPhase(const char *phase) {
    g.phase = phase;
    Log("PHASE", std::string("phase=") + phase + " elapsed_sec=" + std::to_string(ElapsedSec()));
    WriteStatus();
}

const char *CurrentPhase() { return g.phase.c_str(); }

void StatusHeartbeat(const Progress &p) {
    g.has_progress = true;
    g.progress = p;
    WriteStatus();
}

void StatusSetProteins(long long n) { g.proteins = n; }
void StatusSuggest(const std::string &code) { g.suggestion_code = code; }

long long ElapsedSec() {
    if (g.start_tick == 0) return 0;
    return (long long)((plat::TickMs() - g.start_tick) / 1000ULL);
}

int StatusFinish(Result r, const std::string &detail_utf8, const std::string &suggestion_code, const ProcInfo *proc) {
    g.result = r;
    g.exit_code = ResultExitCode(r);
    if (!detail_utf8.empty()) g.detail = detail_utf8;
    if (!suggestion_code.empty()) g.suggestion_code = suggestion_code;
    if (proc && proc->valid) g.failed = *proc;
    g.end_time = NowIso();
    const std::string failed_phase = g.phase;
    if (r == Result::ok) g.phase = "done";   // 失败时保留出错的那一步，GUI 据此显示"在哪一步终止"
    WriteStatus();
    if (r != Result::ok) {
        std::string kv = std::string("result=") + ResultName(r) + " exit_code=" + std::to_string(g.exit_code) +
                         " phase=" + failed_phase;
        if (g.failed.valid) {
            kv += " process=" + g.failed.name + " pid=" + std::to_string(g.failed.pid) +
                  " proc_exit=" + Hex8(g.failed.exit_code);
        }
        if (!g.suggestion_code.empty()) kv += " suggestion_code=" + g.suggestion_code;
        Log("ERROR", kv);
    }
    std::string kv = std::string("result=") + ResultName(r) + " exit_code=" + std::to_string(g.exit_code) +
                     " elapsed_sec=" + std::to_string(ElapsedSec());
    if (r == Result::ok && g.proteins >= 0) kv += " proteins=" + std::to_string(g.proteins);
    Log("DONE", kv);
    return g.exit_code;
}

}  // namespace dbr
