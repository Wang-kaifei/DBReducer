/*
 * @Descripttion:
 * @version:
 * @Author: sueRimn
 * @Date: 2021-12-20 22:30:44
 * @LastEditors: wangkaifei kfwang@stu.xidian.edu.cn
 * @LastEditTime: 2025-01-24 14:27:04
 *
 * 2026-10：加入运行状态报告（退出码 / dbr_status.json / [DBR:事件] 日志行，约定见 docs/status-interface.md），
 *          用进程托管替换 system("Searcher.exe")：pFind 子进程崩溃即终止整棵进程树并报错，
 *          心跳、停滞超时、文件级续跑（resume）、启动前资源预检、产物校验。
 *          各 option 的业务逻辑（读修饰、写 pFind.cfg、调哪个缩库函数）与原来相同。
 */
#include <iostream>
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#endif
#include <string>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <stdio.h>
#include <deque>
#include <chrono>
#include <algorithm>
#include "../utils/IOUtil.h"
#include "../utils/StringUtil.h"
#include "../utils/Platform.h"
#include "Param.hpp"
#include "Rebuild.hpp"
#include "Rebuild-spec.hpp"
#include "Status.hpp"
#include "ProcessGuard.hpp"
#include "Precheck.hpp"

using namespace std;
void ReadMod(const Params &pAnno_param, std::unordered_set<std::string> *mods);
void WriteMod(const std::string &raw_mod_path, const std::unordered_set<std::string> &mods, const std::string &new_mod_path);

namespace {

using dbr::DbrError;
using dbr::Result;
using dbr::AnsiToUtf8;

// 文件大小 / 当前时间（FILETIME 单位）：两个平台的实现在 utils/Platform.h
using plat::FileSize;
using plat::NowFileTime;

#ifdef _WIN32
#define SYS_ERR "Windows 错误码"
#else
#define SYS_ERR "系统错误码"
#endif

std::string Lower(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
	return s;
}

bool EndsWithI(const std::string &s, const std::string &suffix) {
	return s.size() >= suffix.size() && Lower(s.substr(s.size() - suffix.size())) == Lower(suffix);
}

std::string StripSlash(std::string s) {
	while (!s.empty() && (s.back() == '\\' || s.back() == '/')) s.pop_back();
	return s;
}

/// 日志行只允许 ASCII：非 ASCII 字节替换成 '?'
std::string Ascii(const std::string &s) {
	std::string o = s;
	for (auto &c : o) if ((unsigned char)c >= 0x80 || c == ' ') c = (c == ' ') ? '_' : '?';
	return o;
}

std::string Fmt(const char *fmt, double v) {
	char b[64];
	snprintf(b, sizeof(b), fmt, v);
	return b;
}

/// 本次运行里 pFind 搜库的耗时（秒）；-1 表示本次没有搜库（option 1/2/3）
long long g_search_sec = -1;

void MakeDir(const std::string &dir) {
	if (0 != access(dir.c_str(), 0)) {
#ifdef _WIN32
		if (mkdir(dir.c_str()) == -1) {
#else
		if (mkdir(dir.c_str(), 0777) == -1) {
#endif
			cout << "Path: " << dir << endl;
			printf("Build Folder failed!\n");
			throw DbrError(Result::param_error, "无法创建输出目录：" + AnsiToUtf8(dir), "check_input");
		}
	}
}

// ---------------------------------------------------------------- 输出目录扫描（心跳 / 停滞 / 续跑）

struct DirScan {
	std::string newest;              // 最近修改的文件名（ANSI）
	unsigned long long newest_ft = 0;
	int qry = 0;                     // *.qry.res 个数（只数修改时间 >= qry_since 的）
};

void ScanDir(const std::string &dir, int depth, unsigned long long qry_since, DirScan *s) {
	for (const auto &e : plat::ListDir(dir)) {
		const std::string &n = e.name;
		if (e.is_dir) {
			if (depth > 0) ScanDir(dir + PSEP + n, depth - 1, qry_since, s);
			continue;
		}
		if (Lower(n).rfind("dbr_status.json", 0) == 0) continue;   // 自己写的状态文件不算 pFind 的产出
		const unsigned long long ft = e.mtime;
		if (ft > s->newest_ft) {
			s->newest_ft = ft;
			s->newest = n;
		}
		if (EndsWithI(n, ".qry.res") && ft >= qry_since) s->qry++;
	}
}

/// 谱图文件列表：msmsfolder>0 时 msmspath 是目录，取其中扩展名与 msmstype 相同的文件；否则 msmspath 就是文件
std::vector<std::string> ListMsFiles(const Params &p) {
	std::vector<std::string> files;
	if (p.msmsfolder <= 0) {
		for (const auto &raw : p.msmspaths) files.push_back(StripSlash(raw));
		return files;
	}
	const std::string ext = "." + Lower(p.msmstype);
	for (const auto &raw : p.msmspaths) {
		std::string d = StripSlash(raw);
		for (const auto &e : plat::ListDir(d))
			if (!e.is_dir && EndsWithI(e.name, ext)) files.push_back(d + PSEP + e.name);
	}
	return files;
}

/// 谱图文件总数（progress.qry_total）。数不出来返回 -1。
int CountMsFiles(const Params &p) {
	const size_t n = ListMsFiles(p).size();
	return n > 0 ? (int)n : -1;
}

/// 文件是否完整：非空且最后一个字节是换行（pFind 的 qry.res 每行以 \n 结束、没有文件尾标记）
bool EndsWithNewline(const std::string &path) {
	FILE *fp = fopen(path.c_str(), "rb");
	if (!fp) return false;
	bool ok = false;
	if (plat::Seek64(fp, -1, SEEK_END) == 0) ok = (fgetc(fp) == '\n');
	fclose(fp);
	return ok;
}

/// 启动 pFind 前清理输出目录：
///  - 总是删掉残留的 <id>.complete（multip>1 时 Searcher 靠它判断各 worker 是否结束；上次中途被杀留下的
///    会让这次的 Searcher 不等对应 worker 就进入蛋白推断）
///  - resume=1 时删掉写了一半的 *.qry.res（空文件或末字节不是换行），其余的交给 pFind 跳过
void PrepareSearchDir(const Params &p) {
	const std::string dir = StripSlash(p.SearchRes);
	int skip = 0, removed = 0, removed_complete = 0;
	std::string removed_names;
	for (const auto &e : plat::ListDir(dir)) {
		if (e.is_dir) continue;
		const std::string &n = e.name;
		const std::string full = dir + PSEP + n;
		if (EndsWithI(n, ".complete")) {
			if (plat::RemoveFile(full)) ++removed_complete;
		} else if (p.resume == 1 && EndsWithI(n, ".qry.res")) {
			if (!EndsWithNewline(full)) {
				if (plat::RemoveFile(full)) {
					++removed;
					removed_names += (removed_names.empty() ? "" : ",") + Ascii(n);
				}
			} else {
				++skip;
			}
		}
	}
	if (removed_complete > 0) dbr::Log("CLEANUP", "removed_complete=" + std::to_string(removed_complete));
	if (p.resume == 1) {
		std::string kv = "skip=" + std::to_string(skip) + " removed_incomplete=" + std::to_string(removed);
		if (!removed_names.empty()) kv += " file=" + removed_names;
		dbr::Log("RESUME", kv);
	}
}

// ---------------------------------------------------------------- 预检 / 搜库 / 校验

void Precheck(const Params &param) {
	dbr::StatusPhase("precheck");
	dbr::PrecheckResult pc = dbr::RunPrecheck(param.ProteinDatabase, param.maxprolen, param.maxspec, param.multip,
	                                          ListMsFiles(param), param.msmstype);
	if (!pc.ok) {
		dbr::Warn("precheck_skipped");
		return;
	}
	dbr::StatusPrecheck(pc.est_mem_gb, pc.commit_avail_gb, pc.warned, pc.workers, pc.spectra_batch, pc.per_worker_gb);
	printf("Precheck: entries~%lld aa~%lld, spectra/file~%lld -> batch %lld, per_worker~%.2f GB x workers %d "
	       "(multip %d) = %.1f GB, commit available %.1f GB\n",
	       pc.entries_est, pc.aa_est, pc.spectra_per_file_est, pc.spectra_batch, pc.per_worker_gb, pc.workers,
	       param.multip, pc.est_mem_gb, pc.commit_avail_gb);
	if (pc.warned) {
		std::string kv = "est_mem_gb=" + Fmt("%.1f", pc.est_mem_gb) +
		                 " commit_avail_gb=" + Fmt("%.1f", pc.commit_avail_gb) +
		                 " per_worker_gb=" + Fmt("%.2f", pc.per_worker_gb) +
		                 " multip=" + std::to_string(param.multip) +
		                 " workers=" + std::to_string(pc.workers) +
		                 " spectra_batch=" + std::to_string(pc.spectra_batch);
		if (pc.suggest_multip >= 1) kv += " suggest_multip=" + std::to_string(pc.suggest_multip);
		kv += " suggestion_code=" + pc.suggestion_code;
		dbr::Warn("memory", kv);
		dbr::StatusSuggest(pc.suggestion_code);
	}
}

void RunSearcher(const Params &param, const std::string &pFind_path, const std::string &cfg) {
	dbr::StatusPhase("search");
	procguard::Options o;
	o.exe_path = pFind_path + PSEP "Searcher" EXE;
	if (FileSize(o.exe_path) < 0)
		throw DbrError(Result::param_error, "找不到 pFind 的 Searcher" EXE "：" + AnsiToUtf8(o.exe_path), "check_pfind");
	PrepareSearchDir(param);
	// resume=1 时末尾加 " 1"：Searcher 会原样传给 SearchEngine，后者跳过已有 qry.res 的谱图文件
	o.args = "\"" + cfg + "\"" + (param.resume == 1 ? " 1" : "");
	o.work_dir = pFind_path;
	o.fatal_children = {"tagsearcher" EXE, "searchengine" EXE, "pssearcher" EXE, "reranker" EXE, "group" EXE};
	o.tick_ms = 1000;
	printf("Searcher" EXE " %s\n", o.args.c_str());
	fflush(stdout);

	const std::string dir = StripSlash(param.SearchRes);
	const unsigned long long search_start_ft = NowFileTime();
	// resume=1 时上次留下的 qry.res 也算已完成；否则只数本次新写出的
	const unsigned long long qry_since = param.resume == 1 ? 0 : search_start_ft;
	const int qry_total = CountMsFiles(param);
	const unsigned long long hb_ms = (unsigned long long)param.heartbeat_sec * 1000ULL;
	const double stall_sec = param.stall_timeout_min * 60.0;
	unsigned long long scan_ms = hb_ms;
	if (stall_sec > 0) scan_ms = std::min<unsigned long long>(scan_ms, std::max<unsigned long long>(1000ULL, (unsigned long long)(stall_sec * 1000 / 4)));
	scan_ms = std::min<unsigned long long>(scan_ms, 15000ULL);

	unsigned long long newest_ft = search_start_ft;
	std::string newest_name;
	unsigned long long last_hb = plat::TickMs(), last_scan = 0;
	long long stalled_age = -1, stalled_cpu_idle = -1;
	// 停滞 = 输出目录无文件更新 且 整棵进程树 CPU 几乎不增长，两者都持续 stall_sec。
	// 只看文件会误判：ReRanker 等阶段可以满负荷计算很久不写文件
	// 结算窗口取阈值的 1/4（1–10 秒）：满负荷的进程被计为空闲的时长不会超过一个窗口
	const unsigned cpu_window_ms =
	    (unsigned)std::min(10000.0, std::max(1000.0, stall_sec > 0 ? stall_sec * 1000 / 4 : 10000.0));
	procguard::CpuActivity cpu(plat::TickMs(), cpu_window_ms);

	auto make_progress = [&](int qry_done, long long age) {
		dbr::Progress pr;
		pr.qry_done = qry_done;
		pr.qry_total = qry_total;
		pr.newest_file = AnsiToUtf8(newest_name);
		pr.newest_mtime = dbr::FileTimeIso(newest_name.empty() ? 0 : newest_ft);
		pr.newest_age_sec = age;
		return pr;
	};
	{
		DirScan s;
		ScanDir(dir, 2, qry_since, &s);
		dbr::StatusHeartbeat(make_progress(s.qry, 0));   // 先写一次，让读取方立刻看到 qry_total
	}

	o.on_tick = [&](const procguard::TickInfo &ti) -> bool {
		const unsigned long long now = plat::TickMs();
		cpu.Feed(ti, now);
		const bool hb_due = now - last_hb >= hb_ms;
		if (!hb_due && now - last_scan < scan_ms) return false;
		last_scan = now;
		DirScan s;
		ScanDir(dir, 2, qry_since, &s);
		if (s.newest_ft > newest_ft) {
			newest_ft = s.newest_ft;
			newest_name = s.newest;
		}
		const long long age = (long long)((NowFileTime() - newest_ft) / 10000000ULL);
		const long long cpu_idle = cpu.valid() ? (long long)(cpu.idle_ms(now) / 1000ULL) : -1;
		if (hb_due) {
			last_hb = now;
			dbr::Progress pr = make_progress(s.qry, age);
			std::string kv = "elapsed_sec=" + std::to_string(dbr::ElapsedSec()) + " qry_done=" + std::to_string(s.qry) +
			                 " qry_total=" + std::to_string(qry_total) + " newest_age_sec=" + std::to_string(age);
			if (!newest_name.empty()) kv += " newest=" + Ascii(newest_name);
			if (cpu.valid()) kv += " job_cpu_sec=" + Fmt("%.0f", cpu.cpu_sec()) + " cpu_idle_sec=" + std::to_string(cpu_idle);
			dbr::Log("HEARTBEAT", kv);
			dbr::StatusHeartbeat(pr);
		}
		// 取不到 CPU 时间（没能放进 Job）时退回只看文件
		if (stall_sec > 0 && age >= (long long)stall_sec && (cpu_idle < 0 || cpu_idle >= (long long)stall_sec)) {
			stalled_age = age;
			stalled_cpu_idle = cpu_idle;
			dbr::StatusHeartbeat(make_progress(s.qry, age));
			return true;
		}
		return false;
	};
	o.on_other_child_failed = [&](const procguard::ChildExit &c) {
		dbr::Warn("child_exit", "process=" + Ascii(c.name) + " pid=" + std::to_string(c.pid) +
		                            " proc_exit=" + dbr::Hex8(c.exit_code));
	};

	const unsigned long long search_start_tick = plat::TickMs();
	procguard::Result r = procguard::Run(o);
	g_search_sec = (long long)((plat::TickMs() - search_start_tick) / 1000ULL);
	if (!r.job_ok) dbr::Warn("job_unavailable");
	{
		// 结束时再扫一次：心跳之间写出的 qry.res 也要算进最终进度
		DirScan s;
		ScanDir(dir, 2, qry_since, &s);
		if (s.newest_ft > newest_ft) {
			newest_ft = s.newest_ft;
			newest_name = s.newest;
		}
		dbr::StatusHeartbeat(make_progress(s.qry, (long long)((NowFileTime() - newest_ft) / 10000000ULL)));
	}

	switch (r.outcome) {
		case procguard::Outcome::LaunchFailed:
			throw DbrError(Result::pfind_failed,
			               "无法启动 Searcher" EXE "（" SYS_ERR " " + std::to_string(r.launch_error) + "）", "check_pfind");
		case procguard::Outcome::Aborted: {
			char b[400];
			if (stalled_cpu_idle >= 0)
				snprintf(b, sizeof(b),
				         "输出目录已连续 %lld 秒没有文件更新，且 pFind 进程树的 CPU 时间已连续 %lld 秒几乎不增长"
				         "（阈值 stall_timeout_min=%g，即 %.0f 秒），已终止 pFind 进程树",
				         stalled_age, stalled_cpu_idle, param.stall_timeout_min, stall_sec);
			else
				snprintf(b, sizeof(b), "输出目录已连续 %lld 秒没有文件更新（阈值 stall_timeout_min=%g，即 %.0f 秒；"
				         "未能取得进程树 CPU 时间，只按文件判断），已终止 pFind 进程树",
				         stalled_age, param.stall_timeout_min, stall_sec);
			throw DbrError(Result::stalled, b, "see_pfind_log");
		}
		case procguard::Outcome::ChildFailed: {
			dbr::ProcInfo pi;
			pi.valid = true;
			pi.name = r.failed.name;
			pi.pid = r.failed.pid;
			pi.exit_code = r.failed.exit_code;
			std::string detail = "pFind 子进程 " + r.failed.name + " 异常退出（退出码 " + dbr::Hex8(r.failed.exit_code) + "）";
			if (procguard::IsAbortExit(r.failed.exit_code)) detail += "，可能是内存不足，也可能是程序内部错误";
			detail += "，已终止整个 pFind 进程树";
			throw DbrError(Result::worker_crashed, detail, procguard::IsAbortExit(r.failed.exit_code) ? "abort_reduce_multip" : "see_pfind_log", pi);
		}
		case procguard::Outcome::Finished:
			if (r.main_exit_code != 0) {
				dbr::ProcInfo pi;
				pi.valid = true;
				pi.name = "Searcher" EXE;
				pi.pid = r.main_pid;
				pi.exit_code = r.main_exit_code;
				throw DbrError(Result::pfind_failed, "pFind 主进程 Searcher" EXE " 返回非 0（" + dbr::Hex8(r.main_exit_code) + "）",
				               "see_pfind_log", pi);
			}
			break;
	}
}

bool AnyQryRes(const std::string &dir) {
	for (const auto &e : plat::ListDir(StripSlash(dir)))
		if (EndsWithI(e.name, ".qry.res")) return true;
	return false;
}

/// pFind 结果文件必须存在且非 0 字节，否则是故障（14）。
/// 返回值：表头之后是否至少有一行数据。只有表头 = pFind 正常跑完但 0 命中，不算故障，由调用方按 no_hits 处理。
/// （pFind.protein 有 2 行表头，pFind.spectra 有 1 行）
bool RequireDataLines(const std::string &path, const std::string &name, int header_lines, const Params &param) {
	const long long sz = FileSize(path);
	if (sz <= 0) {
		std::string detail = (sz < 0 ? "缺少 pFind 结果文件 " : "pFind 结果文件为空：") + name + "（" + AnsiToUtf8(path) + "）";
		std::string code = "see_pfind_log";
		// 本次刚搜过库、很快就结束、连一个 qry.res 都没写：多半是 pFind 过期（内核按系统日期判断，过期时什么都不写就以 0 退出）
		if (g_search_sec >= 0 && g_search_sec < 300 && !AnyQryRes(param.SearchRes)) {
			detail += "。pFind 搜库仅用时 " + std::to_string(g_search_sec) + " 秒且没有任何中间结果，很可能 pFind 已过期";
			code = "check_pfind_license";
		}
		throw DbrError(Result::output_missing, detail, code);
	}
	std::ifstream in(path);
	std::string line;
	int n = 0;
	while (std::getline(in, line)) {
		if (n++ < header_lines) continue;
		if (!Strip(&line).empty()) return true;
	}
	return false;
}

/// 0 命中：写出空的 DBReducer.fasta，proteins=0，告警 no_hits。不算失败。
void NoHits(const Params &param, const char *reason) {
	FILE *fp = OpenFileWrite(param.RefinedPath);
	fclose(fp);
	dbr::Warn("no_hits", std::string("reason=") + reason);
	dbr::StatusSetProteins(0);
}

/// 数出 DBReducer.fasta 的蛋白条数；0 条时按 no_hits 处理。返回是否有蛋白。
bool CheckReducedFasta(const Params &param) {
	long long n = 0;
	FILE *fp = fopen(param.RefinedPath.c_str(), "rb");
	if (!fp) throw DbrError(Result::output_missing, "没有写出 DBReducer.fasta：" + AnsiToUtf8(param.RefinedPath), "report_bug");
	std::vector<char> buf(1 << 20);
	bool line_start = true;
	size_t k;
	while ((k = fread(buf.data(), 1, buf.size(), fp)) > 0) {
		for (size_t i = 0; i < k; ++i) {
			if (line_start && buf[i] == '>') ++n;
			line_start = (buf[i] == '\n');
		}
	}
	fclose(fp);
	if (n <= 0) {
		dbr::Warn("no_hits", "reason=filtered_out threshold=" + Fmt("%g", param.threshold));
		dbr::StatusSetProteins(0);
		return false;
	}
	dbr::StatusSetProteins(n);
	return true;
}

/// .protein 路径（output_pro_info=1，及 option 2）。返回是否有蛋白。
bool ReduceByProtein(const Params &param) {
	dbr::StatusPhase("verify_search");
	const std::string protein = param.SearchRes + "/pFind.protein";
	if (!RequireDataLines(protein, "pFind.protein", 2, param)) {
		dbr::StatusPhase("write");
		NoHits(param, "header_only");   // 不必为 0 条结果去读整个原始库
		return false;
	}
	dbr::StatusPhase("rebuild");
	DBReducer(protein, param.ProteinDatabase, param.RefinedPath, param.threshold, param.isMaxQuant); // 重建数据库
	return CheckReducedFasta(param);
}

/// .spectra 路径（output_pro_info=0，及 option 1/3）。返回是否有蛋白。
/// full_dict=true（仅 option 3）：整库读进 name_seq，留给后面的 top2 扫全库；否则只读入要写出的蛋白（D9）。
bool ReduceBySpectra(const Params &param, std::unordered_map<std::string, std::string> *name_seq, bool full_dict) {
	dbr::StatusPhase("verify_search");
	const std::string spectra = param.SearchRes + PSEP "pFind.spectra";
	if (!RequireDataLines(spectra, "pFind.spectra", 1, param)) {
		dbr::StatusPhase("write");
		NoHits(param, "header_only");
		return false;
	}
	dbr::StatusPhase("rebuild");
	std::deque<ReduSpec::RawRes> lines; // .spectra文件读取的原始结果
	ReduSpec::ReadRawResbySpec(spectra, &lines, param.threshold);
	std::cout << "read raw end\n";
	std::deque<ReduSpec::Protein> pros; // 将鉴定结果按照protein组织
	ReduSpec::StoreProbyLine(lines, &pros);
	std::cout << "store pro end\n";
	std::unordered_map<std::string, std::string> name_full_name;
	if (full_dict)
		BuildDict(name_seq, param.ProteinDatabase, &name_full_name); // 将原始数据库存成dict形式
	ReduSpec::GroupPros(&pros); // 将蛋白分组（标记subset蛋白）
	std::cout << "group end\n";
	if (!full_dict) {
		// 分组与原始库无关，先分组，再只读入非 subset 蛋白；ReduceDB 按 pros 的顺序写出，结果不变
		std::unordered_set<std::string_view> wanted;   // 指向 pros 里的名字，不复制字符串
		for (const auto &pro : pros)
			if (!pro.is_subset)
				wanted.insert(pro.name);
		LoadSelected(wanted, param.ProteinDatabase, name_seq, &name_full_name);
	}
	ReduSpec::ReduceDB(param.RefinedPath, pros, *name_seq, name_full_name);
	return CheckReducedFasta(param);
}

/// D8b：缩库最后 BuildDict 会把整个原始库读进两张哈希表。系数 4.5 是按数据结构推算的估算值（未实测）：
/// 梨 DNA 库平均每条 125 B，读进内存后每条约 550–600 B。
/// D9 之后只有 option 3（top2 要扫全库）还整库读入；其余 option 流式读库、只留选中的蛋白，不再估算。
void CheckBuildDictMemory(const Params &param, long long db_bytes, int option) {
	if (option != 3) {
		printf("Rebuild reads the database by streaming (only the selected proteins are kept in memory)\n");
		return;
	}
	if (db_bytes <= 0) return;
	const unsigned long long avail_bytes = plat::AvailCommitBytes();
	if (avail_bytes == 0) return;
	const double gb = 1024.0 * 1024 * 1024;
	const double est = db_bytes * 4.5 / gb;
	const double avail = avail_bytes / gb;
	printf("BuildDict memory estimate: %.1f GB (db %.2f GB x 4.5, estimated), commit available %.1f GB\n", est,
	       db_bytes / gb, avail);
	dbr::StatusBuildDict(est, avail, est > avail);
	if (est > avail) {
		dbr::Warn("builddict_memory", "est_gb=" + Fmt("%.1f", est) + " commit_avail_gb=" + Fmt("%.1f", avail) +
		                                  " suggestion_code=shard_db");
		dbr::StatusSuggest("shard_db");
	}
}

void RecordParams(const Params &p, long long db_bytes) {
	dbr::StatusParamNum("multip", p.multip);
	dbr::StatusParamNum("thread", p.thread);
	dbr::StatusParamNum("maxprolen", (double)p.maxprolen);
	dbr::StatusParamNum("maxspec", (double)p.maxspec);
	dbr::StatusParamNum("search_mode", p.search_mode);
	dbr::StatusParamNum("resume", p.resume);
	dbr::StatusParamNum("stall_timeout_min", p.stall_timeout_min);
	dbr::StatusParamNum("heartbeat_sec", p.heartbeat_sec);
	dbr::StatusParamNum("threshold", p.threshold);
	dbr::StatusParamNum("output_pro_info", p.output_pro_info);
	dbr::StatusParamStr("ProteinDatabase", AnsiToUtf8(p.ProteinDatabase));
	dbr::StatusParamNum("db_bytes", (double)db_bytes);
}

// ---------------------------------------------------------------- 主流程

int Run(int argc, char **argv) {
	ReduSpec::GodelInitialize(); // Init godel coding table
	auto start = std::chrono::system_clock::now();
	std::string now = getcwd(NULL, 0);
	auto pos = now.find_last_of("/\\");
	std::string pFind_path = now.substr(0, pos) + PSEP "pFind3" PSEP "bin";
	std::cout << now << std::endl;
	std::string cfgpath = now + PSEP "DBReducer.cfg";
	int option = 0;
	if (argc > 1) {
		cfgpath = argv[1];
	}
	if (argc > 2) {
		try {
			option = std::stoi(argv[2]);
		} catch (const std::exception &) {
			throw DbrError(Result::param_error, "命令行第 2 个参数（option）不是整数：" + AnsiToUtf8(argv[2]), "check_input");
		}
	}
	Params param;
	try {
		param.LoadParam(cfgpath); // 从DBReducer.cfg文件读取参数
	} catch (const FileError &) {
		throw;
	} catch (const std::exception &e) {
		throw DbrError(Result::param_error, std::string("cfg 中有无法解析的参数值（") + e.what() + "）", "check_input");
	}
	param.RefinedPath = param.SearchRes + PSEP "DBReducer.fasta"; // 直接写死输出文件的路径

	dbr::StatusInit(param.SearchRes, option, param.heartbeat_sec);
	const long long db_bytes = FileSize(param.ProteinDatabase);
	RecordParams(param, db_bytes);
	{
		std::string stage = StripSlash(param.SearchRes);
		size_t p = stage.find_last_of("\\/");
		if (p != std::string::npos) stage = stage.substr(p + 1);
		dbr::Log("START", "stage=" + Ascii(stage) + " option=" + std::to_string(option) +
		                      " pid=" + std::to_string(plat::Pid()) + " db_bytes=" + std::to_string(db_bytes) +
		                      " multip=" + std::to_string(param.multip) + " thread=" + std::to_string(param.thread) +
		                      " maxprolen=" + std::to_string(param.maxprolen) + " maxspec=" + std::to_string(param.maxspec) +
		                      " resume=" + std::to_string(param.resume) +
		                      " stall_timeout_min=" + Fmt("%g", param.stall_timeout_min) +
		                      " heartbeat_sec=" + std::to_string(param.heartbeat_sec));
	}
	dbr::StatusPhase("init");
	if (db_bytes < 0)
		throw DbrError(Result::param_error, "找不到原始库 ProteinDatabase：" + AnsiToUtf8(param.ProteinDatabase), "check_input");
	CheckBuildDictMemory(param, db_bytes, option);

	if (option == 0) { // option为0，执行搜索
		MakeDir(param.SearchRes);
		PFindParams pFind_params(param); // 加载pFind参数
		std::unordered_set<std::string> mods; // 总体要考虑的修饰
		ReadMod(param, &mods);
		// 判断路径下是否存在modification.ini文件
		if (0 != access((pFind_path + PSEP "modification.ini").c_str(), 0)) {
			// 原逻辑在这里调 WriteMod 后 return 0；但 WriteMod 读的正是这个不存在的文件，实际以 1 退出
			throw DbrError(Result::param_error, "pFind 目录下缺少 modification.ini：" + AnsiToUtf8(pFind_path), "check_pfind");
		}
		std::string pFind_cfg_path = param.SearchRes + PSEP "pFind.cfg";
		std::cout << pFind_cfg_path << std::endl;
		pFind_params.Write(pFind_cfg_path); // 写出pFind.cfg文件
		Precheck(param);
		RunSearcher(param, pFind_path, pFind_cfg_path);
		if (param.output_pro_info == 1) { // 如果使用pFind蛋白质推断模块
			ReduceByProtein(param);
		}
		else { // 不做蛋白质推断，直接根据.spectra文件建库
			std::unordered_map<std::string, std::string> name_seq;
			ReduceBySpectra(param, &name_seq, false);
		}
	}
	else if (option == 1) { // 如果是1，不做搜索直接根据.spctra文件建库
		std::cout << "No search\n";
		std::unordered_map<std::string, std::string> name_seq;
		ReduceBySpectra(param, &name_seq, false);
	}
	else if (option == 2) { // 如果是2，不做搜索直接根据.protein文件建库
		ReduceByProtein(param);
	}
	else if (option == 3) { // 如果是3，除了根据.spectra文件建库外，还将提取FDR < param.threshold的谱图对应的top2肽段，提取其对应的蛋白质建库
		std::cout << "No search\n";
		std::unordered_map<std::string, std::string> name_seq;
		if (!ReduceBySpectra(param, &name_seq, true)) {
			std::cout << "No hits, skip top2\n";   // 0 命中：没有可信谱图，top2 无从提取
			return dbr::StatusFinish(Result::ok);
		}
		//读取所有经过过滤的spectra
		std::unordered_set<std::string> MS2s;
		ReduSpec::ReadCredibleMS2(param.SearchRes + PSEP "pFind.spectra", &MS2s, param.threshold);
		std::cout << "MS2s size: " << MS2s.size() << std::endl;
		//读取所有top2肽段
		std::vector<std::string> top2seqs;
		ReduSpec::ReadTop2Seq(param.SearchRes, MS2s, &top2seqs);
		//写出所有top2肽段
		auto out_file = OpenFileWrite(param.SearchRes + PSEP "top2.pep");
		for (auto seq : top2seqs) {
			fwrite(seq.c_str(), seq.length(), 1, out_file);
			fwrite("\n", 1, 1, out_file);
		}
		fclose(out_file);
		//根据top2肽段建库
		ReduSpec::GetTop2Fasta(name_seq, top2seqs, param.SearchRes + PSEP "top2.fasta");
	}
	else {
		MakeDir(param.SearchRes);
		PFindParams pFind_params(param); // 加载pFind参数
		std::unordered_set<std::string> mods; // 总体要考虑的修饰
		ReadMod(param, &mods);
		WriteMod(pFind_path + PSEP "modification.ini", mods, pFind_params.modpath); // 写出要考虑的修饰
		std::string pFind_cfg_path = param.SearchRes + PSEP "pFind.cfg";
		std::cout << pFind_cfg_path << std::endl;
		pFind_params.Write(pFind_cfg_path); // 写出pFind.cfg文件
		Precheck(param);
		RunSearcher(param, pFind_path, pFind_cfg_path);
		ReduceByProtein(param);
	}
	auto end = std::chrono::system_clock::now();
	printf("Database rebuild completed in %.2f seconds\n", std::chrono::duration<double>(end - start).count());
	return dbr::StatusFinish(Result::ok);
}

}  // namespace

int main(int argc, char **argv) {
	try {
		return Run(argc, argv);
	} catch (const DbrError &e) {
		return dbr::StatusFinish(e.result, e.what(), e.suggestion_code, &e.proc);
	} catch (const FileError &e) {
		return dbr::StatusFinish(Result::param_error, "无法打开文件：" + AnsiToUtf8(e.path), "check_input");
	} catch (const std::bad_alloc &) {
		return dbr::StatusFinish(Result::internal_error,
		                         std::string("DBReducer 内存分配失败（std::bad_alloc），发生在 ") + dbr::CurrentPhase() + " 阶段",
		                         "free_memory");
	} catch (const std::exception &e) {
		return dbr::StatusFinish(Result::internal_error, "DBReducer 内部异常：" + AnsiToUtf8(e.what()), "report_bug");
	} catch (...) {
		return dbr::StatusFinish(Result::internal_error, "DBReducer 内部出现未知异常", "report_bug");
	}
}

void ReadMod(const Params &pAnno_param, std::unordered_set<std::string> *mods) {
	auto segs1 = Split(pAnno_param.selectmod, ";");
	auto segs2 = Split(pAnno_param.fixmod, ";");
	for (auto seg : segs1) {
		if (seg != "")
			mods->insert(seg);
	}
	for (auto seg : segs2) {
		if (seg != "")
			mods->insert(seg);
	}
}

/**
 * @brief 从原始mod文件中挑选修饰并写出
 *
 * @param mods 待写出的修饰
 * @param raw_mod_path 原始mod文件
 * @param new_mod_path 新mod文件
 */
void WriteMod(const std::string &raw_mod_path, const std::unordered_set<std::string> &mods, const std::string &new_mod_path) {
	std::cout << mods.size() << std::endl;
	for (auto mod : mods) {
		 std::cout << mod << std::endl;
	}
	FILE* out_file = OpenFileWrite(new_mod_path);
	setvbuf ( out_file , NULL , _IOFBF , 102400 );
	std::string title = "@NUMBER_MODIFICATION=" + std::to_string(mods.size()) + "\n";
	fwrite(title.c_str(), title.length(), 1, out_file); // 写文件头
	std::ifstream buffer = ReadFilebyStream(raw_mod_path);
    std::string linestr = "";
	int write_cnt = 0;
    while(!buffer.eof() && getline(buffer, linestr)) {
		if (linestr.substr(0, 4) == "name") {
			std::string mod_name = SplitByChars(linestr, " =\t\n")[1];
			if (find(mods.begin(), mods.end(), mod_name) != mods.end()) {
				write_cnt++;
				std::string name_line = "name" + std::to_string(write_cnt) + "=" + SplitByChars(linestr, "=")[1] + "\n";
				fwrite(name_line.c_str(), name_line.length(), 1, out_file);
				if (buffer.eof()) {
					std::cout << "Please cheack the modification file!\n";
					fclose(out_file);
					throw dbr::DbrError(dbr::Result::param_error,
					                    "modification.ini 格式错误：修饰 " + AnsiToUtf8(mod_name) + " 缺少定义行", "check_pfind");
				}
				getline(buffer, linestr);
				fwrite(linestr.c_str(), linestr.length(), 1, out_file);
				fwrite("\n", 1, 1, out_file);
				if (write_cnt == mods.size())
					break;
			}
		}
    }
	buffer.close();
	fclose(out_file);
}
