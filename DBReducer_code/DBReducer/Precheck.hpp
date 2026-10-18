/**
 * @file Precheck.hpp
 * @brief 启动 pFind 前的资源预检（只告警，不改参数）。
 *
 * 估算依据见 user_questions/已知问题与修复建议.md 的 I-08 与附录（2026-10 实测）：
 *  - 条目数、氨基酸数：采样库文件前 100 MB 外推（外推误差 < 0.3%）
 *  - 可提交余额：GlobalMemoryStatusEx 的 ullAvailPageFile（不是 ullAvailPhys）
 *  - 单个 pFind worker：
 *      蛋白侧 ≈ min(maxprolen, 总氨基酸) × 24 B + 总条目数 × 100 B（首次 reserve）
 *      谱图侧 ≈ min(maxspec, 单个谱图文件的谱图数) × 60 KB（按每张 13–60 KB 的上限算）
 *        —— 一批最多装下当前这一个文件的谱图；文件谱图数取最大的那个文件采样外推
 *  - 实际 worker 数 = min(multip, 谱图文件数)：pFind 按文件把任务分给 worker
 *  - 估算值 × 实际 worker 数超过可提交余额 → 告警
 * 内存系数尚未标定，结果只用于提示。
 */
#ifndef DBR_PRECHECK_HPP
#define DBR_PRECHECK_HPP

#include <string>
#include <vector>

namespace dbr {

struct PrecheckResult {
    bool ok = false;                 ///< false 表示没能完成估算（读不了库等），不告警
    long long db_bytes = 0;
    long long entries_est = 0;
    long long aa_est = 0;
    int workers = 0;                 ///< min(multip, 谱图文件数)
    long long spectra_per_file_est = -1;  ///< 最大谱图文件的谱图数估算；-1 表示估不出（如 PF2），谱图侧按 maxspec 算
    long long spectra_batch = 0;     ///< min(maxspec, spectra_per_file_est)
    double per_worker_gb = 0;
    double est_mem_gb = 0;           ///< × workers 之后
    double commit_avail_gb = 0;
    bool warned = false;
    int suggest_multip = 0;          ///< 告警且单个 worker 装得下时：余额能开的 multip（≥1）；否则为 0，表示不给此项
    std::string suggestion_code;     ///< reduce_multip / reduce_maxprolen / reduce_maxspec
};

PrecheckResult RunPrecheck(const std::string &fasta_path, long long maxprolen, long long maxspec, int multip,
                           const std::vector<std::string> &ms_files, const std::string &ms_type);

}  // namespace dbr

#endif
