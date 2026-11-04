/**
 * @file Precheck.cpp
 * @brief 见 Precheck.hpp。
 */
#include "../utils/Platform.h"

#include "Precheck.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <vector>

namespace dbr {

namespace {
const long long kSampleBytes = 100LL * 1024 * 1024;   // 库采样
const long long kMsSampleBytes = 32LL * 1024 * 1024;  // 谱图文件采样
const double kBytesPerAA = 24.0;             // 蛋白侧，每个氨基酸
const double kBytesPerEntry = 100.0;         // 首次 reserve，每个条目
const double kBytesPerSpectrum = 60.0 * 1024;  // 谱图侧上限
const double kGB = 1024.0 * 1024 * 1024;

long long FileBytes(const std::string &path) {
    plat::Entry e;
    if (!plat::Stat(path, &e)) return -1;
    return e.is_dir ? 0 : e.size;   // 原先用 GetFileAttributesEx，目录的大小是 0
}

/// 统计以 marker 开头的行数，按文件大小外推。失败返回 -1。
/// 在文件里均匀取 8 段采样（谱图大小随液相梯度变化，只采开头会偏：yeast 实测只采开头低估 24%）。
long long EstimateLinesStartingWith(const std::string &path, const char *marker) {
    const long long size = FileBytes(path);
    if (size <= 0) return size == 0 ? 0 : -1;
    FILE *fp = fopen(path.c_str(), "rb");
    if (!fp) return -1;
    const size_t mlen = strlen(marker);
    const int kSegments = 8;
    const long long seg_bytes = kMsSampleBytes / kSegments;
    std::vector<char> buf((size_t)seg_bytes);
    long long sampled = 0, hits = 0;
    for (int s = 0; s < kSegments; ++s) {
        const long long off = size <= kMsSampleBytes ? (long long)s * seg_bytes : size / kSegments * s;
        if (off >= size) break;
        if (plat::Seek64(fp, off, SEEK_SET) != 0) break;
        size_t n = fread(buf.data(), 1, buf.size(), fp);
        if (n == 0) break;
        // 段起点不在文件开头时，从下一个行首开始数
        size_t i = 0;
        if (off > 0) {
            while (i < n && buf[i] != '\n') ++i;
            ++i;
        }
        const size_t begin = i;
        bool line_start = true;
        for (; i + mlen <= n; ++i) {
            if (line_start && memcmp(buf.data() + i, marker, mlen) == 0) ++hits;
            line_start = (buf[i] == '\n');
        }
        if (n > begin) sampled += (long long)(n - begin);
        if (size <= kMsSampleBytes) break;   // 小文件第一段已经读完整个文件
    }
    fclose(fp);
    if (sampled <= 0) return -1;
    return (long long)(hits * ((double)size / (double)sampled));
}

/// 单个谱图文件的谱图数：MGF 数 "BEGIN IONS"，MS2 数 "S\t"；其它格式（PF2 等二进制）估不出
long long EstimateSpectra(const std::string &path, const std::string &ms_type) {
    std::string t = ms_type;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c) { return (char)std::toupper(c); });
    if (t == "MGF") return EstimateLinesStartingWith(path, "BEGIN IONS");
    if (t == "MS2") return EstimateLinesStartingWith(path, "S\t");
    return -1;
}
}  // namespace

PrecheckResult RunPrecheck(const std::string &fasta_path, long long maxprolen, long long maxspec, int multip,
                           const std::vector<std::string> &ms_files, const std::string &ms_type) {
    PrecheckResult r;
    r.db_bytes = FileBytes(fasta_path);
    if (r.db_bytes <= 0) return r;

    // ---- 采样库文件前 100 MB：数条目和氨基酸 ----
    FILE *fp = fopen(fasta_path.c_str(), "rb");
    if (!fp) return r;
    std::vector<char> buf(4 * 1024 * 1024);
    long long sampled = 0, entries = 0, aa = 0;
    bool line_start = true, in_header = false;
    while (sampled < kSampleBytes) {
        size_t want = (size_t)std::min<long long>((long long)buf.size(), kSampleBytes - sampled);
        size_t n = fread(buf.data(), 1, want, fp);
        if (n == 0) break;
        for (size_t i = 0; i < n; ++i) {
            char c = buf[i];
            if (c == '\n') {
                line_start = true;
                in_header = false;
                continue;
            }
            if (line_start && c == '>') {
                ++entries;
                in_header = true;
            } else if (!in_header && ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z'))) {
                ++aa;
            }
            line_start = false;
        }
        sampled += (long long)n;
    }
    fclose(fp);
    if (sampled <= 0) return r;
    const double scale = (double)r.db_bytes / (double)sampled;
    r.entries_est = (long long)(entries * scale);
    r.aa_est = (long long)(aa * scale);

    // ---- 实际 worker 数、单个文件的谱图数（取最大的文件） ----
    if (multip < 1) multip = 1;
    r.workers = ms_files.empty() ? multip : std::min<int>(multip, (int)ms_files.size());
    std::string largest;
    long long largest_bytes = -1;
    for (const auto &f : ms_files) {
        long long b = FileBytes(f);
        if (b > largest_bytes) {
            largest_bytes = b;
            largest = f;
        }
    }
    if (!largest.empty()) r.spectra_per_file_est = EstimateSpectra(largest, ms_type);
    const long long ms = maxspec > 0 ? maxspec : 10000;
    r.spectra_batch = r.spectra_per_file_est >= 0 ? std::min(ms, r.spectra_per_file_est) : ms;

    // ---- 可提交余额 ----
    const unsigned long long avail = plat::AvailCommitBytes();
    if (avail == 0) return r;
    r.commit_avail_gb = avail / kGB;

    // ---- 估算 ----
    const long long batch_aa = maxprolen > 0 ? std::min(maxprolen, r.aa_est) : r.aa_est;
    const double protein_side = batch_aa * kBytesPerAA + r.entries_est * kBytesPerEntry;
    const double spectra_side = r.spectra_batch * kBytesPerSpectrum;
    r.per_worker_gb = (protein_side + spectra_side) / kGB;
    r.est_mem_gb = r.per_worker_gb * r.workers;
    r.ok = true;
    if (r.est_mem_gb > r.commit_avail_gb) {
        r.warned = true;
        if (r.per_worker_gb <= r.commit_avail_gb) {
            // 单个 worker 装得下：降低 multip 即可
            r.suggest_multip = std::max(1, (int)(r.commit_avail_gb / r.per_worker_gb));
            r.suggestion_code = "reduce_multip";
        } else {
            // 一个 worker 都装不下：看哪一侧占大头，不给 suggest_multip
            r.suggest_multip = 0;
            r.suggestion_code = spectra_side > protein_side ? "reduce_maxspec" : "reduce_maxprolen";
        }
    }
    return r;
}

}  // namespace dbr
