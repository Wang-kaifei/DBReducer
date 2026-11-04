/**
 * @file Platform.h
 * @brief Windows / Linux 两个平台不同的地方集中在这里（2026-11，Linux 移植）。只有头文件，不改各 DLL 的导出。
 *
 *  - PSEP / PSEP_C：路径分隔符。源码里的路径一律写成  dir + PSEP "pFind_res1" PSEP "pFind.cfg"  这种形式，
 *    字面量在编译期拼接，Windows 上得到的字符串与原来写死反斜杠时逐字节相同。
 *  - EXE：可执行文件后缀，"Searcher" EXE 在 Windows 上就是 "Searcher.exe"。
 *  - plat::：文件大小 / 修改时间 / 列目录 / 改名覆盖 / 复制 / 计时 / 程序自身位置等。
 *    Windows 版用的就是原先各处直接调用的那些 API，行为不变；Linux 版用 POSIX。
 *  - 文件时间（FileTime）两个平台都用 Windows FILETIME 的单位和起点（100 ns，自 1601-01-01 UTC），
 *    上层比较 / 相减的代码不用分平台。
 *
 * 原件在 pAnno 仓库 utils/，DBReducer 仓库 DBReducer_code/utils/ 里有一份相同的副本，改动请同步。
 */
#ifndef PANNO_PLATFORM_H
#define PANNO_PLATFORM_H

#include <string>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#define PSEP "\\"
#define PSEP_C '\\'
#define EXE ".exe"
#else
#include <dirent.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <time.h>
#include <unistd.h>
#include <cerrno>
#include <climits>
#define PSEP "/"
#define PSEP_C '/'
#define EXE ""
#endif

namespace plat {

using FileTime = unsigned long long;   ///< 100 ns，自 1601-01-01 UTC（即 Windows FILETIME）

struct Entry {
    std::string name;      ///< 只有文件名，不含目录
    bool is_dir = false;
    long long size = 0;
    FileTime mtime = 0;
};

#ifdef _WIN32
// ================================================================ Windows

inline FileTime ToFileTime(const FILETIME &ft) {
    return ((unsigned long long)ft.dwHighDateTime << 32) | ft.dwLowDateTime;
}

/// 文件或目录的属性。不存在返回 false
inline bool Stat(const std::string &path, Entry *e) {
    WIN32_FILE_ATTRIBUTE_DATA d;
    if (!GetFileAttributesExA(path.c_str(), GetFileExInfoStandard, &d)) return false;
    e->is_dir = (d.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
    e->size = ((long long)d.nFileSizeHigh << 32) | d.nFileSizeLow;
    e->mtime = ToFileTime(d.ftLastWriteTime);
    return true;
}

inline FileTime NowFileTime() {
    FILETIME ft;
    GetSystemTimeAsFileTime(&ft);
    return ToFileTime(ft);
}

/// 目录下的条目（不含 . 和 ..），顺序同 FindFirstFile（NTFS 上按名字）
inline std::vector<Entry> ListDir(const std::string &dir) {
    std::vector<Entry> out;
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA((dir + PSEP "*").c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) return out;
    do {
        std::string n = fd.cFileName;
        if (n == "." || n == "..") continue;
        Entry e;
        e.name = n;
        e.is_dir = (fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        e.size = ((long long)fd.nFileSizeHigh << 32) | fd.nFileSizeLow;
        e.mtime = ToFileTime(fd.ftLastWriteTime);
        out.push_back(e);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
    return out;
}

/// 用 src 原子地覆盖 dst（src 一般是同目录的 .tmp）
inline bool ReplaceFile(const std::string &src, const std::string &dst) {
    return MoveFileExA(src.c_str(), dst.c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) != 0;
}
/// 复制并覆盖。与 CopyFile 一样保留源文件的修改时间
inline bool CopyFileOver(const std::string &src, const std::string &dst) {
    return CopyFileA(src.c_str(), dst.c_str(), FALSE) != 0;
}
inline bool RemoveFile(const std::string &path) { return DeleteFileA(path.c_str()) != 0; }

inline unsigned long long TickMs() { return GetTickCount64(); }
inline void SleepMs(unsigned ms) { Sleep(ms); }
inline unsigned long Pid() { return GetCurrentProcessId(); }
inline unsigned long LastError() { return GetLastError(); }

/// 本程序（exe）所在的目录，末尾不带分隔符
inline std::string ExeDir() {
    char buf[MAX_PATH * 4];
    DWORD n = GetModuleFileNameA(nullptr, buf, sizeof(buf));
    std::string p = (n > 0 && n < sizeof(buf)) ? std::string(buf, n) : std::string();
    size_t k = p.find_last_of("\\/");
    return k == std::string::npos ? p : p.substr(0, k);
}

inline bool StdinIsTty() { return _isatty(_fileno(stdin)) != 0; }
/// 双击运行时让窗口停住，看得到结果
inline void PauseConsole() { system("pause"); }
inline void UnsetEnv(const char *name) { SetEnvironmentVariableA(name, nullptr); }

/// 本机物理内存（字节），取不到返回 0
inline unsigned long long TotalPhysBytes() {
    MEMORYSTATUSEX ms;
    ms.dwLength = sizeof(ms);
    return GlobalMemoryStatusEx(&ms) ? ms.ullTotalPhys : 0;
}

/// 还能提交的内存（字节）：可提交上限减已提交（ullAvailPageFile）。取不到返回 0
inline unsigned long long AvailCommitBytes() {
    MEMORYSTATUSEX ms;
    ms.dwLength = sizeof(ms);
    return GlobalMemoryStatusEx(&ms) ? ms.ullAvailPageFile : 0;
}

/// 64 位偏移的 fseek（大于 2 GB 的文件）
inline int Seek64(FILE *fp, long long off, int whence) { return _fseeki64(fp, off, whence); }

/// 绝对路径（不要求存在）
inline std::string FullPath(const std::string &p) {
    char buf[MAX_PATH * 4];
    DWORD n = GetFullPathNameA(p.c_str(), sizeof(buf), buf, nullptr);
    return (n > 0 && n < sizeof(buf)) ? std::string(buf, n) : p;
}
constexpr bool kPathsCaseInsensitive = true;

#else
// ================================================================ Linux / POSIX

constexpr unsigned long long kUnixEpochAsFileTime = 116444736000000000ULL;   // 1970-01-01 距 1601-01-01，单位 100 ns

inline FileTime ToFileTime(const struct timespec &ts) {
    return kUnixEpochAsFileTime + (unsigned long long)ts.tv_sec * 10000000ULL + (unsigned long long)ts.tv_nsec / 100ULL;
}

inline bool Stat(const std::string &path, Entry *e) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    e->is_dir = S_ISDIR(st.st_mode);
    e->size = (long long)st.st_size;
    e->mtime = ToFileTime(st.st_mtim);
    return true;
}

inline FileTime NowFileTime() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ToFileTime(ts);
}

/// 目录下的条目（不含 . 和 ..），按名字排序（readdir 本身无序；Windows 上 NTFS 返回的也是按名字）
inline std::vector<Entry> ListDir(const std::string &dir) {
    std::vector<Entry> out;
    DIR *d = opendir(dir.c_str());
    if (!d) return out;
    while (struct dirent *de = readdir(d)) {
        std::string n = de->d_name;
        if (n == "." || n == "..") continue;
        Entry e;
        if (!Stat(dir + PSEP + n, &e)) continue;   // 刚好被删掉
        e.name = n;
        out.push_back(e);
    }
    closedir(d);
    std::sort(out.begin(), out.end(), [](const Entry &a, const Entry &b) { return a.name < b.name; });
    return out;
}

inline bool ReplaceFile(const std::string &src, const std::string &dst) {
    return rename(src.c_str(), dst.c_str()) == 0;
}
/// 复制并覆盖，保留源文件的修改时间（与 Windows 的 CopyFile 一致）。只用 POSIX，DBReducer 的 utils 按 C++14 编也能用
inline bool CopyFileOver(const std::string &src, const std::string &dst) {
    const int in = open(src.c_str(), O_RDONLY);
    if (in < 0) return false;
    struct stat st;
    if (fstat(in, &st) != 0) { close(in); return false; }
    const int out = open(dst.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (out < 0) { close(in); return false; }
    std::vector<char> buf(1 << 20);
    bool ok = true;
    while (true) {
        const ssize_t n = read(in, buf.data(), buf.size());
        if (n < 0) { if (errno == EINTR) continue; ok = false; break; }
        if (n == 0) break;
        for (ssize_t done = 0; done < n;) {
            const ssize_t w = write(out, buf.data() + done, (size_t)(n - done));
            if (w < 0) { if (errno == EINTR) continue; ok = false; break; }
            done += w;
        }
        if (!ok) break;
    }
    close(in);
    if (close(out) != 0) ok = false;
    if (ok) {
        struct timespec t[2] = {st.st_atim, st.st_mtim};
        utimensat(AT_FDCWD, dst.c_str(), t, 0);
    }
    return ok;
}
inline bool RemoveFile(const std::string &path) { return unlink(path.c_str()) == 0; }

inline unsigned long long TickMs() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (unsigned long long)ts.tv_sec * 1000ULL + (unsigned long long)ts.tv_nsec / 1000000ULL;
}
inline void SleepMs(unsigned ms) { usleep((useconds_t)ms * 1000); }
inline unsigned long Pid() { return (unsigned long)getpid(); }
inline unsigned long LastError() { return (unsigned long)errno; }

/// 本程序所在的目录（/proc/self/exe，已解开符号链接），末尾不带分隔符
inline std::string ExeDir() {
    char buf[PATH_MAX];
    ssize_t n = readlink("/proc/self/exe", buf, sizeof(buf) - 1);
    std::string p = n > 0 ? std::string(buf, (size_t)n) : std::string();
    size_t k = p.find_last_of('/');
    return k == std::string::npos ? p : (k == 0 ? std::string("/") : p.substr(0, k));
}

inline bool StdinIsTty() { return isatty(STDIN_FILENO) != 0; }
/// Linux 终端里运行完不停住（没有 pause，也不需要）
inline void PauseConsole() {}
inline void UnsetEnv(const char *name) { unsetenv(name); }

inline unsigned long long TotalPhysBytes() {
    long pages = sysconf(_SC_PHYS_PAGES), page = sysconf(_SC_PAGESIZE);
    return pages > 0 && page > 0 ? (unsigned long long)pages * (unsigned long long)page : 0;
}

/// 还能申请的内存（字节）。Linux 默认允许超额申请（overcommit），Windows 意义上的提交上限不适用，
/// 取 /proc/meminfo 的 MemAvailable + SwapFree 作对应量。取不到返回 0
inline unsigned long long AvailCommitBytes() {
    FILE *fp = fopen("/proc/meminfo", "r");
    if (!fp) return 0;
    char line[256];
    unsigned long long avail_kb = 0, swap_kb = 0;
    while (fgets(line, sizeof(line), fp)) {
        unsigned long long v;
        if (sscanf(line, "MemAvailable: %llu kB", &v) == 1) avail_kb = v;
        else if (sscanf(line, "SwapFree: %llu kB", &v) == 1) swap_kb = v;
    }
    fclose(fp);
    return (avail_kb + swap_kb) * 1024ULL;
}

inline int Seek64(FILE *fp, long long off, int whence) { return fseeko(fp, (off_t)off, whence); }

/// 绝对路径（不要求存在，不解开符号链接，与 Windows 的 GetFullPathName 一样只做字面上的规整）
inline std::string FullPath(const std::string &p) {
    std::string s = p;
    if (s.empty() || s[0] != '/') {
        char buf[PATH_MAX];
        if (!getcwd(buf, sizeof(buf))) return p;
        s = std::string(buf) + "/" + s;
    }
    std::vector<std::string> parts;
    size_t i = 0;
    while (i <= s.size()) {
        size_t j = s.find('/', i);
        if (j == std::string::npos) j = s.size();
        const std::string seg = s.substr(i, j - i);
        if (seg == "..") { if (!parts.empty()) parts.pop_back(); }
        else if (!seg.empty() && seg != ".") parts.push_back(seg);
        i = j + 1;
    }
    std::string out;
    for (const auto &seg : parts) out += "/" + seg;
    return out.empty() ? "/" : out;
}
constexpr bool kPathsCaseInsensitive = false;

#endif
// ================================================================ 两个平台共用

inline bool Exists(const std::string &p) { Entry e; return Stat(p, &e); }
inline bool IsDir(const std::string &p) { Entry e; return Stat(p, &e) && e.is_dir; }
inline bool IsFile(const std::string &p) { Entry e; return Stat(p, &e) && !e.is_dir; }
/// 普通文件的大小；不存在或是目录返回 -1
inline long long FileSize(const std::string &p) { Entry e; return Stat(p, &e) && !e.is_dir ? e.size : -1; }
/// 修改时间；不存在返回 0
inline FileTime FileMTime(const std::string &p) { Entry e; return Stat(p, &e) ? e.mtime : 0; }

}  // namespace plat

#endif
