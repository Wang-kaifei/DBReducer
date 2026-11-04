/**
 * @file IOUtil.cpp
 * @author Wang Kaifei (wangkaifei20g@ict.ac.cn)
 * @brief 定义了一些文件操作函数：读取文件全部内容；最终结果写出
 * @version 0.1
 * @date 2021-11-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include "IOUtil.h"
#include "Platform.h"
#include <dirent.h>
#include <iostream>
#include <memory.h>
#include <sys/stat.h>

/**
 * @brief 读取全部文件到内存
 * 
 * @param fp 文件指针
 * @param flen 存储缓存区大小
 * @return char* 返回缓存指针
 */
char *ReadFullFile(const std::string & filename, long *flen) {
	FILE *fp = OpenFileRead(filename);
	fseek(fp, 0L, SEEK_END);
	*flen = ftell(fp);
	char *buffer = (char *)malloc(*flen + 1);
	memset(buffer, 0, *flen + 1);
	if (!buffer) {
		std::cout << "NO enough memory!";
		fclose(fp);
		return NULL;
	}
	fseek(fp, 0L, SEEK_SET);
	fread(buffer, *flen, 1, fp);
	fclose(fp);
	return buffer;
}

/**
 * @brief Get the Files List object
 * 
 * @param dirpath 读取的文件夹路径
 * @return vector<string> 文件夹中的所有文件
 */
std::vector<std::string> getFilesList(const std::string &dirpath) {
	std::cout << "OPEN: " << dirpath << std::endl;
	std::vector<std::string> allPath;
	plat::Entry self;
	if (plat::Stat(dirpath, &self) && !self.is_dir) { // 如果本身就是文件，直接返回
		allPath.emplace_back(dirpath);
		return allPath;
	}
	// 2026-11（Linux 移植，与 pAnno 仓库 utils 同步）：原先先 readdir 一次再进循环（丢掉第一个条目，Windows 上恰好总是 "."），
	// 判断子目录时 stat 只传文件名（相对当前目录）。Linux 上 readdir 无序，会丢掉真实文件；改用 plat::ListDir。
	for (const auto &e : plat::ListDir(dirpath)) {
		const std::string p = dirpath + PSEP + e.name;
		if (e.is_dir) { // 如果是文件夹
			std::vector<std::string> tempPath = getFilesList(p);
			allPath.insert(allPath.end(), tempPath.begin(), tempPath.end());
		}
		else { // 如果是文件
			allPath.push_back(p);
		}
	}
	return allPath;
}