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
	struct stat pre_test;
	if (stat(dirpath.c_str(), &pre_test) == 0) { // 如果本身就是文件，直接返回
		if (pre_test.st_mode & S_IFREG) {
			allPath.emplace_back(dirpath);
			return allPath;
		}
	}
	DIR *dir = opendir(dirpath.c_str());
	struct dirent *entry = readdir(dir);
	while ((entry = readdir(dir)) != NULL) // 读取成功则返回下一个目录点
	{
		struct stat s;
		stat(entry->d_name, &s);
		if (s.st_mode & S_IFDIR) { // 如果是文件夹
			if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
				continue;
			std::string dirNew = dirpath + "\\" + entry->d_name;
			std::vector<std::string> tempPath = getFilesList(dirNew);
			allPath.insert(allPath.end(), tempPath.begin(), tempPath.end());
		}
		else { // 如果是文件
			std::string name = entry->d_name;
			std::string imgdir = dirpath + "\\" + name;
			allPath.push_back(imgdir);
		}
	}
	closedir(dir);
	return allPath;
}