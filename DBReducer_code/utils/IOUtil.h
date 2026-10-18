/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2022-06-27 13:26:04
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-06-27 20:10:19
 */
/**
 * @file IOUtil.h
 * @author Wang Kaifei (wangkaifei20g@ict.ac.cn)
 * @brief 声明了一些文件操作函数：打开文件：读/写；读取文件全部内容；最终结果写出
 * @version 0.1
 * @date 2021-11-23
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#ifndef IOUtil_NAME
#define IOUtil_NAME
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

// 2026-10：打不开文件时原先直接 exit(1)，调用方没有机会写状态文件。
// 现改为抛 FileError，由 main 统一收尾（退出码仍为 1 = param_error）。提示信息照旧打印。
class FileError : public std::runtime_error {
public:
	explicit FileError(const std::string & path) : std::runtime_error(path), path(path) {}
	std::string path;
};

inline FILE * OpenFileRead(const std::string & filename) {
	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == NULL) {
		std::cout << "<DBReducer> " << filename << ": Cannot find that file!\n";
		throw FileError(filename);
	}
	return fp;
}

inline std::ifstream ReadFilebyStream(const std::string & filename) {
	std::ifstream buffer(filename);
	if (!buffer.is_open()) {
		printf("Invalid path: %s\n", filename.c_str());
		throw FileError(filename);
	}
	return buffer;
}


inline FILE * OpenFileWrite(const std::string & filename) {
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		std::cout << "<DBReducer> " << filename << ": Cannot find that file!\n";
		throw FileError(filename);
	}
	return fp;
}

inline FILE * OpenFileAdd(const std::string & filename) {
	FILE *fp = fopen(filename.c_str(), "ab");
	if (fp == NULL) {
		std::cout << "<DBReducer> " << filename << ": Cannot find that file!\n";
		throw FileError(filename);
	}
	return fp;
}

inline void ClearFile(const std::string & filename) {
	FILE *fp = OpenFileWrite(filename);
	fclose(fp);
}

char *ReadFullFile(const std::string & filename, long *flen);
std::vector<std::string> getFilesList(const std::string &dirpath);

#endif