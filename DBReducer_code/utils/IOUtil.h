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

inline FILE * OpenFileRead(const std::string & filename) {
	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == NULL) {
		std::cout << "<DBReducer> " << filename << ": Cannot find that file!\n";
		exit(1);
	}
	return fp;
}

inline std::ifstream ReadFilebyStream(const std::string & filename) {
	std::ifstream buffer(filename);
	if (!buffer.is_open()) {
		printf("Invalid path: %s\n", filename.c_str());
        exit(1);
	}
	return buffer;
}


inline FILE * OpenFileWrite(const std::string & filename) {
	FILE *fp = fopen(filename.c_str(), "wb");
	if (fp == NULL) {
		std::cout << "<DBReducer> " << filename << ": Cannot find that file!\n";
		exit(1);
	}
	return fp;
}

inline FILE * OpenFileAdd(const std::string & filename) {
	FILE *fp = fopen(filename.c_str(), "ab");
	if (fp == NULL) {
		std::cout << "<DBReducer> " << filename << ": Cannot find that file!\n";
		exit(1);
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