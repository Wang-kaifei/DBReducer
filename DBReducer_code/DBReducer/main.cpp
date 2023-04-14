/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2021-12-20 22:30:44
 * @LastEditors: Kaifei
 * @LastEditTime: 2023-04-14 10:38:18
 */
#include <iostream>
#include <windows.h>
#include <string>
#include <vector>
#include <string>
#include <stdlib.h>
#include <direct.h>
#include <sstream>
#include <stdio.h>
#include <direct.h>
#include <chrono>   
#include "../utils/IOUtil.h"
#include "../utils/StringUtil.h"
#include "Param.hpp"
#include "Rebuild.hpp"

using namespace std;
void ReadMod(const Params &pAnno_param, std::set<std::string> *mods);
void WriteMod(const std::string &all_mod_path, const std::set<std::string> &sub_mods, const std::string &out_path);

int main(int argc, char **argv) {
	auto start = std::chrono::system_clock::now();
	std::string now = getcwd(NULL, 0);
	auto pos = now.find_last_of("/\\");
	std::string pFind_path = now.substr(0, pos) + "\\pFind3\\bin";
	std::cout << now << std::endl;
	std::string cfgpath = now + "\\DBReducer.cfg";
	if (argc > 1) {
		cfgpath = argv[1];
	}
	Params param;
	param.LoadParam(cfgpath); // 从DBReducer.cfg文件读取参数
	param.RefinedPath = param.SearchRes + "\\DBReducer.fasta"; // 直接写死输出文件的路径
	param.PrintParam();
	if (0 != access(param.SearchRes.c_str(), 0)) {
		if (mkdir(param.SearchRes.c_str()) == -1) {
			cout << "Path: " << param.SearchRes << endl;
			printf("Build Folder failed!\n");
			return 0;
		}
    }
	PFindParams pFind_params(param); // 加载pFind参数
	std::set<std::string> mods; // 总体要考虑的修饰
	ReadMod(param, &mods);
	WriteMod(pFind_path + "\\modification.ini", mods, pFind_params.modpath); // 写出要考虑的修饰
	std::string pFind_cfg_path = param.SearchRes + "\\pFind.cfg";
	pFind_params.Write(pFind_cfg_path); // 写出pFind.cfg文件
	std::string search_cmd = "Searcher.exe " + pFind_cfg_path;
	if (chdir(pFind_path.c_str()) == -1) // 切换到pFind bin目录
		printf("Please cheack your pFind installation path :\n%s\n", pFind_path.c_str());
	printf("%s\n", search_cmd.c_str());
	system(search_cmd.c_str());  // 执行命令，搜索	
	if (chdir(now.c_str()) == -1) // 切换到pAnno目录
		printf("Error: Please cheack the path :\n%s\n", now);
	DBReducer(param.SearchRes + "/pFind.protein", param.ProteinDatabase, param.RefinedPath, param.threshold, param.isMaxQuant); // 重建数据库
	auto end   = std::chrono::system_clock::now();
	printf("Database rebuild completed in %.2f seconds\n", std::chrono::duration<double>{end-start});
	return 0;
}

void ReadMod(const Params &pAnno_param, std::set<std::string> *mods) {
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

void WriteMod(const std::string &all_mod_path, const std::set<std::string> &sub_mods, const std::string &out_path) {
	int32_t cnt = 0;
	std::ifstream buffer = ReadFilebyStream(all_mod_path);
	FILE * refined_mod_file = OpenFileWrite(out_path);
    setvbuf ( refined_mod_file , NULL , _IOFBF , 1024000 );
	std::string linestr;
    while (!buffer.eof()) {
        getline(buffer, linestr);
		if ("name" == linestr.substr(0, 4)) {
			auto segs = SplitByChars(linestr, "= ");
			if (segs.size() < 3) {
				std::cout << "Please cheack the modification: " << linestr << std::endl;
				continue;
			}
			if (sub_mods.find(segs[1]) != sub_mods.end()) { // 应该被写入
				fwrite(linestr.c_str(), linestr.length(), 1, refined_mod_file); // 写入修饰名
				fwrite("\n", 1, 1, refined_mod_file);
				if (buffer.eof()) {
					std::cout << "Please cheack the modification: " << linestr << std::endl;
					break;
				}
				getline(buffer, linestr);
				fwrite(linestr.c_str(), linestr.length(), 1, refined_mod_file); // 写入修饰内容
				fwrite("\n", 1, 1, refined_mod_file);
				if (++cnt >= sub_mods.size())
					break;
			}
		}
    }
	fclose(refined_mod_file);
}