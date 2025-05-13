/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2021-12-20 22:30:44
 * @LastEditors: wangkaifei kfwang@stu.xidian.edu.cn
 * @LastEditTime: 2025-01-24 14:27:04
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
#include <deque>
#include <chrono>   
#include "../utils/IOUtil.h"
#include "../utils/StringUtil.h"
#include "Param.hpp"
#include "Rebuild.hpp"
#include "Rebuild-spec.hpp"

using namespace std;
void ReadMod(const Params &pAnno_param, std::unordered_set<std::string> *mods);
void WriteMod(const std::string &raw_mod_path, const std::unordered_set<std::string> &mods, const std::string &new_mod_path);

int main(int argc, char **argv) {
	ReduSpec::GodelInitialize(); // Init godel coding table
	auto start = std::chrono::system_clock::now();
	std::string now = getcwd(NULL, 0);
	auto pos = now.find_last_of("/\\");
	std::string pFind_path = now.substr(0, pos) + "\\pFind3\\bin";
	std::cout << now << std::endl;
	std::string cfgpath = now + "\\DBReducer.cfg";
	int option = 0;
	if (argc > 1) {
		cfgpath = argv[1];
	}
	if (argc > 2)
		option = std::stoi(argv[2]);
	Params param;
	param.LoadParam(cfgpath); // 从DBReducer.cfg文件读取参数
	param.RefinedPath = param.SearchRes + "\\DBReducer.fasta"; // 直接写死输出文件的路径
	// param.PrintParam();
	// std::cout << "option: " << option << std::endl;
	if (option == 0) { // option为0，执行搜索
		if (0 != access(param.SearchRes.c_str(), 0)) {
			if (mkdir(param.SearchRes.c_str()) == -1) {
				cout << "Path: " << param.SearchRes << endl;
				printf("Build Folder failed!\n");
				return 0;
			}
		}
		PFindParams pFind_params(param); // 加载pFind参数
		std::unordered_set<std::string> mods; // 总体要考虑的修饰
		ReadMod(param, &mods);
		// 判断路径下是否存在modification.ini文件
		if (0 != access((pFind_path + "\\modification.ini").c_str(), 0)) {
			WriteMod(pFind_path + "\\modification.ini", mods, pFind_params.modpath); // 写出要考虑的修饰
			return 0;
		}		
		std::string pFind_cfg_path = param.SearchRes + "\\pFind.cfg";
		std::cout << pFind_cfg_path << std::endl;
		pFind_params.Write(pFind_cfg_path); // 写出pFind.cfg文件
		std::string search_cmd = "Searcher.exe " + pFind_cfg_path;
		if (chdir(pFind_path.c_str()) == -1) // 切换到pFind bin目录
			printf("Please cheack your pFind installation path :\n%s\n", pFind_path.c_str());
		printf("%s\n", search_cmd.c_str());
		system(search_cmd.c_str());  // 执行命令，搜索	
		if (chdir(now.c_str()) == -1) // 切换到DBReducer目录
			printf("Error: Please cheack the path :\n%s\n", now);
		if (param.output_pro_info == 1) { // 如果使用pFind蛋白质推断模块
			DBReducer(param.SearchRes + "/pFind.protein", param.ProteinDatabase, param.RefinedPath, param.threshold, param.isMaxQuant); // 重建数据库
		}
		else { // 不做蛋白质推断，直接根据.spectra文件建库
			std::deque<ReduSpec::RawRes> lines; // .spectra文件读取的原始结果
			ReduSpec::ReadRawResbySpec(param.SearchRes + "\\pFind.spectra", &lines, param.threshold);
			std::deque<ReduSpec::Protein> pros; // 将鉴定结果按照protein组织
			ReduSpec::StoreProbyLine(lines, &pros);
			std::unordered_map<std::string, std::string> name_seq;
			std::unordered_map<std::string, std::string> name_full_name;
			BuildDict(&name_seq, param.ProteinDatabase, &name_full_name); // 将原始数据库存成dict形式
			ReduSpec::GroupPros(&pros); // 将蛋白分组（标记subset蛋白）
			ReduSpec::ReduceDB(param.RefinedPath, pros, name_seq, name_full_name);
		}
	}
	else if (option == 1) { // 如果是1，不做搜索直接根据.spctra文件建库
		std::cout << "No search\n";
		std::deque<ReduSpec::RawRes> lines; // .spectra文件读取的原始结果
		ReduSpec::ReadRawResbySpec(param.SearchRes + "\\pFind.spectra", &lines, param.threshold);
		std::cout << "read raw end\n";
		std::deque<ReduSpec::Protein> pros; // 将鉴定结果按照protein组织
		ReduSpec::StoreProbyLine(lines, &pros);
		std::cout << "store pro end\n";
		std::unordered_map<std::string, std::string> name_seq;
		std::unordered_map<std::string, std::string> name_full_name;
		BuildDict(&name_seq, param.ProteinDatabase, &name_full_name); // 将原始数据库存成dict形式
		ReduSpec::GroupPros(&pros); // 将蛋白分组（标记subset蛋白）
		std::cout << "group end\n";
		ReduSpec::ReduceDB(param.RefinedPath, pros, name_seq, name_full_name);
		return 0;
	}
	else if (option == 2) { // 如果是2，不做搜索直接根据.protein文件建库
		DBReducer(param.SearchRes + "/pFind.protein", param.ProteinDatabase, param.RefinedPath, param.threshold, param.isMaxQuant); // 重建数据库
	}
	else if (option == 3) { // 如果是3，除了根据.spectra文件建库外，还将提取FDR < param.threshold的谱图对应的top2肽段，提取其对应的蛋白质建库 
		std::cout << "No search\n";
		std::deque<ReduSpec::RawRes> lines; // .spectra文件读取的原始结果
		ReduSpec::ReadRawResbySpec(param.SearchRes + "\\pFind.spectra", &lines, param.threshold);
		std::cout << "read raw end\n";
		std::deque<ReduSpec::Protein> pros; // 将鉴定结果按照protein组织
		ReduSpec::StoreProbyLine(lines, &pros);
		std::cout << "store pro end\n";
		std::unordered_map<std::string, std::string> name_seq;
		std::unordered_map<std::string, std::string> name_full_name;
		BuildDict(&name_seq, param.ProteinDatabase, &name_full_name); // 将原始数据库存成dict形式
		ReduSpec::GroupPros(&pros); // 将蛋白分组（标记subset蛋白）
		std::cout << "group end\n";
		ReduSpec::ReduceDB(param.RefinedPath, pros, name_seq, name_full_name);
		//读取所有经过过滤的spectra
		std::unordered_set<std::string> MS2s;
		ReduSpec::ReadCredibleMS2(param.SearchRes + "\\pFind.spectra", &MS2s, param.threshold);
		std::cout << "MS2s size: " << MS2s.size() << std::endl;
		//读取所有top2肽段
		std::vector<std::string> top2seqs;
		ReduSpec::ReadTop2Seq(param.SearchRes, MS2s, &top2seqs);
		//写出所有top2肽段
		auto out_file = OpenFileWrite(param.SearchRes + "\\top2.pep");
		for (auto seq : top2seqs) {
			fwrite(seq.c_str(), seq.length(), 1, out_file);
			fwrite("\n", 1, 1, out_file);
		}
		fclose(out_file);
		//根据top2肽段建库
		ReduSpec::GetTop2Fasta(name_seq, top2seqs, param.SearchRes + "\\top2.fasta");
	}
	else {
		if (0 != access(param.SearchRes.c_str(), 0)) {
			if (mkdir(param.SearchRes.c_str()) == -1) {
				cout << "Path: " << param.SearchRes << endl;
				printf("Build Folder failed!\n");
				return 0;
			}
		}
		PFindParams pFind_params(param); // 加载pFind参数
		std::unordered_set<std::string> mods; // 总体要考虑的修饰
		ReadMod(param, &mods);
		WriteMod(pFind_path + "\\modification.ini", mods, pFind_params.modpath); // 写出要考虑的修饰
		std::string pFind_cfg_path = param.SearchRes + "\\pFind.cfg";
		std::cout << pFind_cfg_path << std::endl;
		pFind_params.Write(pFind_cfg_path); // 写出pFind.cfg文件
		std::string search_cmd = "Searcher.exe " + pFind_cfg_path;
		if (chdir(pFind_path.c_str()) == -1) // 切换到pFind bin目录
			printf("Please cheack your pFind installation path :\n%s\n", pFind_path.c_str());
		printf("%s\n", search_cmd.c_str());
		system(search_cmd.c_str());  // 执行命令，搜索	
		if (chdir(now.c_str()) == -1) // 切换到DBReducer目录
			printf("Error: Please cheack the path :\n%s\n", now);
		DBReducer(param.SearchRes + "/pFind.protein", param.ProteinDatabase, param.RefinedPath, param.threshold, param.isMaxQuant); // 重建数据库
	}
	auto end = std::chrono::system_clock::now();
	printf("Database rebuild completed in %.2f seconds\n", std::chrono::duration<double>{end-start});
	return 0;
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
					exit(1);
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