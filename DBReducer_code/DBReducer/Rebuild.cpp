/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2022-06-27 13:24:46
 * @LastEditors: wangkaifei kfwang@stu.xidian.edu.cn
 * @LastEditTime: 2025-03-12 11:21:01
 */

#include "Param.hpp"
#include <typeinfo>
#include <iostream>
#include <set>
#include <algorithm>
#include <exception>
#include <unordered_map>
#include <unordered_set>
#include <string_view>
#include "../utils/IOUtil.h"
#include "../utils/StringUtil.h"
#include "Status.hpp"
#define VNAME(value) (#value)

void ReadProteinName(std::set<std::string> *protein_names, const std::string &protein_path, float threshold) {
    std::ifstream buffer = ReadFilebyStream(protein_path);
    std::string linestr;
    while (!buffer.eof()) {
        getline(buffer, linestr);
        auto segs = SplitRule(linestr, "\t");
        if (segs.size() < 4)
            continue;
        if (isFloat(Strip(&segs[0]))) { // 存储 lead protein
            if (std::stof(segs[3]) < threshold) {
                if ("REV" != Strip(&segs[1]).substr(0, 3)) {
                    protein_names->insert(segs[1]);
                }
            }
            else
                break;
        }
        else if (isFloat(Strip(&segs[3])) && "SameSet" == Strip(&segs[1]) && "REV" != Strip(&segs[2]).substr(0, 3)) {
            protein_names->insert(Strip(&segs[2]));
        }
    }
}

void BuildDict(std::unordered_map<std::string, std::string> *name_seq, const std::string &fasta_path, std::unordered_map<std::string, std::string> *name_full_name) {
    std::ifstream buffer = ReadFilebyStream(fasta_path);
    std::string linestr;
    std::string name = "", seq = "";
    while (!buffer.eof()) {
        getline(buffer, linestr);
        Strip(&linestr);
        if (linestr.empty()) continue; // 避免空行处理
        if (linestr[0] == '>') { // 如果是蛋白名
            if (seq != "")
                name_seq->insert(std::make_pair(name, seq));
            int pos = linestr.size();
            if (linestr.find_first_of(" \t") != std::string::npos) {
                pos = linestr.find_first_of(" \t");
            }
            auto pos2 = linestr.find_first_of("|");
            if (pos2 != std::string::npos && pos2 > 15)
                pos = pos2;
            name = linestr.substr(1, pos - 1);
            seq = "";
            name_full_name->insert(std::make_pair(name, linestr)); // 存储mini name 2 full name的字典
        }
        else
            seq += linestr;
    }
    // 处理最后一个蛋白质序列
    if (!name.empty() && !seq.empty()) {
        name_seq->insert(std::make_pair(name, seq));
    }
    std::cout << "fasta size: " << name_seq->size() << std::endl;
}

/// D9：同 BuildDict 的解析规则，但只保留 wanted 里的蛋白，其余条目读过即丢，内存只与选中的蛋白有关。
/// 两张表的内容与"BuildDict 整库读入后再按 wanted 取"完全相同（同样 first-wins：全名取第一次出现的表头，
/// 序列取第一个非空序列），因此下游按名字取出的结果不变。
/// wanted 只用来查名字是否在名单里（count），直接用调用方已有的名单，不另复制一份。
template <class NameSet>
void LoadSelected(const NameSet &wanted, const std::string &fasta_path,
                  std::unordered_map<std::string, std::string> *name_seq, std::unordered_map<std::string, std::string> *name_full_name) {
    std::ifstream buffer = ReadFilebyStream(fasta_path);
    std::string linestr;
    std::string name = "", seq = "";
    bool keep = wanted.count(name) > 0; // 首个表头之前的序列行归在空名下，同 BuildDict
    size_t scanned = 0;
    while (!buffer.eof()) {
        getline(buffer, linestr);
        Strip(&linestr);
        if (linestr.empty()) continue;
        if (linestr[0] == '>') {
            if (keep && seq != "")
                name_seq->insert(std::make_pair(name, seq));
            int pos = linestr.size();
            if (linestr.find_first_of(" \t") != std::string::npos) {
                pos = linestr.find_first_of(" \t");
            }
            auto pos2 = linestr.find_first_of("|");
            if (pos2 != std::string::npos && pos2 > 15)
                pos = pos2;
            name = linestr.substr(1, pos - 1);
            seq = "";
            ++scanned;
            keep = wanted.count(name) > 0;
            if (keep)
                name_full_name->insert(std::make_pair(name, linestr));
        }
        else if (keep)
            seq += linestr;
    }
    if (keep && !name.empty() && !seq.empty()) {
        name_seq->insert(std::make_pair(name, seq));
    }
    std::cout << "fasta scanned: " << scanned << " entries, selected " << wanted.size() << " proteins, loaded "
              << name_seq->size() << std::endl;
}
template void LoadSelected(const std::set<std::string> &, const std::string &, std::unordered_map<std::string, std::string> *,
                           std::unordered_map<std::string, std::string> *);
template void LoadSelected(const std::unordered_set<std::string_view> &, const std::string &,
                           std::unordered_map<std::string, std::string> *, std::unordered_map<std::string, std::string> *);

void RebuildDtabaseMaxQuant(const std::set<std::string> &protein_names, const std::unordered_map<std::string, std::string> &name_seq, const std::unordered_map<std::string, std::string> &name_full_name, const std::string &out_path) {
    FILE * refined_file = OpenFileWrite(out_path);
    setvbuf ( refined_file , NULL , _IOFBF , 1024000 );
    int32_t con_cnt = 0;
    for (auto protein_name = protein_names.begin(), end_ = protein_names.end(); protein_name != end_; protein_name++) {
        if (name_seq.find(*protein_name) == name_seq.end())
            std::cout << "The Protein " << *protein_name << " can not be found, please check the format!\n";
        if (name_full_name.find(*protein_name) == name_full_name.end())
            std::cout << "The full name of protein " << *protein_name << " can not be found, please check the format!\n";
        std::string fullname = name_full_name.at(*protein_name);
        if (protein_name->substr(0, 3) == "CON") {
            fwrite(">CON|", 5, 1, refined_file);
            fwrite(std::to_string(con_cnt).c_str(), std::to_string(con_cnt).length(), 1, refined_file);
            con_cnt++;
        }
        else
            fwrite(name_full_name.at(*protein_name).c_str(), name_full_name.at(*protein_name).length(), 1, refined_file);
        fwrite("\n", 1, 1, refined_file);
        fwrite((name_seq.at(*protein_name)).c_str(), (name_seq.at(*protein_name)).length(), 1, refined_file);
        fwrite("\n", 1, 1, refined_file);
    }
    fclose(refined_file);
}

void RebuildDtabase(const std::set<std::string> &protein_names, const std::unordered_map<std::string, std::string> &name_seq, const std::unordered_map<std::string, std::string> &name_full_name, const std::string &out_path) {
    FILE * refined_file = OpenFileWrite(out_path);
    setvbuf ( refined_file , NULL , _IOFBF , 1024000 );
    for (auto protein_name = protein_names.begin(), end_ = protein_names.end(); protein_name != end_; protein_name++) {
        if (name_seq.find(*protein_name) == name_seq.end())
            std::cout << "The Protein " << *protein_name << " can not be found, please check the format!\n";
        if (name_full_name.find(*protein_name) == name_full_name.end())
            std::cout << "The full name of protein " << *protein_name << " can not be found, please check the format!\n";
        fwrite((name_full_name.at(*protein_name)).c_str(), (name_full_name.at(*protein_name)).length(), 1, refined_file);
        fwrite("\n", 1, 1, refined_file);
        fwrite((name_seq.at(*protein_name)).c_str(), (name_seq.at(*protein_name)).length(), 1, refined_file);
        fwrite("\n", 1, 1, refined_file);
    }
    fclose(refined_file);
}

void RenameCON(std::set<std::string> *protein_names) {
    // 删除污染蛋白
    for (auto protein_name = protein_names->begin(), end_ = protein_names->end(); protein_name != end_; protein_name++) {
        if ("CON" == protein_name->substr(0, 3))
            protein_names->erase(protein_name);
    }
}

void DBReducer(const std::string &protein_path, const std::string &fasta_path, const std::string &out_path, float threshold, bool maxquant_flag) {
    std::set<std::string> protein_names;
    std::unordered_map<std::string, std::string> name_seq;
    std::unordered_map<std::string, std::string> name_full_name;
    ReadProteinName(&protein_names, protein_path, threshold); // 读取被鉴定到的蛋白质
    // D9：原先先 BuildDict 把整个原始库读进内存，现在只读入被鉴定到的蛋白（输出按名字排序，与库的顺序无关）
    LoadSelected(protein_names, fasta_path, &name_seq, &name_full_name);
    dbr::StatusPhase("write");
    if (maxquant_flag) // 如果下游是maxquant，则需要重命名CON
        RebuildDtabaseMaxQuant(protein_names, name_seq, name_full_name, out_path);
    else
        RebuildDtabase(protein_names, name_seq, name_full_name, out_path); // 数据库重建
}