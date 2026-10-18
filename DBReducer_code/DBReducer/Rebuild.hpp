/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2022-06-27 13:24:59
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-06-27 16:17:11
 */
#ifndef Rebuild_NAME
#define Rebuild_NAME

#include <string>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string_view>

void DBReducer(const std::string &protein_path, const std::string &fasta_path, const std::string &out_path, float threshold, bool maxquant_flag);
void BuildDict(std::unordered_map<std::string, std::string> *name_seq, const std::string &fasta_path, std::unordered_map<std::string, std::string> *name_full_name);
template <class NameSet>  // 已实例化：std::set<std::string>、std::unordered_set<std::string_view>
void LoadSelected(const NameSet &wanted, const std::string &fasta_path,
                  std::unordered_map<std::string, std::string> *name_seq, std::unordered_map<std::string, std::string> *name_full_name);

#endif