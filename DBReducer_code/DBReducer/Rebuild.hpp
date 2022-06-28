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

void DBReducer(const std::string &protein_path, const std::string &fasta_path, const std::string &out_path, float threshold, bool maxquant_flag);

#endif