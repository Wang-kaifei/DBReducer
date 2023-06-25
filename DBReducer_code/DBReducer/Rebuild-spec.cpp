/*
 * @Descripttion: 从pFind-Filtered.spectra文件去建库
 * @version: 
 * @Author: sueRimn
 * @Date: 2022-06-27 13:24:46
 * @LastEditors: Kaifei
 * @LastEditTime: 2023-05-31 16:48:37
 */

#include "Param.hpp"
#include <typeinfo>
#include <omp.h>
#include <iostream>
#include <set>
#include <algorithm>
#include <exception>
#include <unordered_map>
#include <math.h>
#include "../utils/IOUtil.h"
#include "../utils/StringUtil.h"
#include "Rebuild-spec.hpp"
#include "Rebuild.hpp"
#include <functional>


#define VNAME(value) (#value)

namespace ReduSpec {

    double m_lfCode[256][PRIME_SIZE];

    bool GodelInitialize() // 编码打表
    {
        for(int i = 0;i < 256;++i) // 遍历的应该是ASCII编码，避免aa - 'A'的操作
        {
            for(int j = 0;j < PRIME_SIZE;++j)
            {
                m_lfCode[i][j] = (i + 1) * log(prime[j]);
            }
        }
        return true;
    }

    inline bool CMPPepNum(const Protein &a, const Protein &b) {
        return a.godel_codes.size() > b.godel_codes.size();
    }

    /**
     * @brief Get godel code for peptide
     * 
     * @param peptide aa sequence
     * @return double godel code
     */
    double GetGodel(const std::string &peptide) {
        double lfCode = 0.0;
        lfCode -= lfCode;
        for(size_t i = 0; i < peptide.length(); ++i)
        {
            lfCode += m_lfCode[aacode[peptide[i] - 'A']][i];
            // lfCode += m_lfCode[(unsigned char)peptide[i]][i];
        }
        return lfCode;
    }

/**
 * @brief 并行处理读取的segs，存到pros中
 * 
 * @param lines spectra文件中的每一行原始结果
 * @param pros 蛋白质容器
 */
void StoreProbyLine(const std::deque<RawRes> &lines, std::deque<Protein> *pros) {
    omp_lock_t pro_lock; // 存储蛋白质的锁
    omp_init_lock(&pro_lock);
    std::hash<std::string> hasher; // 散列器
    int i = 0;
    std::unordered_map<size_t, std::vector<Protein>> hash_pro_table; // key = hash code of name, value = the pros have this key hash name
    // #pragma omp parallel for num_threads(20)
    for (auto line = lines.begin(); line != lines.end(); line++) {
        auto names = SplitRule(line->pros, "/"); // 肽段对应到的蛋白质
        double code = GetGodel(line->pep); // 肽段编码
        for (auto name = names.begin(); name != names.end(); name++) { // 遍历所有蛋白名
            if (name->size() == 0 || (name->size() >= 3 && name->substr(0, 3) == "REV")) {
                continue;
            }
            Protein tmp(*name);
            std::size_t hash_code = hasher(*name);
            omp_set_lock(&pro_lock);
            auto stored_pros = hash_pro_table.find(hash_code);
            if (stored_pros != hash_pro_table.end()) { // 如果该hash code 被存储过
                auto iter = std::find(stored_pros->second.begin(), stored_pros->second.end(), tmp);
                if (iter != stored_pros->second.end()) { // 该name已被存储
                    iter->godel_codes.insert(code);
                }
                else { // 该name未被存储（即hash冲突）
                    std::cout << "Hash collision!\n";
                    tmp.godel_codes.insert(code);
                    stored_pros->second.emplace_back(tmp);
                }
            }
            else { // 如果该hash code未被存储
                tmp.godel_codes.insert(code);
                hash_pro_table.insert(std::make_pair(hash_code, std::vector<Protein>{tmp}));
            }
            omp_unset_lock(&pro_lock);
        }
    }
    omp_destroy_lock (&pro_lock);
    for (auto iter = hash_pro_table.begin(); iter != hash_pro_table.end(); iter++) {
        for (auto pro = iter->second.begin(); pro != iter->second.end(); pro++) {
            pros->emplace_back(*pro);
        }
    }
}    

    /**
     * @brief 从.spectra读取PSM并存成RawRes形式
     * 
     * @param file_name .spectra文件路径
     * @param lines 存储RawRes的容器
     * @param threshold q-value阈值
     */
    void ReadRawResbySpec(const std::string &file_name, std::deque<RawRes> *lines, double threshold) {
        std::ifstream buffer = ReadFilebyStream(file_name);
        std::string linestr= "";
        getline(buffer, linestr); // 过滤掉title行
        while (!buffer.eof()) { 
            getline(buffer, linestr);// 按行读取缓冲区内容
            auto segs = SplitRule(linestr, "\t");
            if (std::stof(segs[4]) > threshold)
                break;
            lines->emplace_back(RawRes(segs[12], segs[5]));
        }
    }

    /**
     * @brief 判断test下标对应的protein是否为lead所对应的subset
     * 
     * @param pros 存储蛋白质的容器（按照肽段数从多到少排序）
     * @param lead 
     * @param test 
     * @return true 
     * @return false 
     */
    bool IsContain(const std::deque<Protein> &pros, int lead, int test) {
        for (auto code = pros[test].godel_codes.begin(); code != pros[test].godel_codes.end(); code++) {
            if (pros[lead].godel_codes.find(*code) == pros[lead].godel_codes.end()) // 如果有没有被包含的code
                return false;
        }
        return true;
    }

    /**
     * @brief 判断第index个蛋白质是否为subset
     * 
     * @param pros 存储蛋白质的容器（按照肽段数从多到少排序）
     * @param lead_pros 存储已经被判定为非subset的蛋白索引
     * @param index 
     * @return true 被判定为subset
     * @return false 
     */
    bool IsSubset(const std::deque<Protein> &pros, const std::vector<int> &lead_pros, int index) {
        for (auto iter = lead_pros.begin(); iter != lead_pros.end(); iter++) { // 遍历可能为lead的index
            if (pros[*iter].godel_codes.size() == pros[index].godel_codes.size()) {
                break;
            }
            if (IsContain(pros, *iter, index))
                return true;
        }
        return false;
    }

    /**
     * @brief 给蛋白质分组
     * 
     * @param pros 存储蛋白质的容器
     */
    void GroupPros(std::deque<Protein> *pros) {
        std::sort(pros->begin(), pros->end(), CMPPepNum); // 按照包含肽段数目从高到低排序
        std::vector<int> lead_pros; // 存储非subset protein的index 
        for (int i = 0; i < pros->size(); i++) {
            if (IsSubset(*pros, lead_pros, i)) {
                pros->at(i).is_subset = true;
            }
            else {
                lead_pros.emplace_back(i);
            }
        }
    }

    /**
     * @brief 写出非subset
     * 
     * @param input_path 原始fasta文件
     * @param output_path 输出新的refined fasta文件
     * @param pros 鉴定到的蛋白容器
     * @param name_seq key = pro name, value = pro seq
     */
    void ReduceDB(const std::string &output_path, const std::deque<Protein> &pros, const std::unordered_map<std::string, std::string> &name_seq, const std::unordered_map<std::string, std::string> &name_full_name) {
        FILE *out_file = OpenFileWrite(output_path); // 打开输出文件
        setvbuf(out_file, NULL, _IOFBF, g_out_buffer_size); // 设置缓冲区
        for (auto pro = pros.begin(); pro != pros.end(); pro++) {
            if (!pro->is_subset) {
                // fwrite(">", 1, 1, out_file);
                if (name_full_name.find(pro->name) == name_full_name.end()) {
                    std::cout << "Not found name: " << pro->name << "  &&"<< std::endl;
                    continue;
                }
                if (name_seq.find(pro->name) == name_seq.end()) {
                    std::cout << "Not found sequence: " << pro->name << "  &&" << std::endl;
                    continue;
                }
                fwrite((name_full_name.at(pro->name)).c_str(), (name_full_name.at(pro->name)).length(), 1, out_file);
                fwrite("\n", 1, 1, out_file);
                fwrite(name_seq.at(pro->name).c_str(), name_seq.at(pro->name).size(), 1, out_file);
                fwrite("\n", 1, 1, out_file);
            }
        }
        fclose(out_file); // 关闭文件
    }
}