/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2022-04-06 17:58:06
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-06-27 14:30:27
 */
#ifndef StringUtil_NAME
#define StringUtil_NAME


#include <vector>
#include <string>
#include <algorithm>

bool isFloat(const std::string &myString);
std::vector<std::string> Split(const std::string& s, const std::string & delimiters = " ");
std::vector<std::string> SplitByChars(const std::string& input, const std::string& delimiters = " ");
std::string& Strip(std::string *s, const std::string& chars = " \n\t\r");
std::vector<std::string> SplitRule(const std::string& s, const std::string & delimiter);
inline void capitalizeString(std::string *s) {
    std::transform(s->begin(), s->end(), s->begin(),
                   [](unsigned char c){ return std::toupper(c); });
}

#endif