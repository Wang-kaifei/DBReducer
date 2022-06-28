/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2022-04-06 17:58:06
 * @LastEditors: sueRimn
 * @LastEditTime: 2022-06-27 14:29:45
 */
#include "StringUtil.h"
#include <iostream>
#include <sstream>

bool isFloat(const std::string &myString) {
    std::istringstream iss(myString);
    float f;
    iss >> f;
    return iss.eof() && !iss.fail(); // Check the entire string was consumed and if either failbit or badbit is set
}

std::vector<std::string> Split(const std::string& s, const std::string & delimiter) {
	std::vector<std::string> res;
	auto pos = s.find(delimiter, 0);
	size_t start = 0;
	while(std::string::npos != pos) {
		if (pos > start)
			res.emplace_back(s.substr(start, pos - start));
		start = pos + delimiter.length();
		pos = s.find(delimiter, start);
	}
	if (start < s.length()) 
		res.emplace_back(s.substr(start));
	return res;
}

/**
 * @brief Split,多个分隔符也会存空字符串
 * 
 * @param s 
 * @param delimiter 
 * @return std::vector<std::string> 
 */
std::vector<std::string> SplitRule(const std::string& s, const std::string & delimiter) {
	std::vector<std::string> res;
	auto pos = s.find(delimiter, 0);
	size_t start = 0;
	while(std::string::npos != pos) {
		res.emplace_back(s.substr(start, pos - start));
		start = pos + delimiter.length();
		pos = s.find(delimiter, start);
	}
	res.emplace_back(s.substr(start));
	return res;
}

std::vector<std::string> SplitByChars(const std::string& s, const std::string & delimiters) {
	std::vector<std::string> res;
	auto lastPos = s.find_first_not_of(delimiters, 0); // 第一个不在delimiters中的字符
	auto pos = s.find_first_of(delimiters, lastPos); // lastpos之后的第一个在delimiters中的字符
	while (std::string::npos != pos || std::string::npos != lastPos) {
		res.emplace_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delimiters, pos);
		pos = s.find_first_of(delimiters, lastPos);
    }
	return res;
}

std::string& Strip(std::string *s, const std::string& chars) {
	s->erase(0, s->find_first_not_of(chars.c_str()));
	s->erase(s->find_last_not_of(chars.c_str()) + 1);
	return *s;
}