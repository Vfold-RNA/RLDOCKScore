#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <cctype> // for isspace
// #include <boost/lexical_cast.hpp>
#include <iomanip>//std::stew

namespace rlscore {

inline std::vector<std::string> string2vector(const std::string sline) {
	std::istringstream ss(sline);
	std::string buf;
	std::vector<std::string> token;
	while(ss >> buf) token.push_back(buf);
	return token;
}

// inline bool starts_with(const std::string& str, const std::string& start) {
// 	return str.size() >= start.size() && str.substr(0, start.size()) == start;
// }

inline std::string trim(const std::string &str, const std::string &whitespace = " \t\n\r") {
	const auto str_begin = str.find_first_not_of(whitespace);
	if (str_begin == std::string::npos){
		return ""; // no content
	}
	const auto str_end = str.find_last_not_of(whitespace);
	const auto str_range = str_end - str_begin + 1;
	return str.substr(str_begin, str_range);
}

// inline std::string left_trim(const std::string &str, const std::string &whitespace = " \t") {
// 	const auto str_begin = str.find_first_not_of(whitespace);
// 	if (str_begin == std::string::npos){
// 		return ""; // no content
// 	}
// 	const auto str_end = str.size();
// 	const auto str_range = str_end - str_begin;
// 	return str.substr(str_begin, str_range);
// }

// inline std::string right_trim(const std::string &str, const std::string &whitespace = " \t") {
// 	std::string::int str_begin = 0;
// 	const auto str_end = str.find_last_not_of(whitespace);
// 	if (str_end == std::string::npos){
// 		return ""; // no content
// 	}
// 	const auto str_range = str_end - str_begin + 1;
// 	return str.substr(str_begin, str_range);
// }

// inline std::string reduce(const std::string& str, const std::string& fill = " ", const std::string& whitespace = " \t") {
// 	// Trim first
// 	auto result = trim(str, whitespace);
// 	// replace sub ranges
// 	auto beginSpace = result.find_first_of(whitespace);
// 	while (beginSpace != std::string::npos){
// 		const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
// 		const auto range = endSpace - beginSpace;
// 		result.replace(beginSpace, range, fill);
// 		const auto newStart = beginSpace + fill.length();
// 		beginSpace = result.find_first_of(whitespace, newStart);
// 	}
// 	return result;
// }


// template<typename T>
// std::string to_string(const T& x, std::streamsize width = 0, char fill = ' ') { // default 0 width means no restrictions on width
// 	std::ostringstream out;
// 	out.fill(fill);
// 	if(width > 0)
// 		out << std::setw(width);
// 	out << x;
// 	return out.str();
// }

// struct bad_conversion {};

// template<typename T>
// T convert_substring(const std::string& str, int i, int j) { // indexes are 1-based, the substring should be non-null
// 	if(i < 1 || i > j+1 || j > str.size()) throw bad_conversion();

// 	// omit leading whitespace
// 	while(i <= j && std::isspace(str[i-1]))
// 		++i;

// 	T tmp;
// 	try {
// 		tmp = boost::lexical_cast<T>(str.substr(i-1, j-i+1));
// 	}
// 	catch(...) {
// 		throw bad_conversion();
// 	}
// 	return tmp;
// }

// inline bool substring_is_blank(const std::string& str, int i, int j) { // indexes are 1-based, the substring should be non-null
// 	if(i < 1 || i > j+1 || j > str.size()) throw bad_conversion();
// 	for(int k = i-1; k < j; ++k) {
// 		if(!std::isspace(str[k]))
// 			return false;
// 	}
// 	return true;
// }

// // when this was written, lexical cast to unsigned didn't work very well with "-123", etc.
// template<>
// inline unsigned convert_substring<unsigned>(const std::string& str, int i, int j) { // indexes are 1-based, the substring should be non-null
// 	int tmp = convert_substring<int>(str, i, j);
// 	if(tmp < 0) throw bad_conversion();
// 	return static_cast<unsigned>(tmp);
// }

}