#ifndef _TOOLS_HPP_
#define _TOOLS_HPP_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>



extern long long hashing_value(double x, int digits);

extern std::string toHashStr(const std::vector<double>& x,int digit);
extern std::string toHashStrFit(double fit, int digit);

extern bool loadConfigurationFile(const std::string name, std::vector<std::string>& func,
    std::vector<int>& nvar, std::vector<int>& bounded,
    std::vector<double>& best,
    std::string& step_mode, std::vector<double>& beta,
    std::vector<double>& factor, int& opt_digits,
    int& hash_digits, int& iter, int& runs,
    std::string& init_mode, std::string& output_dir);
extern bool displayConfiguration(const std::vector<std::string>& func,
    const std::vector<int>& nvar, const std::vector<int>& bounded,
    const std::vector<double>& best,
    const std::string& step_mode, const std::vector<double>& beta,
    const std::vector<double>& factor, const int& opt_digits,
    const int& hash_digits, const int& iter, const int& runs,
    const std::string& init_mode, const std::string& output_dir);
#endif
