#ifndef UTILS_PREPROCESS_H
#define UTILS_PREPROCESS_H

#include "string"

bool pick_option(std::vector<std::string>& args, const std::string& option);
std::string pick_option(std::vector<std::string>& args, const std::string& option, const std::string& default_value);
Parameters init_params(const std::string& file_params, int step_alg);
std::string path_abs2rel(std::string& option);
#endif // UTILS_PREPROCESS_H
