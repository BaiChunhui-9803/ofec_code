#ifndef LON_H
#define LON_H
#include <string>
#include <vector>
#include "../sampling/lon_con_sampling.h"

namespace ofec::lon {


	extern int Convergence(const std::vector<std::string>& argv);
	extern int FilterData(const std::vector<std::string>& argv);

	extern int GenerateResults(const std::string& dir_path, const std::vector<std::string>& argv);
	//extern int GenerateResults(const std::string& dir_path,const std::vector<std::string>& argv);


	extern int GenerateResultsEigenSampling(const std::string& rpath, const std::string& dir_path, const ParameterMap& v);




	extern int GenerateResultsEigenSamplingNode(
		const std::string& dir_path, 
		 LonSamplingPar& v,
		LonStruct &lon);
	extern void EvaluateLONbyR(const std::string& rpath, const std::string& dir_path, const LonSamplingPar& v);
	

	extern int Hashing(const std::vector<std::string>& argv);



}

#endif