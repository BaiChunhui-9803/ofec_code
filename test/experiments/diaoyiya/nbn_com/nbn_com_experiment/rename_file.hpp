
#include "../core/global.h"
#include "../core/exception.h"
#include <fstream>
#include <list>
#include <iostream>
#include <chrono>
#include <thread>
//#include <nbnFigureOut.h>
#include <iostream>
//#include "testOFECmatlab.h"
//#include "drawMeshFun.h"



#ifdef  OFEC_PYTHON
#include "../utility/python/python_caller.h"
#include "Python.h"
#endif //  OFEC_PYTHON



#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../utility/nbn_visualization/nbn_onemax_simple.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"

#include "../instance/problem/combination/wmodel/wmodel.h"
#include "../core/problem/continuous/continuous.h"


#include "../utility/general_multithread/general_multithread.h"
#include "../utility/ioh_experimenter/include/ioh/problem/wmodel/wmodel_one_max.hpp"
#include "../utility/ioh_experimenter/include/ioh/problem/wmodel/wmodel_leading_ones.hpp"







#include <filesystem>
#include <string>
#include <vector>
#include <iostream>


void readRemoveDir() {

	std::string path = "F:/diaoyiya/fla_nbn_combinatorial/one_sub_structure/";
	for (const auto& entry : fs::directory_iterator(path)) {
		std::cout << entry.path().filename() << std::endl;
		auto filename = entry.path().filename().string();
		int from = filename.find_last_of("-");
		if (from == -1) {
			fs::remove(entry.path());
			continue;
		}
		auto to = filename.find_last_of(".");
		auto prefix = filename.substr(from + 1, to - from - 1);

		if (prefix == "top" || prefix == "front") {
			//
			//std::cout << "prefix\t" << prefix << std::endl;
			fs::remove(entry.path());
		}
	}
}


void readDirRename() {


	std::vector<int> nameIds = { 0,1,2,7,8,9,10,11,12,15,16 };

	std::string path = "F:/diaoyiya/fla_nbn_combinatorial/one_sub_structure/";
	for (const auto& entry : fs::directory_iterator(path)) {
		std::cout << entry.path().filename() << std::endl;
		auto filename = entry.path().filename().string();

		auto filenameSet = UTILITY::split(filename, "_");

		std::string newName;
		for (auto& it : nameIds) {
			newName += filenameSet[it] + "_";
		}
		newName += ".png";
		fs::rename(entry.path(), path + newName);
		//		int stop = -1;
	}
}


namespace ofec {

	void registerParamAbbr();
	void customizeFileName();
	void run() {
		registerParamAbbr();
		customizeFileName();
		using namespace std;


		readDirRename();

	}

	void registerParamAbbr() {}

	void customizeFileName() {}
}
