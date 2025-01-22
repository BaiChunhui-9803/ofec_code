
#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <limits>

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "../utility/hnsw/hnsw_solution_multithread/hnsw_sol_model.h"
#include "../core/global.h"
#include "../core/problem/optima.h"
#include "../core/definition.h"
#include "../core/global.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../utility/nbn_visualization/nbn_grid_tree_division.h"
#include "../utility/nbn_visualization/nbn_grid_tree_division_allSize.h"
#include  "../utility/function/custom_function.h"
#include "../utility/function/general_function.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"
//#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "interface.h"
//#include "../utility/function/custom_function.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../core/algorithm/population.h"
//#include "../utility/nondominated_sorting/filter_sort.h"
#include "../instance/algorithm/visualize/sampling/instance/sampling_eax_tsp.h"
#include "../instance/algorithm/visualize/sampling/sampling_multiThread.h"
#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include <queue>
#include "../core/definition.h"
#include <deque>
#include <chrono>   
//#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator_multiThead.h"

#include "../utility/hnsw/hnsw_nbn/neighbor_info.h"
#include "../utility/hnsw/hnsw_nbn/select_policy.h"
#include "../utility/hnsw/hnsw_nbn/node.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_model.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_nbn.h"
#include "../core/environment/environment.h"
#include "../core/parameter/variants_to_stream.h"

//#include "../instance/algorithm/combination/lkh/lkh_algorithm.h"
#include<filesystem>
#include<queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_test/eax_tsp_alg.h"
#include "../instance/problem/combination/travelling_salesman/tsp_offline_data/travelling_salesman_offline_data.h"

#include "../utility/nbn_visualization/nbn_fla/nbn_fla.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_basin/eax_tsp_basin_multipop.h"
#include "../utility/nbn_visualization/nbn_fla/nbn_fla_utility.h"
#include "../instance/algorithm/combination/LKH_origin/INCLUDE/LKH.h"


#include <iostream>
#include <iomanip>
#include <ctime>
#include <sstream>


#include "../utility/visualization/idee_trait.h"


int getSolSize(const std::string& fileatt) {
	int inputId;
	std::string line;

	int solSize(0);
	std::ifstream in(fileatt);
	while (in >> inputId) {
		std::getline(in, line);
		//sols.push_back(line);

		++solSize;
	}

	in.close();

	return solSize;
}




void runEXE(const std::string& cmd) {
	std::system(cmd.c_str());
}
std::string LDEE_path;

void neighborSols(const std::string& path, const std::string& filename) {
	int solSize = getSolSize(path + "/" + filename + "_solutions.txt");
	int num_task = std::thread::hardware_concurrency();
	int batchNum = std::ceil(double(solSize) / double(num_task));
	//std::string LDEE_path = "E:/Diao_Yiya/paper/fitness_landscape/IDEE/ldee_tools/";

	std::vector<std::string> cmds;

	for (int idx(0); idx < num_task; ++idx) {
		std::string cmd = LDEE_path + "NeighbourProb.exe --problem-type=perm --batch-size=" + std::to_string(batchNum) + std::string(" --batch-number=") + std::to_string(idx) + " " + path + "/" + filename + "_solutions.txt";

		//	std::string tmp_path = LDEE_path + "NeighbourProb.exe" + " E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data/" + filename + "_solutions.txt";
		//	std::cout << tmp_path << std::endl;
		//	std::system(tmp_path.c_str());
		//	std::cout << cmd << std::endl;
		//	std::system(cmd.c_str());


		cmds.push_back(cmd);
	}

	std::vector<std::thread> thrds;
	for (int idx(0); idx < num_task; ++idx) {
		//std::thread t3(f2, std::ref(n));

		thrds.push_back(std::thread(
			runEXE, std::cref(cmds[idx])));
	}

	for (auto& thrd : thrds)
		thrd.join();


	//" CombineFiles.exe moffp_nearest_neighbour_probabilities.txt 255";
	std::string cmd = LDEE_path + "CombineFiles.exe " + filename + "_nearest_neighbour_probabilities.txt " + std::to_string(num_task - 1);

	std::system(cmd.c_str());


	//std::string cmd = 100 0 moffp_solutions.txt";


}
int getNumDataLine(const std::string& path) {
	std::ifstream in(path);
	std::string inputStr;
	int inputNum;
	double inputDouble;
	int maxNum(0);
	//in >> inputStr >> inputStr >> inputStr >> inputStr;
	while (in) {
		in >> inputNum;
		maxNum = std::max(maxNum, inputNum);
		in >> inputDouble >> inputDouble;
	}
	in.close();

	return maxNum;
}


void runLDEE(const std::string& path, const std::string& filename) {
	std::vector<std::string> subfix = { "_solutions.txt",
		"_nearest_neighbour_probabilities.txt.combined",
		"_.txt_embedded_tsne.txt",
		"_.txt_embedded_ve.txt",
		"_attributes.txt" };


	std::string cmd;


	cmd = LDEE_path + "NeighbourProb.exe --problem-type=perm " + path + "/" + filename + subfix[0];
	std::system(cmd.c_str());

	neighborSols(path, filename);


	return;
	if (!std::filesystem::exists(path + "/" + filename + subfix[1])) {
		return;
	}


	int num_task = std::thread::hardware_concurrency();
	// rename file
//	cmd = LDEE_path + "CombineFiles.exe " + filename + "_nearest_neighbour_probabilities.txt " + std::to_string(num_task - 1);

//	std::system(cmd.c_str());

//	std::filesystem::copy(filename + "_nearest_neighbour_probabilities.txt.combined", filename + subfix[2]);



//	cmd = LDEE_path + "TSNEEmbedding.exe " + "E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data_nearest_neighbor/" + filename + subfix[1];


	if (std::filesystem::exists(path + "/" + filename + subfix[2])) {
	}
	else {

		cmd = LDEE_path + "TSNEEmbedding.exe " + path + "/" + filename + subfix[1];
		std::cout << "cmd\t" << cmd << std::endl;
		std::system(cmd.c_str());
	}

	if (!std::filesystem::exists(path + "/" + filename + subfix[2])) {
		return;
	}



	if (std::filesystem::exists(path + "/" + filename + subfix[3])) {
		return;
	}

	//int numData = getNumDataLine(filename + subfix[2]);

	//int runningData = int(std::sqrt(double(numData)));
	//int skipData = numData - runningData * runningData;
	cmd = LDEE_path + "VacuumEmbedding.exe " + path + "/" + filename + subfix[2];


	//cmd = LDEE_path + "VacuumEmbedding.exe --skip-rows=" + std::to_string(skipData) + " " + filename + subfix[2];
	std::system(cmd.c_str());


}



void runIDEE(const std::string& filename) {

	auto tspname = filename.substr(0, filename.find_first_of("."));
	tspname = tspname + "IDEE";

	auto tool_path = "F:/code/ldee_tools/";
	tool_path = "E:/DiaoYiya/code/ldee_tools/";
	auto data_dir = "F:/code/ofec_data/idee_data4";
	data_dir = "E:/DiaoYiya/code/data/ofec-data/paper_com_experiment_data/totalTsp/idee_eax_lon_data";
	LDEE_path = tool_path;
	runLDEE(data_dir, tspname);



}




void runTask() {

	using namespace ofec;
	using namespace std;
	//ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
	////		ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	////		ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
	//	//ofec::g_working_directory = "E:/DiaoYiya/code/data/ofec-data/";
	////		ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/NBN_datacode_total/data";

	//std::string saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_remote/";
	//saveDir = ofec::g_working_directory + "/paper_com_experiment_data/tsp_typical_comparison_exp3/";
	//saveDir = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp/idee_datas/";

	//std::filesystem::create_directories(saveDir);


	//std::string dir1 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
	//std::string dir2 = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman_origin/";

	std::vector<std::string> filenames =
	{ "2202.tsp","u574.tsp" ,  "5955.tsp","1281.tsp" ,"6702.tsp" ,"6717.tsp" ,  "7310.tsp", "9225.tsp", };

	
	filenames =
	{ "5955.tsp", "u574.tsp" ,  "1281.tsp" };

	
	for (auto& it : filenames) {
		runIDEE(it);
	}


	//	std::reverse(filenames.begin(), filenames.end());

	//for (auto& it : filenames) {
	//	if (it == "u574.tsp") {
	//		calTask(dir2, saveDir, it);
	//	}
	//	else {
	//		calTask(dir1, saveDir, it);
	//	}
	//}

}



namespace ofec {

	void registerParamAbbr() {}
	void customizeFileName() {}
	void run() {

		using namespace ofec;
		using namespace std;

		registerInstance();


		runTask();


		cin.get();

	}


}