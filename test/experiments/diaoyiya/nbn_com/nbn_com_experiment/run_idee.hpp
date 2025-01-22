//
//#include "../core/global.h"
//#include <iostream>
//#include <filesystem>
//#include <memory>
//
//
//#include "custom_method.h"
//#include "prime_method.h"
//#include "../core/global.h"
//
#include <fstream>
#include <list>
#include <iostream>
#include <chrono>
#include <thread>



#include <iostream>
#include <filesystem>
#include <string>
#include<set>
#include <mutex>
#include <thread>


void runAlg() {
	//system("cd E:/Diao_Yiya/paper/fitness_landscape/IDEE/result");
	system("E:/Diao_Yiya/paper/fitness_landscape/IDEE/ldee_tools/NeighbourProb.exe --problem-type=perm  E:/Diao_Yiya/paper/fitness_landscape/IDEE/data/tspPop_1_solutions.txt");

}

void readFile(const std::string& path, const std::string& targetPath) {
	std::vector<std::string> filepaths;
	const std::filesystem::path sandbox{ path };
	std::string filename;
	std::string dirname;

	std::string filefix1 = "_attributes.txt";
	std::string filefix2 = "_solutions.txt";

	// directory_iterator can be iterated using a range-for loop
	for (auto const& dir_entry : std::filesystem::directory_iterator{ sandbox })
	{
		//std::cout << "dir_entry\t" << dir_entry.path() << std::endl;
		//const std::filesystem::path subdir{ dir_entry.path() };
		const std::filesystem::path innerPath{ dir_entry.path() };

		auto dirPath = dir_entry.path().string();
		dirname = dir_entry.path().string();
		dirname = dirname.substr(dirname.find_last_of("\\") + 1, dirname.size());

		//bool firstFlag(false);
		for (auto const& dir_entry2 : std::filesystem::directory_iterator{ innerPath })
		{
			filename = dir_entry2.path().string();
			break;
			//	std::cout << "dir_entry\t" << dir_entry2.path() << std::endl;
		}
		filename = filename.substr(0, filename.find_last_of("_"));
		filename = filename.substr(filename.find_last_of("\\") + 1, filename.size());


		std::filesystem::copy(dirPath + "/" + filename + filefix1, targetPath + "/" + dirname + "_" + filename + filefix1);
		std::filesystem::copy(dirPath + "/" + filename + filefix2, targetPath + "/" + dirname + "_" + filename + filefix2);


		//	int stop = -1;
	}
}

#include <fstream>

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
#include<array>

template<typename T>
inline T mapReal(T value, T input_min, T input_max, T output_min, T output_max) {
	return ((value - input_min) / (input_max - input_min) * (output_max - output_min) + output_min);
}



std::vector<std::string> split(std::string s, std::string delimiter) {
	size_t pos_start = 0, pos_end, delim_len = delimiter.length();
	std::string token;
	std::vector<std::string> res;

	while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
		token = s.substr(pos_start, pos_end - pos_start);
		pos_start = pos_end + delim_len;
		res.push_back(token);
	}

	res.push_back(s.substr(pos_start));
	return res;
}



void filterData(const std::string& ve_file,
	const std::string& origin_file,
	const std::string& newFile, std::vector<std::string>& head) {

	int size = 5;

	std::vector<int> ve_data;


	{
		int lineId(0), a, b;
		std::ifstream in(ve_file);
		while (in >> lineId >> a >> b) {
			ve_data.push_back(lineId);
		}
		in.close();
	}
	if (ve_data.empty())return;

	//std::array<std::string, 4> head;
	//std::vector<std::vector<int>> infos;
	std::vector<std::string> infos;
	std::string curline;
	int lineId = 0;
	std::string headline;
	infos.resize(ve_data.back() + 1);
	{

		std::ifstream in(origin_file);
		std::getline(in, headline);
		//in >> head[0] >> head[1] >> head[2] >> head[3];
		//std::vector<int> info(size);
		while (in >> lineId) {
			std::getline(in, curline);
			if (infos.size() <= lineId) {
				infos.resize(lineId + 1);
			}
			infos[lineId] = curline;
		}
		in.close();
	}

	//// normalize data
	//std::vector<std::array<double, 4>> filter_data;
	//std::array<double, 4> cur_data;
	//{
	//	for (auto& it : ve_data) {
	//		auto& curInt = infos[it];
	//		for (int idx(0); idx < 4; ++idx) {
	//			cur_data[idx] = curInt[idx];
	//		}
	//		filter_data.push_back(cur_data);
	//	}

	//	std::array<double, 4> min_val;
	//	std::fill(min_val.begin(), min_val.end(), std::numeric_limits<double>::max());
	//	std::array<double, 4> max_val;
	//	std::fill(max_val.begin(), max_val.end(), std::numeric_limits<double>::min());
	//	
	//	for (auto& it : filter_data) {

	//		for (int idx(1); idx < 4; ++idx) {
	//			min_val[idx] = std::min(min_val[idx], it[idx]);
	//			max_val[idx] = std::max(max_val[idx], it[idx]);
	//		//	it[idx] = mapReal<double>(it[idx],min_va)
	//		}
	//	}

	//	for (auto& it : filter_data) {

	//		for (int idx(1); idx < 4; ++idx) {
	//			it[idx] = mapReal<double>(it[idx], min_val[idx], max_val[idx],0,1);
	//		}
	//	}

	//}


	{
		head = split(headline, "\t");
		head.pop_back();
		std::ofstream out(newFile);
		out << headline << std::endl;
		//for (auto& it : filter_data) {
		//	out << int(it[0]) << "\t" << it[1] << "\t" << it[2] << "\t" << it[3] << std::endl;
		//}

		for (auto& iter : ve_data) {
			auto it = infos[iter];
			out << iter << it << std::endl;

		}
		out.close();
	}




}



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
	//std::string title;
	//int solSize(0);
	//int num;
	//std::ifstream in(fileatt);
	//in >> title >> title >> title >> title;
	//while (in >> num >> num >> num >> num) {
	//	++solSize;
	//}

	//in.close();
	//return solSize;
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


std::string result_path;

void runLDEE(const std::string& path, const std::string& filename) {
	std::vector<std::string> subfix = { "_solutions.txt",
		"_nearest_neighbour_probabilities.txt.combined",
		"_.txt_embedded_tsne.txt",
		"_.txt_embedded_ve.txt",
		"_attributes.txt" };

	std::string LDEE_path = LDEE_path;

	std::string cmd;


	//cmd = LDEE_path + "NeighbourProb.exe --problem-type=perm " + path + "/" + filename + subfix[0];
	//std::system(cmd.c_str());

	neighborSols(path, filename);
	if (!std::filesystem::exists(filename + subfix[1])) {
		return;
	}


	int num_task = std::thread::hardware_concurrency();
	// rename file
//	cmd = LDEE_path + "CombineFiles.exe " + filename + "_nearest_neighbour_probabilities.txt " + std::to_string(num_task - 1);

//	std::system(cmd.c_str());

//	std::filesystem::copy(filename + "_nearest_neighbour_probabilities.txt.combined", filename + subfix[2]);



//	cmd = LDEE_path + "TSNEEmbedding.exe " + "E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data_nearest_neighbor/" + filename + subfix[1];

	cmd = LDEE_path + "TSNEEmbedding.exe " + filename + subfix[1];

	std::system(cmd.c_str());


	if (!std::filesystem::exists(filename + subfix[2])) {
		return;
	}

	int numData = getNumDataLine(filename + subfix[2]);

	int runningData = int(std::sqrt(double(numData)));
	int skipData = numData - runningData * runningData;

	cmd = LDEE_path + "VacuumEmbedding.exe --skip-rows=" + std::to_string(skipData) + " " + filename + subfix[2];
	std::system(cmd.c_str());

	if (!std::filesystem::exists(filename + subfix[3])) {
		return;
	}

	std::vector<std::string> headers;
	filterData(filename + subfix[3], path + "/" + filename + subfix[4], filename + subfix[4], headers);
	if (!std::filesystem::exists(filename + subfix[4])) {
		return;
	}



	cmd = LDEE_path + "HeightMap.exe --height-map-palette=rainbow " +
		filename + subfix[3] + " " + filename + subfix[4] + " 1 " + std::to_string(headers.size() - 1);
	std::system(cmd.c_str());

	//if (!std::filesystem::exists("height_map_1.png")) {
	//	return;
	//}
	//if (std::filesystem::exists(result_path + "/" + filename + "_" + headers[1] + ".png")) {
	//	return;
	//}
	//if (std::filesystem::exists(result_path + "/" + filename + "_" + headers[2] + ".png")) {
	//	return;
	//}
	//if (std::filesystem::exists(result_path + "/" + filename + "_" + headers[3] + ".png")) {
	//	return;
	//}

	std::cout << "copy file" << std::endl;

	for (int idx(1); idx < headers.size(); ++idx) {
		if (std::filesystem::exists(result_path + "/" + filename + "_" + headers[idx] + ".png")) {
			std::filesystem::remove(result_path + "/" + filename + "_" + headers[idx] + ".png");
		}
	}

	for (int idx(1); idx < headers.size(); ++idx) {

		std::string filename = "height_map_" + std::to_string(idx) + ".png";

		std::filesystem::copy(filename, result_path + "/" + filename + "_" + headers[1] + ".png");

	}

	//std::filesystem::copy("height_map_1.png", result_path + "/" + filename + "_" + headers[1] + ".png");
	//std::filesystem::copy("height_map_2.png", result_path + "/" + filename + "_" + headers[2] + ".png");
	//std::filesystem::copy("height_map_3.png", result_path + "/" + filename + "_" + headers[3] + ".png");

	//if (headers.size() > 4) {
	//	std::filesystem::copy("height_map_4.png", result_path + "/" + filename + "_" + headers[4] + ".png");
	//}

	//	int stop = -1;
}


std::vector<std::string> filenames;
std::mutex total_data_mtx;

void getFilenames(const std::string& path) {
	const std::filesystem::path sandbox{ path };
	std::string filename;
	std::set<std::string> setfilenames;
	// directory_iterator can be iterated using a range-for loop
	for (auto const& dir_entry : std::filesystem::directory_iterator{ sandbox })
	{

		filename = dir_entry.path().string();
		filename = filename.substr(0, filename.find_last_of("_"));
		filename = filename.substr(filename.find_last_of("\\") + 1, filename.size());

		setfilenames.insert(filename);
	}

	for (auto& it : setfilenames) {
		filenames.push_back(it);
	}

	//stop = -1;
}


//void getFilenames2()


std::string filedir;

void runThread() {
	bool flagInfo(false);
	//std::string filedir = "E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data";
	std::string filename;
	while (true) {
		flagInfo = false;
		{
			std::lock_guard<std::mutex> guard(total_data_mtx);
			if (!filenames.empty()) {
				filename = filenames.back();
				filenames.pop_back();
				flagInfo = true;
			}
		}
		if (!flagInfo)break;

		//runGA(info.first, info.second);
		//runEAX_GA(curInfo[0], curInfo[1], curInfo[2], double(curInfo[2] + 1) / double(maxRun + 1));
		runLDEE(filedir, filename);

		//{
		//    std::lock_guard<std::mutex> guard(m_out_alg_mtx);
		//    m_out_alg << info.first[1] << "\t" << info.second << std::endl;
		//}
	}
}




void mergeData(const std::string& sourceDir, const std::string& targetDir) {
	std::vector<std::vector<int>> arributes;
	std::vector<std::string> sols;
	std::array<std::string, 5> headers = { "solution_id","run_id","gen_num","objective", "popNum" };
	std::vector<std::string> filepaths;
	const std::filesystem::path sandbox{ sourceDir };
	std::string filename;
	std::string dirname;

	std::string filefix1 = "_attributes.txt";
	std::string filefix2 = "_solutions.txt";

	// directory_iterator can be iterated using a range-for loop
	for (auto const& dir_entry : std::filesystem::directory_iterator{ sandbox })
	{
		//std::cout << "dir_entry\t" << dir_entry.path() << std::endl;
		//const std::filesystem::path subdir{ dir_entry.path() };
		const std::filesystem::path innerPath{ dir_entry.path() };

		auto dirPath = dir_entry.path().string();
		dirname = dir_entry.path().string();
		dirname = dirname.substr(dirname.find_last_of("\\") + 1, dirname.size());

		std::set<std::string> filenames;

		//bool firstFlag(false);
		for (auto const& dir_entry2 : std::filesystem::directory_iterator{ innerPath })
		{
			filename = dir_entry2.path().string();

			filename = filename.substr(0, filename.find_last_of("_"));
			filename = filename.substr(filename.find_last_of("\\") + 1, filename.size());

			filenames.insert(filename);
			//break;
			//	std::cout << "dir_entry\t" << dir_entry2.path() << std::endl;
		}


		std::string fileAtrri = dirname + filefix1;
		std::string fileSols = dirname + filefix2;

		arributes.clear();
		sols.clear();

		int solId = 0;

		std::string title;
		std::vector<int> data(5);
		std::string line;
		int inputId = 0;
		int popId = 0;

		std::string popStr;
		for (auto& filenameIter : filenames) {
			popStr = filenameIter.substr(filenameIter.find_first_of("_") + 1, filenameIter.size());
			popId = std::stoi(popStr);


			//	filename = filename.substr(filename.find_last_of("\\") + 1, filename.size());


			solId = arributes.size() + 1;
			{
				std::ifstream in(dirPath + "/" + filenameIter + filefix1);
				in >> title >> title >> title >> title;
				while (in >> data[0]) {
					in >> data[1] >> data[2] >> data[3];
					data[0] = solId++;
					data[4] = popId;
					arributes.push_back(data);
				}
				in.close();
			}

			solId = sols.size() + 1;

			//int stop = -1;
			{
				std::ifstream in(dirPath + "/" + filenameIter + filefix2);
				while (in >> inputId) {
					std::getline(in, line);
					sols.push_back(line);
				}

				in.close();
			}



		}



		{
			std::ofstream out(targetDir + "/" + dirname + filefix1);

			for (auto& iter : headers) {
				out << iter << "\t";
			}
			out << std::endl;
			for (auto& it : arributes) {
				for (auto& it2 : it) {
					out << it2 << "\t";
				}
				out << std::endl;
			}

			out.close();
		}


		{
			std::ofstream out(targetDir + "/" + dirname + filefix2);
			int solId(0);
			for (auto& it : sols) {
				out << ++solId << it << std::endl;
			}

			out.close();
		}
	}

}





int main() {

	//readFile("//172.24.24.151/share/student/2018/diaoyiya/code/EAX_GEN/data","E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data");

	//std::cout << "yes" << std::endl;
	//runAlg();


	//mergeData("//172.24.24.151/share/student/2018/diaoyiya/code/EAX_GEN/data", "E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data");


	LDEE_path = "E:/Diao_Yiya/paper/fitness_landscape/IDEE/ldee_tools/";
	filedir = "E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data";

	getFilenames(filedir);

	//result_path = "E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_visualization_data";
	//runLDEE("E:/Diao_Yiya/paper/fitness_landscape/IDEE/idee_data", "TSP_att532_NPOP_100_NCH_30_tspPop_0");




	//std::vector<std::thread> thrds;();
	//int num_task = std::thread::hardware_concurrency();
	////   curInfo = { 1,1,11 };
 //  //    runEAX_GA(curInfo[0], curInfo[1], curInfo[2], double(curInfo[2] + 1) / double(maxRun + 1));




	//for (int idx(0); idx < filenames.size() / 2; ++idx) {
	//	filenames[idx] = filenames[idx * 2];
	//}
	//filenames.resize(filenames.size() / 2);
	//swap(filenames.front(), filenames.back());
	//filenames.pop_back();
	runThread();

	//for (size_t i = 0; i < num_task; ++i) {
	//	thrds.push_back(std::thread(
	//		runThread));}
	//for (auto& thrd : thrds)
	//	thrd.join();
	return 0;
}