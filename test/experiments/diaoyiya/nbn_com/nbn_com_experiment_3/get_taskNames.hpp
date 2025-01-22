

struct Info {
	int m_suc_lkh = 0;
	int m_suc_eax = 0;
	std::string m_tspname;
	// ���� < ��������ṩ�ȽϷ���
	bool operator<(const Info& other) const {
		// ���ȱȽ� m_suc_lkh
		if (m_suc_lkh < other.m_suc_lkh) return true;
		if (m_suc_lkh > other.m_suc_lkh) return false;

		// m_suc_lkh ��ȣ���Ƚ� m_suc_eax
		if (m_suc_eax < other.m_suc_eax) return true;
		if (m_suc_eax > other.m_suc_eax) return false;

		// m_suc_eax Ҳ��ȣ���Ƚ� m_tspname
		return m_tspname < other.m_tspname;
	}


};
void readTxtData(const std::string& filepath, std::map<std::pair<int, int>, Info>& totalInfos) {


	std::string tspname;
	int a(0), b(0);

	std::ifstream in(filepath);
	std::pair<int, int> key;
	Info curinfo;
	while (in >> tspname >> a >> b) {
		//	tspname = tspname.substr(0, tspname.find_first_of("."));
		//	std::cout << tspname << "\t" << a << "\t" << b << std::endl;
		key.first = a;
		key.second = b;
		curinfo.m_suc_lkh = a;
		curinfo.m_suc_eax = b;
		curinfo.m_tspname = tspname;
		totalInfos[key] = curinfo;
	}
	in.close();

	//int stop = -1;
}
void readTspInstance1(std::vector<std::string>& tspnames) {

	std::string dirpath = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp_expriment_numeric_methods/tsp_result2";


	std::map<std::pair<int, int>, Info> totalInfos;

	namespace fs = std::filesystem;
	// Iterate over the directory
	for (const auto& entry : fs::recursive_directory_iterator(dirpath)) {
		if (fs::is_regular_file(entry)) {
			std::cout << "File found: " << entry.path() << std::endl;

			readTxtData(entry.path().string(), totalInfos);
		}
	}

	//std::map<std::pair<int, int>, Info> totalInfos2;
	//readTxtData(dirpath + "/experiments.txt", totalInfos2);

	//for (auto& it : totalInfos2) {
	//	totalInfos[it.first] = it.second;
	//}


	//std::set<Info> sortedInfos;
	//std::vector<Info> tasks;




	//std::map<int, Info> firstInfos;
	//for (auto& it : totalInfos) {
	//	firstInfos[it.first.first] = it.second;
	//}

	//for (auto& it : firstInfos) {
	//	sortedInfos.insert(it.second);
	//	tasks.push_back(it.second);
	//	//	tspnames.push_back(it.second.m_tspname);
	//}

	//firstInfos.clear();

	//for (auto& it : totalInfos) {
	//	firstInfos[it.first.second] = it.second;
	//}

	//for (auto& it : firstInfos) {
	//	sortedInfos.insert(it.second);
	//	tasks.push_back(it.second);
	//	//	tspnames.push_back(it.second.m_tspname);
	//}


	//std::set<Info> setA;
	//std::set<Info> setB;

	//for (auto& it : totalInfos) {
	//	setA.insert(it.second);
	//}
	//for (auto& it : totalInfos2) {
	//	setB.insert(it.second);
	//}



	//std::set<Info> difference; // ���ڴ洢��Ľ��

	//// �������� setA �е����� setB �е�����Ԫ�ز��뵽 difference ��
	//std::set_difference(setA.begin(), setA.end(),
	//	setB.begin(), setB.end(),
	//	std::inserter(difference, difference.begin()));


	//std::set<Info> difference2; // ���ڴ洢��Ľ��

	//// �������� setA �е����� setB �е�����Ԫ�ز��뵽 difference ��
	//std::set_difference(difference.begin(), difference.end(),
	//	sortedInfos.begin(), sortedInfos.end(),
	//	std::inserter(difference2, difference2.begin()));



	//for (auto& it : difference2) {
	//	tasks.push_back(it);
	//}



	tspnames.clear();
	//std::ofstream out(dirpath + "/experiments2.txt");
	for (auto& it : totalInfos) {
		//out << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		std::cout << it.second.m_tspname << "\t" << it.second.m_suc_lkh << "\t" << it.second.m_suc_eax << std::endl;
		tspnames.push_back(it.second.m_tspname);
	}
	//out.close();




}
void readTspInstance(std::vector<std::string>& tspnames) {

	std::string dirpath = ofec::g_working_directory + "/paper_com_experiment_data/totalTsp_expriment_numeric_methods/tsp_result2";


	std::map<std::pair<int, int>, Info> totalInfos;

	namespace fs = std::filesystem;
	// Iterate over the directory
	for (const auto& entry : fs::recursive_directory_iterator(dirpath)) {
		if (fs::is_regular_file(entry)) {
			std::cout << "File found: " << entry.path() << std::endl;

			readTxtData(entry.path().string(), totalInfos);
		}
	}

	std::map<std::pair<int, int>, Info> totalInfos2;
	readTxtData(dirpath + "/experiments.txt", totalInfos2);

	for (auto& it : totalInfos2) {
		totalInfos[it.first] = it.second;
	}


	std::set<Info> sortedInfos;
	std::vector<Info> tasks;




	std::map<int, Info> firstInfos;
	for (auto& it : totalInfos) {
		firstInfos[it.first.first] = it.second;
	}

	for (auto& it : firstInfos) {
		sortedInfos.insert(it.second);
		tasks.push_back(it.second);
		//	tspnames.push_back(it.second.m_tspname);
	}

	firstInfos.clear();

	for (auto& it : totalInfos) {
		firstInfos[it.first.second] = it.second;
	}

	for (auto& it : firstInfos) {
		sortedInfos.insert(it.second);
		tasks.push_back(it.second);
		//	tspnames.push_back(it.second.m_tspname);
	}


	std::set<Info> setA;
	std::set<Info> setB;

	for (auto& it : totalInfos) {
		setA.insert(it.second);
	}
	for (auto& it : totalInfos2) {
		setB.insert(it.second);
	}



	std::set<Info> difference; // ���ڴ洢��Ľ��

	// �������� setA �е����� setB �е�����Ԫ�ز��뵽 difference ��
	std::set_difference(setA.begin(), setA.end(),
		setB.begin(), setB.end(),
		std::inserter(difference, difference.begin()));


	std::set<Info> difference2; // ���ڴ洢��Ľ��

	// �������� setA �е����� setB �е�����Ԫ�ز��뵽 difference ��
	std::set_difference(difference.begin(), difference.end(),
		sortedInfos.begin(), sortedInfos.end(),
		std::inserter(difference2, difference2.begin()));



	for (auto& it : difference2) {
		tasks.push_back(it);
	}



	tspnames.clear();
	//std::ofstream out(dirpath + "/experiments2.txt");
	for (auto& it : tasks) {
		//out << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		std::cout << it.m_tspname << "\t" << it.m_suc_lkh << "\t" << it.m_suc_eax << std::endl;
		tspnames.push_back(it.m_tspname);
	}
	//out.close();




}


