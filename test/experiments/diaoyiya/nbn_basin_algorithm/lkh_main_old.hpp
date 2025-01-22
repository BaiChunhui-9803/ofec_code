#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Genetic.h"
#include "./INCLUDE/BIT.h"
#include "./INCLUDE/Heap.h"
#include "./INCLUDE/Segment.h"
#include <string.h>
#include "../../../../../utility/hash_table/hash_table.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include <list>
#include <map>

/*
 * This file contains the main function of the program.
 */

namespace LKH {
	//thread_local LKHAlg *lkh_ptr = 0;



	char* LKHAlg::lkh_strtok(char* str, const char* delim, char** saveptr) {
#ifdef _WIN32
		return strtok_s(str, delim, saveptr);
		//define something for Windows (32-bit and 64-bit, this part is common)
#elif linux
		// linux
		return strtok_r(str, delim, saveptr);
#elif unix // all unices not caught above
		return strtok_r(str, delim, saveptr);
		// Unix
#endif
		return 0;
	}

	void LKHAlg::readProblem(const std::string& tspDir,
		const std::string& tspname) {
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();
		lkh_ptr = this;
		this->m_filename = tspname;
		//ReadParameters();
		setDefaultParaments();
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();
		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		ReadProblem(tspDir + tspname + ".tsp");
	}


	void LKHAlg::initBestInsideRecorder(const std::vector<int>& initSol, GainType cost) {
		m_bestCostInside = std::numeric_limits<GainType>::max();
		m_bestSolSameEdge = DimensionSaved;
		m_curSameEdge = DimensionSaved;
		m_curCost = cost;
		mvv_diff_matrix.resize(DimensionSaved + 1);
		for (auto& it : mvv_diff_matrix) it.resize(DimensionSaved + 1);
		for (auto& it : mvv_diff_matrix) std::fill(it.begin(), it.end(), 0);
		initDiffMatrix(initSol);
		initDiffMatrix(initSol);
		std::vector<int> initSolVec(initSol);
		for (auto& it : initSolVec) --it;
		ofec::TravellingSalesman::transferEdgeSol(initSolVec, m_curSolRecord);

	}

	void LKHAlg::initDiffMatrix(const std::vector<int>& initSol) {
		for (int idx(0); idx < initSol.size(); ++idx) {
			++mvv_diff_matrix[initSol[idx]][initSol[(idx + 1) % DimensionSaved]];
			++mvv_diff_matrix[initSol[(idx + 1) % DimensionSaved]][initSol[idx]];

		}
	}



	//void getBestInsideSol(const std::vector<int>& initSol, std::vector<int>& bestSol) {
	//	for (int idx(0); idx < initSol.size(); ++idx) {
	//		--mvv_best_diff_matrix[initSol[idx]][initSol[(idx + 1) % DimensionSaved]];
	//		--mvv_best_diff_matrix[initSol[(idx + 1) % DimensionSaved]][initSol[idx]];
	//	}
	//	bestSol.clear();
	//	//bestSol.resize(DimensionSaved);
	//	bestSol.push_back(1);
	//	std::vector<bool> visited(DimensionSaved+1, false);
	//	visited[bestSol.front()] = true;
	//	int from = bestSol.front();
	//	for (int idx(1); idx < DimensionSaved; ++idx) {
	//		for (int idy(1); idy <= DimensionSaved; ++idy) {
	//			if (mvv_best_diff_matrix[from][idy] == 1 && !visited[idy]) {
	//				bestSol.push_back(idy);
	//				visited[idy] = true;
	//				from = idy;
	//				break;
	//			}
	//		}
	//	}
	//	
	//}

	//void getCurInsideSol(const std::vector<int>& initSol, std::vector<int>& bestSol) {
	//	for (int idx(0); idx < initSol.size(); ++idx) {
	//		--mvv_diff_matrix[initSol[idx]][initSol[(idx + 1) % DimensionSaved]];
	//		--mvv_diff_matrix[initSol[(idx + 1) % DimensionSaved]][initSol[idx]];
	//	}
	//	bestSol.clear();
	//	//bestSol.resize(DimensionSaved);
	//	bestSol.push_back(1);
	//	std::vector<bool> visited(DimensionSaved + 1, false);
	//	visited[bestSol.front()] = true;
	//	int from = bestSol.front();
	//	for (int idx(1); idx < DimensionSaved; ++idx) {
	//		for (int idy(1); idy <= DimensionSaved; ++idy) {
	//			if (mvv_diff_matrix[from][idy] == 1 && !visited[idy]) {
	//				bestSol.push_back(idy);
	//				visited[idy] = true;
	//				from = idy;
	//				break;
	//			}
	//		}
	//	}

	//}



	/*
			* Edges (t1,t2) and (t3,t4)
* are exchanged with edges (t2,t3) and (t4,t1).*/
	void LKHAlg::updateCurSolRecord(int t1, int t2, int t3, int t4) {

		std::vector<std::array<int, 2>> fnode_discard_edges;
		//std::vector<std::array<Node*, 2>> fnode_connect_edges;
		fnode_discard_edges.push_back({ t1,t2 });
		fnode_discard_edges.push_back({ t2,t1 });
		fnode_discard_edges.push_back({ t3,t4 });
		fnode_discard_edges.push_back({ t4,t3 });
		int fromId(0), toId(0);
		for (auto& it : fnode_discard_edges) {
			fromId = it.front() - 1;
			toId = it.back() - 1;
			if (m_curSolRecord[fromId].front() == toId) {
				swap(m_curSolRecord[fromId].front(), m_curSolRecord[fromId].back());
			}
		}
		fnode_discard_edges.clear();
		fnode_discard_edges.push_back({ t2,t3 });
		fnode_discard_edges.push_back({ t3,t2 });
		fnode_discard_edges.push_back({ t4,t1 });
		fnode_discard_edges.push_back({ t1,t4 });

		for (auto& it : fnode_discard_edges) {
			fromId = it.front() - 1;
			toId = it.back() - 1;
			m_curSolRecord[fromId].back() = toId;

		}
	}

	/*
			* Edges (t1,t2) and (t3,t4)
* are exchanged with edges (t2,t3) and (t4,t1).*/

	void LKHAlg::updateDiffDis(int t1, int t2, int t3, int t4) {
		int diff = 0;
		std::array<int, 4> fnode = { t1,t2,t3,t4 };
		for (int idx(0); idx < 2; ++idx) {
			if (--mvv_diff_matrix[fnode[idx * 2]][fnode[idx * 2 + 1]] == 1) {
				--diff;
			}
			if (--mvv_diff_matrix[fnode[idx * 2 + 1]][fnode[idx * 2]] == 1) {
				--diff;
			}
		}
		for (int idx(0); idx < 2; ++idx) {
			if (++mvv_diff_matrix[fnode[idx * 2 + 1]][fnode[(idx * 2 + 2) % 4]] == 2) {
				++diff;
			}
			if (++mvv_diff_matrix[fnode[(idx * 2 + 2) % 4]][fnode[idx * 2 + 1]] == 2) {
				++diff;
			}
		}
		diff /= 2;
		m_curSameEdge += diff;
	}


	void LKHAlg::set2optLocalSearchParamenters(const std::string& tspfile) {

		//ReadParameters();
		setDefaultParaments();
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();
		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		ReadProblem(tspfile);
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;
		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();

		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		//	Precision = 1;
	}


	int LKHAlg::randomSamplingWith2optHyperSamplingSpaceDivision(
		const std::string& saveDir,
		const std::string& tspDir,
		const std::string& tspname,
		double radius,
		int totalSample) {
		using namespace ofec;
		this->m_filename = tspname;

		lkh_ptr = this;
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();

		set2optLocalSearchParamenters(tspDir + tspname + ".tsp");



		std::vector<std::vector<int>> hyperSamplingVisited(DimensionSaved + 1, std::vector<int>(DimensionSaved + 1));
		std::vector<bool> curVisited(DimensionSaved + 1);

		AllocateStructures();
		CreateCandidateSet();
		InitializeStatistics();

		BestCost = std::numeric_limits<GainType>::max();

		utility::HashTable<std::vector<int>> lon_set;
		utility::HashTable<std::vector<int>> lon_opt_set;
		lon_set.initiazlize();
		lon_opt_set.initiazlize();

		nodeInfo sampledSol;
		nodeInfo opt2Sol;

		ofstream nodeInfoOut(saveDir + tspname + ".sampled2optInfo");
		sampledSol.outputNodeInfoHeader(nodeInfoOut);

		ofstream nodeRandomSolOut(saveDir + tspname + "_random.sols");
		opt2Sol.outputSolHeader(nodeRandomSolOut);

		ofstream nodeOptSolOut(saveDir + tspname + "_optInside.sols");
		opt2Sol.outputSolHeader(nodeOptSolOut);


		ofec::ParameterMap params;
		params["problem name"] = std::string("TSP");
		params["dataFile1"] = tspname;

		m_problem = ofec::Problem::generateByFactory("TSP");
		m_problem->getInputParameters().input(params);
		m_problem->initialize(0.1);

		//nodeSolOut << "ID\tSOLUTION" << std::endl;
		int totalSamples = totalSample;
		int curId(0);
		//	std::vector<std::vector<int>> isRandom(DimensionSaved, std::vector<int>(DimensionSaved, 0));

		std::vector<int> initSol(DimensionSaved, 0);
		for (int idx(1); idx <= DimensionSaved; ++idx) {
			initSol[idx - 1] = idx;
		}
		std::shared_ptr<ofec::Random> rnd = std::make_shared<ofec::Random>(0.5);


		m_max_radius = DimensionSaved * radius;
		m_flag_divideBySol = true;

		/* Find a specified number (Runs) of local optima */
		for (Run = 1; Run <= totalSamples; Run++) {

			if ((Run - 1) % 1000 == 0) {
				std::cout << "tspname\t" << tspname << "\t" << Run << std::endl;
			}

			//rnd->uniform.shuffle(initSol.begin(), initSol.end());
			sampledSol.m_sampleID = Run;
			generateRandomSolHyperSampling(initSol,
				hyperSamplingVisited,
				curVisited, rnd.get(), DimensionSaved);
			sampledSol.cost = generateSol(initSol);
			initBestInsideRecorder(initSol, sampledSol.cost);
			//sampledSol.cost = generateSolRandom();
			RecordBetterTour();
			transferTourType(BetterTour, sampledSol.m_sol);
			m_test_initSol = initSol;

			unsigned curhash = calHash(sampledSol.m_sol);
			sampledSol.is_exist = lon_set.insertValue(sampledSol.m_sol, curhash, curId);
			sampledSol.m_nodeId = curId + 1;

			opt2Sol.cost = LinKernighanLocalSearch();
			ofec::TravellingSalesman::transferEdgeSol(m_bestSolInsideRecord, opt2Sol.m_sol);
			opt2Sol.cost = m_bestCostInside;
			{
				int cost = 0;
				//for (auto& it : opt2Sol.m_sol) ++it;
				for (int idx(0); idx < opt2Sol.m_sol.size(); ++idx) {
					cost += CAST_TSP(m_problem.get())->cost()[opt2Sol.m_sol[idx]][opt2Sol.m_sol[(idx + 1) % opt2Sol.m_sol.size()]];
				}

				//auto cost = generateSol(opt2Sol.m_sol);
				if (cost != opt2Sol.cost) {
					std::cout << "differ at cost \t Run" << Run << "\tat\tProname\t" << m_filename << std::endl;
				}
			}
			{
				for (auto& it : sampledSol.m_sol) --it;
				//	for (auto& it : opt2Sol.m_sol) --it;
					//for(auto& it:sampe)
				std::vector<std::array<int, 2>> edgex;
				std::vector<std::array<int, 2>> edgey;
				ofec::TravellingSalesman::transferEdgeSol(sampledSol.m_sol, edgex);
				ofec::TravellingSalesman::transferEdgeSol(opt2Sol.m_sol, edgey);


				if (int(ofec::TravellingSalesman::tspDis(edgex, edgey)) + m_bestSolSameEdge != DimensionSaved) {
					std::cout << "error at tspDis \t" << Run << std::endl;
				}
			}

			curhash = calHash(opt2Sol.m_sol);
			opt2Sol.is_exist = lon_opt_set.insertValue(opt2Sol.m_sol, curhash, curId);
			opt2Sol.m_nodeId = curId + 1;

			sampledSol.m_OptId = opt2Sol.m_OptId = opt2Sol.m_nodeId;

			if (!sampledSol.is_exist) {
				sampledSol.outputNodeInfo(nodeInfoOut);
				sampledSol.outputSol(nodeRandomSolOut);
			}
			if (!opt2Sol.is_exist) {
				opt2Sol.outputNodeInfo(nodeInfoOut);
				opt2Sol.outputSol(nodeOptSolOut);
			}
		}

		nodeInfoOut.close();
		nodeRandomSolOut.close();
		nodeOptSolOut.close();
		//out.close();
		freeAll();


		delete ParameterFileName;
		ParameterFileName = nullptr;
		return 0;
	}



	int LKHAlg::randomSamplingWith2optHyperSampling(
		const std::string& saveDir,
		const std::string& tspDir,
		const std::string& tspname, int totalSample) {
		using namespace ofec;
		this->m_filename = tspname;

		lkh_ptr = this;
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();

		set2optLocalSearchParamenters(tspDir + tspname + ".tsp");
		std::vector<std::vector<int>> hyperSamplingVisited(DimensionSaved + 1, std::vector<int>(DimensionSaved + 1));
		std::vector<bool> curVisited(DimensionSaved + 1);

		AllocateStructures();
		CreateCandidateSet();
		InitializeStatistics();

		BestCost = std::numeric_limits<GainType>::max();

		utility::HashTable<std::vector<int>> lon_set;
		lon_set.initiazlize();

		nodeInfo sampledSol;
		nodeInfo opt2Sol;

		ofstream nodeInfoOut(saveDir + tspname + ".sampled2optInfo");
		sampledSol.outputNodeInfoHeader(nodeInfoOut);

		ofstream nodeRandomSolOut(saveDir + tspname + "_random.sols");
		opt2Sol.outputSolHeader(nodeRandomSolOut);

		ofstream nodeOptSolOut(saveDir + tspname + "_opt.sols");
		opt2Sol.outputSolHeader(nodeOptSolOut);


		//nodeSolOut << "ID\tSOLUTION" << std::endl;
		int totalSamples = totalSample;
		int curId(0);
		//	std::vector<std::vector<int>> isRandom(DimensionSaved, std::vector<int>(DimensionSaved, 0));

		std::vector<int> initSol(DimensionSaved, 0);
		for (int idx(1); idx <= DimensionSaved; ++idx) {
			initSol[idx - 1] = idx;
		}
		std::shared_ptr<ofec::Random> rnd = std::make_shared<ofec::Random>(0.5);

		/* Find a specified number (Runs) of local optima */
		for (Run = 1; Run <= totalSamples; Run++) {

			if ((Run - 1) % 1000 == 0) {
				std::cout << "tspname\t" << tspname << "\t" << Run << std::endl;
			}

			rnd->uniform.shuffle(initSol.begin(), initSol.end());
			sampledSol.m_sampleID = Run;

			generateRandomSolHyperSampling(initSol,
				hyperSamplingVisited,
				curVisited, rnd.get(), DimensionSaved);
			sampledSol.cost = generateSol(initSol);
			//sampledSol.cost = generateSolRandom();
			RecordBetterTour();
			transferTourType(BetterTour, sampledSol.m_sol);
			/*		transferSol(initSol);
					if (initSol != sampledSol.m_sol) {
						std::cout << "error" << std::endl;
					}*/
					//for (int idx(0); idx < sampledSol.m_sol.size(); ++idx) {
					//	int curCity = sampledSol.m_sol[idx] - 1;
					//	isRandom[curCity][sampledSol.m_sol[(idx + 1) % sampledSol.m_sol.size()] - 1]++;
					//	isRandom[sampledSol.m_sol[(idx + 1) % sampledSol.m_sol.size()] - 1][curCity]++;
					//	isRandom[curCity][sampledSol.m_sol[(idx - 1 + sampledSol.m_sol.size()) % sampledSol.m_sol.size()] - 1]++;
					//	isRandom[sampledSol.m_sol[(idx - 1 + sampledSol.m_sol.size()) % sampledSol.m_sol.size()] - 1][curCity]++;

					//}


			unsigned curhash = calHash(sampledSol.m_sol);
			sampledSol.is_exist = lon_set.insertValue(sampledSol.m_sol, curhash, curId);
			sampledSol.m_nodeId = curId + 1;

			opt2Sol.cost = LinKernighanLocalSearch();
			RecordBetterTour();
			transferTourType(BetterTour, opt2Sol.m_sol);
			curhash = calHash(opt2Sol.m_sol);
			opt2Sol.is_exist = lon_set.insertValue(opt2Sol.m_sol, curhash, curId);
			opt2Sol.m_nodeId = curId + 1;

			sampledSol.m_OptId = opt2Sol.m_OptId = opt2Sol.m_nodeId;

			if (!sampledSol.is_exist) {
				sampledSol.outputNodeInfo(nodeInfoOut);
				sampledSol.outputSol(nodeRandomSolOut);
			}
			if (!opt2Sol.is_exist) {
				opt2Sol.outputNodeInfo(nodeInfoOut);
				opt2Sol.outputSol(nodeOptSolOut);
			}
		}
		//for (int idx(0); idx < isRandom.size(); ++idx) {
		//	std::cout << "idx\t" << idx << "\t";
		//	double mean(0);
		//	double standard_derivation(0);
		//	ofec::calMeanAndStd<int, double>(isRandom[idx], mean, standard_derivation);

		//	double minVal(1e9);
		//	double maxVal(0);
		//	for (auto& it : isRandom[idx]) {
		//		minVal = std::min<double>(it, minVal);
		//		maxVal = std::max<double>(it, maxVal);

		//	}
		//	std::cout << "max\t" << maxVal << "\t min\t" << minVal << "\t";
		//	std::cout << "mean\t" << mean << "\t" << standard_derivation << std::endl;
		//}

		nodeInfoOut.close();
		nodeRandomSolOut.close();
		nodeOptSolOut.close();
		//out.close();
		freeAll();
		return 0;
	}

	int LKHAlg::randomSamplingWith2opt(const std::string& saveDir,
		const std::string& tspDir,
		const std::string& tspname, int totalSample) {
		using namespace ofec;
		this->m_filename = tspname;



		lkh_ptr = this;
		//ParameterFileName = new char[parFile.size() + 1];
		//strcpy(ParameterFileName, parFile.c_str());
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();

		//ReadParameters();
		setDefaultParaments();
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();
		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		ReadProblem(tspDir + tspname + ".tsp");
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();

		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		//	Precision = 1;






		AllocateStructures();
		CreateCandidateSet();
		InitializeStatistics();

		BestCost = std::numeric_limits<GainType>::max();




		utility::HashTable<std::vector<int>> lon_set;
		lon_set.initiazlize();

		nodeInfo sampledSol;
		nodeInfo opt2Sol;

		ofstream nodeInfoOut(saveDir + tspname + ".sampled2optInfo");
		sampledSol.outputNodeInfoHeader(nodeInfoOut);

		ofstream nodeRandomSolOut(saveDir + tspname + "_random.sols");
		opt2Sol.outputSolHeader(nodeRandomSolOut);

		ofstream nodeOptSolOut(saveDir + tspname + "_opt.sols");
		opt2Sol.outputSolHeader(nodeOptSolOut);


		//nodeSolOut << "ID\tSOLUTION" << std::endl;
		int totalSamples = totalSample;
		int curId(0);
		std::vector<std::vector<int>> isRandom(DimensionSaved, std::vector<int>(DimensionSaved, 0));
		std::vector<int> initSol(DimensionSaved, 0);
		for (int idx(1); idx <= DimensionSaved; ++idx) {
			initSol[idx - 1] = idx;
		}
		std::shared_ptr<ofec::Random> rnd = std::make_shared<ofec::Random>(0.5);

		/* Find a specified number (Runs) of local optima */
		for (Run = 1; Run <= totalSamples; Run++) {

			if ((Run - 1) % 1000 == 0) {
				std::cout << "tspname\t" << tspname << "\t" << Run << std::endl;
			}

			rnd->uniform.shuffle(initSol.begin(), initSol.end());
			sampledSol.m_sampleID = Run;
			sampledSol.cost = generateSol(initSol);

			RecordBetterTour();
			transferTourType(BetterTour, sampledSol.m_sol);
			/*		transferSol(initSol);
					if (initSol != sampledSol.m_sol) {
						std::cout << "error" << std::endl;
					}*/

			for (int idx(0); idx < sampledSol.m_sol.size(); ++idx) {
				int curCity = sampledSol.m_sol[idx] - 1;
				isRandom[curCity][sampledSol.m_sol[(idx + 1) % sampledSol.m_sol.size()] - 1]++;
				isRandom[curCity][sampledSol.m_sol[(idx - 1 + sampledSol.m_sol.size()) % sampledSol.m_sol.size()] - 1]++;
			}



			unsigned curhash = calHash(sampledSol.m_sol);
			sampledSol.is_exist = lon_set.insertValue(sampledSol.m_sol, curhash, curId);
			sampledSol.m_nodeId = curId + 1;

			opt2Sol.cost = LinKernighanLocalSearch();
			RecordBetterTour();
			transferTourType(BetterTour, opt2Sol.m_sol);
			curhash = calHash(opt2Sol.m_sol);
			opt2Sol.is_exist = lon_set.insertValue(opt2Sol.m_sol, curhash, curId);
			opt2Sol.m_nodeId = curId + 1;

			sampledSol.m_OptId = opt2Sol.m_OptId = opt2Sol.m_nodeId;

			if (!sampledSol.is_exist) {
				sampledSol.outputNodeInfo(nodeInfoOut);
				sampledSol.outputSol(nodeRandomSolOut);
			}


			if (!opt2Sol.is_exist) {
				opt2Sol.outputNodeInfo(nodeInfoOut);
				opt2Sol.outputSol(nodeOptSolOut);
			}



		}

		for (int idx(0); idx < isRandom.size(); ++idx) {
			std::cout << "idx\t" << idx << "\t";
			double mean(0);
			double standard_derivation(0);
			ofec::calMeanAndStd<int, double>(isRandom[idx], mean, standard_derivation);

			double minVal(1e9);
			double maxVal(0);
			for (auto& it : isRandom[idx]) {
				minVal = std::min<double>(it, minVal);
				maxVal = std::max<double>(it, maxVal);

			}
			std::cout << "max\t" << maxVal << "\t min\t" << minVal << "\t";
			std::cout << "mean\t" << mean << "\t" << standard_derivation << std::endl;

		}


		//	SRandom(++Seed);

		nodeInfoOut.close();
		nodeRandomSolOut.close();
		nodeOptSolOut.close();
		//out.close();
		freeAll();
		return 0;
	}


	int LKHAlg::lon_sampling_2opt_3kickSave(
		std::vector<int>& nodeFit,
		std::vector<TraceEdge>& edge_history,
		//	std::vector<TraceEdge>& edges,
		std::vector<std::vector<int>>& sols,
		//	const std::string& saveDir,
		const std::string& tspDir,
		const std::string& tspname) {
		using namespace ofec;
		lkh_ptr = this;
		this->m_filename = tspname;
		//ParameterFileName = new char[parFile.size() + 1];
		//strcpy(ParameterFileName, parFile.c_str());
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();

		//ReadParameters();
		setDefaultParaments();
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;
		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();
		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		ReadProblem(tspDir + tspname + ".tsp");
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();

		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		//	Precision = 1;


		AllocateStructures();
		CreateCandidateSet();
		InitializeStatistics();

		BestCost = std::numeric_limits<GainType>::max();
		utility::HashTable<std::vector<int>> lon_set;
		lon_set.initiazlize();
		int chainedLKLoop = 1e4;


		std::vector<std::vector<int>> hyperSamplingVisited(DimensionSaved + 1, std::vector<int>(DimensionSaved + 1));
		std::vector<bool> curVisited(DimensionSaved + 1);

		//std::vector<TraceEdge> edge_history;
		//std::vector<long long> nodeFit;
		//int maxId(0);


		//ofstream nodeFitOut(saveDir + tspname + ".nodes");
		//nodeFitOut << "ID\tFITNESS" << std::endl;

		//ofstream nodeSol(saveDir + tspname + ".sols");
		//nodeSol << "ID\tSOLUTION" << std::endl;

		nodeFit.push_back(0);
		sols.emplace_back(std::vector<int>());



		GainType fromSolCost = 0;
		GainType toSolCost = 0;
		int fromSolId = 0;
		int toSolId = 0;
		std::vector<int> sol(DimensionSaved, 0);
		std::vector<int> initSol(DimensionSaved, 0);
		std::shared_ptr<ofec::Random> rnd = std::make_shared<ofec::Random>(0.5);


		TraceEdge curEdge;
		/* Find a specified number (Runs) of local optima */

		for (Run = 1; Run <= Runs; Run++) {
			curEdge.run = Run;
			std::cout << tspname << "\t" << Run << std::endl;
			//curInfo.runId = Run;
			if (Run % 2) {
				InitialTourAlgorithm = QUICK_BORUVKA;
				fromSolCost = GreedyTour();
			}
			else {


				//	sampledSol.m_sampleID = Run;

				generateRandomSolHyperSampling(initSol,
					hyperSamplingVisited,
					curVisited, rnd.get(), DimensionSaved);
				fromSolCost = generateSol(initSol);
				//fromSolCost = generateSolRandom();
			}
			fromSolCost = LinKernighanLocalSearch();
			RecordBetterTour();
			transferTourType(BetterTour, sol);
			unsigned curhash = calHash(sol);
			int curId(0);
			bool flag_exist = lon_set.insertValue(sol, curhash, curId);
			fromSolId = curId + 1;
			if (!flag_exist) {
				//nodeFitOut << fromSolId << "\t" << fromSolCost << std::endl;
				//	std::cout << fromSolId << "\t" << fromSolCost << std::endl;
			//	std::cout << "new sol\t" << fromSolId << "\t" << fromSolCost << std::endl;

				nodeFit.push_back(fromSolCost);
				sols.push_back(sol);
			}
			//maxId = std::max(curId, maxId);
			//bool flag_exist(false);
			for (int j = 0; j < chainedLKLoop; ++j) {
				curEdge.iter = j;
				Node* N = 0, * NextN = 0, * FirstAlternative = 0, * Last = 0;
				for (Last = FirstNode; (N = Last->BestSuc) != FirstNode; Last = N)
					Follow(N, Last);
				KSwapKick(KickType);
				toSolCost = LinKernighanLocalSearch();
				if (fromSolCost >= toSolCost) {
					if (fromSolCost > toSolCost) j = 0;
					RecordBetterTour();
					transferTourType(BetterTour, sol);

					curhash = calHash(sol);
					flag_exist = lon_set.insertValue(sol, curhash, curId);
					toSolId = curId + 1;

					curEdge.id_from = fromSolId;
					curEdge.id_to = toSolId;
					curEdge.num_kick = 1;

					//	if(curEdge.run>=134)
					//	std::cout << curEdge.run << "\t" << curEdge.iter << "\t" << fromSolId << "\t" << fromSolCost << "\t" << toSolId << "\t" << toSolCost << std::endl;


					fromSolId = toSolId;
					fromSolCost = toSolCost;
					if (!flag_exist) {
						nodeFit.push_back(fromSolCost);
						sols.push_back(sol);
						//		std::cout << "new sol\t" << fromSolId << "\t" << fromSolCost << std::endl;
					}
				}
				else {
					curEdge.id_from = fromSolId;
					curEdge.id_to = -1;
					curEdge.num_kick = 1;
				}


				edge_history.push_back(curEdge);
			}
			//	SRandom(++Seed);
		}
		//out.close();
		freeAll();




		//	nodeFitOut.close();
		//	nodeSol.close();


			//ofstream edgeHistoryOut(saveDir + tspname + ".edge_history");
			//edgeHistoryOut << "RUN\tITER\tID_START\tID_END\tNUM_KICKS" << std::endl;

			//for (auto& it : edge_history) {
			//	edgeHistoryOut << it.run << "\t" << it.iter << "\t" << it.id_from << "\t" << it.id_to << "\t" << it.num_kick << std::endl;
			//}
			//edgeHistoryOut.close();
			// 
			//Random *rnd = ADD_RND(0.5);

		return true;
	}


	int LKHAlg::lon_sampling_2opt_3kick(
		const std::string& saveDir,
		const std::string& tspDir, const std::string& tspname) {
		struct TraceEdge {
			int run = 0;
			int iter = 0;
			int id_from = 0;
			int id_to = 0;
			int num_kick = 0;
			int count = 0;
		};
		lkh_ptr = this;
		//ParameterFileName = new char[parFile.size() + 1];
		//strcpy(ParameterFileName, parFile.c_str());


		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();

		//ReadParameters();
		setDefaultParaments();
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();
		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		ReadProblem(tspDir + tspname + ".tsp");
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();

		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		//	Precision = 1;


		AllocateStructures();
		CreateCandidateSet();
		InitializeStatistics();

		BestCost = std::numeric_limits<GainType>::max();
		utility::HashTable<std::vector<int>> lon_set;
		lon_set.initiazlize();



		int chainedLKLoop = 1e4;

		std::vector<TraceEdge> edge_history;

		//std::vector<long long> nodeFit;
		//int maxId(0);


		ofstream nodeFitOut(saveDir + tspname + ".nodes");
		nodeFitOut << "ID\tFITNESS" << std::endl;

		ofstream nodeSol(saveDir + tspname + ".sols");
		nodeSol << "ID\tSOLUTION" << std::endl;


		GainType fromSolCost = 0;
		GainType toSolCost = 0;
		int fromSolId = 0;
		int toSolId = 0;
		std::vector<int> sol(DimensionSaved, 0);
		TraceEdge curEdge;
		/* Find a specified number (Runs) of local optima */

		for (Run = 1; Run <= Runs; Run++) {
			curEdge.run = Run;
			std::cout << tspname << "\t" << Run << std::endl;
			//curInfo.runId = Run;
			if (Run % 2) {
				InitialTourAlgorithm = QUICK_BORUVKA;
				fromSolCost = GreedyTour();
			}
			else {
				fromSolCost = generateSolRandom();
			}
			fromSolCost = LinKernighanLocalSearch();
			RecordBetterTour();
			transferTourType(BetterTour, sol);
			unsigned curhash = calHash(sol);
			int curId(0);
			bool flag_exist = lon_set.insertValue(sol, curhash, curId);
			fromSolId = curId + 1;
			if (!flag_exist) {
				nodeFitOut << fromSolId << "\t" << fromSolCost << std::endl;
				std::cout << fromSolId << "\t" << fromSolCost << std::endl;
				nodeSol << fromSolId;
				for (int idx(0); idx < sol.size(); ++idx) {
					nodeSol << "\t" << sol[idx];
				}
				nodeSol << std::endl;
			}
			//maxId = std::max(curId, maxId);
			//bool flag_exist(false);
			for (int j = 0; j < chainedLKLoop; ++j) {
				curEdge.iter = j;
				Node* N = 0, * NextN = 0, * FirstAlternative = 0, * Last = 0;
				for (Last = FirstNode; (N = Last->BestSuc) != FirstNode; Last = N)
					Follow(N, Last);
				KSwapKick(KickType);
				toSolCost = LinKernighanLocalSearch();
				if (fromSolCost >= toSolCost) {
					if (fromSolCost > toSolCost) j = 0;
					RecordBetterTour();
					transferTourType(BetterTour, sol);

					curhash = calHash(sol);
					flag_exist = lon_set.insertValue(sol, curhash, curId);
					int toSolId = curId + 1;

					curEdge.id_from = fromSolId;
					curEdge.id_to = toSolId;
					curEdge.num_kick = 1;

					//std::cout <<curEdge.run<<"\t"<<curEdge.iter<<"\t" << fromSolId << "\t" << fromSolCost << "\t" << toSolId << "\t" << toSolCost << std::endl;


					fromSolId = toSolId;
					fromSolCost = toSolCost;
					if (!flag_exist) {
						nodeFitOut << fromSolId << "\t" << fromSolCost << std::endl;
						std::cout << fromSolId << "\t" << fromSolCost << std::endl;
						nodeSol << fromSolId;
						for (int idx(0); idx < sol.size(); ++idx) {
							nodeSol << "\t" << sol[idx];
						}
						nodeSol << std::endl;
					}


				}
				else {
					curEdge.id_from = fromSolId;
					curEdge.id_to = -1;
					curEdge.num_kick = 1;
				}
				edge_history.push_back(curEdge);
			}
			SRandom(++Seed);
		}
		//out.close();
		freeAll();




		nodeFitOut.close();
		nodeSol.close();


		ofstream edgeHistoryOut(saveDir + tspname + ".edge_history");
		edgeHistoryOut << "RUN\tITER\tID_START\tID_END\tNUM_KICKS" << std::endl;

		for (auto& it : edge_history) {
			edgeHistoryOut << it.run << "\t" << it.iter << "\t" << it.id_from << "\t" << it.id_to << "\t" << it.num_kick << std::endl;
		}
		edgeHistoryOut.close();
		std::vector<TraceEdge> edges;
		std::map<std::pair<int, int>, int> edgeInfo2EdgeId;
		std::vector<int> edgeVisited;
		int maxId = 0;
		int curId = 0;
		std::pair<int, int> edgeInfo;
		std::vector<std::pair<int, int>> edgeId2EdgeInfo;
		for (auto& it : edge_history) {
			edgeInfo.first = it.id_from;
			edgeInfo.second = it.id_to;
			if (edgeInfo2EdgeId.find(edgeInfo) == edgeInfo2EdgeId.end()) {
				curId = edgeInfo2EdgeId[edgeInfo] = maxId++;
				edgeVisited.push_back(0);
				edgeId2EdgeInfo.push_back(edgeInfo);
			}
			else {
				curId = edgeInfo2EdgeId[edgeInfo];
			}
			++edgeVisited[curId];
			//int curId= edgeInfo2EdgeId[]
		}
		ofstream edgeOut(saveDir + tspname + ".edges");
		edgeOut << "ID_START\tID_END\tCOUNT" << std::endl;
		for (int idx(0); idx < edgeInfo2EdgeId.size(); ++idx) {
			edgeOut << edgeId2EdgeInfo[idx].first << "\t" << edgeId2EdgeInfo[idx].second << "\t";
			edgeOut << edgeVisited[idx] << std::endl;

		}
		//int primeNumber = 
		edgeOut.close();

		return true;
	}


	int LKHAlg::lon_sampling(const std::string& tspFile) {

		lkh_ptr = this;
		//ParameterFileName = new char[parFile.size() + 1];
		//strcpy(ParameterFileName, parFile.c_str());


		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();

		//ReadParameters();
		setDefaultParaments();
		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();
		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;


		ReadProblem(tspFile);


		MaxMatrixDimension = 20000;
		//	MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
		//		MergeWithTourGPX2;


		Excess = 1;
		ExtraCandidates = 0;
		InitialTourAlgorithm = WALK;
		MaxCandidates = 5;
		//MaxSwaps = std::numeric_limits<int>::max();
		MaxSwaps = DimensionSaved;
		MaxTrials = 10000;
		MoveType = 2;

		Kicks = 1;
		KickType = 4;
		MaxBreadth = std::numeric_limits<int>::max();

		SubsequentMoveType = MoveType;
		RestrictedSearch = 0;
		Runs = 1000;
		//	Precision = 1;
		GainType curSolCost = 0;

		AllocateStructures();
		CreateCandidateSet();
		InitializeStatistics();

		BestCost = std::numeric_limits<GainType>::max();


		// for test
		std::vector<int> sol_27686;
		long long sol_test_cost;

		struct RecordInfo {
			int runId = 0;
			long long beforeFit = 0;
			int beforeId = 0;

			long long afterFit = 0;
			int afterId = 0;

			//std::vector<int> sol;

		};


		struct node {
			int nodeId = 0;
			unsigned hashval;
			std::vector<int> sol;
			GainType solCost;
		};

		std::vector<node> m_nodes;
		node curNode;

		std::list<RecordInfo> lon_trait;
		RecordInfo curInfo;
		int chainedLKLoop = 1e4;
		std::vector<int> sol(DimensionSaved, 0);

		utility::HashTable<std::vector<int>> lon_set;
		lon_set.initiazlize();


		int numX = 1e4;

		int maxId(0);

		//	std::ofstream NodeOut("node.txt");
			//HashInitialize(HTable);


		ofstream out("result.txt");
		out << "Run\tFitness\tSolution\tFitness\tNext_Solution" << std::endl;


		/* Find a specified number (Runs) of local optima */

		for (Run = 1; Run <= Runs; Run++) {

			//int Run = 0;
		//	while(maxId< numX){
				//Run++;
			curInfo.runId = Run;
			if (Run % 2) {
				InitialTourAlgorithm = QUICK_BORUVKA;
				curSolCost = GreedyTour();
			}
			else {

				curSolCost = generateSolRandom();
			}
			curSolCost = LinKernighanLocalSearch();
			RecordBetterTour();
			transferTourType(BetterTour, sol);
			unsigned curhash = calHash(sol);


			int curId(0);

			bool flag = lon_set.insertValue(sol, curhash, curId);

			maxId = std::max(curId, maxId);

			//if (!flag) {

			//	curNode.nodeId = curId;
			//	curNode.hashval = curhash;
			//	curNode.sol = sol;
			//	curNode.solCost = curSolCost;

			//	m_nodes.resize(curNode.nodeId + 1);
			//	m_nodes[curNode.nodeId] = curNode;

			//	NodeOut << "curId\t" << curId << "\tsolFit\t" << curSolCost << "\thash\t" << curhash << "\t";
			//	NodeOut << " sol\t";
			//	for (int idx(0); idx < DimensionSaved; ++idx)
			//		NodeOut << sol[idx] << "\t";
			//	NodeOut << std::endl;
			//	
			//}

	//		std::cout << "solFit\t" << curSolCost << "\thash\t" << curhash << "\tcurId\t" << curId << std::endl;

	/*if (flag) {
				Run--;
				continue;
			}*/


			curInfo.beforeFit = curSolCost;
			curInfo.beforeId = curId;

			bool flag_exist(false);
			for (int j = 0; j < chainedLKLoop; ++j) {
				Node* N = 0, * NextN = 0, * FirstAlternative = 0, * Last = 0;
				for (Last = FirstNode; (N = Last->BestSuc) != FirstNode; Last = N)
					Follow(N, Last);
				KSwapKick(KickType);
				curSolCost = LinKernighanLocalSearch();
				if (curInfo.beforeFit >= curSolCost) {
					RecordBetterTour();
					transferTourType(BetterTour, sol);
					curhash = calHash(sol);
					flag_exist = lon_set.insertValue(sol, curhash, curId);
					maxId = std::max(maxId, curId);
					/*				std::cout << "solFit\t" << curSolCost << "\thash\t" << curhash<< "\tcurId\t" << curId << std::endl;
									if (!flag) {


										curNode.nodeId = curId;
										curNode.hashval = curhash;
										curNode.sol = sol;
										curNode.solCost = curSolCost;

										m_nodes.resize(curNode.nodeId + 1);
										m_nodes[curNode.nodeId] = curNode;


										NodeOut << "curId\t" << curId << "\tsolFit\t" << curSolCost << "\thash\t" << curhash << "\t";
										NodeOut << " sol\t";
										for (int idx(0); idx < DimensionSaved; ++idx)
											NodeOut << sol[idx] << "\t";
										NodeOut << std::endl;



										if (m_nodes.size() >= 3) {
											if (m_nodes.front().solCost == m_nodes.back().solCost) {

												bool equal(true);
												int noequalId(0);
												for (int idx(0); idx < DimensionSaved; ++idx) {
													if (m_nodes.front().sol[idx] != m_nodes.back().sol[idx]) {
														equal = false;
														noequalId = idx;
														break;
													}
												}

												int stop = -1;

												if (m_nodes.front().sol == m_nodes.back().sol) {
													if (m_nodes.front().hashval != m_nodes.back().hashval) {

														std::cout << "error" << std::endl;
													}
												}
											}
										}
									}*/

					curInfo.afterFit = curSolCost;
					curInfo.afterId = curId;
					lon_trait.push_back(curInfo);


					{
						auto& it(lon_trait.back());
						out << it.runId << "\t" << it.beforeFit << "\t" << it.beforeId << "\t" << it.afterFit << "\t" << it.afterId << std::endl;

					}
					if (curInfo.beforeFit > curSolCost) {
						j = 0;
						{
							cout << curInfo.runId << "\t" << curInfo.beforeId << "\t" << curInfo.beforeFit << "\t" << curInfo.afterId << "\t" << curInfo.afterFit << std::endl;
						}
					}
					// output info 


		//			if (flag_exist) break;

					curInfo.beforeFit = curInfo.afterFit;
					curInfo.beforeId = curInfo.afterId;


					//	curInfo.sol = sol;
				}


			}

			/*	if (sol_test_cost == curInfo.beforeFit&&sol_27686 == curInfo.sol) {
					std::cout << "error at hash" << std::endl;
				}*/
				//sol_27686 = curInfo.sol;
				//sol_test_cost = curInfo.beforeFit;
			SRandom(++Seed);
		}
		//PrintStatistics();


		//for (auto& it : lon_trait) {
		//	out << it.runId << "\t" << it.beforeFit << "\t" << it.beforeId << "\t" << it.afterFit << "\t" << it.afterId << std::endl;
		//}
		out.close();


		freeAll();

		delete ParameterFileName;
		ParameterFileName = nullptr;

		return 0;
	}

	unsigned LKHAlg::calHash(const std::vector<int>& sol) {

		unsigned Hash(0);
		for (int idx(1); idx < sol.size(); ++idx) {
			Hash ^= Rand[sol[idx - 1]] * Rand[sol[idx]];
		}
		Hash ^= Rand[sol.back()] * Rand[sol.front()];
		return Hash;
	}


	void LKHAlg::transferSol(std::vector<int>& sol) {
		int firstIdx(0);
		for (; firstIdx < sol.size(); ++firstIdx) {
			if (sol[firstIdx] == 1) {
				break;
			}
		}
		int offset = sol[(firstIdx + 1) % sol.size()] < sol[(firstIdx - 1 + sol.size()) % sol.size()] ? 1 : -1;
		std::vector<int> newSol;
		for (int idx(0); idx < sol.size(); ++idx) {
			newSol.push_back(sol[firstIdx]);
			firstIdx = (firstIdx + offset + sol.size()) % sol.size();
		}
		std::swap(sol, newSol);
	}

	void generateRandomSolHyperSampling(std::vector<int>& sol,
		std::vector<std::vector<int>>& hyperSamplingVisited,
		std::vector<bool>& visited,
		ofec::Random* rnd, int DimensionSaved) {
		using namespace ofec;
		int curIdx(0);
		std::fill(visited.begin(), visited.end(), false);

		std::vector<int> minVal(DimensionSaved + 1, 0);
		for (int idx(1); idx <= DimensionSaved; ++idx) {
			for (auto& it : hyperSamplingVisited[idx]) {
				minVal[idx] += std::min(minVal[idx], it);
			}
		}
		//curCity = 1;
		std::vector<int> shuffle_idxs(DimensionSaved);
		for (int idx(0); idx < shuffle_idxs.size(); ++idx) {
			shuffle_idxs[idx] = idx + 1;
		}
		rnd->uniform.shuffle(shuffle_idxs.begin(), shuffle_idxs.end());

		int curCity = shuffle_idxs.front();
		//int maxSum = std::numeric_limits<int>::max();
		//for (auto& cur : shuffle_idxs) {
		//	if (maxSum > minVal[cur]) {
		//		maxSum = minVal[cur];
		//		curCity = cur;
		//	}
		//}



		sol[curIdx] = curCity;
		visited[curCity] = true;
		++curIdx;
		while (curIdx < sol.size()) {

			int nextCity = -1;
			int totalVisited = std::numeric_limits<int>::max();
			for (int idx(0); idx < shuffle_idxs.size(); ++idx) {
				int next = shuffle_idxs[idx];
				if (!visited[next]) {
					if (totalVisited > hyperSamplingVisited[curCity][next]) {
						nextCity = next;
						totalVisited = hyperSamplingVisited[curCity][next];
					}
				}
			}
			sol[curIdx] = nextCity;
			//sol[] = nextCity;
			++hyperSamplingVisited[curCity][nextCity];
			//	++hyperSamplingVisited[nextCity][curCity];
			curCity = nextCity;
			visited[nextCity] = true;
			++curIdx;
		}
	}


	GainType LKHAlg::generateSol(const std::vector<int>& sol) {

		std::vector<int> sol_next(sol.size() + 1, 0);
		for (int idx(0); idx < sol.size(); ++idx) {
			sol_next[sol[idx]] = sol[(idx + 1) % sol.size()];
		}
		int firstIdx = 1;
		Node* Last = FirstNode = &NodeSet[firstIdx];
		Node* N = 0;
		for (int idx(1); idx <= sol.size(); ++idx) {
			N = &NodeSet[sol_next[firstIdx]];
			Follow(N, Last);
			Last = N;
			firstIdx = sol_next[firstIdx];
		}

		GainType Cost = 0;
		N = FirstNode;
		do
			Cost += (this->*C)(N, N->Suc) - N->Pi - N->Suc->Pi;
		while ((N = N->Suc) != FirstNode);
		Cost /= Precision;
		return Cost;
	}


	GainType LKHAlg::generateSolRandom() {
		Node* N = 0, * NextN = 0, * FirstAlternative = 0, * Last = 0;
		Candidate* NN = 0;
		int Alternatives, Count, i;

	Start:
		/* Mark all nodes as "not chosen" by setting their V field to zero */
		N = FirstNode;
		do
			N->V = 0;
		while ((N = N->Suc) != FirstNode);
		Count = 0;

		/* Choose FirstNode without two incident fixed or common candidate edges */
		do {
			if (FixedOrCommonCandidates(N) < 2)
				break;
		} while ((N = N->Suc) != FirstNode);
		if (ProblemType == ATSP && N->Id <= DimensionSaved)
			N += DimensionSaved;
		FirstNode = N;

		/* Move nodes with two incident fixed or common candidate edges in
		   front of FirstNode */
		for (Last = FirstNode->Pred; N != Last; N = NextN) {
			NextN = N->Suc;
			if (FixedOrCommonCandidates(N) == 2)
				Follow(N, Last);
		}

		/* Mark FirstNode as chosen */
		FirstNode->V = 1;
		N = FirstNode;

		/* Loop as long as not all nodes have been chosen */
		while (N->Suc != FirstNode) {
			FirstAlternative = 0;
			Alternatives = 0;
			Count++;

			/* Case A */
			for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
				if (!NextN->V && Fixed(N, NextN)) {
					Alternatives++;
					NextN->Next = FirstAlternative;
					FirstAlternative = NextN;
				}
			}
			if (Alternatives == 0 && MergeTourFiles > 1) {
				for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
					if (!NextN->V && IsCommonEdge(N, NextN)) {
						Alternatives++;
						NextN->Next = FirstAlternative;
						FirstAlternative = NextN;
					}
				}
			}
			if (Alternatives == 0 && FirstNode->InitialSuc && Trial == 1 &&
				Count <= InitialTourFraction * Dimension) {
				/* Case B */
				for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
					if (!NextN->V && InInitialTour(N, NextN)) {
						Alternatives++;
						NextN->Next = FirstAlternative;
						FirstAlternative = NextN;
					}
				}
			}
			if (Alternatives == 0 && Trial > 1 &&
				ProblemType != HCP && ProblemType != HPP) {
				/* Case C */
				for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
					if (!NextN->V && FixedOrCommonCandidates(NextN) < 2 &&
						NN->Alpha == 0 && (InBestTour(N, NextN) ||
							InNextBestTour(N, NextN))) {
						Alternatives++;
						NextN->Next = FirstAlternative;
						FirstAlternative = NextN;
					}
				}
			}
			if (Alternatives == 0) {
				/* Case D */
				for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
					if (!NextN->V && FixedOrCommonCandidates(NextN) < 2) {
						Alternatives++;
						NextN->Next = FirstAlternative;
						FirstAlternative = NextN;
					}
				}
			}
			if (Alternatives == 0) {
				/* Case E (actually not really a random choice) */
				NextN = N->Suc;
				while ((FixedOrCommonCandidates(NextN) == 2 ||
					Forbidden(N, NextN)) && NextN->Suc != FirstNode)
					NextN = NextN->Suc;
				if (FixedOrCommonCandidates(NextN) == 2 || Forbidden(N, NextN)) {
					FirstNode = N;
					goto Start;
				}
			}
			else {
				NextN = FirstAlternative;
				if (Alternatives > 1) {
					/* Select NextN at random among the alternatives */
					i = Random() % Alternatives;
					while (i--)
						NextN = NextN->Next;
				}
			}
			/* Include NextN as the successor of N */
			Follow(NextN, N);
			N = NextN;
			N->V = 1;
		}
		if (Forbidden(N, N->Suc)) {
			FirstNode = N;
			goto Start;
		}

		GainType Cost = 0;
		N = FirstNode;
		do
			Cost += (this->*C)(N, N->Suc) - N->Pi - N->Suc->Pi;
		while ((N = N->Suc) != FirstNode);
		Cost /= Precision;

		return Cost;

	}


	GainType LKHAlg::LinKernighanLocalSearch()
	{
		GainType Cost = 0, Gain = 0, G0 = 0;
		int X2 = 0, i = 0, it = 0;
		Node* t1 = 0, * t2 = 0, * SUCt1 = 0;
		Candidate* Nt1 = 0;
		Segment* S = 0;
		SSegment* SS = 0;
		double EntryTime = GetTime();

		Cost = 0;
		Reversed = 0;
		S = FirstSegment;
		i = 0;
		do {
			S->Size = 0;
			S->Rank = ++i;
			S->Reversed = 0;
			S->First = S->Last = 0;
		} while ((S = S->Suc) != FirstSegment);
		SS = FirstSSegment;
		i = 0;
		do {
			SS->Size = 0;
			SS->Rank = ++i;
			SS->Reversed = 0;
			SS->First = SS->Last = 0;
		} while ((SS = SS->Suc) != FirstSSegment);

		FirstActive = LastActive = 0;
		Swaps = 0;

		/* Compute the cost of the initial tour, Cost.
		   Compute the corresponding hash value, Hash.
		   Initialize the segment list.
		   Make all nodes "active" (so that they can be used as t1). */

		Cost = 0;
		Hash = 0;
		i = 0;
		t1 = FirstNode;
		do {
			t2 = t1->OldSuc = t1->Suc;
			t1->OldPred = t1->Pred;
			t1->Rank = ++i;
			Cost += (t1->SucCost = t2->PredCost = (this->*C)(t1, t2)) - t1->Pi - t2->Pi;
			Hash ^= Rand[t1->Id] * Rand[t2->Id];
			t1->Cost = std::numeric_limits<int>::max();
			for (Nt1 = t1->CandidateSet; (t2 = Nt1->To); Nt1++)
				if (t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
					t1->Cost = Nt1->Cost;
			t1->Parent = S;

			S->Size++;
			if (S->Size == 1)
				S->First = t1;
			S->Last = t1;
			if (SS->Size == 0)
				SS->First = S;
			S->Parent = SS;
			SS->Last = S;
			if (S->Size == GroupSize) {
				S = S->Suc;
				SS->Size++;
				if (SS->Size == SGroupSize)
					SS = SS->Suc;
			}
			t1->OldPredExcluded = t1->OldSucExcluded = 0;
			t1->Next = 0;

			//if (KickType == 0 || Kicks == 0 || Trial == 1 ||
			//	!InBestTour(t1, t1->Pred) || !InBestTour(t1, t1->Suc))
			Activate(t1);

			//if (KickType == 0 || Kicks == 0 || Trial == 1 ||
			//	!InBestTour(t1, t1->Pred) || !InBestTour(t1, t1->Suc))
			//	Activate(t1);
		} while ((t1 = t1->Suc) != FirstNode);
		if (S->Size < GroupSize)
			SS->Size++;
		Cost /= Precision;
		if (TSPTW_Makespan)
			Cost = TSPTW_CurrentMakespanCost = TSPTW_MakespanCost();
		CurrentPenalty = std::numeric_limits<GainType>::max();
		CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
		/*	if (TraceLevel >= 3 ||
				(TraceLevel == 2 &&
				(CurrentPenalty < BetterPenalty ||
					(CurrentPenalty == BetterPenalty && Cost < BetterCost))))*/
					//StatusReport(Cost, EntryTime, "");
		PredSucCostAvailable = 1;
		BIT_Update(this);

		/* Loop as long as improvements are found */
		do {
			/* Choose t1 as the first "active" node */
			while ((t1 = RemoveFirstActive())) {
				/* t1 is now "passive" */
				SUCt1 = SUC(t1);
				/*			if ((TraceLevel >= 3 || (TraceLevel == 2 && Trial == 1)) &&
								++it % (Dimension >= 100000 ? 10000 :
									Dimension >= 10000 ? 1000 : 100) == 0)*/
									/* printff("#%d: Time = %0.2f sec.\n",
											 it, fabs(GetTime() - EntryTime));*/
											 /* Choose t2 as one of t1's two neighbors on the tour */
				for (X2 = 1; X2 <= 2; X2++) {
					t2 = X2 == 1 ? PRED(t1) : SUCt1;
					if (FixedOrCommon(t1, t2) ||
						(RestrictedSearch && Near(t1, t2) &&
							(Trial == 1 ||
								(Trial > BackboneTrials &&
									(KickType == 0 || Kicks == 0)))))
						continue;
					G0 = (this->*C)(t1, t2);
					OldSwaps = Swaps = 0;
					PenaltyGain = Gain = 0;
					/* Try to find a tour-improving chain of moves */
					do
						t2 = Swaps == 0 ? (this->*BestMove)(t1, t2, &G0, &Gain) :
						(this->*BestSubsequentMove)(t1, t2, &G0, &Gain);
					while (t2);
					if (PenaltyGain > 0 || Gain > 0) {
						/* An improvement has been found */
#ifdef HAVE_LONG_LONG
						assert(Gain % Precision == 0);
#else
						assert(fmod(Gain, Precision) == 0);
#endif
						Cost -= Gain / Precision;
						CurrentPenalty -= PenaltyGain;
						StoreTour();
						TSPTW_CurrentMakespanCost = Cost;
						/*			if (TraceLevel >= 3 ||
										(TraceLevel == 2 &&
										(CurrentPenalty < BetterPenalty ||
											(CurrentPenalty == BetterPenalty &&
												Cost < BetterCost))))*/
												//StatusReport(Cost, EntryTime, "");
	//					if (HashSearch(HTable, Hash, Cost))
	//						goto End_LinKernighan;
						/* Make t1 "active" again */
						Activate(t1);
						OldSwaps = 0;
						break;
					}
					OldSwaps = 0;
					RestoreTour();
				}
			}
			//	if (HashSearch(HTable, Hash, Cost))
			//		goto End_LinKernighan;
			//	HashInsert(HTable, Hash, Cost);
				/* Try to find improvements using non-sequential 4/5-opt moves */
			CurrentPenalty = std::numeric_limits<GainType>::max();
			CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
			PenaltyGain = 0;
			//			if (Gain23Used && ((Gain = Gain23()) > 0 || PenaltyGain > 0)) {
			//				/* An improvement has been found */
			//#ifdef HAVE_LONG_LONG
			//				assert(Gain % Precision == 0);
			//#else
			//				assert(fmod(Gain, Precision) == 0);
			//#endif
			//				Cost -= Gain / Precision;
			//				CurrentPenalty -= PenaltyGain;
			//				TSPTW_CurrentMakespanCost = Cost;
			//				StoreTour();
			//				//if (TraceLevel >= 3 ||
			//				//	(TraceLevel == 2 &&
			//				//	(CurrentPenalty < BetterPenalty ||
			//				//		(CurrentPenalty == BetterPenalty && Cost < BetterCost))))
			//				//	//StatusReport(Cost, EntryTime, "+ ");
			//			//	if (HashSearch(HTable, Hash, Cost))
			//			//			goto End_LinKernighan;
			//			}
		} while (PenaltyGain > 0 || Gain > 0);
	End_LinKernighan:
		PredSucCostAvailable = 0;
		NormalizeNodeList();
		NormalizeSegmentList();
		return Cost;
	}


	int LKHAlg::runGA_tsp(const std::string& dataFile,
		const std::string& parFile, std::vector<int>& tour)
	{
		lkh_ptr = this;
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();
		Node* N;
		int i;
		tour.clear();

		/* Read the specification of the problem */
	   /* if (argc >= 2)
			ParameterFileName = argv[1];*/
			//for (int i = 0; i < 5; ++i) {
			//	char para_name[100] = "./DVRP/fcm_twl/cluster";
			//	//para_name =;
			//	char str[1] = { 0 };
			//	itoa(i, str, 10);
			//	strcat(para_name, str);
			//	strcat(para_name, ".par");

		//std::string parFile2 = "./DVRP/temp file/general.par";

		ParameterFileName = new char[parFile.size() + 1];
		strcpy(ParameterFileName, parFile.c_str());

		//ParameterFileName = "./DVRP/shortest_path/ft70.par";
		ReadParameters();
		MaxMatrixDimension = 20000;
		MergeWithTour = Recombination == IPT ? (&LKH::LKHAlg::MergeWithTourIPT) :
			(&LKH::LKHAlg::MergeWithTourGPX2);

		std::string filepath = dataFile;
		ReadProblem(filepath);
		if (SubproblemSize > 0) {
			if (DelaunayPartitioning)
				SolveDelaunaySubproblems();
			else if (KarpPartitioning)
				SolveKarpSubproblems();
			else if (KCenterPartitioning)
				SolveKCenterSubproblems();
			else if (KMeansPartitioning)
				SolveKMeansSubproblems();
			else if (RohePartitioning)
				SolveRoheSubproblems();
			else if (MoorePartitioning || SierpinskiPartitioning)
				SolveSFCSubproblems();
			else
				SolveTourSegmentSubproblems();
			return EXIT_SUCCESS;
		}
		AllocateStructures();
		if (ProblemType == TSPTW)
			TSPTW_Reduce();
		if (ProblemType == VRPB || ProblemType == VRPBTW)
			VRPB_Reduce();
		if (ProblemType == PDPTW)
			PDPTW_Reduce();
		CreateCandidateSet();
		InitializeStatistics();

		if (Norm != 0 || Penalty) {
			Norm = 9999;
			BestCost = std::numeric_limits<GainType>::max();
			BestPenalty = CurrentPenalty = std::numeric_limits<GainType>::max();
		}
		else {
			/* The ascent has solved the problem! */
			Optimum = BestCost = (GainType)LowerBound;
			UpdateStatistics(Optimum, GetTime() - LastTime);
			RecordBetterTour();
			RecordBestTour();
			CurrentPenalty = std::numeric_limits<GainType>::max();
			BestPenalty = CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
			//WriteTour(OutputTourFileName, BestTour, BestCost);
			//WriteTour(TourFileName, BestTour, BestCost);
			Runs = 0;
		}

		/* Find a specified number (Runs) of local optima */

		for (Run = 1; Run <= Runs; Run++) {
			LastTime = GetTime();
			Cost = FindTour();      /* using the Lin-Kernighan heuristic */

			if (MaxPopulationSize > 1 && !TSPTW_Makespan) {
				/* Genetic algorithm */
				int i;
				for (i = 0; i < PopulationSize; i++) {
					GainType OldPenalty = CurrentPenalty;
					GainType OldCost = Cost;
					Cost = MergeTourWithIndividual(i, this);
					if (TraceLevel >= 1 &&
						(CurrentPenalty < OldPenalty ||
							(CurrentPenalty == OldPenalty && Cost < OldCost))) {
						/*if (CurrentPenalty)
							printff("  Merged with %d: Cost = " GainFormat,
								i + 1, Cost);
						else
							printff("  Merged with %d: Cost = " GainFormat "_"
								GainFormat, i + 1, CurrentPenalty, Cost);*/
								/*if (Optimum != MINUS_INFINITY && Optimum != 0) {
									if (ProblemType != CCVRP && ProblemType != TRP &&
										ProblemType != MLP &&
										MTSPObjective != MINMAX &&
										MTSPObjective != MINMAX_SIZE)
										printff(", Gap = %0.4f%%",
											100.0 * (Cost - Optimum) / Optimum);
									else
										printff(", Gap = %0.4f%%",
											100.0 * (CurrentPenalty - Optimum) /
											Optimum);
								}
								printff("\n");*/
					}
				}
				if (!HasFitness(CurrentPenalty, Cost)) {
					if (PopulationSize < MaxPopulationSize) {
						AddToPopulation(CurrentPenalty, Cost, this);
						if (TraceLevel >= 1)
							PrintPopulation(this);
					}
					else if (SmallerFitness(CurrentPenalty, Cost,
						PopulationSize - 1)) {
						i = ReplacementIndividual(CurrentPenalty, Cost, this);
						ReplaceIndividualWithTour(i, CurrentPenalty, Cost, this);
						if (TraceLevel >= 1)
							PrintPopulation(this);
					}
				}
			}
			else if (Run > 1 && !TSPTW_Makespan)
				Cost = MergeTourWithBestTour();
			if (CurrentPenalty < BestPenalty ||
				(CurrentPenalty == BestPenalty && Cost < BestCost)) {
				BestPenalty = CurrentPenalty;
				BestCost = Cost;
				RecordBetterTour();
				RecordBestTour();
				//WriteTour(TourFileName, BestTour, BestCost);
			}
			OldOptimum = Optimum;
			if (!Penalty ||
				(MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE)) {
				if (CurrentPenalty == 0 && Cost < Optimum)
					Optimum = Cost;
			}
			else if (CurrentPenalty < Optimum)
				Optimum = CurrentPenalty;
			if (Optimum < OldOptimum) {
				printff("*** New OPTIMUM = " GainFormat " ***\n\n", Optimum);
				if (FirstNode->InputSuc) {
					Node* N = FirstNode;
					while ((N = N->InputSuc = N->Suc) != FirstNode);
				}
			}
			Time = fabs(GetTime() - LastTime);
			UpdateStatistics(Cost, Time);
			/*if (TraceLevel >= 1 && Cost != std::numeric_limits<GainType>::max()) {
				printff("Run %d: ", Run);
				StatusReport(Cost, LastTime, "");
				printff("\n");
			}*/
			if (StopAtOptimum && MaxPopulationSize >= 1) {
				if (ProblemType != CCVRP && ProblemType != TRP &&
					ProblemType != MLP &&
					MTSPObjective != MINMAX &&
					MTSPObjective != MINMAX_SIZE ?
					CurrentPenalty == 0 && Cost == Optimum :
					CurrentPenalty == Optimum) {
					Runs = Run;
					break;
				}
			}
			if (PopulationSize >= 2 &&
				(PopulationSize == MaxPopulationSize ||
					Run >= 2 * MaxPopulationSize) && Run < Runs) {
				Node* N;
				int Parent1, Parent2;
				Parent1 = LinearSelection(PopulationSize, 1.25, this);
				do
					Parent2 = LinearSelection(PopulationSize, 1.25, this);
				while (Parent2 == Parent1);
				ApplyCrossover(Parent1, Parent2, this);
				N = FirstNode;
				do {
					if (ProblemType != HCP && ProblemType != HPP) {
						int d = (this->*C)(N, N->Suc);
						AddCandidate(N, N->Suc, d, std::numeric_limits<int>::max());
						AddCandidate(N->Suc, N, d, std::numeric_limits<int>::max());
					}
					N = N->InitialSuc = N->Suc;
				} while (N != FirstNode);
			}
			SRandom(++Seed);
		}
		//get tour
		int Forward = 0;
		int n = DimensionSaved;
		for (i = 1; i < n && BestTour[i] != 1; i++);
		Forward = Asymmetric ||
			BestTour[i < n ? i + 1 : 1] < BestTour[i > 1 ? i - 1 : Dimension];
		for (int j = 1; j <= n; j++) {
			if (BestTour[i] <= n)
				tour.push_back(BestTour[i]);
			if (Forward) {
				if (++i > n)
					i = 1;
			}
			else if (--i < 1)
				i = n;
		}
		//WriteTour(TourFileName, BestTour, BestCost);
		//PrintStatistics();
		freeAll();


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (Salesmen > 1) {
			if (Dimension == DimensionSaved) {
				for (i = 1; i <= Dimension; i++) {
					N = &NodeSet[BestTour[i - 1]];
					(N->Suc = &NodeSet[BestTour[i]])->Pred = N;
				}
			}
			else {
				for (i = 1; i <= DimensionSaved; i++) {
					Node* N1 = &NodeSet[BestTour[i - 1]];
					Node* N2 = &NodeSet[BestTour[i]];
					Node* M1 = &NodeSet[N1->Id + DimensionSaved];
					Node* M2 = &NodeSet[N2->Id + DimensionSaved];
					(M1->Suc = N1)->Pred = M1;
					(N1->Suc = M2)->Pred = N1;
					(M2->Suc = N2)->Pred = M2;
				}
			}
			CurrentPenalty = BestPenalty;
			MTSP_Report(BestPenalty, BestCost);
			MTSP_WriteSolution(MTSPSolutionFileName, BestPenalty, BestCost);
			SINTEF_WriteSolution(SINTEFSolutionFileName, BestCost);
		}
		if (ProblemType == ACVRP ||
			ProblemType == BWTSP ||
			ProblemType == CCVRP ||
			ProblemType == CTSP ||
			ProblemType == CVRP ||
			ProblemType == CVRPTW ||
			ProblemType == MLP ||
			ProblemType == M_PDTSP ||
			ProblemType == M1_PDTSP ||
			MTSPObjective != -1 ||
			ProblemType == ONE_PDTSP ||
			ProblemType == OVRP ||
			ProblemType == PDTSP ||
			ProblemType == PDTSPL ||
			ProblemType == PDPTW ||
			ProblemType == RCTVRP ||
			ProblemType == RCTVRPTW ||
			ProblemType == SOP ||
			ProblemType == TRP ||
			ProblemType == TSPTW ||
			ProblemType == VRPB ||
			ProblemType == VRPBTW || ProblemType == VRPPD) {
			printff("Best %s solution:\n", Type);
			CurrentPenalty = BestPenalty;
			SOP_Report(BestCost);
		}
		//printff("\n");
		//system("pause");

		delete ParameterFileName;
		ParameterFileName = nullptr;
		return EXIT_SUCCESS;
	}



	int LKHAlg::run(const string& fileName, vector<size_t>& tour)
	{
		lkh_ptr = this;
		GainType Cost, OldOptimum;
		double Time, LastTime = GetTime();
		Node* N;
		int i;

		/* Read the specification of the problem */
	   /* if (argc >= 2)
			ParameterFileName = argv[1];*/
			//for (int i = 0; i < 5; ++i) {
			//	char para_name[100] = "./DVRP/fcm_twl/cluster";
			//	//para_name =;
			//	char str[1] = { 0 };
			//	itoa(i, str, 10);
			//	strcat(para_name, str);
			//	strcat(para_name, ".par");
		ParameterFileName = new char[fileName.size() + 1];
		strcpy(ParameterFileName, fileName.c_str());

		//ParameterFileName = "./DVRP/shortest_path/ft70.par";
		ReadParameters();
		//setDefaultParaments();
		MaxMatrixDimension = 20000;
		MergeWithTour = Recombination == IPT ? (&LKH::LKHAlg::MergeWithTourIPT) :
			(&LKH::LKHAlg::MergeWithTourGPX2);

		std::string filepath = static_cast<string>(ofec::g_working_directory) + "/instance/algorithm/realworld/DVRP/LKH/temp file/" + ProblemFileName;
		ReadProblem(filepath);
		if (SubproblemSize > 0) {
			if (DelaunayPartitioning)
				SolveDelaunaySubproblems();
			else if (KarpPartitioning)
				SolveKarpSubproblems();
			else if (KCenterPartitioning)
				SolveKCenterSubproblems();
			else if (KMeansPartitioning)
				SolveKMeansSubproblems();
			else if (RohePartitioning)
				SolveRoheSubproblems();
			else if (MoorePartitioning || SierpinskiPartitioning)
				SolveSFCSubproblems();
			else
				SolveTourSegmentSubproblems();
			return EXIT_SUCCESS;
		}
		AllocateStructures();
		if (ProblemType == TSPTW)
			TSPTW_Reduce();
		if (ProblemType == VRPB || ProblemType == VRPBTW)
			VRPB_Reduce();
		if (ProblemType == PDPTW)
			PDPTW_Reduce();
		CreateCandidateSet();
		InitializeStatistics();

		if (Norm != 0 || Penalty) {
			Norm = 9999;
			BestCost = std::numeric_limits<GainType>::max();
			BestPenalty = CurrentPenalty = std::numeric_limits<GainType>::max();
		}
		else {
			/* The ascent has solved the problem! */
			Optimum = BestCost = (GainType)LowerBound;
			UpdateStatistics(Optimum, GetTime() - LastTime);
			RecordBetterTour();
			RecordBestTour();
			CurrentPenalty = std::numeric_limits<GainType>::max();
			BestPenalty = CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
			//WriteTour(OutputTourFileName, BestTour, BestCost);
			//WriteTour(TourFileName, BestTour, BestCost);
			Runs = 0;
		}

		/* Find a specified number (Runs) of local optima */

		for (Run = 1; Run <= Runs; Run++) {
			LastTime = GetTime();
			Cost = FindTour();      /* using the Lin-Kernighan heuristic */
			if (MaxPopulationSize > 1 && !TSPTW_Makespan) {
				/* Genetic algorithm */
				int i;
				for (i = 0; i < PopulationSize; i++) {
					GainType OldPenalty = CurrentPenalty;
					GainType OldCost = Cost;
					Cost = MergeTourWithIndividual(i, this);
					if (TraceLevel >= 1 &&
						(CurrentPenalty < OldPenalty ||
							(CurrentPenalty == OldPenalty && Cost < OldCost))) {
						/*if (CurrentPenalty)
							printff("  Merged with %d: Cost = " GainFormat,
								i + 1, Cost);
						else
							printff("  Merged with %d: Cost = " GainFormat "_"
								GainFormat, i + 1, CurrentPenalty, Cost);*/
								/*if (Optimum != MINUS_INFINITY && Optimum != 0) {
									if (ProblemType != CCVRP && ProblemType != TRP &&
										ProblemType != MLP &&
										MTSPObjective != MINMAX &&
										MTSPObjective != MINMAX_SIZE)
										printff(", Gap = %0.4f%%",
											100.0 * (Cost - Optimum) / Optimum);
									else
										printff(", Gap = %0.4f%%",
											100.0 * (CurrentPenalty - Optimum) /
											Optimum);
								}
								printff("\n");*/
					}
				}
				if (!HasFitness(CurrentPenalty, Cost)) {
					if (PopulationSize < MaxPopulationSize) {
						AddToPopulation(CurrentPenalty, Cost, this);
						if (TraceLevel >= 1)
							PrintPopulation(this);
					}
					else if (SmallerFitness(CurrentPenalty, Cost,
						PopulationSize - 1)) {
						i = ReplacementIndividual(CurrentPenalty, Cost, this);
						ReplaceIndividualWithTour(i, CurrentPenalty, Cost, this);
						if (TraceLevel >= 1)
							PrintPopulation(this);
					}
				}
			}
			else if (Run > 1 && !TSPTW_Makespan)
				Cost = MergeTourWithBestTour();
			if (CurrentPenalty < BestPenalty ||
				(CurrentPenalty == BestPenalty && Cost < BestCost)) {
				BestPenalty = CurrentPenalty;
				BestCost = Cost;
				RecordBetterTour();
				RecordBestTour();
				//WriteTour(TourFileName, BestTour, BestCost);
			}
			OldOptimum = Optimum;
			if (!Penalty ||
				(MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE)) {
				if (CurrentPenalty == 0 && Cost < Optimum)
					Optimum = Cost;
			}
			else if (CurrentPenalty < Optimum)
				Optimum = CurrentPenalty;
			if (Optimum < OldOptimum) {
				//printff("*** New OPTIMUM = " GainFormat " ***\n\n", Optimum);
				if (FirstNode->InputSuc) {
					Node* N = FirstNode;
					while ((N = N->InputSuc = N->Suc) != FirstNode);
				}
			}
			Time = fabs(GetTime() - LastTime);
			UpdateStatistics(Cost, Time);
			/*if (TraceLevel >= 1 && Cost != std::numeric_limits<GainType>::max()) {
				printff("Run %d: ", Run);
				StatusReport(Cost, LastTime, "");
				printff("\n");
			}*/
			if (StopAtOptimum && MaxPopulationSize >= 1) {
				if (ProblemType != CCVRP && ProblemType != TRP &&
					ProblemType != MLP &&
					MTSPObjective != MINMAX &&
					MTSPObjective != MINMAX_SIZE ?
					CurrentPenalty == 0 && Cost == Optimum :
					CurrentPenalty == Optimum) {
					Runs = Run;
					break;
				}
			}
			if (PopulationSize >= 2 &&
				(PopulationSize == MaxPopulationSize ||
					Run >= 2 * MaxPopulationSize) && Run < Runs) {
				Node* N;
				int Parent1, Parent2;
				Parent1 = LinearSelection(PopulationSize, 1.25, this);
				do
					Parent2 = LinearSelection(PopulationSize, 1.25, this);
				while (Parent2 == Parent1);
				ApplyCrossover(Parent1, Parent2, this);
				N = FirstNode;
				do {
					if (ProblemType != HCP && ProblemType != HPP) {
						int d = (this->*C)(N, N->Suc);
						AddCandidate(N, N->Suc, d, std::numeric_limits<int>::max());
						AddCandidate(N->Suc, N, d, std::numeric_limits<int>::max());
					}
					N = N->InitialSuc = N->Suc;
				} while (N != FirstNode);
			}
			SRandom(++Seed);
		}
		//get tour
		int Forward = 0;
		int n = DimensionSaved;
		for (i = 1; i < n && BestTour[i] != 1; i++);
		Forward = Asymmetric ||
			BestTour[i < n ? i + 1 : 1] < BestTour[i > 1 ? i - 1 : Dimension];
		for (int j = 1; j <= n; j++) {
			if (BestTour[i] <= n)
				tour.push_back(BestTour[i]);
			if (Forward) {
				if (++i > n)
					i = 1;
			}
			else if (--i < 1)
				i = n;
		}
		//WriteTour(TourFileName, BestTour, BestCost);
		//PrintStatistics();
		freeAll();


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (Salesmen > 1) {
			if (Dimension == DimensionSaved) {
				for (i = 1; i <= Dimension; i++) {
					N = &NodeSet[BestTour[i - 1]];
					(N->Suc = &NodeSet[BestTour[i]])->Pred = N;
				}
			}
			else {
				for (i = 1; i <= DimensionSaved; i++) {
					Node* N1 = &NodeSet[BestTour[i - 1]];
					Node* N2 = &NodeSet[BestTour[i]];
					Node* M1 = &NodeSet[N1->Id + DimensionSaved];
					Node* M2 = &NodeSet[N2->Id + DimensionSaved];
					(M1->Suc = N1)->Pred = M1;
					(N1->Suc = M2)->Pred = N1;
					(M2->Suc = N2)->Pred = M2;
				}
			}
			CurrentPenalty = BestPenalty;
			MTSP_Report(BestPenalty, BestCost);
			MTSP_WriteSolution(MTSPSolutionFileName, BestPenalty, BestCost);
			SINTEF_WriteSolution(SINTEFSolutionFileName, BestCost);
		}
		if (ProblemType == ACVRP ||
			ProblemType == BWTSP ||
			ProblemType == CCVRP ||
			ProblemType == CTSP ||
			ProblemType == CVRP ||
			ProblemType == CVRPTW ||
			ProblemType == MLP ||
			ProblemType == M_PDTSP ||
			ProblemType == M1_PDTSP ||
			MTSPObjective != -1 ||
			ProblemType == ONE_PDTSP ||
			ProblemType == OVRP ||
			ProblemType == PDTSP ||
			ProblemType == PDTSPL ||
			ProblemType == PDPTW ||
			ProblemType == RCTVRP ||
			ProblemType == RCTVRPTW ||
			ProblemType == SOP ||
			ProblemType == TRP ||
			ProblemType == TSPTW ||
			ProblemType == VRPB ||
			ProblemType == VRPBTW || ProblemType == VRPPD) {
			printff("Best %s solution:\n", Type);
			CurrentPenalty = BestPenalty;
			SOP_Report(BestCost);
		}
		//printff("\n");
		//system("pause");
		return EXIT_SUCCESS;
	}
	void LKH::LKHAlg::freeAll()
	{
		freeERXT();
		freeGain23();
		//	freeSequence();
		freeCreateDelaunay();
		freeCreateQuadrant();
		freeGreedyTour();
		freePatchCycle();
		freeReadPenalties();
	}
}