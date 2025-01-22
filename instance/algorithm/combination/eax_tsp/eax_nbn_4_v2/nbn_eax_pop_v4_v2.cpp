#include "nbn_eax_pop_v4_v2.h"
#include "../../../../../utility/nbn_visualization/nbn_fla.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"



int ofec::PopNBN_EAX_V4_V2::evolveWithHnsw(Problem* pro,
	Algorithm* alg, Random* rnd)
{

	++m_curIter;
	
	
	// generete solutions
	std::vector<int> shuffleIds(m_curPop.size());
	for (int idx(0); idx < shuffleIds.size(); ++idx) {
		shuffleIds[idx] = idx;
	}
	rnd->uniform.shuffle(shuffleIds.begin(), shuffleIds.end());
	std::vector<std::shared_ptr<Indi>> sonSols(m_curPop.size());

	for (int idx(0); idx < sonSols.size(); ++idx) {
		//m_id2inds.push_back(NBN_indi());
		sonSols[idx].reset(new Indi);
		sonSols[idx]->initialize(pro->numberVariables(), idx, m_curIter);
		sonSols[idx]->copySol(*m_curPop[idx]);

	}

	std::vector<Indi*> newsols(sonSols.size(),nullptr);;
	std::vector<bool> activeSons(sonSols.size(), false);
	//std::vector<int> cursolx;
	for (int idx(0); idx < sonSols.size(); ++idx) {
		tCross->setParents(*sonSols[idx], *m_curPop[shuffleIds[idx]], fFlagC, Nch);
		tCross->doIt(*sonSols[idx], *m_curPop[shuffleIds[idx]], Nch, 1, fFlagC, fEdgeFreq);
		//updateFitness(sonSols[idx].get());
		sonSols[idx]->updateSolBase();
		updateFitness(sonSols[idx].get());
		auto& cursol = sonSols[idx];
	//	cursol->transferSol(cursolx);
		Indi* indi = nullptr;
		bool flag = addSolToHistory(cursol, indi);
		if (flag || indi->m_curClassFlag) {
			activeSons[idx] = true;
			newsols[idx] = indi;
			indi->m_stagnation_time = 0;
			//auto& neighbor = m_hnswModel.getNeighbor(indi->m_solId);
			SolutionBase* parent = nullptr;
			getParentTwoNeighbor(indi, parent, pro);
			if (parent != nullptr) {
				auto solPar = dynamic_cast<Indi*>(parent);
				indi->setCurClass(solPar->m_curClassFlag);
			}
			else indi->setCurClass(m_curEvolveClass);

		}


	//	auto hash = m_solHash.calHash(cursolx);
	//	auto& indi = m_solMap[hash];
	//	if (indi == nullptr) {
	//		m_indis.emplace_back(std::move(cursol));
	//		m_indis.back()->setRndId(hash);
	//	//	m_indis.back()->setCurClass(0);
	//		indi = m_indis.back().get();
	////		m_curPop.push_back(indi);
	//		m_indis.back()->setSolId(m_hnswModel.AddData(indi));
	//		//newsols.push_back(indi);
	//	//	int curId= m_indis
	//	//	newsols[idx] = indi;
	//		activeSons[idx] = true;
	//		newsols[idx] = indi;

	//		//auto& neighbor = m_hnswModel.getNeighbor(indi->m_solId);
	//		SolBase* parent = nullptr;
	//		getParentTwoNeighbor(indi, parent, pro);
	//		if (parent != nullptr) {
	//			auto solPar = dynamic_cast<Indi*>(parent);
	//			m_indis.back()->setCurClass(solPar->m_curClassFlag);
	//		}

	//		//double maxDis = std::numeric_limits<double>::max();
	//		//for (auto& it : neighbor) {
	//		//	if (it.GetNode()->getOriginData()->fitness() > indi->fitness()) {
	//		//		if (maxDis < it.GetDistance()) {
	//		//			maxDis = it.GetDistance();
	//		//			parent = it.GetNode()->getOriginData();
	//		//		}
	//		//	}
	//		//}
	//		//if (parent == nullptr) {
	//		//	for (auto& it : neighbor) {
	//		//		auto& neighbor2 = m_hnswModel.getNeighbor(it.GetNode()->GetId());
	//		//		for (auto& it2 : neighbor2) {
	//		//			if (it2.GetNode()->getOriginData()->fitness() > indi->fitness()) {
	//		//				double d = it2.GetNode()->getOriginData()->variableDistance(*indi, pro);

	//		//				if (maxDis < d) {
	//		//					maxDis = d;
	//		//					parent = it.GetNode()->getOriginData();
	//		//				}
	//		//			}
	//		//		}
	//		//	}
	//		//}
	//	}
	}
	if (m_expandStage && m_curPop.size() < m_numPop) {
		std::vector<Indi*> leftSols;
		//m_expandStage = false;
		// tabu search 
		for (int idx(0); idx < newsols.size(); ++idx) {
			if (activeSons[idx]) {
				if (newsols[idx]->m_curClassFlag == 0 || sonSols[idx]->m_curClassFlag == m_curEvolveClass) {
					newsols[idx]->m_curClassFlag = m_curEvolveClass;
					if (newsols[idx]->fitness() > m_curPop[idx]->fitness()) {
						m_curPop[idx] = newsols[idx];
					}
					else {
						leftSols.push_back(newsols[idx]);
					}
				}
				//	auto& neighbor = m_hnswModel.getNeighbor(newsols[idx])
			}
		}

		rnd->uniform.shuffle(leftSols.begin(), leftSols.end());
		for (auto& it : leftSols) {
			m_curPop.push_back(it);
		}
		m_curPop.resize(m_numPop);
	}
	else {
		m_expandStage = false;
		// tabu search 
		for (int idx(0); idx < newsols.size(); ++idx) {
			if (activeSons[idx]) {
				if (newsols[idx]->m_curClassFlag == 0 || newsols[idx]->m_curClassFlag == m_curEvolveClass) {
					newsols[idx]->m_curClassFlag = m_curEvolveClass;
					if (newsols[idx]->fitness() > m_curPop[idx]->fitness()) {
						m_curPop[idx] = newsols[idx];
					}
				}
				//	auto& neighbor = m_hnswModel.getNeighbor(newsols[idx])
			}
		}

	}




	newsols.clear();
	//int maxsize = std::min<int>(m_maxStagnationTime, m_curPop.size() * (m_curPop.size() - 1) / 2.0);
	for (int idx(0); idx < m_curPop.size(); ++idx) {
		if (++m_curPop[idx]->m_stagnation_time <= m_maxStagnationTime) {
			newsols.push_back(m_curPop[idx]);
		}
	}
	
	swap(newsols, m_curPop);

	if (m_curPop.size() <= 1) {
		m_curPop.clear();

		++m_curEvolveClass;
		m_expandStage = true;
		// calculate nbn

		std::vector<int> belong;
		std::vector<double> dis2parent;
		std::vector<SolutionBase*> sols;
		std::vector<double> vFitness;
		for (auto& it : m_indis) {
			sols.push_back(it.get());
			vFitness.push_back(it->fitness());
		}
		std::vector<std::vector<int>> neighbors(m_indis.size());
		for (int idx(0); idx < m_indis.size(); ++idx) {
			auto& nei = m_hnswModel.getNeighbor(idx);
			auto& curnei = neighbors[idx];
			for (auto& neiId : nei) {
				curnei.push_back(neiId.GetNode()->GetId());
			}
		}


		NBN_NearestBetterCalculator::calculate(sols, vFitness, neighbors,
			belong, dis2parent,
			pro, rnd);
		m_dis2parent = dis2parent;
		std::vector<int> solIds(sols.size());
		for (int idx(0); idx < solIds.size(); ++idx) {
			solIds[idx] = idx;
		}

		std::sort(solIds.begin(), solIds.end(), [&](int a,int b) {
			return vFitness[a] > vFitness[b];
		});
		
		int selectedPeak = -1;

		for (auto& solId : solIds) {
			if (judgePeak(m_indis[solId].get(), dis2parent[solId])) {
				selectedPeak = solId;
				break;
			}
			//auto& cursol = m_indis[solId];
			//if (cursol->m_explorationTimes == 0) {
			//	
			//}
		}


		//std::vector<int> peaks;
		//for (int idx(0); idx < m_indis.size(); ++idx) {
		//	auto& nei = m_hnswModel.getNeighbor(idx);
		//	bool flag = true;
		//	auto curdis = dis2parent[idx];
		//	for (auto& neiId : nei) {
		//		if (curdis < neiId.GetDistance()) {
		//			flag = false;
		//			break;
		//		}
		//	}
		//	if (flag) {
		//		peaks.push_back(idx);
		//	}
		//}

		//std::sort(peaks.begin(), peaks.end(), [&](int a,int b) {
		//	return vFitness[a] > vFitness[b];
		//});


		//for (auto& it : peaks) {
		//	auto cursol = m_indis[it].get();
		//	if (cursol->m_explorationTimes == 0) {
		//		selectedPeak = it;
		//		++cursol->m_explorationTimes;
		//	}
		//}

		if (selectedPeak == -1) {

			for (auto& it : m_indis) {
				it->m_curClassFlag = it->m_searchArea;
			}


			generateSolRandom(pro);
		}
		else {
			std::vector<int> curbasin;
			NBN_FLA_info::getBasin(belong, selectedPeak, curbasin);
			for (auto& it : curbasin) {
				auto& cursol = m_indis[it];
				cursol->setCurClass(m_curEvolveClass);
			}
			std::vector<int> filterBasin;
			std::vector<double> basinDis;
			auto peakSol = m_indis[selectedPeak].get();
			m_curPeak = peakSol;
			for (auto& it : curbasin) {
				auto cursol = m_indis[it].get();
				double dis = pro->variableDistance(*peakSol, *cursol);
				if (dis <= m_searchRadius) {
					filterBasin.push_back(it);
					basinDis.push_back(dis);
				}
			}
			for (auto& it : filterBasin) {
				m_indis[it]->m_searchArea = m_curEvolveClass;
			}
			double basinRadius = 0;
			ofec::calMax(basinDis, basinRadius);

			rnd->uniform.shuffle(filterBasin.begin(), filterBasin.end());
			while (m_curPop.size() < m_numPop&&!filterBasin.empty()) {
				m_curPop.push_back(m_indis[filterBasin.back()].get());
				filterBasin.pop_back();
			}
			
			int expandSearchTimes = 1e4;

			std::shared_ptr<Indi> cursol;
			int step = m_searchRadius / 2.0;

			while (m_curPop.size() < m_numPop&&--expandSearchTimes) {

				cursol.reset(new Indi);
				cursol->initialize(pro->numberVariables(), -1, m_curIter);
				cursol->copy(*peakSol);

				int curstep = rnd->uniform.nextNonStd<int>(1, step + 1);

				TravellingSalesman::mutation(*cursol, curstep, rnd);
				cursol->evaluate(pro, alg);
				//updateFitness(cursol.get());
				cursol->updateLinkInfo();
				tKopt->doIt(*cursol);
				Indi* indi = nullptr;
				bool flag = addSolToHistory(cursol, indi);
				if (flag || indi->m_curClassFlag == 0) {
					SolutionBase* parent = nullptr;
					getParentTwoNeighbor(indi, parent, pro);
					if (parent != nullptr) {
						auto solPar = dynamic_cast<Indi*>(parent);
						indi->setCurClass(solPar->m_curClassFlag);
						if (indi->m_curClassFlag == m_curEvolveClass || indi->m_curClassFlag == 0) {
							m_curPop.push_back(indi);
							indi->m_searchArea = m_curEvolveClass;
							indi->m_stagnation_time = 0;
						}
					}
					else {
						indi->setCurClass(m_curEvolveClass);
						m_curPop.push_back(indi);
						indi->m_searchArea = m_curEvolveClass;
						indi->m_stagnation_time = 0;
					}
				}
			}
			
		}
	}
	

	return 0;
}



void ofec::PopNBN_EAX_V4_V2::generateSolRandom(ofec::Problem* pro) {
	std::shared_ptr<Indi> cursol;
	//std::vector<int> cursolx;
	while (m_curPop.size() < m_numPop) {
		cursol.reset(new Indi());
		cursol->initialize(pro->numberVariables(), m_curPop.size(), m_curIter);

		tKopt->makeRandSol(*cursol); /* Make a random tour */
		tKopt->doIt(*cursol);        /* Apply the local search with the 2-opt neighborhood */
		updateFitness(cursol.get());
		cursol->updateSolBase();
		Indi* indi = nullptr;
		bool flag = addSolToHistory(cursol, indi);
		if (flag || indi->m_curClassFlag == 0) {
			//indi->m_stagnation_time = 0;
			SolutionBase* parent = nullptr;
			getParentTwoNeighbor(indi, parent, pro);
			if (parent != nullptr) {
				auto solPar = dynamic_cast<Indi*>(parent);
				indi->setCurClass(solPar->m_curClassFlag);
				if (indi->m_curClassFlag == m_curEvolveClass || indi->m_curClassFlag == 0) {
					m_curPop.push_back(indi);
					indi->m_stagnation_time = 0;
				}
			}
			else {
				indi->setCurClass(m_curEvolveClass);
				m_curPop.push_back(indi);
				indi->m_stagnation_time = 0;
			}
		}

		//const auto& cursolx = cursol->variable().vect();
		////cursol->transferSol(cursolx);
		//auto hash = m_solHash.calHash(cursolx);
		//auto& indi = m_solMap[hash];
		//if (indi == nullptr) {
		//	m_indis.emplace_back(std::move(cursol));
		//	m_indis.back()->setSolId(m_indis.size() - 1);
		//	m_indis.back()->setRndId(hash);
		////	m_indis.back()->setCurClass(m_curEvolveClass);
		//	indi = m_indis.back().get();
		//	
		//	m_hnswModel.AddData(indi);
		//	SolBase* parent = nullptr;
		//	getParentTwoNeighbor(indi, parent, pro);
		//	if (parent != nullptr) {
		//		auto solPar = dynamic_cast<Indi*>(parent);
		//		indi->setCurClass(solPar->m_curClassFlag);
		//		if (indi->m_curClassFlag == m_curEvolveClass || indi->m_curClassFlag == 0) {
		//			m_curPop.push_back(indi);
		//		}
		//	}
		//}

	}
	std::vector< std::vector<std::vector<int>>*> links;
	for (auto& it : m_curPop) {
		links.push_back(&it->fLink);
	}
	calEdgeFreq(links);
//	m_expandStage = true;
}

void ofec::PopNBN_EAX_V4_V2::calNodeState(Indi* indi) {
	auto& neighbor = m_hnswModel.getNeighbor(indi->m_solId);
	bool flag = true;
	for (auto& it : neighbor) {
		if (it.GetNode()->getOriginData()->fitness() > indi->fitness()) {
			flag = false;
			break;
		}
	}

	if (flag) {
		indi->m_state = Indi::NodeState::kPossiblePeak;
	}
}


void ofec::PopNBN_EAX_V4_V2::getParentTwoNeighbor(Indi* indi, SolutionBase*& parent, ofec::Problem* pro) {
	auto& neighbor = m_hnswModel.getNeighbor(indi->m_solId);
	parent = nullptr;
	double maxDis = std::numeric_limits<double>::max();
	for (auto& it : neighbor) {
		if (it.GetNode()->getOriginData()->fitness() > indi->fitness()) {
			if (maxDis < it.GetDistance()) {
				maxDis = it.GetDistance();
				parent = it.GetNode()->getOriginData();
			}
		}
	}
	if (parent == nullptr) {
		for (auto& it : neighbor) {
			auto& neighbor2 = m_hnswModel.getNeighbor(it.GetNode()->GetId());
			for (auto& it2 : neighbor2) {
				if (it2.GetNode()->getOriginData()->fitness() > indi->fitness()) {
					double d = it2.GetNode()->getOriginData()->variableDistance(*indi, pro);

					if (maxDis < d) {
						maxDis = d;
						parent = it.GetNode()->getOriginData();
					}
				}
			}
		}
	}
}

void ofec::PopNBN_EAX_V4_V2::calNBN(
	std::vector<TreeGraphSimple::NBNdata>& nbn_data, 
	std::vector<TreeGraphSimple::NBNdata>& nbn_data_best, 
	std::vector<int>& curPop, 
	std::vector<int>& curBasin, 
	std::vector<int>& colorfulArea,
	std::vector<int>& bestIds, ofec::Problem* pro, ofec::Random* rnd)
{

	std::vector<int> belong;
	std::vector<double> dis2parent;
	std::vector<double> vFitness;
	std::vector<SolutionBase*> sols;

	for (auto& it : m_indis) {
		sols.push_back(it.get());
		vFitness.push_back(it->fitness());
	}
	std::vector<std::vector<int>> neighbors(m_indis.size());
	for (int idx(0); idx < m_indis.size(); ++idx) {
		auto& nei = m_hnswModel.getNeighbor(idx);
		auto& curnei = neighbors[idx];
		for (auto& neiId : nei) {
			curnei.push_back(neiId.GetNode()->GetId());
		}
	}


	NBN_NearestBetterCalculator::calculate(sols, vFitness, neighbors,
		belong, dis2parent,
		pro, rnd);


	std::vector<int> solIds(sols.size());
	for (int idx(0); idx < solIds.size(); ++idx) {
		solIds[idx] = idx;
	}
	
	transferNBNData(nbn_data, solIds, belong, dis2parent, vFitness);

	curPop.resize(m_curPop.size());
	for (int idx(0); idx < m_curPop.size(); ++idx) {
		curPop[idx] = m_curPop[idx]->m_solId;
	}
	curBasin.clear();

	for (auto& it : m_indis) {
		if (it->m_curClassFlag == m_curEvolveClass) {
			curBasin.push_back(it->m_solId);
		}
	}
	{
		int colorId(0);
		std::map<unsigned long long, int> mii;
		mii[0] = colorId++;
		for (auto& it : m_indis) {
			if (mii.find(it->m_curClassFlag) == mii.end()) {
				mii[it->m_curClassFlag] = colorId++;
			}
		}
		colorfulArea.resize(m_indis.size());
		for (auto& it : m_indis) {
			colorfulArea[it->m_solId] = mii[it->m_curClassFlag];
		}
	}
	
	nbn_data_best = nbn_data;
	
	{

		bestIds.clear();
		nbn_data_best.resize(nbn_data.size() + m_bestSols.size());
		for (int idx(nbn_data.size()); idx < nbn_data_best.size(); ++idx) {
			nbn_data_best[idx].m_belong = idx;
			nbn_data_best[idx].m_dis2parent = std::numeric_limits<double>::max();
			nbn_data_best[idx].m_fitness = m_bestSols[idx - nbn_data.size()].fitness();
			bestIds.push_back(idx);
		}

		std::vector<SolutionBase*> totalSols;
		for (auto& it : m_indis) {
			totalSols.push_back(it.get());
		}
		for (auto& it : m_bestSols) {
			totalSols.push_back(&it);
		}

		for (int idx(nbn_data.size()); idx < nbn_data_best.size(); ++idx) {
			for (int idy(0); idy < idx; ++idy) {
				if (nbn_data_best[idx].m_fitness > nbn_data_best[idy].m_fitness) {

					double dis = totalSols[idx]->variableDistance(*totalSols[idy], pro);
					if (dis < nbn_data_best[idy].m_dis2parent) {
						nbn_data_best[idy].m_dis2parent = dis;
						nbn_data_best[idy].m_belong = idx;
					}
				}
				else if (nbn_data_best[idx].m_fitness < nbn_data_best[idy].m_fitness) {
					double dis = totalSols[idx]->variableDistance(*totalSols[idy], pro);

					if (dis < nbn_data_best[idx].m_dis2parent) {
						nbn_data_best[idx].m_dis2parent = dis;
						nbn_data_best[idx].m_belong = idy;
					}
				}
			}
		}
	}

}

void ofec::PopNBN_EAX_V4_V2::initialize(Problem* pro, Random* rnd)
{
//	PopNBN_EAX_V4::initialize(pro, rnd);
	m_hnswModel.initialize(pro, rnd);
	using namespace ofec;
	using namespace eax_tsp;

	m_curIter = 0;
	Npop = 100;
	Nch = 30;
	define(pro, rnd->getSharedPtr());
	init();



	m_lkh_alg.setDefaultParaments();
	auto& tsp_pro = dynamic_cast<TravellingSalesman&>(*pro);
	std::string file_name = tsp_pro.fileName();
	{
		m_lkh_alg.ReadProblem(tsp_pro.filePath());
	}
	m_lkh_alg.assignedMemory();
	int mc_sta = 4;
	mc_StagnationTimes.resize(mc_sta);
	std::fill(mc_StagnationTimes.begin(), mc_StagnationTimes.end(), 1);
	for (int idx(2); idx <= mc_sta; ++idx) {
		mc_StagnationTimes[idx - 1] = mc_StagnationTimes[idx - 2] * idx;
	}


	{
		std::string foundSolDir = "//172.24.207.203/share/2018/diaoyiya/paper_com_experiment_data/tsp_hard_problem_info/";
		foundSolDir = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data_paper_com/eax_tsp/runningdata22/";
		std::string filename = file_name+ std::string("_found_sols.txt");
		insertBestSol(foundSolDir + filename, pro);
	}
	m_root.reset(new Indi);
	m_root->initialize(pro->numberVariables(), -1, m_curIter);
	m_pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;

	m_solHash.initialize(rnd, pro->numberVariables() + 10);
	m_solMap.clear();
	m_indis.clear();
	m_maxCurClassFlag = 0;
	m_curVisitedTime = 0;
	m_curEvolveClass = 1;

	generateSolRandom(pro);
	m_expandStage = true;

	

	//for (int i = 0; i < Npop; ++i)
	//{
	//	m_indis.emplace_back(new Indi());
	//	m_indis[i]->initialize(pro->numberVariables(), m_indis.size() - 1, m_curIter);
	//}
	//m_curPop.resize(m_numPop);
	//for (int idx(0); idx < m_curPop.size(); ++idx) {
	//	m_curPop[idx] = m_indis[idx].get();
	//	tKopt->makeRandSol(*m_curPop[idx]); /* Make a random tour */
	//	tKopt->doIt(*m_curPop[idx]);        /* Apply the local search with the 2-opt neighborhood */
	//	updateFitness(m_curPop[idx]);
	//}
//	updateRelationShipOfV(m_curPop, rnd);
//	updateRoot(m_curPop);
	



}
