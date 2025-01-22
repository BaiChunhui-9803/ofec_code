#include "nbn_eax_pop.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../utility/nbn_visualization/nbn_fla.h"

int ofec::PopNBN_EAX::evolve(Problem* pro, Algorithm* alg, Random* rnd)
{
	using namespace NBN_learn2;
	using namespace eax_tsp;
	
	++m_curIter;

	if (m_net.m_promising_areas.empty()) {
		return ofec::kNormalEval;
	}
	
	auto cur = m_net.m_promising_areas.back();
	m_net.m_promising_areas.pop_back();

	double curMaxDis = 0;
	

	//double minDis = cur->getSolInfo()->minDis();
	
	std::vector<std::shared_ptr<NetworkNode>> parents;
	std::vector<IndiEAX_info*> parentInfos;
	std::shared_ptr<RelationNode> sonIter = (*cur->sonList().getHead())->getSharedPtr();
	while ((*sonIter->getNext())->isActive()) {
		auto& curRNode = (*sonIter->getNext());
		auto cur = curRNode->getCurInfo()->getSolInfo();
		auto curType = dynamic_cast<IndiEAX_info*> (cur);
		if (m_curIter - curType->getTime() < m_maxNumStagation) {
			parentInfos.push_back(curType);
			auto curNode = curRNode->getCurInfo();
			parents.push_back(curNode->getSharedPtr());
		}
		else if (m_curIter - curType->getTime() == m_maxNumStagation) {
			auto curNode = curRNode->getCurInfo();
			m_net.m_promising_areas.push_back(curNode->getSharedPtr());
		}
		else break;

		sonIter = (*sonIter->getNext());
	}

	for (auto& it : parents) {
		curMaxDis = std::max(it->getSolInfo()->distance(*cur->getSolInfo()), curMaxDis);
	}

	if (curMaxDis <= m_minDis) {
		m_net.m_found_peaks.push_back(cur);
		double maxFit = 0;
		for (auto& it : m_net.m_found_peaks) {
			maxFit =  std::max(maxFit, it->getSolInfo()->fitness() + m_bestObj);
		}
		std::cout << "curmaxfit\t" << maxFit << std::endl;
	//	for(auto&it)
		return ofec::kNormalEval;
	}




	parents.push_back(cur);
	std::vector<TIndi> newParents;
	// exploration
	if (parents.size() <Npop) {
		int size = Npop - parents.size();
		newParents.resize(size);
		std::vector<std::vector<std::vector<int>>*> totalLinks(parentInfos.size());
		std::vector<double> weights(parentInfos.size());
		std::vector<double> objs(parentInfos.size());
		
		for (int idx(0); idx < objs.size(); ++idx) {
			auto& curIndi = parentInfos[idx]->indi();
			totalLinks[idx] = &curIndi.fLink;
			objs[idx] = curIndi.fEvaluationValue;
		}
		m_gl_calculator.calWeight(weights, objs, pro);
		m_gl_calculator.udpateProMat(totalLinks, weights);
		
		std::vector<int> tour;
		for (int idx(0); idx < newParents.size(); ++idx) {
			m_gl_calculator.generateTSPsolution(tour, pro, rnd);

			auto& it = newParents[idx];
			it.define(pro->numVariables());
			it.toCurSol(tour);

			tKopt->doIt(it);        /* Apply the local search with the 2-opt neighborhood */

		//	fEvaluator->doIt(it);

		}

	}
	
	// exploitation 

	std::vector<TIndi*> totalParenst;
	for (auto& it : parentInfos) {
		totalParenst.push_back(&it->indi());
	}
	for (auto& it : newParents) {
		totalParenst.push_back(&it);
	}
	
	std::vector<TIndi> sons(totalParenst.size());
	
	for (int idx(0); idx < sons.size(); ++idx) {
		sons[idx] = *totalParenst[idx];
	}
	std::vector<int> shuffleIdxs(totalParenst.size(), -1);

	std::vector<int> sortIds(totalParenst.size());
	for (int idx(0); idx < sortIds.size(); ++idx) {
		sortIds[idx] = idx;
	}

	std::sort(sortIds.begin(), sortIds.end(), [&](
		int a, int b) {
		return totalParenst[a]->fEvaluationValue > totalParenst[b]->fEvaluationValue;
	});


	std::vector<std::vector<int>> parentIds(totalParenst.size());
	std::vector<double> minDis(totalParenst.size(), std::numeric_limits<double>::max());
	std::vector<double> objs(totalParenst.size(), std::numeric_limits<double>::max());
	for (int id(0); id< totalParenst.size(); ++id) {
		int idx = sortIds[id];
		auto& vPId = parentIds[idx];
		auto& cursol = *totalParenst[idx];
		auto& curobj = objs[idx];
		auto& curMinDis = minDis[idx];
		for (int id_y(id + 1); id_y < totalParenst.size(); ++id_y) {

			int idy = sortIds[id_y];
			double curDis = cursol.distanceTo(*totalParenst[idy]);
			if (curDis < curMinDis) {
				curMinDis = curDis;
				vPId.clear();
				vPId.push_back(idy);
				curobj = totalParenst[idy]->fEvaluationValue;
			}
			else if (curDis == curMinDis) {
				
				if (curobj > totalParenst[idy]->fEvaluationValue) {
					curobj = totalParenst[idy]->fEvaluationValue;
					vPId.clear();
					vPId.push_back(idy);
				}
				else if(curobj == totalParenst[idy]->fEvaluationValue)
				vPId.push_back(idy);
			}
		}
	}
	//for (int idx(0); idx < shuffleIdxs.size(); ++idx) {
	//	shuffleIdxs[idx] = idx;
	//}
	//rnd->uniform.shuffle(shuffleIdxs.begin(), shuffleIdxs.end());
	

	shuffleIdxs.resize(parentIds.size());
	for (int idx(0); idx < shuffleIdxs.size(); ++idx) {
		auto& curIds = parentIds[idx];
		if (curIds.empty()) {
			shuffleIdxs[idx] = -1;
		}
		else		
		shuffleIdxs[idx] = curIds[rnd->uniform.nextNonStd<int>(0, curIds.size())];
	}

	static int iter = 0;
	std::cout << "curIter" << iter++ << std::endl;

	for (int idx(0); idx < shuffleIdxs.size(); ++idx) {
		if (shuffleIdxs[idx] != -1) {
			tCross->setParents(sons[idx], *totalParenst[shuffleIdxs[idx]], fFlagC, Nch);
			double dis = sons[idx].distanceTo(*totalParenst[shuffleIdxs[idx]]);
		//	std::cout << "cur dis\t" << dis << std::endl;

			tCross->doItWithoutParent(sons[idx], Nch, 1, fFlagC, fEdgeFreq);
		}
	}

	std::vector<std::array<int,2>> curParents;
	std::array<int, 2> curInfo;

	std::vector<TIndi*> activeSons;
	for (int idx(0); idx < shuffleIdxs.size(); ++idx) {
		double dis1 = sons[idx].distanceTo(*totalParenst[idx]);
		if (dis1 != 0) {
			dis1 = sons[idx].distanceTo(*totalParenst[shuffleIdxs[idx]]);
			if (dis1 != 0) {
				activeSons.push_back(&sons[idx]);

				curInfo.front() = idx;
				curInfo.back() = shuffleIdxs[idx];
				curParents.push_back(curInfo);
			}
		}
	}

	for (auto& it : newParents) {
		activeSons.push_back(&it);
	}
	

	std::vector<std::shared_ptr<NBN_learn2::IndividualInfo>> indiInfos(activeSons.size());
	for (int idx(0); idx < activeSons.size(); ++idx) {
		auto& it = indiInfos[idx];
		indiInfos[idx].reset(new IndiEAX_info);
		auto& indi = dynamic_cast<IndiEAX_info&>(*it);
		indi.initialize(pro);
		indi.setIndi(*activeSons[idx]);
		indi.setTime(m_curIter);
	}

	std::vector<std::shared_ptr<NBN_learn2::NetworkNode>> networkNodes(indiInfos.size());
	for (int idx(0); idx < networkNodes.size(); ++idx) {
		networkNodes[idx].reset(new NBN_learn2::NetworkNode);
		networkNodes[idx]->initialize(indiInfos[idx]);
	}


	//for (int idx(0); idx < curParents.size(); ++idx) {
	//	auto& curP = curParents[idx];
	//	for (auto idP : curP) {
	//		if (idP < parents.size()) {
	//			parents[idP]->updateSons(networkNodes[idx], rnd);
	//		}
	//	}
	//}

	std::vector<std::shared_ptr<NBN_learn2::NetworkNode>> bestNodes;
	m_net.updateNeighbor(parents, networkNodes, bestNodes);


	rnd->uniform.shuffle(networkNodes.begin(), networkNodes.end());
	for (auto& it : networkNodes) {
		it->updateSonsons(rnd);
	}


	for (auto& it : bestNodes) {
		m_net.m_promising_areas.push_back(it);
	}


	

	



	return ofec::kNormalEval;
}

void ofec::PopNBN_EAX::initialize(Problem* pro, Random* rnd)
{
	using namespace ofec;
	using namespace eax_tsp;

	m_curIter = 0;
	Npop = 100;
	Nch = 30;
	define(pro, rnd->getSharedPtr());
	init();


	fFlagC[0] = 1;          /* Diversity preservation: 1:Greedy, 2:--- , 3:Distance, 4:Entropy (see Section 4) */
	fFlagC[1] = 1;

	m_gl_calculator.initialize(pro);
	m_net.initialize(rnd->getSharedPtr());

	for (int i = 0; i < Npop; ++i)
	{
		tKopt->makeRandSol(tCurPop[i]); /* Make a random tour */
		tKopt->doIt(tCurPop[i]);        /* Apply the local search with the 2-opt neighborhood */
	}
	
	std::vector<std::shared_ptr<NBN_learn2::IndividualInfo>> indiInfos(tCurPop.size());
	for (int idx(0); idx < tCurPop.size(); ++idx) {
		indiInfos[idx].reset(new IndiEAX_info);
		auto& it = indiInfos[idx];
		auto& indi = dynamic_cast<IndiEAX_info&>(*it);
		
		indi.initialize(pro);
		indi.setIndi(tCurPop[idx]);
		indi.setTime(m_curIter);
	}
	
	std::vector<std::shared_ptr<NBN_learn2::NetworkNode>> networkNodes(indiInfos.size());
	for (int idx(0); idx < networkNodes.size(); ++idx) {
		networkNodes[idx].reset(new NBN_learn2::NetworkNode);
		networkNodes[idx]->initialize(indiInfos[idx]);	
	}
	std::vector<std::shared_ptr<NBN_learn2::NetworkNode>> bestNodes;
	m_net.updateNeighbor(networkNodes, bestNodes);
	for (auto& it : bestNodes) {
		m_net.m_promising_areas.push_back(it);
	}
	

	{
		// input  bestSols

		std::string foundSolDir = "//172.24.207.203/share/2018/diaoyiya/paper_com_experiment_data/tsp_hard_problem_info/";
		std::string filename = "eax_0091_found_sols.txt";

		std::vector<int> cursolx;
		int solId(0);
		long long randnum(0);
		double obj(0);
		std::string solutiontag;

		std::ifstream in(foundSolDir + filename);

		while (in >> solId >> randnum >> obj >> solutiontag) {
			CAST_TSP(pro)->inputSol(in, cursolx);
			TIndi indi;
			indi.define(pro->numVariables());
			indi.toCurSol(cursolx);
			m_bestSols.push_back(indi);
		}
		in.close();

		m_bestObj = std::numeric_limits<double>::max();
		for (auto& it : m_bestSols) {
			m_bestObj = std::min(it.fEvaluationValue, m_bestObj);
		}

	}
	
}

void ofec::PopNBN_EAX::calNBN(std::vector<int>& belong,
	std::vector<double>& dis2parent, 
	std::vector<double>& fitness, std::vector<int>& optFlag, 
	std::vector<double>& curoptFit,
	ofec::Random* rnd)
{
	using namespace NBN_learn2;

	std::shared_ptr<NetworkNode> root;
	if (m_net.m_found_peaks.empty()) {
		root = m_net.m_promising_areas.back();
	}
	else {
		root = m_net.m_found_peaks.front();
	}

	std::vector<std::shared_ptr<IndividualInfo>> sols;
	
	m_net.getGraph(root, belong, dis2parent, fitness, sols);
	if (!ofec::NBN_FLA_info::judgeNBN(belong, fitness)) {
		std::cout << "nbn error" << std::endl;
	}
	

	std::vector<double> optDis(m_bestSols.size(), std::numeric_limits<double>::max());
	std::vector<double> optFit(m_bestSols.size(), -std::numeric_limits<double>::max());
	optFlag.resize(m_bestSols.size(), -1);

	curoptFit.resize(m_bestSols.size());
	//auto m_pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	for (int idx(0); idx < optFit.size(); ++idx) {
		curoptFit[idx] = m_bestSols[idx].fEvaluationValue*-1;
	}
	for (int idx(0); idx < m_bestSols.size(); ++idx) {
		auto& curMinDis = optDis[idx];
		auto& curFit = optFit[idx];
		auto& cursol = m_bestSols[idx];
		auto& id = optFlag[idx];
		for (int idy(0); idy < sols.size(); ++idy) {
			auto other = dynamic_cast<IndiEAX_info*> (sols[idy].get());
			double curdis = cursol.distanceTo(other->indi());
			if (curMinDis > curdis) {
				curMinDis = curdis;
				curFit = other->fitness();
				id = idy;
			}
			else if (curMinDis == curdis) {
				if (other->fitness() > curFit) {
					curFit = other->fitness();
					id = idy;
				}
				else if (other->fitness() == curFit && rnd->uniform.next() < 0.5) {
					id = idy;
				}
			}
		}
	}


	for (auto& it : optDis) {
		std::cout << it << "\t";
	}
	std::cout<<std::endl;
}
