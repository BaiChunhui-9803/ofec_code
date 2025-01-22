#include "nbn_eax_pop_v4.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../utility/nbn_visualization/nbn_fla.h"


ofec::PopNBN_EAX_V4::Indi* ofec::PopNBN_EAX_V4::RelationNode::getCur() {
	return m_cur;
}
ofec::PopNBN_EAX_V4::Indi* ofec::PopNBN_EAX_V4::RelationNode::getParent() {
	if (m_list == nullptr) {
		return nullptr;
	}
	else {
		return m_list->m_cur;
	}
}


void ofec::PopNBN_EAX_V4::RelationNode::removeFromLink() {
	if (m_list != nullptr) {

		m_before->m_after = m_after;
		m_after->m_before = m_before;
		m_after = m_before = nullptr;
		++m_list->m_num;
		m_list = nullptr;

		m_after = nullptr;
		m_before = nullptr;
	}
}




void ofec::PopNBN_EAX_V4::RelationNode::insertToLink(RelationList* list) {
	
	list->insertNode(this);
	//m_list = list;
}


void ofec::PopNBN_EAX_V4::RelationList::insertNode(RelationNode* cur) {
	
	cur->removeFromLink();
	cur->m_after = m_head->m_after;
	cur->m_before = m_head.get();
	cur->m_list = this;
	m_head->m_after->m_before = cur;
	m_head->m_after = cur;
	++m_num;
	
	//m_sumDis2parent += cur.m_cur->m_dis2parent;
}

void ofec::PopNBN_EAX_V4::RelationList::removeNode(RelationNode* cur)
{
	cur->removeFromLink();
}


void ofec::PopNBN_EAX_V4::RelationList::clearLink() {
	auto head = m_head.get();
	while (head->m_after->m_cur != nullptr) {
		head->m_after->removeFromLink();
	}
	
}


void ofec::PopNBN_EAX_V4::getDiverSols(std::vector<Indi*>& sols) {
	
}

bool ofec::PopNBN_EAX_V4::judgeInside(
	std::queue<Indi*>& que, 
	Indi* curindi, bool flagInside, unsigned long long curClass) {
	while (!que.empty()) {
		auto& curNei = que.front();
		que.pop();
		auto curDis = curNei->distanceTo(*curindi);
		if (curDis <= curNei->m_maxRadius) {
			flagInside = true;
			break;
		}

		RelationNode* head = curNei->m_sons.m_head.get();
		while (head->m_after->m_cur != nullptr) {
			head = head->m_after;
			if (head->m_cur->m_curClassFlag == curClass
				&& head->m_cur->m_curVisitedId != m_curVisitedTime) {
				que.push(head->m_cur);
			}
		}

		auto parent = curNei->m_parentNode->getParent();
		if (parent != nullptr && parent->m_curClassFlag == curClass
			&& head->m_cur->m_curVisitedId != m_curVisitedTime)
			que.push(parent);


	}
	return flagInside;
}


void ofec::PopNBN_EAX_V4::insertBestSol(const std::string& filepath, ofec::Problem* pro) {
	using namespace ofec;
	using namespace eax_tsp;
	{
		// input  bestSols
		std::vector<int> cursolx;
		int solId(0);
		long long randnum(0);
		double obj(0);
		std::string solutiontag;

		std::ifstream in(filepath);
		if (!in.is_open()) {
			std::cout << "no opt solutions" << std::endl;
			return;
		}

		while (in >> solId >> randnum >> obj >> solutiontag) {
			CAST_TSP(pro)->inputSol(in, cursolx);
			TIndi indi;
			indi.define(pro->numberVariables());
			indi.toCurSol(cursolx);
			fEvaluator->doIt(indi);
			updateFitness(&indi);
			indi.updateSolBase();

			m_bestSols.push_back(indi);
		}
		in.close();

		m_bestObj = std::numeric_limits<double>::max();
		for (auto& it : m_bestSols) {
			m_bestObj = std::min(it.fEvaluationValue, m_bestObj);
		}

	}

}

void ofec::PopNBN_EAX_V4::filterNondominatedState(
	std::vector<Indi*>& totalStates) const{
	std::sort(totalStates.begin(), totalStates.end(),
		[&](Indi* a, Indi* b) {
		if (a->fitness() == b->fitness()) {
			return a->m_dis2parent > b->m_dis2parent;
		}
		else return a->fitness() > b->fitness();
	});

	double beforeFit = totalStates.front()->fitness();
	double beforeDis = 0;
	int curIdx = 0;
	bool saveFlag = false;
	for (int idx(0); idx < totalStates.size(); ++idx) {
		if (totalStates[idx]->m_exploitationFlag)continue;
		saveFlag = false;
		if (totalStates[idx]->fitness() == beforeFit) {
			if (totalStates[idx]->m_dis2parent == beforeDis) {
				saveFlag = true;
			}
		}
		else if (totalStates[idx]->m_dis2parent > beforeDis) {
			beforeDis = totalStates[idx]->m_dis2parent;
			beforeFit = totalStates[idx]->fitness();
			saveFlag = true;
		}
		if (saveFlag) {
			++curIdx;
			if (curIdx != idx) {
				totalStates[curIdx] = totalStates[idx];
			}
		}
	}


	totalStates.resize(curIdx + 1);
}

void ofec::PopNBN_EAX_V4::updateNondiminatedStates(ofec::Problem* pro, ofec::Random* rnd) {

	m_nondiminatedStates.clear();
	std::vector<int> cursolx;



	


	for (int idx(0); idx < m_indis.size(); ++idx) {
		auto cur = m_indis[idx].get();
		if (!cur->m_exploitationFlag) {
			m_nondiminatedStates.push_back(cur);
		}
	}


	if (m_nondiminatedStates.size() <=1) {
		return ;
	}
	
	std::sort(m_nondiminatedStates.begin(), m_nondiminatedStates.end(), 
		[&](Indi* a, Indi*b) {
		if (a->fitness() == b->fitness()) {
			return a->m_dis2parent > b->m_dis2parent;
		}
		else return a->fitness() > b->fitness();
	});

	double beforeFit = m_nondiminatedStates.front()->fitness();
	double beforeDis = 0;
	int curIdx = 0;
	bool saveFlag = false;
	for (int idx(0); idx < m_nondiminatedStates.size(); ++idx) {
		if (m_nondiminatedStates[idx]->m_exploitationFlag)continue;
		saveFlag = false;
		if (m_nondiminatedStates[idx]->fitness() == beforeFit) {
			if (m_nondiminatedStates[idx]->m_dis2parent == beforeDis) {
				saveFlag = true;
			}
		}
		else if (m_nondiminatedStates[idx]->m_dis2parent > beforeDis) {
			beforeDis = m_nondiminatedStates[idx]->m_dis2parent;
			beforeFit = m_nondiminatedStates[idx]->fitness();
			saveFlag = true;
		}
		if (saveFlag) {
			++curIdx;
			if (curIdx!= idx) {
				m_nondiminatedStates[curIdx] = m_nondiminatedStates[idx];
			}
		}
	}


	m_nondiminatedStates.resize(curIdx + 1);

	//std::vector<Indi*> totalSols(m_indis.size());
	//for (int idx(0); idx < totalSols.size(); ++idx) {
	//	totalSols[idx] = m_indis[idx].get();
	//}

	//std::vector<Indi*> newSols;
	//std::shared_ptr<Indi> newsol;

	//for (int idx(0); idx < m_nondiminatedStates.size(); ++idx) {
	//	auto cur = m_nondiminatedStates[idx];
	//	/*if (!cur->m_exploitationFlag)*/ {

	//		cur->transferSol(cursolx);
	//		for (auto& it : cursolx) ++it;
	//		auto cost = m_lkh_alg.LKlocalSearch(cursolx, m_moveType);
	//		for (auto& it : cursolx) --it;
	//		auto hash = m_solHash.calHash(cursolx);
	//		if (m_solMap.find(hash) == m_solMap.end()) {

	//			newsol.reset(new Indi);
	//			newsol->initialize(pro->numberVariables(), -1, m_curIter);
	//			newsol->toCurSol(cursolx);
	//			fEvaluator->doIt(*newsol);
	//			updateFitness(newsol.get());
	//			auto& indi = m_solMap[hash];
	//			if (indi == nullptr) {
	//				m_indis.emplace_back(std::move(newsol));
	//				m_indis.back()->setSolId(m_indis.size() - 1);
	//				m_indis.back()->setRndId(hash);
	//				m_indis.back()->setCurClass(m_curClassFlag);
	//				indi = m_indis.back().get();
	//				newSols.push_back(indi);
	//			}


	//		}
	//	//	m_nondiminatedStates.push_back(cur);

	//	}

	//}
	//updateRelationShipOf2V(newSols, totalSols, rnd);
	//m_nondiminatedStates.clear();
	//for (auto& it : newSols) {
	//	m_nondiminatedStates.push_back(it);
	//}

}

int ofec::PopNBN_EAX_V4::exploitationState_v1(ofec::Problem* pro, ofec::Random* rnd) {
	Indi* cur = nullptr;

	//	updateNondiminatedStates();



	while (!m_stagnationIndis.empty()) {
		cur = m_stagnationIndis.back();
		m_stagnationIndis.pop_back();
		if (cur->m_curClassFlag != 0)break;
	}
	if (cur == nullptr) {
		return EvaluationTag::kNormalEval;
	}


	m_curIter = 0;
	//	m_id2curId.resize(m_indis.size());
	increaseClassFlag();
	//m_root->m_curClassFlag = m_curClassFlag;


	std::queue<Indi*> que;
	cur->m_curClassFlag = m_maxCurClassFlag;
	que.push(cur);
	//m_stagnationIndis.pop_back();

	std::cout << "curPeak\t" << que.front()->m_solId << std::endl;
	while (!que.empty()) {
		auto cur = que.front();
		que.pop();
		if (m_curPop.size() < m_numPop)
			m_curPop.push_back(cur);
		//	m_id2curId[cur->m_solId] = m_curPop.size() - 1;


		RelationNode* head = cur->m_sons.m_head.get();
		auto& curMaxDis = cur->m_maxRadius;
		curMaxDis = 0;
		while (head->m_after->m_cur != nullptr) {
			head = head->m_after;

			que.push(head->m_cur);
			curMaxDis = std::max(head->m_cur->m_maxRadius, curMaxDis);
		}

	}
	if (m_curPop.size() < m_numPop) {
		increaseVisitedTimes();

		std::vector<double> radiusRate(m_curPop.size(), 1.0);
		std::vector<std::shared_ptr<Indi>> newRandSols;
		std::vector<std::shared_ptr<Indi>> crossSols;
		int dim = pro->numberVariables();
		std::vector<int> cursol;
		std::shared_ptr<Indi> curindi;
		int maxLoop = 1000;

		std::vector<Indi*> borders;
		getBorder(m_curPop, borders);

		//while (newRandSols.size()+ m_curPop.size() < m_numPop&& maxLoop--) {
		//	curindi.reset(new Indi);
		//	curindi->initialize(dim, -1, m_curIter);
		//	auto& randIndi = m_curPop[rnd->uniform.nextNonStd<int>(0, m_curPop.size())];

		//	randIndi->transferSol(cursol);
		//	int solIdx = m_id2curId[randIndi->m_solId];
		//	int radius = rnd->uniform.nextNonStd<int>(1, randIndi->m_maxRadius/4.0);
		//	if (radius) {

		//		while (radius--) {
		//			auto a = rnd->uniform.nextNonStd<int>(0, dim);
		//			auto b = rnd->uniform.nextNonStd<int>(0, dim);
		//			std::swap(cursol[a], cursol[b]);
		//		}
		//		curindi->toCurSol(cursol);
		//		tKopt->doIt(*curindi);        /* Apply the local search with the 2-opt neighborhood */
		//		updateFitness(curindi.get());
		//		
		//		auto hash = m_solHash.calHash(curindi->fLink);
		//		curindi->setRndId(hash);
		//		if (m_solMap.find(hash) != m_solMap.end()) {
		//			continue;
		//		}
		//		bool flagInside = true;
		//		//std::queue<Indi*> que;
		//		//que.push(randIndi);
		//		//increaseVisitedTimes();
		//		//randIndi->m_curVisitedId = m_curVisitedTime;
		//		////visited[m_id2curId[randIndi->m_solId]] = curVisited;


		//		//judgeInside(que, curindi.get(), flagInside);
		//		double minDis = 0;
		//		
		//		if (flagInside) {
		//			newRandSols.push_back(curindi);
		//			//radiusRate[solIdx] = radiusRate[solIdx] * 0.9 + 1.0 * 0.1;
		//		}

		//	}
		//}

		maxLoop = 1000;

		clearEdgeFreq();
		std::vector<int> cursolx;
		while (crossSols.size() + m_curPop.size() < m_numPop
			&& maxLoop--) {
			auto a = m_curPop[rnd->uniform.nextNonStd<int>(0, m_curPop.size())];
			auto b = borders[rnd->uniform.nextNonStd<int>(0, borders.size())];
			curindi.reset(new Indi);
			curindi->initialize(dim, -1, m_curIter);
			curindi->copySol(*a);

			tCross->setParents(*curindi, *b, fFlagC, Nch);
			tCross->doIt(*curindi, *b, Nch, 1, fFlagC, fEdgeFreq);
			updateFitness(curindi.get());
			curindi->transferSol(cursolx);
			auto hash = m_solHash.calHash(cursolx);
			curindi->setRndId(hash);
			if (m_solMap.find(hash) != m_solMap.end()) {
				continue;
			}


			double minDis = std::numeric_limits<double>::max();
			for (auto& it : m_curPop) {
				minDis = std::min(minDis, it->distanceTo(*curindi));
			}

			bool flagInside = true;
			for (auto& it : borders) {
				double curdis = it->distanceTo(*curindi);
				if (curdis < minDis) {
					flagInside = false;
					break;
				}
			}

			if (flagInside) {
				crossSols.push_back(curindi);
				//radiusRate[solIdx] = radiusRate[solIdx] * 0.9 + 1.0 * 0.1;
			}
		}

		//std::vector<int> cursolx;
		for (int idx(0); idx < crossSols.size(); ++idx) {
			crossSols[idx]->transferSol(cursolx);
			auto hash = m_solHash.calHash(cursolx);
			auto& indi = m_solMap[hash];
			if (indi == nullptr) {
				m_indis.emplace_back(std::move(crossSols[idx]));
				m_indis.back()->setSolId(m_indis.size() - 1);
				m_indis.back()->setRndId(hash);
				m_indis.back()->setCurClass(m_maxCurClassFlag);
				indi = m_indis.back().get();
				m_curPop.push_back(indi);
			}

		}

		//NBN_indi* best = m_curPop.front();
		//double maxFit = m_curPop.front()->fitness();
		//for (auto& it : m_curPop) {
		//	if (maxFit < it->fitness()) {
		//		maxFit = it->fitness();
		//	}
		//}
		for (auto& it : m_curPop) {
			//if (it->fitness() == maxFit) {
			//	it->setRoot(m_root.get());
			//}
			it->updateState();
		}

		updateRelationShipOfV(m_curPop, rnd);


		std::vector< std::vector<std::vector<int>>*> links;
		for (auto& it : m_curPop) {
			links.push_back(&it->fLink);
		}
		calEdgeFreq(links);

	}

}

void ofec::PopNBN_EAX_V4::getSubTree(Indi* cur, std::vector<Indi*>& subtrees) {

	std::queue<Indi*> que;
	que.push(cur);
	//m_stagnationIndis.pop_back();

	std::cout << "curPeak\t" << que.front()->m_solId << std::endl;
	std::cout << "curdis2opt\t" << cur->distanceTo(m_bestSols.front()) << std::endl;
	std::cout << "curdis2parent\t" << cur->m_dis2parent << std::endl;

	
	if (cur->distanceTo(m_bestSols.front()) <= cur->m_dis2parent) {
		int stop = -1;
	}
	
	while (!que.empty()) {
		auto cur = que.front();
		que.pop();
	//	if (subtrees.size() < m_numPop)
		subtrees.push_back(cur);
		//	m_id2curId[cur->m_solId] = m_curPop.size() - 1;
		RelationNode* head = cur->m_sons.m_head.get();
		
		//auto& curMaxDis = cur->m_maxRadius;
		//curMaxDis = 0;
		while (head->m_after->m_cur != nullptr) {
			head = head->m_after;
			if(head->m_cur->m_dis2parent<=cur->m_dis2parent)
			que.push(head->m_cur);
		//	curMaxDis = std::max(head->m_cur->m_maxRadius, curMaxDis);
		}

	}
}

void ofec::PopNBN_EAX_V4::getTotalSubTree(Indi* cur, std::vector<Indi*>& subtrees)const{
	std::queue<Indi*> que;
	que.push(cur);
	//m_stagnationIndis.pop_back();
	while (!que.empty()) {
		auto cur = que.front();
		que.pop();
		//	if (subtrees.size() < m_numPop)
		subtrees.push_back(cur);
		//	m_id2curId[cur->m_solId] = m_curPop.size() - 1;
		RelationNode* head = cur->m_sons.m_head.get();
		//auto& curMaxDis = cur->m_maxRadius;
		//curMaxDis = 0;
		while (head->m_after->m_cur != nullptr) {
			head = head->m_after;
			///if (head->m_cur->m_dis2parent <= cur->m_dis2parent)
			que.push(head->m_cur);
			//	curMaxDis = std::max(head->m_cur->m_maxRadius, curMaxDis);
		}

	}
}


void ofec::PopNBN_EAX_V4::getDiversitySolInBasin(Indi* cur, std::vector<Indi*>& sols)const {
	std::queue<Indi*> que;
	que.push(cur);
	//m_stagnationIndis.pop_back();

	std::vector<double> dis;
	double meanDis(0), stdDis(0);
	while (!que.empty()) {
		auto cur = que.front();
		que.pop();
		RelationNode* head = cur->m_sons.m_head.get();
		dis.clear();
		while (head->m_after->m_cur != nullptr) {
			head = head->m_after;
			///if (head->m_cur->m_dis2parent <= cur->m_dis2parent)
			que.push(head->m_cur);
			dis.push_back(head->m_cur->m_dis2parent);
			//	curMaxDis = std::max(head->m_cur->m_maxRadius, curMaxDis);
		}

		calMeanAndStd(dis, meanDis, stdDis);
	
		double filterDis = meanDis + stdDis * 3;
		head = cur->m_sons.m_head.get();

		while (head->m_after->m_cur != nullptr) {
			head = head->m_after;
			if (head->m_cur->m_dis2parent > filterDis) {
				sols.push_back(head->m_cur);
			}
		}
	}
}


bool ofec::PopNBN_EAX_V4::generateSols(Indi* cur, 
	std::shared_ptr<Indi>& newSol,
	ofec::Problem* pro, ofec::Random* rnd) {
	newSol.reset(new Indi);
	newSol->initialize(pro->numberVariables(), -1, m_curIter);
	int radius = rnd->uniform.nextNonStd<int>(1,  std::max<int>(1,cur->m_dis2parent / 4.0));
	radius = cur->m_dis2parent/4;
	radius = 19 / 4;
	int dim = pro->numberVariables();
	std::vector<int> cursol;
	
	if (radius) {
		cur->transferSol(cursol);
	

		while (radius--) {
			auto a = rnd->uniform.nextNonStd<int>(0, dim);
			auto b = rnd->uniform.nextNonStd<int>(0, dim);
			std::swap(cursol[a], cursol[b]);
		}
		newSol->toCurSol(cursol);
		tKopt->doIt(*newSol);        /* Apply the local search with the 2-opt neighborhood */
	//	fEvaluator->doIt(*newSol);
		updateFitness(newSol.get());
		

		double dis = newSol->distanceTo(*cur);
		/*if (dis < cur->m_dis2parent)*/ {
			newSol->transferSol(cursol);
			auto hash = m_solHash.calHash(cursol);
			newSol->setRndId(hash);
			return true;
			//if (m_solMap.find(hash) == m_solMap.end()) {
			//	return true;
			//}
		}
		//else {
		//	newSol->toCurSol(cursol);
		//	//tKopt->doIt(*newSol);        /* Apply the local search with the 2-opt neighborhood */
		//	updateFitness(newSol.get());
		//	auto hash = m_solHash.calHash(cursol);
		//	newSol->setRndId(hash);
		//	if (m_solMap.find(hash) == m_solMap.end()) {
		//		return true;
		//	}
		//}
	}

	return false;
}


void ofec::PopNBN_EAX_V4::setBasinsCurClass(Indi* cur, int curClass) {
	std::vector<Indi*> basins;
	getTotalSubTree(cur, basins);
	for (auto& it : basins) {
		it->m_curClassFlag = curClass;
	}
}




ofec::PopNBN_EAX_V4::Indi* ofec::PopNBN_EAX_V4::addSolToIndis(
	std::shared_ptr<Indi>& cur,
	unsigned long long hashValue) {
	auto& indi = m_solMap[hashValue];
	if (indi == nullptr) {
		m_indis.emplace_back(std::move(cur));
		m_indis.back()->setSolId(m_indis.size() - 1);
		m_indis.back()->setRndId(hashValue);
		//m_indis.back()->setCurClass(m_curEvolveClass);
		indi = m_indis.back().get();

	//	newsols.push_back(indi);
	}

	return indi;
	
}

void ofec::PopNBN_EAX_V4::generateRandomSols(Indi* cur, 
	std::vector<Indi*>& newsols, int numSamples,int maxSamples, 
	ofec::Problem*pro, ofec::Random*rnd) {
	std::shared_ptr<Indi> newSol;

	while (maxSamples--) {
		if (newsols.size() >= numSamples)break;
		newSol.reset(new Indi);
		newSol->initialize(pro->numberVariables(),-1,m_curIter);
		newSol->copySol(*cur);
		

		int radius = rnd->uniform.nextNonStd<int>(1, std::max<int>(1, cur->m_curSearchDis / 4.0));
		int dim = pro->numberVariables();
		std::vector<int> cursol;

		if (radius) {
			cur->transferSol(cursol);


			while (radius--) {
				auto a = rnd->uniform.nextNonStd<int>(0, dim);
				auto b = rnd->uniform.nextNonStd<int>(0, dim);
				std::swap(cursol[a], cursol[b]);
			}
			newSol->toCurSol(cursol);
			tKopt->doIt(*newSol);        /* Apply the local search with the 2-opt neighborhood */
			updateFitness(newSol.get());

			newSol->transferSol(cursol);
			auto hash = m_solHash.calHash(cursol);
			
			auto& indi = m_solMap[hash];
			if (indi == nullptr) {
				m_indis.emplace_back(std::move(newSol));
				m_indis.back()->setSolId(m_indis.size() - 1);
				m_indis.back()->setRndId(hash);
				m_indis.back()->setCurClass(m_curEvolveClass);
				indi = m_indis.back().get();

				newsols.push_back(indi);
			}

		}
	}
}


void ofec::PopNBN_EAX_V4::anayalizeNewSols(std::vector<Indi*>& newsols, ofec::Random* rnd) {
	updateRelationShipOfV(newsols, rnd);
	std::sort(newsols.begin(), newsols.end(), [&](Indi*a , Indi* b) {
		return a->m_dis2parent > b->m_dis2parent;
	});

	std::cout << "curIndo" << std::endl;
}


void ofec::PopNBN_EAX_V4::expandSearch_v1(Problem* pro, Algorithm* alg, Random* rnd) {

	increaseClassFlag();
	m_curEvolveClass = m_maxCurClassFlag;



	updateNondiminatedStates(pro, rnd);
	m_expanedStage = true;
	//int stop = -1;
	Indi* cur = nullptr;
	for (auto& it : m_nondiminatedStates) {
		if (it->m_dis2parent <= 3)continue;
		if (!it->m_exploitationFlag) {

			cur = it;
			break;
		}
	}

	//cur = m_indis[6506].get();
	//m_peakAround = cur;
	//cur->m_exploitationFlag = true;
	//cur->m_curSearchDis = cur->m_dis2parent;
	//cur->m_dis2parent *= 0.9;
	//std::vector<Indi*> newsols;
	//generateRandomSols(cur,newsols,
	//	1e3, 1e4,pro, rnd);
	//
	//anayalizeNewSols(newsols, rnd);

	//std::vector<Indi*> cursons;
	//{
	//	RelationNode* head = cur->m_sons.m_head.get();
	//	while (head->m_after->m_cur != nullptr) {
	//		head = head->m_after;
	//		cursons.push_back(head->m_cur);
	//	}
	//}
	//anayalizeNewSols(cursons, rnd);





	std::cout << "curPeak\t" << cur->m_solId << std::endl;
	std::cout << "curdis2opt\t" << cur->distanceTo(m_bestSols.front()) << std::endl;
	std::cout << "curdis2parent\t" << cur->m_dis2parent << std::endl;



	//	getSubTree(cur, m_curPop);
		// genereate randomSols;
	std::vector<Indi*> basins;
	getTotalSubTree(cur, basins);

	std::vector<Indi*> randSol;
	std::shared_ptr<Indi> cursol;
	int maxSampleTime = m_maxSampleTime;

	std::vector<Indi*> initsols;
	for (auto& it : basins) {
		//it->m_curClassFlag = m_curEvolveClass;
		//if (it->m_dis2parent >= 4)
		{
			initsols.push_back(it);
		}
	}
	for (auto& it : basins) {
		it->m_curClassFlag = m_curEvolveClass;
		it->m_stagnation_time = 0;
	}

	//for (auto& it : basins) {
	//	double curdis = cur->distanceTo(*it);
	//	if (curdis <= cur->m_dis2parent) {
	//		m_curPop.push_back(it);
	//	}
	//}
	//
	//rnd->uniform.shuffle(m_curPop.begin(), m_curPop.end());
	//if (m_curPop.size() > m_numPop) {
	//	m_curPop.resize(m_numPop);
	//}
	//m_curPop.push_back(cur);
	int leftSize = std::max<int>(m_numPop - basins.size(), m_numPop / 2.0);
	leftSize = m_numPop;
	while (randSol.size() < m_numPop && --maxSampleTime) {
		auto initSol = initsols[rnd->uniform.nextNonStd<int>(0, initsols.size())];
		if (generateSols(initsols.front(), cursol, pro, rnd)) {

			auto hash = cursol->m_randId;
			auto& indi = m_solMap[hash];
			if (indi == nullptr) {
				m_indis.emplace_back(std::move(cursol));
				m_indis.back()->setSolId(m_indis.size() - 1);
				m_indis.back()->setRndId(hash);
				m_indis.back()->setCurClass(m_curEvolveClass);
				indi = m_indis.back().get();
			}
			randSol.push_back(indi);
		}
	}
	//maxSampleTime = m_maxSampleTime;

	//while (randSol.size() < leftSize && --maxSampleTime) {
	//	auto initSol = initsols[rnd->uniform.nextNonStd<int>(0, initsols.size())];
	//	if (generateSols(initSol, cursol, pro, rnd)) {
	//		auto hash = cursol->m_randId;
	//		auto& indi = m_solMap[hash];
	//		if (indi == nullptr) {
	//			m_indis.emplace_back(std::move(cursol));
	//			m_indis.back()->setSolId(m_indis.size() - 1);
	//			m_indis.back()->setRndId(hash);
	//			m_indis.back()->setCurClass(m_curClassFlag);
	//			indi = m_indis.back().get();
	//		}
	//		randSol.push_back(indi);
	//	}
	//}

	std::vector<Indi*> newSols;
	for (auto& it : randSol) {
		newSols.push_back(it);
	}

	//for (int idx(0); idx < randSol.size(); ++idx) {

	//	auto hash = randSol[idx]->m_randId;
	//	auto& indi = m_solMap[hash];
	//	if (indi == nullptr) {
	//		m_indis.emplace_back(std::move(randSol[idx]));
	//		m_indis.back()->setSolId(m_indis.size() - 1);
	//		m_indis.back()->setRndId(hash);
	//		m_indis.back()->setCurClass(m_curClassFlag);
	//		indi = m_indis.back().get();
	//	}
	//	newSols.push_back(indi);
	//}


	//int curId(0);
	//while (newSols.size() < m_numPop&&curId< basins.size()) {
	//	newSols.push_back(basins[curId++]);
	//}

	swap(newSols, m_curPop);

	{

		std::vector<Indi*> neighbors;

		getNeighbors(basins, neighbors);
		for (auto& it : basins) {
			neighbors.push_back(it);
		}
		updateRelationShipOf2V(m_curPop, neighbors, rnd);
		updateRelationShipOfV(m_curPop, rnd);
		updateRoot(m_curPop);



		std::vector< std::vector<std::vector<int>>*> links;
		for (auto& it : m_curPop) {
			links.push_back(&it->fLink);
		}
		calEdgeFreq(links);

		fFlagC[0] = 4;          /* Diversity preservation: 1:Greedy, 2:--- , 3:Distance, 4:Entropy (see Section 4) */
		fFlagC[1] = 1;          /* Eset Type: 1:Single-AB, 2:Block2 (see Section 3) */

	}


	for (auto& it : m_curPop) {
		it->m_stagnation_time = 0;
	}
}


void ofec::PopNBN_EAX_V4::expandSearch_v2(Problem* pro, Algorithm* alg, Random* rnd) {
	m_maxStag = 100;
	increaseClassFlag();
	m_curEvolveClass = m_maxCurClassFlag;
	m_peakAround = m_peaks.back();
//	m_peaks.pop_back();
	std::vector<Indi*>basins;
	getTotalSubTree(m_peakAround, basins);

	std::cout << "curPeak\t" << m_peakAround->m_solId << std::endl;
	std::cout << "curdis2opt\t" << m_peakAround->distanceTo(m_bestSols.front()) << std::endl;
	std::cout << "curdis2parent\t" << m_peakAround->m_dis2parent << std::endl;
	std::vector<Indi*> indis;
	for (auto& it : basins) {
		if (it != m_peakAround) {
			indis.push_back(it);
		}
	}

	for (int idx(0); idx < m_numPop / 2.0; ++idx) {
		m_curPop.push_back(indis[idx]);
		m_curPop.back()->setCurClass(m_curEvolveClass);
	}


	std::sort(indis.begin(), indis.end(), [](Indi* a, Indi*b) {
		return a->m_dis2parent > b->m_dis2parent;
	});
	int curIdx(0);
	while (m_curPop.size() < m_numPop && curIdx < indis.size()) {
		if (indis[curIdx]->m_curClassFlag != m_curEvolveClass) {
			m_curPop.push_back(indis[curIdx]);
			m_curPop.back()->setCurClass(m_curEvolveClass);

		}
		++curIdx;
	}
	//
	//if (indis.size() > m_numPop) {
	//	indis.resize(m_numPop);
	//	swap(m_curPop, indis);
	//}

	for (auto& it : m_curPop) {
		it->m_stagnation_time = 0;
	}


	//filterNondominatedState(indis);
	//
	//std::vector<Indi*> filterSols;
	//getDiversitySolInBasin(m_peakAround, filterSols);

	//int stop = -1;
	
}


void ofec::PopNBN_EAX_V4::expandSearch_test(Problem* pro, Algorithm* alg, Random* rnd) {

	clearEdgeFreq();
	std::vector<Indi*> totalsols;
	for (auto& it : m_indis) {
		totalsols.push_back(it.get());
	}
	std::vector<Indi*> cursols;
	std::vector<int> parentIds;

	std::vector<Indi*> newsols;
	{
		{
			std::cout << "curIter\t" << m_curIter << std::endl;
			std::vector<int> belong;
			std::vector<double> dis2parent;
			std::vector<double> fitness;
			std::vector<int> optFlag;
			std::vector<double> curoptFit;
			calNBN(belong, dis2parent, fitness, optFlag, curoptFit,
				rnd);


			analysis(optFlag);

			
			parentIds.clear();
			int curId = optFlag.front();
			parentIds.push_back(curId);
			while (m_indis[curId]->m_parentNode->getParent() != nullptr) {
				curId = m_indis[curId]->m_parentNode->getParent()->m_solId;
				parentIds.push_back(curId);
			}
		}
	}


	{
		for (auto& it:parentIds) {
			cursols.push_back(m_indis[it].get());
		}

		std::shared_ptr<Indi> cur;
		
		Indi* a = nullptr;
		Indi* b = nullptr;
		std::vector<int> solx;
		while (true) {

			a = cursols[rnd->uniform.nextNonStd<int>(0, cursols.size())];
			b = cursols[rnd->uniform.nextNonStd<int>(0, cursols.size())];


			cur.reset(new Indi);
			cur->initialize(pro->numberVariables(), -1, m_curIter);
			cur->copySol(*a);
			tCross->setParents(*cur, *b, fFlagC, Nch);
			tCross->doIt(*cur, *b, Nch, 1, fFlagC, fEdgeFreq);
			updateFitness(cur.get());

			cur->transferSol(solx);
			auto hash = m_solHash.calHash(solx);
			auto& indi = m_solMap[hash];
			if (indi == nullptr) {

				std::cout << "create newsols" << std::endl;

				m_indis.emplace_back(std::move(cur));
				m_indis.back()->setSolId(m_indis.size() - 1);
				m_indis.back()->setRndId(hash);
				indi = m_indis.back().get();
				indi->m_curClassFlag = m_curEvolveClass;

				newsols.clear();
				
				newsols.push_back(indi);
				updateRelationShipOf2V(
					newsols, totalsols,rnd);

				totalsols.push_back(indi);


				{
					{
						std::cout << "curIter\t" << m_curIter << std::endl;
						std::vector<int> belong;
						std::vector<double> dis2parent;
						std::vector<double> fitness;
						std::vector<int> optFlag;
						std::vector<double> curoptFit;
						calNBN(belong, dis2parent, fitness, optFlag, curoptFit,
							rnd);


						analysis(optFlag);
					}
				}
				
			}


		}
		
	}
	
}

void ofec::PopNBN_EAX_V4::expandSearch_test2(Problem* pro, Algorithm* alg, Random* rnd) {
	
	increaseClassFlag();
	m_curEvolveClass = m_maxCurClassFlag;
	
}


int ofec::PopNBN_EAX_V4::evolve2(Problem* pro, Algorithm* alg, Random* rnd) {


	if (m_curPop.empty()) {
	
		expandSearch_v1(pro, alg, rnd);
	}


	//// for test

	//{
	//	int testSolId = 5906;
	//	double dis = 91;
	//	if (m_indis.size() > testSolId) {
	//		double curdis = m_indis[testSolId]->distanceTo(m_bestSols.front());
	//		std::cout << "testSolId\t" << testSolId << "\tcurDis\t" << curdis << "\tdisToParent\t" << m_indis[testSolId]->m_dis2parent<< std::endl;
	//	}
	//}

	++m_curIter;

	std::cout << "curPopSize\t" << m_curPop.size() << std::endl;

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


	std::vector<bool> activeSons(sonSols.size(), false);

	for (int idx(0); idx < sonSols.size(); ++idx) {
		tCross->setParents(*sonSols[idx], *m_curPop[shuffleIds[idx]], fFlagC, Nch);
		tCross->doIt(*sonSols[idx], *m_curPop[shuffleIds[idx]], Nch, 1, fFlagC, fEdgeFreq);
		updateFitness(sonSols[idx].get());
	}

	
	std::vector<Indi*> sonIndiPtr(sonSols.size());
	std::vector<int> cursolx;
	for (int idx(0); idx < sonSols.size(); ++idx) {
		sonSols[idx]->transferSol(cursolx);
		auto hash = m_solHash.calHash(cursolx);
		auto& indi = m_solMap[hash];
		if (indi == nullptr) {
			m_indis.emplace_back(std::move(sonSols[idx]));
			m_indis.back()->setSolId(m_indis.size() - 1);
			m_indis.back()->setRndId(hash);
			indi = m_indis.back().get();
			indi->m_curClassFlag = m_curEvolveClass;
			activeSons[idx] = true;
		}
		else if (indi->m_curClassFlag != m_curEvolveClass) {
			indi->m_curClassFlag = m_curEvolveClass;
			activeSons[idx] = true;
		}
		sonIndiPtr[idx] = indi;
	}

	{
		std::vector<Indi*> neighbors;
		getNeighbors(m_curPop, neighbors);
		for (auto& it : m_curPop) {
			neighbors.push_back(it);
		}

		updateRelationShipOf2V(sonIndiPtr, neighbors, rnd);
		updateRelationShipOfV(sonIndiPtr, rnd);
		updateRoot(sonIndiPtr);

	}


	if (m_expanedStage && m_curPop.size() == m_numPop) {
		m_expanedStage = false;
	}
	
	if (m_expanedStage) {
		std::vector<Indi*> oldPop;
		std::vector<Indi*> newPop;
		for (int idx(0); idx < sonIndiPtr.size(); ++idx) {
			if (activeSons[idx]) {
				oldPop.push_back(m_curPop[idx]);
				newPop.push_back(sonIndiPtr[idx]);
			}
			else {
				newPop.push_back(m_curPop[idx]);
			}
		}

		int leftSize = std::min(m_numPop - newPop.size(), oldPop.size());
		if (leftSize == 0) {
			m_expanedStage = false;
			
		}
		rnd->uniform.shuffle(oldPop.begin(), oldPop.end());
		for (int idx(0); idx < leftSize; ++idx) {
			newPop.push_back(oldPop[idx]);
		}
		swap(m_curPop, newPop);




		std::vector< std::vector<std::vector<int>>*> links;
		for (auto& it : m_curPop) {
			links.push_back(&it->fLink);
		}
		calEdgeFreq(links);

	}
	else {

		std::vector<Indi*> newPop;
		for (int idx(0); idx < sonIndiPtr.size(); ++idx) {
			if (activeSons[idx]) {
				newPop.push_back(sonIndiPtr[idx]);
			}
			else {
				newPop.push_back(m_curPop[idx]);
			}
		}
		swap(m_curPop, newPop);

		double maxFit = m_curPop.front()->fitness();

		for (auto& it : m_curPop) {
			maxFit = std::max(maxFit, it->fitness());
		}

		bool stop = false;
		int solId = 0;
		for (auto& it : m_curPop) {
			if (it->m_stagnation_time >= m_maxStag && it->fitness() == maxFit) {
				solId = it->m_solId;
				stop = true;
				break;
			}
		}

		if (stop) {
			if (m_peakAround != nullptr) {
				std::cout << "converge solution to init peak\n";
				std::cout << m_indis[solId]->distanceTo(*m_peakAround) << std::endl;
				std::cout << m_indis[solId]->fitness() - m_peakAround->fitness() << std::endl;

				std::cout << "dis to opt\t" << m_indis[solId]->distanceTo(m_bestSols.front()) << std::endl;
				std::cout << "end" << std::endl;
			}

			std::sort(m_curPop.begin(), m_curPop.end(), [&](Indi* a, Indi* b) {
				return a->m_dis2parent > b->m_dis2parent;
			});
			
			std::vector<double> disVec;
			for (auto& it : m_curPop) {
				if(it->m_dis2parent<1e9)
				disVec.push_back(it->m_dis2parent);
			}

			double meanDis(0), stdDis(0);
			calMeanAndStd(disVec, meanDis, stdDis);

			{
				std::cout << "analysis curpop" << std::endl;
				std::cout << "meanDis\t" << meanDis << "stdDis\t" << stdDis << std::endl;
				std::cout << "end" << std::endl;
			}

			
			double filterDis = meanDis + stdDis * 3;
			
			std::vector<Indi*> filterData;
			for (auto& it : m_curPop) {
				if (it->m_dis2parent > filterDis) {
					filterData.push_back(it);
				}
			}

			for (auto& it : filterData) {
				m_peaks.push_back(it);
			}

			
			
			//m_peaks.push_back(m_indis[solId].get());
			m_curPop.clear();
		}



		for (auto& it : m_curPop) {
			++it->m_stagnation_time;
		}

		//bool firstBest = true;
		//
		//std::vector<Indi*> filterIndis;
		//for (auto& it : m_curPop) {
		//	if (it->m_stagnation_time >= m_maxStag) {
		//		if (firstBest && it->fitness() == maxFit) {
		//			firstBest = false;
		//		}
		//		else {
		//			filterIndis.push_back(it);
		//			//m_stagnationIndis.push_back(it);
		//		}
		//	}
		//	//else newPop.push_back(it);
		//}

		//for (auto& it : filterIndis) {
		//	if (it->m_curClassFlag == m_curEvolveClass) {
		//		increaseClassFlag();
		//		setBasinsCurClass(it,m_maxCurClassFlag);
		//	}
		//}

		//newPop.clear();

		//for (auto& it : m_curPop) {
		//	if (it->m_curClassFlag == m_curEvolveClass) {
		//		newPop.push_back(it);
		//	}
		//	//else newPop.push_back(it);
		//}
	
		//swap(m_curPop, newPop);



		//for (int idx(0); idx < sonIndiPtr.size(); ++idx) {
		//	if (activeSons[idx] ||
		//		sonIndiPtr[idx]->m_randId != m_curPop[idx]->m_randId) {
		//		//double dis = sonIndiPtr[idx]->distanceTo(*m_curPop[idx]);
		//		m_curPop[idx]->m_improve = true;
		//	}
		//}
		//for (auto& it : m_curPop) {
		//	++it->m_stagnation_time;
		//}

		//{


		//	double bestFit = m_curPop.front()->fitness();
		//	for (auto& it : m_curPop) {
		//		bestFit = std::max(it->fitness(), bestFit);
		//	}


		//	for (int idx(0); idx < activeSons.size(); ++idx) {
		//		if (activeSons[idx]) {
		//			bestFit = std::max(sonIndiPtr[idx]->fitness(), bestFit);
		//		}
		//	}

		//	std::vector<Indi*> newIndi;
		//	for (int idx(0); idx < activeSons.size(); ++idx) {
		//		if (activeSons[idx]) {
		//			if (sonIndiPtr[idx]->fitness() == bestFit) {
		//				newIndi.push_back(sonIndiPtr[idx]);
		//			}
		//			else {
		//				auto parent = sonIndiPtr[idx]->m_parentNode->m_cur;
		//				if (parent->m_curClassFlag == m_curClassFlag) {
		//					newIndi.push_back(sonIndiPtr[idx]);
		//				}
		//			}
		//		}
		//	}

		//	std::vector<Indi*> newCurPop;
		//	for (int idx(0); idx < m_curPop.size(); ++idx) {
		//		if (!m_curPop[idx]->m_improve) {
		//			newCurPop.push_back(m_curPop[idx]);
		//		}
		//	}
		//	for (auto& it : newIndi) {
		//		newCurPop.push_back(it);
		//	}
		//	swap(newCurPop, m_curPop);
		//}


		//{
		//	std::vector<Indi*> newCurPop;

		//	bool firstFlag = true;
		//	double maxFit = 0;
		//	for (auto& it : m_curPop) {
		//		maxFit = std::max(maxFit, it->fitness());
		//	}

		//	for (auto& it : m_curPop) {

		//		if (it->m_stagnation_time >= m_maxStag) {
		//			if (it->fitness() == maxFit && firstFlag) {
		//				newCurPop.push_back(it);
		//			}
		//			else {
		//				m_stagnationIndis.push_back(it);
		//			}
		//		}
		//		else newCurPop.push_back(it);

		//		if (it->m_parentNode->getParent() == nullptr) {
		//			firstFlag = false;
		//		}
		//	}
		//	swap(newCurPop, m_curPop);
		//}


		if (m_curPop.size() ==1) {

			//m_stagnationIndis.push_back(m_curPop.back());
			m_curPop.back()->m_exploitationFlag = true;
			m_peaks.push_back(m_curPop.back());
			m_curPop.pop_back();
		}
	}

	{
		std::cout << "curIter\t" << m_curIter << std::endl;
		std::vector<int> belong;
		std::vector<double> dis2parent;
		std::vector<double> fitness;
		std::vector<int> optFlag;
		std::vector<double> curoptFit;
		calNBN(belong, dis2parent, fitness, optFlag, curoptFit,
			rnd);


		analysis(optFlag);
	}
	return EvaluationTag::kNormalEval;
}


int ofec::PopNBN_EAX_V4::evolve3(Problem* pro, Algorithm* alg, Random* rnd) {
	++m_curIter;

	std::cout << "curPopSize\t" << m_curPop.size() << std::endl;

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


	for (int idx(0); idx < sonSols.size(); ++idx) {
		tCross->setParents(*sonSols[idx], *m_curPop[shuffleIds[idx]], fFlagC, Nch);
		tCross->doIt(*sonSols[idx], *m_curPop[shuffleIds[idx]], Nch, 1, fFlagC, fEdgeFreq);
		updateFitness(sonSols[idx].get());
	}


	std::vector<Indi*> sonIndiPtr;
	std::vector<int> cursolx;
	for (int idx(0); idx < sonSols.size(); ++idx) {
		sonSols[idx]->transferSol(cursolx);
		auto hash = m_solHash.calHash(cursolx);
		auto& indi = m_solMap[hash];
		if (indi == nullptr) {
			m_indis.emplace_back(std::move(sonSols[idx]));
			m_indis.back()->setSolId(m_indis.size() - 1);
			m_indis.back()->setRndId(hash);
			indi = m_indis.back().get();
			indi->m_curClassFlag = m_curEvolveClass;
			sonIndiPtr.push_back(indi);
		}
		//	sonIndiPtr[idx] = indi;
	}

	{
		std::vector<Indi*> neighbors;
		getNeighbors(m_curPop, neighbors);
		for (auto& it : m_curPop) {
			neighbors.push_back(it);
		}

		updateRelationShipOf2V(sonIndiPtr, neighbors, rnd);
		updateRelationShipOfV(sonIndiPtr, rnd);
		updateRoot(sonIndiPtr);
	}

	for (auto& it : sonIndiPtr) {
		it->m_originDis = it->m_dis2parent;
		m_curPop.push_back(it);
	}

	for (auto& it : m_curPop) {
		++it->m_stagnation_time;
	}

	sonIndiPtr.clear();
	for (auto& it : m_curPop) {
		if (it->m_originDis == it->m_dis2parent) {
			if (it->m_stagnation_time >= 30) {
				m_stagnationIndis.push_back(it);
			}
			else
				sonIndiPtr.push_back(it);
		}
	}

	swap(m_curPop, sonIndiPtr);


	return 0;
}

int ofec::PopNBN_EAX_V4::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
	return evolve2(pro, alg, rnd);
}


void ofec::PopNBN_EAX_V4::initialize(Problem* pro, Random* rnd) {

	using namespace ofec;
	using namespace eax_tsp;

	m_curIter = 0;
	Npop = 100;
	Nch = 30;
	define(pro, rnd->getSharedPtr());
	init();



	m_lkh_alg.setDefaultParaments();
	{

		auto& tsp_pro = dynamic_cast<TravellingSalesman&>(*pro);
	/*	auto& v = pro->getParam();
	    std::string file_name = v->get<std::string>("dataFile1");
		std::string filedir;
		if (v->has("dataFile2")) {
			filedir = v->get<std::string>("dataFile2");
		}
		else {
			filedir = "data";
		}
		std::ostringstream oss1;
		oss1 << static_cast<std::string>(g_working_dir);
		std::cout << "filename\t" << file_name << std::endl;
		oss1 << "/instance/problem/combination/travelling_salesman/" + filedir + "/" << file_name << ".tsp";
	*/	m_lkh_alg.ReadProblem(tsp_pro.filePath());
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
		std::string filename = "2048_found_sols.txt";
		insertBestSol(foundSolDir + filename, pro);
	}
	m_root.reset(new Indi);
	m_root->initialize(pro->numberVariables(), -1, m_curIter);
	m_pos = (pro->optimizeMode(0) == ofec::OptimizeMode::kMaximize) ? 1 : -1;
	
	m_solHash.initialize(rnd, pro->numberVariables()+10);
	m_solMap.clear();
	m_indis.clear();
	m_maxCurClassFlag = 0;
	m_curVisitedTime = 0;
	m_curEvolveClass = 1;


	std::shared_ptr<Indi> cursol;
	std::vector<int> cursolx;
	while (m_curPop.size() < m_numPop) {
		cursol.reset(new Indi());
		cursol->initialize(pro->numberVariables(), m_curPop.size(), m_curIter);

		tKopt->makeRandSol(*cursol); /* Make a random tour */
		tKopt->doIt(*cursol);        /* Apply the local search with the 2-opt neighborhood */
		updateFitness(cursol.get());
		cursol->transferSol(cursolx);
		auto hash = m_solHash.calHash(cursolx);
		auto& indi = m_solMap[hash];
		if (indi == nullptr) {
			m_indis.emplace_back(std::move(cursol));
			m_indis.back()->setSolId(m_indis.size() - 1);
			m_indis.back()->setRndId(hash);
			m_indis.back()->setCurClass(m_curEvolveClass);
			indi = m_indis.back().get();
			m_curPop.push_back(indi);
		}

	}
	
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
	updateRelationShipOfV(m_curPop, rnd);
	updateRoot(m_curPop);
	std::vector< std::vector<std::vector<int>>*> links;
	for (auto& it : m_curPop) {
		links.push_back(&it->fLink);
	}
	calEdgeFreq(links);

	for (auto& it : m_curPop) {
		it->m_originDis = it->m_dis2parent;
	}
}


void ofec::PopNBN_EAX_V4::analysis(std::vector<int>& optNearIds) {
	std::cout << std::endl << std::endl;
	std::cout << "showInfo" << std::endl;

	std::vector<int> parentIds;
	std::vector<double> dis2opt;
	std::vector<double> dis2parent;
	for (int idx(0); idx < optNearIds.size(); ++idx) {
		parentIds.clear();
		int curId = optNearIds[idx];
		parentIds.push_back(curId);
		while (m_indis[curId]->m_parentNode->getParent()!=nullptr) {
			curId = m_indis[curId]->m_parentNode->getParent()->m_solId;
			parentIds.push_back(curId);
		}
		dis2opt.clear();
		dis2parent.clear();
		for (auto& it : parentIds) {
			dis2opt.push_back(m_indis[it]->distanceTo(m_bestSols[idx]));
			dis2parent.push_back(m_indis[it]->m_dis2parent);
		}
		std::cout << "The \t" << idx << "\tsolution" << std::endl;
		
		for (auto& it : parentIds) {
			std::cout << it << "\t";
		}
		std::cout << std::endl;

		for (auto& it : dis2opt) {
			std::cout << it << "\t";
		}
		std::cout << std::endl;
		for (auto& it : dis2parent) {
			std::cout << it << "\t";
		}
		std::cout << std::endl;
		
	//	std::cout << std::endl;
	}
}



void ofec::PopNBN_EAX_V4::getPopIds(std::vector<int>& popIds) {
	popIds.clear();
	for (auto& it : m_curPop) {
		popIds.push_back(it->m_solId);
	}
}

void ofec::PopNBN_EAX_V4::calNBN(
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	std::vector<double>& fitness,
	std::vector<int>& optFlag,
	std::vector<double>& curoptFit,
	ofec::Random* rnd
) {
	auto& totalSols(m_indis);

	std::vector<int> id2curIdxs(totalSols.back()->m_solId + 1, -1);
	for (int idx(0); idx < totalSols.size(); ++idx) {
		id2curIdxs[totalSols[idx]->m_solId] = idx;
	}

	belong.resize(totalSols.size());
	dis2parent.resize(totalSols.size());
	fitness.resize(totalSols.size());
	for (int idx(0); idx < totalSols.size(); ++idx) {
		if (totalSols[idx]->m_parentNode->getParent() ==nullptr) {
			belong[idx] = id2curIdxs[totalSols[idx]->m_solId];
		}
		else belong[idx] = id2curIdxs[totalSols[idx]->m_parentNode->getParent()->m_solId];
		dis2parent[idx] = totalSols[idx]->m_dis2parent;
		fitness[idx] = totalSols[idx]->fitness();
	}

	//m_net.getGraph(root, belong, dis2parent, fitness, sols);
	if (!ofec::NBN_FLA_info::judgeNBN(belong, fitness)) {
		std::cout << "nbn error" << std::endl;
	}




	std::vector<double> optDis(m_bestSols.size(), std::numeric_limits<double>::max());
	std::vector<double> optFit(m_bestSols.size(), -std::numeric_limits<double>::max());
	optFlag.resize(m_bestSols.size(), -1);

	curoptFit.resize(m_bestSols.size());
	//auto m_pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	for (int idx(0); idx < optFit.size(); ++idx) {
		curoptFit[idx] = m_bestSols[idx].fEvaluationValue * -1;
	}
	for (int idx(0); idx < m_bestSols.size(); ++idx) {
		auto& curMinDis = optDis[idx];
		auto& curFit = optFit[idx];
		auto& cursol = m_bestSols[idx];
		auto& id = optFlag[idx];
		for (int idy(0); idy < totalSols.size(); ++idy) {
			auto& other = totalSols[idy];
			double curdis = cursol.distanceTo(*other);
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

	for (auto& it : optFlag) {
		std::cout << it << "\t";
	}
	std::cout << std::endl;
	for (auto& it : optDis) {
		std::cout << it << "\t";
	}	
	std::cout << std::endl;

	for (int idx(0); idx < optFit.size(); ++idx) {
		std::cout << m_bestSols[idx].fEvaluationValue << "\t";
	}

	std::cout << std::endl;


	for (int idx(0); idx < optFit.size(); ++idx) {
		std::cout << m_indis[optFlag[idx]]->fEvaluationValue << "\t";
	}

	std::cout << std::endl;
}