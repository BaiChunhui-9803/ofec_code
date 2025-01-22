#include "nbn_eax_pop_v3.h"
#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../utility/nbn_visualization/nbn_fla.h"
#include <queue>


int ofec::PopNBN_EAX_V3::evolve(Problem* pro, Algorithm* alg, Random* rnd)
{
	using namespace eax_tsp;


	if (m_curPop.empty() && m_stagnationIndis.empty()) {
		return EvalTag::kNormalEval;
	}

	if (m_curPop.empty()) {
		m_id2curId.resize(m_id2inds.size());
		increaseClassFlag();

		std::queue<int> que;
		que.push(m_stagnationIndis.back());
		m_stagnationIndis.pop_back();
		while (!que.empty()) {
			auto cur = que.front();
			que.pop();
			m_curPop.push_back(cur);
			auto& curSol = m_id2inds[cur];
			m_id2curId[curSol->m_solId] = m_curPop.size() - 1;
			curSol->m_curClassFlag = m_curClassFlag;
			if (m_curPop.size() == m_maxPop)break;
			RelationNode* head = curSol->m_sons.m_head.get();
			auto& curMaxDis = curSol->m_maxRadius;
			while (head->m_after->m_curSolId != -1) {
				head = head->m_after;
				que.push(head->m_curSolId);
				curMaxDis = std::max(m_id2inds[head->m_curSolId]->m_maxRadius, curMaxDis);
			}

		}
		if (m_curPop.size() < m_maxPop) {
			unsigned long long curVisited = 0;
			std::vector<unsigned long long> visited(m_curPop.size(), 0);
			std::vector<double> radiusRate(m_curPop.size(), 1.0);
			std::vector<NBN_indi> newSols;

			int dim = pro->numVariables();
			std::vector<int> cursol;
			NBN_indi curindi;
			curindi.initialize(dim, -1);
			while (m_curPop.size() + newSols.size() < m_maxPop) {
				auto& randIndi = *m_id2inds[m_curPop[rnd->uniform.nextNonStd<int>(0, m_curPop.size())]];

				randIndi.transferSol(cursol);
				int solIdx = m_id2curId[randIndi.m_solId];
				int radius = std::round(randIndi.m_maxRadius * radiusRate[solIdx] / 4.0);
				if (radius) {

					while (radius--) {
						auto a = rnd->uniform.nextNonStd<int>(0, dim);
						auto b = rnd->uniform.nextNonStd<int>(0, dim);
						std::swap(cursol[a], cursol[b]);
					}
					curindi.toCurSol(cursol);
					tKopt->doIt(curindi);        /* Apply the local search with the 2-opt neighborhood */
					m_fitnessUpdate(curindi, pro);


					std::queue<int> que;
					que.push(randIndi.m_solId);
					visited[m_id2curId[randIndi.m_solId]] = curVisited;
					bool flagInside = false;
					++curVisited;
					while (!que.empty()) {
						auto& curNei = m_id2inds[que.front()];
						que.pop();
						auto curDis = curNei->distanceTo(curindi);
						if (curDis <= curNei->m_maxRadius) {
							flagInside = true;
							break;
						}

						RelationNode* head = curNei->m_sons.m_head.get();
						while (head->m_after->m_curSolId != -1) {
							head = head->m_after;
							if (m_id2inds[head->m_curSolId]->m_curClassFlag == m_curClassFlag
								&& visited[m_id2curId[head->m_curSolId]] != curVisited) {
								que.push(head->m_curSolId);
							}
						}

						auto parent = curNei->m_parentNode->m_curSolId;
						if (m_id2inds[parent]->m_curClassFlag == m_curClassFlag
							&& visited[m_id2curId[parent]] != curVisited)
							que.push(parent);


					}
					if (flagInside) {
						newSols.push_back(curindi);
					}

					if (flagInside) {
						newSols.push_back(curindi);
						radiusRate[solIdx] = radiusRate[solIdx] * 0.9 + 1.0 * 0.1;
					}

				}
			}

			for (int idx(0); idx < newSols.size(); ++idx) {
				{
				//	m_id2inds.emplace_back(std::move(newSols[idx]));
					m_id2inds.back()->setSolId(m_curSolId++);
					m_curPop.push_back(m_id2inds.size() - 1);
				}
			}

			//NBN_indi* best = m_curPop.front();
			double maxFit = m_id2inds[m_curPop.front()]->m_fitness;
			for (auto& it : m_curPop) {
				if (maxFit < m_id2inds[it]->m_fitness) {
					maxFit = m_id2inds[it]->m_fitness;
				}
			}
			for (auto& it : m_curPop) {
				if (m_id2inds[it]->m_fitness == maxFit) {
					m_id2inds[it]->clearParent();
				}
			}

			updateRelationShipOfV(m_curPop, rnd);


			std::vector< std::vector<std::vector<int>>*> links;
			for (auto& it : m_curPop) {
				links.push_back(&m_id2inds[it]->fLink);
			}
			calEdgeFreq(links);

		}

	}

	++m_curIter;
	std::vector<int> shuffleIds(m_curPop.size());
	for (int idx(0); idx < shuffleIds.size(); ++idx) {
		shuffleIds[idx] = idx;
	}
	rnd->uniform.shuffle(shuffleIds.begin(), shuffleIds.end());
	std::vector<NBN_indi> sonSols(m_curPop.size());
	for (int idx(0); idx < sonSols.size(); ++idx) {
		//m_id2inds.push_back(NBN_indi());
		//sonSols[idx] = &m_id2inds.back();
		sonSols[idx].initialize(pro->numVariables(), -1);
		sonSols[idx].copySol(*m_id2inds[m_curPop[idx]]);
	}

	for (int idx(0); idx < sonSols.size(); ++idx) {
		tCross->setParents(sonSols[idx], *m_id2inds[m_curPop[shuffleIds[idx]]], fFlagC, Nch);
		tCross->doIt(sonSols[idx], *m_id2inds[m_curPop[shuffleIds[idx]]], Nch, 1, fFlagC, fEdgeFreq);
		m_fitnessUpdate(sonSols[idx], pro);
	}



	std::vector<bool> activeSons(sonSols.size(), false);
	//std::vector<double> disMat(sonSols.size(), -1);
	for (int idx(0); idx < sonSols.size(); ++idx) {
		double dis = sonSols[idx].distanceTo(*m_id2inds[m_curPop[idx]]);
	//	disMat[idx] = dis;
		if (dis > 0) {
			double curdis = sonSols[idx].distanceTo(*m_id2inds[m_curPop[shuffleIds[idx]]]);
			if (curdis > 0) {
				activeSons[idx] = true;
				m_id2inds[m_curPop[idx]]->m_improve = true;
				if (dis >= m_id2inds[m_curPop[idx]]->m_dis2parent) {
					++m_id2inds[m_curPop[idx]]->m_stagnation_time;
					m_id2inds[m_curPop[idx]]->m_stag = true;
				}
			}
			else ++m_id2inds[m_curPop[idx]]->m_stagnation_time;
		}
		else {
			++m_id2inds[m_curPop[idx]]->m_stagnation_time;
		}
	}


	{
		std::vector<int> newIndi;
		for (int idx(0); idx < activeSons.size(); ++idx) {
			if (activeSons[idx]) {
			//	m_id2inds.emplace_back(std::move(sonSols[idx]));
				m_id2inds.back()->setSolId(m_curSolId++);
				newIndi.push_back(m_id2inds.size() - 1);
			}
		}

		//std::sort(newIndi.begin(), newIndi.end(), 
		//	[](NBN_indi*a, NBN_indi* b) {
		//	return a->m_fitness < b->m_fitness;
		//});

		{
			updateRelationShipOf2V(newIndi, m_curPop, rnd);
			updateRelationShipOfV(newIndi, rnd);
		}

		//if (m_curPop.size() <= 12) {
		//	for (int idx(0); idx < m_curPop.size(); ++idx) {
		//		std::cout << m_id2inds[m_curPop[idx]]->m_stagnation_time << "\t";
		//	}
		//	std::cout << std::endl;
		//}
		std::vector<int> newCurPop;
		for (int idx(0); idx < m_curPop.size(); ++idx) {
			if (!m_id2inds[m_curPop[idx]]->m_improve) {
				//		newCurPop.push_back(m_curPop[idx]);
				if (m_id2inds[m_curPop[idx]]->m_parentNode->m_curSolId != -1
					&& m_id2inds[m_curPop[idx]]->m_stagnation_time >= m_maxStag) {
					m_stagnationIndis.push_back(m_curPop[idx]);
				}
				else {
					newCurPop.push_back(m_curPop[idx]);
				}

			}
			else if (m_id2inds[m_curPop[idx]]->m_stag) {
				m_stagnationIndis.push_back(m_curPop[idx]);
			}
			//else if (m_id2inds[m_curPop[idx]]->m_stagnation_time) {
			//	if (m_id2inds[m_curPop[idx]]->m_dis2parent > 5) {
			//		newCurPop.push_back(m_curPop[idx]);
			//	}
			//}

		}
		//for (int idx(sonSols.size()); idx < m_curPop.size(); ++idx) {
		//	newCurPop.push_back(m_curPop[idx]);
		//}

		for (auto& it : newIndi) {
			newCurPop.push_back(it);
		}
		swap(newCurPop, m_curPop);
	}

	{
		std::vector<int> newCurPop;
		for (auto& it : m_curPop) {
			if (m_id2inds[it]->m_parentNode->m_curSolId == -1
				&& m_id2inds[it]->m_stagnation_time >= m_maxStag) {
				m_peaks.push_back(it);
			}
			else newCurPop.push_back(it);
		}

		swap(newCurPop, m_curPop);
	}
	

	if (m_curPop.size() == 1) {
		m_peaks.push_back(m_curPop.front());
		m_curPop.pop_back();
	}

	std::cout << "curIter\t" << m_curIter << std::endl;

	{
		std::vector<int> belong;
		std::vector<double> dis2parent;
		std::vector<double> fitness;
		std::vector<int> optFlag;
		std::vector<double> curoptFit;
		calNBN(belong, dis2parent, fitness, optFlag, curoptFit,
			rnd);
	}
	return ofec::kNormalEval;
}

void ofec::PopNBN_EAX_V3::initialize(Problem* pro, Random* rnd)
{
	using namespace ofec;
	using namespace eax_tsp;

	m_curIter = 0;
	Npop = 100;
	Nch = 30;
	define(pro, rnd->getSharedPtr());
	init();
	m_gl_calculator.initialize(pro);
	m_curPop.resize(m_maxPop);

	m_fitnessUpdate =
		[](NBN_indi& sol, Problem* pro) {
		ofec::Real pos = (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		sol.m_fitness = pos * sol.fEvaluationValue;
	};


	for (int i = 0; i < Npop; ++i)
	{
		m_id2inds.emplace_back(new NBN_indi());
		m_curPop[i] = m_id2inds.size() - 1;
		m_id2inds[m_curPop[i]]->initialize(pro->numVariables(), m_curSolId++);
		tKopt->makeRandSol(*m_id2inds[m_curPop[i]]); /* Make a random tour */
		tKopt->doIt(*m_id2inds[m_curPop[i]]);        /* Apply the local search with the 2-opt neighborhood */
		m_fitnessUpdate(*m_id2inds[m_curPop[i]], pro);
	}

	// calculate nbn
	{

		std::sort(m_curPop.begin(), m_curPop.end(), [&](
			int a, int b) {
			return m_id2inds[a]->m_fitness < m_id2inds[b]->m_fitness;
		});

		for (int idx(0); idx < m_curPop.size(); ++idx) {
			//auto& parent = m_curPop[idx]->m_parent;
			auto& cursol = m_id2inds[m_curPop[idx]];
			auto& curMinDis = m_id2inds[m_curPop[idx]]->m_dis2parent;
			for (int idy(idx + 1); idy < m_curPop.size(); ++idy) {
				auto& othersol = m_id2inds[m_curPop[idy]];
				double curDis = cursol->distanceTo(*othersol);

				if (curDis < curMinDis) {
					curMinDis = curDis;
					cursol->updateParent(*othersol);
				}
				else if (curDis == curMinDis && rnd->uniform.next() < 0.5) {
					cursol->updateParent(*othersol);
				}
			}
		}
	}

	std::vector< std::vector<std::vector<int>>*> links;
	for (auto& it : m_id2inds) {
		links.push_back(&it->fLink);
	}
	calEdgeFreq(links);
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

