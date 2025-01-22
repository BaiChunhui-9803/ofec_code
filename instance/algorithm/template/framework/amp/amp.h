#ifndef AMP_H
#define AMP_H

#include "../../../../../core/algorithm/multi_population.h"
#include "../../../../../core/global.h"
#include "../../../../../utility/clustering/hslh.h"
#include <vector>

namespace ofec {
	template <typename TPopulation>
	class AMP : public MultiPopulation<TPopulation> {
	public:
		using PopulationType = TPopulation;
		using MultiPopulation<TPopulation>::m_pops;
		using typename MultiPopulation<TPopulation>::IteratorType;
		using IndividualType = typename TPopulation::IndividualType;

	protected:
		Real m_avgRadius = 0;
		const int mc_init_indiSize;
		int m_next_indiSize = 0;
		int m_pre_indiSize = 0;
		int m_total_indiSize=0;
		int m_converged_indiSize = 0;
		const int mc_max_idxs, mc_min_idxs;
		int m_min_num_inds = 2;
		std::vector<std::unique_ptr<IndividualType>> m_new_indis;	
		int m_pre_subpopSize=0; 
		int m_cur_subpopSize=0;	
		const int mc_init_subpopSize = 1;
		int m_convergedPops;
		const int mc_offPeak, mc_stepIndis;
		bool m_first_update;
		int m_increaseDivIter = 0;
		const int mc_maxIncreaseDivIter = 50;
		std::vector<std::vector<double>> mv_indis;

		virtual void updateMemory();
		virtual void updateNumIdxs(Random *rnd);
		virtual void initSolutions(int Solution_num, Problem *pro, Algorithm *alg, Random *rnd) = 0;
		virtual int createNewSwarms(int Solution_num, Problem *pro, Algorithm *alg, Random *rnd);
		virtual void degradeExploredAreas(Problem *pro, Algorithm *alg, Random *rnd) = 0;
		virtual void removeRedundentHiber(Problem *pro) = 0;
		virtual void updateStagnateState(Problem *pro) = 0;
		void calAvgRadis();
		size_t totalNumInds() const {
			size_t total(0);
			for (int i(0); i < m_pops.size(); ++i) {
				total += m_pops[i]->size();
			}
			return total;
		}
		virtual void removeOverlapping(Problem *pro);
		virtual bool judgeOverlapping(const TPopulation &p1, const TPopulation &p2, Problem *pro) { return false; }
		virtual bool checkIncreaseDiv(Problem *pro) {
			if (m_increaseDivIter >= mc_maxIncreaseDivIter) {
				m_increaseDivIter = 0;
				return true;
			}
			else {
				return false;
			}	
		};
		virtual void hibernatePop(TPopulation &pop, Problem *pro, Algorithm *alg, Random *rnd) { pop.setActive(false); }
	
		IteratorType merge(IteratorType iter1, IteratorType iter2, Problem *pro) {
			int better(0), worse(0);
			for (auto& it : (*iter1)->m_best) {
				for (auto& it2 : (*iter2)->m_best) {
					auto r = objectiveCompare(it->objective(), it2->objective(), pro->optimizeMode());
					if (r == Dominance::kDominant) {
						better++;
					}
					else if (r == Dominance::kDominated) {
						worse++;
					}
				}
			}
			if (better >= worse) {
				(*iter1)->m_initial_radius = (*iter1)->updateRadius((*iter1)->m_initial_radius, (*iter2)->m_initial_radius);
			}
			else {
				(*iter1)->m_initial_radius = (*iter1)->updateRadius((*iter2)->m_initial_radius, (*iter1)->m_initial_radius);
			}
			for (auto& it : (*iter2)->best()) {
				(*iter1)->updateBest(*it, pro);
			}
			size_t pre_size((*iter1)->size());
			(*iter1)->m_individuals.resize((*iter1)->size() + (*iter2)->size());
			for (size_t i = pre_size; i < (*iter1)->size(); ++i) {
				std::swap((*iter1)->m_individuals[i], (*iter2)->m_individuals[i - pre_size]);
			}
			if ((*iter1)->size() > 10/*mc_max_idxs*/) {
				(*iter1)->sort(pro);
				std::sort((*iter1)->begin(), (*iter1)->end(), [](std::unique_ptr<IndividualType>& a, std::unique_ptr<IndividualType>& b) {
					return a->fitness() < b->fitness();
					});
				(*iter1)->m_individuals.resize(10);
			}
			return m_pops.erase(iter2);
		}
		
	public:
		AMP(size_t size_pop) : 
			mc_init_indiSize(size_pop),
			mc_max_idxs(10000), 
			mc_min_idxs(10), 
			mv_indis(100), 
			mc_offPeak(3), 
			mc_stepIndis(5)
		{
			for (int i = 0; i < 100; i++) {
				mv_indis[i].push_back(0);
				mv_indis[i].push_back(0);
			}
		};
		std::vector<std::unique_ptr<IndividualType>> & remainInds() { return m_new_indis; }
		void initialize(Problem *pro, Algorithm *alg, Random *rnd);
		virtual int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
	};

	template<typename TPopulation>
	void AMP<TPopulation>::updateMemory() {
		if (m_cur_subpopSize >= mv_indis.size()) {
			unsigned size = mv_indis.size();
			mv_indis.resize(m_cur_subpopSize + 1);
			for (unsigned i = size; i < mv_indis.size(); i++) {
				mv_indis[i].push_back(0);
				mv_indis[i].push_back(0);
			}
		}
		mv_indis[m_cur_subpopSize].push_back(m_pre_indiSize);
		double mean = mv_indis[m_cur_subpopSize][0];
		mv_indis[m_cur_subpopSize][0] = (mean*(mv_indis[m_cur_subpopSize].size() - 3) + m_pre_indiSize) / (mv_indis[m_cur_subpopSize].size() - 2);
		mv_indis[m_cur_subpopSize][1] = 0;
		for (unsigned k = 2; k < mv_indis[m_cur_subpopSize].size(); k++) {
			mv_indis[m_cur_subpopSize][1] += (mv_indis[m_cur_subpopSize][k] - mv_indis[m_cur_subpopSize][0])*(mv_indis[m_cur_subpopSize][k] - mv_indis[m_cur_subpopSize][0]);
		}
		mv_indis[m_cur_subpopSize][1] = sqrt(mv_indis[m_cur_subpopSize][1] / (mv_indis[m_cur_subpopSize].size() - 2));
	}

	template<typename TPopulation>
	void AMP<TPopulation>::updateNumIdxs(Random *rnd) {
		m_cur_subpopSize = m_convergedPops + m_pops.size();
		m_total_indiSize = totalNumInds();
		if (m_first_update) {
			m_pre_subpopSize = m_cur_subpopSize;
		    m_pre_indiSize = m_convergedPops + m_total_indiSize;
			m_first_update = false;
		}
		updateMemory();
		int indisNum = static_cast<int>(rnd->normal.nextNonStd(mv_indis[m_cur_subpopSize][0], mv_indis[m_cur_subpopSize][1]));
		int dif = m_cur_subpopSize - m_pre_subpopSize;
		double ratio = fabs(dif*1. / mc_offPeak);
		if (ratio >= 1) {
			m_next_indiSize = indisNum + mc_stepIndis * dif;
		}
		else if (ratio == 0) {
			m_next_indiSize = indisNum;
		}
		else {
			double p = rnd->uniform.next();
			if (p <= ratio) m_next_indiSize = indisNum + mc_stepIndis * dif;
			else m_next_indiSize = indisNum;
		}
		m_pre_subpopSize = m_cur_subpopSize;	
		if (m_next_indiSize <= m_total_indiSize) {
			m_next_indiSize = m_total_indiSize + mc_stepIndis * 2;
		}
		if (m_next_indiSize > mc_max_idxs) m_next_indiSize = mc_max_idxs;
		if (m_next_indiSize < mc_min_idxs) m_next_indiSize = mc_min_idxs;	
	}


	template<typename TPopulation>
	int AMP<TPopulation>::createNewSwarms(int Solution_num, Problem *pro, Algorithm *alg, Random *rnd){
		initSolutions(Solution_num, pro, alg, rnd);
		HSLH<IndividualType> cluster(m_new_indis, pro);
		std::vector<std::unique_ptr<IndividualType>> left_indis;
		cluster.clustering(mc_init_subpopSize, pro);
		for (int i(0); i < cluster.size(); ++i) {
			if (cluster[i].size() > m_min_num_inds) {
				TPopulation a(cluster[i].size(), pro);
				this->append(std::move(a));
				auto iter = cluster[i].begin();
				for (int j(0); j < cluster[i].size(); ++j, ++iter) {
					(*m_pops.back())[j] = *(iter->second);
					(*m_pops.back())[j] = *(iter->second);
				}
				m_pops.back()->initRadius(pro);
			}
			else {
				auto iter = cluster[i].begin();
				for (int j(0); j < cluster[i].size(); ++j, ++iter) {
					left_indis.emplace_back(new IndividualType(*(iter->second)));
				}
			}
		}
		m_new_indis.resize(left_indis.size());
		for (size_t i = 0; i < left_indis.size(); i++) {
			m_new_indis[i].reset(left_indis[i].release());
		}
		return cluster.size();
	}

	template<typename TPopulation>
	void AMP<TPopulation>::calAvgRadis() {
		m_avgRadius = 0;
		int totalPop(0);
		for (int i(0); i < m_pops.size(); ++i) {
			if (m_pops[i]->isActive()) {
				++totalPop;
				m_avgRadius += m_pops[i]->curRadius();
			}
		}
		if(totalPop)
		m_avgRadius /= static_cast<double>(totalPop);
	}

	template<typename TPopulation>
	void AMP<TPopulation>::removeOverlapping(Problem *pro) {
		for (auto iter1 = m_pops.begin(); iter1 != m_pops.end(); ++iter1) {
			for (auto iter2 = iter1 + 1; iter2 != m_pops.end();) {
				if (judgeOverlapping(**iter1, **iter2, pro))
					iter2 = merge(iter1, iter2, pro);
				else
					iter2++;
			}
		}
	}

	template<typename TPopulation>
	int AMP<TPopulation>::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int rf = MultiPopulation<TPopulation>::evolve(pro, alg, rnd);
		if (rf != kNormalEval) {
			return rf;
		}
		calAvgRadis();
		removeOverlapping(pro);
		removeRedundentHiber(pro);
		for (int i(0); i < m_pops.size(); ++i) {
			if (m_pops[i]->judgeConverged()) {
				hibernatePop(*m_pops[i], pro, alg, rnd);
			}
		}
		degradeExploredAreas(pro, alg, rnd);
		updateStagnateState(pro);
		if (checkIncreaseDiv(pro)) {
			updateNumIdxs(rnd);
			for (int i(0); i < m_pops.size(); ++i) {
				m_pops[i]->setActive(true);
			}
			int new_idx_size(m_next_indiSize - totalNumInds());
			if (new_idx_size > 0) {
				createNewSwarms(new_idx_size, pro, alg, rnd);
			}
			//std::cout << "Solution size :\t " << totalNumInds() <<"\t"<<"Pop size \t"<< m_pops.size()<< std::endl;
		}
		return rf;
	}

	template<typename TPopulation>
	void AMP<TPopulation>::initialize(Problem *pro, Algorithm *alg, Random *rnd) {
		m_first_update = true;
		m_total_indiSize = 0;
		m_converged_indiSize = 0;
		m_pre_subpopSize = 0;
		m_cur_subpopSize = 0;
		m_convergedPops = 0;
		m_increaseDivIter = 0;
		createNewSwarms(mc_init_indiSize, pro, alg, rnd);
	}
}


#endif