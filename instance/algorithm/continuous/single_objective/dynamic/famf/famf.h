#ifndef OFEC_FAMF_H
#define OFEC_FAMF_H

#include "../metrics_dynamic.h"
#include "famf_pop.h"
#include "famf_Solution.h"
#include "hibernated_zone.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../problem/continuous/single_objective/dynamic/uncertianty_continuous.h"
#include "../../../../../record/dynamic/rcr_vec_real_dynamic.h"
#include "../../../../../../utility/clustering/hslh.h"

namespace ofec {
	template<typename TPopulation>
	class FAMF : public MetricsDynamicConOEA {
	protected:
		HibernatedZone m_hibernated_area;
		Real m_convergThreshold = 0.0001;
		Real m_convFactor = 0.005;

	public:
		FAMF() {}
		using IndividualType = typename TPopulation::IndividualType;
		void initialize_() override;
		void record() override{}
		int create_new_swarms(int Solution_num);
		void init_Solutions(int Solution_num);
		void measureMultiPop(bool);
		int removeOverlapping();
		Real getAvgCurRadius();
		void hibernate_pop(TPopulation& curPop){}
		void removeRedundentHiber(){}
		void degrade_explored_areas();
		void update_stagnate_state();
		int evolve();
		int evolve_subpop();

		bool checkIncreaseDiv();
		void increase_diversity();
		void update_num_Solutions();
		void update_memory();
		int total_indi_size();

		bool hiberstart_Flag;

	#ifdef OFEC_DEMO
		void updateBuffer() override{};
	#endif
	protected:
		void run_() override;

		MultiPopulation<TPopulation> m_sub_pop;
		std::vector<IndividualType> m_new_indis;
		int mc_init_indiSize;
		int m_min_num_inds = 4;
		const int mc_max_idxs = 1000;
		const int mc_min_idxs = 10;
		const int min_subpopSize = 1;

		int m_next_indiSize = 0;
		int m_pre_indiSize = 0;
		int m_total_indiSize = 0;
		int m_converged_indiSize = 0;

		int m_pre_subpopSize = 0;
		int m_cur_subpopSize = 0;
		int m_convergedPops = 0;

		bool m_first_update;
		int mc_offPeak, mc_stepIndis;
		int m_increaseDivIter = 0;
		int mc_maxIncreaseDivIter = 50;

		Real m_avgRadius = 0;

		std::vector<std::vector<double>> mv_indis;	//map


	};

	template<typename TPopulation>
	void FAMF<TPopulation>::initialize_(){
		MetricsDynamicConOEA::initialize_();
		auto& v = *m_param;
		mc_init_indiSize = v.has("population size") ? v.get<int>("population size") : 100;
		mc_offPeak = 3;
		mc_stepIndis = 5;
	}

	template<typename TPopulation>
	void FAMF<TPopulation>::measureMultiPop(bool effictive) {
		for (int i(0); i < m_sub_pop.size(); ++i) {
			for (int j(0); j < m_sub_pop[i].size(); j++){
				m_sub_pop[i][j].evaluate(m_problem.get(), this,effictive);
			}
		}
	}

	template<typename TPopulation>
	int FAMF<TPopulation>::removeOverlapping()
	{
		for (unsigned i = 0; i < m_sub_pop.size(); i++) {
			for (unsigned j = i + 1; j < m_sub_pop.size(); j++) {
				Real dist = m_problem->variableDistance(this->m_sub_pop[i].get_center(), this->m_sub_pop[j].get_center());//dist表示子种群i与子种群j中心的距离
				if (dist < this->m_sub_pop[i].getInitialraidus() || dist < this->m_sub_pop[j].getInitialraidus()) {	//如果dist小于种群i或种群j的种群内距离
					int c1 = 0, c2 = 0;
					for (int k = 0; k < this->m_sub_pop[j].size(); k++) { //k<j种群的Popsize时，dist等于种群i的中心到每种群j第k个代表性个体的距离
						dist = m_problem->variableDistance(this->m_sub_pop[i].get_center(), this->m_sub_pop[j][k]);
						if (dist < this->m_sub_pop[i].getInitialraidus()) c1++;//如果dist小于种群i的种群内距离，C1++;C1表示种群j种有多少个个体少于种群i内距离
					}
					for (int k = 0; k < this->m_sub_pop[i].size(); k++) {
						dist = m_problem->variableDistance(this->m_sub_pop[i].get_center(), this->m_sub_pop[i][k]);
						if (dist < this->m_sub_pop[i].getInitialraidus()) c2++;//C2表示种群i种有多少个个体少于种群j内距离
					}
					if (c1 > 0 && c2 > 0) {
						int idx = -1;
						m_sub_pop[i].merge_pop(m_sub_pop[j], 0.7*(m_sub_pop[i].size() + m_sub_pop[j].size()));
						m_sub_pop - j;
						--j;
						return j;//返回被删除种群的index
						}					//					cout << "Delete Pop id is " << idx << endl;
					}
				}
			}
		return -1;
	}

	template<typename TPopulation>
	Real FAMF<TPopulation>::getAvgCurRadius()
	{
		if (this->m_sub_pop.size() == 0) return 0;
		double r = 0;
		int count = 0;
		for (unsigned int j = 0; j < this->m_sub_pop.size(); j++) {
			if (this->m_sub_pop[j].isHibernating()) continue;
			r += this->m_sub_pop[j].getCurrentraidus();
			count++;
		}
		if (count > 0) 	r /= count;
		return r;
	}

	//template<typename TPopulation>
	//inline void FAMF<TPopulation>::hibernate_pop(TPopulation& curPop){
	//	curPop.setHiberFlag(true);
	//	int idx(0);
	//	Real minDis(0);
	//	for (auto& it : curPop.best()) {
	//		m_hibernated_area.get_nearest_optimum(*it, idx, minDis);
	//		if (minDis > CAST_CONOP(m_problem.get())->variableAccuracy()) {
	//			m_hibernated_area.add_optimum(*it);
	//			if (m_hibernated_area[m_hibernated_area.size() - 1].getReady() == false) {
	//				measureMultiPop(true);
	//				m_hibernated_area.erase(m_hibernated_area.size() - 1);
	//			}

	//		}
	//	}
	//}

	//template<typename TPopulation>
	//inline void FAMF<TPopulation>::removeRedundentHiber() {
	//	for (int j(0); j < m_sub_pop.size(); ++j) {
	//		Real minDis(0);
	//		int idx;
	//		std::vector<int> hiberPeak; hiberPeak.resize(m_hibernated_area.size());
	//		if (m_sub_pop[j].isHibernating()) {
	//			m_hibernated_area.get_nearest_optimum(*m_sub_pop[j].best().front(), idx, minDis);
	//			if(idx != -1) hiberPeak[idx]++;
	//			if (idx != -1 && minDis < CAST_CONOP(m_problem.get())->variableAccuracy()) {
	//				if (hiberPeak[idx] != 1) {
	//					m_sub_pop - j;
	//					--j;
	//				}
	//			}
	//		}
	//	}
	//}

	template<typename TPopulation>
	void FAMF<TPopulation>::degrade_explored_areas(){
		for (int i(0); i < m_sub_pop.size(); ++i) {
			m_sub_pop[i].degrade_explored_areas(m_hibernated_area);
		}
	}

	template<typename TPopulation>
	void FAMF<TPopulation>::update_stagnate_state(){
		for (int i(0); i < m_sub_pop.size(); ++i) {
			if (!m_sub_pop[i].isHibernating()) {
				m_sub_pop[i].isStagnant(m_avgRadius);
			}
		}
	}

	template<typename TPopulation>
	int FAMF<TPopulation>::create_new_swarms(int Solution_num){
		init_Solutions(Solution_num);
		HSLH<IndividualType> cluster(m_new_indis);
		std::vector<IndividualType> left_indis;

		cluster.clustering(min_subpopSize);
		for (int i(0); i < cluster.size(); ++i) {
			if (cluster[i].size() > m_min_num_inds) {
				auto a = new TPopulation(cluster[i].size());
				m_sub_pop + (a);
				auto iter = cluster[i].begin();
				for (int j(0); j < cluster[i].size(); ++j, ++iter) {
					m_sub_pop[m_sub_pop.size() - 1][j] = *(iter->second);
				}
				m_sub_pop[m_sub_pop.size() - 1].estimateInitialRadius();
			}
			else {
				auto iter = cluster[i].begin();
				for (int j(0); j < cluster[i].size(); ++j, ++iter) {
					left_indis.push_back(*(iter->second));
				}
			}
			//m_sub_pop.add(newPop);
		}
		swap(m_new_indis, left_indis);
		return cluster.size();
	}

	template<typename TPopulation>
	void FAMF<TPopulation>::init_Solutions(int Solution_num){
		int rf;
		int originSize = m_new_indis.size();
		m_new_indis.resize(Solution_num + originSize);

		for (int i(0); i < m_new_indis.size(); ++i) {
			m_new_indis[i].initialize(i);
			rf = m_new_indis[i].evaluate();
			if (rf == kChangeNextEval) measureMultiPop(true);
		}
	}

	template<typename TPopulation>
	int FAMF<TPopulation>::evolve()
	{
		std::cout << "pop_size: " << m_sub_pop.size() << '\t' << "total_indi_size: "<< total_indi_size() << std::endl;

		int rf = evolve_subpop();
		if (rf != kNormalEval) {
			return rf;
		}
		m_avgRadius = getAvgCurRadius();
		//cout << "AvgRadius: "<< m_avgRadius << endl;

		while(removeOverlapping()!= -1);
		for (int i(0); i < m_sub_pop.size(); ++i) {
			if (m_sub_pop[i].isConverged()) {
				m_sub_pop[i].setHiberFlag(true);
			}
		}

		//峰屏蔽机制，有问题
		//hiberstart_Flag = false;
		//if (hiberstart_Flag == true) {
		//	for (int i(0); i < m_sub_pop.size(); ++i) {
		//		if (m_sub_pop[i].isConverged()) {
		//			hibernate_pop(m_sub_pop[i]);
		//			cout << "Population " << i << " is converged now" << endl;
		//		}
		//	}
		//	removeRedundentHiber();
		//}
		degrade_explored_areas();
		update_stagnate_state();

		if (checkIncreaseDiv()) {
			if(m_problem->hasTag(ProblemTag::kDOP)) for (int i(0); i < m_sub_pop.size();i++) m_sub_pop[i].wakeUp();
			//m_hibernated_area.clear();
			std::cout << "Increase Diversity Now!" << std::endl;

			increase_diversity();
		}
		return rf;
	}

	template<typename TPopulation>
	int FAMF<TPopulation>::evolve_subpop()
	{
		int rf(kNormalEval);
		for (int i(0); i < m_sub_pop.size(); ++i) {
			if (m_sub_pop[i].isHibernating()) continue;
			rf = m_sub_pop[i].evolve();
			if (rf == kChangeNextEval) {
				measureMultiPop(true);
				return rf;
			}
			if (rf != kNormalEval) {
				return rf;
			}

		}
		return rf;
	}

	template<typename TPopulation>
	bool FAMF<TPopulation>::checkIncreaseDiv(){
		if (m_sub_pop.size() == 0)
			return true;
		double avgRadius = 0;
		int count = 0;
		for (unsigned k = 0; k < m_sub_pop.size(); k++) {
			double r = this->getAvgCurRadius();
			if (!this->m_sub_pop[k].isHibernating() && !this->m_sub_pop[k].isStagnant(r)) {
				avgRadius += this->m_sub_pop[k].getCurrentraidus();
				count++;
			}
		}
		if (count > 0)	avgRadius /= count;
		else return true;
		
		if (avgRadius <= this->m_convFactor * GET_DYNCONOP(m_problem.get())->domainArea()) {//
			//cout << "Touch off Method 1:" << ++explode_times << endl;
			return true;
		}
		else return false;
	}

	template<typename TPopulation>
	void FAMF<TPopulation>::increase_diversity(){
		update_num_Solutions();
		int new_idx_size(m_next_indiSize - total_indi_size());
		std::cout << "new idx size\t" << new_idx_size << std::endl;
		if (new_idx_size > 0) {
			create_new_swarms(new_idx_size);
			//measureMultiPop(false);
		}
	}

	template<typename TPopulation>
	void FAMF<TPopulation>::update_num_Solutions(){
		m_cur_subpopSize = m_convergedPops + m_sub_pop.size();
		m_total_indiSize = 0;
		for (size_t k = 0; k < m_sub_pop.size(); k++) {
			m_total_indiSize += m_sub_pop[k].size();
		}

		if (m_first_update) {
			m_pre_subpopSize = m_cur_subpopSize;
			m_pre_indiSize = m_convergedPops + m_total_indiSize;
			m_first_update = false;
		}
		update_memory();
		int indisNum = static_cast<int>(m_random->normal.nextNonStd(mv_indis[m_cur_subpopSize][0], mv_indis[m_cur_subpopSize][1]));
		int dif = m_cur_subpopSize - m_pre_subpopSize;
		double ratio = fabs(dif * 1. / mc_offPeak);
		if (ratio >= 1) {
			m_next_indiSize = indisNum + mc_stepIndis * dif;
		}
		else if (ratio == 0) {
			m_next_indiSize = indisNum;
		}
		else {
			double p = m_random->uniform.next();
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
	void FAMF<TPopulation>::update_memory(){
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
		mv_indis[m_cur_subpopSize][0] = (mean * (mv_indis[m_cur_subpopSize].size() - 3) + m_pre_indiSize) / (mv_indis[m_cur_subpopSize].size() - 2);
		mv_indis[m_cur_subpopSize][1] = 0;
		for (unsigned k = 2; k < mv_indis[m_cur_subpopSize].size(); k++) {
			mv_indis[m_cur_subpopSize][1] += (mv_indis[m_cur_subpopSize][k] - mv_indis[m_cur_subpopSize][0]) * (mv_indis[m_cur_subpopSize][k] - mv_indis[m_cur_subpopSize][0]);
		}
		mv_indis[m_cur_subpopSize][1] = sqrt(mv_indis[m_cur_subpopSize][1] / (mv_indis[m_cur_subpopSize].size() - 2));
	}

	template<typename TPopulation>
	int FAMF<TPopulation>::total_indi_size()
	{
		int total(0);
		for (int i(0); i < m_sub_pop.size(); ++i) {
			total += m_sub_pop[i].size();
		}
		return total;
	}

	template<typename TPopulation>
	void FAMF<TPopulation>::run_(){
		while (!terminating()) {
			evolve();
	#ifdef OFEC_DEMO
			update_buffer();
	#endif
		}
	}
}
#endif