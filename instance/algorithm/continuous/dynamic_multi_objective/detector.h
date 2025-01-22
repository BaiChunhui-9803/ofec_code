#ifndef OFEC_DETECTOR_H
#define OFEC_DETECTOR_H

#include "../../../../core/algorithm/algorithm.h"
#include "../../../../core/algorithm/population.h"
#include "../../../../core/problem/solution.h"
#include "../../../../core/problem/continuous/continuous.h"
#include "../../template/selection/multi_objective/nsgaii.h"

//detect changes in dynamic environments
namespace ofec {
	template <typename TPopulation>
	bool ifProblemChanged(const TPopulation& pop, Real ratio,Problem *pro,Algorithm *alg, Random *rnd) {
		size_t sample_num = std::floor(pop.size() * ratio);//get the number of re-evaluate Solutions
		std::vector<Solution<>*> sample;//store re-evaluate Solutions
		std::vector<size_t> index;//store the index of the re-evaluate Solutions
		//int tag = kNormalEval;
		for (size_t i = 0; i < sample_num;) {
			size_t ind = std::floor(pop.size() * rnd->uniform.next());//index of the random Solution
			//ensure no repeat inidvidual
			if (index.empty() || find(index.begin(), index.end(), ind) == index.end()) {
				index.push_back(ind);
				Solution<> temp_sol(pro->numberObjectives(), pro->numberConstraints());
				temp_sol.variable().vect()= pop[ind].variable().vect();
				temp_sol.evaluate(pro,alg);//re-evaluate and is invalid
				sample.emplace_back(&temp_sol);
				for (size_t j = 0; j <CAST_CONOP(pro)->numberObjectives(); j++) {//detect whether the objectives have been changed or not
					if (temp_sol.objective()[j] != pop[ind].objective()[j])
						return true;
				}
				i++;
			}
		}
		return false;
	}

	template <typename TPopulation>
	bool ifPopConverged(const TPopulation& pop) {
		for (size_t i = 0; i < pop.size(); i++) {
			if (pop[i]->fitness() == 0) {
				if (i == pop.size() - 1) {
					return true;
				}
			}
			else {
				return false;
			}
		}
	}

	/*template <typename T>
	bool if_pop_stagnation(std::vector<std::unique_ptr<Solution>>& pop) {
		for (size_t i = 0; i < pop.size(); i++) {
			if (pop[i]->fitness() == 0) {
				if (i == pop.size() - 1) {
					return true;
				}
			}
			else {
				return false;
			}
		}
	}*/

	//template <typename T>
	//void update_var_dim(std::vector<std::unique_ptr<T>>& pop) {
	//	std::vector<int> variable_flag = CONTINUOUS_CAST->get_variable_flag();
	//	for (size_t i = 0; i < pop.size(); i++) {
	//		std::vector<real> temp1 = pop[i]->variable().vect();
	//		std::vector<real> temp2;
	//		for (size_t j = 0, m = 0; j < variable_flag.size(); j++) {
	//			if (variable_flag[j] == 0) {
	//				temp2.push_back(temp1[m]);
	//				m++;
	//			}
	//			else if (variable_flag[j] == 1) {
	//				real tmp = CONTINUOUS_CAST->range()[m].limit.second - CONTINUOUS_CAST->range()[m].limit.first;
	//				temp2.push_back(CONTINUOUS_CAST->range()[0].limit.first + tmp * global::ms_global->m_uniform[caller::Problem]->next());
	//			}
	//			else {
	//				m++;
	//			}
	//		}
	//		pop[i]->variable().vect() = temp2;
	//	}

	//	/*for (size_t i = 0; i < m_offspring.size(); i++) {
	//		std::vector<real> temp1 = m_offspring[i].variable().vect();
	//		std::vector<real> temp2;
	//		for (size_t j = 0, m = 0; j < variable_flag.size(); j++) {
	//			if (variable_flag[j] == 0) {
	//				temp2.push_back(temp1[m]);
	//				m++;
	//			}
	//			else if (variable_flag[j] == 1) {
	//				real tmp = CONTINUOUS_CAST->range()[m].limit.second - CONTINUOUS_CAST->range()[m].limit.first;
	//				temp2.push_back(CONTINUOUS_CAST->range()[0].limit.first + tmp * global::ms_global->m_uniform[caller::Problem]->next());
	//			}
	//			else {
	//				m++;
	//			}
	//		}
	//		m_offspring[i].variable().vect() = temp2;
	//	}*/
	//}

	//template <typename T>
	//void update_obj_dim(std::vector<std::unique_ptr<T>>& pop) {
	//	size_t num_obj = pop[0]->objective_size();
	//	size_t new_num_obj = global::ms_global->m_problem->objective_size();
	//	for (size_t i = 0; i < pop.size(); i++) {
	//		if (new_num_obj > num_obj) {
	//			for (size_t j = 0; j < new_num_obj - num_obj; j++) {
	//				pop[i]->objective().push_back(1);
	//			}
	//		}
	//		else {
	//			for (size_t j = 0; j < num_obj - new_num_obj; j++) {
	//				pop[i]->objective().pop_back();
	//			}
	//		}
	//	}
	//	/*for (size_t i = 0; i < m_offspring.size(); i++) {
	//		if (new_num_obj > num_obj) {
	//			for (size_t j = 0; j < new_num_obj - num_obj; j++) {
	//				m_offspring[i].objective().push_back(1);
	//			}
	//		}
	//		else {
	//			for (size_t j = 0; j < num_obj - new_num_obj; j++) {
	//				m_offspring[i].objective().pop_back();
	//			}
	//		}
	//	}*/
	//}
}

#endif // !OFEC_DETECTOR_H