/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Li Zhou
* Email: 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/
// updated Jan 19, 2019 by Li Zhou

#ifndef OFEC_DYNSAMULTIPOP
#define OFEC_DYNSAMULTIPOP

#include "../../../../core/algorithm/multi_population.h"
#include "../../../problem/realworld/csiwdn/csiwdn.h"
#include "../../../../core/exception.h"

namespace ofec {

	template <typename TPopulation>
	class DynSaMultiPop : public MultiPopulation<TPopulation>
	{
		using PAIR = std::pair<VarCSIWDN, Real>;
		struct CmpByValue {
			bool operator()(const PAIR& lhs, const PAIR& rhs) {
				return lhs.second < rhs.second;
			}
		};
	public:
		DynSaMultiPop() = default;
		DynSaMultiPop(size_t n, size_t subsize,Environment* env): MultiPopulation<TPopulation>(n,subsize, env){}
		void populationCompetition(size_t indx1, size_t indx2);
		bool isHighCoverageRate(Environment* env, size_t indx1, size_t indx2, double maxrate);
		bool isStable(Environment* env);
		void addNewPopulation(Environment* env, Random *rnd, size_t num_pop, size_t size_pop, CSIWDN::InitType type, const std::pair<int, int>& source_index);
	protected:
		std::vector<size_t> m_count_last;
	};

	template <typename TPopulation>
	void DynSaMultiPop<TPopulation>::populationCompetition(size_t indx1, size_t indx2) {
		
		size_t sub_size1 = MultiPopulation<TPopulation>::m_populations[0]->m_first_pop->size();
		std::vector<std::pair<VarCSIWDN, Real>> temp1;
		for (size_t i = 0; i < sub_size1; ++i)
			temp1.push_back(std::make_pair((*(MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop))[i].variable(), (*(MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop))[i].objective(0)));
		for (size_t i = 0; i < sub_size1; ++i)
			temp1.push_back(std::make_pair((*(MultiPopulation<TPopulation>::m_populations[indx2]->m_first_pop))[i].variable(), (*(MultiPopulation<TPopulation>::m_populations[indx2]->m_first_pop))[i].objective(0)));
		sort(temp1.begin(), temp1.end(), CmpByValue());
		for (size_t i = 0; i < sub_size1; ++i) {
			(*(MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop))[i].variable() = temp1[i].first;
			(*(MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop))[i].objective()[0] = temp1[i].second;
		}

		size_t sub_size2 = MultiPopulation<TPopulation>::m_populations[0]->m_second_pop->size();
		std::vector<std::pair<VarCSIWDN, Real>> temp2;
		for (size_t i = 0; i < sub_size2; ++i)
			temp2.push_back(std::make_pair((*(MultiPopulation<TPopulation>::m_populations[indx1]->m_second_pop))[i].variable(), (*(MultiPopulation<TPopulation>::m_populations[indx1]->m_second_pop))[i].objective(0)));
		for (size_t i = 0; i < sub_size2; ++i)
			temp2.push_back(std::make_pair((*(MultiPopulation<TPopulation>::m_populations[indx2]->m_second_pop))[i].variable(), (*(MultiPopulation<TPopulation>::m_populations[indx2]->m_second_pop))[i].objective(0)));
		sort(temp2.begin(), temp2.end(), CmpByValue());
		for (size_t i = 0; i < sub_size2; ++i) {
			(*(MultiPopulation<TPopulation>::m_populations[indx1]->m_second_pop))[i].variable() = temp2[i].first;
			(*(MultiPopulation<TPopulation>::m_populations[indx1]->m_second_pop))[i].objective()[0] = temp2[i].second;
		}

	}

	template <typename TPopulation>
	bool DynSaMultiPop<TPopulation>::isHighCoverageRate(Environment* env,size_t indx1, size_t indx2, double maxrate) {
		auto pro = env->problem();
		int size1 = MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop->size();
		int size2 = MultiPopulation<TPopulation>::m_populations[indx2]->m_first_pop->size();
		std::vector<bool> flag(size2, true);
		int count = 0;
		for (size_t i = 0; i < size1; ++i) {
			for (size_t j = 0; j < size2; ++j) {
				int z = 0;
				while ((z + 1) < CAST_CSIWDN(pro)->numberSource() && CAST_CSIWDN(pro)->phase() >= (*(MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop))[i].variable().startTime(z + 1)) {
					z++;
				}
				if (flag[j] && (*(MultiPopulation<TPopulation>::m_populations[indx1]->m_first_pop))[i].variable().index(z) == (*(MultiPopulation<TPopulation>::m_populations[indx2]->m_first_pop))[j].variable().index(z)) {
					++count;
					flag[j] = false;
					break;
				}
			}
		}
		double rate = (double)count / (double)size1;
		return rate > maxrate;
	}

	template <typename TPopulation>
	bool DynSaMultiPop<TPopulation>::isStable(Environment* env) {
		auto pro = env->problem();
		size_t node_size = CAST_CSIWDN(pro)->numberNode();
		std::vector<size_t> count(node_size, 0);
		int num_pop = MultiPopulation<TPopulation>::size();
		for (size_t i = 0; i < num_pop; ++i) {
			std::vector<bool> flag(node_size, false);
			for (size_t j = 0; j < MultiPopulation<TPopulation>::m_populations[i]->m_first_pop->size(); ++j) {
				int z = 0;
				while ((z + 1) < CAST_CSIWDN(pro)->numSource() && CAST_CSIWDN(pro)->phase() >= (*(MultiPopulation<TPopulation>::m_populations[i]->m_first_pop))[j].variable().startTime(z + 1)) {
					z++;
				}
				size_t temp = (*(MultiPopulation<TPopulation>::m_populations[i]->m_first_pop))[j].variable().index(z) - 1;
				if (!flag[temp]) {
					++count[temp];
					flag[temp] = true;
				}
			}
		}
		
		if (m_count_last.size() == 0) {
			m_count_last = count;
			return false;
		}
		for (size_t i = 0; i < count.size(); ++i) {
			if (m_count_last[i] != count[i]) {
				m_count_last = count;
				return false;
			}

		}
		return true;


	
	}

	template <typename TPopulation>
	void DynSaMultiPop<TPopulation>::addNewPopulation(Environment* env, Random *rnd,size_t num_pop, size_t size_pop, CSIWDN::InitType type, const std::pair<int, int> & source_index) {
		auto pro = env->problem();
		size_t flag = MultiPopulation<TPopulation>::size();
		for (size_t i = 0; i < num_pop; ++i)
			MultiPopulation<TPopulation>::m_populations.emplace_back(std::move(std::unique_ptr<TPopulation>(new TPopulation(size_pop,env))));
		for (size_t i = 0; i < flag; i++) {
			this->m_populations[i]->updateBest(env);
		}
		Real best_obj = this->m_populations[0]->m_best->objective()[0];
		size_t idx_pop = 0;
		for (size_t i = 1; i < flag; ++i) {
			if (best_obj > this->m_populations[i]->m_best->objective()[0]) {
				best_obj = this->m_populations[i]->m_best->objective()[0];
				idx_pop = i;
			}
		}
		auto &var_best = this->m_populations[idx_pop]->m_best->variable();
		CAST_CSIWDN(pro)->setInitType(type);
		if (type == CSIWDN::InitType::kRandom || type == CSIWDN::InitType::kDistance) {
			for (size_t i = flag; i < flag + num_pop; ++i) {
				MultiPopulation<TPopulation>::m_populations[i]->initializeNewPop(env, rnd, var_best, source_index);
			}
		}
		else if (type == CSIWDN::InitType::kBeVisited) {
			CAST_CSIWDN(pro)->clustering(num_pop, rnd);
			CAST_CSIWDN(pro)->calProByBeVisited();

			//Real sum = 0.0;
			//for (size_t i = 0; i < flag; ++i) {
			//	sum += this->m_populations[i]->m_best->objective()[0];
			//}
			//std::vector<Real> roulette(flag + 1, 0);
			//for (size_t i = 0; i < roulette.size() - 1; ++i) {
			//	roulette[i + 1] = roulette[i] + (1 - this->m_populations[i]->m_best->objective()[0] / sum);
			//}
			for (size_t i = flag; i < flag + num_pop; ++i) {
				CAST_CSIWDN(pro)->setPopIdentifier(i - flag);
				//Real random_value = rnd->uniform.nextNonStd<Real>(roulette[0], roulette[roulette.size() - 1]);
				//for (size_t j = 0; j < roulette.size() - 1; ++j) {
				//	if (random_value >= roulette[j] && random_value < roulette[j + 1]) {
				//		idx_pop = j;
				//		break;
				//	}
				//}
				//auto &var_best = this->m_populations[idx_pop]->m_best->variable();
				MultiPopulation<TPopulation>::m_populations[i]->initializeNewPop(env, rnd, var_best, source_index);
			}
		}
		else if (type == CSIWDN::InitType::kCluster) {
			CAST_CSIWDN(pro)->clustering(num_pop, rnd);
			for (size_t i = flag; i < flag + num_pop; ++i) {
				CAST_CSIWDN(pro)->setPopIdentifier(i - flag);
				MultiPopulation<TPopulation>::m_populations[i]->initializeNewPop(env, rnd, var_best, source_index);
			}
		}
		else throw Exception("No this initialization type");

	}


}

#endif