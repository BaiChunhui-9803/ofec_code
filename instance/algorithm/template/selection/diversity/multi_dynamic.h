/* Segura C, Coello C A C, Segredo E, et al.
A novel diversity-based replacement strategy for evolutionary algorithms[J].
IEEE transactions on cybernetics, 2015, 46(12): 3233-3246.*/

#ifndef OFEC_SELECT_DIVER_MULTI_DYNAMIC_H
#define OFEC_SELECT_DIVER_MULTI_DYNAMIC_H

#include "../../../../../utility/functional.h"
#include "../../../../../core/environment/environment.h"

namespace ofec::selection::diversity {
	/*  multi dynamic survivor selection scheme */
	template<typename TPopulation, typename TPopAndOff>
	void multiDynamic(TPopulation &pop, TPopAndOff &pop_and_off, Environment *env, Random *rnd, Real elapsed) {
		std::vector<std::vector<Real>> objs(pop_and_off.size(), std::vector<Real>(2));
		std::vector<OptimizeMode> mode = { env->problem()->optimizeMode(0), OptimizeMode::kMaximize };
		for (size_t i = 0; i < pop_and_off.size(); ++i)
			objs[i][0] = pop_and_off[i].objective(0);
		std::set<size_t> current_members;
		for (size_t i = 0; i < pop_and_off.size(); ++i)
			current_members.insert(i);
		size_t best = 0;
		for (size_t i = 1; i < pop_and_off.size(); ++i) {
			if (dominate(pop_and_off[i], pop_and_off[best], env->problem()->optimizeMode()))
				best = i;
		}
		std::set<size_t> new_pop = { best };
		current_members.erase(best);
		std::vector<std::vector<Real>> dist_mat(pop_and_off.size(), std::vector<Real>(pop_and_off.size(), -1));
		Real DCN;
		Real real_max = std::numeric_limits<Real>::max();
		Real real_min = -std::numeric_limits<Real>::max();
		Real mean_range = 0;
		for (size_t j = 0; j < env->problem()->numberVariables(); ++j)
			mean_range += CAST_CONOP(env->problem())->range(j).second - CAST_CONOP(env->problem())->range(j).first;
		mean_range /= env->problem()->numberVariables();
		Real D_I = mean_range / pop.size();
		Real D = D_I * (1 - elapsed); 
		while (new_pop.size() < pop.size()) {
			std::vector<std::vector<Real> *> data;
			for (size_t i : current_members) {
				DCN = real_max;
				for (size_t j : new_pop) {
					if (dist_mat[i][j] == -1)
						dist_mat[i][j] = dist_mat[j][i] = pop_and_off[i].variableDistance(pop_and_off[j], env);
					if (DCN > dist_mat[i][j])
						DCN = dist_mat[i][j];
				}
				objs[i][1] = DCN;
				if (DCN < D)
					objs[i][0] = env->problem()->optimizeMode(0) == OptimizeMode::kMinimize ? real_max : real_min;
				data.push_back(&objs[i]);
			}
			std::vector<int> rank;
			nd_sort::fastSort(data, rank, mode);
			std::vector<size_t> ND;
			auto iter = rank.begin();
			for (size_t i : current_members) {
				if (*iter++ == 0)
					ND.push_back(i);
			}
			size_t selected = *rnd->uniform.nextElem(ND.begin(), ND.end());
			new_pop.insert(selected);
			current_members.erase(selected);
		}
		auto iter = new_pop.begin();
		for (size_t i = 0; i < pop.size(); ++i)
			pop[i] = pop_and_off[*iter++];
	}
}

#endif // !OFEC_SELECT_DIVER_MULTI_DYNAMIC_H