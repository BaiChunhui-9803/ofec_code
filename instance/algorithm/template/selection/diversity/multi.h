/* C. Segura, C. Coello, E. Segredo, G. Miranda, and C. Leon,
¡°Improving the diversity preservation of multi-objective approaches used for single-objective optimization,¡±
in 2013 IEEE Congress on Evolutionary Computation (CEC), June 2013, pp. 3198¨C3205*/

#ifndef OFEC_SELECT_DIVER_MULTI_H
#define OFEC_SELECT_DIVER_MULTI_H

#include <list>
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../core/problem/problem.h"
#include "../../../../../utility/functional.h"

namespace ofec::selection::diversity {
	/* multi survivor selection scheme */
	template<typename TPopulation, typename TPopAndOff>
	void multi(TPopulation &pop, TPopAndOff &pop_and_off, Environment *env, Random *rnd) {
		std::vector<std::vector<Real>> objs(pop_and_off.size(), std::vector<Real>(2));
		std::vector<OptimizeMode> mode = { env->problem()->optimizeMode(0), OptimizeMode::kMaximize };
		for (size_t i = 0; i < pop_and_off.size(); ++i)
			objs[i][0] = pop_and_off[i].objective(0);
		std::list<size_t> current_members;
		for (size_t i = 0; i < pop_and_off.size(); ++i)
			current_members.push_back(i);
		size_t best = 0;
		for (size_t i = 1; i < pop_and_off.size(); ++i) {
			if (dominate(pop_and_off[i], pop_and_off[best], env->problem()->optimizeMode()))
				best = i;
		}
		std::list<size_t> new_pop = { best };
		current_members.erase(std::find(current_members.begin(), current_members.end(), best));
		std::vector<std::vector<Real>> dist_mat(pop_and_off.size(), std::vector<Real>(pop_and_off.size(), -1));
		Real DCN;
		while (new_pop.size() < pop.size()) {	
			std::vector<std::vector<Real>*> data;
			for (size_t i : current_members) {
				DCN = std::numeric_limits<Real>::max();
				for (size_t j : new_pop) {
					if (dist_mat[i][j] == -1)
						dist_mat[i][j] = dist_mat[j][i] = pop_and_off[i].variableDistance(pop_and_off[j], env);
					if (DCN > dist_mat[i][j])
						DCN = dist_mat[i][j];
				}
				objs[i][1] = DCN;
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
			new_pop.push_back(selected);
			current_members.erase(std::find(current_members.begin(), current_members.end(), selected));
		}
		auto iter = new_pop.begin();
		for (size_t i = 0; i < pop.size(); ++i)
			pop[i] = pop_and_off[*iter++];
	}
}

#endif // !OFEC_SELECT_DIVER_MULTI_H
