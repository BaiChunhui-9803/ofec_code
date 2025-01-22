#include "hghec_dynamic_adaptor.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/heap/heap.h"
#include "../../../../../../utility/functional.h"
#include "../../../../../../utility/clustering/nbc.h"
#include "../../../../../../utility/nondominated_sorting/fast_sort.h"
#include <numeric>
#include <algorithm>

namespace ofec {
	AdaptorDynamicHGHEC::AdaptorDynamicHGHEC(Problem *pro, Random *rnd, Algorithm *alg, size_t init_num) :
		AdaptorHGHEC(pro, rnd, alg, init_num) {}

	void AdaptorDynamicHGHEC::updateHills() {
		if (!m_his_sols.empty()) {
			// std::cout << hisSols().size() << std::endl;
			// if (hisSols().size() < 500) return;
			for (size_t id_hill = 0; id_hill < m_num_hills; ++id_hill) {
				std::vector<const SolutionBase*> his_sols;
				std::vector<const SolutionBase*> his_explore_sols;
				getHisSolsInHill(his_sols, id_hill);
				getHisExploreSolsInHill(his_explore_sols, id_hill);
				std::vector<const SolutionBase*> seed_sols;
				identifySeedSolutions(his_sols, seed_sols);
				subdivideSpace(seed_sols, id_hill);
			}
			groupSubspace();
		}
	}

	void AdaptorDynamicHGHEC::archiveSolution(const SolutionBase& sol, TaskHGHE task){
		AdaptorHGHE::archiveSolution(sol, task);
		updateInfoSSP(m_his_sols.back().get(), task);
	}

	void AdaptorDynamicHGHEC::clearPreHis(){
		for (auto & info_ssp : m_info_ssp){
			info_ssp.best_sol = nullptr;
			info_ssp.his_sols.clear();
			info_ssp.his_exploit.clear();
			info_ssp.his_explore.clear();
		}
		m_his_sols.clear();
		m_his_sols_explore.clear();
		m_his_sols_exploit.clear();
	}

	void AdaptorDynamicHGHEC::calculateExplorationPotential(){
		m_exploration_potential.resize(m_num_hills);
		std::vector<size_t> num_sols_each_hill;
		std::vector<Real> sols_density_each_hill;
		std::vector<Real> volume_each_hill;
		num_sols_each_hill.resize(m_num_hills);
		sols_density_each_hill.resize(m_num_hills);
		volume_each_hill.resize(m_num_hills);
		for (size_t id_hill = 0; id_hill < m_num_hills; id_hill++) {
			size_t sols_count = 0;
			Real ssp_volume = 0.;
			for (const auto & id_ssp : m_ssp_hills[id_hill]) {
				sols_count += m_info_ssp[id_ssp].his_sols.size();
			}
			for (const auto & volume : m_volume_ssp_hills[id_hill]) {
				ssp_volume += volume;
			}
			if (sols_count == 0) {
				sols_count = 1;
			}
			num_sols_each_hill[id_hill] = sols_count;
			volume_each_hill[id_hill] = ssp_volume;
			// 单位个体所占的体积
			sols_density_each_hill[id_hill] = (ssp_volume * 1.0) / sols_count;
		}
		// density decided the potential hill has
		for (int i = 0; i < sols_density_each_hill.size(); ++i){
			m_exploration_potential[i] = sols_density_each_hill[i] / (*std::max_element(sols_density_each_hill.begin(), sols_density_each_hill.end()) * 1.0);
		}
		//for (int i = 0; i < sols_density_each_hill.size(); ++i) {
		//	m_exploration_potential[i] = 1.0;
		//}
	}

	std::vector<const SolutionBase*> AdaptorDynamicHGHEC::getBestSolsEachHill(){
		std::vector<const SolutionBase*> seed_sols;
		for (const auto& hill : m_ssp_hills) {
			const SolutionBase* best_sol = nullptr;
			for (size_t id_ssp : hill) {
				if (m_info_ssp[id_ssp].best_sol && (!best_sol || m_info_ssp[id_ssp].best_sol->dominate(*best_sol, m_problem.get()))) {
					best_sol = m_info_ssp[id_ssp].best_sol;
				}
			}
			if (!best_sol)
				throw MyExcept("best solution of hill can not be null");
			seed_sols.push_back(best_sol);
		}
		return seed_sols;
	}

	void AdaptorDynamicHGHEC::getHisSolsInHill(std::vector<const SolutionBase*>& his_sols, int id_hill) {
		for (size_t id_ssp : m_ssp_hills[id_hill]) {
			auto& his_sols_ssp = m_info_ssp[id_ssp].his_sols;
			his_sols.insert(his_sols.end(), his_sols_ssp.begin(), his_sols_ssp.end());
		}
	}

	void AdaptorDynamicHGHEC::getHisExploreSolsInHill(std::vector<const SolutionBase*>& his_explore_sols, int id_hill) {
		for (size_t id_ssp : m_ssp_hills[id_hill]) {
			auto& his_sols_ssp = m_info_ssp[id_ssp].his_explore;
			his_explore_sols.insert(his_explore_sols.end(), his_sols_ssp.begin(), his_sols_ssp.end());
		}
	}

	void AdaptorDynamicHGHEC::assignSubspaceFitness() {
		size_t num_missing = 0;
		for (size_t id_ssp = 0; id_ssp < m_info_ssp.size(); ++id_ssp) {
			if (m_info_ssp[id_ssp].best_sol != nullptr)
				m_info_ssp[id_ssp].obj = m_info_ssp[id_ssp].best_sol->objective(0);
			else
				num_missing++;
		}
		if (num_missing > 0)
			bilinearInterpolation();
	}
}

