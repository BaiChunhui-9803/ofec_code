#ifndef OFEC_TRACER_POP_NBD_H
#define OFEC_TRACER_POP_NBD_H

#include "../../../../../core/algorithm/algorithm.h"
#include <torch/script.h>
#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/algorithm/multi_population.h"
#include "../../../../../utility/clustering/nbc.h"
#include "../../../../../core/problem/continuous/continuous.h"

namespace ofec {
	class TracerPopNBD : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(TracerPopNBD)
	private:
		std::string m_model_file_name;
		int m_max_num_iters;
		torch::jit::script::Module m_img_seg_mdl;
		int m_nbd_resolution;
		std::vector<std::vector<const SolutionBase*>> m_evo_pop;
		std::vector<std::vector<Real>> m_evo_nbd;
		std::vector<std::vector<bool>> m_outlier;
		std::vector<std::vector<Real>> m_evo_nbd2;
		std::list<std::unique_ptr<const SolutionBase>> m_sols;

	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;

		template<typename TIndividual>
		void addLatestPop(const Population<TIndividual> &pop, Environment *env) {
			if (m_evo_pop.size() < m_max_num_iters) {
				m_evo_pop.resize(m_evo_pop.size() + 1);
				for (size_t i = 0; i < pop.size(); ++i) {
					m_sols.emplace_back(new typename TIndividual::SolutionType(pop[i]));
					m_evo_pop.back().push_back(m_sols.back().get());
				}
				NBC nbc;
				nbc.setData(m_evo_pop.back());
				if (CAST_CONOP(env->problem()) != nullptr) {
					nbc.updateNbDistByKDTree(env);
				}
				else {
					nbc.updateNbDistByDistMat(env);
				}
				m_evo_nbd.push_back(nbc.nearestBetterDis());
				m_outlier.push_back(std::vector<bool>(m_evo_nbd.back().size(), false));
				m_evo_nbd2.push_back(nbc.nearestBetterDis2());
			}
		}

		template<typename TPopulation>
		void addLatestPop(const MultiPopulation<TPopulation> &multi_pop, Environment *env) {
			if (m_evo_pop.size() < m_max_num_iters) {
				m_evo_pop.resize(m_evo_pop.size() + 1);
				for (size_t k = 0; k < multi_pop.size(); ++k) {
					for (size_t i = 0; i < multi_pop[k].size(); ++i) {
						m_sols.emplace_back(new typename TPopulation::IndividualType::SolutionType(multi_pop[k][i]));
						m_evo_pop.back().push_back(m_sols.back().get());
					}
				}
				NBC nbc;
				nbc.setData(m_evo_pop.back());
				if (CAST_CONOP(env->problem()) != nullptr) {
					nbc.updateNbDistByKDTree(env);
				}
				else {
					nbc.updateNbDistByDistMat(env);
				}
				m_evo_nbd.push_back(nbc.nearestBetterDis());
				m_outlier.push_back(std::vector<bool>(m_evo_nbd.back().size(), false));
				m_evo_nbd2.push_back(nbc.nearestBetterDis2());
			}
		}

		void clearEvoPop();

		void identifyOutliers(std::list<const Solution<>*> &solutions, Environment *env);

	public:
		const std::vector<std::vector<Real>>& evoNBD() const { return m_evo_nbd; }
		const std::vector<std::vector<bool>>& outlier() const { return m_outlier; }
	};
}

#endif // !OFEC_TRACE_POP_NBD_H
