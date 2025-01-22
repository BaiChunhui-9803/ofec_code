#ifndef OFEC_CLUSTERING_PCC_H
#define OFEC_CLUSTERING_PCC_H

#include <torch/script.h>
#include "../../../core/problem/encoding.h"
#include "../../../core/algorithm/population.h"
#include "../../../core/algorithm/multi_population.h"
#include "../nbc.h"

namespace ofec {
	class PCC {
	private:
		int m_problem.get();
		torch::jit::script::Module m_img_seg_mdl;
		int m_nbd_resolution;
		std::list<std::unique_ptr<const SolutionBase>> m_sols;
		std::vector<std::vector<const SolutionBase*>> m_evo_pop;
		std::vector<std::vector<Real>> m_evo_nbd;
#ifdef OFEC_DEMO
		std::vector<std::vector<bool>> m_outlier;
		std::vector<std::vector<const SolutionBase*>> m_kmeans_m_clusters;
#endif
		std::vector<std::vector<Real>> m_evo_nbd2;
		std::vector<std::vector<const SolutionBase*>> m_clusters;
		std::vector<const SolutionBase*> m_cluster_centers;

	public:
		void setIdPro(Problem *pro) { m_problem.get() = pro; }
		void setModel(const std::string &file_name);
		template<typename TInd>
		void addLatestPop(const Population<TInd> &pop) {
			m_evo_pop.resize(m_evo_pop.size() + 1);
			for (size_t i = 0; i < pop.size(); ++i) {
				m_sols.emplace_back(new TInd(pop[i]));
				m_evo_pop.back().push_back(m_sols.back().get());
			}
			NBC nbc;
			nbc.setData(m_evo_pop.back(), m_problem.get());
			nbc.updateNbDistByKDTree();
			m_evo_nbd.push_back(nbc.nearestBetterDis());
#ifdef OFEC_DEMO
			m_outlier.push_back(std::vector<bool>(m_evo_nbd.back().size(), false));
#endif
			m_evo_nbd2.push_back(nbc.nearestBetterDis2());
		}
		template<typename TPopulation>
		void addLatestPop(const MultiPopulation<TPopulation> &multi_pop) {
			m_evo_pop.resize(m_evo_pop.size() + 1);
			for (size_t k = 0; k < multi_pop.size(); ++k) {
				for (size_t i = 0; i < multi_pop[k].size(); ++i) {
					m_sols.emplace_back(new TPopulation::IndividualType(multi_pop[k][i]));
					m_evo_pop.back().push_back(m_sols.back().get());
				}
			}
			NBC nbc;
			nbc.setData(m_evo_pop.back(), m_problem.get());
			nbc.updateNbDistByKDTree();
			m_evo_nbd.push_back(nbc.nearestBetterDis());
#ifdef OFEC_DEMO
			m_outlier.push_back(std::vector<bool>(m_evo_nbd.back().size(), false));
#endif
			m_evo_nbd2.push_back(nbc.nearestBetterDis2());
		}
		void clustering();
		void clear();

		const std::vector<std::vector<const SolutionBase*>>& evoPop() const { return m_evo_pop; }
		const std::vector<std::vector<Real>>& evoNBD() const { return m_evo_nbd; }
#ifdef OFEC_DEMO
		const std::vector<std::vector<bool>>& outlier() const { return m_outlier; }
		const std::vector<std::vector<const SolutionBase*>>& kMeansClusters() const { return m_kmeans_m_clusters; }
#endif
		const std::vector<std::vector<const SolutionBase*>>& clusters() const { return m_clusters; }
		const std::vector<const SolutionBase*>& clusterCenters() const { return m_cluster_centers; }
	};
}

#endif // !OFEC_CLUSTERING_PCC_H
