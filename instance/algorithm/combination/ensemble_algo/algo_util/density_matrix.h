
#ifndef DENSITY_MATRIX_H
#define DENSITY_MATRIX_H

//#include "../../../../../core/definition.h"
#include"../../../../../core/definition.h"
#include "../../../../../core/problem/encoding.h"
#include "../../../template/combination/sequence/sequence_algorithm.h"
#include"fitness_calculator.h"
#include "../../../../../utility/function/custom_function.h"
#include<vector>
#include<list>

namespace ofec {
	class DensityMatrix {
	protected:
		std::vector<std::vector<double>> m_density_matrix;
		std::vector<std::vector<double>> m_temp_matrix;
		double m_sum_fitness = 0;
		double m_num_var = 0;
		
		enum class DisType { Eula, PopDensity, IndiPop };
		DisType m_dis_type = DisType::IndiPop;
	public:

		void initialize(Problem *pro, Algorithm *alg) {
			//FitnessCalculator::initialize(pro);
			m_num_var = pro->numberVariables();
			auto& matSize(GET_ASeq(alg).interpreter().getMatrixSize());
			UTILITY::assignVVector<double>(m_density_matrix, matSize, 0);
			UTILITY::assignVVector<double>(m_temp_matrix, matSize, 0);
			//m_density_matrix.resize(matSize.size());
			//for (int idx(0); idx < m_density_matrix.size(); ++idx) {
			//	m_density_matrix[idx].resize(matSize[idx]);
			//	std::fill(m_density_matrix[idx].begin(), m_density_matrix.end(), 0);
			//}
		}
		void clear() {
			UTILITY::assignVVector<double>(m_density_matrix, 0);
			m_sum_fitness = 0;
		}

		template<typename TIndi>
		void update(
			const std::vector<std::unique_ptr<TIndi>>& pop,
			const std::vector<double>& fitness) {
			//FitnessCalculator::update(fitness);
			auto& popFit(fitness);
			for (auto& it : popFit) {
				m_sum_fitness += it;
			}
			for (int idx(0); idx < pop.size(); ++idx) {
				for (auto& edge : pop[idx]->edges()) {
					m_density_matrix[edge.first][edge.second] += popFit[idx];
				}
			}
		}


		template<typename TIndi>
		double disBetweenInds(
			const TIndi& sol1, double fitness1,
			const TIndi& sol2, double fitness2)const {
			if (m_dis_type == DisType::Eula) {

			}
			double sumFit(fitness1 + fitness2 + m_sum_fitness);
			UTILITY::assignVVector<double>(m_temp_matrix, 0);
			for (auto& edge : sol1.edges()) {
				m_temp_matrix[edge.first][edge.second] +=
					m_density_matrix[edge.first][edge.second] + fitness1;
			}
			double totalFit(0);
			for (auto& edge : sol2.edges()) {
				totalFit+=m_temp_matrix[edge.first][edge.second];
			}
			return 1.0 - totalFit / (sumFit*m_num_var);
		}



		template<typename TIndi>
		double disToPop(const TIndi & sol, double fitness) const{
			double popFit(fitness);
			double sumFit(m_sum_fitness + popFit);
			double dis(0);
			for (auto& edge : sol.edges()) {
				dis += m_density_matrix[edge.first][edge.second];
			}
			dis += popFit * m_num_var;
			dis = dis / (sumFit * m_num_var);
			return 1.0- dis;
		}


		template<typename TIndi>
		double innerDisToPop(const TIndi& sol, double fitness)const {
			double dis(0);
			double sumFit(m_sum_fitness);
			for (auto& edge : sol.edges()) {
				dis += m_density_matrix[edge.first][edge.second];
			}
			dis = dis / (sumFit * m_num_var);
			return 1.0 - dis;
		}
	};
}


#endif