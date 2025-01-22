#ifndef OFEC_PopGL_NBN_COM_ALG_H
#define OFEC_PopGL_NBN_COM_ALG_H

#include "../../../../core/algorithm/population.h"
#include "../../../../core/problem/solution.h"


#include "../../../problem/combination/travelling_salesman/travelling_salesman.h"

namespace ofec {
	class PopGL_NBN_COM_ALG : public Population<Solution<VarVec<int>>> {
	public:
		PopGL_NBN_COM_ALG() = default;
		PopGL_NBN_COM_ALG(size_t size_pop, Problem* pro) : Population(size_pop, pro) {};
		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override;
		virtual void initialize(Problem* pro, Random* rnd) override; 
		int evaluate(Problem* pro, Algorithm* alg) override;
		virtual ~PopGL_NBN_COM_ALG() = default;

		void udpateProMat(Problem* pro);

		void getSolInfo(std::vector<int>& belong, std::vector<double>& dis2parent, 
			std::vector<double>& fitness) const{
			belong = mv_belong;
			dis2parent = mv_dis2parent;
			fitness = m_fitness;
		}

		const std::vector<std::vector<double>>& getProMat() const{
			return mvv_pro_mat;
		}

	protected:
		std::vector<int> mv_belong;
		std::vector<double> mv_dis2parent;
		std::vector<std::vector<ofec::Real>> mvv_pro_mat;
		std::vector<double> m_fitness;
		std::vector<double> m_weight;
		std::vector<double> m_obj;
		std::vector<std::vector<std::array<int, 2>>>  mv_sol_edges;

		std::shared_ptr<ofec::Random> m_rnd;
	};
}

#endif // !OFEC_PopGL_NBN_COM_ALG_H