/********* Begin Register Information **********
{
	"description": "",
	"identifier": "ESPDE",
	"name": "ESPDE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/********************************************************************
@article{li2023history,
  volume = {631},
  journal = {Information Sciences},
  pages = {15--30},
  year = {2023},
  author = {Yu Li and Lingling Huang and Weifeng Gao and Zhifang Wei
	and Tianqi Huang and Jingwei Xu and Maoguo Gong},
  title = {History information-based hill-valley technique
	for multimodal optimization problems}
}
*********************************************************************/

#ifndef OFEC_ESPDE_H
#define OFEC_ESPDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/differential_evolution/population.h"
#include "../../../../../../core/algorithm/multi_population.h"

namespace ofec {
	class ESPDE :virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(ESPDE)
	protected:
		enum State { kExplore, kExploit };
		
		size_t m_NP, m_t, m_K;
		std::vector<std::shared_ptr<Solution<>>> m_P;
		std::list<std::shared_ptr<const Solution<>>> m_EI;
		MultiPopulation<PopulationDE<>> m_S;
		Real m_c;
		Real m_positive;
		Real m_f_max, m_f_min;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		int initPop(Environment *env);
		int speciation(Environment *env);
		Real sigmaNegInf(const PopulationDE<> &S_i, Environment *env) const;
		Real gamma() const;
		size_t nearest(const PopulationDE<> &S_i, const Solution<> &s, Environment *env) const;
		void generateParam(Real mu_F, Real mu_CR, const PopulationDE<> &S_i,
			std::vector<Real> &F_i, std::vector<Real> &CR_i);
		void updateMuParam(const std::list<Real> & A_F, const std::list<Real> &A_CR,
			Real &mu_F, Real &mu_CR);
		void updateFMinMax(const Solution<> &s);

	public:
		static bool test(const Solution<> &x_a, const Solution<> &x_b,
			const std::list<std::shared_ptr<const Solution<>>> &EI, Environment *env);
	};
}

#endif // ! OFEC_ESPDE_H
