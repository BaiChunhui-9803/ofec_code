/********* Begin Register Information **********
{
	"name": "FEP",
	"identifier": "FEP",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_FEP_H
#define OFEC_FEP_H

#include "../../../../template/classic/evolutionary_programming/population.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class PopFEP : public PopEP<> {
	protected:
		void mutate(Random *rnd, Environment *env) override;
	};

	class FEP : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(FEP)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

	protected:
		Real m_tau, m_tau_prime;
		size_t m_q, m_pop_size;
		PopFEP m_pop;
	};
}

#endif // !OFEC_FEP_H
