/********* Begin Register Information **********
{
	"name": "ConSA",
	"identifier": "ConSA",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#include "../../../../template/classic/simulated_annealing/sa.h"
#include "../../../../../../core/problem/encoding.h"


#ifndef OFEC_CON_SA_H
#define OFEC_CON_SA_H
namespace ofec {
	class ConSA : virtual public SA<VariableVector<Real>> {
		OFEC_CONCRETE_INSTANCE(ConSA)
	protected:
		std::vector<Real> m_step_size;

		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void pickNeighbour(Environment *env) override;
	};
}

#endif