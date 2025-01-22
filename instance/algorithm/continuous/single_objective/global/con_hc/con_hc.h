/********* Begin Register Information **********
{
	"name": "ConHC",
	"identifier": "ConHC",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#include "../../../../template/classic/hill_climbing/hc.h"
#include "../../../../../../core/problem/encoding.h"


#ifndef OFEC_CON_HC_H
#define OFEC_CON_HC_H

namespace ofec {
	class ConHC : virtual public HC<VariableVector<Real>> {
		OFEC_CONCRETE_INSTANCE(ConHC)
	protected:
		std::vector<Real> m_step_size;

		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void pickNeighbour(Environment *env) override;
	};
}

#endif