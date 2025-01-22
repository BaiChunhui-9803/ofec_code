/********* Begin Register Information **********
{
	"name": "NBDE",
	"identifier": "NBDE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_NBDE_H
#define OFEC_NBDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class NBDE : public Algorithm {
		OFEC_CONCRETE_INSTANCE(NBDE)
	protected:
		PopulationDE<> m_pop;
		size_t m_pop_size;

		void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // !OFEC_NBDE_H
