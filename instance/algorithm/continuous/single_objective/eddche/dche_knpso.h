/********* Begin Register Information **********
{
	"description": "",
	"identifier": "DCHE_kNPSO",
	"name": "DCHE-kNPSO",
	"tags": [
		"continuous",
		"single-objective"
	],
	"dependency on libraries": [ "LibTorch" ]
}
*********** End Register Information **********/


#ifndef OFEC_DCHE_KNPSO_H
#define OFEC_DCHE_KNPSO_H

#include "../../../template/framework/dche/dche.h"
#include "../eeec/knpso.h"
#include "../../../template/framework/dche/tracer_pop_nbd.h"

namespace ofec {
	class DCHE_kNPSO : public DCHE, public kNPSO, public TracerPopNBD {
		OFEC_CONCRETE_INSTANCE(DCHE_kNPSO)
	protected:
		size_t m_max_stagant_iters;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		void initializePopulationInHill(const Hill *hill, Environment *env);

		void identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) override;
	};
}

#endif // ! OFEC_DCHE_KNPSO_H
