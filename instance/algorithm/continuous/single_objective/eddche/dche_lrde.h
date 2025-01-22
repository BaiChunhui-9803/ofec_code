/********* Begin Register Information **********
{
	"description": "",
	"identifier": "DCHE_LRDE",
	"name": "DCHE-LRDE",
	"tags": [
		"continuous",
		"single-objective"
	]
}
*********** End Register Information **********/


#ifndef OFEC_DCHE_LRDE_H
#define OFEC_DCHE_LRDE_H

#include "../../../template/framework/dche/dche.h"
#include "../eeec/lrde.h"

namespace ofec {
	class DCHE_LRDE : public DCHE, public LRDE {
		OFEC_CONCRETE_INSTANCE(DCHE_LRDE)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		bool reject(const Solution<> &sol, Environment *env) const override;
		void initializePopulationInHill(const Hill *hill, Environment *env);

		void identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) override;
	};
}

#endif // ! OFEC_DCHE_LRDE_H
