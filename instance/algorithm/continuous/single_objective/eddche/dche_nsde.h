/********* Begin Register Information **********
{
	"description": "",
	"identifier": "DCHE_NSDE",
	"name": "DCHE-NSDE",
	"tags": [
		"continuous",
		"single-objective"
	]
}
*********** End Register Information **********/


#ifndef OFEC_DCHE_NSDE_H
#define OFEC_DCHE_NSDE_H

#include "../../../template/framework/dche/dche.h"
#include "../multi_modal/nmde/nsde.h"

namespace ofec {
	class DCHE_NSDE : public DCHE, public NSDE {
		OFEC_CONCRETE_INSTANCE(DCHE_NSDE)
	protected:
		size_t m_max_stagant_iters;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		void initPopInHill(const Hill *hill, Environment *env);

		void identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) override;
	};
}

#endif // ! OFEC_DCHE_NSDE_H
