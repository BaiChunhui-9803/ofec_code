/********* Begin Register Information **********
{
    "description": "",
    "identifier": "DCHE_HillVallEA",
    "name": "DCHE-HillVallEA",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_DCHE_HILLVALLEA_H
#define OFEC_DCHE_HILLVALLEA_H

#include "../../../template/framework/dche/dche.h"
#include "../multi_modal/hill_vallea/hill_vallea.h"

namespace ofec {
	class DCHE_HillVallEA : public DCHE, public HillVallEA {
		OFEC_CONCRETE_INSTANCE(DCHE_HillVallEA)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		std::vector<std::shared_ptr<Solution<>>> uniformlySampleInHill(const Hill *hill, Environment *env);

		void identifyCandidates(std::list<const Solution<> *> &candidates, Environment *env) override;
	};
}

#endif // ! OFEC_DCHE_HILLVALLEA_H
