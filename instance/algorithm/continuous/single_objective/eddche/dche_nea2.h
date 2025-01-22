/********* Begin Register Information **********
{
    "description": "",
    "identifier": "DCHE_NEA2",
    "name": "DCHE-NEA2",
    "tags": [
        "continuous",
        "single-objective"
    ],
	"dependency on libraries": [ "mlpack", "Armadillo" ]
}
*********** End Register Information **********/


#ifndef OFEC_DCHE_NEA2_H
#define OFEC_DCHE_NEA2_H

#include "../../../template/framework/dche/dche.h"
#include "../multi_modal/nea2/nea2.h"

namespace ofec {
	class DCHE_NEA2 : public DCHE, public NEA2 {
		OFEC_CONCRETE_INSTANCE(DCHE_NEA2)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;

		void addSubpopsInHill(const Hill *hill, Environment *env);

		void identifyCandidates(std::list<const Solution<>*> &candidates, Environment *env) override;
	};
}

#endif // ! OFEC_DCHE_NEA2_H
