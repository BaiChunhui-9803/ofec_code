/********* Begin Register Information **********
{
    "description": "",
    "identifier": "GuidedDE",
    "name": "Guided-DE",
    "tags": [
        "continuous",
        "single-objective"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_GUIDED_DE_H
#define OFEC_GUIDED_DE_H

#include "canonical_de.h"

namespace ofec {
	class GuidedDE : public CanonicalDE{
		OFEC_CONCRETE_INSTANCE(GuidedDE)
	protected:
        void addInputParameters();
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_GUIDED_DE_H
