/********* Begin Register Information **********
{
    "description": "",
    "identifier": "GuidedJDE",
    "name": "guided-jDE",
    "tags": [
        "single-objective",
        "continuous"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_GUIDED_JDE_H
#define OFEC_GUIDED_JDE_H

#include "../global/jde/jde.h"

namespace ofec {
	class GuidedJDE : public jDE {
		OFEC_CONCRETE_INSTANCE(GuidedJDE)
	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
	};

	class PopGuidedJDE final : public PopJDE {
	public:
		PopGuidedJDE(size_t size_pop, Environment* env);
		void initialize(Environment *env, Random *rnd) override;
	};
}

#endif // ! OFEC_GUIDED_JDE_H
