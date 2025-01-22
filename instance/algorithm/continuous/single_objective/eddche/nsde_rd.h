/********* Begin Register Information **********
{
    "description": "NSDE with restart-and-double strategy",
    "identifier": "NSDE_RD",
    "name": "NSDE-RD",
    "tags": [
        "single-objective",
        "continuous"
    ]
}
*********** End Register Information **********/


#ifndef OFEC_NSDE_RD_H
#define OFEC_NSDE_RD_H

#include "../multi_modal/nmde/nsde.h"

namespace ofec {
	class NSDE_RD : public NSDE {
		OFEC_CONCRETE_INSTANCE(NSDE_RD)
	protected:
		size_t m_max_stagant_iters;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}

#endif // ! OFEC_NSDE_RD_H
