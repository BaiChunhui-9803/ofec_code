/********* Begin Register Information **********
{
	"name": "MOP_WFG6",
	"identifier": "WFG6",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_WFG6_H
#define OFEC_WFG6_H

#include "wfg.h"

namespace ofec {
	class WFG6 final : public WFG {
	protected:
		void t1(std::vector<Real> &y) override;
		void t2(std::vector<Real> &y) override;
		void shape(std::vector<Real> &y) override;
	};
}


#endif
