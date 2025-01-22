/********* Begin Register Information **********
{
	"name": "MOP_WFG9",
	"identifier": "WFG9",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_WFG9_H
#define OFEC_WFG9_H

#include "wfg.h"

namespace ofec {
	class WFG9 final : public WFG {
	protected:
		void t1(std::vector<Real> &y) override;
		void t2(std::vector<Real> &y) override;
		void t3(std::vector<Real> &y) override;
		void shape(std::vector<Real> &y) override;
	};
}


#endif
