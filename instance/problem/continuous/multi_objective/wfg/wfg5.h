/********* Begin Register Information **********
{
	"name": "MOP_WFG5",
	"identifier": "WFG5",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_WFG5_H
#define OFEC_WFG5_H

#include "wfg.h"

namespace ofec {
	class WFG5 final : public WFG {
	protected:
		void t1(std::vector<Real> &y) override;
		void t2(std::vector<Real> &y) override;
		void shape(std::vector<Real> &y) override;
	};
}


#endif
