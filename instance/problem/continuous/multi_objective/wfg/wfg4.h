/********* Begin Register Information **********
{
	"name": "MOP_WFG4",
	"identifier": "WFG4",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_WFG4_H
#define OFEC_WFG4_H

#include "wfg.h"

namespace ofec {
	class WFG4 final : public WFG {
	protected:
		void t1(std::vector<Real> &y) override;
		void t2(std::vector<Real> &y) override;
		void shape(std::vector<Real> &y) override;
	};
}


#endif
