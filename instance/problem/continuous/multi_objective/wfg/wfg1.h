/********* Begin Register Information **********
{
	"name": "MOP_WFG1",
	"identifier": "WFG1",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_WFG1_H
#define OFEC_WFG1_H

#include "wfg.h"

namespace ofec {
	class WFG1 final : public WFG {
	protected:
		void t1(std::vector<Real> &y) override;
		void t2(std::vector<Real> &y) override;
		void t3(std::vector<Real> &y) override;
		void t4(std::vector<Real> &y) override;
		void shape(std::vector<Real> &y) override;
	};
}


#endif