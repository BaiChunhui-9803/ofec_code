/********* Begin Register Information **********
{
	"name": "MOP_UF06",
	"identifier": "UF06",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF6_H
#define OFEC_UF6_H

#include "uf.h"

namespace ofec {
	class UF06 : public UF {
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif

