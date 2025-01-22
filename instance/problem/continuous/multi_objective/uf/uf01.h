/********* Begin Register Information **********
{
	"name": "MOP_UF01",
	"identifier": "UF01",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF1_H
#define OFEC_UF1_H

#include "uf.h"

namespace ofec {
	class UF01 : public UF {
	protected:
		void evaluateObjective(Real *x, std::vector<Real> &obj) override;
	};
}

#endif
