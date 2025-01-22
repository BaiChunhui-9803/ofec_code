/********* Begin Register Information **********
{
	"name": "MOP_UF02",
	"identifier": "UF02",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF2_H
#define OFEC_UF2_H

#include "uf.h"

namespace ofec {
	class UF02 : public UF{
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif
