/********* Begin Register Information **********
{
	"name": "MOP_UF07",
	"identifier": "UF07",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF7_H
#define OFEC_UF7_H

#include "uf.h"

namespace ofec {
	class UF07 : public UF {
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif
