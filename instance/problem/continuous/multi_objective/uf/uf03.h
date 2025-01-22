/********* Begin Register Information **********
{
	"name": "MOP_UF03",
	"identifier": "UF03",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF3_H
#define OFEC_UF3_H

#include "uf.h"

namespace ofec {
	class UF03 : public UF{
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif

