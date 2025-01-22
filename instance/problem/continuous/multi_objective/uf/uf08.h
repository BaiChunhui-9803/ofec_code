/********* Begin Register Information **********
{
	"name": "MOP_UF08",
	"identifier": "UF08",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF8_H
#define OFEC_UF8_H

#include "uf.h"

namespace ofec {
	class UF08 : public UF {
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif
