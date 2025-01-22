/********* Begin Register Information **********
{
	"name": "MOP_UF09",
	"identifier": "UF09",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UF9_H
#define OFEC_UF9_H

#include "uf.h"

namespace ofec {
	class UF09 : public UF {
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif
