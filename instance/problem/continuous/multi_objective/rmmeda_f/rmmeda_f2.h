/********* Begin Register Information **********
{
	"name": "MOP_RMMEDA_F2",
	"identifier": "RMMEDA_F2",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_RMMEDA_F2_H
#define OFEC_RMMEDA_F2_H

#include "rmmeda_f.h"

namespace ofec {
	class RMMEDA_F2 : public RMMEDA_F {
	public:
		void updateOptima() override;
		void sampleParetoSols(size_t sample_num);

		/*void sampleParetoFront(size_t sample_num);
		void loadParetoFront(size_t sample_num);*/
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif