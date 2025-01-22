/********* Begin Register Information **********
{
	"name": "MOP_RMMEDA_F9",
	"identifier": "RMMEDA_F9",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_RMMEDA_F9_H
#define OFEC_RMMEDA_F9_H

#include "rmmeda_f.h"

namespace ofec {
	class RMMEDA_F9 : public RMMEDA_F {
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