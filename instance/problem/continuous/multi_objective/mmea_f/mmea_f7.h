/********* Begin Register Information **********
{
	"name": "MOP_MMEA_F7",
	"identifier": "MMEA_F7",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_MMEA_F7_H
#define OFEC_MMEA_F7_H

#include "mmea_f.h"

namespace ofec {
	class MMEA_F7 : public MMEA_F {
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