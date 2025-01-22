/********* Begin Register Information **********
{
	"name": "MOP_MMEA_F4",
	"identifier": "MMEA_F4",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_MMEA_F4_H
#define OFEC_MMEA_F4_H

#include "mmea_f.h"

namespace ofec {
	class MMEA_F4 : public MMEA_F {
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