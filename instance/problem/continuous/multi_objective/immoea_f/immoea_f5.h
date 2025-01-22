/********* Begin Register Information **********
{
	"name": "MOP_IMMOEA_F5",
	"identifier": "IMMOEA_F5",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_IMMOEA_F5_H
#define OFEC_IMMOEA_F5_H

#include "immoea_f.h"

namespace ofec {
	class IMMOEA_F5 : public IMMOEA_F {
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