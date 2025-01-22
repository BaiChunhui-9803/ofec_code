/********* Begin Register Information **********
{
	"name": "MOP_MOEADDE_F8",
	"identifier": "MOEADDE_F8",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_MOEADDE_F8_H
#define OFEC_MOEADDE_F8_H

#include "moeadde_f.h"

namespace ofec {
	class MOEADDE_F8 : public MOEADDE_F {
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