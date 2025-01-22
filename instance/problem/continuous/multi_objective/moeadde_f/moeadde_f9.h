/********* Begin Register Information **********
{
	"name": "MOP_MOEADDE_F9",
	"identifier": "MOEADDE_F9",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_MOEADDE_F9_H
#define OFEC_MOEADDE_F9_H

#include "moeadde_f.h"

namespace ofec {
	class MOEADDE_F9 : public MOEADDE_F {
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