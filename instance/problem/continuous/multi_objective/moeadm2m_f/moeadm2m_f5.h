/********* Begin Register Information **********
{
	"name": "MOP_MOEADM2M_F5",
	"identifier": "MOEADM2M_F5",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_MOEADM2M_F5_H
#define OFEC_MOEADM2M_F5_H

#include "moeadm2m_f.h"

namespace ofec {
	class MOEADM2M_F5 : public MOEADM2M_F {
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