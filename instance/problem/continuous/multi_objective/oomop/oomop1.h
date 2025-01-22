/********* Begin Register Information **********
{
	"name": "MOP_OOMOP1",
	"identifier": "OOMOP1",
	"problem tags": [ "OOMOP", "MOP", "ConOP" ]
}
*********** End Register Information **********/

/*********************************************************************
* This MOP is constructed by object-oriented construction method
**********************************************************************/

// created by tanqingshan on June 20th, 2022

#ifndef OFEC_OOMOP1_H
#define OFEC_OOMOP1_H

#include"oomop.h"

namespace ofec {
	class OOMOP1:public OOMOP {
	private:
		std::vector<std::vector<std::vector<Real>>> m_global_anchor_points;
	public:
		void initialize_();
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
		Real evaluate_opt_pub_sol(std::vector<std::vector<Real>>& sols,size_t num_global_ps);

		void calPubValue(std::vector<Real>& x, std::vector<Real>& pub_objs) override;
		std::vector<std::vector<Real>> samplePubOptima(size_t dim_num, size_t num_g_ps) override;

		void generatePF() override;

	};
}

#endif // !OFEC_OOMOP1_H
