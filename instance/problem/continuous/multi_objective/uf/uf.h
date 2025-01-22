/****************************************************************************
* ------------------------------ - Reference--------------------------------
* Q.Zhang, A.Zhou, S.Zhao, P.N.Suganthan, W.Liu, and S.Tiwari,
* Multiobjective optimization test instances for the CEC 2009 special
* session and competition, School of CS& EE, University of Essex, Working
* Report CES - 487, 2009.
****************************************************************************/


#ifndef OFEC_UF_H
#define OFEC_UF_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../multi_objective/metrics_mop.h"

namespace ofec {
	class UF : public Continuous, public MetricsMOP  {
	protected:
		void initialize_() override;
		void updateOptima() override;
		void loadParetoFront(size_t sample_num);
		void sampleParetoSets(size_t sample_num);
		void sampleParetoFront(size_t sample_num);
		//void loadParetoFront();
	};
}

#endif
