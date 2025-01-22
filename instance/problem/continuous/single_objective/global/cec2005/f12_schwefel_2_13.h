/********* Begin Register Information **********
{
	"name": "GOP_CEC2005_F12",
	"identifier": "GOP_CEC2005_F12",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_F12_SCHWEFEL_2_13_H
#define OFEC_F12_SCHWEFEL_2_13_H

#include "../../function.h"
#include "../../../../single_objective/metrics_gop.h"

namespace ofec {
	namespace cec2005 {
		class Schwefel_2_13 : public Function, public MetricsGOP {
			OFEC_CONCRETE_INSTANCE(Schwefel_2_13)
		protected:
			void addInputParameters();
			void initialize_(Environment *env) override;
			void loadData(const std::string & path);
			void evaluateOriginalObj(Real *x, std::vector<Real>& obj) const override;
		private:
			std::vector<std::vector<Real>> m_a;
			std::vector<std::vector<Real>> m_b;
			std::vector<Real> m_alpha;
		};
	}
	using GOP_CEC2005_F12 = cec2005::Schwefel_2_13;
}
#endif // !OFEC_F12_SCHWEFEL_2_13_H
