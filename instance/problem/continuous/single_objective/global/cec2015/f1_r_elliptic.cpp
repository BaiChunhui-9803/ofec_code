#include "f1_r_elliptic.h"

namespace ofec {
	namespace cec2015 {
		void F1_R_elliptic::initialize_(Environment *env) {
			Elliptic::initialize_();
			setConditionNumber(1.);
			setBias(100.);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
			setGlobalOpt(m_translation.data());
		}
	}
}