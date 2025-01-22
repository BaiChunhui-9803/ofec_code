#include "f5_sr_schwefel.h"

namespace ofec {
	namespace cec2015 {
		void F5_SR_schwefel::initialize_(Environment *env) {
			Schwefel::initialize_();
			setDomain(-100., 100.);
			setConditionNumber(1.);
			setBias(500.);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
			setScale(1. / 10.);
			setGlobalOpt(m_translation.data());
		}
	}
}