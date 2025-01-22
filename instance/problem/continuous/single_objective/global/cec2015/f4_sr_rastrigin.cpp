#include "f4_sr_rastrigin.h"

namespace ofec {
	namespace cec2015 {
		void F4_SR_rastrigin::initialize_(Environment *env) {
			Rastrigin::initialize_();
			setDomain(-100., 100.);
			setConditionNumber(1.);
			setBias(400.);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
			setScale(100. / 5.12);
			setGlobalOpt(m_translation.data());
		}
	}
}