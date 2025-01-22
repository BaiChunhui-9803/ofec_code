#include "f3_sr_ackley.h"

namespace ofec {
	namespace cec2015 {
		void F3_SR_ackley::initialize_(Environment *env) {
			Ackley::initialize_();
			setDomain(-100., 100.);
			setConditionNumber(1.);
			setBias(300.);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");
			setGlobalOpt(m_translation.data());
		}
	}
}