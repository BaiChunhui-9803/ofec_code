#include "f2_r_cigar.h"

namespace ofec {
	namespace cec2015 {
		void F2_R_cigar::initialize_(Environment *env) {
			BentCigar::initialize_();
			setConditionNumber(1.);
			setBias(200.);
			loadTranslation("/instance/problem/continuous/single_objective/global/cec2015/");
			loadRotation("/instance/problem/continuous/single_objective/global/cec2015/");	
			setGlobalOpt(m_translation.data());
		}
	}
}