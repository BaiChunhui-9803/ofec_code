#include "dynamic_noisy_strategy.h"

bool ofec::DynamicNoisyStartegy::judgeChange(Problem *pro, const std::vector<Real>& originObjs, const std::vector<Real>& curObjs) const {
	if (pro->hasTag(ProblemTag::kNoisyOP)) {
		if (curObjs.size() != originObjs.size())return true;
		for (int objId(0); objId < originObjs.size(); ++objId) {
			if (abs(originObjs[objId] - curObjs[objId]) / abs(originObjs[objId]) > m_uncertianty_thresHold) {
				return true;
			}
		}
		return false;
	}
	else return originObjs != curObjs;
}
