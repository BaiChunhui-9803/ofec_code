#include "nl_shade_lbc.h"
#include "optimizer.h"
#include "../../../../../../core/problem/continuous/continuous.h"


namespace ofec {


	void NL_SHADE_LBC::initialize_(Environment* env) {
		Algorithm::initialize_(env);
	}
	void NL_SHADE_LBC::run_(Environment* env) {
		
		ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
		int GNVars = con_pro->numberVariables();


		nl_shade_lbc::Optimizer OptZ;
		////////////////NInds     NVars  func     Run  memory    arch size
		OptZ.Initialize(GNVars * 23, GNVars, 20 * GNVars, 1, m_random.get());
		OptZ.MainCycle(env, m_random.get());
		OptZ.Clean();
	}
}