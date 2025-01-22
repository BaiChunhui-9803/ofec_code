#include"ttFun.h"
#include"ttFunPar.h"

const int ofec::sp::FunTt::m_numFun = 3;
const double ofec::sp::FunsTt::Three_Sigma_Limit = 3.0;
double ofec::sp::FunsTt::getRandomWithin3Sigma(double rand_num) {
	if (rand_num > Three_Sigma_Limit) rand_num = Three_Sigma_Limit;
	else if (rand_num < -Three_Sigma_Limit) rand_num = -Three_Sigma_Limit;
	return rand_num;
}

void ofec::sp::FunsTt::initialize(Random *rnd, const FunsTtPar& par) {
	for (int id_fun(0); id_fun < 2; ++id_fun) {
		if (par.m_fun_exist[id_fun]) {
			m_funs[id_fun].reset(new FunTt());
			for (int idx(0); idx < m_funs[id_fun]->m_funs.size(); ++idx) {
				par.m_pars[idx].resetFun(m_funs[id_fun]->m_funs[idx]);
				m_funs[id_fun]->m_funs[idx]->initialize(rnd->uniform);
			}
		}
		else {
			m_funs[id_fun].reset(nullptr);
		}
		m_funs[id_fun]->m_type = par.m_fun_type;
	}
	m_noisy_ratio = par.m_obj_noisy_ratio;
	m_random_type = par.m_random_type;
	m_range_type = par.m_range_type;
}
void ofec::sp::FeasibleDistributionFunTt::initialize(Random *rnd, const FunsTtPar& par) {
	DistributionFunTt::initialize(rnd, par);
	m_feasible_threshold = par.m_feasible_threshold;
}

void ofec::sp::Distribution3DFunTt::initialize(Random *rnd, const FunsTtPar& par) {
	FunsTt::initialize(rnd, par);
	m_scale_xy = par.m_scale_xy;
	m_offset_xy = par.m_offset_xy;
}

