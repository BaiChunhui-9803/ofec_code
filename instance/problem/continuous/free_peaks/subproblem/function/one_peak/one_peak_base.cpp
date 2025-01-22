#include "one_peak_base.h"
#include "../../../type.h"
#include "../../../free_peaks.h"
#include <fstream>

namespace ofec::free_peaks {
	OnePeakBase::OnePeakBase(Problem *pro, const std::string &subspace_name, const ParameterMap &param, Random *rnd) :
		m_pro(pro), m_subspace_name(subspace_name), m_param(param),
		m_height(param.get<Real>("height")) 
	{
		if (m_pro->numberObjectives() == 1) {
			if (m_param.get<std::string>("center") == "default")
				m_center.assign(CAST_CONOP(m_pro)->numberVariables(), 0);
			else if (m_param.get<std::string>("center") == "random") {
				m_center.resize(CAST_CONOP(m_pro)->numberVariables());
				Real max = param.get<Real>("max_center");
				Real min = param.get<Real>("min_center");
				for (size_t i = 0; i < m_center.size(); ++i) {
					m_center[i] = rnd->uniform.nextNonStd<Real>(min, max);// todo
				}
			}
		}
		else if (m_pro->numberObjectives() > 1) {
			if (m_param.get<std::string>("center") == "default")
				m_center.assign(CAST_CONOP(m_pro)->numberVariables(), 0);
			else if (m_param.get<std::string>("center") == "random") {
				m_center.resize(CAST_CONOP(m_pro)->numberVariables());
				for (size_t i = 0; i < m_center.size(); ++i) {
					m_center[i] = rnd->uniform.nextNonStd<Real>(-50, 50);// todo
				}
			}
			else if (m_param.get<std::string>("center") == "designated") {
				m_center.resize(CAST_CONOP(m_pro)->numberVariables());
				for (size_t i = 0; i < m_center.size(); ++i) {
					m_center[i] = 10;// todo
				}
			}
		}
	}

	void OnePeakBase::bindData() {

	}
}