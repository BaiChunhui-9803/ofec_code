#include "map_Weierstrass_rugeness.h"
#include "../../free_peaks.h"
#include "../../../../../../utility/functional.h"

namespace ofec::free_peaks {
	WeierstrassRugenessMap::WeierstrassRugenessMap(Problem *pro, const std::string& subspace_name, const ParameterMap& param) :
		TransformBase(pro, subspace_name, param) {
		if (param.has("alpha")) {
			m_alpha = param.get<double>("alpha");
		}
		}

	void WeierstrassRugenessMap::transfer(std::vector<Real>& obj, const std::vector<Real>& var) {
		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables,0);
		auto& curOpt = m_sub_pro->optima();
		if (curOpt.numberSolutions()>0) {
			double dis = 0;
			/* TRANSFORMATION IN SEARCH SPACE*/
			for (int i = 0; i < m_number_variables; ++i) {
				tmpvect[i] = 0.;
				for (int j = 0; j < m_number_variables; ++j) {
					dis += (var[j] - curOpt.solution(0).variable().vect()[j]);
					tmpvect[i] += m_rot[i][j] * ((var[j] - curOpt.solution(0).variable().vect()[j]) * 5);
				}
			}
		}
		irregularize(tmpvect.data());
		for (int i = 0; i < m_number_variables; ++i) {
			trasX[i] = 0.;
			for (int j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_linearTF[i][j] * tmpvect[j];
			}
		}
		double tmp;
		/* COMPUTATION core*/
	//	obj[0] = 0.;

		double value(0);
		for (int i = 0; i < m_number_variables; ++i) {
			tmp = 0.;
			for (int j = 0; j < 12; ++j) {
				tmp += cos(2 * OFEC_PI * (trasX[i] + 0.5) * m_bK[j]) * m_aK[j];
			}
			value += tmp;

			//for (int idO(0); idO < m_number_objectives; ++idO) {
			//	obj[idO] += tmp;
			//}
		}
		for (int idO(0); idO < m_number_objectives; ++idO) {
			obj[idO] -= 5* pow(value / (Real)m_number_variables - m_F0, 3.);
		}
		for (int idO(0); idO < m_number_objectives; ++idO) {
			ofec::mapReal(obj[idO], m_from[idO].first, m_from[idO].second, 
				m_to[idO].first , m_to[idO].second
				);
			//obj[idO] -= 5 * pow(value / (Real)m_number_variables - m_F0, 3.);
		}

		//double curOffset= 5 * pow( - m_F0, 3.);

		//int stop = -1;
	}

	void WeierstrassRugenessMap::bindData() {
		auto& subpro = CAST_FPs(TransformBase::m_pro)->subspaceTree().name_box_subproblem.at(m_subspace_name).second;
		m_random = CAST_FPs(TransformBase::m_pro)->random();
		
		m_from = subpro->function()->objRanges();

		for (size_t j = 0; j < m_from.size(); ++j)
			m_from[j].first = 0;

		m_to = m_from;
		double curOffset = 5 * pow(-m_F0, 3.);
		for (size_t j = 0; j < m_from.size(); ++j)
			m_from[j].first -= curOffset;
		



		for (size_t j = 0; j < m_to.size(); ++j)
			m_to[j].first = 0;


		m_number_variables = CAST_FPs(TransformBase::m_pro)->numberVariables();
		m_condition_number = 100.0;
		m_F0 = 0.;
		m_aK.resize(12);
		m_bK.resize(12);



		//m_alpha = 0.5;
		loadRotation(1. / sqrt(m_condition_number));
		for (size_t i = 0; i < 12; ++i) /// number of summands, 20 in CEC2005, 10/12 saves 30% of time
		{
			m_aK[i] = pow(m_alpha, (Real)i);
			m_bK[i] = pow(3., (Real) i);
			m_F0 += m_aK[i] * cos(2 * OFEC_PI * m_bK[i] * m_alpha);
		}



	}
}