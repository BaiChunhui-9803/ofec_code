#include "moea_f.h"
#include "../../../../../core/global.h"
#include <fstream>
#include <sstream>

namespace ofec {
	void MOEA_F::initialize_() {
		Continuous::initialize_();
		auto& v = *m_param;;

		resizeObjective(v.get<int>("number of objectives"));
		if (m_name == "MOEA_F6") {
			if (m_number_objectives != 3)
				throw MyExcept("The number of objectives must be equal to 3.");
		}
		else {
			if (m_number_objectives != 2)
				throw MyExcept("The number of objectives must be equal to 2.");
		}
		for (size_t j = 0; j < m_number_objectives; j++)
			setOptMode(OptimizeMode::kMinimize, j);

		resizeVariable(v.get<int>("number of variables")); //recomend n=10
		if (m_number_variables < m_number_objectives)
			throw MyExcept("The number of variables must be no less than the number of objectives.");
		if (m_name == "MOEA_F1") {
			setDomain(0., 1.);
			m_dtype = 1;
			m_ptype = 21;
			m_ltype = 21;
		}
		else if (m_name == "MOEA_F2") {
			m_domain.setRange(0., 1., 0);
			for (size_t i = 1; i < m_number_variables; ++i)
				m_domain.setRange(-1., 1., i);
			m_dtype = 1;
			m_ptype = 21;
			m_ltype = 22;
		}
		else if (m_name == "MOEA_F3") {
			m_domain.setRange(0., 1., 0);
			for (size_t i = 1; i < m_number_variables; ++i)
				m_domain.setRange(-1., 1., i);
			m_dtype = 1;
			m_ptype = 21;
			m_ltype = 23;
		}
		else if (m_name == "MOEA_F4") {
			m_domain.setRange(0., 1., 0);
			for (size_t i = 1; i < m_number_variables; ++i)
				m_domain.setRange(-1., 1., i);
			m_dtype = 1;
			m_ptype = 21;
			m_ltype = 24;
		}
		else if (m_name == "MOEA_F5") {
			m_domain.setRange(0., 1., 0);
			for (size_t i = 1; i < m_number_variables; ++i)
				m_domain.setRange(-1., 1., i);
			m_dtype = 1;
			m_ptype = 21;
			m_ltype = 26;
		}
		else if (m_name == "MOEA_F6") {
			for (size_t i = 0; i < m_number_variables; ++i) {
				if (i == 0 || i == 1)
					m_domain.setRange(0., 1., i);
				else
					m_domain.setRange(-2., 2., i);
			}
			m_dtype = 1;
			m_ptype = 31;
			m_ltype = 32;
		}
		else if (m_name == "MOEA_F7") {
			setDomain(0., 1.);
			m_dtype = 3;
			m_ptype = 21;
			m_ltype = 21;
		}
		else if (m_name == "MOEA_F8") {
			setDomain(0., 1.);
			m_dtype = 4;
			m_ptype = 21;
			m_ltype = 21;
		}
		else if (m_name == "MOEA_F9") {
			m_domain.setRange(0., 1., 0);
			for (size_t i = 1; i < m_number_variables; ++i)
				m_domain.setRange(-1., 1., i);
			m_dtype = 1;
			m_ptype = 22;
			m_ltype = 22;
		}
		else
			throw MyExcept("The name of MOEA_F is wrong.");

		LoadParetoFront();
	}

	void MOEA_F::alphaFunction(Real alpha[], const Real *x, int dim, int type) {
		if (dim == 2)
		{
			if (type == 21) {
				alpha[0] = x[0];
				alpha[1] = 1 - std::sqrt(x[0]);
			}

			if (type == 22) {
				alpha[0] = x[0];
				alpha[1] = 1 - x[0] * x[0];
			}

			if (type == 23) {
				alpha[0] = x[0];
				alpha[1] = 1 - std::sqrt(alpha[0]) - alpha[0] * sin(10 * alpha[0] * alpha[0] * OFEC_PI);
			}

			if (type == 24) {
				alpha[0] = x[0];
				alpha[1] = 1 - x[0] - 0.05 * sin(4 * OFEC_PI * x[0]);
			}

		}
		else
		{

			if (type == 31) {
				alpha[0] = cos(x[0] * OFEC_PI / 2) * cos(x[1] * OFEC_PI / 2);
				alpha[1] = cos(x[0] * OFEC_PI / 2) * sin(x[1] * OFEC_PI / 2);
				alpha[2] = sin(x[0] * OFEC_PI / 2);
			}

			if (type == 32) {
				alpha[0] = 1 - cos(x[0] * OFEC_PI / 2) * cos(x[1] * OFEC_PI / 2);
				alpha[1] = 1 - cos(x[0] * OFEC_PI / 2) * sin(x[1] * OFEC_PI / 2);
				alpha[2] = 1 - sin(x[0] * OFEC_PI / 2);
			}

			if (type == 33) {
				alpha[0] = x[0];
				alpha[1] = x[1];
				alpha[2] = 3 - (sin(3 * OFEC_PI * x[0]) + sin(3 * OFEC_PI * x[1])) - 2 * (x[0] + x[1]);
			}

			if (type == 34) {
				alpha[0] = x[0] * x[1];
				alpha[1] = x[0] * (1 - x[1]);
				alpha[2] = (1 - x[0]);
			}
		}
	}

	Real MOEA_F::betaFunction(const std::vector<Real> &x, int type) {
		Real beta=0.;
		size_t dim = x.size();

		if (dim == 0) {
			// a bug here when dim=0
			//beta = 0;
			return beta;
		}

		if (type == 1) {
			for (size_t i = 0; i < dim; i++) {
				beta += x[i] * x[i];
			}
			beta = 2.0 * beta / dim;
		}

		if (type == 2) {
			for (size_t i = 0; i < dim; i++) {
				beta += std::sqrt(i + 1) * x[i] * x[i];
			}
			beta = 2.0 * beta / dim;
		}

		if (type == 3) {
			Real sum = 0, xx;
			for (size_t i = 0; i < dim; i++) {
				xx = 2 * x[i];
				sum += (xx * xx - cos(4 * OFEC_PI * xx) + 1);
			}
			beta = 2.0 * sum / dim;
		}

		if (type == 4) {
			Real sum = 0, prod = 1, xx;
			for (size_t i = 0; i < dim; i++) {
				xx = 2 * x[i];
				sum += xx * xx;
				prod *= cos(10 * OFEC_PI * xx / std::sqrt(i + 1));
			}
			beta = 2.0 * (sum - 2 * prod + 2) / dim;
		}

		return beta;
	}

	Real MOEA_F::paretoSetFunc2(const Real &x, const Real &t1, size_t dim, int type, int css) {
		// type:  the type of curve 
		// css:   the class of index
		Real beta = 0.;
		size_t numDim = numberVariables();

		dim++;

		if (type == 21) {
			//Real xy = 2 * (x - 0.5);
			// a bug here when numDim=2
			if (numDim == 2) 
				beta = x - pow(t1, 2.0);
			else	
				beta = x - pow(t1, 0.5 * (1.+3.*(dim-2)/(numDim-2.)));

		}

		if (type == 22) {
			Real theta = 6 * OFEC_PI * t1 + dim * OFEC_PI / numDim;
			//Real xy = 2 * (x - 0.5);
			beta = t1 - sin(theta);
		}

		if (type == 23) {
			Real theta = 6 * OFEC_PI * t1 + dim * OFEC_PI / numDim;
			Real ra = 0.8 * t1;
			//Real xy = 2 * (x - 0.5);
			if (css == 1)
				beta = x - ra * cos(theta);
			else {
				beta = x - ra * sin(theta);
			}
		}

		if (type == 24) {
			Real theta = 6 * OFEC_PI * t1 + dim * OFEC_PI / numDim;
			//Real xy = 2 * (x - 0.5);
			Real ra = 0.8 * t1;
			if (css == 1)
				beta = x - ra * cos(theta / 3);
			else {
				beta = x - ra * sin(theta);
			}
		}

		if (type == 25) {
			Real rho = 0.8;
			Real phi = OFEC_PI * t1;
			Real theta = 6 * OFEC_PI * t1 + dim * OFEC_PI / numDim;
			//Real xy = 2 * (x - 0.5);
			if (css == 1)
				beta = x - rho * sin(phi) * sin(theta);
			else if (css == 2)
				beta = x - rho * sin(phi) * cos(theta);
			else
				beta = x - rho * cos(phi);
		}

		if (type == 26) {
			Real theta = 6 * OFEC_PI * t1 + dim * OFEC_PI / numDim;
			Real ra = 0.3 * t1 * (t1 * cos(4 * theta) + 2);
			//Real xy = 2 * (x - 0.5);
			if (css == 1)
				beta = x - ra * cos(theta);
			else {
				beta = x - ra * sin(theta);
			}
		}

		return beta;
	}

	Real MOEA_F::paretoSetFunc3(const Real &x, const Real &t1, const Real &t2, int dim, int type) {
		// type:  the type of curve 
		// css:   the class of index
		Real beta=0.;
		size_t numDim = numberVariables();

		dim++;

		if (type == 31) {
			Real xy = 4 * (x - 0.5);
			Real rate = 1.0 * dim / numDim;
			beta = xy - 4 * (t1 * t1 * rate + t2 * (1.0 - rate)) + 2;
		}

		if (type == 32) {
			Real theta = 2 * OFEC_PI * t1 + dim * OFEC_PI / numDim;
			Real xy = 4 * (x - 0.5);
			beta = xy - 2 * t2 * sin(theta);
		}

		return beta;
	}

	void MOEA_F::evaluateObjective(Real *x, std::vector<Real> &y_obj) {
		// 2-objective case
		Real *x_var = x;
		size_t nobj = numberObjectives();
		size_t nDim = numberVariables();
		if (nobj == 2)
		{
			if (m_ltype == 21 || m_ltype == 22 || m_ltype == 23 || m_ltype == 24 || m_ltype == 26)	//using in F1-F5
			{
				Real g = 0, h = 0, a, b;
				std::vector <Real> aa;
				std::vector <Real> bb;
				for (size_t n = 1; n < nDim; n++)
				{

					if (n % 2 == 0) {
						a = paretoSetFunc2(x_var[n], x_var[0], n, m_ltype, 1);  // linkage
						aa.push_back(a);
					}
					else
					{
						b = paretoSetFunc2(x_var[n], x_var[0], n, m_ltype, 2);
						bb.push_back(b);
					}

				}

				g = betaFunction(aa, m_dtype);
				h = betaFunction(bb, m_dtype);

				Real alpha[2];
				alphaFunction(alpha, x_var, 2, m_ptype);  // shape function
				y_obj[0] = alpha[0] + g;
				y_obj[1] = alpha[1] + h;
				aa.clear();
				bb.clear();
			}

			if (m_ltype == 25)
			{
				Real g = 0, h = 0, a, b;
				Real e = 0, c;
				std::vector <Real> aa;
				std::vector <Real> bb;
				for (size_t n = 1; n < nDim; n++) {
					if (n % 3 == 0) {
						a = paretoSetFunc2(x_var[n], x_var[0], n, m_ltype, 1);
						aa.push_back(a);
					}
					else if (n % 3 == 1)
					{
						b = paretoSetFunc2(x_var[n], x_var[0], n, m_ltype, 2);
						bb.push_back(b);
					}
					else {
						c = paretoSetFunc2(x_var[n], x_var[0], n, m_ltype, 3);
						if (n % 2 == 0)    aa.push_back(c);
						else          bb.push_back(c);
					}
				}
				g = betaFunction(aa, m_dtype);          // distance function
				h = betaFunction(bb, m_dtype);
				Real alpha[2];
				alphaFunction(alpha, x_var, 2, m_ptype);  // shape function
				y_obj[0] = alpha[0] + h;
				y_obj[1] = alpha[1] + g;
				aa.clear();
				bb.clear();
			}
		}


		// 3-objective case
		if (nobj == 3)
		{
			if (m_ltype == 31 || m_ltype == 32)
			{
				Real g = 0, h = 0, e = 0, a;
				std::vector <Real> aa;
				std::vector <Real> bb;
				std::vector <Real> cc;
				for (int n = 2; n < nDim; n++)
				{
					a = paretoSetFunc3(x_var[n], x_var[0], x_var[1], n, m_ltype);
					if (n % 3 == 0)	    aa.push_back(a);
					else if (n % 3 == 1)	bb.push_back(a);
					else            cc.push_back(a);
				}

				g = betaFunction(aa, m_dtype);
				h = betaFunction(bb, m_dtype);
				e = betaFunction(cc, m_dtype);

				Real alpha[3];
				alphaFunction(alpha, x_var, 3, m_ptype);  // shape function
				y_obj[0] = alpha[0] + h;
				y_obj[1] = alpha[1] + g;
				y_obj[2] = alpha[2] + e;
				aa.clear();
				bb.clear();
				cc.clear();
			}
		}
	}

	void MOEA_F::LoadParetoFront() {
		std::ifstream infile;
		std::stringstream os;
		os << g_working_dir << "/instance/problem/continuous/multi_objective/moea_f/data/pf_P" << m_ptype << "D" << m_dtype << "L" << m_ltype << ".dat";
		infile.open(os.str());
		if (!infile)
			return;
		int lines = 0;
		m_optima.reset(new Optima<>());
		std::string str;
		while (getline(infile, str))
			++lines;
		m_optima.reset(new Optima<>);
		m_optima->resizeObjectiveSet(lines);
		infile.close();
		infile.clear();
		infile.open(os.str());
		std::vector<Real> temp_obj(m_number_objectives);
		for (size_t i = 0; i < lines; i++) {
			for (size_t j = 0; j < m_number_objectives; j++)
				infile >> temp_obj[j];
			m_optima->setObjective(temp_obj, i);
		}
		m_optima->setObjectiveGiven(true);
		infile.close();
	}
}