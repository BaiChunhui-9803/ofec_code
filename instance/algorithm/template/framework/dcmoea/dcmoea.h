#ifndef DCMOEA_H
#define DCMOEA_H
#include "../../../../../core/algorithm/Solution.h"

namespace ofec {
	template<typename TSolution>
	class DCMOEA_ind : public TSolution {
	public:
		template<typename ... Args>
		DCMOEA_ind(size_t num_obj, size_t num_con, Args&& ... args) :
			TSolution(num_obj, num_con, std::forward<Args>(args)...) {}
		void resize_vio_obj(size_t size_vio) { m_violation_objectives.clear(); m_violation_objectives.resize(size_vio); }
		void set_efeasible(bool flag) { m_efeasible = flag; }
		bool get_efeasible() { return m_efeasible; }
		void set_vio_obj(size_t id, Real val) { m_violation_objectives[id] = val; }
		std::vector<Real> &get_vio_obj() { return m_violation_objectives; }
		DCMOEA_ind &operator=(const DCMOEA_ind &rhs) {
			if (this == &rhs) return *this;
			TSolution::operator=(rhs);
			m_violation_objectives = rhs.m_violation_objectives;
			m_efeasible = rhs.m_efeasible;
		}
	protected:
		std::vector<Real> m_violation_objectives;
		bool m_efeasible;
	};


	template<typename TSolution>
	class DCMOEA {
	public:
		void initialize(size_t num_cons);
	protected:
		void calculate_initial_max_violation(std::vector<std::unique_ptr<TSolution>> &pop);
		void mark_Solution_efeasible(std::vector<std::unique_ptr<TSolution>> &pop);
		bool judge_population_efeasible(std::vector<std::unique_ptr<TSolution>> &pop);
		void reduce_boundary();
		virtual void calculate_violation_objective(std::vector<std::unique_ptr<TSolution>> &pop);
	protected:
		Real m_z = 1.0e-08;
		Real m_Nearzero = 1.0e-15;
		std::vector<Real> m_max_G, m_e;
		int m_k = 0, m_max_K = 100;
		size_t m_number_constraints;
	};

	template<typename TSolution>
	void DCMOEA<TSolution>::initialize(size_t num_cons) {
		m_number_constraints = num_cons;
		m_max_G.assign(num_cons, 1);
		m_e = m_max_G;
	}

	template<typename TSolution>
	void DCMOEA<TSolution>::calculate_initial_max_violation(std::vector<std::unique_ptr<TSolution>> &pop) {
		for (int i = 0; i < pop.size(); ++i) {
			for (int k = 0; k < pop[0]->constraint().size(); ++k) {
				if (pop[i]->constraint()[k] > m_max_G[k]) {
					m_max_G[k] = pop[i]->constraint()[k];
				}
			}
		}
	}

	template<typename TSolution>
	void DCMOEA<TSolution>::mark_Solution_efeasible(std::vector<std::unique_ptr<TSolution>> &pop) {
		for (int i = 0; i < pop.size(); ++i) {
			bool test = true;
			for (int j = 0; j < m_e.size(); ++j) {
				if (pop[i]->constraint()[j] > m_e[j]) {
					test = false;
					break;
				}
			}
			pop[i]->set_efeasible(test);
		}
	}

	template<typename TSolution>
	bool DCMOEA<TSolution>::judge_population_efeasible(std::vector<std::unique_ptr<TSolution>> &pop) {
		bool flag = true;
		for (int i = 0; i < pop.size(); ++i) {
			if (pop[i]->get_efeasible() == false) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	template<typename TSolution>
	void DCMOEA<TSolution>::reduce_boundary() {
		for (int i = 0; i < m_max_G.size(); ++i) {
			auto temp = sqrt(log2((m_max_G[i] + m_z) / m_z) / log2(exp(1)));
			auto B = m_max_K / temp;
			if (B == 0.0) {
				B = m_Nearzero;
			}
			auto A = m_max_G[i] + m_z;
			auto f = A * exp(-pow(m_k / B, 2)) - m_z;
			if (abs(f - m_z) < m_Nearzero)
				f = m_z;
			if (abs(f) < m_Nearzero)
				f = 0.0;
			m_e[i] = f;
		}
	}

	template<typename TSolution>
	void DCMOEA<TSolution>::calculate_violation_objective(std::vector<std::unique_ptr<TSolution>> &pop) {
		for (int i = 0; i < pop.size(); ++i) {
			Real vObj = 0.0;
			for (int h = 0; h < m_number_constraints; ++h) {
				vObj += (pop[i]->constraint()[h] / m_max_G[h]);
			}
			pop[i]->set_vio_obj(0, vObj / m_number_constraints);
		}
	}
}




#endif // !DCO_H

