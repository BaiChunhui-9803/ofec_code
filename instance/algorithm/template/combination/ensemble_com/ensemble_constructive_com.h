#ifndef ENSEMBLE_CONSTRUCTIVE_COM_H
#define ENSEMBLE_CONSTRUCTIVE_COM_H

#include"constructive_com_operator.h"
#include<functional>
namespace ofec {

	template<typename TInterpreter>
	class ensemble_constructive_com {
	public:
		using IndividualType = typename TInterpreter::SolutionType;

	protected:
		std::vector<std::unique_ptr<ConstructiveComOp<TInterpreter>>> m_ops;
		std::function<Real(int cur_iter)> m_init_pop_pro_fun;
		std::function<Real(int cur_iter)> m_step_init_pro_fun;
		std::function<Real(int cur_iter)> m_step_next_pro_fun;
		std::function<Real(int cur_iter)> m_step_final_pro_fun;
		std::function<Real(int cur_iter)> m_update_final_pro_fun;

		std::vector<int> m_updateFinal_op_idxs;
	//	std::vector<int> m_
	public:

		ensemble_constructive_com() = default;
		virtual ~ensemble_constructive_com() = default;

		virtual void initialize() {
			for (auto& it : m_ops) {
				it->initialize();
			}
			m_init_pop_pro_fun = [&](int cur) {
				if (cur < 0 || cur >= m_ops.size())return 0;
				else {
					return m_ops[cur]->getInitPopPro();
				}
			};
			m_step_init_pro_fun = [&](int cur) {
				if (cur < 0 || cur >= m_ops.size())return 0.0;
				else {
					return m_ops[cur]->getStepInitPro();
				}
			};
			m_step_next_pro_fun = [&](int cur) {
				if (cur < 0 || cur >= m_ops.size())return 0.0;
				else {
					return m_ops[cur]->getStepNextPro();
				}
			};
			m_step_final_pro_fun = [&](int cur) {
				if (cur < 0 || cur >= m_ops.size())return 0.0;
				else {
					return m_ops[cur]->getStepFinalPro();
				}
			};
			m_update_final_pro_fun = [&](int cur) {
				if (cur < 0 || cur >= m_ops.size())return 0.0;
				else {
					return m_ops[cur]->getUpdateFinalPro();
				}
			};
		}

		virtual void initialize_solutions(std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) {
			std::vector<int> select_order;
			m_random->uniform.sequencePro<int>(0, m_ops.size(), m_init_pop_pro_fun);
			for (auto iter : select_order) {
				if (m_ops[iter]->initializePop(pops, interperter, alg, pro, rnd)) {
					break;
				}
			}

			for (auto& iter_op : m_ops) {
				iter_op->constructiveFinal(pops, interperter, alg, pro, rnd);
			}
			
		}
		virtual void construct_solutions(std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter,Algorithm *alg, Problem *pro, Random *rnd) {
			std::vector<int> select_order;
			for (int pop_id(0); pop_id < pops.size(); ++pop_id) {
				interperter.stepInit(*pops[pop_id]);
				m_random->uniform.sequencePro<int>(0, m_ops.size(), m_step_init_pro_fun,select_order);
				for (auto iter : select_order) {
					if (m_ops[iter]->stepInit(*pops[pop_id], pop_id, interperter, alg, pro, rnd))break;
				}
				
				while (!interperter.stepFinish(*pops[pop_id])) {
					if (interperter.stepFeasible(*pops[pop_id])) {

						m_random->uniform.sequencePro<int>(0, m_ops.size(), m_step_next_pro_fun, select_order);
						for (auto iter : select_order) {
							if (m_ops[iter]->stepNext(*pops[pop_id], pop_id, interperter, alg, pro, rnd))break;
						}
					}
					else {
						interperter.stepBack(*pops[pop_id]);
					}
				}
				interperter.stepFinal(*pops[pop_id]);

			}
			
		}
		template<
			template <class> class Op
		>
		void addOp() {
			m_ops.emplace_back(new Op<TInterpreter>())
		}


		
	};

}



#endif