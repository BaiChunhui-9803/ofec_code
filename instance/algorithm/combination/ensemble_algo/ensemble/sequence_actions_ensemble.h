/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* Framework of genetic learning (PopGL) algorithm
*
*********************************************************************************/

#ifndef SEQUENCE_ACTIONS_ENSEMBLE_H
#define SEQUENCE_ACTIONS_ENSEMBLE_H
#include"../../../template/combination/sequence/sequence_action.h"
#include"../gl_actions/sequence_frequency.h"
#include"../gl_actions/sequence_genetic_matrix.h"
#include"../ga_actions/sequence_learn_other.h"
#include"../ga_actions/sequence_self_learn.h"
#include "../global_learn/sequence_local_search.h"
#include<memory>

#include "../../sequence/SP/sp_interpreter.h"
#include "../../sequence/SP/sp_operator.h"

namespace ofec {
	enum class SequenceAlgType {
		kSequenceFrequency = 0,
		kSequenceGeneticMatrix,
		kSequenceSelfLearn,
		kSequenceLearnOther,
		kLocalSearch
	};

	class SequenceActionEnsembleParametersBase:public SequenceParametersBase {
	public:
		std::vector<std::unique_ptr<SequenceParametersBase>> m_pars;
		std::vector<double> m_action_pros;
		std::vector<SequenceAlgType> m_action_types;
		std::vector<bool> m_global_act;

		void initGL_GA(Problem *pro,Algorithm *alg);
		void initGL_GA_LS(Problem *pro, Algorithm *alg);
		void initGA(Problem *pro, Algorithm *alg);

		void initGL_GlobalNone(Problem *pro, Algorithm *alg);
	};
	

	template<typename TSequenceOperator>
	class SequenceActionEnsemble :public SequenceAction<TSequenceOperator> {
	public:
		using OpType = TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		using ParType = typename  SequenceAction<TSequenceOperator>::ParType;


	public:

		virtual void initialize(const SequenceParametersBase& par) override;
		virtual void resize(const std::vector<int>& mat_size, double* initVal = nullptr) {
			for (auto& it : m_actions) {
				it->resize(mat_size, initVal);
			}
		}

		virtual void setPopId(int popId)override {
			m_id_pop = popId;
			for (auto& it : m_actions) {
				it->setPopId(popId);
			}
		}
		virtual void setEvalStrategy(const std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra)override {
			m_eval_stra = eval_stra;
			for (auto& it : m_actions) {
				it->setEvalStrategy(m_eval_stra);
			}
		}

		// return -1 means failed
		virtual int learn(
			SolutionType& cur,
			Random *rnd, Problem *pro,
			const InterpreterType& interpreter,
			const OpType& op,
			const std::vector<SolutionType*>& pop,
			int id_indiv
		)const override {
			std::vector<double> temp_pro(m_actions_pro);
			int left = temp_pro.size();
			int next(-1);
			int action_idx(-1);
			while (left-- || next != -1) {
				action_idx = rnd->uniform.spinWheel(temp_pro.begin(), temp_pro.end()) - temp_pro.begin();
				if (action_idx != temp_pro.size())
				{
					next = m_actions[action_idx]->learn(cur,rnd, pro, interpreter, op, pop, id_indiv);
					if (next == -1) temp_pro[action_idx] = 0;
					else return next;
				}
				else {
					next = -1;
				}
			}
			return next;
		}

		virtual void globalUpdate(Problem *pro, Random *rnd,
			const std::vector<SolutionType*>& pop) override {
			for (auto& it : m_actions) {
				it->globalUpdate(pro,rnd,pop);
			}
		}
		virtual void localUpdate(Problem *pro,
			const SolutionType& cur) override{
			for (auto& it : m_actions) {
				it->localUpdate(pro, cur);
			}
		}
		//std::unique_ptr<SequenceAction>& getActions(int idx) {
		//	return m_actions[idx];
		//}

		virtual int createSolution(Problem *pro, Algorithm *alg, Random *rnd,
			std::vector<SolutionType*>& pop,
			std::vector<SolutionType>& offspring);


		virtual void setPro(std::vector<double>& pro) {
			m_actions_pro = pro;
		}

	protected:
		std::vector<std::unique_ptr<SequenceAction>> m_actions;
		std::vector<double> m_actions_pro;
		std::vector<double> m_temp_pro;
		std::vector<bool> m_global_act;
		
	};


	template<typename TSequenceOperator>
	inline int SequenceActionEnsemble<TSequenceOperator>::createSolution(Problem *pro, Algorithm *alg, Random *rnd,
		std::vector<SolutionType*>& pop,
		std::vector<SolutionType>& offspring) {
	
		if (m_flagStepLearn) {
			auto& op(GET_ASeq(alg).getOp<TSequenceOperator>());
			auto& interpreter(GET_ASeq(alg).getInterpreter<InterpreterType>());
			for (int offId(0); offId < offspring.size(); ++offId)
			{
				auto& it = offspring[offId];
				it.reset();
				interpreter.stepInit(pro, it);
				while (!interpreter.stepFinish(pro, it)) {
					int next = learn(it, rnd, pro, interpreter, op, pop, offId);
					if (next != -1) {
						interpreter.stepNext(pro, it, next);
						localUpdate(pro, it);
					}
					else {
						interpreter.stepBack(pro, it);
					}
				}
				//	op.learn_from_other(it, rnd, pro, interpreter, m_center, m_radius);
				interpreter.stepFinal(pro, it);
			}

		}


		int rf = 0;

		for (int idx(0); idx < m_actions.size(); ++idx) {
			if (m_global_act[idx]) {
				rf |= m_actions[idx]->globalLearn(pro, alg, rnd,
					pop, offspring);
			}
		}
		return rf;
	}



	template<typename TSequenceInterpreter>
	inline void SequenceActionEnsemble<TSequenceInterpreter>::initialize(const SequenceParametersBase& par){
		auto& cur_par = dynamic_cast<const SequenceActionEnsembleParametersBase&>(par);
		m_actions_pro = cur_par.m_action_pros;
		m_actions.resize(cur_par.m_action_types.size());
		m_global_act = cur_par.m_global_act;
		for (int idx(0); idx < m_actions.size(); ++idx) {
			switch (cur_par.m_action_types[idx])
			{
			case ofec::SequenceAlgType::kSequenceFrequency:
				m_actions[idx].reset(new SequenceFrequency<TSequenceInterpreter>());
				break;
			case ofec::SequenceAlgType::kSequenceGeneticMatrix:
				m_actions[idx].reset(new SequenceGeneticMatrix<TSequenceInterpreter>());

				break;
			case ofec::SequenceAlgType::kSequenceSelfLearn:
				m_actions[idx].reset(new SequenceSelfLearn<TSequenceInterpreter>());

				break;
			case ofec::SequenceAlgType::kSequenceLearnOther:
				m_actions[idx].reset(new SequenceLearnOther<TSequenceInterpreter>());

				break;

			case ofec::SequenceAlgType::kLocalSearch:
				m_actions[idx].reset(new SequenceLocalSearch<TSequenceInterpreter>());

				break;
			default:
				break;
			}

			m_actions[idx]->initialize(*cur_par.m_pars[idx].get());

		}

		m_flagStepLearn = false;
		for (auto& it : m_actions) {
			m_flagStepLearn |= it->flagStepLearn();
		}

		//m_actions.resize(3);

		//m_actions[0].reset(new SequenceSelfLearn());
		//m_actions[0].reset(new SequenceSelfLearn());
	}

	//template<typename TSequenceOperator>
	//inline int SequenceActionEnsemble<TSequenceOperator>::updateSolution(Problem *pro, Algorithm *alg, Random *rnd, std::vector<SolutionType*>& pop, std::vector<SolutionType>& offspring, int& num_improve)
	//{
	//	int rf = 0;
	//	num_improve = 0;
	//	for (auto& it : m_actions) {
	//		rf|= it->updateSolution(pro, alg, rnd,
	//			pop,offspring, num_improve);
	//	}
	//	return rf;
	//}


	//extern void getParameters(std::vector)
}

#endif