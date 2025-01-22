#include "sequence_actions_ensemble.h"

void ofec::SequenceActionEnsembleParametersBase::initGL_GA(Problem *pro,Algorithm *alg)
{

	this = alg;
	m_problem.get() = pro;
	int totalA(2);
	m_pars.resize(totalA);
	m_action_pros.resize(totalA);
	m_action_types.resize(totalA);
	m_global_act.resize(totalA);
	for (auto& it : m_global_act) it = true;

	


	{
	//		int m_problem.get() = -1;
	//	int this = -1;
		std::unique_ptr<SequenceGeneticMatrixParameter> cur_ptr(new SequenceGeneticMatrixParameter);
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;
		
		m_pars[0].reset(cur_ptr.release());
		m_action_pros[0] = 0.2;
		m_action_types[0] = SequenceAlgType::kSequenceGeneticMatrix;
		
	}

	{
		//		int m_problem.get() = -1;
		//	int this = -1;
		std::unique_ptr<SequenceParametersBase> cur_ptr(new SequenceParametersBase());
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;

		m_pars[1].reset(cur_ptr.release());
		m_action_pros[1] = 0.8;
		m_action_types[1] = SequenceAlgType::kSequenceSelfLearn;

	}   
}

void ofec::SequenceActionEnsembleParametersBase::initGL_GA_LS(Problem *pro, Algorithm *alg)
{

	this = alg;
	m_problem.get() = pro;
	int totalA(3);
	m_pars.resize(totalA);
	m_action_pros.resize(totalA);
	m_action_types.resize(totalA);

	m_global_act.resize(totalA);
	for (auto& it : m_global_act) it = true;

	int action_idx = 0;

	{
		std::unique_ptr<SequenceGeneticMatrixParameter> cur_ptr(new SequenceGeneticMatrixParameter);
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;

		m_pars[action_idx].reset(cur_ptr.release());
		m_action_pros[action_idx] = 0.2;
		m_action_types[action_idx] = SequenceAlgType::kSequenceGeneticMatrix;

	}
	action_idx++;
	{
		//		int m_problem.get() = -1;
		//	int this = -1;
		std::unique_ptr<SequenceParametersBase> cur_ptr(new SequenceParametersBase());
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;

		m_pars[action_idx].reset(cur_ptr.release());
		m_action_pros[action_idx] = 0.8;
		m_action_types[action_idx] = SequenceAlgType::kSequenceSelfLearn;

	}


	action_idx++;
	{
		//		int m_problem.get() = -1;
		//	int this = -1;
		std::unique_ptr<SequenceParametersBase> cur_ptr(new SequenceParametersBase());
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;

		m_pars[action_idx].reset(cur_ptr.release());
		m_action_pros[action_idx] = 0;
		m_action_types[action_idx] = SequenceAlgType::kSequenceSelfLearn;

	}
}


void ofec::SequenceActionEnsembleParametersBase::initGA(Problem *pro, Algorithm *alg) {


	this = alg;
	m_problem.get() = pro;
	int totalA(1);
	m_pars.resize(totalA);
	m_action_pros.resize(totalA);
	m_action_types.resize(totalA);

	m_global_act.resize(totalA);
	for (auto& it : m_global_act) it = true;

	{
		//		int m_problem.get() = -1;
		//	int this = -1;
		std::unique_ptr<SequenceParametersBase> cur_ptr(new SequenceParametersBase());
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;

		m_pars[0].reset(cur_ptr.release());
		m_action_pros[0] = 1.0;
		m_action_types[0] = SequenceAlgType::kSequenceSelfLearn;

	}
}

void ofec::SequenceActionEnsembleParametersBase::initGL_GlobalNone(Problem *pro, Algorithm *alg)
{

	this = alg;
	m_problem.get() = pro;
	int totalA(1);
	m_pars.resize(totalA);
	m_action_pros.resize(totalA);
	m_action_types.resize(totalA);
	m_global_act.resize(totalA);
	for (auto& it : m_global_act) it = true;




	{
		//		int m_problem.get() = -1;
		//	int this = -1;
		std::unique_ptr<SequenceGeneticMatrixParameter> cur_ptr(new SequenceGeneticMatrixParameter);
		cur_ptr->m_problem.get() = pro;
		cur_ptr->this = alg;
		cur_ptr->m_updateTag = SequenceGeneticMatrixParameter::MatrixUpdateTag::kNone;

		m_pars[0].reset(cur_ptr.release());
		m_action_pros[0] = 1.0;
		m_action_types[0] = SequenceAlgType::kSequenceGeneticMatrix;

	}

}
