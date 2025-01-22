/********* Begin Register Information **********
[
	{ "name": "GLU_SP", "identifier": "GLU_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP_GLU_SP", "identifier": "MP_GLU_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "AMP_GLU_SP", "identifier": "AMP_GLU_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP_GLU_WD_SP", "identifier": "MP_GLU_WD_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP_GLU_WM_SP", "identifier": "MP_GLU_WM_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "GLU_WD_SP", "identifier": "GLU_WD_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "GLU_WM_SP", "identifier": "GLU_WM_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP_GLU_SWM_SP", "identifier": "MP_GLU_SWM_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "GLU_SWM_SP", "identifier": "GLU_SWM_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "AMP_GLU_SWM_SP", "identifier": "AMP_GLU_SWM_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP2E_GLU_SP", "identifier": "MP2E_GLU_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP2E_GLU_SWM_SP", "identifier": "MP2E_GLU_SWM_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "GL_SP", "identifier": "GL_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "MP_GL_SP", "identifier": "MP_GL_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },	{ "name": "AS_SP", "identifier": "AS_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] },
	{ "name": "PACO_SP", "identifier": "PACO_SP", "problem tags": [ "ComOP", "SP", "SOP", "MMOP", "NoisyOP", "DOP" ] }
]
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (ofec)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@163.com
* Language: C++
*************************************************************************
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

Dorigo, M. (1996). "Ant system optimization by a colony of cooperating agents."
IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS.

*************************************************************************/
// Created: 7 Oct 2014
// Last modified:
// updated: by Yiya Diao in July  2021

#ifndef SET_COM_ALG_SP_H
#define SET_COM_ALG_SP_H

//#include "aco/aco.h"
//#include "aco/as/as.h"
//#include "aco/acs/acs.h"
//#include "aco/mmas/mmas.h"
//#include "aco/paco/paco.h"
//#include "sequence/TSP/TSP_interpreter.h"
#include "sequence/SP/SP_interpreter.h"
#include"gl/gl_adaptor_com.h"
#include "gl/gl_com.h"
#include"gl/mp/mp_gl.h"
#include"sequence/SP/SP_operator.h"

#include"ensemble_algo/ensemble/sequence_actions_ensemble.h"
//#include"../uncertianty/noisy/evaluation_strategy_matrix.h"
#include"ensemble_algo/gl/adaptor_gl_seq.h"

#include "../template/framework/uncertianty/algorithm_uncertianty.h"
#include "../template/framework/uncertianty/evaluation_strategy.h"
#include "../template/framework/uncertianty/pop_uncertianty.h"
#include "ensemble_algo/adaptor/adaptor_seq.h"
#include "uncertianty_algorithm/pop_uncertianty_seq.h"
#include "uncertianty_algorithm/pop_uncertainty_seq_ensemble.h"

#include "uncertianty_algorithm/multi_pop/pop_mp_diversity_main_uncertianty_seq.h"


//#include "uncertianty_algorithm/adaptor_gl_uncertianty_seq.h"
//#include "uncertianty_algorithm/gl_pop_uncertainty_seq.h"
//#include "uncertianty_algorithm/gl_uncertianty_seq.h"

#include"uncertianty_algorithm/gl_uncertianty_seq.h"
#include"uncertianty_algorithm/multi_pop/mp_E2_uncertianty_seq.h"


#include "../template/framework/uncertianty/evaluation_strategy.h"


#include "uncertianty_algorithm/multi_pop/mp_uncertianty_seq.h"
#include "uncertianty_algorithm/multi_pop/amp_uncertianty_seq.h"


// distance calculator
#include "../template/framework/uncertianty/distance_calculator.h"
#include "../uncertianty/distance_calculator/distance_calculator_weight.h"
#include "../uncertianty/distance_calculator/distance_calculator_seq_weight.h"
#include "../uncertianty/distance_calculator/distance_calculator_seq.h"



namespace ofec {
	//using AS_TSP = ACO<AS<TSP_Interpreter>>;
	//using ACS_TSP = ACO<ACS<TSP_Interpreter>>;
	//using MMAS_TSP = ACO<MMAS<TSP_Interpreter>>;
	//using PACO_TSP = ACO<PACO<TSP_Interpreter>>;
	//using AS_SP = ACO<AS<SP_Interpreter>>;
	//using GL_SP = GL_Seq<AdaptorGLSeq<SelectionProblemOperator>>;
	//using GL_SP = GL_Seq<AdaptorGLSeqEval<SelectionProblemOperator>>;



	using GLU_SP = GeneticLearningUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorBase,
		SelectionProblemOperator>;


	using GLU_WD_SP = GeneticLearningUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorWeight,
		SelectionProblemOperator>;


	using GLU_WM_SP = GeneticLearningUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorSeqWeight,
		SelectionProblemOperator>;


	using GLU_SWM_SP = GeneticLearningUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorSeqMat,
		SelectionProblemOperator>;



	//using PACO_SP = ACO<AS<SP_Interpreter>>;
	//using MP_GL_SP = MP_GL<AdaptorGLSeq<SelectionProblemOperator>>;
	using MP_GLU_SP = MultiPopUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorBase,
		SelectionProblemOperator>;

	using MP_GLU_WD_SP = MultiPopUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorWeight,
		SelectionProblemOperator>;


	using MP_GLU_WM_SP = MultiPopUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorSeqWeight,
		SelectionProblemOperator>;

	using MP_GLU_SWM_SP = MultiPopUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorSeqMat,
		SelectionProblemOperator>;



	using AMP_GLU_SP = AdaptiveMultiPopUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorBase,
		SelectionProblemOperator>;


	using AMP_GLU_SWM_SP = AdaptiveMultiPopUncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorSeqMat,
		SelectionProblemOperator>;


	using MP2E_GLU_SWM_SP = MultiPopE2UncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorSeqMat,
		SelectionProblemOperator>;



	using MP2E_GLU_SP = MultiPopE2UncertiantySeq<
		EvaluationStrategyBase,
		DistanceCalculatorBase,
		SelectionProblemOperator>;




}

#endif