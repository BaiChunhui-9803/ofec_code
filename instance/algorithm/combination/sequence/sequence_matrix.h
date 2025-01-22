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

#ifndef SEQUENCE_MATRIX_H
#define SEQUENCE_MATRIX_H

/*
* 	template<typename TSequenceOperator>
	class AdaptorGLSeq: public AdaptorGL<typename TSequenceOperator::SolutionType>{

	public:
		using OpType =  TSequenceOperator;
		using SolutionType = typename TSequenceOperator::SolutionType;
		using InterpreterType = typename TSequenceOperator::InterpreterType;
		
*/

namespace ofec {
	template<typename TSequenceInterpreter>
	class SequenceMatrix  {
	public:
		enum class UpdateType = { EDensity,EFre,EACO };
		using SolutionType = typename TSequenceInterpreter::SolutionType;
		using InterpreterType = typename TSequenceInterpreter;

	public:

		void initialize(Problem *pro, Algorithm *alg, Random *rnd,double initVal=1.0) ;
		void resize(const std::vector<int>& mat_size);

	protected:
		double m_evaporate_rate = 0.0;
		
		std::vector<std::vector<double>> m_matrix;
	};

}

#endif