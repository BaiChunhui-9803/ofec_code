/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@163.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

Dorigo, M. (1996). "Ant system optimization by a colony of cooperating agents."
IEEE TRANSACTIONS ON SYSTEMS, MAN, AND CYBERNETICS.

*************************************************************************/
// Created: 7 Oct 2014
// Last modified:
// updated: by Yiya Diao in July  2021

#ifndef OFEC_ACO_H
#define OFEC_ACO_H

#include "aco_pop.h"
#include "../../../../core/algorithm/algorithm.h"
#include <iostream>

#ifdef  OFEC_DEMO
#include <custom/buffer/algorithm/combination/density_matrix/buffer_density_matrix.h>
#include <core/global_ui.h>
#endif //  OFEC_DEMO


namespace ofec {
	template<typename T_ACO>
	class ACO : public Algorithm{

	public:
	public:
		using SolutionType = typename T_ACO::SolutionType;
		using InterpreterType = typename T_ACO::InterpreterType;
	public:
		using SolutionType = typename T_ACO::SolutionType;
	protected:
		PopACO<T_ACO> m_pops;
	protected:
	void updateFitness(Problem *pro, SolutionType& sol) const;
	public:
		virtual void initialize_()override;
		virtual void run_()override;
		virtual void record() override {}

#ifdef OFEC_DEMO
	public:
		virtual void updateBuffer();
#endif
	};

#ifdef  OFEC_DEMO
	template<typename T_ACO>
	inline void ACO<T_ACO>::updateBuffer() {
		{
			auto& mat(ofec_demo::BufferDensityMatrix::ms_density_matrix[1] = m_pops.getMatrix());
			auto& range(ofec_demo::BufferDensityMatrix::ms_density_range[1] = { std::numeric_limits<Real>::max(),0 });
			for (auto& it : mat) {
				for (auto& it2 : it) {
					range.first = std::min(it2, range.first);
					range.second = std::max(it2, range.second);
				}
			}
			for (auto& it : mat) {
				sort(it.begin(), it.end());
			}
		}
		{
			auto& mat(ofec_demo::BufferDensityMatrix::ms_density_matrix[0]);
			auto& mat_size(m_pops.getProInterpreter().getMatrixSize());
			mat.resize(mat_size.size());
			for (int idx(0); idx < mat.size(); ++idx) {
				mat[idx].resize(mat_size[idx]);
			}
			for (int idx(0); idx < m_pops.size(); ++idx) {
				m_pops.getProInterpreter().updateEdges(m_problem.get(), m_pops[idx], true);
			}
			calculateDensityMatrix(m_pops, mat);
			auto& range(ofec_demo::BufferDensityMatrix::ms_density_range[0] = { std::numeric_limits<Real>::max(),0 });
			for (auto& it : mat) {
				for (auto& it2 : it) {
					range.first = std::min(it2, range.first);
					range.second = std::max(it2, range.second);
				}
			}
			for (auto& it : mat) {
				sort(it.begin(), it.end());
			}

		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif

	template<typename T_ACO>
	inline void ACO<T_ACO>::initialize_()
	{

		Algorithm::initialize_();
		m_pops.initializeParameters(m_id_param, m_random.get(),m_problem.get(), this);
		m_keep_candidates_updated = true;

	}

	template<typename T_ACO>
	inline void ACO<T_ACO>::run_() {
		m_pops.initialize(m_random.get(),m_problem.get(), this);
#ifdef OFEC_DEMO
		updateBuffer();
#endif
		while (!terminating()) {
			m_pops.evolve(m_random.get(), m_problem.get(), this);
#ifdef OFEC_DEMO
			updateBuffer();
#endif
			m_problem->showInfomations(this);	
		}
	}
}

#endif //!OFEC_ACO_H