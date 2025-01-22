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

#ifndef OFEC_POPULATION_ACO_H
#define OFEC_POPULATION_ACO_H

#include<vector>
#include"../../../../core/algorithm/population.h"
#include<functional>
#include<memory>

namespace ofec {
	template<typename T_ACO>
	class PopACO : public Population<typename T_ACO::SolutionType> {
	public:
		using AdaptorType = typename T_ACO::AdaptorType;
		using OpType = typename T_ACO::OpType;
		using SolutionType = typename T_ACO::SolutionType;
		using InterpreterType = typename T_ACO::InterpreterType;

	protected:
		InterpreterType m_pro_interperter;
		T_ACO m_ACO_model;

	protected:
		virtual inline Real calFitness(Problem *pro, SolutionType& sol) const {
			if (pro->optimizeMode()[0] == OptimizeMode::kMaximize) {
				return   sol.objective(0);
			}
			else {
				return (1.0 / (1e-9 + sol.objective(0)));
			}
		}
		virtual int updateSolutions(Random *rnd,Problem *pro, Algorithm *alg);
		virtual int localSearch(Random *rnd,Problem *pro, Algorithm *alg) {
			int rf;
			for (auto& it : m_individuals) {
				rf |= m_pro_interperter.evaluate(pro, alg, rnd, *it);
			//	if (rf != kNormalEval) return rf;
			}
		//	updateFitness(pro);
			return rf;
		}

	public:
		const IntrepreterType& getProInterpreter() {
			return m_pro_interperter;
		}
		const std::vector<std::vector<ofec::Real>>& getMatrix() {
			return m_ACO_model.getMatrix();
		}
		
		virtual int evolve(Random *rnd,Problem *pro, Algorithm *alg);
		virtual void initializeParameters(int id_param, Random *rnd, Problem *pro, Algorithm *alg);
		virtual void initialize(Random *rnd,Problem *pro,Algorithm *alg); //a uniformly distributed initialization by default

	};

	template<typename T_ACO>
	inline int PopACO<T_ACO>::updateSolutions(Random *rnd,Problem *pro, Algorithm *alg) {
		for (auto& it : m_individuals) {
			m_pro_interperter.stepInit(pro,*it);
			while (!m_pro_interperter.stepFinish(pro,*it)) {
				if (m_pro_interperter.stepFeasible(pro,*it)) {
					if (m_ACO_model.stepNext(rnd, pro, alg, m_pro_interperter, *it)) {
						m_ACO_model.localUpdatePheromone(rnd, pro, alg, m_pro_interperter, *it);
					}
					else {
						m_pro_interperter.stepBack(pro,*it);
					}
				}
				else {
					m_pro_interperter.stepBack(pro,*it);
				}
			}
			m_pro_interperter.stepFinal(pro,*it);
		}
		return localSearch(rnd,pro, alg);
	}

	template<typename T_ACO>
	inline int PopACO<T_ACO>::evolve(Random *rnd,Problem *pro, Algorithm *alg) {
		int tag = updateSolutions(pro, alg, rnd);
		if (tag != kNormalEval) {
			handleEvaluationTag(tag);
		}
		m_ACO_model.globalUpdatePheromone(rnd, pro, alg, m_pro_interperter, m_individuals);
		return tag;
	}

	template<typename T_ACO>
	inline void PopACO<T_ACO>::initializeParameters(int id_param, Random *rnd, Problem *pro,Algorithm *alg) {
		Population<SolutionType>::clear();
		Population<SolutionType>::resize(GET_PARAM(id_param).get<int>("population size"), pro, pro->numberVariables());
		m_pro_interperter.initializeByProblem(pro);
		auto heuristic_sol(m_pro_interperter.heuristicSol(pro, alg, rnd));
//		updateFitness(pro, heuristic_sol);
		m_ACO_model.fitnessUpdate() =
			std::bind(&PopACO::calFitness, this,  std::placeholders::_1, std::placeholders::_2);
		m_ACO_model.initialize(id_param, rnd, pro, alg, m_pro_interperter);

		//m_ACO_model.setHeuristicFitness(heuristic_sol.fitness());
	}

	template<typename T_ACO>
	inline void PopACO<T_ACO>::initialize(Random *rnd,Problem *pro,Algorithm *alg)
	{
		
		Population::initialize(pro, rnd);
		for (auto& it : m_individuals) {
			m_pro_interperter.stepFinal(pro,*it);
			m_pro_interperter.evaluate(pro, alg, rnd, *it);
			//it->evaluate(pro, alg);
		}
	}
}

#endif