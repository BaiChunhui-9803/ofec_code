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
// updated: Yiya Diao() 3 Dec 2018
//
#ifndef OFEC_MMAS_H
#define OFEC_MMAS_H

#include "../../../template/classic/ant_colony_optimization/aco_model.h"
#include "../../../../../utility/functional.h"
#include <cmath>

namespace ofec {
	template<typename TInterpreter>
	class MMAS : public ModelACO<TInterpreter> {
	public:
		using SolutionType = typename TInterpreter::SolutionType;
		using IntrepreterType = TInterpreter;

	protected:
		long long m_iteration = 0;
		Real m_alpha = 1.0;
		Real m_beta = 2.0;
		Real m_rho = 0.98;
		Real m_pheroMax;
		Real m_pheroMin;
		Real m_branchFrc = 1.00001;
		Real m_branchFactor;
		long long m_uGB = 25;
		Real m_lambda = 0.05;
		int m_length = 20; //len of candidate list
		long long m_restartFoundBest;
		std::vector<std::unique_ptr<SolutionType>> m_globalBest;//the globally best tour from the beginning of the trial
		std::vector<std::unique_ptr<SolutionType>> m_restartBest;//the restart best tour from the beginning of the trial
		std::vector<std::vector<int>> mvv_neighbors;

	private:
		std::vector<std::vector<bool>> m_best_sols_edges;

	protected:
		Real getPro(const TInterpreter &interpreter, Problem *pro, const SolutionType &cur, int to) override {
			return  pow(mvv_phero[interpreter.curPositionIdx(pro, cur)][to], m_alpha) * pow(interpreter.heursticInfo(pro, cur, to), m_beta);
		}
		void pheroSmoothing();
		void updatePheromonePolicy(const std::vector<std::unique_ptr<SolutionType>> &update_sols);
		void checkPheromoneTrailLimits();
		void updatePheroMinAndMax();
		Real nodeBranching();

		bool updateSolutions(Problem *pro, SolutionType &cur_indiv, std::vector<std::unique_ptr<SolutionType>> &curPops);

	public:
		void initialize(int id_param, Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter &interpreter) override {
			ModelACO::initialize(id_param, rnd, pro, alg, interpreter);
			auto &v = GET_PARAM(id_param);
			/*	if (v.find("ACS_beta") != v.end()) {
					m_beta = v.get<Real>("ACS_beta");
				}
				else m_beta = 5.0;
				if (v.find("ACS_Q") != v.end()) {
					m_Q = v.get<Real>("ACS_Q");
				}
				else m_Q = 100.;
				if (v.find("ACS_rho") != v.end()) {
					m_rho = v.get<Real>("ACS_rho");
				}
				else m_rho = 0.5;*/
			auto &mat(interpreter.getMatrixSize());

			m_best_sols_edges.resize(mat.size());
			for (int idx(0); idx < m_best_sols_edges.size(); ++idx) {
				m_best_sols_edges[idx].resize(mat[idx]);
			}
			// find the nearest candidates
			mvv_neighbors.resize(mat.size());
			for (int idx(0); idx < mat.size(); ++idx) {
				interpreter.getNearestNeighbors(pro, m_length, idx, mvv_neighbors[idx]);
			}

			auto sol(interpreter.heuristicSol(pro, alg, rnd));

			m_pheroMax = 1. / (m_rho * m_fitness_update(pro, sol));
			m_pheroMin = m_pheroMax / (2. * mat.size());

			mvv_phero.resize(mat.size());
			for (int idx(0); idx < mat.size(); ++idx) {
				mvv_phero[idx].resize(mat[idx]);
				std::fill(mvv_phero[idx].begin(), mvv_phero[idx].end(), m_pheroMax);
			}
			m_globalBest.clear();
			m_restartBest.clear();

		}

		void globalUpdatePheromone(Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter &interpreter, std::vector<std::unique_ptr<SolutionType>> &pops) override {
			++m_iteration;
			if (pops.empty())return;
			for (auto &it : pops) {
				interpreter.updateEdges(pro, *it);
			}
			bool isUpdate(false);
			for (auto &it : pops) {
				isUpdate |= updateSolutions(pro, *it, m_globalBest);
			}
			if (isUpdate) {
				updatePheroMinAndMax();
			}
			isUpdate = false;
			for (auto &it : pops) {
				isUpdate |= updateSolutions(pro, *it, m_restartBest);
			}
			if (isUpdate) {
				m_restartFoundBest = m_iteration;
			}
			if (m_iteration % m_uGB) updatePheromonePolicy(m_globalBest);
			else updatePheromonePolicy(m_restartBest);
			checkPheromoneTrailLimits();
			pheroSmoothing();
		}

		bool stepNext(Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter &interpreter, SolutionType &cur) override {
			std::function<Real(int to)> nei_weight_fun = std::bind(&MMAS::getPro, this, std::cref(interpreter), pro, std::cref(cur), std::placeholders::_1);
			bool finish(false);
			std::vector<int> feasible(mvv_neighbors[interpreter.curPositionIdx(pro, cur)]);
			interpreter.fillterFeasible(pro, cur, feasible);
			if (!feasible.empty()) {
				auto next_step_iter = (rnd->uniform.greedyRandom(feasible.begin(), feasible.end(), nei_weight_fun));
				finish = interpreter.stepNext(pro, cur, *next_step_iter);
			}
			if (!finish) {
				interpreter.calNextFeasible(pro, cur);
				auto next_step_iter = (rnd->uniform.spinWheelIdx(cur.feasibleNeighbors().begin(), cur.feasibleNeighbors().end(), nei_weight_fun));
				if (next_step_iter == cur.feasibleNeighbors().end()) return false;
				return interpreter.stepNext(pro, cur, *next_step_iter);;
			}
			return finish;

		}
	};
	template<typename TInterpreter>
	inline void MMAS<TInterpreter>::pheroSmoothing() {
		if (!(m_iteration % 100))
		{
			m_branchFactor = nodeBranching();
			if (m_branchFactor < m_branchFrc && ((m_iteration - m_restartFoundBest) > 250))
			{
				m_restartBest.clear();
				for (auto &it : mvv_phero) {
					std::fill(it.begin(), it.end(), m_pheroMax);
				}

			}
		}
	}
	template<typename TInterpreter>
	inline void MMAS<TInterpreter>::updatePheromonePolicy(const std::vector<std::unique_ptr<SolutionType>> &update_sols) {
		for (auto &it : m_best_sols_edges) {
			std::fill(it.begin(), it.end(), false);
		}
		for (auto &best_iter : update_sols) {
			for (auto &edge : best_iter->edges()) {
				m_best_sols_edges[edge.first][edge.second] = true;
			}
		}

		for (int dim1(0); dim1 < mvv_phero.size(); ++dim1) {
			for (int dim2(0); dim2 < mvv_phero[dim1].size(); ++dim2) {
				if (m_best_sols_edges[dim1][dim2])
					mvv_phero[dim1][dim2] = m_rho * mvv_phero[dim1][dim2] + update_sols.front()->fitness();
				else {
					mvv_phero[dim1][dim2] *= m_rho;
				}
			}
		}
	}

	template<typename TInterpreter>
	inline bool MMAS<TInterpreter>::updateSolutions(Problem *pro, SolutionType &cur_indiv, std::vector<std::unique_ptr<SolutionType>> &curPops) {
		bool insertSol(false);
		{
			Real fitness(m_fitness_update(pro, cur_indiv));
			if (curPops.empty()) {
				insertSol = true;
			}
			else if (curPops.front()->fitness() < fitness) {
				curPops.clear();
				insertSol = true;
			}
			else if (curPops.front()->fitness() == fitness) {
				insertSol = true;
			}
			if (insertSol) {
				curPops.emplace_back(new SolutionType(cur_indiv));
				curPops.back()->setFitness(fitness);
			}
			insertSol = false;
		}
		return insertSol;
	}

	template<typename TInterpreter>
	inline void MMAS<TInterpreter>::checkPheromoneTrailLimits() {
		Real curPheroMin(0);
		for (int i = 0; i < mvv_phero.size(); i++)
		{
			curPheroMin = m_pheroMin / (mvv_neighbors[i].size() - 1);
			for (int j(0); j < mvv_neighbors[i].size(); ++j) {
				if (mvv_phero[i][j] < curPheroMin)
				{
					mvv_phero[i][j] = curPheroMin;
				}
				else if (mvv_phero[i][j] > m_pheroMax) {
					mvv_phero[i][j] = m_pheroMax;
				}
			}
		}
	}


	template<typename TInterpreter>
	inline void MMAS<TInterpreter>::updatePheroMinAndMax() {
		int n = mvv_phero.size();
		Real p_x = exp(log(0.05) / n);
		m_pheroMin = 1. * (1. - p_x) / (p_x / 2);
		m_pheroMax = 1. / (1.0 - m_rho) * m_globalBest.front()->fitness();
		m_pheroMin = m_pheroMax * m_pheroMin;
	}
	template<typename TInterpreter>
	inline Real MMAS<TInterpreter>::nodeBranching() {
		int n = mvv_phero.size();
		Real minTemp, maxTemp, cutoff, avg = 0.;
		int num_branches(0);

		for (int i(0); i < n; ++i) {
			maxTemp = minTemp = mvv_phero[i][mvv_neighbors[i].front()];

			for (auto &iter : mvv_neighbors[i]) {
				maxTemp = std::max(mvv_phero[i][iter], maxTemp);
				minTemp = std::min(mvv_phero[i][iter], minTemp);
			}
			cutoff = minTemp + m_lambda * (maxTemp - minTemp);
			num_branches = 0;

			for (auto &iter : mvv_neighbors[i]) {
				if (mvv_phero[i][iter] > cutoff) {
					num_branches++;
				}
			}
			avg += Real(num_branches) / static_cast<Real>(mvv_neighbors[i].size());
		}
		return avg / 2.0;
	}
}

#endif