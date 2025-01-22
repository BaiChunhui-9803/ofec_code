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

#ifndef OFEC_PACO_H
#define OFEC_PACO_H

#include "../../../template/classic/ant_colony_optimization/aco_model.h"
#include "../../../../../utility/functional.h"
#include <cmath>

#include"../../../template/classic/ant_colony_optimization/aco_model.h"
#include"../../../../../utility/functional.h"
#include<cmath>

namespace ofec {
	template<typename 
		TInterpreter
	>
	class PACO : public ACO_Model<TInterpreter> {
	public:
		enum class ReplacementStrategy { Age, Quality, Prob, Prob_Age};
	public:
		using SolutionType = typename TInterpreter::SolutionType;
		using IntrepreterType = TInterpreter;
		//using SolutionIter = typename std::vector<std::unique_ptr<SolutionType>>::iterator;
	protected:

	protected:
		ReplacementStrategy m_replace_strategy;
		bool m_elite_update = false;
		long long m_iteration = 0;
		Real m_alpha = 1.0;
		Real m_beta = 5.0;
		Real m_pheroMax = 1.0;
		Real m_pheroInit;
		Real m_pheroDelta = 1.0;
		int m_k = 3;
		int m_del_ant_idx = -1;
		std::vector<std::unique_ptr<SolutionType>> m_memory;
		std::unique_ptr<SolutionType> m_elite;
	protected:
		virtual Real getPro(const TInterpreter& interpreter, Problem *pro, const SolutionType& cur, int to) {
			return  pow(mvv_phero[interpreter.curPositionIdx(pro, cur)][to], m_alpha) * pow(interpreter.heursticInfo(pro, cur, to), m_beta);
		}

	private:
		std::vector<int> m_idxs;
		std::vector<Real> m_pro;

	public:
		virtual void initialize(int id_param, Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter& interpreter) override {
			ACO_Model::initialize(id_param, rnd, pro, alg, interpreter);
			auto& v = GET_PARAM(id_param);
			m_replace_strategy = ReplacementStrategy::Age;
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
			m_pheroInit = m_pheroMax / static_cast<Real>(mat.size());
			if (m_elite_update)
				m_pheroDelta = (m_pheroMax - m_pheroInit) / static_cast<Real>(m_k + 1);
			else
				m_pheroDelta = (m_pheroMax - m_pheroInit) / static_cast<Real>(m_k);
			mvv_phero.resize(mat.size());
			for (int idx(0); idx < mat.size(); ++idx) {
				mvv_phero[idx].resize(mat[idx]);
				std::fill(mvv_phero[idx].begin(), mvv_phero[idx].end(), m_pheroInit);
			}
			m_memory.clear();
			m_elite = nullptr;
			m_del_ant_idx = -1;
			m_idxs.resize(m_k + 1);
			for (int idx(0); idx < m_idxs.size(); ++idx) {
				m_idxs[idx] = idx;
			}
			m_pro.resize(m_k + 1);

		}

		virtual void globalUpdatePheromone(Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter& interpreter, std::vector<std::unique_ptr<SolutionType>>& pops) {
			++m_iteration;
			if (pops.empty())return;
			for (auto &it : pops)  it->setFitness(m_fitness_update(pro, *it));
			//std::function<Real(const std::vector<std::unique_ptr<SolutionType>>::iterator& cur_iter)> pro_fun = [](const std::vector<std::unique_ptr<SolutionType>>::iterator& iter) {
			//	return (*iter)->fitness();
			//	//return iter->fitness();
			//};
			//rnd->uniform.greedyRandom<std::vector<std::unique_ptr<SolutionType>>::iterator>(pops.begin(), pops.end(),pro_fun);
			auto sol_iter = rnd->uniform.greedyRandom<decltype(pops.begin()), Real>(pops.begin(), pops.end(),
				[](const decltype(pops.begin()) &iter) { return (*iter)->fitness(); });
			if (sol_iter != pops.end())return;
			else {
				if (m_memory.size() >= m_k) {
					switch (m_replace_strategy) {
					case ReplacementStrategy::Age: {
						m_del_ant_idx = (m_del_ant_idx + 1) % m_k;
						break;
					}
					case ReplacementStrategy::Prob: {
						Real min_pro(std::numeric_limits<Real>::max()), max_pro(-min_pro);
						for (auto &it : m_memory) {
							max_pro = std::max(max_pro, it->fitness());
							min_pro = std::min(min_pro, it->fitness());
						}
						for (int idx(0); idx < m_k; ++idx) {
							m_pro[idx] = 1.0 - mapReal<Real>(m_memory[idx]->fitness(), min_pro, max_pro, 0, 1);
						}
						m_pro.back() = 1.0 - mapReal<Real>((*sol_iter)->fitness(), min_pro, max_pro, 0., 1.0);
						m_del_ant_idx = rnd->uniform.spinWheel(m_pro.begin(), m_pro.end()) - m_pro.begin();
						if (m_del_ant_idx == m_k) m_del_ant_idx = -1;
						break;
					}
					case ReplacementStrategy::Prob_Age: {
						Real min_pro(std::numeric_limits<Real>::max()), max_pro(-min_pro);
						for (auto &it : m_memory) {
							max_pro = std::max(max_pro, it->fitness());
							min_pro = std::min(min_pro, it->fitness());
						}
						for (int idx(0); idx < m_k; ++idx) {
							m_pro[idx] = 1.0 - mapReal<Real>(m_memory[idx]->fitness(), min_pro, max_pro, 0, 1);
						}
						m_del_ant_idx = rnd->uniform.spinWheel(m_pro.begin(), m_pro.begin() + m_k) - m_pro.begin();
						break;
					}
					case ReplacementStrategy::Quality: {
						m_del_ant_idx = rnd->uniform.spinWheel<decltype(m_idxs.begin()), Real>(m_idxs.begin(), m_idxs.begin() + m_k,
							[&](const decltype(m_idxs.begin()) &iter) {return double(-pops[*iter]->fitness()); }) - m_idxs.begin();
						if ((*sol_iter)->fitness() < pops[m_del_ant_idx]->fitness()
							|| ((*sol_iter)->fitness() == pops[m_del_ant_idx]->fitness()
								&& rnd->uniform.next() < 0.5)) {
							m_del_ant_idx = -1;
						}
						break;
					}
					}
				}
				if (m_memory.size() < m_k) {
					m_memory.emplace_back(new SolutionType(**sol_iter));
					interpreter.updateEdges(pro, *m_memory.back());
					for (auto &it : m_memory.back()->edges()) {
						mvv_phero[it.first][it.second] += m_pheroDelta;
					}
				}
				else {
					if (m_del_ant_idx != -1) {
						interpreter.updateEdges(pro, *(*sol_iter));
						for (auto &it : m_memory[m_del_ant_idx]->edges()) {
							mvv_phero[it.first][it.second] -= m_pheroDelta;
						}
						m_memory[m_del_ant_idx].reset(new SolutionType(**sol_iter));
						for (auto &it : m_memory[m_del_ant_idx]->edges()) {
							mvv_phero[it.first][it.second] += m_pheroDelta;
						}
					}
				}
				if (m_elite_update) {
					if (m_elite == nullptr || (*sol_iter)->fitness() > m_elite->fitness()
						|| ((*sol_iter)->fitness() == m_elite->fitness()
							&& rnd->uniform.next() < 0.5)) {
						if (m_elite != nullptr) {
							for (auto &it : m_elite->edges()) {
								mvv_phero[it.first][it.second] -= m_pheroDelta;
							}
						}
						m_elite.reset(new SolutionType(**sol_iter));
						for (auto &it : m_elite->edges()) {
							mvv_phero[it.first][it.second] += m_pheroDelta;
						}
					}
				}
			}
		}
		virtual bool stepNext(Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter& interpreter, SolutionType& cur) {

		bool stepNext(Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter &interpreter, SolutionType &cur) override {
			std::function<Real(int to)> nei_weight_fun = std::bind(&PACO::getPro, this, std::cref(interpreter), pro, std::cref(cur), std::placeholders::_1);
			{
				interpreter.calNextFeasible(pro, cur);
				auto next_step_iter = (rnd->uniform.spinWheelIdx(cur.feasibleNeighbors().begin(), cur.feasibleNeighbors().end(), nei_weight_fun));
				if (next_step_iter == cur.feasibleNeighbors().end()) return false;
				return interpreter.stepNext(pro, cur, *next_step_iter);;
			}
		}
	};

}

#endif