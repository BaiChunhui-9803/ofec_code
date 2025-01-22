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
#ifndef OFEC_ANT_SYSTEM_H
#define OFEC_ANT_SYSTEM_H

#include "../../../template/classic/ant_colony_optimization/aco_model.h"
#include <cmath>

namespace ofec {
	template<typename TInterpreter>
	class AS : public ModelACO<TInterpreter> {
	public:
		using SolutionType = typename TInterpreter::SolutionType;
		using IntrepreterType = TInterpreter;
	protected:
		Real m_alpha = 1.0;
		Real m_beta = 5.0;
		Real m_Q = 100.;
		Real m_rho = 0.5;
		Real getPro(const TInterpreter &interpreter, int pro_id, const SolutionType &cur, int to) override {
			return  pow(mvv_phero[interpreter.curPositionIdx(pro_id, cur)][to], m_alpha) * pow(interpreter.heursticInfo(pro_id, cur, to), m_beta);
		}

	public:
		void initialize(int id_param, Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter &interpreter) override {
			ModelACO::initialize(id_param, rnd, pro, alg, interpreter);
			auto &v = GET_PARAM(id_param);
			if (v.has("AS_alpha")) {
				m_alpha = v.get<Real>("AS_alpha");
			}
			else m_alpha = 1.0;
			if (v.has("AS_beta")) {
				m_beta = v.get<Real>("AS_beta");
			}
			else m_beta = 5.0;
			if (v.has("AS_Q")) {
				m_Q = v.get<Real>("AS_Q");
			}
			else m_Q = 100.;
			if (v.has("AS_rho")) {
				m_rho = v.get<Real>("AS_rho");
			}
			else m_rho = 0.5;
		}

		void initializePheromone() override {
			for (auto &it : mvv_phero) {
				for (auto &it2 : it) {
					it2 = 1.0;
				}
			}
		}

		void globalUpdatePheromone(Random *rnd, Problem *pro, Algorithm *alg, const TInterpreter &interpreter, std::vector<std::unique_ptr<SolutionType>> &pops) override {
			for (auto &it : mvv_phero) {
				for (auto &it2 : it) {
					it2 *= m_rho;
				}
			}
			Real fitness(0);
			for (int idx(0); idx < pops.size(); ++idx) {
				fitness = m_fitness_update(pro, *pops[idx]);
				interpreter.updateEdges(pro, *pops[idx]);
				for (auto &cur_edge : pops[idx]->edges()) {
					mvv_phero[cur_edge.first][cur_edge.second] += m_Q * fitness;
				}
			}
		}
	};
}

#endif // ! OFEC_ANT_SYSTEM_H

