#include"sp_interpreter.h"
#include "../../../../../instance/problem/combination/selection_problem/selection_problem.h"
using namespace ofec;
using namespace sp;

void ofec::SP_Interpreter::initializeByProblem(Problem *pro)
{
	m_start_state_idx = 0;
	auto& sp(GET_DYN_SP(pro));
	m_dim_size = sp.numberVariables();
	m_cur2idx.resize(m_dim_size);
	m_pos_size = 1;
	for (int dim(0); dim < m_dim_size; ++dim) {
		m_cur2idx[dim].resize(sp.numPoints(dim));
	}
	for (int idx(0); idx+1 < m_dim_size; ++idx) {
		m_pos_size += sp.numPoints(idx);
	}

	m_matrix_size.resize(m_pos_size);
	m_idx2cur.resize(m_pos_size);
	
	int idx(0);
	m_matrix_size[idx] = sp.numPoints(idx);
	m_idx2cur[idx] = std::make_pair(-1, -1);
	++idx;
	for (int x(0); x+1 < m_dim_size; ++x) {
		for (int y(0); y < sp.numPoints(x);  ++y) {
			m_idx2cur[idx] = std::make_pair(x, y);
			m_cur2idx[x][y] = idx;
			m_matrix_size[idx] = sp.numPoints(x+1);
			++idx;
		}
	}
	for (auto& it : m_cur2idx.back()) it = -1;
	
}

ofec::Real ofec::SP_Interpreter::heursticInfo(Problem *pro, int from, int to, int obj_idx) const
{
	return 1.0/(GET_DYN_SP(pro)->getHeuristicInfo(m_idx2cur[from].first, m_idx2cur[from].second, to)+m_eps);
}

ofec::Real ofec::SP_Interpreter::heursticInfo(Problem *pro, const SolutionType& indi, int to, int obj_idx) const
{
	Real heu_info(m_eps);
	int cur_pos(curPositionIdx(pro,indi));
	if (cur_pos >= 0) {
		heu_info += heursticInfo(pro,indi.getLoc(),to);
	}
	return heu_info;
}

ofec::SP_Interpreter::SolutionType ofec::SP_Interpreter::heuristicSol(Problem *pro, Algorithm *alg, Random *rnd) const
{
	SolutionType sol;
	stepInit(pro, sol);
	while (!stepFinish(pro,sol)) {                                                                          (-1);
		Real maxH(-1e9), curH(0);
		std::vector<int> nei;
		int curPos(curPositionIdx(pro,sol));

		for (int to(0); to < m_matrix_size[curPos]; ++to) {
			curH = heursticInfo(pro, sol, to);
			if (maxH < curH) {
				maxH = curH;
				nei.clear();
				nei.push_back(maxH);
			}
			else if (maxH == curH) {
				nei.push_back(to);
			}
		}
		if (nei.empty()) {
			break;
		}
		stepNext(pro, sol, nei[rnd->uniform.nextNonStd<int>(0,nei.size())]);
	}
	return std::move(sol);
}

int ofec::SP_Interpreter::curPositionIdx(Problem *pro, const SolutionType& sol) const
{
	if (sol.getLoc() == 0)return m_start_state_idx;
	else return m_cur2idx[sol.getLoc()-1][sol.variable()[sol.getLoc()-1]];
}

void ofec::SP_Interpreter::stepInit(Problem *pro, SolutionType& sol) const
{
	SequenceInterpreter::stepInit(pro, sol);
	sol.getLoc() = 0;
	sol.variable().resize(m_dim_size);
	sol.flagCurrentEdgesUpdated() = false;
	sol.flagEdgesUpdated() = false;
	sol.flagFeasibleNeighborsUpdated() = false;
	
}

void ofec::SP_Interpreter::calNextFeasible(Problem *pro, SolutionType& sol) const
{
	if (!sol.flagFeasibleNeighborsUpdated()) {
		int cur_idx(curPositionIdx(pro, sol));
		if (cur_idx >= 0) {
			sol.feasibleNeighbors().resize(m_matrix_size[cur_idx]);
			for (int idx(0); idx < sol.feasibleNeighbors().size(); ++idx) {
				sol.feasibleNeighbors()[idx] = idx;
			}
		}
		else {
			sol.feasibleNeighbors().clear();
		}
		sol.flagFeasibleNeighborsUpdated() = true;
	}

}

bool ofec::SP_Interpreter::stepFeasible(Problem *pro, SolutionType& sol) const
{
	return sol.getLoc() != m_dim_size;
}

void ofec::SP_Interpreter::stepBack(Problem *pro, SolutionType& sol) const
{
	if (sol.getLoc()) {
		--sol.getLoc();
		sol.flagCurrentEdgesUpdated() = false;
		sol.flagEdgesUpdated() = false;
		sol.flagFeasibleNeighborsUpdated() = false;
	}
}

bool ofec::SP_Interpreter::stepNext(Problem *pro, SolutionType& sol, int next) const
{
	if (sol.getLoc() < m_dim_size) {
		sol.variable()[sol.getLoc()++] = next;
		sol.flagCurrentEdgesUpdated() = false;
		sol.flagEdgesUpdated() = false;
		sol.flagFeasibleNeighborsUpdated() = false;
		return true;
	}
	return false;
}

void ofec::SP_Interpreter::updateCurEdges(Problem *pro, SolutionType& sol) const
{
	if (!sol.flagCurrentEdgesUpdated()) {
		int from(m_start_state_idx);
		if (sol.getLoc()) from = m_cur2idx[sol.getLoc() - 1][sol.variable()[sol.getLoc() - 1]];
		int to = sol.variable()[sol.getLoc()];
		sol.currentEdges().resize(1);
		sol.currentEdges().front() = std::make_pair(from, to);
		sol.flagCurrentEdgesUpdated() = true;
	}
}

void ofec::SP_Interpreter::updateEdges(Problem *pro, SolutionType& sol, bool all ) const
{
	if (!sol.flagEdgesUpdated()) {
		int loc = sol.getLoc();
		if (all) {
			loc = sol.variable().size();
		}

		int from(m_start_state_idx), to(-1);
		sol.edges().resize(loc);
		int idx(0);
		while (true) {
			sol.edges()[idx] = std::make_pair(from, sol.variable()[idx]);
			from = m_cur2idx[idx][sol.variable()[idx]];
			if (++idx == loc)break;
		}
		sol.flagEdgesUpdated() = true;
	}
}

bool ofec::SP_Interpreter::stepFinish(Problem *pro,const SolutionType& sol) const
{
	return sol.getLoc() == m_dim_size;
}

void ofec::SP_Interpreter::stepFinal(Problem *pro, SolutionType& sol) const
{
	sol.getLoc() = m_dim_size;
	sol.flagCurrentEdgesUpdated() = false;
	sol.flagEdgesUpdated() = false;
	sol.flagFeasibleNeighborsUpdated() = false;
}

void ofec::SP_Interpreter::getNearestNeighbors(Problem *pro, int len, int from_idx, std::vector<int>& neighbors) const
{
	std::vector<std::pair<int, double>> idx2dis(m_matrix_size[from_idx]);
	for (int to(0); to < idx2dis.size(); ++to) {
		idx2dis[to].first = to;
		idx2dis[to].second = heursticInfo(pro, from_idx, to);
	}
	std::sort(idx2dis.begin(), idx2dis.end(), [&](const std::pair<int,double>& a,const std::pair<int,double>& b) {
		return a.second > b.second;
	});
	len = std::min<int>(len, idx2dis.size());
	neighbors.resize(len);
	for (int idx(0); idx < len; ++idx) {
		neighbors[idx] = idx2dis[idx].first;
	}
}

void ofec::SP_Interpreter::fillterFeasible(Problem *pro, SolutionType& sol, std::vector<int>& feasible) const
{
	int cur_idx(curPositionIdx(pro, sol));
	if (cur_idx >= 0) {
		feasible.resize(m_matrix_size[cur_idx]);
		for (int idx(0); idx < feasible.size(); ++idx) {
			feasible[idx] = idx;
		}
	}
	else feasible.clear();
}


