#include "sp_operator.h"
#include"../../../../problem/combination/selection_problem/selection_problem.h"

void ofec::SelectionProblemOperator::initialize(Problem *pro, Algorithm *alg) 
{
	
	auto & par(GET_PARAM(alg.idParam()));
	getTypeVal<int, EFitnessCalType>(par, "combinatorial fitness type", m_fitness_type,EFitnessCalType::Objective);
	getTypeVal<int, EMemoryType>(par, "combinatorial memory type", m_memory_type, EMemoryType::Change);
	//m_fitness_type = EFitnessCalType::HLS;
	//getTypeVal<bool,bool>(par,"",m_noisy_flag,)
	
	
	//if (m_fitness_type == EFitnessCalType::HLS) {
	//	cout << "HLS" << std::endl;
	//}
	//else cout << "Single" << std::endl;
}
int ofec::SelectionProblemOperator::learn_from_local(const SolutionType& cur, Random *rnd, Problem *pro, const InterpreterType& interpreter) const
{
	if (interpreter.stepFinish(pro, cur)) {
		return -1;
	}
	else {
		return cur.variable()[cur.getLoc()];
	}
}

int ofec::SelectionProblemOperator::learn_from_global(const SolutionType& cur, Random *rnd, Problem *pro, const InterpreterType& interpreter, const std::function<Real(const SolutionType& cur, int to)>& proFun) const
{
	if (interpreter.stepFinish(pro, cur)) {
		return -1;
	}	
	int matSize=interpreter.getMatrixSize()[interpreter.curPositionIdx(pro,cur)];
	std::function<Real(int to)> nei_weight_fun = [&](int to) {
		return proFun(cur, to);
	};
	auto next_step_iter(rnd->uniform.spinWheel<int,Real>(0, matSize, nei_weight_fun));
	if (next_step_iter == matSize)next_step_iter = -1;
	return next_step_iter;
}

int ofec::SelectionProblemOperator::learn_from_other(const SolutionType& cur, Random *rnd, Problem *pro, const InterpreterType& interpreter, const SolutionType& other) const
{
	if (interpreter.stepFinish(pro, cur)) {
		return -1;
	}
	else {
		return other.variable()[cur.getLoc()];
	}
	
}


bool ofec::SelectionProblemOperator::learn_from_other(SolutionType& cur, Random *rnd, Problem *pro, const InterpreterType& interpreter, const SolutionType& other, Real radius) const
{
	Real dis(cur.variableDistance(other, pro));
	
	if (dis > radius) {
		cur.reset();
		int numVar = std::ceil((dis - radius) * cur.variable().size());
		std::vector<int> idxs(cur.variable().size());
		for (int idx(0); idx < cur.variable().size();++idx) {
			idxs[idx] = idx;
		}
		rnd->uniform.shuffle(idxs.begin(), idxs.end());
		for (auto& dim : idxs) {
			if (cur.variable()[dim] != other.variable()[dim]) {
				cur.variable()[dim] = other.variable()[dim];
				if (dim + 1 < cur.variable().size()) {
					cur.variable()[dim+1] = other.variable()[dim+1];
					--numVar;
					if (numVar == 0)break;
				}

			}
		}
		return true;
	}
	return false;
}

void ofec::SelectionProblemOperator::localSearch(SolutionType& cur, Random *rnd, Problem *pro, Algorithm *alg, int totalEvals, int curType)
{

//	std::cout << "localSearch" << std::endl;
	thread_local std::vector<std::pair<int, int>> neis;
	SolutionType neiSolution(cur);
	neiSolution.reset();
	neis.clear();
	auto& sp_object(GET_DYN_SP(pro));
	for (int dimIdx(0); dimIdx < sp_object.numberVariables(); ++dimIdx) {
		for (int posIdx(0); posIdx < sp_object.numPoints(dimIdx); ++posIdx) {
			neis.push_back(std::make_pair(dimIdx, posIdx));
		}
	}
	rnd->uniform.shuffle(neis.begin(), neis.end());
	int neiIdx(0);
	while (totalEvals--) {

		auto& nei(neis[neiIdx]);
		neiSolution.variable()[nei.first] = nei.second;
		neiSolution.reset();
		evaluate(pro, alg, rnd, neiSolution, true);
		if (better(pro, neiSolution, cur)) {
			cur = neiSolution;
		}
		else neiSolution = cur;
		while (neiSolution.variable()[neis[neiIdx].first] != neis[neiIdx].second) {
			++neiIdx;
			if (neiIdx == neis.size()) {
				rnd->uniform.shuffle(neis.begin(), neis.end());
				neiIdx = 0;
			}
		}
	}	
}

void ofec::SelectionProblemOperator::localSearch(SolutionType& cur,
	int totalEvals, int curType,
	const std::function<void(SolutionType& cur)>& fitnessCal,
	Random *rnd, Problem *pro, Algorithm *alg
) {
	//	std::cout << "localSearch" << std::endl;
	thread_local std::vector<std::pair<int, int>> neis;
	SolutionType neiSolution(cur);
    neiSolution.reset();
	neis.clear();
	auto& sp_object(GET_DYN_SP(pro));
	for (int dimIdx(0); dimIdx < sp_object.numberVariables(); ++dimIdx) {
		for (int posIdx(0); posIdx < sp_object.numPoints(dimIdx); ++posIdx) {
			neis.push_back(std::make_pair(dimIdx, posIdx));
		}
	}
	rnd->uniform.shuffle(neis.begin(), neis.end());
	int neiIdx(0);
	while (totalEvals--) {
		auto& nei(neis[neiIdx]);
		neiSolution.reset();
		neiSolution.variable()[nei.first] = nei.second;
		fitnessCal(neiSolution);
		if (neiSolution.fitness() > cur.fitness()) {
			cur = neiSolution;
		}
		else {
			neiSolution = cur;
		}
		while (neiSolution.variable()[neis[neiIdx].first] != neis[neiIdx].second) {
			++neiIdx;
			if (neiIdx == neis.size()) {
				rnd->uniform.shuffle(neis.begin(), neis.end());
				neiIdx = 0;
			}
		}
	}

}


int ofec::SelectionProblemOperator::localSearch(SolutionType& cur,
	int total_iters,
	Problem *pro, Algorithm *alg, Random *rnd,
	std::shared_ptr<EvaluationStrategyBase<SolutionType>>& eval_stra,
	int pop_idx 
) {
	int rf(0);
	//	std::cout << "localSearch" << std::endl;
	thread_local std::vector<std::pair<int, int>> neis;
	SolutionType neiSolution(cur);
	neiSolution.reset();
	neis.clear();
	auto& sp_object(GET_DYN_SP(pro));
	for (int dimIdx(0); dimIdx < sp_object.numberVariables(); ++dimIdx) {
		for (int posIdx(0); posIdx < sp_object.numPoints(dimIdx); ++posIdx) {
			neis.push_back(std::make_pair(dimIdx, posIdx));
		}
	}
	rnd->uniform.shuffle(neis.begin(), neis.end());
	int neiIdx(0);
	while (total_iters--) {
		auto& nei(neis[neiIdx]);
		neiSolution.reset();
		neiSolution.variable()[nei.first] = nei.second;

		rf|=eval_stra->calFitness(neiSolution, pop_idx);
		if (neiSolution.fitness() > cur.fitness()) {
			cur = neiSolution;
		}
		else {
			neiSolution = cur;
		}
		while (neiSolution.variable()[neis[neiIdx].first] != neis[neiIdx].second) {
			++neiIdx;
			if (neiIdx == neis.size()) {
				rnd->uniform.shuffle(neis.begin(), neis.end());
				neiIdx = 0;
			}
		}
	}
	return rf;
}


int ofec::SelectionProblemOperator::evaluate(Problem *pro, Algorithm *alg, Random *rnd, SolutionType& curSol, bool effective) const
{
	std::vector<Real> origin(curSol.objective());
	
	int rf = 0;
	
	switch (m_fitness_type) {
	case EFitnessCalType::Objective: {
		int rf = curSol.evaluate(pro, alg, effective);
		break;
	}
	case EFitnessCalType::HLS: {

		std::vector< std::unique_ptr<SolutionBase>> samples(m_sample_num);
		GET_DYN_SP(pro)->generate_HLS_samples<SP_Solution>(curSol, samples, rnd);
		//generateSamplesSolutions(pro, rnd, curSol, samples);
		Real fitness(0);
		for (auto& it : samples) {
			rf |= it->evaluate(pro, alg, effective);
			fitness += it->objective()[0];
		}
		rf |= curSol.evaluate(pro, alg, effective);
		fitness += curSol.objective()[0];
		fitness /= (m_sample_num + 1.0);
		curSol.objective()[0] = fitness;

		break;
	}
	}
	if (m_memory_type == EMemoryType::Decrease) {
		int curTimes = ++curSol.evaluateTimes();
		curSol.objective()[0] = (curSol.objective()[0] + origin.front() * (curTimes - 1)) / Real(curTimes);

	}
	return rf;
}




