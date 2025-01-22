#include "ghc_onemax.h"
#include "../../../../problem/combination/one_max/one_max.h"


void ofec::GHC_OneMax::evolve() {
	
	m_cur_sol->copy(*m_best_sol);
	//auto it = dynamic_cast<decltype(*m_cur_sol)&>(*m_cur_sol);


	auto& curSolX = dynamic_cast<ofec::OneMax::solutionType&>(*m_cur_sol);
	// Filp the bit at flip_index.
	curSolX.variable()[m_flit_idx] = (curSolX.variable()[m_flit_idx] + 1) % 2;
	if (++m_flit_idx >= m_dim)
	{
		m_flit_idx = 0;
	}

	m_cur_sol->evaluate(m_problem.get(), this);
	//double pos = std
	ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	m_cur_sol->setFitness(pos * m_cur_sol->objective()[0]);


	if (m_cur_sol->fitness() >= m_best_sol->fitness())
	{
		m_best_sol->copy(*m_cur_sol);
	}
}

void ofec::GHC_OneMax::run_()
{
	while (!terminating()) {
		evolve();
#ifdef OFEC_DEMO
	//	updateBuffer();
#endif
	}
}

void ofec::GHC_OneMax::initialize_()
{
	Algorithm::initialize_();
	m_flit_idx = 0;
	m_dim = m_problem->numVariables();
	m_cur_sol.reset(m_problem->createSolution());
	m_best_sol.reset(m_problem->createSolution());


	m_cur_sol->initialize(m_problem.get(), m_random.get());
	m_cur_sol->evaluate(m_problem.get(), this);
	ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	m_cur_sol->setFitness(pos * m_cur_sol->objective()[0]);

	m_best_sol->copy(*m_cur_sol);
	

}
