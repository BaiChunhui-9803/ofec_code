#include "sars_algorithm.h"
#include "../../../../problem/combination/one_max/one_max.h"


void ofec::SARS_OneMax::evolve() {

	if (step == 1) {
		innerBudget += innerBudget;
		innerBudgetParam = (double)innerBudget;
		++stepMain;
		Tend = modularGA::sars::T_from_DeltaE_and_P(1.0, 1.0 / sqrt(innerBudget));
		epsilon = modularGA::sars::epsilon_from_T_and_step(Tstart, Tend, innerBudget);
		m_new_sol->initialize(m_problem.get(), m_random.get());
		m_new_sol->evaluate(m_problem.get(), this);
		ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
		m_new_sol->setFitness(pos * m_new_sol->objective()[0]);
		m_cur_sol->copy(*m_new_sol);


		if (m_new_sol->fitness() > m_best_sol->fitness()) {
			m_best_sol->copy(*m_new_sol);
		}
	}


	bool unchanged = true;
	auto& solx = dynamic_cast<ofec::OneMax::solutionType&>(*m_new_sol);
	// until the solution changes, repeat
	do
	{
		// flip each bit with the independent probability of 1/n
		for (int i = m_nvar; (--i) >= 0;)
		{
			if (m_random->uniform.next() < p)
			{
				unchanged = false; // there was a change

				solx.variable()[i] ^= 1;      // flip the bit
			}
		}
	} while (unchanged); // repeat until at least one change

	// evaluate the new candidate solution
	m_new_sol->evaluate(m_problem.get(), this);


	ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	m_new_sol->setFitness(pos * m_new_sol->objective()[0]);

	if (m_new_sol->fitness() > m_best_sol->fitness()) {
		*m_best_sol = *m_new_sol;
	}

	// if new solution is at least as good as current one, accept it
	// otherwise: if check if it is acceptable at the current temperature
	if ((m_new_sol->fitness() >= m_cur_sol->fitness()) ||
		(m_random->uniform.next() < modularGA::sars::p_accept(m_cur_sol->fitness() - m_new_sol->fitness(), modularGA::sars::temperature(Tstart, epsilon, step))))
	{
		*m_cur_sol = *m_new_sol;
	}


	std::cout << "best fitness\t" << m_best_sol->fitness() << std::endl;

	step = (step + 1) % innerBudget + 1;
}
void ofec::SARS_OneMax::run_()
{
	while (!terminating()) {
		evolve();
#ifdef OFEC_DEMO
//		updateBuffer();
#endif
	}
	
}

void ofec::SARS_OneMax::initialize_()
{
	Algorithm::initialize_();
	m_nvar = m_problem->numVariables();
	Tend = 0;
	epsilon = 0;
	innerBudgetParam = 0;
	step = 1;
	innerBudget = 512;
	stepMain = 1;
	p = 1.0 / ((double)m_nvar);
	Tstart = modularGA::sars::T_from_DeltaE_and_P(std::max(1.0, m_nvar / 4.0), 0.1);





	m_cur_sol.reset(m_problem->createSolution());
	m_best_sol.reset(m_problem->createSolution());
	m_new_sol.reset(m_problem->createSolution());

	//m_best_sol->objective()[0]= 


	m_best_sol->initialize(m_problem.get(), m_random.get());
	m_best_sol->evaluate(m_problem.get(), this);
	ofec::Real pos = (m_problem->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
	m_best_sol->setFitness(pos * m_best_sol->objective()[0]);

	//m_best_sol->copy(m_best_sol);
	
}
