#include "iterated_local_search.h"
#include "../../../../../../core/problem/continuous/continuous.h"

void ofec::IteratedLocalSearch::run_()
{
    using namespace ofec;


    static int id = 0;
    ProFunction::vector_t solution_x(m_problem->numberVariables());
    int dim = m_problem->numberVariables();
    for (int idx(0); idx < dim; ++idx) {
        solution_x(idx) = m_local_best.variable()[idx];
    }
    while (!terminating()) {
  
        auto [solutionBest, solver_state] = m_proSolver->Minimize(m_proFun, solution_x);
        
        if (solutionBest.value* m_minimize_flag < m_local_best.objective()[0] ) {
            m_local_best.objective()[0] = solutionBest.value * m_minimize_flag;
       //     std::cout << "cur improve id\t" << ++id <<"curOjb\t"<< solutionBest.value<< std::endl;
    
            for (int idx(0); idx < dim; ++idx) {
                m_local_best.variable()[idx] = solutionBest.x(idx);
            }
            solution_x = solutionBest.x;
        }
        solution_x = ProFunction::BoundedPerturbation(m_random.get(), solution_x, m_perturbation_range, m_domain);
        
        //m_pop.evolve(m_problem.get(), this, m_random.get());
    }
}


void ofec::IteratedLocalSearch::initialize_()
{
	using namespace ofec;
	Algorithm::initialize_();


    auto curPro = CAST_CONOP(m_problem.get());
	m_domain = curPro->domain();
	m_perturbation_range.resize(m_domain.size());
	for (int idx(0); idx < m_perturbation_range.size(); ++idx) {
		m_perturbation_range[idx]  = m_perturbation_ratio  * (m_domain[idx].limit.second - m_domain[idx].limit.first);
	}

    if (m_problem->optimizeMode()[0] == OptimizeMode::kMaximize) {
        m_minimize_flag = -1.0;
    }
    else m_minimize_flag = 1.0;
	m_proFun.initialize(m_problem.get(), this);


    {
        auto& domain(m_domain);
        ProFunction::vector_t lower_bound(m_problem->numberVariables());
        ProFunction::vector_t upper_bound(m_problem->numberVariables());
        for (int idx(0); idx < m_problem->numberVariables(); ++idx) {
            lower_bound(idx) = domain[idx].limit.first;
            upper_bound(idx) = domain[idx].limit.second;
        }

        cppoptlib::solver::State<ProFunction::scalar_t> stop_state;
        stop_state.num_iterations = 10000;
        double accu(1e-5);
        //if (curPro->variableAccuracy() > 0) accu = curPro->variableAccuracy();
        stop_state.x_delta = ProFunction::scalar_t{ accu };
        stop_state.x_delta_violations = 5;
        accu = 1e-5;
        //if (curPro->objectiveAccuracy() > 0) accu = curPro->objectiveAccuracy();
        stop_state.f_delta = ProFunction::scalar_t{ accu};
        stop_state.f_delta_violations = 5;
        stop_state.gradient_norm = ProFunction::scalar_t{ 1e-4 };
        stop_state.condition_hessian = ProFunction::scalar_t{ 0 };
        stop_state.status = cppoptlib::solver::Status::NotStarted;

        m_proSolver.reset(new ProSolver(lower_bound, upper_bound, false, stop_state));
    }

    m_local_best = solutionType(curPro->numberObjectives(), curPro->numberConstraints(), curPro->numberVariables());
    m_local_best.initialize(m_problem.get(), m_random.get());
    m_local_best.evaluate(m_problem.get(), this);
}
