/********* Begin Register Information **********
{
	"name": "ConILS",
	"identifier": "IteratedLocalSearch",
	"problem tags": [ "ConOP", "GOP", "MMOP", "SOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_ITERATED_LOCAL_SEARCH_H
#define OFEC_ITERATED_LOCAL_SEARCH_H
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../../../core/problem/solution.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/numerical_solvers/solver/lbfgsb.h"
#include "../../../../../../utility/numerical_solvers/function.h"


/*
* #include "../../../../../utility/numerical_solvers/solver/lbfgsb.h"
#include "../../../../../utility/numerical_solvers/function.h"
*/

namespace ofec {
	class IteratedLocalSearch : public Algorithm {
    public:
        class ProFunction : public cppoptlib::function::Function<double> {

        protected:
            Problem* m_pro=nullptr;
            Algorithm* m_alg = nullptr;
            int m_dim = -1;
            double m_minimize_flag = 1.0;
            // ofec::Solution<> m_cur_sol;

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using FunctionXd = cppoptlib::function::Function<double>;
            using FunctionXd::hessian_t;
            using FunctionXd::vector_t;

            Problem* idPro()const {
                return m_pro;
            }

            void initialize(Environment *env,Algorithm *alg= nullptr) {
                using namespace ofec;
                m_pro = pro;
                m_alg = alg;
                
                m_dim = m_pro->numberVariables();
                if (m_pro->optimizeMode()[0] == OptimizeMode::kMaximize) {
                    m_minimize_flag = -1.0;
                }
            }



            scalar_t operator()(const vector_t& x) const override {
                //using namespace ofec;
                auto pro(m_pro);
                ofec::Solution<>  cur_sol(pro->numberObjectives(), pro->numberConstraints(), pro->numberVariables());

                for (int idx(0); idx < m_dim; ++idx) {
                    cur_sol.variable()[idx] = x(idx);
                }
                cur_sol.evaluate(m_pro, m_alg);

                return cur_sol.objective()[0] * m_minimize_flag;
            }



            static vector_t BoundedPerturbation(Random *rnd, const vector_t& x,
                const std::vector<double>& p, const ofec::Domain<ofec::Real>& domain) {
                vector_t y(p.size());
                for (int idx(0); idx < p.size(); ++idx) {
                    y(idx) = x(idx) + rnd->uniform.nextNonStd<double>(-p[idx], p[idx]);
                    if (y(idx) < domain[idx].limit.first) y(idx) = domain[idx].limit.first;
                    else if (y(idx) > domain[idx].limit.second) y(idx) = domain[idx].limit.second;
                }
                return std::move(y);
            }

            //void Gradient(const vector_t& x, vector_t* grad) const override {
            //    (*grad)[0] = 2 * 5 * x[0];
            //    (*grad)[1] = 2 * 100 * x[1];
            //}

            //void Hessian(const vector_t& x, hessian_t* hessian) const override {
            //    (*hessian)(0, 0) = 10;
            //    (*hessian)(0, 1) = 0;
            //    (*hessian)(1, 0) = 0;
            //    (*hessian)(1, 1) = 200;
            //}
        };

        using  ProSolver = cppoptlib::solver::Lbfgsb<ProFunction>;

	public:
		using variableType = VariableVector<>;
		using solutionType = Solution<variableType>;

	protected:            
        double m_minimize_flag = 1.0;
        solutionType m_local_best ;
		double m_perturbation_ratio = 0.1;
        ofec::Domain<ofec::Real> m_domain;
        std::vector<double> m_perturbation_range;
        std::unique_ptr<ProSolver> m_proSolver = nullptr;
        ProFunction m_proFun;
        
        //std::vector<std::unique_ptr<ofec::SolutionBase>> m_solution_traits;

	protected:

		// Í¨¹ý Algorithm ¼Ì³Ð
		virtual void run_() override;
        virtual void initialize_()override;


   public:
       //void getBestSolution(std::unique_ptr<ofec::SolutionBase>& bestSol)const {
       //    bestSol.reset(nullptr);
       //}
       

	};
}


#endif