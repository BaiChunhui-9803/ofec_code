#ifndef LON_CON_SAMPLING_H
#define LON_CON_SAMPLING_H

#include "../../../../../utility/parameter/param_map.h"
#include "../../../../../utility/numerical_solvers/solver/lbfgsb.h"
#include "../../../../../utility/numerical_solvers/function.h"
#include "../../../../../core/problem/solution.h"
#include "../../../../../core/problem/domain.h"
#include "../../../../../core/problem/problem.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include <vector>

namespace ofec::lon {


    struct LonSamplingPar {
        
        std::shared_ptr<Problem> m_problem = nullptr;
        std::shared_ptr<const ParameterMap> m_param = nullptr;
        std::shared_ptr<Random> m_rnd = nullptr;
        double seed = 0.5;
        std::string ofile = "output.txt";
        std::string proname = "Classic_Ackley";
        int numVar = 2;
        
      //  double fopt = 0.0;
      //  std::string  tstep = "per";
     //  
     //   std::string tinit = "uni";
        int Iterations = 1000;//  ## iter, using R because iter is already used in python
        //int RUN = 300;
        int RUN = 300;
        
        //decimals of the system
        int system_digit = -3;
        int hash_digit = 3;
        double m_step = 0.05;
        

        double m_global_fitness = 0;

        ofec::Domain<ofec::Real> curDomain;
        ofec::Domain<ofec::Real> originDomain;




        void setDefault() {
            seed = 0.5;
            ofile = "output.txt";
            proname = "Classic_Ackley";
            numVar = 2;
        //    fopt = 0.0;
        //    tstep = "per";
       //     step = 0.01;
       //     tinit = "uni";
            Iterations = 100;//  ## iter, using R because iter is already used in python
            RUN = 10;
            system_digit = -1;
            hash_digit = 1;
            m_step = 0.01;
        }

        void setParameters(const ofec::ParameterMap& par) {
            {
                seed = par.get<double>("--seed",seed);
                ofile = par.get<std::string>("--ofile", ofile);
                proname = par.get<std::string>("problem name", proname);
                numVar = par.get<int>("number of variables", numVar);
                //  ofec::getParValue(par, "--fopt", fopt, fopt);
                //ofec::getParValue(par, "--tstep", tstep, tstep);
                //ofec::getParValue(par, "--step", step, step);
                Iterations = par.get<int>("--iter", Iterations);
                RUN = par.get<int>("--runs", RUN);
                system_digit = par.get<int>("--prec", system_digit);
                //      ofec::getParValue(par, "--bounded", prec, prec);
            }
        }

        void initialize(const ofec::ParameterMap& v) {
            
            m_param = std::shared_ptr<const ParameterMap>(new ParameterMap(v));
            m_problem = Problem::generateByFactory(m_param, 0.1);
            m_problem->initialize();
            m_rnd = std::shared_ptr<ofec::Random>(new Random(0.5));;

            int dim = m_problem->numberVariables();
            curDomain.resize(dim);
            for (int idx(0); idx < dim; ++idx) {
                curDomain.setRange(-1, 1, idx);
            };
            originDomain = CAST_CONOP(m_problem.get())->domain();
        }

        void release() {
            m_param = nullptr;
            m_problem = nullptr;
            m_rnd = nullptr;
        }
    };


  




    class LonFunction : public cppoptlib::function::Function<double> {

    protected:
        std::shared_ptr<Problem> m_problem = nullptr;
        int m_dim = -1;
        double m_minimize_flag = 1.0;
        // ofec::Solution<> m_cur_sol;
        ofec::Domain<ofec::Real> m_cur_domain;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using FunctionXd = cppoptlib::function::Function<double>;
        using FunctionXd::hessian_t;
        using FunctionXd::vector_t;

        Problem* getPro()const {
            return m_problem.get();
        }

        int dim()const {
            return m_dim;
        }

        virtual void initialize(const std::shared_ptr<Problem>& pro) {
            using namespace ofec;
            m_problem = pro;
            m_dim = pro->numberVariables();
            if (pro->optimizeMode()[0] == OptimizeMode::kMaximize) {
                m_minimize_flag = -1.0;
            }
            else m_minimize_flag = 1.0;

            m_cur_domain = CAST_CONOP(pro.get())->domain();
        }
        virtual const ofec::Domain<ofec::Real>& curDomain()const {
            return m_cur_domain;
        }



        scalar_t operator()(const vector_t& x) const override {
            using namespace ofec;
            auto pro(m_problem.get());
            ofec::Solution<>  cur_sol(pro->numberObjectives(), pro->numberConstraints(), pro->numberVariables());

            for (int idx(0); idx < m_dim; ++idx) {
                cur_sol.variable()[idx] = x(idx);
            }
            cur_sol.evaluate(m_problem.get());

            return cur_sol.objective()[0] * m_minimize_flag;
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


    extern std::vector<double> m_xx;

    class LonNormalizeFunction : public LonFunction {

    protected:
        ofec::Domain<ofec::Real> m_domain;



    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        using FunctionXd = cppoptlib::function::Function<double>;
        using FunctionXd::hessian_t;
        using FunctionXd::vector_t;




        virtual void initialize(const std::shared_ptr<Problem>& pro) override{
            LonFunction::initialize(pro);
            using namespace ofec;
            m_domain = CAST_CONOP(pro.get())->domain();
            m_cur_domain.resize(m_dim);
            for (int idx(0); idx < m_dim; ++idx) {
                m_cur_domain.setRange(-1, 1, idx);
            }
        }

        void transfer(const vector_t& x, const ofec::VariableBase& xx) {

        }

        scalar_t operator()(const vector_t& x) const override {
            using namespace ofec;
            auto pro(m_problem);

            m_xx.resize(pro->numberVariables());
            ofec::Solution<>  cur_sol(pro->numberObjectives(), pro->numberConstraints(), pro->numberVariables());

            for (int idx(0); idx < m_dim; ++idx) {

                m_xx[idx] = x(idx);
                cur_sol.variable()[idx] = ofec::mapReal(x(idx), m_cur_domain.range(idx).limit.first, m_cur_domain.range(idx).limit.second, 
                    m_domain[idx].limit.first, m_domain[idx].limit.second);
     
            }
            cur_sol.evaluate(m_problem.get());

            return cur_sol.objective()[0] * m_minimize_flag;
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


    struct OptimaSta {
        int m_numOpts = 0;
        int m_numFoundOpts = 0;
        int m_numFoundLocalOpts = 0;

        void clear() {
            m_numOpts = 0;
            m_numFoundOpts = 0;
            m_numFoundLocalOpts = 0;
        }
    };

    struct OptimaInfo {
        OptimaSta m_globalOpt;
        OptimaSta m_localOpt;

        OptimaSta m_Opts;

        int m_outerOpts = 0;

        double m_maxDis = 0;
        //double m_minDis = 0;
        double m_meanDis = 0;
        double m_stdDis = 0;

        int m_numOptInGrid = 0;


        std::pair<double, double> m_firstNotOptimaDis;
        std::pair<double, double> m_lastOptDis;
        std::pair<double, double> m_dis_threadhold;
        

        void clear() {
            m_globalOpt.clear();
            m_localOpt.clear();
            m_Opts.clear();
            m_outerOpts = 0;
            m_maxDis = 0;
            //m_minDis = 0;
            m_meanDis = 0;
            m_stdDis = 0;
            m_numOptInGrid = 0;
            m_firstNotOptimaDis.first = 0;
            m_firstNotOptimaDis.second = 0;
            m_lastOptDis.first = 0;
            m_lastOptDis.second = 0;
            m_dis_threadhold.first = m_dis_threadhold.second = 0;
            
        }

        static void outputHead(std::ostream& out) {
            out << "NumOpts,\tNumFoundOpts,\tNumFoundLOpt,\t";
            out << "NumLocalOpts,\tNumFoundOpts,\tNumFoundLopt,\t";
            out << "NumOpts,\tNumFoundOpts,\t";
            out << "outerOpts,\t";
            out << "maxDis,\tmeanDis,\tstdDis,\t";
        }
        void output(std::ostream& out) {
            out << m_globalOpt.m_numOpts << ",\t" << m_globalOpt.m_numFoundOpts << ",\t" << m_globalOpt.m_numFoundLocalOpts << ",\t";
            out << m_localOpt.m_numOpts << ",\t" << m_localOpt.m_numFoundOpts << ",\t" << m_localOpt.m_numFoundLocalOpts << ",\t";
            out << m_Opts.m_numOpts << ",\t" << m_Opts.m_numFoundOpts << ",\t";
            out << m_outerOpts << ",\t";
            out << m_maxDis << ",\t" << m_meanDis << ",\t" << m_stdDis << ",\t";
        }
    };


    struct LonNode {
        std::vector<double> m_variables;
        std::string m_x_hash;
        double m_fit = 0;
        std::string m_fit_hash;
        
        void updateState(LonFunction::state_t solutionBest,int dim) {
            m_variables.resize(dim);
            for (int idx(0); idx < dim; ++idx) {
                m_variables[idx]= solutionBest.x(idx) ;
            }
            m_fit= solutionBest.value ;
        }
    };

    

    struct LonStruct {
        std::vector<std::vector<LonNode>> m_runs;

        void clear() {
            m_runs.clear();
        }
    };

  
    using LonSolver = cppoptlib::solver::Lbfgsb<LonFunction>;
    using LonNorSolver = cppoptlib::solver::Lbfgsb<LonNormalizeFunction>;
    extern LonFunction::vector_t BoundedPerturbation(Random *rnd, const LonFunction::vector_t& x,
        const std::vector<double>& p, const ofec::Domain<ofec::Real> & domain);

    extern double lon_round(double var, int digit);
    extern double AverageEscapes(Random *rnd,const LonFunction& fun, LonSolver& solver, double step, int digits);
    extern double getTunningStep(Random *rnd, const LonFunction& fun, LonSolver& solver);
    
    extern void testStepTunning(Random *rnd, const LonFunction& fun, LonSolver& solver);
   
    // normalize the solution in cubic of [0,1] 
    extern void outputSolution(std::ofstream& out, const LonSolver::function_state_t& solution, const ofec::Domain<ofec::Real>& domain, int  system_digit );

    extern void BasicHoppingSamplingFun(const LonSamplingPar& lonPar, Random *rnd, const LonFunction& fun, LonSolver& solver);
    
    extern void localOptimaNetworkConSampling(const ofec::ParameterMap& par, LonSamplingPar&lonpar);

    extern void BasicHoppingSamplingFunPro(const LonSamplingPar& lonPar, Random *rnd,
        const LonNormalizeFunction& fun, LonNorSolver& solver, OptimaInfo& info);

    extern void localOptimaNetworkConSamplinPro(const ofec::ParameterMap& par, OptimaInfo& info);


  //  extern void BasicHoppingSamplingFunNormalize(const LonSamplingPar& lonPar, Random *rnd, const LonFunction& fun, LonSolver& solver);

    extern void localOptimaNetworkConSamplingNormalize(const ofec::ParameterMap& par, LonSamplingPar& lonpar);
   

    extern void BasicHoppingSamplingFun(const LonSamplingPar& lonPar, Random *rnd, const LonFunction& fun, LonSolver& solver, LonStruct& lon);

    extern void localOptimaNetworkConSamplingLON(const LonSamplingPar& lonpar, LonStruct & lon);

    extern std::vector<double> mapX(const std::vector<double>& x,
        const  ofec::Domain<ofec::Real>& originDomain, const   ofec::Domain<ofec::Real>& curDomain);
    

    extern void filterLON( LonStruct& lon, const LonSamplingPar& par);
    extern void outputLON(std::ostream& out, LonStruct& lon);

    extern void analyzeLON(const LonSamplingPar& par, const LonStruct& lon, OptimaInfo & optInfo);
    
    
}

#endif