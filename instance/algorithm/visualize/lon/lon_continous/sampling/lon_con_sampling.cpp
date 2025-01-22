#include "lon_con_sampling.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include <iostream>
#include <fstream>
#include "../experiment/include/Tools.h"
#include <set>

namespace ofec::lon {


        std::vector<double> mapX(const std::vector<double>& x, 
        const  ofec::Domain<ofec::Real>& originDomain,const   ofec::Domain<ofec::Real>& curDomain) {
        std::vector<double> newX(x);
        for (int idx(0); idx < newX.size(); ++idx) {
            newX[idx] = ofec::mapReal(x[idx], originDomain.range(idx).limit.first, originDomain.range(idx).limit.second,
                curDomain[idx].limit.first, curDomain[idx].limit.second);

        }
        return std::move(newX);
    }

    std::vector<double> m_xx;

    LonFunction::vector_t BoundedPerturbation(Random *rnd, const LonFunction::vector_t& x,
        const std::vector<double>& p, const ofec::Domain<ofec::Real>& domain) {
        LonFunction::vector_t y(p.size());
        for (int idx(0); idx < p.size(); ++idx) {
            y(idx) = x(idx) + rnd->uniform.nextNonStd<double>(-p[idx], p[idx]);
            if (y(idx) < domain[idx].limit.first) y(idx) = domain[idx].limit.first;
            else if (y(idx) > domain[idx].limit.second)y(idx) = domain[idx].limit.second;
        }
        return std::move(y);
    }

    double lon_round(double var, int digit)
    {

        double digit10 = pow(10, digit);
        double value = (int)(var * digit10 + .5);
        return (value / digit10);
    }

    double AverageEscapes(Random *rnd, const LonFunction& fun, LonSolver& solver, double step, int digits) {

        using namespace ofec;
        int RUN = 100;
        int SUBS = 100;
        //  int my_digits = 4;


        auto pro(fun.getPro());
        ofec::Solution<>  cur_sol(pro->numberObjectives(), pro->numberConstraints(), pro->numberVariables());

        int dim = pro->numberVariables();
        std::vector<double> escapes(RUN, 0);
        LonFunction::vector_t x0(dim);


        //auto conPro = GET_CONOP(fun.idPro());
        auto& curDomain(fun.curDomain());
        std::vector<double> p(dim, 0);
        for (int idx(0); idx < dim; ++idx) {
            p[idx] = step * (curDomain[idx].limit.second - curDomain[idx].limit.first);
            p[idx] = lon_round(p[idx], digits);
        }

        for (int r(0); r < RUN; ++r) {
            cur_sol.initialize(fun.getPro(), rnd);
            for (int idx(0); idx < dim; ++idx) {
                x0(idx) = cur_sol.variable()[idx];
            }
            auto [solution, solver_state] = solver.Minimize(fun, x0);
            auto x1 = solution.x;
            auto x1_c = x1;
            for (int idx(0); idx < dim; ++idx) {
                x1_c(idx) = lon_round(x1(idx), digits);
            }
            double count = 0;

            for (int s(0); s < SUBS; ++s) {
                auto xt = BoundedPerturbation(rnd, x1, p, curDomain);
                auto [solution, solver_state] = solver.Minimize(fun, xt);
                auto xs = solution.x;
                auto xs_c = xs;
                for (int idx(0); idx < dim; ++idx) {
                    xs_c(idx) = lon_round(xs_c(idx), digits);
                }
                bool flag_unequal(true);
                for (int idx(0); idx < dim; ++idx) {
                    if (x1_c(idx) == xs_c(idx)) {
                        flag_unequal = false;
                        break;
                    }
                }
                if (flag_unequal) {
                    count = count + 1;
                }
            }
            escapes[r] = count / SUBS;
        }
        double avg(0);
        ofec::calMean(escapes, avg);

        return avg;
    }

    double getTunningStep(Random *rnd, const LonFunction& fun, LonSolver& solver) {
        using namespace ofec;
        int digits = 4;
        //  std::cout << GET_PRO(fun.idPro()).name() << "\t" << "D\t" << GET_PRO(fun.idPro()).numberVariables() << std::endl;
        double inc = 1.0;
        double old_step = 0.0;
        double test_step = 0.1;
        double new_step = -1.0;
        double avg_old = -1.0;
        double avg_new = -1.0;
        double avg_test = -1.0;
        bool best_found = false;
        int dig = 1;
        double avg(0);
        int dim = fun.getPro()->numberVariables();
        std::vector<double> size(dim, 0);
    //    auto conPro = GET_CONOP(fun.idPro());
        auto& curDomain(fun.curDomain());
        while (dig <= digits && best_found == false) {
            inc = inc / 10.0;
            while (true) {
                avg = AverageEscapes(rnd, fun, solver, test_step, digits);
                if (lon_round(test_step, digits) == 1.0) {
                    best_found = true;
                    break;
                }
                if (avg < 0.5) {
                    old_step = test_step;
                    avg_old = avg;
                    test_step = test_step + inc;
                    if (lon_round(test_step, digits) == lon_round(new_step, digits)) {
                        if (digits == dig) {
                            best_found = true;
                            break;
                        }
                        else {
                            test_step = test_step - inc + inc / 10.0;
                            break;
                        }
                    }
                }
                else if (avg == 0.5) {
                    best_found = true;
                    break;
                }
                else {
                    if (dig == digits) {
                        best_found = true;
                        break;
                    }
                    new_step = test_step;
                    avg_new = avg;
                    if (dig < digits) {
                        test_step = old_step + inc / 10.0;
                    }
                    break;
                }
            }
            dig += 1;
        }


        double a = abs(avg_old - 0.5);
        double b = abs(avg - 0.5);
        double c = abs(avg_new - 0.5);

        double step(0), escapes(0), error(0);

        if (a < b && a < c) {
            step = old_step;
            escapes = avg_old;
            error = a;
        }
        else if (c < a && c < b) {
            step = new_step;
            escapes = avg_new;
            error = c;
        }
        else {
            step = test_step;
            escapes = avg;
            error = b;
        }

        return step;
    }

    void testStepTunning(Random *rnd, const LonFunction& fun, LonSolver& solver) {
        using namespace ofec;
        int digits = 4;
        std::cout << fun.getPro()->name() << "\t" << "D\t" << fun.getPro()->numberVariables() << std::endl;
        double inc = 1.0;
        double old_step = 0.0;
        double test_step = 0.1;
        double new_step = -1.0;
        double avg_old = -1.0;
        double avg_new = -1.0;
        double avg_test = -1.0;
        bool best_found = false;
        int dig = 1;
        double avg(0);

        int dim =fun.getPro()->numberVariables();
        std::vector<double> size(dim, 0);
        auto conPro = CAST_CONOP(fun.getPro());
        auto& curDomain(conPro->domain());

        int outterLoop = 0;
        int innerLoop = 0;

        while (dig <= digits && best_found == false) {
            inc = inc / 10.0;
            outterLoop++;
            while (true) {
                innerLoop++;

                if (outterLoop == 3 && innerLoop == 9) {
                    int stop = -1;
                }
                avg = AverageEscapes(rnd, fun, solver, test_step, digits);
                std::cout << "outerLoopt \t " << outterLoop << "::: inner loop \t" << innerLoop << std::endl;
                std::cout << "old step " << old_step << " :::avg " << avg_old << std::endl;
                std::cout << "new step " << new_step << " :::avg " << avg_new << std::endl;
                std::cout << "test step " << test_step << " :::avg " << avg << std::endl;

                for (int idx(0); idx < dim; ++idx) {
                    size[idx] = test_step * (curDomain[idx].limit.second - curDomain[idx].limit.first);

                    std::cout << size[idx] << "\t";
                }
                std::cout << std::endl;
                std::cout << std::endl << std::endl;


                if (lon_round(test_step, digits) == 1.0) {
                    best_found = true;
                    break;
                }

                if (avg < 0.5) {
                    old_step = test_step;
                    avg_old = avg;
                    test_step = test_step + inc;
                    if (lon_round(test_step, digits) == lon_round(new_step, digits)) {
                        if (digits == dig) {
                            best_found = true;
                            break;
                        }
                        else {
                            test_step = test_step - inc + inc / 10.0;
                            break;
                        }
                    }
                }
                else if (avg == 0.5) {
                    best_found = true;
                    break;
                }
                else {
                    if (dig == digits) {
                        best_found = true;
                        break;
                    }
                    new_step = test_step;
                    avg_new = avg;
                    if (dig < digits) {
                        test_step = old_step + inc / 10.0;
                    }
                    break;
                }
            }
            dig += 1;
        }


        double a = abs(avg_old - 0.5);
        double b = abs(avg - 0.5);
        double c = abs(avg_new - 0.5);

        double step(0), escapes(0), error(0);

        if (a < b && a < c) {
            step = old_step;
            escapes = avg_old;
            error = a;
        }
        else if (c < a && c < b) {
            step = new_step;
            escapes = avg_new;
            error = c;
        }
        else {
            step = test_step;
            escapes = avg;
            error = b;
        }


        std::vector<double> p(dim, 0);
        for (int idx(0); idx < dim; ++idx) {
            p[idx] = step * (curDomain[idx].limit.second - curDomain[idx].limit.first);
            p[idx] = lon_round(p[idx], digits);
        }

        std::cout << "Selected step size: " << std::endl;
        std::cout << "per :\t" << step << std::endl;
        std::cout << "escapes :\t" << escapes << std::endl;
        std::cout << "error :\t" << error << std::endl;
        std::cout << "beta :\t";
        for (int idx(0); idx < dim; ++idx) {
            std::cout << p[idx] << "\t";
        }
        std::cout << std::endl;
    }

    void outputSolution(std::ofstream& out, const LonSolver::function_state_t& solution, const ofec::Domain<ofec::Real>& domain, int  system_digit) {

        double curVal(0);
        out << solution.value << "\t";
        for (int idx(0); idx < domain.size(); ++idx) {
            curVal = ofec::mapReal<double>(solution.x(idx), domain[idx].limit.first, domain[idx].limit.second, 0, 1);


            if (system_digit > 0) {
                curVal = lon_round(curVal, system_digit);
            }
            out << curVal << " ";
        }
        // out << "\t";

    }


    void BasicHoppingSamplingFun(const LonSamplingPar& lonPar, Random *rnd,
        const LonFunction& fun, LonSolver& solver) {
        using namespace ofec;
        auto conPro = CAST_CONOP(fun.getPro());
        auto& curDomain(fun.curDomain());
        auto& proDomain(conPro->domain());
        int dim = conPro->numberVariables();

        std::vector<double> p(dim, 0);
        for (int idx(0); idx < dim; ++idx) {
            p[idx] = lonPar.m_step * (curDomain[idx].limit.second - curDomain[idx].limit.first);
            //    p[idx] = lon_round(p[idx], digits);
         //       p[idx] = 1.63;
        }


        std::cout << "lon con sampleing path\t" << lonPar.ofile << std::endl;
        std::ofstream out(lonPar.ofile);


        out << "Run\tValue_Current_Best\tCurrent_Best " << dim << "\tValue_Solution_At_Iteration\tSolution_At_Iteration\n";
        ofec::Solution<>  cur_sol(conPro->numberObjectives(), conPro->numberConstraints(), conPro->numberVariables());

        LonFunction::vector_t x0(dim);

        bool flag_improvement(false);
        for (int run(0); run < lonPar.RUN; ++run) {
            flag_improvement = false;
            cur_sol.initialize(fun.getPro(), rnd);
            for (int idx(0); idx < dim; ++idx) {

                x0(idx) = ofec::mapReal(cur_sol.variable()[idx], proDomain.range(idx).limit.first, proDomain.range(idx).limit.second,
                    curDomain[idx].limit.first, curDomain[idx].limit.second);

            }
            auto [solutionBest, solver_state] = solver.Minimize(fun, x0);


            if (lonPar.system_digit > 0) {
                for (int idx(0); idx < dim; ++idx) {
                    solutionBest.x(idx) = lon_round(solutionBest.x(idx), lonPar.system_digit);
                }
                solutionBest.value = lon_round(solutionBest.value, lonPar.system_digit);
            }

            int r(1);
            while (r <= lonPar.Iterations) {
                auto x = BoundedPerturbation(rnd, solutionBest.x, p, curDomain);
                auto [solution, solver_state] = solver.Minimize(fun, x);


                if (lonPar.system_digit > 0) {
                    for (int idx(0); idx < dim; ++idx) {
                        solution.x(idx) = lon_round(solution.x(idx), lonPar.system_digit);
                    }
                    solution.value = lon_round(solution.value, lonPar.system_digit);
                }

                if (solution.value <= solutionBest.value) {

                    out << run + 1 << "\t";

                    //                    outputSolution(out, solutionBest, curDomain, lonPar.system_digit);

                    out << solutionBest.value << "\t";
                    for (int idx(0); idx < dim; ++idx) {
                        out << solutionBest.x(idx) << " ";
                    }
                    out << "\t";

                    solutionBest = solution;

                    //    outputSolution(out, solutionBest, curDomain, lonPar.system_digit);


                    out << solutionBest.value << "\t";
                    for (int idx(0); idx < dim; ++idx) {
                        out << solutionBest.x(idx) << " ";
                    }
                    out << std::endl;



                    flag_improvement = true;


                }

                r++;
            }


            if (!flag_improvement) {
                out << run + 1 << "\t";
                // outputSolution(out, solutionBest, curDomain, lonPar.system_digit);

                out << solutionBest.value << "\t";
                for (int idx(0); idx < dim; ++idx) {
                    out << solutionBest.x(idx) << " ";
                }
                out << "\t";

                //  outputSolution(out, solutionBest, curDomain, lonPar.system_digit);


                out << solutionBest.value << "\t";
                for (int idx(0); idx < dim; ++idx) {
                    out << solutionBest.x(idx) << " ";
                }
                out << std::endl;
            }
        }

        out.close();

    }


    void localOptimaNetworkConSampling(const ofec::ParameterMap& par, LonSamplingPar& lonPar) {
        using namespace ofec;

        //// set paramers
        //LonSamplingPar lonPar;
        //lonPar.setDefault();
        //lonPar.setParameters(par);

       // auto id_param = std::make_shared<const ParameterMap>(par);
        auto id_param = std::shared_ptr<const ofec::ParameterMap>(new ofec::ParameterMap(par));

        auto pro = Problem::generateByFactory(id_param, 0.1);
        pro->initialize();
        double m_minimize_flag = 1.0;
        if (pro->optimizeMode()[0] == OptimizeMode::kMaximize) {
            m_minimize_flag = -1.0;
        }
        auto optima = CAST_CONOP(pro.get())->optima();
        if (optima->isObjectiveGiven()) {
         lonPar.m_global_fitness=   optima->objective().at(0)* m_minimize_flag;
        }
        
        LonFunction fun;
        fun.initialize(pro);
        using Solver = cppoptlib::solver::Lbfgsb<LonFunction>;
        auto& domain(fun.curDomain());
        LonFunction::vector_t lower_bound(pro->numberVariables());
        LonFunction::vector_t upper_bound(pro->numberVariables());
        for (int idx(0); idx < pro->numberVariables(); ++idx) {
            lower_bound(idx) = domain[idx].limit.first;
            upper_bound(idx) = domain[idx].limit.second;
        }

        cppoptlib::solver::State<LonFunction::scalar_t> stop_state;
        stop_state.num_iterations = 10000;
        stop_state.x_delta = LonFunction::scalar_t{ 1e-7 };
        //stop_state.x_delta_violations = 5;
        stop_state.f_delta = LonFunction::scalar_t{ 1e-5 };
        //stop_state.f_delta_violations = 5;
        stop_state.gradient_norm = LonFunction::scalar_t{ 1e-5 };
        stop_state.condition_hessian = LonFunction::scalar_t{ 0 };
        stop_state.status = cppoptlib::solver::Status::NotStarted;

        Solver solver(lower_bound, upper_bound, false, stop_state);

        Random rnd (lonPar.seed);
        // int digits = 100;
        // AverageEscapes(rnd, fun, solver, step, digits);
 //        testStepTunning(rnd, fun, solver);

      //  lonPar.m_step = getTunningStep(rnd, fun, solver);
        BasicHoppingSamplingFun(lonPar, &rnd, fun, solver);


        {


        }



    }




    void BasicHoppingSamplingFunPro(const LonSamplingPar& lonPar, Random *rnd,
        const LonNormalizeFunction& fun, LonNorSolver& solver, OptimaInfo& info) {



        using namespace ofec;
        auto& curDomain = fun.curDomain();
        int dim = fun.dim();

        std::vector<double> p(dim, 0);
        for (int idx(0); idx < dim; ++idx) {
            p[idx] = lonPar.m_step * (curDomain[idx].limit.second - curDomain[idx].limit.first);
            //    p[idx] = lon_round(p[idx], digits);
         //       p[idx] = 1.63;
        }

        auto conPro = CAST_CONOP(fun.getPro());
        int hashDigit = 1;


        struct NodeInfo {


            enum Tag {None,Global,Local,Funnel};
            int m_id = -1;
            Solution<> curSol;
            std::string hash;
            bool m_found = false;
            // m_tag= 0 global,tag = 1, local optima, tag =2 funnel
            Tag m_tag = None;
            std::vector<double> dis;

            NodeInfo() = default;
            

            NodeInfo(const NodeInfo& info) :m_id(info.m_id),curSol(info.curSol),hash(info.hash),
                m_found(info.m_found),m_tag(info.m_tag),dis(info.dis)
            {

            }
        };
        

        //std::vector<double> dis;
        //std::vector<bool> flag;
        

        std::map<std::string, int> hashToId;
        std::vector<NodeInfo> totalInfos;

        //std::vector<Solution<>> globalOpt;
        //std::vector<std::string> globalOptHash;
        //std::vector<Solution<>> localOpt;
        //std::vector<std::string> localOptHash;

      //  int numVar = conPro->optima()->numberVariables();
        std::list<Solution<>> totalOpts;
         double maxFit(-std::numeric_limits<double>::max());
        double modeValue = 1.0;
        if (conPro->optimizeMode(0) == OptimizeMode::kMinimize) {
            modeValue = -1.0;
        }

        for (int idx(0); idx < conPro->optima()->numberVariables(); ++idx) {
            Solution<> curSol(conPro->numberObjectives(), conPro->numberConstraints(), conPro->numberVariables());

            curSol.variable() = dynamic_cast<const VariableVector<Real>&>(conPro->optima()->variableBase(idx));
            curSol.evaluate(fun.getPro());
            curSol.setFitness(curSol.objective()[0] * modeValue);
            totalOpts.push_back(curSol);
            maxFit = std::max(maxFit, curSol.fitness());
        }

        NodeInfo curInfo;

        while (!totalOpts.empty()) {


            if (abs(maxFit - totalOpts.back().fitness()) <= conPro->objectiveAccuracy()) {
                curInfo.m_tag = NodeInfo::Global;
            }
            else {
                curInfo.m_tag = NodeInfo::Local;
            }
            auto hashValue = toHashStr(totalOpts.back().variable().vect(), hashDigit);
            curInfo.hash = hashValue;
            curInfo.curSol = std::move(totalOpts.back());
   
            if (hashToId.find(curInfo.hash) == hashToId.end()) {
                curInfo.m_id = totalInfos.size();
                hashToId[curInfo.hash] = totalInfos.size();
            //    totalInfos.push_back(NodeInfo());
             //   totalInfos.back() = curInfo;
                totalInfos.push_back(curInfo);
            }
            else {
                curInfo.m_id = hashToId.at(curInfo.hash);
            }


            if (totalInfos.back().curSol.variablePointer() == nullptr) {
                int stop = -1;
            }
            
            totalOpts.pop_back();
        }

       ofec::Solution<>  cur_sol(conPro->numberObjectives(), conPro->numberConstraints(), conPro->numberVariables());
        LonFunction::vector_t x0(dim);
      //  std::vector<double> vecx(mat.data(), mat.data() + mat.rows() * mat.cols());
        std::set<std::string> lon_node;
        std::vector<NodeInfo> funnel_infos;
        bool flag_improvement(false);

        
        for (int run(0); run < lonPar.RUN; ++run) {

            std::set<std::string> hash_pathset;
            std::vector<std::string> path;
            flag_improvement = false;
            cur_sol.initialize(fun.getPro(), rnd);
            for (int idx(0); idx < dim; ++idx) {
                x0(idx) = cur_sol.variable()[idx];
            }
            auto [solutionBest, solver_state] = solver.Minimize(fun, x0);


            if (lonPar.system_digit > 0) {
                for (int idx(0); idx < dim; ++idx) {
                    solutionBest.x(idx) = lon_round(solutionBest.x(idx), lonPar.system_digit);
                }
                solutionBest.value = lon_round(solutionBest.value, lonPar.system_digit);
            }



            int r(1);
            while (r <= lonPar.Iterations) {
                auto x = BoundedPerturbation(rnd, solutionBest.x, p, curDomain);
                auto [solution, solver_state] = solver.Minimize(fun, x);


                if (lonPar.system_digit > 0) {
                    for (int idx(0); idx < dim; ++idx) {
                        solution.x(idx) = lon_round(solution.x(idx), lonPar.system_digit);
                    }
                    solution.value = lon_round(solution.value, lonPar.system_digit);
                }

                if (solution.value <= solutionBest.value) {

                    for (int idx(0); idx < dim; ++idx) {
                        cur_sol.variable()[idx] = solutionBest.x(idx);
                    }
                    cur_sol.evaluate(fun.getPro());
                    curInfo.hash = toHashStr(cur_sol.variable().vect(), hashDigit);
                    lon_node.insert(curInfo.hash);

                    solutionBest = solution;
                    flag_improvement = true;

                }

                r++;
            }
      //      ofec::Solution<>  cur_sol(pro->numberObjectives(), pro->numberConstraints(), pro->numberVariables());
            for (int idx(0); idx < dim; ++idx) {
                cur_sol.variable()[idx] = solutionBest.x(idx);
            }
            cur_sol.evaluate(fun.getPro());
            curInfo.hash = toHashStr(cur_sol.variable().vect(), hashDigit);
            curInfo.curSol = cur_sol;
            if (curInfo.curSol.variablePointer()==nullptr) {
                int stop = -1;
            }
      //      funnel_infos.push_back(NodeInfo());
       //     funnel_infos.back() = curInfo;

           funnel_infos.push_back(curInfo);
            if (funnel_infos.back().curSol.variablePointer() == nullptr) {
                int stop = -1;
            }
        }

        for (auto& it : funnel_infos) {
            if (lon_node.find(it.hash) == lon_node.end()) {
                if (hashToId.find(it.hash) == hashToId.end()) {
                    curInfo = it;
                    curInfo.m_id = totalInfos.size();
                    hashToId[curInfo.hash] = totalInfos.size();
                    curInfo.m_tag = NodeInfo::Funnel;
                    totalInfos.emplace_back(curInfo);
                }
                else {
                    curInfo.m_id = hashToId.at(curInfo.hash);
                    auto& center = totalInfos[curInfo.m_id];
                    curInfo.dis.push_back(center.curSol.variableDistance(it.curSol, fun.getPro()));
                }
            }
        }


        std::vector<double> dis2opt;
        std::vector<double> minDis;
        for (auto& it : totalInfos) {
            if (it.m_tag == NodeInfo::Global) {
                if (!it.dis.empty()) {
                    info.m_globalOpt.m_numFoundOpts++;
                }
                info.m_globalOpt.m_numOpts++;
                for (auto& curdis : it.dis) dis2opt.push_back(curdis);
            }
            else if (it.m_tag == NodeInfo::Local) {
                if (!it.dis.empty()) info.m_localOpt.m_numFoundOpts++;
                ++info.m_localOpt.m_numOpts;
                for (auto& curdis : it.dis) dis2opt.push_back(curdis);
            }
            else {
                ++info.m_outerOpts;
                double dis(0);
                
            }
        }
        if(!dis2opt.empty())
        calMeanAndStd(dis2opt, info.m_meanDis, info.m_stdDis);
        info.m_maxDis = 0;
        for (auto& it : dis2opt) {
            info.m_maxDis = std::max(it, info.m_maxDis);
        }
        

    }

    void localOptimaNetworkConSamplinPro(const ofec::ParameterMap& par, OptimaInfo& info) {
        using namespace ofec;
        // set paramers
        LonSamplingPar lonPar;
        lonPar.setDefault();
        lonPar.setParameters(par);

        auto id_param = std::shared_ptr<const ParameterMap>(new ParameterMap(par));

        auto pro = Problem::generateByFactory(id_param, 0.1);
        pro->initialize();
        LonNormalizeFunction fun;
        fun.initialize(pro);
     //   using Solver = cppoptlib::solver::Lbfgsb<LonNormalizeFunction>;
        LonNormalizeFunction::vector_t lower_bound(pro->numberVariables());
        LonNormalizeFunction::vector_t upper_bound(pro->numberVariables());
        for (int idx(0); idx < pro->numberVariables(); ++idx) {
            lower_bound(idx) = fun.curDomain()[idx].limit.first;
            upper_bound(idx) = fun.curDomain()[idx].limit.second;
        }

        //cppoptlib::solver::State<LonFunction::scalar_t> stop_state;

        cppoptlib::solver::State<LonFunction::scalar_t> stop_state = cppoptlib::solver::DefaultStoppingSolverState<LonFunction::scalar_t>();

     //   stop_state.num_iterations = 10000;
        stop_state.x_delta = LonFunction::scalar_t{ 1e-7 };
      //  stop_state.x_delta_violations = 5;
        stop_state.f_delta = LonFunction::scalar_t{ 1e-5 };
       // stop_state.f_delta_violations = 5;
        stop_state.gradient_norm = LonFunction::scalar_t{ 1e-5 };
        stop_state.condition_hessian = LonFunction::scalar_t{ 0 };
        stop_state.status = cppoptlib::solver::Status::NotStarted;

        LonNorSolver solver(lower_bound, upper_bound, false, stop_state);

        Random rnd(lonPar.seed);
        // int digits = 100;
        // AverageEscapes(rnd, fun, solver, step, digits);
 //        testStepTunning(rnd, fun, solver);

     //   lonPar.m_step = getTunningStep(rnd, fun, solver);
        BasicHoppingSamplingFunPro(lonPar, &rnd, fun, solver, info);





    }



    void localOptimaNetworkConSamplingNormalize(const ofec::ParameterMap& par, LonSamplingPar& lonPar) {
        using namespace ofec;

        //// set paramers
        //LonSamplingPar lonPar;
        //lonPar.setDefault();
        //lonPar.setParameters(par);

        auto id_param = std::shared_ptr<const ParameterMap>(new ParameterMap(par));
        
        auto pro = Problem::generateByFactory(id_param, 0.1);
        pro->initialize();
        double m_minimize_flag = 1.0;
        if (pro->optimizeMode()[0] == OptimizeMode::kMaximize) {
            m_minimize_flag = -1.0;
        }
        auto optima = CAST_CONOP(pro.get())->optima();
        if (optima->isObjectiveGiven()) {
            lonPar.m_global_fitness = optima->objective().at(0) * m_minimize_flag;
        }

        LonNormalizeFunction fun;
        fun.initialize(pro);
        using Solver = cppoptlib::solver::Lbfgsb<LonFunction>;
        auto& domain(fun.curDomain());
        LonNormalizeFunction::vector_t lower_bound(pro->numberVariables());
        LonNormalizeFunction::vector_t upper_bound(pro->numberVariables());
        for (int idx(0); idx < pro->numberVariables(); ++idx) {
            lower_bound(idx) = domain[idx].limit.first;
            upper_bound(idx) = domain[idx].limit.second;
        }


        cppoptlib::solver::State<LonNormalizeFunction::scalar_t> stop_state = cppoptlib::solver::DefaultStoppingSolverState<LonNormalizeFunction::scalar_t>();

      //  stop_state.num_iterations = 10000;
        stop_state.x_delta = LonNormalizeFunction::scalar_t{ 1e-7 };
        //stop_state.x_delta_violations = 5;
        stop_state.f_delta = LonNormalizeFunction::scalar_t{ 1e-5 };
        //stop_state.f_delta_violations = 5;
        stop_state.gradient_norm = LonNormalizeFunction::scalar_t{ 1e-5 };
        stop_state.condition_hessian = LonNormalizeFunction::scalar_t{ 0 };
        stop_state.status = cppoptlib::solver::Status::NotStarted;

        Solver solver(lower_bound, upper_bound, false, stop_state);

        Random rnd (lonPar.seed);
        // int digits = 100;
        // AverageEscapes(rnd, fun, solver, step, digits);
 //        testStepTunning(rnd, fun, solver);

      //  lonPar.m_step = getTunningStep(rnd, fun, solver);
        BasicHoppingSamplingFun(lonPar, &rnd, fun, solver);


        {


        }

    }

    void BasicHoppingSamplingFun(const LonSamplingPar& lonPar, Random *rnd, const LonFunction& fun, LonSolver& solver, LonStruct& lon) {
        using namespace ofec;
        auto conPro = CAST_CONOP(fun.getPro());
        auto& curDomain(fun.curDomain());
        auto& proDomain(conPro->domain());
        int dim = conPro->numberVariables();

        std::vector<double> p(dim, 0);
        for (int idx(0); idx < dim; ++idx) {
            p[idx] = lonPar.m_step * (curDomain[idx].limit.second - curDomain[idx].limit.first);
            //    p[idx] = lon_round(p[idx], digits);
         //       p[idx] = 1.63;
        }
        LonNode curNode;
        lon.m_runs.clear();

      //  out << "Run\tValue_Current_Best\tCurrent_Best " << dim << "\tValue_Solution_At_Iteration\tSolution_At_Iteration\n";
        ofec::Solution<>  cur_sol(conPro->numberObjectives(), conPro->numberConstraints(), conPro->numberVariables());

        LonFunction::vector_t x0(dim);
        LonFunction::state_t xx;
        
        bool flag_improvement(false);
        for (int run(0); run < lonPar.RUN; ++run) {
            lon.m_runs.push_back(std::vector<LonNode>());
            auto& cur_trait = lon.m_runs.back();
            flag_improvement = false;
            cur_sol.initialize(fun.getPro(), rnd);
            for (int idx(0); idx < dim; ++idx) {

                x0(idx) = ofec::mapReal(cur_sol.variable()[idx], proDomain.range(idx).limit.first, proDomain.range(idx).limit.second,
                    curDomain[idx].limit.first, curDomain[idx].limit.second);

            }
            auto [solutionBest, solver_state] = solver.Minimize(fun, x0);


            if (lonPar.system_digit > 0) {
                for (int idx(0); idx < dim; ++idx) {
                    solutionBest.x(idx) = lon_round(solutionBest.x(idx), lonPar.system_digit);
                }
                solutionBest.value = lon_round(solutionBest.value, lonPar.system_digit);
            }

            curNode.updateState(solutionBest, dim);
            cur_trait.push_back(curNode);

            int r(1);
            while (r <= lonPar.Iterations) {
                auto x = BoundedPerturbation(rnd, solutionBest.x, p, curDomain);
                auto [solution, solver_state] = solver.Minimize(fun, x);


                if (lonPar.system_digit > 0) {
                    for (int idx(0); idx < dim; ++idx) {
                        solution.x(idx) = lon_round(solution.x(idx), lonPar.system_digit);
                    }
                    solution.value = lon_round(solution.value, lonPar.system_digit);
                }

                if (solution.value <= solutionBest.value) {

               //     out << run + 1 << "\t";
                    //                    outputSolution(out, solutionBest, curDomain, lonPar.system_digit);
                  //  out << solutionBest.value << "\t";
                    //for (int idx(0); idx < dim; ++idx) {
                    //    out << solutionBest.x(idx) << " ";
                    //}
                    //out << "\t";

                    solutionBest = solution;

                    //    outputSolution(out, solutionBest, curDomain, lonPar.system_digit);

                    curNode.updateState(solutionBest, dim);
                    cur_trait.push_back(curNode);

                    //out << solutionBest.value << "\t";
                    //for (int idx(0); idx < dim; ++idx) {
                    //    out << solutionBest.x(idx) << " ";
                    //}
                    //out << std::endl;



                    flag_improvement = true;


                }

                r++;
            }


          }

      //  out.close();
    }

    void localOptimaNetworkConSamplingLON(
        const LonSamplingPar& lonPar, LonStruct& lon) {
        using namespace ofec;

        Problem *pro = lonPar.m_problem.get();
        double m_minimize_flag = 1.0;
        if (pro->optimizeMode()[0] == OptimizeMode::kMaximize) {
            m_minimize_flag = -1.0;
        }
        //auto optima = CAST_CONOP(pro)->optima();
        //if (optima->isObjectiveGiven()) {
        //    lonPar.m_global_fitness = optima->objective().at(0) * m_minimize_flag;
        //}

        LonFunction fun;
        fun.initialize(lonPar.m_problem);
        using Solver = cppoptlib::solver::Lbfgsb<LonFunction>;
        auto& domain(fun.curDomain());
        LonFunction::vector_t lower_bound(pro->numberVariables());
        LonFunction::vector_t upper_bound(pro->numberVariables());
        for (int idx(0); idx < pro->numberVariables(); ++idx) {
            lower_bound(idx) = domain[idx].limit.first;
            upper_bound(idx) = domain[idx].limit.second;
        }

        cppoptlib::solver::State<LonFunction::scalar_t> stop_state = cppoptlib::solver::DefaultStoppingSolverState<LonFunction::scalar_t>();
      //  stop_state.num_iterations = 10000;
        stop_state.x_delta = LonFunction::scalar_t{ 1e-7 };
        //stop_state.x_delta_violations = 5;
        stop_state.f_delta = LonFunction::scalar_t{ 1e-5 };
        //stop_state.f_delta_violations = 5;
        stop_state.gradient_norm = LonFunction::scalar_t{ 1e-5 };
        stop_state.condition_hessian = LonFunction::scalar_t{ 0 };
        stop_state.status = cppoptlib::solver::Status::NotStarted;

        Solver solver(lower_bound, upper_bound, false, stop_state);


        // int digits = 100;
        // AverageEscapes(rnd, fun, solver, step, digits);
 //        testStepTunning(rnd, fun, solver);

      //  lonPar.m_step = getTunningStep(rnd, fun, solver);
        BasicHoppingSamplingFun(lonPar, lonPar.m_rnd.get(), fun, solver, lon);


        {


        }
    }

    void filterLON(LonStruct& lon, const LonSamplingPar& par) {

        using namespace ofec;
        int digit = par.hash_digit;
        LonStruct newLon;
        newLon.m_runs.resize(lon.m_runs.size());
       // std::pair<std::string, std::string> cur_hash;

        std::map<std::string, LonNode> node_map;

        for (auto& it : lon.m_runs) {
            for (auto& curNode : it) {

                curNode.m_fit_hash = toHashStrFit(curNode.m_fit, digit);
                curNode.m_x_hash = toHashStr(curNode.m_variables, digit);

                if (node_map.find(curNode.m_x_hash) == node_map.end()) {
                    node_map[curNode.m_x_hash] = curNode;
                }
                else {
                    auto& map_node = node_map[curNode.m_x_hash];
                    if (map_node.m_fit > curNode.m_fit) {
                        map_node = curNode;
                    }
                }
            }
        }


        for (auto& it : lon.m_runs) {
            for (auto& curNode : it) {
                curNode = node_map[curNode.m_x_hash];
            }
        }

        for (int idRun(0); idRun < lon.m_runs.size(); ++idRun) {
            auto& cur_trait = newLon.m_runs[idRun];
            //std::set<std::pair<std::string, std::string>> cur_trait_hash;
            std::map< std::string, int> map_id;
            std::vector<int> curId(lon.m_runs[idRun].size(), 1e9);
            for (int id(0); id < lon.m_runs[idRun].size(); ++id) {
                auto& curNode = lon.m_runs[idRun][id];

                curNode.m_fit_hash = toHashStrFit(curNode.m_fit, digit);
                curNode.m_x_hash = toHashStr(mapX(curNode.m_variables,par.originDomain,par.curDomain), digit);
                //  cur_hash.first = curNode.m_fit_hash;
                //  cur_hash.second = curNode.m_x_hash;
                if (map_id.find(curNode.m_x_hash) == map_id.end()) {
                    map_id[curNode.m_x_hash] = id;
                }
            }
            //for (auto& curNode : lon.m_runs[idRun]) {

            //    //if (cur_trait_hash.find(cur_hash) == cur_trait_hash.end()) {
            //    //    cur_trait.push_back(curNode);
            //    //    
            //    //}
            //}


            int id = lon.m_runs[idRun].size() - 1;
           // auto& old_trait = lon.m_runs[idRun];
            while (id >= 0) {
                auto& curNode = lon.m_runs[idRun][id];
                cur_trait.push_back(curNode);
                id = map_id[curNode.m_x_hash];
                --id;
            }


            std::reverse(std::begin(cur_trait), std::end(cur_trait));

        }



        std::swap(newLon, lon);
        
    }

    void outputLON(std::ostream& fout, LonStruct& lon){
        using namespace std;
        fout << "Run" << "\t";
        fout << "Fitness" << "\t";
        fout << "Solution" << "\t";
        fout << "Fitness" << "\t";
        fout << "Next_Solution" << "\n";


        
        for (int idRun(0); idRun < lon.m_runs.size(); ++idRun) {
            if (lon.m_runs[idRun].size() > 1) {
                for (int idNode(0); idNode + 1 < lon.m_runs[idRun].size(); ++idNode) {
                    auto& curNode = lon.m_runs[idRun][idNode];
                    auto& nextNode = lon.m_runs[idRun][idNode + 1];
                    fout << idRun + 1 << "\t" << curNode.m_fit_hash << "\t" << curNode.m_x_hash << "\t";
                    fout << nextNode.m_fit_hash << "\t" << nextNode.m_x_hash << std::endl;
                }
            }
            else {
                auto& curNode = lon.m_runs[idRun].front();
                auto& nextNode = lon.m_runs[idRun].back();
                fout << idRun + 1 << "\t" << curNode.m_fit_hash << "\t" << curNode.m_x_hash << "\t";
                fout << nextNode.m_fit_hash << "\t" << nextNode.m_x_hash << std::endl;
            }
        }
    }



    void analyzeLON(const LonSamplingPar& par, const LonStruct& lon, OptimaInfo& optInfo) {
        using namespace ofec;
        auto conPro = CAST_CONOP(par.m_problem.get());

        double modeValue = 1.0;
        if (conPro->optimizeMode(0) == OptimizeMode::kMinimize) {
            modeValue = -1.0;
        }

        double maxFit(-std::numeric_limits<double>::max());
        std::vector<Solution<>> totalOpts;
        for (int idx(0); idx < conPro->optima()->numberVariables(); ++idx) {
            Solution<> curSol(conPro->numberObjectives(), conPro->numberConstraints(), conPro->numberVariables());

            curSol.variable() = dynamic_cast<const VariableVector<Real>&>(conPro->optima()->variableBase(idx));
            curSol.evaluate(par.m_problem.get());
            curSol.setFitness(curSol.objective()[0] * modeValue);
            totalOpts.push_back(curSol);
            maxFit = std::max(maxFit, curSol.fitness());

        }

        struct solInfo {
            LonNode m_cur;
            bool flag = false;
            bool m_found = false;

        };

        std::map<std::string, solInfo> node_state;
        for (auto& curtrait : lon.m_runs) {
            for (int idx(0); idx + 1 < curtrait.size(); ++idx) {
                auto& curNode = curtrait[idx];
                auto& curInfo = node_state[curNode.m_x_hash];
                curInfo.m_cur = curNode;
                curInfo.flag = false;
                // node_state[curNode.m_x_hash] = false;
            }
            auto& lastNode = curtrait.back();
            if (node_state.find(lastNode.m_x_hash) == node_state.end()) {
                auto& curInfo = node_state[lastNode.m_x_hash];

                curInfo.m_cur = lastNode;
                curInfo.flag = true;
            }
        }

        optInfo.clear();
        optInfo.m_Opts.m_numOpts = totalOpts.size();
        
        std::vector<double> dis;
        bool optFlag = false;
        bool foundFlag = false;
        double curDis = 0;

        double funnelFlag = false;
        std::string curHash;

        VariableVector<Real> curx;
        for (auto& optSol : totalOpts) {
            if (abs(optSol.fitness()- maxFit) < par.m_problem.get()->objectiveAccuracy()) {
                optFlag = true;
            }
            else optFlag = false;

            if (optFlag) ++optInfo.m_globalOpt.m_numOpts;
            else ++optInfo.m_localOpt.m_numOpts;

            curHash = toHashStr(mapX(optSol.variable().vect(),par.originDomain,par.curDomain), par.hash_digit);
            foundFlag = false;
            funnelFlag = false;
            if (node_state.find(curHash) != node_state.end()) {
                auto& curInfo = node_state[curHash];
                if (!curInfo.m_found) {
                    curInfo.m_found = true;
                    foundFlag = true;
                    funnelFlag = curInfo.flag;
                    curx = curInfo.m_cur.m_variables;
                }

            }
            if (foundFlag) {
                OptimaSta* optSta = nullptr;
                if (optFlag) {
                    optSta = &optInfo.m_globalOpt;
                }
                else {
                    optSta = &optInfo.m_localOpt;
                }
                if (funnelFlag) ++(optSta->m_numFoundOpts);
                ++(optSta->m_numFoundLocalOpts);
                
                ++optInfo.m_Opts.m_numFoundOpts;
      
                curDis = par.m_problem->normalizedVariableDistance(curx, optSol.variable());
                dis.push_back(curDis);
            }
        }

        
        ofec::calMeanAndStd<double, double>(dis, optInfo.m_meanDis, optInfo.m_stdDis);
        ofec::calMax(dis, optInfo.m_maxDis);

        
        for (auto& curNode: node_state) {
            if (!curNode.second.m_found) {
                ++optInfo.m_outerOpts;
            }
        }


    }
}