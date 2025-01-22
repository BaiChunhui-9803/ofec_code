
#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <limits>

#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "../utility/hnsw/hnsw_solution_multithread/hnsw_sol_model.h"
#include "../core/global.h"
#include "../core/problem/optima.h"
#include "../core/definition.h"
#include "../core/global.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../utility/nbn_visualization/nbn_grid_tree_division.h"
#include "../utility/nbn_visualization/nbn_grid_tree_division_allSize.h"
#include  "../utility/function/custom_function.h"
#include "../utility/function/general_function.h"
#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator.h"

//#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "interface.h"

//#include "../utility/function/custom_function.h"
#include "../utility/nbn_visualization/tree_graph_simple.h"
#include "../core/algorithm/population.h"
//#include "../utility/nondominated_sorting/filter_sort.h"
#include "../instance/algorithm/sampling/instance/sampling_eax_tsp.h"
#include "../instance/algorithm/sampling/sampling_multiThread.h"
#include "../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include <queue>
#include "../core/definition.h"
#include <deque>
#include <chrono>   
//#include "../utility/nbn_visualization/nbn_calculator/nearest_better_calculator_multiThead.h"

#include "../utility/hnsw/hnsw_nbn/neighbor_info.h"
#include "../utility/hnsw/hnsw_nbn/select_policy.h"
#include "../utility/hnsw/hnsw_nbn/node.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_model.h"
#include "../utility/hnsw/hnsw_nbn/hnsw_nbn.h"

#include "../core/environment/environment.h"
#include "../core/parameter/variants_to_stream.h"

//#include "../instance/algorithm/combination/lkh/lkh_algorithm.h"
#include<filesystem>
#include<queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_test/eax_tsp_alg.h"



using namespace std;




bool compareTwoRealVariant(const ofec::Real& data1, const ofec::Real& data2) {
    std::string info1, info2;
    std::stringstream tmp;
    tmp << data1;
    info1 = tmp.str();
    tmp.str("");

    //ofec::Real tmp;

    tmp << data2;
    info2 = tmp.str();
    tmp.clear();

    return info1 == info2;
}

bool compareParameterString(const ofec::ParameterVariant& data1, const ofec::ParameterVariant& data2) {
    if (data1.index() != data2.index())return false;
    auto type = static_cast<ofec::ParameterType>(data1.index());
    if (type == ofec::ParameterType::kReal) {
        return compareTwoRealVariant(std::get<ofec::Real>(data1), std::get<ofec::Real>(data2));
    }
    else if (type == ofec::ParameterType::kVectorReal) {
        auto& vdata1 = std::get<std::vector<ofec::Real>>(data1);
        auto& vdata2 = std::get<std::vector<ofec::Real>>(data2);
        if (vdata1.size() != vdata2.size())return false;
        for (int idx(0); idx < vdata1.size(); ++idx) {
            if (!compareTwoRealVariant(vdata1[idx], vdata2[idx])) {
                return false;
            }
        }
        return true;

    }
    else {
        return data1 == data2;
    }
}

void compare(ofec::ParameterVariantStream& paramsStream1, ofec::ParameterVariantStream& paramsStream2) {
    std::vector<ofec::ParameterVariantStream::Type> vectorParam1;
    std::vector<ofec::ParameterVariantStream::Type> vectorParam2;
    while (!paramsStream1.empty()) {
        vectorParam1.emplace_back(paramsStream1.popVariant());
    }
    while (!paramsStream2.empty()) {
        vectorParam2.emplace_back(paramsStream2.popVariant());
    }

    auto size = std::min(vectorParam1.size(), vectorParam2.size());
    for (int idx(0); idx < size; ++idx) {
        if (idx == 10) {
            int stop = -1;
        }
        if (!compareParameterString(vectorParam1[idx].first, vectorParam2[idx].first)) {
            //std::cout << "from\t" << ofec::variants_stream::m_errorInfos[idx].from;
            //std::cout << "\tto\t" << ofec::variants_stream::m_errorInfos[idx].to;
            //std::cout << "\tThreadId\t" << ofec::variants_stream::m_errorInfos[idx].InfoIdx << std::endl;

            std::cout << "error idx\t" << idx << "\tparameters1 \t";
            ofec::variants_stream::output(std::cout, vectorParam1[idx].first);
            std::cout << std::endl;
            // << vectorParam1[idx].first << std::endl;

            std::cout << "error idx\t" << idx << "\tparameters2 \t";
            //<< vectorParam2[idx].first << std::endl;

            ofec::variants_stream::output(std::cout, vectorParam2[idx].first);
            std::cout << std::endl;

            //   break;
        }
    }
}



void testBufReader(std::stringstream& buf, int from, int to) {
    using namespace ofec;
    ofec::variants_stream::ThreadLocalInfo info;
    ofec::variants_stream::GlobalThreadInfo globalInfo;
    globalInfo.m_info = buf.str();
    info.from = from;
    info.to = to;
    variants_stream::str2DataOneTaskLine(globalInfo, info);
}


void filterUniqueSols(std::vector<std::vector<std::vector<std::shared_ptr<ofec::SolutionBase>>>>& sol,
    std::vector<std::vector<std::vector<int>>>& solIds,
    std::vector<ofec::SolutionBase*>& solbases,
    ofec::TravellingSalesman* tsp_pro) {
    solIds.resize(sol.size());
    for (int idx(0); idx < solIds.size(); ++idx) {
        solIds[idx].resize(sol[idx].size());
        for (int idy(0); idy < solIds[idx].size(); ++idy) {
            solIds[idx][idy].resize(sol[idx][idy].size());
        }
    }
    //int totalSols = 0;
    std::map<unsigned long long, size_t> solHashToId;

    for (int idx(0); idx < sol.size(); ++idx) {
        for (int idy(0); idy < sol[idx].size(); ++idy) {
            for (int idz(0); idz < sol[idx][idy].size(); ++idz) {
                auto hash = tsp_pro->calHash(*sol[idx][idy][idz]);
                if (solHashToId.find(hash) == solHashToId.end()) {

                    solHashToId[hash] = solbases.size();
                    solbases.push_back(sol[idx][idy][idz].get());

                }
                solIds[idx][idy][idz] = solHashToId[hash];
            }
        }
    }



}

void calTaskMultiThreadInputOutputTest(const std::string& dirpath,
    const std::string& tspfilename,
    std::ostream& error_out, int numRun = 30) {
    using namespace ofec;
    using namespace std;
    using namespace chrono;
    auto filename = tspfilename.substr(0, tspfilename.find_first_of("."));

    std::cout << "calculating filename\t" << filename << std::endl;

    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");
    std::string proname = "TSP";


    std::shared_ptr<Environment> env(Environment::create());
    env->recordInputParameters();
    env->initialize();

    env->setProblem(ofec::Factory<ofec::Problem>::produce(proname));
    env->problem()->inputParameters().input(params);
    env->problem()->recordInputParameters();
    env->initializeProblem(0.5);

    auto pro = env->problem();


    //std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce(proname));
    //pro->inputParameters().input(params);
    //pro->recordInputParameters();
    ////  pro->initialize(0.5);


    //auto env = std::make_shared<ofec::Environment>();
    //env->initialize();
    //env->setProblem(pro);
    //
    //env->initializeProblem(0.5);


    if (pro->numberVariables() > 2000)return;
    auto rnd = std::make_shared<ofec::Random>(0.5);
    std::string algName = "EAX-TSP-sampling";


    auto sols = ofec::SamplingMutithread::runAlgMultiTask(0.5, 0.5, proname, params, algName, params, 1);
    // auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, numRun);


     //std::vector<int> popSols;
    std::vector<ofec::SolutionBase*> solbases;
    auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
    std::vector<std::vector<std::vector<int>>> solIds;

    filterUniqueSols(sols, solIds, solbases, tsp_pro);
    std::cout << "total solutions\t" << solbases.size() << std::endl;
    {
        ofec::ParameterVariantStream paramsStream;
        paramsStream << solIds.size();
        for (auto& it : solIds) {
            paramsStream << it.size();
            for (auto& it2 : it) {
                paramsStream << it2;
            }
        }
        std::stringstream buf;
        ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "solIds_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }




    //{
    //    int totalSols = 0;
    //    std::set<unsigned long long> setSols;
    //    for (auto& it : sols) {
    //        for (auto& it2 : it) {
    //            for (auto& it3 : it2) {

    //                ++totalSols;
    //                auto hash = tsp_pro->getHash(*it3);
    //                if (setSols.find(hash) == setSols.end()) {
    //                    solbases.push_back(it3.get());
    //                    setSols.insert(hash);

    //                }
    //            }
    //        }
    //    }
    //}

    auto eval_fun =
        [](ofec::SolutionBase& sol, ofec::Environment* env) {
        using namespace ofec;
        sol.evaluate(env);
        ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
        sol.setFitness(pos * sol.objective(0));
        };
    UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);




    {
        ofec::ParameterVariantStream paramsStream;
        paramsStream << solIds.size();
        for (auto& it : solIds) {
            paramsStream << it.size();
            for (auto& it2 : it) {
                paramsStream << it2;
            }
        }
        std::stringstream buf;
        ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "solIds_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();
    }


    {
        ofec::ParameterVariantStream paramsStream;
        paramsStream << solbases.size();
        for (auto& it : solbases) {
            auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
            paramsStream << cursol.variable().vect();
        }
        std::stringstream buf;
        ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "solVariables_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();
    }


    {
        // for test copy solutions


        std::vector<std::unique_ptr<ofec::SolutionBase>> solbasesCopy;
        auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro);
        std::vector<std::vector<std::vector<int>>> solIdsCopy;
        {

            ofec::ParameterVariantStream paramsStream;

            std::stringstream buf;
            std::ifstream in(dirpath + "solIds_" + tspfilename + ".txt");
            if (!in) {
                std::cout << "not open file" << std::endl;
                std::cout << dirpath + "solIds_" + tspfilename + ".txt" << std::endl;
            }
            buf << in.rdbuf();
            in.close();
            ofec::variants_stream::stream2ParametersMutithreadLine(buf, paramsStream);

            size_t dataSize;
            paramsStream >> dataSize;
            solIdsCopy.resize(dataSize);
            for (auto& it : solIdsCopy) {
                paramsStream >> dataSize;
                it.resize(dataSize);
                for (auto& it2 : it) {
                    paramsStream >> it2;
                }
            }

        }
        {


            ofec::ParameterVariantStream paramsStream;

            std::stringstream buf;
            std::ifstream in(dirpath + "solVariables_" + tspfilename + ".txt");
            if (!in) {
                std::cout << "not open file" << std::endl;
                std::cout << dirpath + "solVariables_" + tspfilename + ".txt" << std::endl;
            }
            buf << in.rdbuf();
            in.close();
            ofec::variants_stream::stream2ParametersMutithreadLine(buf, paramsStream);


            size_t dataSize;
            paramsStream >> dataSize;
            solbasesCopy.resize(dataSize);
            for (auto& it : solbasesCopy) {
                it.reset(tsp_pro->createSolution());
                auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
                paramsStream >> cursol.variable().vect();
            }
        }

        // compare two set of solutions

        if (solIds != solIdsCopy) {
            error_out << "error at solIds input \t" << tspfilename << std::endl;
            std::cout << "error at solIds input " << std::endl;
            return;
        }

        if (solbasesCopy.size() != solbases.size()) {
            error_out << "error at solbases size input \t" << tspfilename << std::endl;
            std::cout << "error at solbases size input " << std::endl;
            return;
        }

        for (int idx(0); idx < solbases.size(); ++idx) {
            if (!tsp_pro->same(solbases[idx]->variableBase(), solbasesCopy[idx]->variableBase())) {
                error_out << "error at solbases variable compare input \t" << tspfilename << std::endl;
                std::cout << "error at solbases variable compare input " << std::endl;
                return;
            }
        }

    }






    nbn::HnswModel hnswModel;
    hnswModel.initialize(env.get(), rnd.get(), std::thread::hardware_concurrency());
    {
        std::cout << "begin hnsw" << std::endl;
        auto start = system_clock::now();
        std::vector<int> solIds;
        //	new_datum->m_hnswModel2.setNumberThreads(1);
        hnswModel.addDataMutliThread(solbases, solIds);
        auto end = system_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        std::cout << "hnsw calculate costs"
            << double(duration.count()) * microseconds::period::num / microseconds::period::den
            << "seconds" << endl;
    }
    //for (auto& it : solbases) {
    //	new_datum->m_hnswModel.addDataSingleThread(it);
    //}


    std::vector<int> belong;
    std::vector<double> dis2parent;
    std::vector<double> vFitness;

    for (auto& it : solbases) {
        vFitness.push_back(it->fitness());
    }
    {
        std::cout << "begin nbn" << std::endl;
        auto start = system_clock::now();
        nbn::HnswModelNBN model2;
        model2.copy(hnswModel);
        //model2.setNumberThreads(1);
        model2.updateFitness(vFitness);
        model2.calculateNBN(belong, dis2parent);
        auto end = system_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        std::cout << "nbn calculate costs"
            << double(duration.count()) * microseconds::period::num / microseconds::period::den
            << "seconds" << endl;

    }

    nbn::HnswModel hnswModel_copy;
    hnswModel_copy.initialize(env.get(), rnd.get(), std::thread::hardware_concurrency());
    hnswModel_copy.setSolutions(solbases);

    {
        ofec::ParameterVariantStream paramsStream;

        hnswModel.configsToParameters(paramsStream);

        std::stringstream buf;
        ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "hnsw_parameters_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }

    ofec::ParameterVariantStream paramsStream_out;
    {
        ofec::ParameterVariantStream paramsStream;

        hnswModel.datasToParameters(paramsStream);
        auto  paramsStream2 = paramsStream;
        //   hnswModel_copy.datasFromParameters(paramsStream2);
        paramsStream_out = paramsStream;
        std::cout << "parameter size" << paramsStream.size() << std::endl;
        std::stringstream buf;
        ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "hnsw_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }

    {
        ofec::ParameterVariantStream paramsStream;
        std::stringstream buf;
        std::ifstream in(dirpath + "hnsw_parameters_" + tspfilename + ".txt");
        if (!in) {
            std::cout << "not open file" << std::endl;
            std::cout << dirpath + "hnsw_parameters_" + tspfilename + ".txt" << std::endl;
        }
        buf << in.rdbuf();
        in.close();

        std::cout << "buf size\t" << buf.str().size() << std::endl;
        ofec::variants_stream::stream2ParametersMutithreadLine(buf, paramsStream);

        hnswModel_copy.configsfromParameters(paramsStream);

    }

    ofec::ParameterVariantStream paramsStream_in;
    {

        // ofec::variants_stream::m_out.clear();


        ofec::ParameterVariantStream paramsStream;
        std::stringstream buf;
        std::ifstream in(dirpath + "hnsw_data_" + tspfilename + ".txt");
        buf << in.rdbuf();
        in.close();

        //   testBufReader(buf);
        ofec::variants_stream::stream2ParametersMutithreadLine(buf, paramsStream);
        paramsStream_in = paramsStream;
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::cout << "parameter size" << paramsStream.size() << std::endl;


        //    compare(paramsStream_in, paramsStream_out);
        hnswModel_copy.datasFromParameters(paramsStream);



        // testBufReader(buf,0,14017);

    }
    // compare(paramsStream_in, paramsStream_out);
    // 
    //  compare each node in hnswModel;
    if (hnswModel.numberNodes() != hnswModel_copy.numberNodes()) {
        error_out << tspfilename << std::endl;
        std::cout << "error" << std::endl;
        return;
        //   int stop = -1;
    }
    for (int idx(0); idx < hnswModel.numberNodes(); ++idx) {
        if (hnswModel.node(idx) != hnswModel_copy.node(idx)) {
            error_out << tspfilename << std::endl;
            std::cout << "error" << std::endl;
            return;
            //     int stop = -1;
        }
    }


    //{

    //    ofec::ParameterVariantStream paramsStream;
    //    paramsStream << belong;
    //    std::stringstream buf;
    //    ofec::FileReader::data2StrMutithread(buf, paramsStream);
    //    std::ofstream out(dirpath + "nbn_data_" + tspfilename + ".txt");
    //    out << buf.rdbuf();
    //    out.close();

    //}

}

void getTask(const std::string& path, std::set<std::string>& setNames) {


    // std::set<std::string> setNames;
   //  std::string path = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        auto filepath = entry.path().filename().string();
        int index = filepath.find_first_of("_") + 1;
        auto filename = filepath.substr(index, filepath.size() - index);
        index = filename.find_first_of("_") + 1;
        filename = filename.substr(index, filename.size() - index);
        filename = filename.substr(0, filename.find_first_of("."));
        filename = filename + ".tsp";
        setNames.insert(filename);
    }

    for (auto& it : setNames) {
        //     std::cout << it << std::endl;
    }
}

void calTask(const std::string& saveDir, std::ostream& error_out) {


    std::vector<std::string> filenames;
    std::string path = ofec::g_working_directory + "/instance/problem/combination/travelling_salesman/";
    for (const auto& entry : std::filesystem::directory_iterator(path)) {
        auto filepath = entry.path().filename().string();
        int index = filepath.find_first_of(".");
        if (index != -1) {
            auto fix = filepath.substr(index, filepath.size() - index);
            if (fix == ".tsp") {
                // std::cout << filepath << std::endl;
                filenames.push_back(filepath);
            }
        }

    }


    std::queue<std::string> tasks;
    std::set<std::string> finishTasks;
    getTask(saveDir, finishTasks);
    finishTasks.insert("XMC10150.tsp");
    finishTasks.insert("XSC6880.tsp");
    finishTasks.insert("pm8079.tsp");
    //  finishTasks.insert("f13795.tsp");
    std::cout << "finish Tasks" << std::endl;
    for (auto& it : finishTasks) {
        std::cout << it << std::endl;
    }

    for (auto& it : filenames) {
        if (finishTasks.find(it) == finishTasks.end()) {
            tasks.push(it);
        }
    }

    std::cout << "running task" << std::endl;
    //   calTaskMultiThread(saveDir, "EIL101.tsp");
    while (!tasks.empty()) {
        auto curfilename = tasks.front();
        std::cout << curfilename << std::endl;
        tasks.pop();
        calTaskMultiThreadInputOutputTest(saveDir, curfilename, error_out);
    }

}


void readRue500filenames(std::string& filepath, std::vector<std::string>& filenames) {
    std::ifstream in(filepath);
    std::string tag, tspname;
    double errData;
    int a, b;
    while (in >> tag >> tag >> tspname >> errData >> a >> b) {
        filenames.push_back(tspname + ".tsp");
        std::cout << tspname << std::endl;
    }
    in.close();
}
void calTask2(const std::string& saveDir, std::vector<std::string>& filenames, std::ostream& error_out) {

    std::queue<std::string> tasks;
    std::set<std::string> finishTasks;
    getTask(saveDir, finishTasks);
    //  finishTasks.insert("f13795.tsp");
    std::cout << "finish Tasks" << std::endl;
    for (auto& it : finishTasks) {
        std::cout << it << std::endl;
    }

    for (auto& it : filenames) {
        if (finishTasks.find(it) == finishTasks.end()) {
            tasks.push(it);
        }
    }

    std::cout << "running task" << std::endl;
    //   calTaskMultiThread(saveDir, "EIL101.tsp");
    while (!tasks.empty()) {
        auto curfilename = tasks.front();
        std::cout << curfilename << std::endl;
        tasks.pop();
        calTaskMultiThreadInputOutputTest(saveDir, curfilename, error_out);
    }
}


void runEax(const std::string& tspfilename, double seed) {
    using namespace ofec;
    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");

 //   std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce("TSP"));
 //   pro->inputParameters().input(params);
 //   pro->recordInputParameters();
    // pro->initialize(0.5);
   // if (pro->numberVariables() > 1000)return;

    std::shared_ptr<Environment> env(Environment::create());
    //auto env = std::make_shared<ofec::Environment>();
    env->initialize();
    env->setProblem(ofec::Factory<ofec::Problem>::produce("TSP"));
    env->problem()->inputParameters().input(params);
    env->problem()->recordInputParameters();
    env->initializeProblem(0.5);
    env->setAlgorithm(Factory<Algorithm>::produce("EAX-TSP-sampling"));
    env->initializeAlgorithm(seed);

    if (env->problem()->numberVariables() > 2000)return;

    auto samplingAlg = CAST_SAMPLING_ALG(env->algorithm());
    //   std::cout << "first Random\t" << env->algorithm()->randomNumber() << std::endl;
    samplingAlg->samplingDuringRun(env.get());
    //   std::cout << "last Random\t" << env->algorithm()->randomNumber() << std::endl;

    auto& sols = samplingAlg->getSols();
    int totalSol = 0;
    for (auto& it : sols) {
        for (auto& it2 : it) {
            ++totalSol;
        }
    }

    std::cout << "sols size\t" << totalSol << std::endl;

}



bool testEAX_test(const std::string& tspfilename, double seed,
    std::vector<unsigned long long>& rndSeq, std::vector<std::string>& errName) {

    bool normalRet = true;
    using namespace ofec;
    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");
    //auto env = std::make_shared<ofec::Environment>();

    //std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce("TSP"));
    //pro->inputParameters().input(params);
    //pro->recordInputParameters();
    //if (pro->numberVariables() > 1000)return true;

    std::shared_ptr<Environment> env(Environment::create());
    env->initialize();
    env->setProblem(ofec::Factory<ofec::Problem>::produce("TSP"));
    env->inputParameters().input(params);
    env->recordInputParameters();
    env->initializeProblem(0.5);
    env->setAlgorithm(Factory<Algorithm>::produce("EAX-TSP2"));
    env->initializeAlgorithm(seed);

    auto eaxAlg = CAST_EAX_TSP2(env->algorithm());

    auto  newRnd = std::make_shared<ofec::Random>(seed);
    eaxAlg->setRandom(newRnd);
    eaxAlg->environment().randomSequence = rndSeq;

    std::cout << "first Random\t" << eaxAlg->randomNumber() << "\t random idx\t" << eaxAlg->getRandom()->uniform.index() << std::endl;


    int numRun = 200;
    //  try{
    while (numRun--) {


        ++eaxAlg->environment().fCurNumOfGen;
        eaxAlg->environment().selectForMating(eaxAlg->getRandom());

        for (int s = 0; s < eaxAlg->environment().Npop; ++s)
        {
            eaxAlg->environment().generateKids(s, eaxAlg->getRandom());
            //  this->selectForSurvival(s);
        }

        eaxAlg->environment().setAverageBest();

    }
    //  }
    //  catch (const ofec::Exception& a) {
    //      std::cout <<"exption Info" << a.what() << std::endl;
     //     normalRet = false;
    //  }
    std::cout << "last Random\t" << eaxAlg->randomNumber() << "\t random idx\t" << eaxAlg->getRandom()->uniform.index() << std::endl;

    rndSeq = eaxAlg->environment().randomSequence;

    errName = eaxAlg->environment().error_infos;
    return normalRet;
}



void generateRandomTest() {
    ofec::Random rnd(0.5);
    unsigned size = 40000;
    std::vector<int> vi(1000);
    while (size--) {
        //rnd.uniform.next();
        //rnd.uniform.nextNonStd<double>(0.0, rnd.uniform.next() * 1000.0);
        rnd.uniform.shuffle(vi.begin(), vi.end());
    }
    std::cout << "final seed\n" << std::setprecision(7) << rnd.uniform.nextNonStd<double>(0.0, rnd.uniform.next() * 1000.0) << std::endl;
}

void testEax() {
    using namespace ofec;
    int curTime = 0;
    std::vector<unsigned long long> rndSeq;
    std::vector<std::string> errInfos;

    Real rnda = 0.5;

    testEAX_test("E500-25.tsp", rnda, rndSeq, errInfos);

    std::cout << "running time\t" << curTime++ << std::endl;


    int test = 1000;
    while (test--) {
        if (!testEAX_test("E500-25.tsp", rnda, rndSeq, errInfos)) {

        }
    }
}


void otherTask() {
    //int test = 1000;
   //while (test--) {
   //    runEax("E500-25.tsp", 0.5);
   //    
   //}


   //while (true) {
   //    testEax();
   //    //generateRandomTest();
   //}

   //calTask(saveDir, err_out);
   //ofec::variants_stream::m_out.open("test.txt");
   //calTaskMultiThreadInputOutputTest(saveDir, "E500-25.tsp", std::cout,3);
   //ofec::variants_stream::m_out.close();
   // calTask(saveDir);

   //stringstream str;
   //std::string info;
   //std::getline(str, info);
   //testParametersStringInputOutput();
}

void calRUE500task() {
    ofec::g_working_directory = "E:/Diao_Yiya/code/OFEC_v2_data";
    ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
    ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
    // ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";


    std::string saveDir = ofec::g_working_directory + "/nbn_tsp_rue500_sols_onerun/";
    std::string errorInfoPath = ofec::g_working_directory + "/nbn_data/error_data.txt";

    std::filesystem::create_directories(saveDir);

    //std::string filepath = "//172.24.242.8/share/Student/2018/YiyaDiao/NBN_data_paper_com/eax_tsp/runningdata22/problem.txt";

    auto filepath = ofec::g_working_directory + "/nbn_task/" + "problem.txt";
    std::vector<std::string> filenames;
    readRue500filenames(filepath, filenames);
    std::ofstream err_out(errorInfoPath);
    calTask2(saveDir, filenames, err_out);
    err_out.close();

}

void test() {

    ofec::g_working_directory = "E:/Diao_Yiya/code/OFEC_v2_data";
    ofec::g_working_directory = "//172.24.207.203/share/2018/diaoyiya/ofec-data";
    std::string saveDir = ofec::g_working_directory + "/nbn_tsp_rue500_sols_onerunTest/";
    std::filesystem::create_directories(saveDir);
    calTaskMultiThreadInputOutputTest(saveDir, "3141.tsp", std::cout, 1);
}


namespace ofec {

    void registerParamAbbr() {}
    void customizeFileName() {}
    void run() {

        using namespace ofec;
        using namespace std;

        registerInstance();
        // test();
         // test();
         // calRUE500task();

    }


}