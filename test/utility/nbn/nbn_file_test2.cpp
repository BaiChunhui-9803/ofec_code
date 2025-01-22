
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
#include "../core/parameter/file_reader.h"
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


void calTaskMultiThread(const std::string& dirpath, const std::string& tspfilename) {
    using namespace std;
    using namespace chrono;
    auto filename = tspfilename.substr(0, tspfilename.find_first_of("."));

    std::cout << "calculating filename\t" << filename << std::endl;

    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");

    std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce("TSP"));
    pro->getInputParameters().input(params);
    pro->recordInputParameters();
    pro->initialize(0.5);

    auto env = std::make_shared<ofec::Environment>();
    env->initialize();
    env->setProblem(pro);

    auto rnd = std::make_shared<ofec::Random>(0.5);
    std::string algName = "EAX-TSP-sampling";
    auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, 30);


    auto eval_fun =
        [](ofec::SolutionBase& sol, ofec::Environment* env) {
        using namespace ofec;
        sol.evaluate(env);
        ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
        sol.setFitness(pos * sol.objective(0));
        };
    std::vector<ofec::SolutionBase*> solbases;
    auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro.get());
    std::set<unsigned long long> setSols;
    for (auto& it : sols) {
        for (auto& it2 : it) {
            for (auto& it3 : it2) {
                auto hash = tsp_pro->getHash(*it3);
                if (setSols.find(hash) == setSols.end()) {
                    solbases.push_back(it3.get());
                    setSols.insert(hash);
                }
            }
        }
    }


    ////  solbases.resize(1000);
    //for (int idx(0); idx < solbases.size(); ++idx) {
    //    solbases[idx]->solId() = idx;
    //}

    UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);

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



    {
        ofec::ParameterVariantStream paramsStream;

        hnswModel.configsToParameters(paramsStream);

        std::stringstream buf;
        ofec::FileReader::data2StrMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "hnsw_parameters_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }


    {
        ofec::ParameterVariantStream paramsStream;

        hnswModel.datasToParameters(paramsStream);

        std::stringstream buf;
        ofec::FileReader::data2StrMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "hnsw_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }


    {

        ofec::ParameterVariantStream paramsStream;
        paramsStream << belong;
        std::stringstream buf;
        ofec::FileReader::data2StrMutithread(buf, paramsStream);
        std::ofstream out(dirpath + "nbn_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }

}



void calTask(const std::string& dirpath, const std::string& tspfilename) {
    using namespace std;
    using namespace chrono;
    auto filename = tspfilename.substr(0, tspfilename.find_first_of("."));

    std::cout << "calculating filename\t" << filename << std::endl;

    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");


    std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce("TSP"));
    pro->getInputParameters().input(params);
    pro->recordInputParameters();
    pro->initialize(0.5);
    std::shared_ptr<ofec::Environment> env = std::make_shared<ofec::Environment>();
    env->initialize();
    env->setProblem(pro);


    auto rnd = std::make_shared<ofec::Random>(0.5);
    std::string algName = "EAX-TSP-sampling";
    auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, 30);
    auto eval_fun =
        [](ofec::SolutionBase& sol, ofec::Environment* env) {
        using namespace ofec;
        sol.evaluate(env);
        ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
        sol.setFitness(pos * sol.objective(0));
        };
    std::vector<ofec::SolutionBase*> solbases;
    auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro.get());
    std::set<unsigned long long> setSols;
    for (auto& it : sols) {
        for (auto& it2 : it) {
            for (auto& it3 : it2) {
                auto hash = tsp_pro->getHash(*it3);
                if (setSols.find(hash) == setSols.end()) {
                    solbases.push_back(it3.get());
                    setSols.insert(hash);
                }
            }
        }
    }

    UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);
    n2::HnswSolModel model;
    model.initialize(env.get(), rnd.get());

    {
        std::cout << "begin hnsw" << std::endl;
        auto start = system_clock::now();
        std::vector<int> solIds;
        model.addDataMutliThread(solbases, solIds);
        auto end = system_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        std::cout << "hnsw calculate costs"
            << double(duration.count()) * microseconds::period::num / microseconds::period::den
            << "seconds" << endl;

    }


    {
        ofec::ParameterVariantStream paramsStream;
        std::vector<std::vector<int>> neighbors;
        paramsStream << model.numberNodes();
        for (int idx(0); idx < model.numberNodes(); ++idx) {
            model.getTotalNeighbors(idx, neighbors);
            paramsStream << neighbors.size();
            for (auto& it : neighbors) {
                paramsStream << it;
            }
        }
        std::stringstream buf;
        ofec::FileReader::data2StrMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "hnsw_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }


    std::vector<int> belong;
    std::vector<double> dis2parent;
    std::vector<double> vFitness;

    for (auto& it : solbases) {
        vFitness.push_back(it->fitness());
    }

    std::vector<std::vector<int>> neighbors(solbases.size());
    for (int idx(0); idx < solbases.size(); ++idx) {
        auto& nei = model.getNeighbor(idx);
        auto& curnei = neighbors[idx];
        for (auto& neiId : nei) {
            curnei.push_back(neiId.GetNode()->GetId());
        }
    }

    {
        std::cout << "begin nbn" << std::endl;
        auto start = system_clock::now();
        ofec::NBN_NearestBetterCalculator::calculate(solbases, vFitness, neighbors,
            belong, dis2parent,
            pro.get(), rnd.get());
        auto end = system_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        std::cout << "nbn calculate costs"
            << double(duration.count()) * microseconds::period::num / microseconds::period::den
            << "seconds" << endl;


    }





    {

        ofec::ParameterVariantStream paramsStream;
        paramsStream << belong;
        std::stringstream buf;
        ofec::FileReader::data2StrMutithread(buf, paramsStream);
        std::ofstream out(dirpath + "nbn_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }

}


void calculateNBNwithNeighborFromData(const std::string& dirpath, const std::string& tspfilename) {
    using namespace std;
    using namespace chrono;
    auto filename = tspfilename.substr(0, tspfilename.find_first_of("."));

    std::cout << "calculating filename\t" << filename << std::endl;

    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");

    std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce("TSP"));
    pro->getInputParameters().input(params);
    pro->recordInputParameters();
    pro->initialize(0.5);
    ofec::Environment env;
    env.initialize();
    env.setProblem(pro);


    auto rnd = std::make_shared<ofec::Random>(0.5);
    std::string algName = "EAX-TSP-sampling";
    auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, 30);
    auto eval_fun =
        [](ofec::SolutionBase& sol, ofec::Environment* env) {
        using namespace ofec;
        sol.evaluate(env);
        ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
        sol.setFitness(pos * sol.objective(0));
        };
    std::vector<ofec::SolutionBase*> solbases;
    auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro.get());
    std::set<unsigned long long> setSols;
    for (auto& it : sols) {
        for (auto& it2 : it) {
            for (auto& it3 : it2) {
                auto hash = tsp_pro->getHash(*it3);
                if (setSols.find(hash) == setSols.end()) {
                    solbases.push_back(it3.get());
                }
            }
        }
    }

    UTILITY::evaluateRandomSolutionsMultiThreads(solbases, &env, eval_fun);
    //n2::HnswSolModel model;
    //model.initialize(pro.get(), rnd.get());

    //{
    //    std::cout << "begin hnsw" << std::endl;
    //    auto start = system_clock::now();
    //    std::vector<int> solIds;
    //    model.addDataMutliThread(solbases, solIds);
    //    auto end = system_clock::now();
    //    auto duration = duration_cast<microseconds>(end - start);
    //    std::cout << "hnsw calculate costs"
    //        << double(duration.count()) * microseconds::period::num / microseconds::period::den
    //        << "seconds" << endl;

    //}


    std::vector<std::vector<std::vector<int>>> model_neighbors;
    {
        ofec::ParameterVariantStream paramsStream;
        std::stringstream buf;
        std::ifstream in(dirpath + "hnsw_data_" + tspfilename + ".txt");
        buf << in.rdbuf();
        in.close();

        ofec::FileReader::stream2ParametersMutithread(buf, paramsStream);

        int numNode(0);
        size_t numSize;
        paramsStream >> numNode;
        model_neighbors.resize(numNode);
        for (int idx(0); idx < numNode; ++idx) {
            paramsStream >> numSize;
            auto& curNeighbor = model_neighbors[idx];
            curNeighbor.resize(numSize);
            for (auto& it : curNeighbor) {
                paramsStream >> it;
            }
        }

    }


    std::vector<int> belong;
    std::vector<double> dis2parent;
    std::vector<double> vFitness;

    for (auto& it : solbases) {
        vFitness.push_back(it->fitness());
    }

    std::vector<std::vector<int>> neighbors(solbases.size());
    for (int idx(0); idx < solbases.size(); ++idx) {
        auto& nei = model_neighbors[idx][0];
        auto& curnei = neighbors[idx];
        curnei = nei;
    }

    {
        std::cout << "begin nbn" << std::endl;
        auto start = system_clock::now();
        ofec::NBN_NearestBetterCalculator::calculate(solbases, vFitness, neighbors,
            belong, dis2parent,
            pro.get(), rnd.get());
        auto end = system_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        std::cout << "nbn calculate costs"
            << double(duration.count()) * microseconds::period::num / microseconds::period::den
            << "seconds" << endl;


    }

    {

        ofec::ParameterVariantStream paramsStream;
        paramsStream << belong;
        std::stringstream buf;
        ofec::FileReader::data2StrMutithread(buf, paramsStream);
        std::ofstream out(dirpath + "nbn_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }

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
        //if (filename.front() == 'e') {
        //    int stop = -1;
        //}
        filename = filename + ".tsp";
        setNames.insert(filename);

        //if (index != -1) {
        //    auto fix = filepath.substr(index, filepath.size() - index);
        //    if (fix == ".tsp") {
        //        std::cout << filepath << std::endl;

        //        filenames.push(filepath);
        //    }
        //}

    }

    for (auto& it : setNames) {
        //     std::cout << it << std::endl;
    }
}
void calTask(const std::string& saveDir) {


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
        calTaskMultiThread(saveDir, curfilename);
    }

}


void calTask_hnsw_nbn(const std::string& saveDir) {
    std::vector<std::string> filenames;

    calTask(saveDir);
    // auto path = "//172.24.207.203/share/2018/diaoyiya/ofec-data/nbn_data/";
    // getTask(path, filenames);
}

void testBufReader(std::stringstream& buf) {
    using namespace ofec;
    ofec::variants_stream::ThreadLocalInfo info;
    ofec::variants_stream::GlobalThreadInfo globalInfo;
    globalInfo.m_info = buf.str();
    info.from = 0;
    info.to = 14191;

    variants_stream::str2DataOneTask(globalInfo, info);
}

// input HNSW_NBN from data

void input() {

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
        if (vectorParam1[idx] != vectorParam2[idx]) {
            std::cout << "from\t" << ofec::variants_stream::m_errorInfos[idx].from;
            std::cout << "\tto\t" << ofec::variants_stream::m_errorInfos[idx].to;
            std::cout << "\tThreadId\t" << ofec::variants_stream::m_errorInfos[idx].InfoIdx << std::endl;

            std::cout << "error idx\t" << idx << "\tparameters1 \t";
            ofec::output(std::cout, vectorParam1[idx].first);
            std::cout << std::endl;
            // << vectorParam1[idx].first << std::endl;

            std::cout << "error idx\t" << idx << "\tparameters2 \t";
            //<< vectorParam2[idx].first << std::endl;

            ofec::output(std::cout, vectorParam2[idx].first);
            std::cout << std::endl;

            break;
        }
    }


}

void calTaskMultiThreadInputOutputTest(const std::string& dirpath, const std::string& tspfilename) {
    using namespace std;
    using namespace chrono;
    auto filename = tspfilename.substr(0, tspfilename.find_first_of("."));

    std::cout << "calculating filename\t" << filename << std::endl;

    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");


    std::shared_ptr<ofec::Problem> pro(ofec::Factory<ofec::Problem>::produce("TSP"));
    pro->getInputParameters().input(params);
    pro->recordInputParameters();
    pro->initialize(0.5);

    auto env = std::make_shared<ofec::Environment>();
    env->initialize();
    env->setProblem(pro);


    auto rnd = std::make_shared<ofec::Random>(0.5);
    std::string algName = "EAX-TSP-sampling";
    auto sols = ofec::SamplingMutithread::runAlgMultiTask(pro.get(), algName, params, 1);


    auto eval_fun =
        [](ofec::SolutionBase& sol, ofec::Environment* env) {
        using namespace ofec;
        sol.evaluate(env);
        ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
        sol.setFitness(pos * sol.objective(0));
        };
    std::vector<ofec::SolutionBase*> solbases;
    auto tsp_pro = dynamic_cast<ofec::TravellingSalesman*>(pro.get());
    std::set<unsigned long long> setSols;
    for (auto& it : sols) {
        for (auto& it2 : it) {
            for (auto& it3 : it2) {
                auto hash = tsp_pro->getHash(*it3);
                if (setSols.find(hash) == setSols.end()) {
                    solbases.push_back(it3.get());
                    setSols.insert(hash);
                }
            }
        }
    }



    UTILITY::evaluateRandomSolutionsMultiThreads(solbases, env.get(), eval_fun);

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



    //{
    //    ofec::ParameterVariantStream paramsStream;

    //    hnswModel.configsToParameters(paramsStream);

    //    std::stringstream buf;
    //    ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
    //    std::cout << "buf size\t" << buf.str().size() << std::endl;
    //    std::ofstream out(dirpath + "hnsw_parameters_" + tspfilename + ".txt");
    //    out << buf.rdbuf();
    //    out.close();

    //}

    ofec::ParameterVariantStream paramsStream_out;
    {
        ofec::ParameterVariantStream paramsStream;

        hnswModel.datasToParameters(paramsStream);
        paramsStream_out = paramsStream;
        std::cout << "parameter size" << paramsStream.size() << std::endl;
        std::stringstream buf;
        ofec::variants_stream::parameters2StreamMutithread(buf, paramsStream);
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::ofstream out(dirpath + "hnsw_data_" + tspfilename + ".txt");
        out << buf.rdbuf();
        out.close();

    }

    nbn::HnswModel hnswModel_copy;
    hnswModel_copy.initialize(env.get(), rnd.get(), std::thread::hardware_concurrency());

    //{
    //    ofec::ParameterVariantStream paramsStream;
    //    std::stringstream buf;
    //    std::ifstream in(dirpath + "hnsw_parameters_" + tspfilename + ".txt");
    //    if (!in) {
    //        std::cout << "not open file" << std::endl;
    //        std::cout << dirpath + "hnsw_parameters_" + tspfilename + ".txt" << std::endl;
    //    }
    //    buf << in.rdbuf();
    //    in.close();

    //    std::cout << "buf size\t" << buf.str().size() << std::endl;
    //    ofec::variants_stream::stream2ParametersMutithread(buf, paramsStream);

    //    hnswModel_copy.configsfromParameters(paramsStream);

    //}
    hnswModel_copy.setSolutions(solbases);
    ofec::ParameterVariantStream paramsStream_in;
    {

        ofec::variants_stream::m_out.clear();


        ofec::ParameterVariantStream paramsStream;
        std::stringstream buf;
        std::ifstream in(dirpath + "hnsw_data_" + tspfilename + ".txt");
        buf << in.rdbuf();
        in.close();

        //   testBufReader(buf);
        ofec::variants_stream::stream2ParametersMutithread(buf, paramsStream);
        paramsStream_in = paramsStream;
        std::cout << "buf size\t" << buf.str().size() << std::endl;
        std::cout << "parameter size" << paramsStream.size() << std::endl;
        //  hnswModel_copy.datasFromParameters(paramsStream);

    }


    compare(paramsStream_in, paramsStream_out);

    // compare each node in hnswModel;
    //if (hnswModel.numberNodes() != hnswModel_copy.numberNodes()) {

    //    std::cout << "error" << std::endl;
    //    int stop = -1;
    //}
    //for (int idx(0); idx < hnswModel.numberNodes(); ++idx) {
    //    if (hnswModel.node(idx) != hnswModel_copy.node(idx)) {
    //        std::cout << "error" << std::endl;

    //        int stop = -1;
    //    }
    //}


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

using namespace std;


// 采样字符集
static constexpr char CCH[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

// sz: 字符串的长度
// printable：是否可打印。如果用作key，可以使用不可打印的字符哟
string makeRandStr(int sz, ofec::Random* rnd, bool printable)
{
    string ret;
    ret.resize(sz);
    // std::mt19937 rng(std::random_device{}());
    for (int i = 0; i < sz; ++i)
    {
        auto rndNum = rnd->uniform.nextNonStd<size_t>(std::numeric_limits<size_t>::lowest(), std::numeric_limits<size_t>::max());
        if (printable)
        {

            uint32_t x = rndNum % (sizeof(CCH) - 1);
            ret[i] = CCH[x];
        }
        else
        {
            ret[i] = rndNum % 0xFF;
        }
    }

    return ret;
}

template<typename T>
void setValue(ofec::ParameterVariant& data, ofec::Random* rnd) {
    T value = rnd->uniform.nextNonStd<T>(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
    data = value;
}


// 函数模板
template <>
void setValue<bool>(ofec::ParameterVariant& data, ofec::Random* rnd) {
    bool value = rnd->uniform.nextNonStd<int>(0, 2);
    data = value;
}

void setString(ofec::ParameterVariant& data, ofec::Random* rnd, size_t maxSize = 100) {
    std::string value = makeRandStr(rnd->uniform.nextNonStd<size_t>(0, maxSize), rnd, false);
    data = value;
}

template<typename T>
void setVectorValue(ofec::ParameterVariant& data, ofec::Random* rnd, size_t maxSize = 100) {
    std::vector<T> values(rnd->uniform.nextNonStd<size_t>(0, maxSize));
    for (auto& value : values) {
        value = rnd->uniform.nextNonStd<T>(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max());
    }
    data = values;
}

void generateVariants(ofec::ParameterVariant& data, ofec::Random* rnd, size_t maxSize = 100) {
    using namespace ofec;
    int numType = 10;
    auto type = static_cast<ofec::ParameterType>(rnd->uniform.nextNonStd<int>(0, numType));
    switch (type)
    {

    case ParameterType::kBool: {
        setValue<bool>(data, rnd);
        break;
    }
    case ParameterType::kInt: {
        setValue<int>(data, rnd);
        break;
    }
    case ParameterType::kSizeT: {
        setValue<size_t>(data, rnd);
        //info << std::get<size_t>(data) << "\t";
        break;
    }
    case ParameterType::kChar: {
        setValue<char>(data, rnd);
        //info << std::get<char>(data) << "\t";
        break;
    }
    case ParameterType::kReal: {
        setValue<ofec::Real>(data, rnd);
        break;
    }

    case ParameterType::kVectorBool: {
        setVectorValue<bool>(data, rnd, maxSize);

        break;
    }

    case ParameterType::kVectorInt: {
        setVectorValue<int>(data, rnd, maxSize);
        break;
    }

    case ParameterType::kVectorSizeT: {
        setVectorValue<size_t>(data, rnd, maxSize);
        break;
    }

    case ParameterType::kString: {
        setString(data, rnd, maxSize);
        break;
    }

    case ParameterType::kVectorReal: {
        setVectorValue<Real>(data, rnd, maxSize);
        break;
    }
    default:
        break;
    }
}

void generateParameterString(ofec::Random* rnd, ofec::ParameterVariantStream& paramsStream) {
    //    ofec::ParameterVariantStream paramsStream;
    ofec::ParameterVariant data;
    int numData = 100;
    while (numData--) {
        generateVariants(data, rnd);
        paramsStream << data;
    }
}

void testParametersStringInputOutput() {
    auto rnd = std::make_shared<ofec::Random>(0.5);
    ofec::ParameterVariantStream paramsStreamOut;
    ofec::ParameterVariantStream paramsStreamIn;
    generateParameterString(rnd.get(), paramsStreamOut);
    std::stringstream buf;
    ofec::variants_stream::parameters2StreamMutithread(buf, paramsStreamOut);
    ofec::variants_stream::stream2ParametersMutithreadLine(buf, paramsStreamIn);





    //{

    //    std::ofstream out("testParamters.txt");
    //    out << buf.rdbuf();
    //    out.close();
    //}
}





#include<filesystem>
#include<queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"

namespace ofec {

    void registerParamAbbr() {}
    void customizeFileName() {}
    void run() {

        using namespace ofec;
        using namespace std;

        ofec::g_working_directory = "E:/Diao_Yiya/code/OFEC_v2_data";
        // ofec::g_working_directory = "/home/lab408/share/2018/diaoyiya/ofec-data";
        // ofec::g_working_directory = "/mnt/Data/Student/2018/YiyaDiao/code_total/data";


        std::string saveDir = ofec::g_working_directory + "/nbn_data/";
        //    saveDir = "//172.24.207.203/share/2018/diaoyiya/ofec-data/nbn_data/";

            //std::filesystem::create_directories(saveDir);

        registerInstance();

        //ofec::variants_stream::m_out.open("test.txt");
        //calTaskMultiThreadInputOutputTest(saveDir, "A280.tsp");
        //ofec::variants_stream::m_out.close();
        // calTask(saveDir);

        //stringstream str;
        //std::string info;
        //std::getline(str, info);
        testParametersStringInputOutput();

    }


}