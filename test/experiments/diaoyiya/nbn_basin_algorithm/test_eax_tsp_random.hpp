
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
#include<filesystem>
#include<queue>
#include <set>

#include "../core/parameter/variants_to_stream.h"
#include "../instance/algorithm/combination/eax_tsp/eax_tsp_alg.h"

using namespace std;



bool testEAX_test(const std::string& tspfilename, double seed,
    std::vector<unsigned long long>& rndSeq, std::vector<std::string>& errName) {

    bool normalRet = true;
    using namespace ofec;
    ofec::ParameterMap params;
    //params["problem name"] = std::string("TSP");
    params["dataFile1"] = tspfilename;
    //params["algorithm name"] = std::string("EAX-TSP-sampling");
    auto env = std::make_shared<ofec::Environment>();


    //std::shared_ptr<ofec::Problem> pro();
    //pro->inputParameters().input(params);
    //pro->recordInputParameters();
    //if (pro->numberVariables() > 1000)return true;


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


void test() {

    ofec::g_working_directory = "E:/Diao_Yiya/code/OFEC_v2_data";
    std::string saveDir = ofec::g_working_directory + "/nbn_tsp_rue500_sols";
    //calTaskMultiThreadInputOutputTest(saveDir, "E500-25.tsp", std::cout, 3);
}



void copyTSPdata() {

}

namespace ofec {

    void registerParamAbbr() {}
    void customizeFileName() {}
    void run() {

        using namespace ofec;
        using namespace std;

        registerInstance();
        ofec::g_working_directory = "E:/Diao_Yiya/code/OFEC_v2_data";
        //  test();
        std::vector<unsigned long long> rndSeq;
        std::vector<std::string> errName;
        testEAX_test("A280.tsp", 0.01,
            rndSeq, errName);

        // EAX_TSP2 ALG;


        // calRUE500task();

    }


}