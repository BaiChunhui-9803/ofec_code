#include "../utility/catch.hpp"
#include "../run/custom_method.h"
#include"../instance/problem/combination/airport_planing/airport_planning_problem.h"

#include "../core/global.h"

#include<string>
#include<iostream>
#include <iomanip>      
#include <ctime>  
#include <regex>
#include <string>
using namespace ofec;
using namespace std;


void read_time(std::tm& t, const std::string & time_s, const std::string & format) {
    std::istringstream ss(time_s);
    ss.imbue(std::locale("de_DE.utf-8"));
    
    ss >> std::get_time(&t, format.c_str());
//    ss >> std::get_time(&t, "%d%b%y %h%m");
    if (ss.fail()) {
        std::cout << "Parse failed\n";
    }
    else {
        std::cout << std::put_time(&t, "%d-%m-%Y-%H:%M") << '\n';
    }
}

void try_get_date(const std::string& s)
{
    std::cout << "Parsing the date out of '" << s <<
        "' in the locale " << std::locale().name() << '\n';
    std::istringstream str(s);
    std::ios_base::iostate err = std::ios_base::goodbit;

    std::tm t;
    std::istreambuf_iterator<char> ret =
        std::use_facet<std::time_get<char>>(str.getloc()).get_date(
            std::istreambuf_iterator<char>(str),
            std::istreambuf_iterator<char>(),
            str, err, &t
        );
    str.setstate(err);
    if (str) {
        std::cout << "Day: " << t.tm_mday << ' '
            << "Month: " << t.tm_mon + 1 << ' '
            << "Year: " << t.tm_year + 1900 << '\n';
    }
    else {
        std::cout << "Parse failed. Unparsed string: ";
        std::copy(ret, std::istreambuf_iterator<char>(),
            std::ostreambuf_iterator<char>(std::cout));
        std::cout << '\n';
    }
}

#include<fstream>

void test_file() {
    const std::string m_filepath = "E:/Diao_Yiya/code/OFEC/instance/problem/combination/airport_planing/data";
    std::string memberFileName2 = m_filepath + "data2";
    ofstream test(memberFileName2);
    test << "yes" << std::endl;
    test.close();
}


void test_file(const std::string& filenames) {

    ofstream test(filenames);
    test << "yes" << std::endl;
    test.close();
}



void print(Problem *pro,Solution<VariableVector<std::vector<int>>>& best,int id) {
    const std::string m_filepath = "E:/Diao_Yiya/code/OFEC/instance/problem/combination/airport_planing/data/";

    {
        ofstream outSol(m_filepath + "solution" + std::to_string(id) + ".csv");
        GET_APP(pro)->printSolReal(outSol, best);
        GET_APP(pro)->printSolReal(std::cout, best);
        outSol.close();
    }
    {
        ofstream outObj(m_filepath + "objective" + std::to_string(id) + ".csv");
        GET_APP(pro)->printSolObj(outObj, best);
        GET_APP(pro)->printSolObj(std::cout, best);
        outObj.close();
    }
}
void testRandomSolutions() {


    Solution s1;
    ParameterMap params;

    params["problem name"] = std::string("ComOP_APP");
    params["numDim"] = "5";

    Random *rnd = ADD_RND(0.5);
    int id_param = std::make_shared<const ParameterMap>(params);
    Problem *pro = Problem::generateByFactory(id_param, 0.1);

    pro->initialize();
  
    Solution<VariableVector<std::vector<int>>> best;
    GET_APP(pro)->initSolutionRandomRandMem(best, rnd);
    GET_APP(pro)->calObj(best);
    int id = 0;
    print(pro, best, id);
    Solution<VariableVector<std::vector<int>>> cur;


    while (true) {
        GET_APP(pro)->initSolutionRandomRandMem(cur, rnd);
        GET_APP(pro)->calObj(cur);
        if (GET_APP(pro)->better2Obj(cur, best)) {
            best = cur;
            print(pro, best, ++id);
        }
    }

}


void testrandomDFSSolutions() {


    Solution s1;
    ParameterMap params;

    params["problem name"] = std::string("ComOP_APP");
    params["numDim"] = "5";

    Random *rnd = ADD_RND(0.5);
    int id_param = std::make_shared<const ParameterMap>(params);
    Problem *pro = Problem::generateByFactory(id_param, 0.1);

    pro->initialize();

    Solution<VariableVector<std::vector<int>>> best;
    GET_APP(pro)->initSolutionRandomNetwork(best, rnd);
    GET_APP(pro)->calObj(best);
    int id = 0;
    print(pro, best, id);
    Solution<VariableVector<std::vector<int>>> cur;


    while (true) {
        GET_APP(pro)->initSolutionRandomNetwork(cur, rnd);
        GET_APP(pro)->calObj(cur);
        if (GET_APP(pro)->better2Obj(cur, best)) {
            best = cur;
            print(pro, best, ++id);
        }
    }

}



void testNetwork() {


    Solution s1;
    ParameterMap params;

    params["problem name"] = std::string("ComOP_APP");
    params["numDim"] = "5";

    Random *rnd = ADD_RND(0.5);
    int id_param = std::make_shared<const ParameterMap>(params);
    Problem *pro = Problem::generateByFactory(id_param, 0.1);

    pro->initialize();
    
  //  GET_APP(pro)
}

TEST_CASE("function", "[time_algorithm]") {

    testrandomDFSSolutions();
  //  testRandomSolutions();


//
//    Solution s1;
//    ParameterMap params;
//
//    params["problem name"] = std::string("ComOP_APP");
//    params["numDim"] = "5";
//
//    Random *rnd = ADD_RND(0.5);
//    int id_param = std::make_shared<const ParameterMap>(params);
//    Problem *pro = Problem::generateByFactory(id_param, 0.1);
//
//    pro->initialize();
//    const std::string m_filepath = "E:/Diao_Yiya/code/OFEC/instance/problem/combination/airport_planing/data/";
////std::string memberFileName2 = m_filepath + "data2";
//    Solution<VariableVector<std::vector<int>>> cur;
//    GET_APP(pro)->initSolutionRandom2Mem(cur, rnd);
//    {
//        ofstream outSol(m_filepath + "sol.csv");
//        GET_APP(pro)->printSolReal(outSol, cur);
//        GET_APP(pro)->printSolReal(std::cout, cur);
//        outSol.close();
//    }
//    {
//        ofstream outObj(m_filepath + "obj.csv");
//        GET_APP(pro)->printSolObj(outObj, cur);
//        GET_APP(pro)->printSolObj(std::cout, cur);
//        outObj.close();
//    }




    //string input = "C1F1";
    //std::string comma = "[CF]";
    //auto res_str = APP::split(input, comma);
    //    for (auto str : res_str)
    //    std::cout << "str is: " << str << std::endl;


    //const std::string m_filepath = "E:/Diao_Yiya/code/OFEC/instance/problem/combination/airport_planing/data/";
    //std::string memberFileName2 = m_filepath + "data2";



    //std::string memberFileName = m_filepath + "Data A-Crew.csv";
    //fstream memberFile(memberFileName);
    //std::string input;
    //while (memberFile) {
    //    std::getline(memberFile, input);
    //    std::cout << input << std::endl;
    //    std::string comma = ",";
    //    auto res_str = APP::split(input, comma);
    //    for (auto str : res_str)
    //        std::cout << "str is: " << str << std::endl;
    //}
	//std::string long_str = "FA2,8/12/2021,10:10,PGX,8/12/2021,11:40,NKX,C1F1";
	//std::string comma = ",";
	//auto res_str = APP::split(long_str, comma);
	//for (auto str : res_str)
	//	std::cout << "str is: " << str << std::endl;

 //   struct std::tm now {};
 //   const std::string& format = "%d/%m/%Y-%H:%M";
 //   auto datetime = res_str[1] + "-" + res_str[2];
 //   std::cout<<datetime << std::endl;
 //   read_time(now, datetime, format);




 //   struct std::tm next {};
 //   datetime = res_str[4] + "-" + res_str[5];
 //   std::cout << datetime << std::endl;
 //   read_time(next, datetime, format);


 //   time_t now_t = std::mktime(&now);
 //   time_t next_t = std::mktime(&next);

 //   std::cout << "difftime\t" << std::difftime(now_t, next_t) << std::endl;
}



TEST_CASE("problem", "[airport_plan]") {
	
}
