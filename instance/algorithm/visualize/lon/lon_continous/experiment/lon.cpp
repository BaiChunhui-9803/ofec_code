#include "lon.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <cstdlib>
#include <ratio>
#include <cmath>
#include <fstream>
#include "include/Tools.h"
#include <filesystem>
//#include <Python.h>
//#include "../python_c/call_python.h"
#include "../R_lan/callR.h"
//#include "../sampling/lon_con_sampling.h"

namespace ofec::lon {

    int Convergence(const std::vector<std::string>& argv) {

        using namespace std;
        // std::vector<std::string>  argv = { "0","cfg/SpreadSpectrumRadarPollyPhase.cfg", "10","10000" };

         //char* argv[] = { "0",
         // "rat575.tsp",
         //  "100",
         //  "30",
         //  "-1",
         //  "3600",
         //        "5"
         //};

        int argc = argv.size();


        string exe = "BasinHopping/BasinHoppingIter.py";
        vector < string > func;
        vector < double > best;
        vector < int > nvar, bounded;
        int opt_digits, hash_digits, runs, iter;
        vector< double > beta;
        vector< double > factor;
        string  init_mode, step_mode, name, output_dir;
        double step;

        string cmd, cmd_out, path;
        double cost;
        std::chrono::high_resolution_clock::time_point t1;
        std::chrono::high_resolution_clock::time_point t2;
        std::chrono::duration<double> time_span;
        ifstream fin;
        ofstream fout;
        ofstream fsta;
        int count;
        vector< long int > s;
        int RUNS = 10;
        int R = 10000;

        if (argc != 4) {
            std::cout << "./Convergence <cfg>  runs iter\n";
            std::cout << "./Convergence cfg/SpreadSpectrumRadarPollyPhase.cfg 10 10000\n";
            return 1;
        }
        name = argv[1];
        RUNS = std::stoi(argv[2]);
        R = std::stoi(argv[3]);

        cout << "RUNS " << RUNS << ", R " << R << "\n";

        loadConfigurationFile(name, func, nvar, bounded, best, step_mode, beta, factor,
            opt_digits, hash_digits, iter, runs,
            init_mode, output_dir);

        displayConfiguration(func, nvar, bounded, best, step_mode, beta, factor,
            opt_digits, hash_digits, iter, runs,
            init_mode, output_dir);

        cmd = "mkdir ./../Result";
        // system(cmd.c_str());
        std::cout << "cmd code\t" << cmd.c_str() << std::endl;
        std::string pathname = "./../Result";
        std::filesystem::create_directories(pathname);



        cmd = "mkdir " + output_dir;
        //  system(cmd.c_str());
        std::cout << "cmd code\t" << cmd.c_str() << std::endl;
        pathname = output_dir;
        std::filesystem::create_directories(pathname);



        path = output_dir + "/convergence/";
        cmd = "mkdir " + path;
        // system(cmd.c_str());
        std::cout << "cmd code\t" << cmd.c_str() << std::endl;
        pathname = path;
        std::filesystem::create_directories(pathname);




        // USE THE SAME SEED PER TRAINING INSTANCE ACROSS DIFFERENT STEPS
        s.resize(RUNS * name.size(), 0);
        s[0] = time(nullptr); // initial seed
        for (unsigned int i = 1; i < s.size(); i++)
            s[i] = s[i - 1] + 1;

        cmd = path + "/sta.txt";
        fsta.open(cmd, ios::out);
        count = 0;


        if (!fsta.is_open()) {
            cout << cmd << " not found !!!\n";
            return false;
        }

        for (unsigned int f = 0; f < func.size(); f++) {
            for (unsigned int i = 0; i < factor.size(); i++) {
                step = factor[i] * beta[f];
                cmd_out = path + "/data_" + func[f] + "n" + to_string(nvar[f]) + "_p" + to_string(step) + ".txt";
                fout.open(cmd_out, ios::out);

                for (int run = 1; run <= RUNS; run++) {
                    fsta << s[count] << " " << func[f] << " " << run << " ";

                    //! CREATE THE COMMAND
                    cmd = exe + " --seed " + to_string(s[count]) + " -i " + func[f];
                    cmd += " --nvar " + to_string(nvar[f]) + " --tstep " + step_mode;
                    cmd += " --step " + to_string(step) + " --iter " + to_string(R);
                    cmd += " --bounded " + to_string(bounded[f]);
                    cmd += " > out.txt";
                    std::cout << "cmd: " << cmd << std::endl;


                    ofec::ParameterMap par;
                    //par[]

                    


                    //std::string pythonFileName = "BasinHoppingIterFun";
                    //std::string pythonFunName = "BasinHoppingIterFun";

                    //std::vector<std::string> pars = { "python_path",
                    // "--seed " , to_string(s[count]) , "-i" , func[f],
                    // "--nvar" , to_string(nvar[f]) ,"--tstep" , step_mode,
                    // "--step" , to_string(step) , "--iter" , to_string(R),
                    // "--bounded" , to_string(bounded[f]) };

                    //PythonCaller::callFunction(pythonFileName, pythonFunName, pars);


                    //   std::vector<std::string> args;
                    //   args.push_back("--seed");
                    //   args.push_back(to_string(s[count]));
                    //  
                    //   



                    //   
                    //   



                    //   //! EXECUTE THE PROGRAM
                    //   t1 = high_resolution_clock::now();
                    // //  system(cmd.c_str());
                    //   
                    //   
                    //   {
                    //       PyObject* pModule = NULL;//声明变量
                    //       PyObject* pFunc = NULL;// 声明变量

                    //       pModule = PyImport_ImportModule(pythonFileName.c_str());//这里是要调用的文件名  pythonFileName.py
                    //       pFunc = PyObject_GetAttrString(pModule, pythonFunName.c_str());//这里是要调用的函数名

                    //       if (pModule == NULL)
                    //       {
                    //           std::cout << "没找到该Python文件:\t" << pythonFileName << std::endl;
                    //       }
                    //       else {


                    //          // PyObject* args = Py_BuildValue("(ss)", "a","--help");
                    //           //给python函数参数赋值

                    //           //paraments_type+= 
                    //           PyObject* args = Py_BuildValue("(sssssssssssssss)", "a",
                    //               "--seed ", to_string(s[count]),
                    //               "-i" , func[f],
                    //               "--nvar", to_string(nvar[f]), "--tstep" , step_mode,
                    //               "--step", to_string(step), "--iter" ,to_string(R),
                    //               "--bounded" , to_string(bounded[f])
                    //               );//给python函数参数赋值
                    //           PyObject* pRet = PyObject_CallObject(pFunc, args);//调用函数

                    //       }
                    //      
                    //   }
                    //   
                    //   t2 = high_resolution_clock::now();
                    //   time_span = duration_cast<duration< double >>(t2 - t1);
                    //   fsta << time_span.count() << endl;

                    //   //! READ THE COST
                    //   fin.open("out.txt", ios::in);
                    //   if (!fin.is_open()) {
                    //       cout << "out.txt" << " not found !!!\n";
                    //       return false;
                    //   }


                    //   while (!fin.eof()) {
                    //       fin >> cost;
                    //       if (fin.eof())
                    //           break;
                    //       fout << cost << " ";
                    //   }
                    //   fout << "\n";
                    //   fin.close();

                    // //  std::filesystem::remove("out.txt");
                    ////   system("rm -r out.txt");
                    //   count++;
                }
                fout.close();
            }
        }
        cout << " " << cmd_out << " is ready!!!\n";
        fsta.close();

        return 0;
    }



    int FilterData(const std::vector<std::string>& argv) {
        using namespace std;
        string input, output, report;
        ifstream fin;
        ofstream fout, frep;
        vector<int> run, fit, fit_next;
        vector<string> node, node_next;
        int aux_run, aux_fit, aux_fit_next;
        string aux_node, aux_node_next, head;
        vector<bool> flag_node, flag_next, ready_node, ready_next;
        bool dup;
        int min;
        int argc = argv.size();
        if (argc != 4) {
            cout << "./FilterData <input> <output> <report>\n";
            cout << "./FilterData data_hash/rastrigin_n5.txt data_clean/rastrigin_n5.txt report.txt\n";
            return 1;
        }
        input = argv[1];
        output = argv[2];
        report = argv[3];
        fin.open(input, ios::in);
        if (!fin.is_open()) {
            cout << "File: " << input << " not found!!!\n";
            return 1;
        }
        std::getline(fin, head);
        while (!fin.eof()) {
            fin >> aux_run >> aux_fit >> aux_node >> aux_fit_next >> aux_node_next;
            if (fin.eof())
                break;
            run.push_back(aux_run);
            fit.push_back(aux_fit);
            node.push_back(aux_node);
            fit_next.push_back(aux_fit_next);
            node_next.push_back(aux_node_next);
        }
        fin.close();
        cout << " " << input << " was loaded!!!\n";

        flag_node.resize(run.size(), false);
        flag_next.resize(run.size(), false);
        ready_node.resize(run.size(), false);
        ready_next.resize(run.size(), false);
        frep.open(report, ios::out);
        if (!frep.is_open()) {
            cout << "File: " << input << " not found!!!\n";
            return 1;
        }
        for (unsigned int i = 0; i < run.size() - 1; i++) {
            if (!ready_node[i]) {
                min = fit[i];
                dup = false;
                if (node[i] == node_next[i]) {
                    min = (min < fit_next[i]) ? min : fit_next[i];
                    flag_next[i] = true;
                    if (fit[i] != fit_next[i])
                        dup = true;
                }
                for (unsigned int j = i + 1; j < run.size(); j++) {
                    if (node[i] == node[j]) {
                        min = (min < fit[j]) ? min : fit[j];
                        flag_node[j] = true;
                        if (fit[i] != fit[j])
                            dup = true;
                    }
                    if (node[i] == node_next[j]) {
                        min = (min < fit_next[j]) ? min : fit_next[j];
                        flag_next[j] = true;
                        if (fit[i] != fit_next[j])
                            dup = true;
                    }
                }
                if (dup) {
                    frep << "node     " << node[i] << " is duplicated in line " << i + 2 << " current_fit: " << fit[i] << ", new fitness: " << min << endl;
                    fit[i] = min;
                    ready_node[i] = true;
                }
                if (flag_next[i]) {
                    flag_next[i] = false;
                    ready_next[i] = true;
                    if (dup) {
                        frep << "node_next " << node_next[i] << " is duplicated in line " << i + 2 << " current_fit: " << fit_next[i] << ", new fitness: " << min << endl;
                        fit_next[i] = min;
                    }
                }
                for (unsigned int j = i + 1; j < run.size(); j++) {
                    if (flag_node[j]) {
                        flag_node[j] = false;
                        ready_node[j] = true;
                        if (dup) {
                            frep << "node     " << node[j] << " is duplicated in line " << j + 2 << " current_fit: " << fit[j] << ", new fitness: " << min << endl;
                            fit[j] = min;
                        }
                    }
                    if (flag_next[j]) {
                        flag_next[j] = false;
                        ready_next[j] = true;
                        if (dup) {
                            frep << "node_next " << node_next[j] << " is duplicated in line " << j + 2 << " current_fit: " << fit_next[j] << ", new fitness: " << min << endl;
                            fit_next[j] = min;
                        }
                    }
                }
            }

            if (!ready_next[i]) {
                min = fit_next[i];
                dup = false;
                for (unsigned int j = i + 1; j < run.size(); j++) {
                    if (node_next[i] == node[j]) {
                        min = (min < fit[j]) ? min : fit[j];
                        flag_node[j] = true;
                        if (fit_next[i] != fit[j])
                            dup = true;
                    }
                    if (node_next[i] == node_next[j]) {
                        min = (min < fit_next[j]) ? min : fit_next[j];
                        flag_next[j] = true;
                        if (fit_next[i] != fit_next[j])
                            dup = true;
                    }
                }
                if (dup) {
                    frep << "node_next " << node_next[i] << " is duplicated in line " << i + 2 << " current_fit: " << fit_next[i] << ", new fitness: " << min << endl;
                    fit_next[i] = min;
                    ready_next[i] = true;
                }
                for (unsigned int j = i + 1; j < run.size(); j++) {
                    if (flag_node[j]) {
                        flag_node[j] = false;
                        ready_node[j] = true;
                        if (dup) {
                            frep << "node      " << node[j] << " is duplicated in line " << j + 2 << " current_fit: " << fit[j] << ", new fitness: " << min << endl;
                            fit[j] = min;
                        }
                    }
                    if (flag_next[j]) {
                        flag_next[j] = false;
                        ready_next[j] = true;
                        if (dup) {
                            frep << "node_next " << node_next[j] << " is duplicated in line " << j + 2 << " current_fit: " << fit_next[j] << ", new fitness: " << min << endl;
                            fit_next[j] = min;
                        }
                    }
                }
            }
        }
        frep.close();
        cout << " " << report << " is ready !!!\n";

        fout.open(output, ios::out);
        if (!fout.is_open()) {
            cout << "File: " << output << " not found!!!\n";
            return 1;
        }
        fout << head << endl;
        for (unsigned int i = 0; i < run.size(); i++) {
            fout << setfill(' ') << setw(4) << std::dec << run[i];
            fout << setw(11) << dec << fit[i] << "\t\t" << node[i];
            fout << setw(11) << dec << fit_next[i] << "\t\t" << node_next[i] << "\n";
        }

        fout.close();
        cout << " " << output << " is ready !!!\n";

        return 0;
    }


    

    int GenerateResults(const std::string& dir_path, const std::vector<std::string>& argv) {
        using namespace std;
        auto   CURRNET_DIR = dir_path;
            
          string program = "./../BasinHopping/BasinHoppingSampling.py";
        
          std::string pythonFileName = "BasicHoppingSamplingFun";
          std::string pythonFunName =  "BasicHoppingSamplingFun";
        
        
          // PARAMETERS
          vector < string > func;
          vector < double > best;
          vector < int > nvar, bounded;
          int opt_digits, hash_digits, runs, iter;
          vector< double > beta;
          vector< double > factor;
          string  init_mode, step_mode;
          double prec,pvalue;
          // OUTPUTS
          string cmd, name, ipath, opath, ofiltered, output_dir, table_name;
          int noptima, nfunnels, ngfunnels;
          double neutral, success, strength, deviation, sddeviation;
          ifstream fin;
          ofstream fout, fout2;
          int argc = argv.size();
          if( argc != 2 ) {
            std::cout << "Template: ./GenerateResults  <problem_config>\n";
            std::cout << "Example: ./GenerateResults cfg/SpreadSpectrumRadarPollyPhase.cfg\n";
            return 1;
          }
          name = argv[ 1 ];
          table_name = "table";
        
          loadConfigurationFile( name, func, nvar, bounded, best, step_mode, beta, factor,
                                 opt_digits, hash_digits, iter, runs,
                                  init_mode, output_dir );
        
          prec = pow( 10.0, -hash_digits );
        
          std::string resultPath =  std::string(CURRNET_DIR) + "/" + output_dir;
          std::filesystem::create_directories(resultPath);
          std::vector<std::string> paths = {"data_raw/","data_hash/","data_filtered/"};
        
          //for (auto& it : paths) {
          //    std::cout << "create path\t" << resultPath + it << std::endl;
          //    std::filesystem::create_directories(resultPath + it);
          //}
          {
              ipath = resultPath + "data_raw/";
              std::cout << "create path\t" << ipath << std::endl;
              std::filesystem::create_directories(ipath);
          }
          {
              opath = resultPath + "data_hash/";
              std::cout << "create path\t" << opath << std::endl;
              std::filesystem::create_directories(opath);
          }
        
          {
              ofiltered = resultPath + "data_filtered/";
              std::cout << "create path\t" << ofiltered << std::endl;
              std::filesystem::create_directories(ofiltered);
          }
        
        //
        //  cmd = "mkdir ./../Result";
        //  std::cout << cmd << std::endl;
        //  std::string pathname = "./../Result";
        //  std::filesystem::create_directories(pathname);
        //
        //  cmd = "mkdir " + output_dir;
        // // system(cmd.c_str());
        //
        //  std::cout << cmd << std::endl;
        //  pathname = output_dir;
        //  std::filesystem::create_directories(pathname);
        //
        //  ipath =  output_dir + "/data_raw/";
        //  cmd = "mkdir " + ipath;
        //  //system( cmd.c_str() );
        //  std::cout << cmd << std::endl;
        //  pathname = ipath;
        //  std::filesystem::create_directories(pathname);
        //
        //
        //  opath = output_dir + "/data_hash/";
        //  cmd = "mkdir " + opath;
        ////  system( cmd.c_str() );
        //  std::cout << cmd << std::endl;
        //  pathname = opath;
        //  std::filesystem::create_directories(pathname);
        //
        //  ofiltered = output_dir + "/data_filtered/";
        //  cmd = "mkdir " + ofiltered;
        //
        //  std::cout << cmd << std::endl;
        //  pathname = ofiltered;
        //  std::filesystem::create_directories(pathname);
        
          cmd = std::string(CURRNET_DIR) + "/" + output_dir + table_name + "_LON.txt";
          std::cout << "LON_TXT\t" << cmd << std::endl;
          fout.open( cmd, ios::out );
          if( !fout.is_open() ) {
            cout << cmd << "_LON.txt is not found!!!\n";
            return 1;
          }
        
          fout << setfill(' ') << setw(25) << "Instance" << "\t";
          fout << setfill(' ') << setw(10) << "Beta" << "\t";
          fout << setfill(' ') << setw(10) << "pstep" << "\t";
          fout << setfill(' ') << setw(10) << "step" << "\t";
          fout << setfill(' ') << setw(10) << "noptima" << "\t";
          fout << setfill(' ') << setw(10) << "nfunnels" << "\t";
          fout << setfill(' ') << setw(10) << "ngfunnels" << "\t";
          fout << setfill(' ') << setw(10) << "neutral" << "\t";
          fout << setfill(' ') << setw(10) << "strength" << "\t";
          fout << setfill(' ') << setw(10) << "success" << "\t";
          fout << setfill(' ') << setw(10) << "deviation" << "\t";
          fout << setfill(' ') << setw(10) << "sd_deviation" << "\n";
        
          std::flush(fout);
        
          cmd = std::string(CURRNET_DIR) + "/" +output_dir + "/" + table_name + "_CMLON.txt";
        
          std::cout << "CLON_TXT\t" << cmd << std::endl;
          fout2.open( cmd, ios::out );
          if( !fout2.is_open() ) {
            std::cout << cmd << "_CMLON.txt is not found!!!\n";
            return 1;
          }
        
          fout2 << setfill(' ') << setw(25) << "Instance" << "\t";
          fout2 << setfill(' ') << setw(10) << "Beta" << "\t";
          fout2 << setfill(' ') << setw(10) << "pstep" << "\t";
          fout2 << setfill(' ') << setw(10) << "step" << "\t";
          fout2 << setfill(' ') << setw(10) << "noptima" << "\t";
          fout2 << setfill(' ') << setw(10) << "nfunnels" << "\t";
          fout2 << setfill(' ') << setw(10) << "ngfunnels" << "\t"; 
          fout2 << setfill(' ') << setw(10) << "neutral" << "\t";
          fout2 << setfill(' ') << setw(10) << "strength" << "\t";
          fout2 << setfill(' ') << setw(10) << "success" << "\t";
          fout2 << setfill(' ') << setw(10) << "deviation" << "\t";
          fout2 << setfill(' ') << setw(10) << "sd_deviation" << "\n";
        
        
          std::flush(fout2);
          int s = time( nullptr ); // initial seed
        
        
          std::string outputFilepath = std::string(CURRNET_DIR) + "/" + output_dir + "output.txt";
        
          std::vector<std::string> pars;
          for( unsigned int f = 0; f < func.size(); f++ ) {
              for (unsigned int i = 0; i < factor.size(); i++) {
                  pvalue = beta[f] * factor[i];
        
                  name = "data_" + func[f] + "n" + to_string(nvar[f]) + "_p" + to_string(pvalue) + ".txt";
        
                  cmd = program + " --seed " + to_string(s) + " --ofile " + ipath + name;
                  cmd += " -i " + func[f] + " --nvar " + to_string(nvar[f]);
                  cmd += " --fopt " + to_string(best[f]);
                  cmd += " --tstep " + step_mode + " --step " + to_string(pvalue);
                  cmd += " --tinit " + init_mode + " --iter " + to_string(iter);
                  cmd += " --runs " + to_string(runs);
                  cmd += " --prec " + to_string(opt_digits);
                  cmd += " --bounded " + to_string(bounded[f]);
                  std::cout << "cmd: " << cmd << endl;
        
                  pars = { "python_path",
                      "--seed" ,to_string(s) , "--ofile" ,  ipath + name,
                      "-i" , func[f] ,"--nvar" ,to_string(nvar[f]),
                      "--fopt" ,to_string(best[f]),
                      "--tstep" , step_mode , "--step" ,to_string(pvalue),
                      "--tinit" , init_mode , "--iter" , to_string(iter),
                      "--runs" , to_string(runs),
                      "--prec" , to_string(opt_digits),
                      "--bounded" , to_string(bounded[f])
                  };

                 // ofec::getParValue(par, "problem name", ofile, ofile);
                 // ofec::getParValue(par, "number of variables", numVar, numVar);

                  ofec::ParameterMap par;
                  par["problem name"] = std::string("BBOB_F01");
                  par["number of variables"] = 2;

                  par["problem name"] = std::string("Classic_Ackley");
                  par["number of variables"] = 2;
                  par["--seed"] = double(0.5);
                  par["--ofile"] = ipath + name;
                  par["--iter"] = iter;
                  par["--runs"] = runs;
                  par["--prec"] = opt_digits;
        

            //      localOptimaNetworkConSampling(par);
            //      //hash_digits = 1;

            ////      PythonCaller::callFunction(pythonFileName, pythonFunName, pars);
            //      s++;
            //  //    //   system( cmd.c_str() );
        
        
            //         // Generate the hashing version
            //      cmd = "./Hashing " + ipath + name + " " + opath + name + " " + to_string(hash_digits);
            //      std::cout << "cmd: " << cmd << endl;
            //      pars = { "Hashing", ipath + name ,opath + name ,to_string(hash_digits) };
            //      Hashing(pars);
        
        
            //  //// system( cmd.c_str() );
        
            //    // Filter repeated nodes
            //    cmd = "./FilterData " + opath + name + " " + ofiltered + name + " " + ofiltered + "rep_" + name;
            //    std::cout << "cmd: " << cmd << endl;
            //    pars = { "FilterData", opath + name ,ofiltered + name ,ofiltered + "rep_" + name };
            //    FilterData(pars);
               // system( cmd.c_str() );
        
                // Compute the metrics
                //cmd = "./../Graph/GraphMetrics.r " +  ofiltered + name + " ";
                //cmd += to_string(hashing_value(best[f], hash_digits)) + " > out.txt";
                //std::cout << "cmd: " << cmd << endl;




                cmd = std::string(CURRNET_DIR)+"/Graph/EvalGraph.R " + ofiltered + name + " ";
                cmd += to_string(hashing_value(best[f], hash_digits));
                cmd += " " + outputFilepath;
        
                std::cout << "cmd\t" << cmd << std::endl;
  //              RCaller::runR( cmd.c_str() );
        
        
                
                cmd = std::string(CURRNET_DIR) + "/Graph/Data2GraphVizPic.r " + ofiltered + name + " ";
                cmd += to_string(hashing_value(best[f], hash_digits));
                cmd += " " + resultPath;

                std::cout << "cmd\t" << cmd << std::endl;
                RCaller::runR(cmd.c_str());
        
        
                fin.open( "out.txt", ios::in );
        
                fin >> cmd >> noptima;
                fin >> cmd >> nfunnels;
                fin >> cmd >> ngfunnels;
                fin >> cmd >> neutral;
                fin >> cmd >> strength;
                fin >> cmd >> success;
                fin >> cmd >> deviation;
                fin >> cmd >> sddeviation;
        
                fout << setfill(' ') << setw(25) << func[ f ] + "n" + to_string( nvar[ f ] ) << "\t";
                fout << setfill(' ') << setw(10) << beta[ f ] << "\t";
                fout << setfill(' ') << setw(10) << factor[ i ] << "\t";
                fout << setfill(' ') << setw(10) << pvalue << "\t";
                fout << setfill(' ') << setw(10) << noptima << "\t";
                fout << setfill(' ') << setw(10) << nfunnels << "\t";
                fout << setfill(' ') << setw(10) << ngfunnels << "\t";
                fout << setfill(' ') << setw(10) << neutral << "\t";
                fout << setfill(' ') << setw(10) << strength << "\t";
                fout << setfill(' ') << setw(10) << success << "\t";
                fout << setfill(' ') << setw(10) << deviation*prec << "\t";
                fout << setfill(' ') << setw(10) << sddeviation*prec << "\n";
        
                fin >> cmd >> noptima;
                fin >> cmd >> nfunnels;
                fin >> cmd >> ngfunnels;
                fin >> cmd >> neutral;
                fin >> cmd >> strength;
        
                fout2 << setfill(' ') << setw(25) << func[ f ] + "n" + to_string( nvar[ f ] ) << "\t";
                fout2 << setfill(' ') << setw(10) << beta[ f ] << "\t";
                fout2 << setfill(' ') << setw(10) << factor[ i ] << "\t";
                fout2 << setfill(' ') << setw(10) << pvalue << "\t";
                fout2 << setfill(' ') << setw(10) << noptima << "\t";
                fout2 << setfill(' ') << setw(10) << nfunnels << "\t";
                fout2 << setfill(' ') << setw(10) << ngfunnels << "\t";
                fout2 << setfill(' ') << setw(10) << neutral << "\t";
                fout2 << setfill(' ') << setw(10) << strength << "\t";
                fout2 << setfill(' ') << setw(10) << success << "\t";
                fout2 << setfill(' ') << setw(10) << deviation*prec << "\t";
                fout2 << setfill(' ') << setw(10) << sddeviation*prec << "\n";
        
                fin.close();
               // system( "rm -r out.txt " );
        
        
                std::flush(fout);
                std::flush(fout2);
            }
          }
          fout.close();
          fout2.close();
          cout << table_name << "_LON.txt is ready!!!\n";
          cout << table_name << "_CMLON.txt is ready!!!\n";

        return 0;
    }


    int GenerateResultsEigenSampling(const std::string& rpath, const std::string& dir_path, const ParameterMap& v)
    {

        LonSamplingPar par;
        par.setParameters(v);
        std::string output_dir = "proname_" + par.proname + "_dim_" + std::to_string(par.numVar);
        std::string name = output_dir + ".txt";

        std::string resultPath = dir_path + "/" + output_dir + "/";
        
            std::string ipath = resultPath + "data_raw/";
            std::cout << "create path\t" << ipath << std::endl;
            std::filesystem::create_directories(ipath);

            par.ofile = ipath + name;
        
        
            std::string opath = resultPath + "data_hash/";
            std::cout << "create path\t" << opath << std::endl;
            std::filesystem::create_directories(opath);
        

        
            std::string ofiltered = resultPath + "data_filtered/";
            std::cout << "create path\t" << ofiltered << std::endl;
            std::filesystem::create_directories(ofiltered);
        

            localOptimaNetworkConSampling(v,par);


      std::vector<std::string> pars;
      pars = { "Hashing", ipath + name ,opath + name ,std::to_string(par.hash_digit) };
      Hashing(pars);

      pars = { "FilterData", opath + name ,ofiltered + name ,ofiltered + "rep_" + name };
      FilterData(pars);



      resultPath = resultPath + "result.txt";
      std::string cmd = std::string(rpath) + "/Graph/EvalGraph.R " + ofiltered + name + " ";
      cmd += std::to_string(hashing_value(par.m_global_fitness, par.hash_digit));
      cmd += " " + resultPath;


      std::cout << "cmd\t" << cmd << std::endl;
      RCaller::runR(cmd.c_str());


      cmd = std::string(rpath) + "/Graph/Data2GraphViz.r " + ofiltered + name + " ";
      cmd += std::to_string(hashing_value(par.m_global_fitness, par.hash_digit));
      cmd += " " + resultPath;


      std::cout << "cmd\t" << cmd << std::endl;
      RCaller::runR(cmd.c_str());




      std::ifstream fin;
      fin.open(resultPath, std::ios::in);
      {

          int noptima, nfunnels, ngfunnels;
          double neutral, success, strength, deviation, sddeviation;


          fin >> cmd >> noptima;
          fin >> cmd >> nfunnels;
          fin >> cmd >> ngfunnels;
          fin >> cmd >> neutral;
          fin >> cmd >> strength;
          fin >> cmd >> success;
          fin >> cmd >> deviation;
          fin >> cmd >> sddeviation;



          fin >> cmd >> noptima;
          fin >> cmd >> nfunnels;
          fin >> cmd >> ngfunnels;
          fin >> cmd >> neutral;
          fin >> cmd >> strength;


          std::cout << noptima << "\t" << nfunnels << "\t" << ngfunnels << "\t" << neutral << "\t" << strength << std::endl;
      }
      fin.close(); 
      return 0;
    }



    int GenerateResultsEigenSamplingNode(
        const std::string& dir_path,  LonSamplingPar& par, LonStruct &lon) {
        using namespace ofec;

      //  par.RUN = 1;
        auto&  curparam = *par.m_param;
        std::string proname = curparam.get<std::string>("problem name");
        if (proname == "free_peaks") {

            //par["taskName"] = curcurInfo.m_taskname;
            proname = curparam.get<std::string>("taskName");

        }
        std::string output_dir = "proname_" + proname + "_dim_" + std::to_string(par.numVar) + "_numStep_" + std::to_string(int(par.m_step*100));
        std::string name = output_dir + ".txt";

        std::string resultPath = dir_path + "/" + output_dir + "/";

        std::string ipath = resultPath + "data_raw/";
    //    std::cout << "create path\t" << ipath << std::endl;
        std::filesystem::create_directories(ipath);
        par.ofile = ipath + name;
        std::string opath = resultPath + "data_hash/";
    //    std::cout << "create path\t" << opath << std::endl;
        std::filesystem::create_directories(opath);
        std::string ofiltered = resultPath + "data_filtered/";
    //    std::cout << "create path\t" << ofiltered << std::endl;
        std::filesystem::create_directories(ofiltered);


      // LonStruct lon;




        localOptimaNetworkConSamplingLON(par,lon);


 /*       std::vector<std::string> pars;
        pars = { "Hashing", ipath + name ,opath + name ,std::to_string(par.hash_digit) };
        Hashing(pars);

        pars = { "FilterData", opath + name ,ofiltered + name ,ofiltered + "rep_" + name };
        FilterData(pars);*/

        filterLON(lon, par);
      //  std::cout << "filepath\t" << ofiltered + name << std::endl;
        std::ofstream fout(ofiltered + name, std::ios::out);
        outputLON(fout,lon);
        fout.close();


 //       par.release();
        return 0;

    }


    void EvaluateLONbyR(const std::string& rpath, const std::string& dir_path, const LonSamplingPar& par) {

        std::string output_dir = "proname_" + par.proname + "_dim_" + std::to_string(par.numVar);
        std::string name = output_dir + ".txt";

        std::string resultPath = dir_path + "/" + output_dir + "/";

        std::string ipath = resultPath + "data_raw/";
        std::string opath = resultPath + "data_hash/";
        std::string ofiltered = resultPath + "data_filtered/";


        resultPath = resultPath + "result.txt";
      std::string cmd = std::string(rpath) + "/Graph/EvalGraph.R " + ofiltered + name + " ";

      
      cmd += std::to_string(hashing_value(par.m_global_fitness, par.hash_digit));
      cmd += " " + resultPath;


      std::cout << "cmd\t" << cmd << std::endl;
      RCaller::runR(cmd.c_str());


      cmd = std::string(rpath) + "/Graph/Data2GraphViz.r " + ofiltered + name + " ";
      cmd += std::to_string(hashing_value(par.m_global_fitness, par.hash_digit));
      cmd += " " + resultPath;


      std::cout << "cmd\t" << cmd << std::endl;
      RCaller::runR(cmd.c_str());




      std::ifstream fin;
      fin.open(resultPath, std::ios::in);
      {

          int noptima, nfunnels, ngfunnels;
          double neutral, success, strength, deviation, sddeviation;

          fin >> cmd >> noptima;
          fin >> cmd >> nfunnels;
          fin >> cmd >> ngfunnels;
          fin >> cmd >> neutral;
          fin >> cmd >> strength;
          fin >> cmd >> success;
          fin >> cmd >> deviation;
          fin >> cmd >> sddeviation;

          fin >> cmd >> noptima;
          fin >> cmd >> nfunnels;
          fin >> cmd >> ngfunnels;
          fin >> cmd >> neutral;
          fin >> cmd >> strength;


          std::cout << noptima << "\t" << nfunnels << "\t" << ngfunnels << "\t" << neutral << "\t" << strength << std::endl;
      }
      fin.close();
    }

    int Hashing(const std::vector<std::string>& argv) {
        using namespace std;
        string input, output, cmd;
        double pre, fit, temp;
        int D, run;
        ifstream fin;
        ofstream fout;
        int argc = argv.size();
        if (argc != 4) {
            cout << "./hashing <input> <output> <resolution>\n";
            cout << "./hashing SpreadSpectrumRadarPollyPhase_n5.txt SpreadSpectrumRadarPollyPhase_n5.txt 4 \n";
            return 1;
        }
        input = argv[1];
        output = argv[2];
        pre = std::stoi(argv[3]);
        fin.open(input, ios::in);
        fout.open(output, ios::out);
        if (!fin.is_open()) {
            cout << "File: " << input << " not found!!!\n";
            return 1;
        }
        fin >> cmd;
        fout << setfill(' ') << setw(4) << cmd << "\t";
        fin >> cmd;
        fout << setfill(' ') << setw(11) << "Fitness" << "\t";
        fin >> cmd >> D;
        fout << setfill(' ') << setw(8 * D) << "Solution" << "\t";
        fin >> cmd;
        fout << setfill(' ') << setw(11) << "Fitness" << "\t";
        fin >> cmd;
        fout << setfill(' ') << setw(8 * D) << "Next_Solution" << "\n";
        while (true) {
            fin >> run;
            if (fin.eof())
                break;
            fin >> fit;
            fout << setfill(' ') << setw(4) << std::dec << run << "\t" << setw(11) << dec << hashing_value(fit, pre) << "\t";
            for (int i = 0; i < D; i++) {
                fin >> temp;
                fout << std::setfill('0') << setw(8) << hex << hashing_value(temp, pre);
            }
            fin >> fit;
            fout << "\t" << setfill(' ') << setw(11) << dec << hashing_value(fit, pre) << "\t";
            for (int i = 0; i < D; i++) {
                fin >> temp;
                fout << setfill('0') << setw(8) << hex << hashing_value(temp, pre);
            }
            fout << "\n";
        }
        fin.close();
        fout.close();
        cout << " " << input << " was loaded!!!\n";
        cout << " " << output << " is ready!!!\n";

        return 0;
    }

}