#ifndef EAX_TSP_TENVIRONMENT_H
#define EAX_TSP_TENVIRONMENT_H

#include <fstream> 
#include <string>
#include <vector>

#include "kopt.h"
#include "cross.h"
#include "../../../../../core/environment/environment.h"

namespace ofec {
    namespace eax_tsp {


        class TEnvironment
        {
        protected:

            //std::ofstream m_output_popBest;
            //std::string m_popBest_filename;

            //std::ofstream m_output_pop;
            //std::string m_pop_filename;

            //std::ofstream m_output_popChildren;
            //std::string m_popChildren_filename;


            //std::ofstream m_output_popSol;
            //std::ofstream m_output_solAttributes;
            //std::vector<int> m_outputSol;
            //std::string m_filepath;
            //int m_generation = 0;
            //int m_RunId = 0;
          //  int m_solId = 0;

        public:
            TEnvironment();
            ~TEnvironment();

            int getNpop()const {
                return Npop;
            }
            const std::vector<TIndi>& getCurPop()const {
                return tCurPop;
            }
            std::vector<TIndi>& getCurPop() {
                return tCurPop;
            }

            void setPopulation();
            

            //  void initialize(Random *rnd);
              //void setFilePath(const std::string& filepath) {
              //    m_filepath = filepath;
              //}
              //void setRunId(int runId) {
              //    m_RunId = runId;
              //    m_generation = 0;
              //    m_solId = 0;
              //}

            void setMaxGeneration(int gen) {
                m_maxGeneration = gen;
            }
            void setMaxStagnationIter(int gen) {
                m_maxStagnationIter = gen;

            }
            //void updateLEDDFile();
            //void outputLEDDInfo();

            void initialize() {
                this->init();
                this->initPop();

                this->getEdgeFreq();
                this->setAverageBest();
            }


            void evolve() {
                ++fCurNumOfGen;
                this->selectForMating();
                for (int s = 0; s < Npop; ++s)
                {
                    this->generateKids(s);
                    //  this->selectForSurvival(s);
                }
                this->setAverageBest();


                //if (fCurNumOfGen % 50 == 0)
                //{
                //    printf("%d:\t%lf\t%lf\n", fCurNumOfGen, fBestValue, fAverageValue);
                //    // record time every 50 gens

                //}

            }


            void define(ofec::Environment *env, const std::shared_ptr<Random> &rnd);   /* Define the variables */
            void doIt();   /* Main procedure of the GA */
            void doItOrigin();/* Origin version of GA*/

            void init();   /* Initialization of the GA */
            bool terminationCondition(); /* Decide whether to proceed to next stage (or treminate the GA) */
            void setAverageBest();       /* Compute average and best tour lengths of the population */

            void initPop();        /* Create an initial population */
            void selectForMating();    /* Determine a set of pairs of parents at each generation */
            void generateKids(int s);  /* Generate offspring solutions from a selected pair of parents. Selection for survival is also performed here. */
            void getEdgeFreq();        /* Compute the frequency of the edges of the population */
            void calEdgeFreq(std::vector< std::vector<std::vector<int>>*> & links);
            void clearEdgeFreq();

            void printOn();  /* Display and write summary of results */
            //void writeBest(); /* Write the best tour */


            // set parementer

            void setNPop(int npop) {
                Npop = npop;
            }
            void setNch(int nch) {
                Nch = nch;
            }
            void setTMax(int ctmax) {
                tmax = ctmax;
            }
            void setTeminate(bool flag) {
                terminate = flag;
            }


            double getRandom() {
                return tCross->tRand->m_random->uniform.next();
            }

            double getfBestValue()const  {
                return fBestValue;
            }

         
        protected:
            
            std::shared_ptr<TEvaluator> fEvaluator; /* Distance of the edges */
            std::shared_ptr<TCross> tCross;         /* Eede assembly crossover */
            std::shared_ptr<TKopt> tKopt;           /* Local search with the 2-opt neighborhood */
            std::string fFileNameTSP;     /* File name of an TSP instance */
            double optimum;            /* best known optimum cost */
            int tmax;               /* maximum running time in seconds*/
            bool terminate;         /* if terminate immediately */

            int Npop; /* Number of population members (N_pop in the paper) */
            int Nch;  /* Number of offspring solutions (N_ch in the paper) */
            std::vector<TIndi> tCurPop; /* Current population members */
            TIndi tBest;    /* Best solution in the current population */
            TIndi gBest;    /* Best solution so far */

            double gBestValue = -1;

            int m_maxGeneration = -1;
            int m_maxStagnationIter = -1;
            int m_stagnation_iteration = 0;
            double m_curBestCost = 0;



            int fCurNumOfGen;         /* The current number of generations */
            long int fAccumurateNumCh; /* The accumulated number of offspring solutions */

            int fBestNumOfGen;             /* The number of generations at which the current best solution was found */
            long int fBestAccumeratedNumCh; /* The accumulated number of offspring solutions at which the current best solution was found */
            std::vector<std::vector<int>> fEdgeFreq;                /* The frequency of the edges of the population */
            double fAverageValue;           /* The average tour lengths of the population */
            double fBestValue;                 /* The tour lenght of the best tour in the population */
            int fBestIndex;                 /* Index of the best population member */

            std::vector<int> fIndexForMating; /* Mating list (r[] in the paper) */
            int fStagBest;        /* The number of generations during which no improvement is found in the best tour */
            std::vector<int> fFlagC;       /* Specify configurations of EAX and selection strategy */
            int fStage;           /* Current stage */
            int fMaxStagBest;     /* If fStagBest = fMaxStagBest, proceed to the next stage */
            int fCurNumOfGen1;    /* Number of generations at which Stage I is terminated */

          //  clock_t fTimeStart, fTimeInit, fTimeEnd; /* Use them to measure the execution time */
        };

    }
}

#endif