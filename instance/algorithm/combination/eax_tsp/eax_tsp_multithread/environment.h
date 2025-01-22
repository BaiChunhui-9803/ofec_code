#ifndef EAX_TSP_MULTITHREAD_TENVIRONMENT_H
#define EAX_TSP_MULTITHREAD_TENVIRONMENT_H

#include <array>
#include <fstream> 
#include <string>
#include <vector>

#include "kopt.h"
#include "cross.h"
#include "../../../../../core/environment/environment.h"

namespace ofec {
    namespace eax_tsp_mt {


        class TEnvironment
        {
        public:
            TEnvironment() :fFlagC(10, -1) {}
            ~TEnvironment() = default;


            void define(ofec::Environment *env, const std::shared_ptr<Random> &rnd);   /* Define the variables */
           // void doIt();   /* Main procedure of the GA */
           // void doItOrigin();/* Origin version of GA*/

            void init();   /* Initialization of the GA */
            //bool terminationCondition(); /* Decide whether to proceed to next stage (or treminate the GA) */
            bool terminationCondition(const std::vector<TIndi*>& curPop);/* Decide whether to proceed to next stage (or treminate the GA) */
            void setAverageBest(const std::vector<TIndi*>& curPop);       /* Compute average and best tour lengths of the population */

           // void initPop();        /* Create an initial population */
          //  void selectForMating();    /* Determine a set of pairs of parents at each generation */
           // void generateKids(int s);  /* Generate offspring solutions from a selected pair of parents. Selection for survival is also performed here. */
           // void getEdgeFreq();        /* Compute the frequency of the edges of the population */
            void calEdgeFreq(std::vector< std::vector<std::vector<int>>*> & links);
            void clearEdgeFreq();


            // set parementer
            void setMaxGeneration(int gen) {
                m_maxGeneration = gen;
            }
            void setMaxStagnationIter(int gen) {
                m_maxStagnationIter = gen;
            }
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


            int getNpop()const {
                return Npop;
            }
         
        protected:
            std::shared_ptr<Random> m_random;
            std::shared_ptr<TEvaluator> fEvaluator; /* Distance of the edges */
            std::shared_ptr<TCross> tCross;         /* Eede assembly crossover */
            std::shared_ptr<TKopt> tKopt;           /* Local search with the 2-opt neighborhood */
            std::string fFileNameTSP;     /* File name of an TSP instance */

          //  int m_N = 0;
            double optimum;            /* best known optimum cost */
            int tmax;               /* maximum running time in seconds*/
            bool terminate;         /* if terminate immediately */

            int Npop; /* Number of population members (N_pop in the paper) */
            int Nch;  /* Number of offspring solutions (N_ch in the paper) */
           // std::vector<TIndi> tCurPop; /* Current population members */
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

          //  std::vector<int> fIndexForMating; /* Mating list (r[] in the paper) */
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