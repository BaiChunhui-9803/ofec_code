#ifndef EAX_TSP2_CROSS_H
#define EAX_TSP2_CROSS_H


#include <vector>
#include <array>
#include "indi.h"
#include "evaluator.h"

#include "randomize.h"
#include "sort.h"

#include <memory>

namespace ofec {
    namespace eax_tsp2 {

        class TEnvironment;

        class TCross
        {

            TEnvironment* m_tmp_environemt = nullptr;

        public:

            void bindEnvironment(TEnvironment* env);

            TCross(int N);
            //   ~TCross();
            void doIt(TIndi& tKid, TIndi& tPa2, int numOfKids, int flagP, 
                const std::vector<int>& flagC, 
                std::vector<std::vector<int>>& fEdgeFreq, Random *rnd);/* Main procedure of EAX */
            void doItWithoutParent(TIndi& tKid, int numOfKids, int flagP,
                const std::vector<int>& flagC, std::vector<std::vector<int>>& fEdgeFreq, Random* rnd);/* Main procedure of EAX */
            
            
            void setParents(const TIndi& tPa1, const TIndi& tPa2, const std::vector<int>& flagC,
                int numOfKids, Random* rnd);/* Set information of the parent tours */
            void setABcycle(const TIndi& parent1, const TIndi& parent2, 
                const std::vector<int>& flagC, int numOfKidsm, Random* rnd); /* Step 2 of EAX */

            void swap(int& x, int& y); /* Swap */
            void formABcycle();        /* Store an AB-cycle found */
            void changeSol(TIndi& tKid, int ABnum, int type); /* Apply an AB-cycle to an intermediate solution */

            void makeCompleteSol(TIndi& tKid, Random* rnd); /* Step 5 of EAX */
            void makeUnit(Random* rnd);                   /* Step 5-1 of EAX */
            void backToPa1(TIndi& tKid);       /* Undo the parent p_A */
            void goToBest(TIndi& tKid);        /* Modify tKid to the best offspring solution */

            void incrementEdgeFreq(std::vector<std::vector<int>>& fEdgeFreq); /* Increment fEdgeFreq[][] */
            int calAdpLoss(std::vector<std::vector<int>>& fEdgeFreq);         /* Compute the difference in the averate distance */
            double calEntLoss(std::vector<std::vector<int>>& fEdgeFreq);      /* Compute the difference in the edge entropy */

            /* Block2 */
            void setWeight(const TIndi& parent1, const TIndi& parent2);
            int calCNaive();
            void searchEset(int num, Random* rnd);
            void addAB(int num);
            void deleteAB(int num);


            //const std::vector<TIndi>& Localsearch_indis() {
            //    return m_localsearch_indis;
            //}

            int fNumOfGeneratedCh;
            std::shared_ptr<TEvaluator> eval;
            int Npop;

            std::shared_ptr<TRandom> tRand = nullptr;
            std::shared_ptr<TSort> tSort = nullptr;


        private:
            // tmp opt-2 
           //std::vector<std::array<int, 4>> m_opt2_circles;
            //std::vector<int> m_opt2_circles;
         //   std::vector<TIndi> m_localsearch_indis;


            TIndi tBestTmp;
            int fFlagImp;
            int fN;
            int r;
            int exam;
            int examFlag;
            int flagSt;
            int flagCycle;
            int prType;
            int chDis;
            int koritsuMany;
            int bunkiMany;
            int st;
            int ci;
            int pr;
            int stock;
            int stAppear;
            int fEvalType;
            int fEsetType;
            int fNumOfABcycleInESet;
            int fNumOfABcycle;
            int fPosiCurr;
            int fMaxNumOfABcycle;

            std::vector<int> koritsu;
            std::vector<int> bunki;
            std::vector<int> koriInv;
            std::vector<int> bunInv;
            std::vector<int> checkKoritsu;
            std::vector<int> fRoute;
            std::vector<int> fPermu;
            std::vector<int> fC;
            std::vector<int> fJun;
            std::vector<int> fOrd1;
            std::vector<int> fOrd2;

            std::vector<std::vector<int>> nearData;
            std::vector<std::vector<int>> fABcycle;

            // Speed Up Start
            int fNumOfUnit;
            int fNumOfSeg;
            //     int fNumOfSPL;
            int fNumOfElementInCU;
            int fNumOfSegForCenter;
            double fGainModi;
            int fNumOfModiEdge;
            int fNumOfBestModiEdge;
            int fNumOfAppliedCycle;
            int fNumOfBestAppliedCycle;

            std::vector<int> fOrder;
            std::vector<int> fInv;
            std::vector<int> fSegUnit;
            std::vector<int> fSegPosiList;
            std::vector<int> LinkAPosi;
            std::vector<int> fPosiSeg;
            std::vector<int> fNumOfElementInUnit;
            std::vector<int> fCenterUnit;
            std::vector<int> fListOfCenterUnit;
            std::vector<int> fSegForCenter;
            std::vector<double> fGainAB;
            std::vector<int> fAppliedCylce;
            std::vector<int> fBestAppliedCylce;

            std::vector<std::vector<int>> fSegment;
            std::vector<std::vector<int>> LinkBPosi;
            std::vector<std::vector<int>> fModiEdge;
            std::vector<std::vector<int>> fBestModiEdge;
            // Speed Up End

            // Block2
            int fNumOfUsedAB;
            int fNumC;
            int fNumE;
            int fTmax;
            int fMaxStag;
            int fNumOfABcycleInEset;
            int fDisAB;
            int fBestNumC;
            int fBestNumE;

            std::vector<int> fNumOfElementINAB;
            std::vector<int> fWeightSR;
            std::vector<int> fWeightC;
            std::vector<int> fUsedAB;
            std::vector<int> fMovedAB;
            std::vector<int> fABcycleInEset;

            std::vector<std::vector<int>> fInEffectNode;
            std::vector<std::vector<int>> fWeightRR;
        };

    }

}
#endif
