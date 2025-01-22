#ifndef EAX_TSP_TKOPT_H
#define EAX_TSP_TKOPT_H

#include "indi.h"
#include "evaluator.h"
#include "randomize.h"
#include <memory>

namespace ofec {
    namespace eax_tsp {
        class TKopt {
        public:
            TKopt(int N);
            ~TKopt();
            void setInvNearList();
            void transIndiToTree(TIndi& indi);
            void transTreeToIndi(TIndi& indi);
            void doIt(TIndi& tIndi); /* Apply a local search with the 2-opt neighborhood */

            int getNext(int t);
            int getPrev(int t);
            int turn(int& orient);

            void sub();
            void incrementImp(int flagRev);
            void combineSeg(int segL, int segS);

            void checkDetail();
            void checkValid();
            void swap(int& x, int& y);
            void makeRandSol(TIndi& indi); /* Set a random tour */
            void setSol(TIndi& indi, const std::vector<int>& sol);

            std::shared_ptr<TEvaluator> eval;
            std::shared_ptr<TRandom> tRand;
        private:
            int fN;
            int fFixNumOfSeg;
            int fNumOfSeg;
            int fFlagRev;
            double fTourLength;

            std::vector<std::vector<int>> fLink;
            std::vector<std::vector<int>> fLinkSeg;
            std::vector<std::vector<int>> fCitySeg;
            std::vector<std::vector<int>> fInvNearList;

            std::vector<int> fT;
            std::vector<int> fB;
            std::vector<int> fSegCity;
            std::vector<int> fOrdCity;
            std::vector<int> fOrdSeg;
            std::vector<int> fOrient;
            std::vector<int> fSizeSeg;
            std::vector<int> fActiveV;
            std::vector<int> fNumOfINL;
            std::vector<int> fArray;
            std::vector<int> fCheckN;
            std::vector<int> fGene;
        };

    }
}

#endif
