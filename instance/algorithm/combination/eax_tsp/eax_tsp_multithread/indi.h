#ifndef EAX_TSP_MULTITHREAD_INDI_H
#define EAX_TSP_MULTITHREAD_INDI_H

#include <vector>
#include "../../../../../core/problem/solution.h"

namespace ofec::eax_tsp_mt {

    class TIndi : public Solution<VariableVector<int>>
    {

    public:
        int fN = 0;     /* Number of cities */
        std::vector<std::vector<int>> fLink;
        // int **fLink;         /* fLink[i][]: two vertices adjacent to i */
        double fEvaluationValue; /* Tour length of */


        TIndi();
        ~TIndi() = default;


        TIndi& operator=(const TIndi& src);      /* Copy */
        bool operator==(const TIndi& indi2);     /* Return true if two tours are the same, false otherwise */


        void copy(const TIndi& other) {
            Solution::copy(other);
            fN = other.fN;
            fLink = other.fLink;
            fEvaluationValue = other.fEvaluationValue;
        }

        void copy(const SolutionBase& sol) {
            Solution::copy(sol);
            updateLinkInfo();
        }

        void define(int N);


        void transferSol(std::vector<int>& sol)const;
        void toCurSol(const std::vector<int>& sol);

        double distanceTo(const TIndi& indi2)const;


        void updateSolBase() {
            transferSol(variable().vect());
            resizeObjective(1);
            objective(0) = fEvaluationValue;
        }

        void updateLinkInfo()
        {
            toCurSol(variable().vect());
            fEvaluationValue = objective(0);
        }
    };

}

#endif
