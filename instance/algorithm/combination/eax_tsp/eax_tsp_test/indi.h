#ifndef EAX_TSP2_INDI_H
#define EAX_TSP2_INDI_H

#include <vector>
#include "../../../../../core/problem/solution.h"

namespace ofec::eax_tsp2 {

    class TIndi : public Solution<VariableVector<int>>
    {
    public:


        TIndi();
        ~TIndi() = default;
        void define(int N);
        void transferSol(std::vector<int>& sol)const;
        void toCurSol(const std::vector<int>& sol);

        double distanceTo(const TIndi& indi2)const;

        TIndi& operator=(const TIndi& src);      /* Copy */
        bool operator==(const TIndi& indi2);     /* Return true if two tours are the same, false otherwise */
        int fN = 0;     /* Number of cities */
        std::vector<std::vector<int>> fLink;
        // int **fLink;         /* fLink[i][]: two vertices adjacent to i */
        double fEvaluationValue; /* Tour length of */


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
