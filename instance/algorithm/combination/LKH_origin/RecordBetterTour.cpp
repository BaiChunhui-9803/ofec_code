#include "./INCLUDE/LKH.h"

/*
 * The RecordBetterTour function is called by FindTour each time
 * the LinKernighan function has returned a better tour.
 *
 * The function records the tour in the BetterTour array and in the
 * BestSuc field of each node. Furthermore, for each node the previous
 * value of BestSuc is saved in the NextBestSuc field.
 *
 * Recording a better tour in the BetterTour array when the problem is
 * asymmetric requires special treatment since the number of nodes has
 * been doubled.
 */

void LKH::LKHAlg::RecordBetterTour()
{
    Node *N = Asymmetric ? Depot : FirstNode;
    Node *Stop = N;

    if (!Asymmetric) {
        int i = 1;
        do
            BetterTour[i++] = N->Id;
        while ((N = N->Suc) != Stop);
    } else if (N->Suc->Id != DimensionSaved + N->Id) {
        int i = 1;
        do
            if (N->Id <= DimensionSaved)
                BetterTour[i++] = N->Id;
        while ((N = N->Suc) != Stop);
    } else {
        int i = DimensionSaved;
        do
            if (N->Id <= DimensionSaved)
                BetterTour[i--] = N->Id;
        while ((N = N->Suc) != Stop);
    }
    BetterTour[0] = BetterTour[DimensionSaved];
    do {
        N->NextBestSuc = N->BestSuc;
        N->BestSuc = N->Suc;
    }
    while ((N = N->Suc) != Stop);
}


void  LKH::LKHAlg::RecordCurTour(std::vector<int>& sol)const {
    sol.resize(DimensionSaved + 1);
    Node* N = Asymmetric ? Depot : FirstNode;
    Node* Stop = N;
    if (!Asymmetric) {
        int i = 1;
        do
            sol[i++] = N->Id;
        while ((N = N->Suc) != Stop);
    }
    else if (N->Suc->Id != DimensionSaved + N->Id) {
        int i = 1;
        do
            if (N->Id <= DimensionSaved)
                sol[i++] = N->Id;
        while ((N = N->Suc) != Stop);
    }
    else {
        int i = DimensionSaved;
        do
            if (N->Id <= DimensionSaved)
                sol[i--] = N->Id;
        while ((N = N->Suc) != Stop);
    }
    sol[0] = sol[DimensionSaved];
    sol.pop_back();
 //   for (auto& it : sol) --it;
}
