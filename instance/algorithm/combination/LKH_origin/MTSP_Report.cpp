#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"

void LKH::LKHAlg::MTSP_Report(GainType Penalty, GainType Cost)
{
    GainType MinCost = std::numeric_limits<GainType>::max(), MaxCost = std::numeric_limits<GainType>::min();
    int MinSize = std::numeric_limits<int>::max(), MaxSize = 0;
    Node *N = Depot, *NextN;
    int Forward = SUCC(Depot)->Id != Depot->Id + DimensionSaved;
    int SalesmenUsed = 0;

    do {
        GainType Cost = 0;
        int Size = -1;
        do {
            NextN = Forward ? SUCC(N) : PREDD(N);
            Cost += (this->*C)(N, NextN) - N->Pi - NextN->Pi;
            if (NextN->Id > DimensionSaved)
                NextN = Forward ? SUCC(NextN) : PREDD(NextN);
            Size++;
        } while ((N = NextN)->DepotId == 0);
        Cost /= Precision;
        if (Cost < MinCost)
            MinCost = Cost;
        if (Cost > MaxCost)
            MaxCost = Cost;
        if (Size < MinSize)
            MinSize = Size;
        if (Size > MaxSize)
            MaxSize = Size;
        if (Size > 0)
            SalesmenUsed++;
    } while (N != Depot);
    if (MTSPMinSize == 0)
        printff("  Salesmen/vehicles used = %d\n", SalesmenUsed);
    printff("  Size.min = %d, Size.max = %d, Penalty = " GainFormat "\n",
            MinSize, MaxSize, Penalty);
    printff("  Cost.min = " GainFormat ", Cost.max = " GainFormat,
            MinCost, MaxCost);
    printff(", Cost.sum = " GainFormat "\n", Cost);
}
