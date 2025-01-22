#include "./INCLUDE/LKH.h"

void LKH::LKHAlg::StatusReport(GainType Cost, double EntryTime, char *Suffix)
{
    if (Penalty) {
        printff("Cost = " GainFormat "_" GainFormat, CurrentPenalty, Cost);
        if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0) {
            if (ProblemType != CCVRP && ProblemType != TRP &&
                ProblemType != MLP &&
                MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE)
                printff(", Gap = %0.4f%%",
                        100.0 * (Cost - Optimum) / Optimum);
            else
                printff(", Gap = %0.4f%%",
                        100.0 * (CurrentPenalty - Optimum) / Optimum);
        }
        printff(", Time = %0.2f sec. %s",
                fabs(GetTime() - EntryTime), Suffix);
    } else {
        printff("Cost = " GainFormat, Cost);
        if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
            printff(", Gap = %0.4f%%", 100.0 * (Cost - Optimum) / Optimum);
        printff(", Time = %0.2f sec. %s%s",
                fabs(GetTime() - EntryTime), Suffix,
                Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
    }
    printff("\n");
}
