#include "./INCLUDE/LKH.h"

void LKH::LKHAlg::SOP_Report(GainType Cost)
{
    std::cout << " Cost = " << Cost << "_" << CurrentPenalty << std::endl;
    //printff("  Cost = " GainFormat "_" GainFormat "\n",
    //        CurrentPenalty, Cost);
}
