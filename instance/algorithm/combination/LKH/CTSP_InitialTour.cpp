#include "./INCLUDE/LKH.h"
  
/* The CTSP_InitialTour function computes an initial tour for a
 * colored TSP.
 */
namespace LKH {
	GainType LKHAlg::CTSP_InitialTour()
	{
		Node *N = 0, *NextN = 0;
		GainType Cost;
		int Set;
		double EntryTime = GetTime();

		if (TraceLevel >= 1)
			printff("CTSP = ");
		assert(!Asymmetric);
		for (Set = 2; Set <= Salesmen; Set++)
			Follow(&NodeSet[Dim + Set - 1],
				Set == 2 ? Depot : &NodeSet[Dim + Set - 2]);
		N = Depot;
		do {
			NextN = N->Suc;
			if (N->DepotId == 0) {
				Set = N->Color != 0 ? N->Color : Random() % Salesmen + 1;
				Follow(N, Set == 1 ? Depot : &NodeSet[Dim + Set - 1]);
			}
		} while ((N = NextN) != Depot);
		Cost = 0;
		N = FirstNode;
		do
			Cost += (this->*C)(N, N->Suc) - N->Pi - N->Suc->Pi;
		while ((N = N->Suc) != FirstNode);
		Cost /= Precision;
		CurrentPenalty = std::numeric_limits<GainType>::max();
		CurrentPenalty = (this->*Penalty)();
		if (TraceLevel >= 1) {
			if (Salesmen > 1 || ProblemType == SOP)
				printff("%lld" "_" "%lld", CurrentPenalty, Cost);
			else
				printff("%lld", Cost);
			if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
				printff(", Gap = %0.2f%%", 100.0 * (Cost - Optimum) / Optimum);
			printff(", Time = %0.2f sec.\n", fabs(GetTime() - EntryTime));
		}
		return Cost;
	}
}