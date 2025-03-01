#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
namespace LKH {
	GainType LKHAlg::Penalty_M_PDTSP()
	{
		Node *N;
		GainType P = 0;
		int Forward = SUCC(Depot)->Id != Depot->Id + DimensionSaved;
		int Load = 0, MaxLoad = std::numeric_limits<int>::min(), MinLoad = std::numeric_limits<int>::max(), k;
		static int *ComLoad = 0, *Min = 0;

		if (!ComLoad) {
			assert(ComLoad = (int *)malloc(DemandDimension * sizeof(int)));
			assert(Min = (int *)malloc(DemandDimension * sizeof(int)));
		}
		for (k = 0; k < DemandDimension; k++) {
			ComLoad[k] = 0;
			Min[k] = std::numeric_limits<int>::max();
		}
		N = Depot;
		do {
			if (N->Id <= Dim) {
				Load += N->Demand;
				MinLoad = 0;
				for (k = 0; k < DemandDimension; k++) {
					ComLoad[k] += N->M_Demand[k];
					Load += N->M_Demand[k];
					if (ComLoad[k] < Min[k])
						Min[k] = ComLoad[k];
					MinLoad += Min[k];
					if (ProblemType == M1_PDTSP && ComLoad[k] < 0)
						P -= ComLoad[k];
				}
				if (Load > MaxLoad)
					MaxLoad = Load;
				if (MaxLoad - MinLoad > Capacity)
					P += MaxLoad - MinLoad - Capacity;
				if (P > CurrentPenalty)
					return CurrentPenalty + 1;
			}
			N = Forward ? SUCC(N) : PREDD(N);
		} while (N != Depot);
		return P;
	}
}