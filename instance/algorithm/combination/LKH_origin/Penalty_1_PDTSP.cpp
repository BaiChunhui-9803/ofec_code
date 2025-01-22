#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
namespace LKH {
	GainType LKHAlg::Penalty_1_PDTSP()
	{
		Node *N;
		GainType P = 0;
		int Forward = SUCC(Depot)->Id != Depot->Id + DimensionSaved;
		int Load = 0, MaxLoad = std::numeric_limits<int>::min(), MinLoad = std::numeric_limits<int>::max();

		N = Depot;
		do {
			if (N->Id <= Dim) {
				Load += N->Demand;
				if (Load > MaxLoad)
					MaxLoad = Load;
				if (Load < MinLoad)
					MinLoad = Load;
				if (MaxLoad - MinLoad > Capacity) {
					P += MaxLoad - MinLoad - Capacity;
					if (P > CurrentPenalty ||
						(P == CurrentPenalty && CurrentGain <= 0))
						return CurrentPenalty + (CurrentGain > 0);
				}
			}
			N = Forward ? SUCC(N) : PREDD(N);
		} while (N != Depot);
		return P;
	}
}