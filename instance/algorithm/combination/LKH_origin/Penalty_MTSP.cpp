#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
namespace LKH {
	GainType LKHAlg::Penalty_MTSP_MINSUM()
	{
		int Forward = SUCC(Depot)->Id != Depot->Id + DimensionSaved;
		Node *N = Depot, *NextN;
		GainType P = 0, DistanceSum;
		int MinSize = std::numeric_limits<int>::max(), MaxSize = 0;

		do {
			int Size = -1;
			do {
				Size++;
				NextN = Forward ? SUCC(N) : PREDD(N);
				if (NextN->Id > DimensionSaved)
					NextN = Forward ? SUCC(NextN) : PREDD(NextN);
			} while ((N = NextN)->DepotId == 0);
			if (Size < MinSize)
				MinSize = Size;
			if (Size > MaxSize)
				MaxSize = Size;
			if (Size > 0)
				Salesmen++;
		} while (N != Depot);
		if (MTSPMaxSize < Dimension - Salesmen && MaxSize > MTSPMaxSize)
			P += MaxSize - MTSPMaxSize;
		if (MTSPMinSize >= 1 && MinSize < MTSPMinSize)
			P += MTSPMinSize - MinSize;
		if (DistanceLimit != std::numeric_limits<double>::max()) {
			do {
				if (P > CurrentPenalty ||
					(P == CurrentPenalty && CurrentGain <= 0))
					return CurrentPenalty + (CurrentGain > 0);
				DistanceSum = 0;
				do {
					NextN = Forward ? SUCC(N) : PREDD(N);
					DistanceSum += ((this->*C)(N, NextN) - N->Pi - NextN->Pi) /
						Precision;
					if (NextN->Id > DimensionSaved)
						NextN = Forward ? SUCC(NextN) : PREDD(NextN);
				} while ((N = NextN)->DepotId == 0);
				if (DistanceSum > DistanceLimit)
					P += DistanceSum - DistanceLimit;
			} while (N != Depot);
		}
		return P;
	}

	GainType LKHAlg::Penalty_MTSP_MINMAX()
	{
		int Forward = SUCC(Depot)->Id != Depot->Id + DimensionSaved;
		static Node *StartRoute = 0;
		Node *N, *NextN, *CurrentRoute;
		GainType Cost, MaxCost = std::numeric_limits<GainType>::min();

		if (!StartRoute)
			StartRoute = Depot;
		if (StartRoute->Id > DimensionSaved)
			StartRoute -= DimensionSaved;
		N = StartRoute;
		do {
			Cost = 0;
			CurrentRoute = N;
			do {
				NextN = Forward ? SUCC(N) : PREDD(N);
				Cost += (this->*C)(N, NextN) - N->Pi - NextN->Pi;
				if (NextN->Id > DimensionSaved)
					NextN = Forward ? SUCC(NextN) : PREDD(NextN);
			} while ((N = NextN)->DepotId == 0);
			Cost /= Precision;
			if (Cost > MaxCost) {
				if (Cost > CurrentPenalty ||
					(Cost == CurrentPenalty && CurrentGain <= 0)) {
					StartRoute = CurrentRoute;
					return CurrentPenalty + (CurrentGain > 0);
				}
				MaxCost = Cost;
			}
		} while (N != StartRoute);
		return MaxCost;
	}

	GainType LKHAlg::Penalty_MTSP_MINMAX_SIZE()
	{
		int Forward = SUCC(Depot)->Id != Depot->Id + DimensionSaved;
		static Node *StartRoute = 0;
		Node *N, *NextN, *CurrentRoute;
		int Size, MaxSize = std::numeric_limits<int>::min();

		if (!StartRoute)
			StartRoute = Depot;
		if (StartRoute->Id > DimensionSaved)
			StartRoute -= DimensionSaved;
		N = StartRoute;
		do {
			Size = 0;
			CurrentRoute = N;
			do {
				NextN = Forward ? SUCC(N) : PREDD(N);
				Size++;
				if (NextN->Id > DimensionSaved)
					NextN = Forward ? SUCC(NextN) : PREDD(NextN);
			} while ((N = NextN)->DepotId == 0);
			if (Size > MaxSize) {
				if (Size > CurrentPenalty ||
					(Size == CurrentPenalty && CurrentGain <= 0)) {
					StartRoute = CurrentRoute;
					return CurrentPenalty + (CurrentGain > 0);
				}
				MaxSize = Size;
			}
		} while (N != StartRoute);
		return MaxSize;
	}
}