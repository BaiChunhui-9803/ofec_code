#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
namespace LKH {
	GainType LKHAlg::Penalty_CVRP()
	{
		static Node *StartRoute = 0;
		Node *N, *CurrentRoute;
		GainType DemandSum, DistanceSum, P = 0;

		if (!StartRoute)
			StartRoute = Depot;
		if (StartRoute->Id > DimensionSaved)
			StartRoute -= DimensionSaved;
		N = StartRoute;
		do {
			CurrentRoute = N;
			DemandSum = 0;
			do
				DemandSum += N->Demand;
			while ((N = SUCC(N))->DepotId == 0);
			if (DemandSum > Capacity &&
				((P += DemandSum - Capacity) > CurrentPenalty ||
				(P == CurrentPenalty && CurrentGain <= 0))) {
				StartRoute = CurrentRoute;
				return CurrentPenalty + (CurrentGain > 0);
			}
			if (DistanceLimit != std::numeric_limits<double>::max()) {
				DistanceSum = 0;
				N = CurrentRoute;
				do {
					DistanceSum += ((this->*C)(N, SUCC(N)) - N->Pi - SUCC(N)->Pi) /
						Precision;
					if (!N->DepotId)
						DistanceSum += N->ServiceTime;
				} while ((N = SUCC(N))->DepotId == 0);
				if (DistanceSum > DistanceLimit &&
					((P += DistanceSum - DistanceLimit) > CurrentPenalty ||
					(P == CurrentPenalty && CurrentGain <= 0))) {
					StartRoute = CurrentRoute;
					return CurrentPenalty + (CurrentGain > 0);
				}
			}
		} while (N != StartRoute);
		return P;
	}
}