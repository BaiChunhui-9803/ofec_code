#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
namespace LKH {
	GainType Penalty_MVRPB(LKHAlg *Alg)
	{
		static LKHAlg::Node *StartRoute = 0;
		LKHAlg::Node *N, *NextN, *CurrentRoute;
		GainType Duration, P = 0;
		int Load, DeliverySum, PickupSum;
		int Forward = ((Alg->Depot)->Parent ? (Alg->Reversed == (Alg->Depot)->Parent->Reversed ? (Alg->Depot)->Suc : (Alg->Depot)->Pred) : \
			Alg->Reversed ? (Alg->Depot)->Pred : (Alg->Depot)->Suc)->Id != Alg->Depot->Id + Alg->DimensionSaved;

		if (!StartRoute)
			StartRoute = Alg->Depot;
		if (StartRoute->Id > Alg->DimensionSaved)
			StartRoute -= Alg->DimensionSaved;
		N = StartRoute;
		do {
			CurrentRoute = N;
			DeliverySum = PickupSum = 0;
			do {
				if (N->Id <= Alg->Dim && N != Alg->Depot) {
					DeliverySum += N->Delivery;
					PickupSum += N->Pickup;
				}
				NextN = Forward ? ((N)->Parent ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Suc : (N)->Pred) : \
					Alg->Reversed ? (N)->Pred : (N)->Suc) : ((N)->Parent ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Pred : (N)->Suc) : Alg->Reversed ? (N)->Suc : (N)->Pred);
				N = Forward ? ((NextN)->Parent ? (Alg->Reversed == (NextN)->Parent->Reversed ? (NextN)->Suc : (NextN)->Pred) : \
					Alg->Reversed ? (NextN)->Pred : (NextN)->Suc) : ((NextN)->Parent ? (Alg->Reversed == (NextN)->Parent->Reversed ? (NextN)->Pred : (NextN)->Suc) : Alg->Reversed ? (NextN)->Suc : (NextN)->Pred);
			} while (N->DepotId == 0);
			if (DeliverySum > Alg->Capacity)
				P += DeliverySum - Alg->Capacity;
			if (PickupSum > Alg->Capacity)
				P += PickupSum - Alg->Capacity;
			if (P > Alg->CurrentPenalty ||
				(P == Alg->CurrentPenalty && Alg->CurrentGain <= 0)) {
				StartRoute = CurrentRoute;
				return Alg->CurrentPenalty + (Alg->CurrentGain > 0);
			}
			Duration = Alg->Depot->Earliest;
			Load = DeliverySum;
			N = CurrentRoute;
			do {
				if (N->Id <= Alg->Dim && N != Alg->Depot) {
					Load -= N->Delivery;
					Load += N->Pickup;
					if (Load > Alg->Capacity)
						P += Load - Alg->Capacity;
					if (Duration < N->Earliest)
						Duration = N->Earliest;
					if (Duration > N->Latest)
						P += Duration - N->Latest;
					if (P > Alg->CurrentPenalty) {
						StartRoute = CurrentRoute;
						return Alg->CurrentPenalty + 1;
					}
					Duration += N->ServiceTime;
				}
				NextN = Forward ? ((N)->Parent ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Suc : (N)->Pred) : \
					Alg->Reversed ? (N)->Pred : (N)->Suc) : ((N)->Parent ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Pred : (N)->Suc) : Alg->Reversed ? (N)->Suc : (N)->Pred);
				Duration += ((Alg->*(Alg->C))(N, NextN) - N->Pi - NextN->Pi) / Alg->Precision;
				if (Alg->DistanceLimit != 0 && Duration > Alg->DistanceLimit)
					P += Duration - Alg->DistanceLimit;
				N = Forward ? ((NextN)->Parent ? (Alg->Reversed == (NextN)->Parent->Reversed ? (NextN)->Suc : (NextN)->Pred) : \
					Alg->Reversed ? (NextN)->Pred : (NextN)->Suc) : ((NextN)->Parent ? (Alg->Reversed == (NextN)->Parent->Reversed ? (NextN)->Pred : (NextN)->Suc) : Alg->Reversed ? (NextN)->Suc : (NextN)->Pred);
			} while (N->DepotId == 0);
			if (Duration > Alg->Depot->Latest &&
				(P += Duration - Alg->Depot->Latest) > Alg->CurrentPenalty) {
				StartRoute = CurrentRoute;
				return Alg->CurrentPenalty + 1;
			}
		} while (N != StartRoute);
		return P;
	}
}