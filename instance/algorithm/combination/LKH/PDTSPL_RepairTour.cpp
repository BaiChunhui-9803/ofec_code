#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"

/* The PDTSPL_RepairTour function creates a feasible tour from a given tour.
 */
namespace LKH {
	GainType PDTSPL_RepairTour(LKHAlg *Alg)
	{
		LKHAlg::Node **Fringe, *First, *Last, *N;
		int FringeNodes = 0, Min, i, j;
		GainType Cost;
		int *Stack, StackTop = -1;
		int Forward = ((Alg->Depot)->Parent ? (Alg->Reversed == (Alg->Depot)->Parent->Reversed ? (Alg->Depot)->Suc : (Alg->Depot)->Pred) : \
			Alg->Reversed ? (Alg->Depot)->Pred : (Alg->Depot)->Suc)->Id != Alg->Depot->Id + Alg->DimensionSaved;//((Depot)->Parent ? SUC(Depot) :Reversed ? (Depot)->Pred : (Depot)->Suc)
		//SUC(a) (Reversed == (Depot)->Parent->Reversed ? (Depot)->Suc : (Depot)->Pred)
		assert(Stack = (int *)malloc(Alg->Dim * sizeof(int)));
		assert(Fringe = (LKHAlg::Node **) malloc(Alg->Dim * sizeof(LKHAlg::Node *)));
		First = Last = Alg->Depot;
		First->Prev = First->Next = First;
		FringeNodes = 0;
		N = Alg->Depot;
		i = 0;
		do {
			if (N->Delivery)
				Fringe[FringeNodes++] = N;
			if (N->Delivery || N->Pickup)
				N->Rank = ++i;
			N = Forward ? ((N)->Parent ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Suc : (N)->Pred) : \
				Alg->Reversed ? (N)->Pred : (N)->Suc) : ((N)->Parent ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Pred : (N)->Suc) : Alg->Reversed ? (N)->Suc : (N)->Pred);//((N)->Parent ? PRED(N) :Reversed ? (N)->Suc : (N)->Pred)
			//PRED(N) (Reversed == (N)->Parent->Reversed ? (N)->Pred : (N)->Suc)
		} while (N != Alg->Depot);
		while (FringeNodes > 0) {
			Min = std::numeric_limits<int>::max();
			for (j = 0; j < FringeNodes; j++) {
				N = Fringe[j];
				if ((N->Delivery ||
					(StackTop >= 0 && N->Pickup == Stack[StackTop])) &&
					N->Rank < Min) {
					Min = N->Rank;
					i = j;
				}
			}
			N = Fringe[i];
			Fringe[i] = Fringe[--FringeNodes];
			if (N->Delivery) {
				Stack[++StackTop] = N->Id;
				Fringe[FringeNodes++] = &Alg->NodeSet[N->Delivery];
			}
			else
				StackTop--;
			N->Prev = Last;
			N->Next = First;
			First->Prev = Last->Next = N;
			Last = N;
		}
		free(Fringe);
		free(Stack);
		N = First;
		Alg->Follow(N, N);
		do {
			Alg->Follow(N->Next, N);
		} while ((N = N->Next) != First);
		do {
			Alg->Precede(&Alg->NodeSet[N->Id + Alg->DimensionSaved], N);
		} while ((N = N->Next) != First);
		Cost = 0;
		N = First;
		do
			Cost += (Alg->*(Alg->C))(N, N->Suc) - N->Pi - N->Suc->Pi;
		while ((N = N->Suc) != First);
		Alg->CurrentPenalty = 0;
		return Cost / Alg->Precision;
	}
}