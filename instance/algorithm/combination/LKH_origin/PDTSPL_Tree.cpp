#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
namespace LKH {
#define Dad Nearest
#define FirstSon OldPred
#define LastSon OldSuc
#define NextTourNode Mark

	void Perturb(LKHAlg  *Alg);
	void InsertSubTree(LKHAlg::Node *N, LKHAlg::Node *Prev, LKHAlg::Node *Dad);
	void RemoveSubTree(LKHAlg::Node *N, LKHAlg  *Alg);
	GainType Tree2Tour(LKHAlg *Alg);
	void Tour2Tree(LKHAlg *Alg);

	GainType TreeCost(LKHAlg *Alg)
	{
		GainType Cost = 0;
		LKHAlg::Node *N = Alg->Depot;
		do
			Cost += (Alg->*(Alg->C))(N, N->Suc) - N->Pi - N->Suc->Pi;
		while ((N = N->Suc) != Alg->Depot);
		return Cost / Alg->Precision;
	}

	void Tour2Tree(LKHAlg *Alg)
	{
		LKHAlg::Node *Current = Alg->Depot, *N;
		int *Stack, StackTop = -1;
		int Forward = ((Alg->Depot)->Parent ? (Alg->Reversed == (Alg->Depot)->Parent->Reversed ? (Alg->Depot)->Suc : (Alg->Depot)->Pred) : \
			Alg->Reversed ? (Alg->Depot)->Pred : (Alg->Depot)->Suc)->Id != Alg->Depot->Id + Alg->DimensionSaved;
		assert(Stack = (int *)malloc(Alg->Dim * sizeof(int)));
		N = Alg->FirstNode;
		do
			N->Dad = N->FirstSon = N->LastSon = N->Next = 0;
		while ((N = N->Suc) != Alg->FirstNode);
		Stack[++StackTop] = Alg->Depot->Id;
		do {
			if (Current->Delivery) {
				N = &Alg->NodeSet[Stack[StackTop]];
				Current->Dad = N;
				Current->Prev = N->LastSon;
				if (!N->FirstSon)
					N->FirstSon = N->LastSon = Current;
				else {
					N->LastSon->Next = Current;
					N->LastSon = Current;
				}
				Stack[++StackTop] = Current->Id;
			}
			else if (Current->Pickup)
				StackTop--;
			Current = Forward ? ((Current)->Parent ? (Alg->Reversed == (Current)->Parent->Reversed ? (Current)->Suc : (Current)->Pred) : \
				Alg->Reversed ? (Current)->Pred : (Current)->Suc) : ((N)->Parent ? (Alg->Reversed == (Current)->Parent->Reversed ? (Current)->Pred : (Current)->Suc) : Alg->Reversed ? (Current)->Suc : (Current)->Pred);
		} while (Current != Alg->Depot);
	}

	void RemoveSubTree(LKHAlg::Node *N, LKHAlg *Alg)
	{
		LKHAlg::Node *Prev = N->Prev;
		if (Prev)
			Prev->Next = N->Next;
		if (N->Next)
			N->Next->Prev = Prev;
		if (N->Dad->FirstSon == N)
			N->Dad->FirstSon = N->Next;
		if (N->Dad->LastSon == N)
			N->Dad->LastSon = Prev;
		N->Dad = N->Prev = N->Next = 0;
	}

	thread_local LKHAlg::Node *BestPrev, *BestDad;
	thread_local GainType BestMin;

	void FindBestPrevDad(LKHAlg::Node *N, LKHAlg::Node *Root, LKHAlg *Alg) {
		LKHAlg::Node *Son, *Prev;
		GainType d;
		if (Root == Alg->Depot) {
			if (!Alg->Depot->FirstSon) {
				BestPrev = 0;
				BestDad = Alg->Depot;
				return;
			}
			BestMin = std::numeric_limits<int>::max();
		}
		Prev = Root->LastSon;
		if (Prev) {
			InsertSubTree(N, Prev, Root);
			if ((d = Alg->Random() /* Tree2Tour() */) < BestMin) {
				BestMin = d;
				BestPrev = Prev;
				BestDad = Root;
			}
			RemoveSubTree(N, Alg);
		}
		Son = Root->FirstSon;
		while (1) {
			Prev = Son ? Son->Prev : 0;
			InsertSubTree(N, Prev, Root);
			if ((d = Alg->Random()/* Tree2Tour() */) < BestMin) {
				BestMin = d;
				BestPrev = Prev;
				BestDad = Root;
			}
			RemoveSubTree(N, Alg);
			if (!Son)
				break;
			FindBestPrevDad(N, Son, Alg);
			Son = Son->Next;
		}
	}

	void InsertSubTree(LKHAlg::Node *N, LKHAlg::Node *Prev, LKHAlg::Node *Dad) {
		N->Prev = Prev;
		if (Prev) {
			N->Next = Prev->Next;
			if (Prev->Next)
				Prev->Next->Prev = N;
			Prev->Next = N;
		}
		else
			N->Next = Dad->FirstSon;
		if (N->Next)
			N->Next->Prev = N;
		if (!Dad->FirstSon)
			Dad->FirstSon = Dad->LastSon = N;
		else if (!N->Prev)
			Dad->FirstSon = N;
		else if (!N->Next)
			Dad->LastSon = N;
		N->Dad = Dad;
	}

	void Perturb(LKHAlg  *Alg)
	{
		int Count = 0, i, j;
		LKHAlg::Node **TreeSet, *N;
		assert(TreeSet = (LKHAlg::Node **) malloc(Alg->Dim * sizeof(LKHAlg::Node *)));
		N = Alg->Depot;
		do {
			if (N->Delivery)
				TreeSet[Count++] = N;
		} while ((N = N->Suc) != Alg->Depot);
		for (i = 1; i < Count; i++) {
			N = TreeSet[i];
			TreeSet[i] = TreeSet[j = Alg->Random() % (i + 1)];
			TreeSet[j] = N;
		}
		Tour2Tree(Alg);
		for (i = 0; i < 40; i++) {
			N = TreeSet[i];
			RemoveSubTree(N, Alg);
			FindBestPrevDad(N, Alg->Depot, Alg);
			InsertSubTree(N, BestPrev, BestDad);
		}
		Tree2Tour(Alg);
	}

	static LKHAlg::Node *Last;

	void DFS(LKHAlg::Node *Current, LKHAlg *Alg)
	{
		LKHAlg::Node *N = Current == Alg->Depot || Current->Delivery ?
			Current : &Alg->NodeSet[Current->Pickup];
		if (Last)
			Last = Last->NextTourNode = N;
		else
			Last = N;
		N = Current->FirstSon;
		while (N) {
			DFS(N, Alg);
			N = N->Next;
		}
		if (Current != Alg->Depot) {
			N = Current->Delivery ? &Alg->NodeSet[Current->Delivery] : Current;
			Last = Last->NextTourNode = N;
		}
	}

	GainType Tree2Tour(LKHAlg *Alg)
	{
		LKHAlg::Node *N;
		GainType Cost;
		Last = 0;
		DFS(Alg->Depot, Alg);
		Last->NextTourNode = Alg->Depot;
		N = Alg->Depot;
		LKH::LKHAlg::Follow(N, N);
		do {
			LKH::LKHAlg::Follow(N->NextTourNode, N);
		} while ((N = N->NextTourNode) != Alg->Depot);
		do {
			LKH::LKHAlg::Precede(&Alg->NodeSet[N->Id + Alg->DimensionSaved], N);
		} while ((N = N->NextTourNode) != Alg->Depot);
		Cost = 0;
		N = Alg->Depot;
		do
			Cost += (Alg->*(Alg->C))(N, N->Suc) - N->Pi - N->Suc->Pi;
		while ((N = N->Suc) != Alg->Depot);
		//    printff("Cost = %lld\n", Cost / Precision);
		assert((Alg->*(Alg->Penalty))() == 0);
		return Cost / Alg->Precision;
	}
}