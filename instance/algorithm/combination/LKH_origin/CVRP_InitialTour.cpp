#include "./INCLUDE/LKH.h"

/* The CVRP_InitialTour function computes an initial tour using the 
 * Clarke and Wright savings algorithm.
 */
namespace LKH {
#define ITERATIONS 10

#define Degree V
#define Size LastV
#define Load Loc

	static LKHAlg::Node *Find(LKHAlg::Node * v);
	static void Union(LKHAlg::Node * x, LKHAlg::Node * y, LKHAlg *Alg);
	static void MakeSets(LKHAlg *Alg);
	static void Distribute(int Constrained, double R, LKHAlg *Alg);
	static int compareValue(const void *s1, const void *s2);
	static void CreateS(LKHAlg *Alg);

	typedef struct Saving {
		int Value, i, j;
	} Saving;

	static thread_local Saving *S = 0;
	static thread_local int Sets, SSize;

	GainType LKHAlg::CVRP_InitialTour()
	{
		Node *N, *Last, *Next, *Tour;
		int s, Dim = Dimension - Salesmen + 1, it;
		GainType Cost, BestCost = std::numeric_limits<GainType>::max(), BestPenalty = std::numeric_limits<GainType>::max();
		double EntryTime = GetTime();

		if (TraceLevel >= 1)
			printff("CVRP = ");
		assert(!Asymmetric);
		if (!S)
			CreateS(this);
		for (it = 0; it < ITERATIONS; it++) {
			for (s = 0; s < Salesmen; s++) {
				Tour = s == 0 ? Depot : &NodeSet[Dim + s];
				if (Tour == FirstNode)
					FirstNode = Tour->Suc;
				Follow(Tour, Depot);
				Tour->Dad = Tour;
				Tour->Prev = Tour->Next = 0;
				Tour->Degree = 0;
				Tour->Size = 1;
			}
			MakeSets(this);
			if (it > 0)
				Distribute(1, 0.01, this);
			if (Sets > Salesmen)
				Distribute(1, 0, this);
			if (Sets > Salesmen) {
				if (BestPenalty == 0)
					continue;
				Distribute(0, 0, this);
			}
			Sets += Salesmen;
			N = FirstNode;
			do {
				if (N->Degree <= 1 && (s = N->Special) > 0) {
					Tour = s == 1 ? Depot : &NodeSet[Dim + s - 1];
					Union(N, Tour, this);
					s = s == Salesmen ? 1 : s + 1;
					if (N->Degree <= 1) {
						Tour = s == 1 ? Depot : &NodeSet[Dim + s - 1];
						Union(N, Tour, this);
					}
				}
			} while ((N = N->Suc) != FirstNode);
			while (Sets > 0) {
				if (N->Degree <= 1) {
					for (s = 1; s <= Salesmen; s++) {
						Tour = s == 1 ? Depot : &NodeSet[Dim + s - 1];
						if (Tour->Degree <= 1 &&
							(Sets == 1 || Find(N) != Find(Tour))) {
							Union(N, Tour, this);
							if (N->Degree <= 1)
								N = N->Pred;
							break;
						}
					}
				}
				N = N->Suc;
			}
			Last = N = FirstNode = Depot;
			do {
				Next = N->Next != Last ? N->Next : N->Prev;
				Follow(N, Last);
				Last = N;
			} while ((N = Next) != FirstNode);
			N = FirstNode = Depot;
			Cost = 0;
			do
				Cost += (this->*C)(N, N->Suc) - N->Pi - N->Suc->Pi;
			while ((N = N->Suc) != FirstNode);
			Cost /= Precision;
			Cost += ServiceTime * (Dim - 1);
			CurrentPenalty = std::numeric_limits<GainType>::max();
			CurrentPenalty = (this->*Penalty)();
			if (CurrentPenalty < BestPenalty ||
				(CurrentPenalty == BestPenalty && Cost < BestCost)) {
				N = FirstNode;
				while ((N = N->OldSuc = N->Suc) != FirstNode);
				BestCost = Cost;
				BestPenalty = CurrentPenalty;
			}
		}
		N = FirstNode;
		do
			(N->Suc = N->OldSuc)->Pred = N;
		while ((N = N->Suc) != FirstNode);
		Cost = BestCost;
		CurrentPenalty = BestPenalty;
		if (TraceLevel >= 1) {
			if (Salesmen > 1 || ProblemType == SOP)
				std::cout << CurrentPenalty << "\t" << Cost;
			//printff(GainFormat "_" GainFormat, CurrentPenalty, Cost);
			else
				std::cout << Cost << std::endl;
				//printff(GainFormat, Cost);
			if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
				std::cout << " Gap = " << 100.0 * (Cost - Optimum) / Optimum;
		//		printff(", Gap = %0.2f%%", 100.0 * (Cost - Optimum) / Optimum);

			std::cout << ", Time = " << fabs(GetTime() - EntryTime) << " sec.\n";
//			printff(", Time = %0.2f sec.\n", fabs(GetTime() - EntryTime));
		}
		if (Run == Runs) {
			free(S);
			S = 0;
		}
		return Cost;
	}

	void CreateS(LKHAlg *Alg)
	{
		int Dim = Alg->Dimension - Alg->Salesmen + 1, i, j;
		LKHAlg::Node *Ni, *Nj;
		SSize = 0;
		assert(S =
			(Saving *)malloc((Dim - 2) * (Dim - 1) / 2 * sizeof(Saving)));
		/* Compute savings */
		for (i = 1; i < Dim; i++) {
			Ni = &Alg->NodeSet[i];
			if (Ni == Alg->Depot)
				continue;
			for (j = i + 1; j <= Dim; j++) {
				Nj = &Alg->NodeSet[j];
				if (Nj == Alg->Depot || Alg->Forbidden(Ni, Nj))
					continue;
				S[SSize].Value = (LKH::LKHAlg::Fixed(Ni, Nj) || Alg->IsCommonEdge(Ni, Nj)) ? std::numeric_limits<int>::max() ://(Fixed(Ni, Nj) || IsCommonEdge(Ni, Nj))
					(Alg->*(Alg->OldDistance))(Ni, Alg->Depot) +
					(Alg->*(Alg->OldDistance))(Alg->Depot, Nj) - (Alg->*(Alg->OldDistance))(Ni, Nj);
				S[SSize].i = i;
				S[SSize].j = j;
				SSize++;
			}
		}
		/* Rank the savings in descending order */
		qsort(S, SSize, sizeof(Saving), compareValue);
	}

	static LKHAlg::Node *Find(LKHAlg::Node * v)
	{
		if (v != v->Dad)
			v->Dad = Find(v->Dad);
		return v->Dad;
	}

	static void Union(LKHAlg::Node * x, LKHAlg::Node * y, LKHAlg *Alg)
	{
		LKHAlg::Node *u = Find(x), *v = Find(y);
		if (u->Size < v->Size) {
			u->Dad = v;
			v->Size += u->Size;
			v->Load += u->Load;
			v->Cost += u->Cost + (Alg->*(Alg->OldDistance))(x, y) -
				(Alg->*(Alg->OldDistance))(x, Alg->Depot) - (Alg->*(Alg->OldDistance))(y, Alg->Depot);
		}
		else {
			v->Dad = u;
			u->Size += v->Size;
			u->Load += v->Load;
			u->Cost += v->Cost + (Alg->*(Alg->OldDistance))(x, y) -
				(Alg->*(Alg->OldDistance))(x, Alg->Depot) - (Alg->*(Alg->OldDistance))(y, Alg->Depot);
		}
		if (x->Degree++ == 0)
			x->Prev = y;
		else
			x->Next = y;
		if (y->Degree++ == 0)
			y->Prev = x;
		else
			y->Next = x;
		Sets--;
	}

	static void MakeSets(LKHAlg *Alg)
	{
		LKHAlg::Node *N = Alg->FirstNode;
		Sets = 0;
		do {
			N->Dad = N;
			N->Prev = N->Next = 0;
			N->Degree = 0;
			N->Size = 1;
			N->Load = N->Demand;
			N->Cost = 2 * (Alg->*(Alg->OldDistance))(N, Alg->Depot);
			Sets++;
		} while ((N = N->Suc) != Alg->FirstNode);
	}

	static void Distribute(int Constrained, double R, LKHAlg *Alg)
	{
		LKHAlg::Node *Ni, *Nj, *u, *v;
		int i;

		for (i = 0; i < SSize && Sets > Alg->Salesmen; i++) {
			if (R > 0 && Alg->Random() % 1000 <= 1000 * R)
				continue;
			Ni = &Alg->NodeSet[S[i].i];
			Nj = &Alg->NodeSet[S[i].j];
			if (Ni->Degree < 2 && Nj->Degree < 2) {
				u = Find(Ni);
				v = Find(Nj);
				if (u == v)
					continue;
				if (!Constrained ||
					(u->Load + v->Load <= Alg->Capacity &&
					(Alg->DistanceLimit == std::numeric_limits<double>::max() ||
						u->Cost + v->Cost + (Alg->*(Alg->OldDistance))(Ni, Nj) -
						(Alg->*(Alg->OldDistance))(Ni, Alg->Depot) - (Alg->*(Alg->OldDistance))(Nj, Alg->Depot) +
						(u->Size + v->Size) * Alg->ServiceTime <= Alg->DistanceLimit)))
					Union(Ni, Nj, Alg);
			}
		}
	}

	static int compareValue(const void *s1, const void *s2)
	{
		int v1 = ((Saving *)s1)->Value;
		int v2 = ((Saving *)s2)->Value;
		return v1 > v2 ? -1 : v1 == v2 ? 0 : 1;
	}
}