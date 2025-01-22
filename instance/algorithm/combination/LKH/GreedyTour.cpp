#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Heap.h"
namespace LKH {
#define Tail SucSaved

	/*
	 * The GreedyTour function computes a tour using either
	 *
	 *      (1) Nearest Neighbor (NEAREST-NEIGHBOR),
	 *      (2) Bentley's multiple fragment heuristic (GREEDY),
	 *      (3) Boruvka (BORUVKA), or
	 *      (4) Applegate, Cook and Rohe's Quick-Boruvka heuristic
	 *          (QUICK-BORUVKA).
	 *
	 *    J. L. Bentley,
	 *    Fast Algorithms for Geometric Traveling Salesman Problems,
	 *    ORSA Journal on Computing, Vol. 4, 4, pp. 387-411 (1992).
	 *
	 *    D. Applegate, W. Cook, and A. Rohe,
	 *    Chained Lin-Kernighan for large traveling salesman problems.
	 *    Technical Report No. 99887, Forschungsinstitut
	 *    fuer Diskrete Mathematik, University of Bonn, Germany (1999).
	 *
	 * The function returns the cost of the resulting tour.
	 *
	 * Nearest neighbor searching is performed by using the candidate graph.
	 * To provide different tours a randomized nearest neighbor search is used.
	 * The first found feasible neighbor is returned with probability 2/3.
	 *
	 * A heap contains the nearest neighbor edges when the multiple fragment
	 * heuristic is used.
	 */

#define Degree V        /* Number of edges currently adjacent to a node */
#define Mark LastV      /* Mark of a node during breadth-first search   */
#define Level BestPi    /* Search level */

	struct GreedyTourInfo {
		int EdgesInFragments;
		GainType Cost;
	};

	static LKHAlg::Node* NearestNeighbor(LKHAlg::Node* From, LKHAlg* Alg, GreedyTourInfo& info);
	static LKHAlg::Node* NearestInList(LKHAlg::Node* From, LKHAlg::Node* First, LKHAlg* Alg, GreedyTourInfo& info);
	static int MayBeAddedToFragments(LKHAlg::Node* From, LKHAlg::Node* To, LKHAlg* Alg, GreedyTourInfo& info);
	static void AddEdgeToFragments(LKHAlg::Node* From, LKHAlg::Node* To, LKHAlg* Alg, GreedyTourInfo& info);
	static void RemoveFromList(LKHAlg::Node* N, LKHAlg::Node** First, LKHAlg* Alg);
	static int compareX(const void* Na, const void* Nb);
	static int compareCost(const void* Na, const void* Nb);




	GainType LKHAlg::GreedyTour()
	{
		GreedyTourInfo info;


		Node *From, *To, *First, *Last = 0, **Perm;
		int Count, i;
		double EntryTime = GetTime();
		//if (TraceLevel >= 1) {
		//	if (InitialTourAlgorithm == BORUVKA)
		//		printff("Boruvka = ");
		//	else if (InitialTourAlgorithm == GREEDY)
		//		printff("Greedy = ");
		//	else if (InitialTourAlgorithm == NEAREST_NEIGHBOR)
		//		printff("Nearest-Neighbor = ");
		//	else if (InitialTourAlgorithm == QUICK_BORUVKA)
		//		printff("Quick-Boruvka = ");
		//}
		info.Cost = 0;
		info.EdgesInFragments = 0;
		From = FirstNode;
		do {
			From->Degree = 0;
			From->Tail = From;
			From->Mark = 0;
			From->Next = From->Suc;
			From->Pred = From->Suc = 0;
		} while ((From = From->Next) != FirstNode);
		Count = 0;
		for (;;) {
			if ((To = NearestNeighbor(From, this, info)) && FixedOrCommon(From, To))
				AddEdgeToFragments(From, To, this,info);
			else {
				if ((From->Nearest = To)) {
					if (InitialTourAlgorithm == GREEDY) {
						From->Rank = From->Cost;
						HeapLazyInsert(From, this);
					}
					else
						Count++;
				}
				if ((From = From->Next) == FirstNode)
					break;
			}
		}
		if (InitialTourAlgorithm == NEAREST_NEIGHBOR) {
			if (info.EdgesInFragments < Dimension) {
				while (From->Degree == 2)
					From = From->Tail;
				for (;;) {
					int Min = std::numeric_limits<int>::max(), d;
					Node *Nearest = 0;
					while ((To = NearestNeighbor(From, this, info))) {
						AddEdgeToFragments(From, To, this,info);
						if (info.EdgesInFragments == Dimension)
							break;
						From = To->Degree < 2 ? To : To->Tail;
						while (From->Degree == 2)
							From = From->Tail;
					}
					if (info.EdgesInFragments == Dimension)
						break;
					To = FirstNode;
					do {
						if (MayBeAddedToFragments(From, To, this,info) &&
							(!c || (this->*c)(From, To) < Min)
							&& (d = (this->*C)(From, To)) < Min) {
							Min = From->Cost = d;
							Nearest = To;
						}
					} while ((To = To->Next) != FirstNode);
					assert(Nearest);
					To = Nearest;
					AddEdgeToFragments(From, To, this,info);
					if (info.EdgesInFragments == Dimension)
						break;
					while (From->Degree == 2)
						From = From->Tail;
					From = To->Degree < 2 ? To : To->Tail;
				}
			}
		}
		else {
			if (InitialTourAlgorithm == GREEDY) {
				Heapify(this);
				while ((From = HeapDeleteMin(this))) {
					To = From->Nearest;
					if (MayBeAddedToFragments(From, To, this,info))
						AddEdgeToFragments(From, To, this,info);
					if ((From->Nearest = NearestNeighbor(From, this, info))) {
						From->Rank = From->Cost;
						HeapInsert(From, this);
					}
				}
			}
			else {
				assert(Perm = (Node **)malloc(Count * sizeof(Node *)));
				for (From = FirstNode, i = 0; i < Count; From = From->Next)
					if (From->Nearest)
						Perm[i++] = From;
				if (InitialTourAlgorithm == QUICK_BORUVKA) {
					qsort(Perm, Count, sizeof(Node *), compareX);
					for (i = 0; i < Count; i++) {
						From = Perm[i];
						if ((To = NearestNeighbor(From, this, info))) {
							AddEdgeToFragments(From, To, this,info);
							i--;
						}
					}
				}
				else if (InitialTourAlgorithm == BORUVKA) {
					while (Count > 0) {
						qsort(Perm, Count, sizeof(Node *), compareCost);
						for (i = 0; i < Count; i++) {
							From = Perm[i];
							To = From->Nearest;
							if (MayBeAddedToFragments(From, To, this,info))
								AddEdgeToFragments(From, To, this, info);
						}
						for (i = 0; i < Count;) {
							From = Perm[i];
							if (!(To = NearestNeighbor(From, this, info)))
								Perm[i] = Perm[--Count];
							else {
								From->Nearest = To;
								i++;
							}
						}
					}
				}
				free(Perm);
			}
			if (info.EdgesInFragments < Dimension) {
				/* Create a list of the degree-0 and degree-1 nodes */
				First = 0;
				From = FirstNode;
				do {
					if (From->Degree != 2) {
						From->OldSuc = First;
						if (First)
							First->OldPred = From;
						else
							Last = From;
						First = From;
					}
				} while ((From = From->Next) != FirstNode);
				First->OldPred = Last;
				Last->OldSuc = First;
				From = First;
				do {
					if ((From->Nearest = NearestInList(From, First, this,info))) {
						From->Rank = From->Cost;
						HeapLazyInsert(From, this);
					}
				} while ((From = From->OldSuc) != First);
				Heapify(this);
				/* Find the remaining fragments */
				while ((From = HeapDeleteMin(this))) {
					To = From->Nearest;
					if (MayBeAddedToFragments(From, To, this, info)) {
						AddEdgeToFragments(From, To, this, info);
						if (From->Degree == 2)
							RemoveFromList(From, &First, this);
						if (To->Degree == 2)
							RemoveFromList(To, &First, this);
					}
					if (From->Degree != 2
						&& (From->Nearest = NearestInList(From, First, this, info))) {
						From->Rank = From->Cost;
						HeapInsert(From, this);
					}
				}
				assert(info.EdgesInFragments == Dimension);
			}
		}
		/* Orient Pred and Suc so that the list of nodes represents a tour */
		To = FirstNode;
		From = To->Pred;
		do {
			if (To->Suc == From) {
				To->Suc = To->Pred;
				To->Pred = From;
			}
			From = To;
		} while ((To = From->Suc) != FirstNode);
		To->Pred = From;
		info.Cost /= Precision;
		CurrentPenalty = std::numeric_limits<GainType>::max();
		CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
	/*	if (TraceLevel >= 1) {
			if (Salesmen > 1 || ProblemType == SOP)
				printff(GainFormat "_" GainFormat, CurrentPenalty, info.Cost);
			else
				printff(GainFormat, info.Cost);
			if (Optimum != MINUS_INFINITY && Optimum != 0) {
				if (MTSPObjective == MINMAX || MTSPObjective == MINMAX_SIZE)
					printff(", Gap = %0.1f%%",
						100.0 * (CurrentPenalty - Optimum) / Optimum);
				else
					printff(", Gap = %0.1f%%",
						100.0 * (info.Cost - Optimum) / Optimum);
			}
			printff(", Time = %0.2f sec.\n", fabs(GetTime() - EntryTime));
		}*/
		return info.Cost;
	}

	/*
	 * The NearestNeighbor function returns for a given node, From, the "nearest"
	 * neighbor, To, for which the addition of the edge (From, To) will not make
	 * it impossible to complete a tour.
	 *
	 * The neighbor is determined by searching the candidate graph breadth-first,
	 * starting at the node From. If possible, the nearest node on level 1 is
	 * chosen. Otherwise, on level 2. And so on. A queue is used for this search.
	 *
	 * Note that this algorithm finds an approximation to the nearest neighbor.
	 * However, by giving candidate edges precedence over non-candidate edges the
	 * algorithm will often produce better results than an algorithm that
	 * determines the true nearest neighbors (for example, by using a KD-tree).
	 */
	static thread_local int mark = 0;

	void LKHAlg::freeGreedyTour()
	{
		mark = 0;
	}
	static LKHAlg::Node *NearestNeighbor(LKHAlg::Node * From, LKHAlg *Alg, GreedyTourInfo& info)
	{

		LKHAlg::Candidate *NN;
		LKHAlg::Node *To, *N, *First = 0, *Last = 0, *Nearest = 0;
		int MaxLevel = Alg->Dimension, Min = std::numeric_limits<int>::max(), d;

		if (From->Degree == 2)
			return 0;
		for (NN = From->CandidateSet; (To = NN->To); NN++) {
			if ((Alg->Fixed(From, To) || Alg->IsCommonEdge(From, To)) && MayBeAddedToFragments(From, To, Alg, info)) {
				From->Cost = NN->Cost;
				return To;
			}
		}
		From->Level = 0;
		if (++mark == 0)
			mark = 1;
		From->Mark = mark;
		/* Insert From into an empty queue */
		First = Last = From;
		From->OldSuc = 0;

		while ((N = First) && N->Level < MaxLevel) {
			/* Remove the first node from the queue */
			if (N == Last)
				First = Last = 0;
			else
				First = N->OldSuc;
			for (NN = N->CandidateSet; (To = NN->To); NN++) {
				if (To->Mark != mark) {
					To->Mark = mark;
					To->Level = N->Level + 1;
					if (MayBeAddedToFragments(From, To, Alg, info) &&
						(N == From ? (d = NN->Cost) < Min :
						(!Alg->c || (Alg->*(Alg->c))(From, To) < Min)
							&& (d = (Alg->*(Alg->C))(From, To)) < Min)) {
						Min = From->Cost = d;
						/* Randomization */
						if (!Nearest && Alg->Random() % 3 != 0)
							return To;
						Nearest = To;
						MaxLevel = To->Level;
					}
					else if (To->Level < MaxLevel) {
						/* Insert To as the last element of the queue */
						if (Last)
							Last->OldSuc = To;
						else
							First = To;
						Last = To;
						Last->OldSuc = 0;
					}
				}
			}
		}
		return Nearest;
	}

	/*
	 * The NearestInList function returns for a given node, From, the "nearest"
	 * neighbor, To, for which the addition of the edge (From, To) will not make
	 * it impossible to complete a tour.
	 *
	 * The neighbor is determined by searching the list of nodes that are not in
	 * any fragment.
	 */

	static LKHAlg::Node *NearestInList(LKHAlg::Node * From, LKHAlg::Node * First, LKHAlg *Alg, GreedyTourInfo& info)
	{
		LKHAlg::Node *To, *Nearest = 0;
		int Min = std::numeric_limits<int>::max(), d;

		To = First;
		do {
			if (MayBeAddedToFragments(From, To, Alg, info) &&
				(!Alg->c || (Alg->*(Alg->c))(From, To) < Min) && (d = (Alg->*(Alg->C))(From, To)) < Min) {
				Min = From->Cost = d;
				Nearest = To;
			}
		} while ((To = To->OldSuc) != First);
		return Nearest;
	}

	/*
	 * The MayBeAddedToFragments is used to test if the addition of a given edge,
	 * (From, To), makes it impossible to complete a tur.
	 * If the edge may be added, the function returns 1; otherwise 0.
	 */

	static int MayBeAddedToFragments(LKHAlg::Node * From, LKHAlg::Node * To, LKHAlg *Alg, GreedyTourInfo& info)
	{
		return From != To && From->Degree != 2 && To->Degree != 2 &&
			(From->Tail != To || info.EdgesInFragments == Alg->Dimension - 1) &&
			!Alg->Forbidden(From, To);
	}

	/*
	 * The AddEdgeToFragments function adds a given edge, (From, To), to the
	 * current set of fragments.
	 */

	static void AddEdgeToFragments(LKHAlg::Node * From, LKHAlg::Node * To, LKHAlg *Alg, GreedyTourInfo& info)
	{
		LKHAlg::Node *Temp;

		if (!From->Pred)
			From->Pred = To;
		else
			From->Suc = To;
		if (!To->Pred)
			To->Pred = From;
		else
			To->Suc = From;
		From->Degree++;
		To->Degree++;
		Temp = From->Tail;
		Temp->Tail = To->Tail;
		To->Tail->Tail = Temp;
		info.EdgesInFragments++;
		info.Cost += From->Cost - From->Pi - To->Pi;
	}

	/*
	 * The RemoveFromList function removes a given node, N, from the list of
	 * non-fragment nodes. The parameter FIrst is a reference to the first node
	 * of this list.
	 */

	static void RemoveFromList(LKHAlg::Node * N, LKHAlg::Node ** First, LKHAlg *Alg)
	{
		if (*First == N)
			*First = N->OldSuc;
		N->OldPred->OldSuc = N->OldSuc;
		N->OldSuc->OldPred = N->OldPred;
	}

	static int compareX(const void *Na, const void *Nb)
	{
		double x1 = (*(LKHAlg::Node **) Na)->X;
		double y1 = (*(LKHAlg::Node **) Na)->Y;
		double x2 = (*(LKHAlg::Node **) Nb)->X;
		double y2 = (*(LKHAlg::Node **) Nb)->Y;
		return x1 < x2 ? -1 : x1 > x2 ? 1 : y1 < y2 ? -1 : y1 > y2 ? 1 : 0;
	}

	static int compareCost(const void *Na, const void *Nb)
	{
		return (*(LKHAlg::Node **) Na)->Cost - (*(LKHAlg::Node **) Nb)->Cost;
	}
}