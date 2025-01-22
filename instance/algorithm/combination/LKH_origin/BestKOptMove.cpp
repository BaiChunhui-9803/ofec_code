#include "./INCLUDE/Segment.h"
#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Sequence.h"

#include <cstring>

/*
 * The BestKOptMove function makes edge exchanges. If possible, it makes a 
 * r-opt move (r >= 2) that improves the tour. Otherwise, it makes the most 
 * promising sequential K-opt move that fulfils the positive gain criterion. 
 * To prevent an infinity chain of moves the last edge in a K-opt move must
 * not previously have been included in the chain. 
 *
 * The edge (t[1],t[2]) is the first edge to be exchanged. G0 is a pointer to 
 * the accumulated gain.
 *
 * In case a K-opt move is found that improves the tour, the improvement of 
 * the cost is made available to the caller through the parameter Gain. 
 * If *Gain > 0, an improvement of the current tour has been found. In this
 * case the function returns 0.
 *
 * Otherwise, the best K-opt move is made, and a pointer to the node that was 
 * connected to t[1] (in order to close the tour) is returned. The new 
 * accumulated gain is made available to the caller through the parameter G0. 
 *
 * The function is called from the LinKernighan function. 
 */
//static thread_local unique_ptr<GainType> BestG2;
namespace LKH {
	static thread_local GainType BestG2 = 0;

	static GainType BestKOptMoveRec(int k, GainType G0, LKHAlg *Alg);

	LKHAlg::Node *LKHAlg::BestKOptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
	{
		/*if (!BestG2.get())
			BestG2.reset(new GainType(0));*/
		OldSwaps = Swaps;
		PenaltyGain = 0;
		K = Swaps == 0 ? MoveType : SubsequentMoveType;
		*Gain = 0;
		t[1] = t1;
		t[2] = t2;
		T[2 * K] = 0;
		BestG2 = std::numeric_limits<GainType>::min();

		/*
		 * Determine (T[3],T[4], ..., T[2K]) = (t[3],t[4], ..., t[2K])
		 * such that
		 *
		 *     G[2 * K] = *G0 - C(t[2],T[3]) + C(T[3],T[4])
		 *                    - C(T[4],T[5]) + C(T[5],T[6])
		 *                      ...
		 *                    - C(T[2K-3],T[2K-2]) + C(T[2K-1],T[2K])
		 *
		 * is maximum, and (T[2K-1],T[2K]) has not previously been included.
		 * If during this process a legal move with *Gain > 0 is found, then
		 * make the move and exit BestKOptMove immediately.
		 */

		MarkDeleted(t1, t2);
		*Gain = BestKOptMoveRec(2, *G0, this);
		UnmarkDeleted(t1, t2);

		if (PenaltyGain <= 0 && *Gain <= 0 && T[2 * K]) {
			int i;
			std::memcpy(t + 1, T + 1, 2 * K * sizeof(Node *));
			for (i = 2; i < 2 * K; i += 2)
				incl[incl[i] = i + 1] = i;
			incl[incl[1] = 2 * K] = 1;
			MakeKOptMove(K);
			for (i = 1; i < 2 * K; i += 2)
				Exclude(T[i], T[i + 1]);
			*G0 = BestG2;
			return T[2 * K];
		}
		return 0;
	}

	static GainType BestKOptMoveRec(int k, GainType G0, LKHAlg *Alg)
	{
		LKHAlg::Candidate *Nt2;
		LKHAlg::Node *t1, *t2, *t3, *t4, *SUCt1;
		GainType G1, G2, G3, Gain;
		int X4, i;
		int Breadth2 = 0;

		t1 = t[1];
		t2 = t[i = 2 * k - 2];
		incl[incl[i] = i + 1] = i;
		incl[incl[1] = i + 2] = 1;
		SUCt1 = (Alg->Reversed == (t1)->Parent->Reversed ? (t1)->Suc : (t1)->Pred);
		/* Choose (t2,t3) as a candidate edge emanating from t2 */
		for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
			if (t3 == t2->Pred || t3 == t2->Suc ||
				((G1 = G0 - Nt2->Cost) <= 0 &&
					Alg->GainCriterionUsed && (k > 2 ||
					(Alg->ProblemType != LKH::HCP
						&& Alg->ProblemType != LKH::HPP)))
				|| Added(t2, t3))
				continue;
			if (++Breadth2 > Alg->MaxBreadth)
				break;
			MarkAdded(t2, t3);
			t[2 * k - 1] = t3;
			G[2 * k - 2] = G1 + t3->Pi;
			/* Choose t4 as one of t3's two neighbors on the tour */
			for (X4 = 1; X4 <= 2; X4++) {
				t4 = X4 == 1 ? (Alg->Reversed == (t3)->Parent->Reversed ? (t3)->Pred : (t3)->Suc) : (Alg->Reversed == (t3)->Parent->Reversed ? (t3)->Suc : (t3)->Pred);
				if ((LKH::LKHAlg::Fixed(t3, t4) || Alg->IsCommonEdge(t3, t4)) || Deleted(t3, t4))
					continue;
				t[2 * k] = t4;
				G2 = G1 + (Alg->*(Alg->C))(t3, t4);
				G3 = std::numeric_limits<GainType>::min();
				if (t4 != t1 && !Alg->Forbidden(t4, t1) && !Added(t4, t1) &&
					(Alg->CurrentPenalty > 0 || Alg->TSPTW_Makespan ||
					((!Alg->c || G2 - (Alg->*(Alg->c))(t4, t1) > 0) &&
						(G3 = G2 - (Alg->*(Alg->C))(t4, t1)) > 0)) && FeasibleKOptMove(k)) {
					if (Alg->CurrentPenalty > 0 || Alg->TSPTW_Makespan)
						G3 = G2 - (Alg->*(Alg->C))(t4, t1);
					if (Alg->CurrentPenalty || Alg->TSPTW_Makespan || G3 > 0) {
						Alg->MakeKOptMove(k);
						if (Alg->Improvement(&G3, t1, SUCt1)) {
							UnmarkAdded(t2, t3);
							return G3;
						}
					}
				}
				if (Alg->Backtracking && !Alg->Excludable(t3, t4))
					continue;
				MarkDeleted(t3, t4);
				G[2 * k - 1] = G2 - t4->Pi;
				if (k < (K)) {
					Gain = BestKOptMoveRec(k + 1, G2, Alg);
					if (Alg->PenaltyGain > 0 || Gain > 0) {
						UnmarkAdded(t2, t3);
						UnmarkDeleted(t3, t4);
						return Gain;
					}
					incl[incl[1] = 2 * k] = 1;
				}
				if (t4 != t1 && !Alg->Forbidden(t4, t1) &&
					k + 1 < Alg->NonsequentialMoveType &&
					Alg->PatchingC >= 2 && Alg->PatchingA >= 1 &&
					(Alg->Swaps == 0 || Alg->SubsequentPatching)) {
					if (G3 == std::numeric_limits<GainType>::min())
						G3 = G2 - (Alg->*(Alg->C))(t4, t1);
					if ((Alg->PatchingCRestricted ? G3 > 0 && Alg->IsCandidate(t4, t1) :
						Alg->PatchingCExtended ? G3 > 0
						|| Alg->IsCandidate(t4, t1) : G3 > 0)) {
						Gain = Alg->PatchCycles(k, G3);
						if (Alg->PenaltyGain > 0 || Gain > 0) {
							UnmarkAdded(t2, t3);
							UnmarkDeleted(t3, t4);
							return Gain;
						}
					}
				}
				UnmarkDeleted(t3, t4);
				if (k == (K) && t4 != t1 && t3 != t1 && G3 <= 0 &&
					!Added(t4, t1) &&
					(!Alg->GainCriterionUsed || G2 - Alg->Precision >= t4->Cost)) {
					if (!Alg->Backtracking || Alg->Swaps > 0) {
						if ((G2 > BestG2 ||
							(G2 == BestG2 && !LKH::LKHAlg::Near(t3, t4) &&
								LKH::LKHAlg::Near(T[2 * (K)-1], T[2 * (K)]))) &&
							Alg->Swaps < Alg->MaxSwaps &&
							Alg->Excludable(t3, t4) && !LKH::LKHAlg::InInputTour(t3, t4)) {
							if (Alg->RestrictedSearch && (K) > 2 &&
								Alg->ProblemType != LKH::HCP && Alg->ProblemType != LKH::HPP) {
								/* Ignore the move if the gain does not vary */
								G[0] = G[2 * (K)-2];
								G[1] = G[2 * (K)-1];
								for (i = 2 * (K)-3; i >= 2; i--)
									if (G[i] != G[i % 2])
										break;
								if (i < 2)
									continue;
							}
							if (FeasibleKOptMove(K)) {
								BestG2 = G2;
								std::memcpy(T + 1, t + 1, 2 * (K) * sizeof(LKHAlg::Node *));
							}
						}
					}
					else if (Alg->MaxSwaps > 0 && FeasibleKOptMove(K)) {
						LKHAlg::Node *SUCt1 = (Alg->Reversed == (t1)->Parent->Reversed ? (t1)->Suc : (t1)->Pred);
						Alg->MakeKOptMove(K);
						for (i = 1; i < 2 * k; i += 2) {
							Alg->Exclude(t[i], t[i + 1]);
							UnmarkDeleted(t[i], t[i + 1]);
						}
						for (i = 2; i < 2 * k; i += 2)
							UnmarkAdded(t[i], t[i + 1]);
						std::memcpy(tSaved + 1, t + 1, 2 * k * sizeof(LKHAlg::Node *));
						while ((t4 = (Alg->*(Alg->BestSubsequentMove))(t1, t4, &G2, &Gain)));
						if (Alg->PenaltyGain > 0 || Gain > 0) {
							UnmarkAdded(t2, t3);
							return Gain;
						}
						Alg->OldSwaps = 0;
						Alg->RestoreTour();
						K = k;
						memcpy(t + 1, tSaved + 1, 2 * (K) * sizeof(LKHAlg::Node *));
						for (i = 1; i < 2 * (K)-2; i += 2)
							MarkDeleted(t[i], t[i + 1]);
						for (i = 2; i < 2 * (K); i += 2)
							MarkAdded(t[i], t[i + 1]);
						for (i = 2; i < 2 * (K); i += 2)
							incl[incl[i] = i + 1] = i;
						incl[incl[1] = 2 * (K)] = 1;
						if (SUCt1 != (Alg->Reversed == (t1)->Parent->Reversed ? (t1)->Suc : (t1)->Pred))
							Alg->Reversed ^= 1;
						T[2 * (K)] = 0;
					}
				}
			}
			UnmarkAdded(t2, t3);
			if (t3 == t1)
				continue;
			/* Try to delete an added edge, (_,t3) or (t3,_) */
			for (i = 2 * k - 4; i >= 2; i--) {
				if (t3 == t[i]) {
					t4 = t[i ^ 1];
					if (t4 == t1 || Alg->Forbidden(t4, t1) || (LKH::LKHAlg::Fixed(t3, t4) || Alg->IsCommonEdge(t3, t4))
						|| Added(t4, t1))
						continue;
					G2 = G1 + (Alg->*(Alg->C))(t3, t4);
					if ((!(Alg->c) || G2 - (Alg->*(Alg->c))(t4, t1) > 0)
						&& (Gain = G2 - (Alg->*(Alg->C))(t4, t1)) > 0) {
						incl[incl[i ^ 1] = 1] = i ^ 1;
						incl[incl[i] = 2 * k - 2] = i;
						if (FeasibleKOptMove(k - 1)) {
							Gain = G1 + (Alg->*(Alg->C))(t3, t4) - (Alg->*(Alg->C))(t4, t1);
							Alg->MakeKOptMove(k - 1);
							if (Alg->Improvement(&Gain, t1, SUCt1))
								return Gain;
						}
						incl[incl[i ^ 1] = i] = i ^ 1;
					}
				}
			}
			incl[1] = 2 * k;
			incl[2 * k - 2] = 2 * k - 1;
		}
		return 0;
	}
}