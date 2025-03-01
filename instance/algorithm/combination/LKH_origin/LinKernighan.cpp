#include "./INCLUDE/Segment.h"
#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Hashing.h"
#include "./INCLUDE/Sequence.h"
#include "./INCLUDE/BIT.h"

/*
 * The LinKernighan function seeks to improve a tour by sequential
 * and non-sequential edge exchanges.
 *
 * The function returns the cost of the resulting tour.
 */
namespace LKH {
	GainType LKHAlg::LinKernighan()
	{
		GainType Cost, Gain, G0;
		int X2, i, it = 0;
		Node *t1 = 0, *t2 = 0, *SUCt1 = 0;
		Candidate *Nt1 = 0;
		Segment *S = 0;
		SSegment *SS = 0;
		double EntryTime = GetTime();

		Cost = 0;
		Reversed = 0;
		S = FirstSegment;
		i = 0;
		do {
			S->Size = 0;
			S->Rank = ++i;
			S->Reversed = 0;
			S->First = S->Last = 0;
		} while ((S = S->Suc) != FirstSegment);
		SS = FirstSSegment;
		i = 0;
		do {
			SS->Size = 0;
			SS->Rank = ++i;
			SS->Reversed = 0;
			SS->First = SS->Last = 0;
		} while ((SS = SS->Suc) != FirstSSegment);

		FirstActive = LastActive = 0;
		Swaps = 0;

		/* Compute the cost of the initial tour, Cost.
		   Compute the corresponding hash value, Hash.
		   Initialize the segment list.
		   Make all nodes "active" (so that they can be used as t1). */

		Cost = 0;
		Hash = 0;
		i = 0;
		t1 = FirstNode;
		do {
			t2 = t1->OldSuc = t1->Suc;
			t1->OldPred = t1->Pred;
			t1->Rank = ++i;
			Cost += (t1->SucCost = t2->PredCost = (this->*C)(t1, t2)) - t1->Pi - t2->Pi;
			Hash ^= Rand[t1->Id] * Rand[t2->Id];
			t1->Cost = std::numeric_limits<int>::max();
			for (Nt1 = t1->CandidateSet; (t2 = Nt1->To); Nt1++)
				if (t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
					t1->Cost = Nt1->Cost;
			t1->Parent = S;
			S->Size++;
			if (S->Size == 1)
				S->First = t1;
			S->Last = t1;
			if (SS->Size == 0)
				SS->First = S;
			S->Parent = SS;
			SS->Last = S;
			if (S->Size == GroupSize) {
				S = S->Suc;
				SS->Size++;
				if (SS->Size == SGroupSize)
					SS = SS->Suc;
			}
			t1->OldPredExcluded = t1->OldSucExcluded = 0;
			t1->Next = 0;
			if (KickType == 0 || Kicks == 0 || Trial == 1 ||
				!InBestTour(t1, t1->Pred) || !InBestTour(t1, t1->Suc))
				Activate(t1);
		} while ((t1 = t1->Suc) != FirstNode);
		if (S->Size < GroupSize)
			SS->Size++;
		Cost /= Precision;
		if (TSPTW_Makespan)
			Cost = TSPTW_CurrentMakespanCost = TSPTW_MakespanCost();
		CurrentPenalty = std::numeric_limits<GainType>::max();
		CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
	/*	if (TraceLevel >= 3 ||
			(TraceLevel == 2 &&
			(CurrentPenalty < BetterPenalty ||
				(CurrentPenalty == BetterPenalty && Cost < BetterCost))))*/
			//StatusReport(Cost, EntryTime, "");
			PredSucCostAvailable = 1;
		BIT_Update(this);

		/* Loop as long as improvements are found */
		do {
			/* Choose t1 as the first "active" node */
			while ((t1 = RemoveFirstActive())) {
				/* t1 is now "passive" */
				SUCt1 = SUC(t1);
	/*			if ((TraceLevel >= 3 || (TraceLevel == 2 && Trial == 1)) &&
					++it % (Dimension >= 100000 ? 10000 :
						Dimension >= 10000 ? 1000 : 100) == 0)*/
					/* printff("#%d: Time = %0.2f sec.\n",
							 it, fabs(GetTime() - EntryTime));*/
							 /* Choose t2 as one of t1's two neighbors on the tour */
					for (X2 = 1; X2 <= 2; X2++) {
						t2 = X2 == 1 ? PRED(t1) : SUCt1;
						if (FixedOrCommon(t1, t2) ||
							(RestrictedSearch && Near(t1, t2) &&
							(Trial == 1 ||
								(Trial > BackboneTrials &&
								(KickType == 0 || Kicks == 0)))))
							continue;
						G0 = (this->*C)(t1, t2);
						OldSwaps = Swaps = 0;
						PenaltyGain = Gain = 0;
						/* Try to find a tour-improving chain of moves */
						do
							t2 = Swaps == 0 ? (this->*BestMove)(t1, t2, &G0, &Gain) :
							(this->*BestSubsequentMove)(t1, t2, &G0, &Gain);
						while (t2);
						if (PenaltyGain > 0 || Gain > 0) {
							/* An improvement has been found */
#ifdef HAVE_LONG_LONG
							assert(Gain % Precision == 0);
#else
							assert(fmod(Gain, Precision) == 0);
#endif
							Cost -= Gain / Precision;
							CurrentPenalty -= PenaltyGain;
							StoreTour();
							TSPTW_CurrentMakespanCost = Cost;
				/*			if (TraceLevel >= 3 ||
								(TraceLevel == 2 &&
								(CurrentPenalty < BetterPenalty ||
									(CurrentPenalty == BetterPenalty &&
										Cost < BetterCost))))*/
								//StatusReport(Cost, EntryTime, "");
								if (HashSearch(HTable, Hash, Cost))
									goto End_LinKernighan;
							/* Make t1 "active" again */
							Activate(t1);
							OldSwaps = 0;
							break;
						}
						OldSwaps = 0;
						RestoreTour();
					}
			}
			if (HashSearch(HTable, Hash, Cost))
				goto End_LinKernighan;
			HashInsert(HTable, Hash, Cost);
			/* Try to find improvements using non-sequential 4/5-opt moves */
			CurrentPenalty = std::numeric_limits<GainType>::max();
			CurrentPenalty = Penalty ? (this->*Penalty)() : 0;
			PenaltyGain = 0;
			if (Gain23Used && ((Gain = Gain23()) > 0 || PenaltyGain > 0)) {
				/* An improvement has been found */
#ifdef HAVE_LONG_LONG
				assert(Gain % Precision == 0);
#else
				assert(fmod(Gain, Precision) == 0);
#endif
				Cost -= Gain / Precision;
				CurrentPenalty -= PenaltyGain;
				TSPTW_CurrentMakespanCost = Cost;
				StoreTour();
				//if (TraceLevel >= 3 ||
				//	(TraceLevel == 2 &&
				//	(CurrentPenalty < BetterPenalty ||
				//		(CurrentPenalty == BetterPenalty && Cost < BetterCost))))
				//	//StatusReport(Cost, EntryTime, "+ ");
					if (HashSearch(HTable, Hash, Cost))
						goto End_LinKernighan;
			}
		} while (PenaltyGain > 0 || Gain > 0);
	End_LinKernighan:
		PredSucCostAvailable = 0;
		NormalizeNodeList();
		NormalizeSegmentList();
		return Cost;
	}
}