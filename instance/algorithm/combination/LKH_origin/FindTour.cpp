#include "./INCLUDE/LKH.h"

/*
 * After the candidate set has been created the FindTour function is called
 * a predetermined number of times (Runs).
 *
 * FindTour performs a number of trials, where in each trial it attempts
 * to improve a chosen initial tour using the modified Lin-Kernighan edge
 * exchange heuristics.
 *
 * Each time a better tour is found, the tour is recorded, and the candidates
 * are reorderded by the AdjustCandidateSet function. Precedence is given to
 * edges that are common to two currently best tours. The candidate set is
 * extended with those tour edges that are not present in the current set.
 * The original candidate set is re-established at exit from FindTour.
 */
namespace LKH {
	static void SwapCandidateSets(LKHAlg *Alg);
	static thread_local GainType OrdinalTourCost;




	GainType LKHAlg::FindTourVer2() {
		GainType Cost;
		Node* t;
		int i;
		double EntryTime = GetTime();

		t = FirstNode;
		do
			t->OldPred = t->OldSuc = t->NextBestSuc = t->BestSuc = 0;
		while ((t = t->Suc) != FirstNode);
		if (Run == 1 && Dimension == DimensionSaved) {
			OrdinalTourCost = 0;
			for (i = 1; i < Dimension; i++)
				OrdinalTourCost += (this->*C)(&NodeSet[i], &NodeSet[i + 1])
				- NodeSet[i].Pi - NodeSet[i + 1].Pi;
			OrdinalTourCost += (this->*C)(&NodeSet[Dimension], &NodeSet[1])
				- NodeSet[Dimension].Pi - NodeSet[1].Pi;
			OrdinalTourCost /= Precision;
		}
		BetterCost = std::numeric_limits<GainType>::max();
		BetterPenalty = CurrentPenalty = std::numeric_limits<GainType>::max();
		if (MaxTrials > 0)
			HashInitialize(HTable);
		else {
			Trial = 1;
			ChooseInitialTour();
			CurrentPenalty = std::numeric_limits<GainType>::max();
			CurrentPenalty = BetterPenalty = Penalty ? (this->*Penalty)() : 0;
		}
		for (Trial = 1; Trial <= MaxTrials; Trial++) {
			//if (GetTime() - EntryTime >= TimeLimit ||
			//	GetTime() - StartTime >= TotalTimeLimit) {
			//	//if (TraceLevel >= 1)
			//	//	printff("*** Time limit exceeded ***\n");
			//	Trial--;
			//	break;
			//}
			/* Choose FirstNode at random */
			if (Dimension == DimensionSaved)
				FirstNode = &NodeSet[1 + Random() % Dimension];
			else
				for (i = Random() % Dimension; i > 0; i--)
					FirstNode = FirstNode->Suc;
			ChooseInitialTour();
			if ((ProblemType == SOP || ProblemType == M1_PDTSP) &&
				(InitialTourAlgorithm != SOP_ALG || Trial > 1))
				SOP_RepairTour();
			Cost = LinKernighanLocalSearchVer2();
			if (GetTime() - EntryTime < TimeLimit /*&&
				GetTime() - StartTime < TotalTimeLimit*/) {
				if (FirstNode->BestSuc && !TSPTW_Makespan) {
					/* Merge tour with current best tour */
					t = FirstNode;
					while ((t = t->Next = t->BestSuc) != FirstNode);
					Cost = (this->*MergeWithTour)();
				}
				if (Dimension == DimensionSaved && Cost >= OrdinalTourCost &&
					BetterCost > OrdinalTourCost && !TSPTW_Makespan) {
					/* Merge tour with ordinal tour */
					for (i = 1; i < Dimension; i++)
						NodeSet[i].Next = &NodeSet[i + 1];
					NodeSet[Dimension].Next = &NodeSet[1];
					Cost = (this->*MergeWithTour)();
				}
			}
			if (CurrentPenalty < BetterPenalty ||
				(CurrentPenalty == BetterPenalty && Cost < BetterCost)) {
				if (TraceLevel >= 1) {
					//printff("* %d: ", Trial);
					//StatusReport(Cost, EntryTime, "");
				}
				BetterCost = Cost;
				BetterPenalty = CurrentPenalty;
				RecordBetterTour();
				//if ((BetterPenalty < BestPenalty ||
				//	(BetterPenalty == BestPenalty && BetterCost < BestCost)) &&
				//	SubproblemSize == 0)
				////	WriteTour(OutputTourFileName, BetterTour, BetterCost);
				//if (StopAtOptimum) {
				//	if (!Penalty ||
				//		(ProblemType != CCVRP &&
				//			ProblemType != MLP &&
				//			ProblemType != TRP &&
				//			Penalty != Penalty_MTSP_MINMAX &&
				//			Penalty != Penalty_MTSP_MINMAX_SIZE) ?
				//		CurrentPenalty == 0 && Cost == Optimum :
				//		CurrentPenalty == Optimum)
				//		break;
				//}
				AdjustCandidateSet();
				HashInitialize(HTable);
				HashInsert(HTable, Hash, Cost);
			}
			else if (TraceLevel >= 2) {
				//printff("  %d: ", Trial);
				//StatusReport(Cost, EntryTime, "");
			}
			/* Record backbones if wanted */
			if (Trial <= BackboneTrials && BackboneTrials < MaxTrials) {
				SwapCandidateSets(this);
				AdjustCandidateSet();
				if (Trial == BackboneTrials) {
	/*				if (TraceLevel >= 1) {
						printff("# %d: Backbone candidates ->\n", Trial);
						CandidateReport();
					}*/
				}
				else
					SwapCandidateSets(this);
			}
		}
		if (BackboneTrials > 0 && BackboneTrials < MaxTrials) {
			if (Trial > BackboneTrials ||
				(Trial == BackboneTrials &&
					(!StopAtOptimum || BetterCost != Optimum)))
				SwapCandidateSets(this);
			t = FirstNode;
			do {
				free(t->BackboneCandidateSet);
				t->BackboneCandidateSet = 0;
			} while ((t = t->Suc) != FirstNode);
		}
		t = FirstNode;
		if (Norm == 0 || MaxTrials == 0 || !t->BestSuc) {
			do
				t = t->BestSuc = t->Suc;
			while (t != FirstNode);
		}
		Hash = 0;
		do {
			(t->Suc = t->BestSuc)->Pred = t;
			Hash ^= Rand[t->Id] * Rand[t->Suc->Id];
		} while ((t = t->BestSuc) != FirstNode);
		if (Trial > MaxTrials)
			Trial = MaxTrials;
		ResetCandidateSet();
		CurrentPenalty = BetterPenalty;
		return BetterCost;
	}


	GainType LKHAlg::FindTour()
	{
		GainType Cost;
		Node *t = 0;
		int i;
		double EntryTime = GetTime();

		t = FirstNode;
		do
			t->OldPred = t->OldSuc = t->NextBestSuc = t->BestSuc = 0;
		while ((t = t->Suc) != FirstNode);
		if (Run == 1 && Dimension == DimensionSaved) {
			OrdinalTourCost = 0;
			for (i = 1; i < Dimension; i++)
				OrdinalTourCost += (this->*C)(&NodeSet[i], &NodeSet[i + 1])
				- NodeSet[i].Pi - NodeSet[i + 1].Pi;
			OrdinalTourCost += (this->*C)(&NodeSet[Dimension], &NodeSet[1])
				- NodeSet[Dimension].Pi - NodeSet[1].Pi;
			OrdinalTourCost /= Precision;
		}
		BetterCost = std::numeric_limits<GainType>::max();
		BetterPenalty = CurrentPenalty = std::numeric_limits<GainType>::max();
		if (MaxTrials > 0)
			HashInitialize(HTable);
		else {
			Trial = 1;
			ChooseInitialTour();
			CurrentPenalty = std::numeric_limits<GainType>::max();
			CurrentPenalty = BetterPenalty = Penalty ? (this->*Penalty)() : 0;
		}
		for (Trial = 1; Trial <= MaxTrials; Trial++) {
			if (GetTime() - EntryTime >= TimeLimit) {
				/*if (TraceLevel >= 1)
					printff("*** Time limit exceeded ***\n");*/
				break;
			}
			/* Choose FirstNode at random */
			if (Dimension == DimensionSaved)
				FirstNode = &NodeSet[1 + Random() % Dimension];
			else
				for (i = Random() % Dimension; i > 0; i--)
					FirstNode = FirstNode->Suc;
			ChooseInitialTour();


			//if (m_recordedGrails) {
			//	std::vector<int> curSol;
			//	RecordCurTour(curSol);
			//	m_trails.push_back(curSol);
			//}


			if ((ProblemType == SOP || ProblemType == M1_PDTSP) &&
				InitialTourAlgorithm != SOP_ALG)
				SOP_RepairTour();
			Cost = LinKernighan();

			if (m_recordedGrails) {
				std::vector<int> curSol;
				RecordCurTour(curSol);
				m_trails.push_back(curSol);
			}
			if (FirstNode->BestSuc && !TSPTW_Makespan) {
				/* Merge tour with current best tour */
				t = FirstNode;
				while ((t = t->Next = t->BestSuc) != FirstNode);
				Cost = (this->*MergeWithTour)();

				if (m_recordedGrails) {
					std::vector<int> curSol;
					RecordCurTour(curSol);
					m_trails.push_back(curSol);
				}
			}
			if (Dimension == DimensionSaved && Cost >= OrdinalTourCost &&
				BetterCost > OrdinalTourCost && !TSPTW_Makespan) {
				/* Merge tour with ordinal tour */
				for (i = 1; i < Dimension; i++)
					NodeSet[i].Next = &NodeSet[i + 1];
				NodeSet[Dimension].Next = &NodeSet[1];
				Cost = (this->*MergeWithTour)();


				if (m_recordedGrails) {
					std::vector<int> curSol;
					RecordCurTour(curSol);
					m_trails.push_back(curSol);
				}
			}
			if (CurrentPenalty < BetterPenalty ||
				(CurrentPenalty == BetterPenalty && Cost < BetterCost)) {
				if (TraceLevel >= 1) {
					//printff("* %d: ", Trial);
					//StatusReport(Cost, EntryTime, "");
				}
				BetterCost = Cost;
				BetterPenalty = CurrentPenalty;
				RecordBetterTour();
				if (BetterPenalty < BestPenalty ||
					(BetterPenalty == BestPenalty && BetterCost < BestCost))
					WriteTour(OutputTourFileName, BetterTour, BetterCost);
				if (StopAtOptimum) {
					if (ProblemType != CCVRP && ProblemType != TRP &&
						MTSPObjective != MINMAX &&
						MTSPObjective != MINMAX_SIZE ?
						CurrentPenalty == 0 && Cost == Optimum :
						CurrentPenalty == Optimum)
						break;
				}
				AdjustCandidateSet();
				HashInitialize(HTable);
				HashInsert(HTable, Hash, Cost);
			}
			else if (TraceLevel >= 2) {
				//printff("  %d: ", Trial);
				//StatusReport(Cost, EntryTime, "");
			}
			/* Record backbones if wanted */
			if (Trial <= BackboneTrials && BackboneTrials < MaxTrials) {
				SwapCandidateSets(this);
				AdjustCandidateSet();
				if (Trial == BackboneTrials) {
					if (TraceLevel >= 1) {
						//printff("# %d: Backbone candidates ->\n", Trial);
						CandidateReport();
					}
				}
				else
					SwapCandidateSets(this);
			}
		}
		if (BackboneTrials > 0 && BackboneTrials < MaxTrials) {
			if (Trial > BackboneTrials ||
				(Trial == BackboneTrials &&
				(!StopAtOptimum || BetterCost != Optimum)))
				SwapCandidateSets(this);
			t = FirstNode;
			do {
				free(t->BackboneCandidateSet);
				t->BackboneCandidateSet = 0;
			} while ((t = t->Suc) != FirstNode);
		}
		t = FirstNode;
		if (Norm == 0) {
			do
				t = t->BestSuc = t->Suc;
			while (t != FirstNode);
		}
		Hash = 0;
		do {
			(t->Suc = t->BestSuc)->Pred = t;
			Hash ^= Rand[t->Id] * Rand[t->Suc->Id];
		} while ((t = t->BestSuc) != FirstNode);
		if (Trial > MaxTrials)
			Trial = MaxTrials;
		ResetCandidateSet();
		CurrentPenalty = BetterPenalty;
		return BetterCost;
	}

	/*
	 * The SwapCandidateSets function swaps the normal and backbone candidate sets.
	 */

	static void SwapCandidateSets(LKHAlg *Alg)
	{
		LKHAlg::Node *t = Alg->FirstNode;
		do {
			LKHAlg::Candidate *Temp = t->CandidateSet;
			t->CandidateSet = t->BackboneCandidateSet;
			t->BackboneCandidateSet = Temp;
		} while ((t = t->Suc) != Alg->FirstNode);
	}
}