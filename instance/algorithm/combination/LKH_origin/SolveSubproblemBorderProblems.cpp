#include "./INCLUDE/LKH.h"
#include "./INCLUDE/GainType.h"
#include "./INCLUDE/GeoConversion.h"

/*
 * Given a tour and a partitioning of the problem into subproblems, the 
 * SolveSubproblemBorderProblems function attempt to improve the tour by 
 * means of a new partitioning given by the borders of these subproblems.
 *
 * For each of the original subproblems a new subproblem is defined by
 * the SubproblemSize points that are closest to its windowing border.
 * These border points together with the given tour induces a subproblem 
 * consisting of all border points, and with edges fixed between points 
 * that are connected by tour segments whose interior points are outside 
 * the border.  
 *  
 * If an improvement is found, the new tour is written to TourFile. 
 * The original tour is given by the SubproblemSuc references of the nodes.
 *
 * The parameter Subproblems specifies the number of subproblems.
 * The parameter GlobalBestCost references the current best cost of the
 * whole problem.
 */
namespace LKH {
	static void MarkBorderPoints(int CurrentSubproblem, LKHAlg *Alg);
	static void QuickSelect(LKHAlg::Node ** A, int n, int k, LKHAlg *Alg);

	void LKHAlg::SolveSubproblemBorderProblems(int Subproblems,
		GainType * GlobalBestCost)
	{
		Node *N;
		GainType OldGlobalBestCost;
		int CurrentSubproblem;
		int *SubproblemSaved;
		double EntryTime = GetTime();

		assert(SubproblemSaved =
			(int *)malloc((DimensionSaved + 1) * sizeof(int)));
		/* Compute upper bound for the original problem */
		N = FirstNode;
		do {
			N->Suc = N->SubproblemSuc;
			N->Suc->Pred = N;
			if (N->Subproblem > Subproblems)
				N->Subproblem -= Subproblems;
			SubproblemSaved[N->Id] = N->Subproblem;
			N->FixedTo1Saved = N->FixedTo2Saved = 0;
			N->SubBestPred = N->SubBestSuc = 0;
		} while ((N = N->SubproblemSuc) != FirstNode);
		if (TraceLevel >= 1)
			printff("\n*** Solve subproblem border problems *** [" GainFormat
				"]\n", *GlobalBestCost);
		for (CurrentSubproblem = 1;
			CurrentSubproblem <= Subproblems; CurrentSubproblem++) {
			MarkBorderPoints(CurrentSubproblem, this);
			OldGlobalBestCost = *GlobalBestCost;
			SolveSubproblem(CurrentSubproblem, Subproblems, GlobalBestCost);
			if (SubproblemsCompressed && *GlobalBestCost == OldGlobalBestCost)
				SolveCompressedSubproblem(CurrentSubproblem, Subproblems,
					GlobalBestCost);
			N = FirstNode;
			do
				N->Subproblem = SubproblemSaved[N->Id];
			while ((N = N->SubproblemSuc) != FirstNode);
		}
		free(SubproblemSaved);
		printff("\nCost = " GainFormat, *GlobalBestCost);
		if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
			printff(", Gap = %0.4f%%",
				100.0 * (*GlobalBestCost - Optimum) / Optimum);
		printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
			*GlobalBestCost < Optimum ? "<" : *GlobalBestCost ==
			Optimum ? "=" : "");
	}

#define Coord(N, axis) (axis == 0 ? (N)->X : axis == 1 ? (N)->Y : (N)->Z)

	/*
	 * The MarkBorderPoints function marks the border points of a given
	 * subproblem (CurrentSubproblem >= 1) by setting their Subproblem value to
	 * CurrentSubproblem. All other points are given a Subproblem value of 0.
	 */

	static void MarkBorderPoints(int CurrentSubproblem, LKHAlg *Alg)
	{
		double Min[3], Max[3];
		int dMin, dMax, d, i, axis, ActualSubproblemSize = 0, Size = 0;
		LKHAlg::Node **A, *N;

		assert(A = (LKHAlg::Node **) malloc(Alg->DimensionSaved * sizeof(LKHAlg::Node *)));
		Min[0] = Min[1] = Min[2] = std::numeric_limits<double>::max();
		Max[0] = Max[1] = Max[2] = -std::numeric_limits<double>::max();
		if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEOM ||
			Alg->WeightType == LKH::GEO_MEEUS || Alg->WeightType == LKH::GEOM_MEEUS) {
			N = Alg->FirstNode;
			do {
				N->Xc = N->X;
				N->Yc = N->Y;
				N->Zc = N->Z;
				if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEO_MEEUS)
					GEO2XYZ(N->Xc, N->Yc, &N->X, &N->Y, &N->Z);
				else
					GEOM2XYZ(N->Xc, N->Yc, &N->X, &N->Y, &N->Z);
			} while ((N = N->SubproblemSuc) != Alg->FirstNode);
			Alg->CoordType = LKH::THREED_COORDS;
		}
		N = Alg->FirstNode;
		do {
			if (N->Subproblem == CurrentSubproblem) {
				for (i = Alg->CoordType == LKH::THREED_COORDS ? 2 : 1; i >= 0; i--) {
					if (Coord(N, i) < Min[i])
						Min[i] = Coord(N, i);
					if (Coord(N, i) > Max[i])
						Max[i] = Coord(N, i);
				}
				ActualSubproblemSize++;
			}
		} while ((N = N->SubproblemSuc) != Alg->FirstNode);
		do {
			if (N->Subproblem == CurrentSubproblem ||
				(N->X >= Min[0] && N->X <= Max[0] &&
					N->Y >= Min[1] && N->Y <= Max[1] &&
					(Alg->CoordType == LKH::TWOD_COORDS ||
					(N->Z >= Min[2] && N->Z <= Max[2])))) {
				N->Rank = std::numeric_limits<int>::max();
				for (i = Alg->CoordType == LKH::THREED_COORDS ? 2 : 1; i >= 0; i--) {
					dMin = (int)(fabs(Coord(N, i) - Min[i]) + 0.5);
					dMax = (int)(fabs(Coord(N, i) - Max[i]) + 0.5);
					d = dMin < dMax ? dMin : dMax;
					if (d < N->Rank)
						N->Rank = d;
				}
			}
			else {
				axis = -1;
				if (Alg->CoordType == LKH::TWOD_COORDS) {
					if (N->X >= Min[0] && N->X <= Max[0])
						axis = 1;
					else if (N->Y >= Min[1] && N->Y <= Max[1])
						axis = 0;
				}
				else if (N->X >= Min[0] && N->X <= Max[0]) {
					if (N->Y >= Min[1] && N->Y <= Max[1])
						axis = 2;
					else if (N->Z >= Min[2] && N->Z <= Max[2])
						axis = 1;
				}
				else if (N->Y >= Min[1] && N->Y <= Max[1] &&
					N->Z >= Min[2] && N->Z <= Max[2])
					axis = 0;
				if (axis != -1) {
					dMin = (int)(fabs(Coord(N, axis) - Min[axis]) + 0.5);
					dMax = (int)(fabs(Coord(N, axis) - Max[axis]) + 0.5);
					N->Rank = dMin < dMax ? dMin : dMax;
				}
				else {
					N->Rank = 0;
					for (i = Alg->CoordType == LKH::THREED_COORDS ? 2 : 1; i >= 0; i--) {
						dMin = (int)(fabs(Coord(N, i) - Min[i]) + 0.5);
						dMax = (int)(fabs(Coord(N, i) - Max[i]) + 0.5);
						d = dMin < dMax ? dMin : dMax;
						if (d > N->Rank)
							N->Rank = d;
					}
				}
			}
			N->Subproblem = 0;
			if (!Alg->SubproblemsCompressed ||
				((N->SubproblemPred != N->SubBestPred ||
					N->SubproblemSuc != N->SubBestSuc) &&
					(N->SubproblemPred != N->SubBestSuc ||
						N->SubproblemSuc != N->SubBestPred)))
				A[Size++] = N;
		} while ((N = N->SubproblemSuc) != Alg->FirstNode);
		if (ActualSubproblemSize > Size)
			ActualSubproblemSize = Size;
		else
			QuickSelect(A, Size, ActualSubproblemSize, Alg);
		for (Size = 0; Size < ActualSubproblemSize; Size++)
			A[Size]->Subproblem = CurrentSubproblem;
		free(A);
		if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEOM ||
			Alg->WeightType == LKH::GEO_MEEUS || Alg->WeightType == LKH::GEOM_MEEUS) {
			N = Alg->FirstNode;
			do {
				N->X = N->Xc;
				N->Y = N->Yc;
				N->Z = N->Zc;
			} while ((N = N->SubproblemSuc) != Alg->FirstNode);
			Alg->CoordType = LKH::TWOD_COORDS;
		}
	}

#define SWAP(a, b) { temp = (a); (a) = (b); (b) = temp; }

	/*
	 * The QuickSelect function rearranges the array A[0..n-1] such that
	 * A[0..k-1]->Rank are less than or equal to A[k]->Rank.
	 */

	static void QuickSelect(LKHAlg::Node ** A, int n, int k, LKHAlg *Alg)
	{
		int i, j, l = 0, r = n - 1, mid, pivot;
		LKHAlg::Node *temp, *v;

		for (;;) {
			if (r <= l + 1) {
				if (r == l + 1 && A[r]->Rank < A[l]->Rank)
					SWAP(A[l], A[r]);
				return;
			}
			else {
				mid = (l + r) / 2;
				SWAP(A[mid], A[l + 1]);
				if (A[l]->Rank > A[r]->Rank)
					SWAP(A[l], A[r]);
				if (A[l + 1]->Rank > A[r]->Rank)
					SWAP(A[l + 1], A[r]);
				if (A[l]->Rank > A[l + 1]->Rank)
					SWAP(A[l], A[l + 1]);
				i = l + 1;
				j = r;
				v = A[l + 1];
				pivot = v->Rank;
				for (;;) {
					do
						i++;
					while (A[i]->Rank < pivot);
					do
						j--;
					while (A[j]->Rank > pivot);
					if (j < i)
						break;
					SWAP(A[i], A[j]);
				}
				A[l + 1] = A[j];
				A[j] = v;
				if (j >= k)
					r = j - 1;
				if (j <= k)
					l = i;
			}
		}
	}
}