#include "./INCLUDE/LKH.h"
#include "./INCLUDE/GeoConversion.h"

/*
 * The SolveKarpSubproblems function attempts to improve a given tour 
 * by means of Karp's partitioning scheme. 
 *
 * The overall region containing the nodes is subdivided into rectangles
 * with SubproblemSize nodes in each rectangle. Each rectangle together 
 * with the given tour induces a subproblem consisting of all nodes inside 
 * the rectangle, and with edges fixed between nodes that are connected 
 * by tour segments whose interior points are outside the rectangle.  
 *  
 * If an improvement is found, the new tour is written to TourFile. 
 * The original tour is given by the SubproblemSuc references of the nodes.
 */
namespace LKH {

	typedef struct KarpInfo {
		LKHAlg::Node** KDTree;
		GainType GlobalBestCost, OldGlobalBestCost;
		int CurrentSubproblem, Subproblems;
	};


	static void KarpPartition(int start, int end, LKHAlg *Alg, KarpInfo& info);
	static void CalculateSubproblems(int start, int end, LKHAlg *Alg, KarpInfo& info);





	void LKHAlg::SolveKarpSubproblems()
	{

		KarpInfo info;
		Node *N;
		double EntryTime = GetTime();

		AllocateStructures();
		ReadPenalties();

		/* Compute upper bound for the original problem */
		info.GlobalBestCost = 0;
		N = FirstNode;
		do {
			if (!Fixed(N, N->SubproblemSuc))
				info.GlobalBestCost += (this->*(this->Distance))(N, N->SubproblemSuc);
			N->Subproblem = 0;
		} while ((N = N->SubproblemSuc) != FirstNode);
		if (TraceLevel >= 1) {
			if (TraceLevel >= 2)
				printff("\n");
			printff("*** Karp partitioning *** [Cost = " GainFormat "]\n",
				info.GlobalBestCost);
		}
		if (WeightType == GEO || WeightType == GEOM ||
			WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS) {
			N = FirstNode;
			do {
				N->Xc = N->X;
				N->Yc = N->Y;
				N->Zc = N->Z;
				if (WeightType == GEO || WeightType == GEO_MEEUS)
					GEO2XYZ(N->Xc, N->Yc, &N->X, &N->Y, &N->Z);
				else
					GEOM2XYZ(N->Xc, N->Yc, &N->X, &N->Y, &N->Z);
			} while ((N = N->SubproblemSuc) != FirstNode);
			CoordType = THREED_COORDS;
		}
		info.KDTree = BuildKDTree(SubproblemSize);
		if (WeightType == GEO || WeightType == GEOM ||
			WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS) {
			N = FirstNode;
			do {
				N->X = N->Xc;
				N->Y = N->Yc;
				N->Z = N->Zc;
			} while ((N = N->SubproblemSuc) != FirstNode);
			CoordType = TWOD_COORDS;
		}

		info.Subproblems = 0;
		CalculateSubproblems(0, Dimension - 1, this,info);
		info.CurrentSubproblem = 0;
		KarpPartition(0, Dimension - 1, this,info);
		free(info.KDTree);
		printff("\nCost = " GainFormat, info.GlobalBestCost);
		if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
			printff(", Gap = %0.4f%%",
				100.0 * (info.GlobalBestCost - Optimum) / Optimum);
		printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
			info.GlobalBestCost < Optimum ? "<" : info.GlobalBestCost ==
			Optimum ? "=" : "");
		if (SubproblemBorders && info.Subproblems > 1)
			SolveSubproblemBorderProblems(info.Subproblems, &info.GlobalBestCost);
	}

	/*
	 * The KarpPartition function subidivides the overall region into
	 * rectangles and attempts to solve the induced subproblems.
	 */

	static void KarpPartition(int start, int end, LKHAlg *Alg, KarpInfo& info)
	{
		if (end - start + 1 <= Alg->SubproblemSize) {
			int i;
			info.CurrentSubproblem++;
			for (i = start; i <= end; i++)
				info.KDTree[i]->Subproblem = info.CurrentSubproblem;
			info.OldGlobalBestCost = info.GlobalBestCost;
			Alg->SolveSubproblem(info.CurrentSubproblem, info.Subproblems, &info.GlobalBestCost);
			if (Alg->SubproblemsCompressed && info.GlobalBestCost == info.OldGlobalBestCost)
				Alg->SolveCompressedSubproblem(info.CurrentSubproblem, info.Subproblems,
					&info.GlobalBestCost);
		}
		else {
			int mid = (start + end) / 2;
			KarpPartition(start, mid, Alg,info);
			KarpPartition(mid + 1, end, Alg,info);
		}
	}

	/*
	 * The CalculateSubproblems function is used to calculate the number of
	 * subproblems (info.Subproblems) created by the KarpPartition function.
	 */

	static void CalculateSubproblems(int start, int end, LKHAlg *Alg, KarpInfo& info)
	{
		if (end - start + 1 <= Alg->SubproblemSize)
			info.Subproblems++;
		else {
			int mid = (start + end) / 2;
			CalculateSubproblems(start, mid, Alg,info);
			CalculateSubproblems(mid + 1, end, Alg,info);
		}
	}
}