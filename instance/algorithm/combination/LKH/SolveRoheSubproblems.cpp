#include "./INCLUDE/LKH.h"
#include "./INCLUDE/GeoConversion.h"

/*
 * The SolveRoheSubproblems function attempts to improve a given tour by means 
 * of Rohe's random rectangle/box partitioning scheme. If an improvement is 
 * found, the new tour is written to TourFile.
 * 
 * The original tour is given by the SubproblemSuc references of the nodes. 
 * The size of each segment is SubproblemSize.
 *
 * The algorithm is described on page 81 in
 *  
 *   A Rohe, 
 *   Parallele Heuristiken fur sehr grosse Traveling Salesman Probleme, 
 *   Diplomarbeit, Research Institute for Discrete Mathematics, 
 *   Universitat Bonn (1997).
 */
namespace LKH {
#define PRANDMAX std::numeric_limits<int>::max()


	typedef struct RoheSubInfo {
		int Size;
		LKHAlg::Node** KDTree;
	};


	static void WindowSize(double XMin, double XMax, double YMin, double YMax,
		double ZMin, double ZMax, int start, int end, LKHAlg *Alg, RoheSubInfo& info);
	static void MakeSubproblem(double XMin, double XMax, double YMin,
		double YMax, double ZMin, double ZMax,
		int Subproblem, int start, int end, LKHAlg *Alg, RoheSubInfo& info);


	void LKHAlg::SolveRoheSubproblems()
	{
		Node *N;
		int CurrentSubproblem, Subproblems, Remaining, i;
		GainType GlobalBestCost, OldGlobalBestCost;
		double XMin, XMax, YMin, YMax, ZMin, ZMax, DX, DY, DZ, CLow, CMid,
			CHigh;
		double EntryTime = GetTime();

		AllocateStructures();
		ReadPenalties();
		Subproblems = 0;

		RoheSubInfo info;

		/* Compute upper bound for the original problem */
		GlobalBestCost = 0;
		N = FirstNode;
		do {
			if (!Fixed(N, N->SubproblemSuc))
				GlobalBestCost += (this->*(this->Distance))(N, N->SubproblemSuc);
			N->Subproblem = 0;
		} while ((N = N->SubproblemSuc) != FirstNode);
		if (TraceLevel >= 1) {
			if (TraceLevel >= 2)
				printff("\n");
			printff("*** Rohe partitioning *** [Cost = " GainFormat "]\n",
				GlobalBestCost);
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
		N = FirstNode;
		XMin = XMax = N->X;
		YMin = YMax = N->Y;
		ZMin = ZMax = N->Z;
		while ((N = N->SubproblemSuc) != FirstNode) {
			if (N->X < XMin)
				XMin = N->X;
			else if (N->X > XMax)
				XMax = N->X;
			if (N->Y < YMin)
				YMin = N->Y;
			else if (N->Y > YMax)
				YMax = N->Y;
			if (N->Z < ZMin)
				ZMin = N->Z;
			else if (N->Z > ZMax)
				ZMax = N->Z;
		}
		info.KDTree = BuildKDTree(SubproblemSize);
		Remaining = Dimension;
		while (Remaining > SubproblemSize) {
			N = FirstNode;
			i = Random() % Remaining;
			while (i--)
				N = N->Suc;
			DX = (0.5 + 0.5 * Random() / (PRANDMAX - 1)) * (XMax - XMin);
			DY = (0.5 + 0.5 * Random() / (PRANDMAX - 1)) * (YMax - YMin);
			DZ = (0.5 + 0.5 * Random() / (PRANDMAX - 1)) * (ZMax - ZMin);
			CLow = 0;
			CHigh = 2;
			/* Binary search */
			do {
				CMid = (CLow + CHigh) / 2;
				info.Size = 0;
				WindowSize(N->X - CMid * DX, N->X + CMid * DX,
					N->Y - CMid * DY, N->Y + CMid * DY,
					N->Z - CMid * DZ, N->Z + CMid * DZ,
					0, Dimension - 1, this, info);
				if (info.Size >= 2.0 / 3 * SubproblemSize && info.Size <= SubproblemSize)
					break;
				if (info.Size < 2.0 / 3 * SubproblemSize)
					CLow = CMid;
				else
					CHigh = CMid;
			} while (CHigh - CLow > std::numeric_limits<double>::min());
			MakeSubproblem(N->X - CMid * DX, N->X + CMid * DX,
				N->Y - CMid * DY, N->Y + CMid * DY,
				N->Z - CMid * DZ, N->Z + CMid * DZ,
				++Subproblems, 0, Dimension - 1, this,info);
			Remaining -= info.Size;
		}
		if (Remaining > 3) {
			Subproblems++;
			N = FirstNode;
			do
				N->Subproblem = Subproblems;
			while ((N = N->Suc) != FirstNode);
		}
		if (WeightType == GEO || WeightType == GEOM || WeightType == GEO_MEEUS
			|| WeightType == GEOM_MEEUS) {
			N = FirstNode;
			do {
				N->X = N->Xc;
				N->Y = N->Yc;
				N->Z = N->Zc;
			} while ((N = N->SubproblemSuc) != FirstNode);
			CoordType = TWOD_COORDS;
		}
		free(info.KDTree);
		for (CurrentSubproblem = 1;
			CurrentSubproblem <= Subproblems; CurrentSubproblem++) {
			OldGlobalBestCost = GlobalBestCost;
			SolveSubproblem(CurrentSubproblem, Subproblems, &GlobalBestCost);
			if (SubproblemsCompressed && GlobalBestCost == OldGlobalBestCost)
				SolveCompressedSubproblem(CurrentSubproblem, Subproblems,
					&GlobalBestCost);
		}
		printff("\nCost = " GainFormat, GlobalBestCost);
		if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
			printff(", Gap = %0.4f%%",
				100.0 * (GlobalBestCost - Optimum) / Optimum);
		printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
			GlobalBestCost < Optimum ? "<" : GlobalBestCost ==
			Optimum ? "=" : "");
		if (SubproblemBorders && Subproblems > 1)
			SolveSubproblemBorderProblems(Subproblems, &GlobalBestCost);
	}

	/*
	 * The WindowSize function computes the number of unused nodes that belong
	 * to a given box.
	 */

	static void WindowSize(double XMin, double XMax, double YMin, double YMax,
		double ZMin, double ZMax, int start, int end, LKHAlg *Alg, RoheSubInfo&info)
	{
		LKHAlg::Node *N;

		if (end - start + 1 <= Alg->SubproblemSize) {
			int i;
			for (i = start; i <= end; i++) {
				N = info.KDTree[i];
				if (N->Subproblem == 0 && N->X >= XMin && N->X <= XMax &&
					N->Y >= YMin && N->Y <= YMax &&
					N->Z >= ZMin && N->Z <= ZMax)
					if (++info.Size > Alg->SubproblemSize)
						return;
			}
		}
		else {
			int mid = (start + end) / 2;
			N = info.KDTree[mid];
			if (N->Subproblem == 0 && N->X >= XMin && N->X <= XMax &&
				N->Y >= YMin && N->Y <= YMax && N->Z >= ZMin && N->Z <= ZMax)
				info.Size++;
			if (info.Size <= Alg->SubproblemSize &&
				(N->Axis == 0 ? N->X >= XMin : N->Axis == 1 ? N->Y >= YMin :
					N->Z >= ZMin))
				WindowSize(XMin, XMax, YMin, YMax, ZMin, ZMax, start, mid - 1, Alg, info);
			if (info.Size <= Alg->SubproblemSize &&
				(N->Axis == 0 ? N->X <= XMax : N->Axis == 1 ? N->Y <= YMax :
					N->Z <= ZMax))
				WindowSize(XMin, XMax, YMin, YMax, ZMin, ZMax, mid + 1, end, Alg, info);
		}
	}

	/*
	 * The MakeSubproblem function associates each unused node belonging to
	 * a given box with a subproblem number. Each such node is removed
	 * from the list of unused nodes.
	 */

	static void MakeSubproblem(double XMin, double XMax, double YMin,
		double YMax, double ZMin, double ZMax,
		int Subproblem, int start, int end, LKHAlg *Alg, RoheSubInfo& info)
	{
		LKHAlg::Node *N, *Next;

		if (end - start + 1 <= Alg->SubproblemSize) {
			int i;
			for (i = start; i <= end; i++) {
				N = info.KDTree[i];
				if (N->Subproblem == 0 && N->X >= XMin && N->X <= XMax &&
					N->Y >= YMin && N->Y <= YMax &&
					N->Z >= ZMin && N->Z <= ZMax) {
					N->Subproblem = Subproblem;
					if (N == Alg->FirstNode) {
						Next = N->Suc;
						Alg->Follow(N, N);
						Alg->FirstNode = Next;
					}
					else
						Alg->Follow(N, N);
				}
			}
		}
		else {
			int mid = (start + end) / 2;
			N = info.KDTree[mid];
			if (N->Subproblem == 0 && N->X >= XMin && N->X <= XMax &&
				N->Y >= YMin && N->Y <= YMax && N->Z >= ZMin && N->Z <= ZMax) {
				N->Subproblem = Subproblem;
				if (N == Alg->FirstNode) {
					Next = N->Suc;
					Alg->Follow(N, N);
					Alg->FirstNode = Next;
				}
				else
					Alg->Follow(N, N);
			}
			if (N->Axis == 0 ? N->X >= XMin : N->Axis == 1 ? N->Y >=
				YMin : N->Z >= ZMin)
				MakeSubproblem(XMin, XMax, YMin, YMax, ZMin, ZMax, Subproblem,
					start, mid - 1, Alg, info);
			if (N->Axis == 0 ? N->X <= XMax : N->Axis == 1 ? N->Y <=
				YMax : N->Z <= ZMax)
				MakeSubproblem(XMin, XMax, YMin, YMax, ZMin, ZMax, Subproblem,
					mid + 1, end, Alg, info);
		}
	}
}