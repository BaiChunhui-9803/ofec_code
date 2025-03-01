#include "./INCLUDE/LKH.h"
#include "./INCLUDE/GeoConversion.h"

/*
 * The SolveKMeansSubproblems function attempts to improve a given tour
 * by means of a partitioning scheme based on K-means clustering.
 *
 * The overall region containing the nodes is subdivided into K clusters,
 * where K = ceil(Dimension/SubproblemSize). Each cluster together with
 * the given tour induces a subproblem consisting of all nodes in the
 * cluster and with edges fixed between nodes that are connected by tour
 * segments whose interior points do not belong to the cluster.
 *  
 * If an improvement is found, the new tour is written to TourFile.
 * The original tour is given by the SubproblemSuc references of the nodes.
 */
namespace LKH {
	static void KMeansClustering(int K, LKHAlg *Alg);

	void LKHAlg::SolveKMeansSubproblems()
	{
		Node *N;
		GainType GlobalBestCost, OldGlobalBestCost;
		double EntryTime = GetTime();
		int CurrentSubproblem, Subproblems;

		AllocateStructures();
		ReadPenalties();

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
			printff("*** K-means partitioning *** [Cost = " GainFormat "]\n",
				GlobalBestCost);
		}

		Subproblems = (int)ceil((double)Dimension / SubproblemSize);
		KMeansClustering(Subproblems, this);
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

#define M Beta

	/*
	 * The KMeansClustering function performs K-means clustering. Each node
	 * is given a unique cluster number in its Subproblem field.
	 *
	 * The algorithm is accellerated using ideas from
	 *
	 *     Dan Judd, Philip K. McKinley, Anil K. Jain:
	 *     Large-Scale Parallel Data Clustering.
	 *     IEEE Transactions on Pattern Analysis and
	 *     Machine Intelligence 20(8):871-876 (1998)
	 */

	static void KMeansClustering(int K, LKHAlg *Alg)
	{
		LKHAlg::Node *Center, **CenterRef, **Perm, *N, Old;
		int *Count, i, j, d, OldSubproblem;
		double *SumXc, *SumYc, *SumZc, Xc, Yc, Zc;
		int *Movement, *MMax, Max;
		int Moving = 0;
		LKHAlg::CostFunction OldDistance = Alg->Distance;
		//Alg->OldDistance = Alg->Distance;////

		assert(Center = (LKHAlg::Node *) calloc((K + 1), sizeof(LKHAlg::Node)));
		assert(CenterRef = (LKHAlg::Node **) calloc((K + 1), sizeof(LKHAlg::Node *)));
		assert(SumXc = (double *)calloc(K + 1, sizeof(double)));
		assert(SumYc = (double *)calloc(K + 1, sizeof(double)));
		assert(SumZc = (double *)calloc(K + 1, sizeof(double)));
		assert(Count = (int *)calloc(K + 1, sizeof(int)));
		assert(Movement = (int *)calloc(K + 1, sizeof(int)));
		assert(MMax = (int *)calloc(K + 1, sizeof(int)));
		assert(Perm = (LKHAlg::Node **) malloc(Alg->Dimension * sizeof(LKHAlg::Node *)));

		/* Pick random initial centers */
		for (i = 0; i < Alg->Dimension; i++)
			Perm[i] = &Alg->NodeSet[i + 1];
		for (j = 1; j < Alg->Dimension; j++) {
			i = Alg->Random() % j;
			N = Perm[j];
			Perm[j] = Perm[i];
			Perm[i] = N;
		}
		for (i = 1; i <= K; i++) {
			Center[i].X = Perm[i - 1]->X;
			Center[i].Y = Perm[i - 1]->Y;
			Center[i].Z = Perm[i - 1]->Z;
		}

		/* Assign each node to center 0 (a ghost cluster) */
		N = Alg->FirstNode;
		do {
			N->BestPi = N->Pi;
			N->Pi = 0;
			if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEO_MEEUS)
				GEO2XYZ(N->X, N->Y, &N->Xc, &N->Yc, &N->Zc);
			else if (Alg->WeightType == LKH::GEOM || Alg->WeightType == LKH::GEOM_MEEUS)
				GEOM2XYZ(N->X, N->Y, &N->Xc, &N->Yc, &N->Zc);
			else {
				N->Xc = N->X;
				N->Yc = N->Y;
				N->Zc = N->Z;
			}
			N->Cost = std::numeric_limits<int>::max() / 2;
			N->M = std::numeric_limits<int>::min();
			N->Subproblem = N->LastV = 0;
			SumXc[0] += N->Xc;
			SumYc[0] += N->Yc;
			SumZc[0] += N->Zc;
			Count[0]++;
		} while ((N = N->Suc) != Alg->FirstNode);
		Xc = Center[0].Xc = SumXc[0] / Count[0];
		Yc = Center[0].Yc = SumYc[0] / Count[0];
		Zc = Center[0].Zc = SumZc[0] / Count[0];
		if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEO_MEEUS)
			XYZ2GEO(Xc, Yc, Zc, &Center[0].X, &Center[0].Y);
		if (Alg->WeightType == LKH::GEOM || Alg->WeightType == LKH::GEOM_MEEUS)
			XYZ2GEOM(Xc, Yc, Zc, &Center[0].X, &Center[0].Y);
		else {
			Center[0].X = Xc;
			Center[0].Y = Yc;
			Center[0].Z = Zc;
		}
		if (Alg->Distance == &LKHAlg::Distance_TOR_2D)
			Alg->Distance = &LKHAlg::Distance_EUC_2D;
		else if (Alg->Distance == &LKHAlg::Distance_TOR_3D)
			Alg->Distance = &LKHAlg::Distance_EUC_3D;
		for (;;) {
			/* Assign each node to its closest center */
			Moving = 0;
			do {
				N->M -= Movement[N->Subproblem] + MMax[N->Subproblem];
				if (N->M < 0) {
					OldSubproblem = N->Subproblem;
					if (OldSubproblem == 0)
						N->Cost = N->NextCost = std::numeric_limits<int>::max();
					else {
						N->Cost =
							(Alg->*(Alg->Distance))(N, &Center[N->Subproblem]) * Alg->Precision;
						N->NextCost =
							(Alg->*(Alg->Distance))(N, &Center[N->LastV]) * Alg->Precision;
						if (N->Cost > N->NextCost) {
							i = N->LastV;
							N->LastV = N->Subproblem;
							N->Subproblem = i;
							d = N->NextCost;
							N->NextCost = N->Cost;
							N->Cost = d;
						}
					}
					for (i = 1; i <= K; i++) {
						if (i == N->Subproblem || i == N->LastV)
							continue;
						d = std::numeric_limits<int>::min();
						if ((!Alg->c || (Alg->*(Alg->c))(N, &Center[i]) <= N->Cost) &&
							(d =
							(Alg->*(Alg->Distance))(N, &Center[i]) * Alg->Precision) <= N->Cost) {
							N->NextCost = N->Cost;
							N->Cost = d;
							N->LastV = N->Subproblem;
							if (d < N->NextCost)
								N->Subproblem = i;
						}
						else if (d < N->NextCost &&
							(!Alg->c || (Alg->*(Alg->c))(N, &Center[i]) < N->NextCost) &&
							(d != std::numeric_limits<int>::min() || (d =
							(Alg->*(Alg->Distance))(N,
								&Center[i]) *
								Alg->Precision) <
								N->NextCost)) {
							N->NextCost = d;
							N->LastV = i;
						}
					}
					N->M = N->NextCost - N->Cost;
					if (N->Subproblem != OldSubproblem) {
						Moving++;
						SumXc[OldSubproblem] -= N->Xc;
						SumYc[OldSubproblem] -= N->Yc;
						SumZc[OldSubproblem] -= N->Zc;
						Count[OldSubproblem]--;
						SumXc[N->Subproblem] += N->Xc;
						SumYc[N->Subproblem] += N->Yc;
						SumZc[N->Subproblem] += N->Zc;
						Count[N->Subproblem]++;
					}
				}
			} while ((N = N->Suc) != Alg->FirstNode);
			if (!Moving)
				break;
			if (Alg->TraceLevel >= 2)
				Alg->printff("Moving %d %s\n", Moving,
					Moving > 1 ? "points" : "point");
			/* Move centers */
			Max = std::numeric_limits<int>::min();
			for (i = 1; i <= K; i++) {
				if (Count[i] > 0) {
					Old.X = Center[i].X;
					Old.Y = Center[i].Y;
					Old.Z = Center[i].Z;
					Xc = Center[i].Xc = SumXc[i] / Count[i];
					Yc = Center[i].Yc = SumYc[i] / Count[i];
					Zc = Center[i].Zc = SumZc[i] / Count[i];
					if (Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEO_MEEUS)
						XYZ2GEO(Xc, Yc, Zc, &Center[i].X, &Center[i].Y);
					else if (Alg->WeightType == LKH::GEOM || Alg->WeightType == LKH::GEOM_MEEUS)
						XYZ2GEOM(Xc, Yc, Zc, &Center[i].X, &Center[i].Y);
					else {
						Center[i].X = Xc;
						Center[i].Y = Yc;
						Center[i].Z = Zc;
					}
					Movement[i] = (Alg->*(Alg->Distance))(&Old, &Center[i]) * Alg->Precision;
					if (Movement[i] > Max)
						Max = Movement[i];
				}
				else
					Movement[i] = 0;
			}
			for (i = 1; i <= K; i++) {
				if (Movement[i] != Max)
					MMax[i] = Max;
				else {
					MMax[i] = std::numeric_limits<int>::min();
					for (j = 1; j <= K; j++)
						if (j != i && Movement[j] > MMax[i])
							MMax[i] = Movement[j];
				}
			}
		}
		Alg->Distance = Alg->OldDistance;
		N = Alg->FirstNode;
		do
			N->Pi = N->BestPi;
		while ((N = N->Suc) != Alg->FirstNode);
		for (i = 1; i <= K; i++)
			CenterRef[i] = &Center[i];
		Alg->AdjustClusters(K, CenterRef);
		delete(Center); Center = nullptr;
		delete(CenterRef); CenterRef = nullptr;
		delete(SumXc); SumXc = nullptr;
		delete(SumYc); SumYc = nullptr;
		delete(SumZc); SumZc = nullptr;
		delete(Count); Count = nullptr;
		delete(Perm); Perm = nullptr;
		delete(Movement); Movement = nullptr;
		delete(MMax); MMax = nullptr;
	}
}