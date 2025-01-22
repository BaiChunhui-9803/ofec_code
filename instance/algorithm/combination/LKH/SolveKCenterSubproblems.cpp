#include "./INCLUDE/LKH.h"

/*
 * The SolveKCenterSubproblems function attempts to improve a given tour
 * by means of a partitioning scheme based on approximate K-center clustering. 
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
	static void KCenterClustering(int K, LKHAlg *Alg);

	void LKHAlg::SolveKCenterSubproblems()
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
			printff("*** K-center partitioning *** [Cost = " GainFormat "]\n",
				GlobalBestCost);
		}

		Subproblems = (int)ceil((double)Dimension / SubproblemSize);
		KCenterClustering(Subproblems, this);
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
	 * The KCenterClustering function performs K-center clustering. Each node
	 * is given a unique cluster number in its Subproblem field.
	 */

	static void KCenterClustering(int K, LKHAlg *Alg)
	{
		LKHAlg::Node **Center, *N;
		int d, dMax, i;

		assert(Center = (LKHAlg::Node **) calloc((K + 1), sizeof(LKHAlg::Node *)));

		/* Pick first cluster arbitrarily */
		Center[1] = &Alg->NodeSet[Alg->Random() % Alg->Dimension + 1];
		/* Assign all cities to cluster1 */
		N = Alg->FirstNode;
		do {
			N->Subproblem = 1;
			N->Cost = (Alg->*(Alg->Distance))(N, Center[1]);
		} while ((N = N->Suc) != Alg->FirstNode);
		for (i = 2; i <= K; i++) {
			/* Take as the cluster center Center[i] a city furthest from
			 * Center[1..i-1]. */
			dMax = std::numeric_limits<int>::min();
			N = Alg->FirstNode;
			do {
				if ((d = N->Cost) > dMax) {
					Center[i] = N;
					dMax = d;
				}
			} while ((N = N->Suc) != Alg->FirstNode);
			do {
				if ((d = (Alg->*(Alg->Distance))(N, Center[i])) < N->Cost) {
					N->Cost = d;
					N->Subproblem = i;
				}
			} while ((N = N->Suc) != Alg->FirstNode);
		}
		Alg->AdjustClusters(K, Center);
		free(Center);
	}
}