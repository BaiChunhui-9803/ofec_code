#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Genetic.h"

/*
 * The AddToPopulation function adds the current tour as an individual to 
 * the population. The fitness of the individual is set equal to the cost
 * of the tour. The population is kept sorted in increasing fitness order.
 */

//thread_local unique_ptr<int> MaxPopulationSize; /* The maximum size of the population */
//thread_local unique_ptr<int> PopulationSize;    /* The current size of the population */
//thread_local unique_ptr<CrossoverFunction> Crossover;
//thread_local unique_ptr<int*> Population;      /* Array of individuals (solution tours) */
//thread_local unique_ptr<GainType> Fitness;     /* The fitness (tour cost) of each individual */
//thread_local unique_ptr<GainType*> PenaltyFitness;     /* The fitness (tour cost) of each individual */
namespace LKH {
	thread_local int MaxPopulationSize; /* The maximum size of the population */
	thread_local int PopulationSize;    /* The current size of the population */

	thread_local CrossoverFunction Crossover;
	thread_local int **Population = 0;      /* Array of individuals (solution tours) */
	thread_local GainType *PenaltyFitness = 0;  /* The f itnesslty  (tour penalty) of each
	i                             individual */
	thread_local GainType *Fitness = 0;     /* The fitness (tour cost) of each individual */

	void AddToPopulation(GainType Penalty, GainType Cost, LKHAlg *Alg)
	{
		int i, *P;
		LKHAlg::Node *N;
		/*if (!MaxPopulationSize.get()) {
			MaxPopulationSize.reset(new int(0));
		}*/
		if (!Population) {
			assert(Population =
				(int **)malloc(MaxPopulationSize * sizeof(int *)));
			for (i = 0; i < MaxPopulationSize; i++)
				assert(Population[i] =
				(int *)malloc((1 + Alg->Dimension) * sizeof(int)));
			assert(PenaltyFitness =
				(GainType *)malloc(MaxPopulationSize * sizeof(GainType)));
			assert(Fitness =
				(GainType *)malloc(MaxPopulationSize * sizeof(GainType)));
		}
		for (i = PopulationSize;
			i >= 1 && SmallerFitness(Penalty, Cost, i - 1); i--) {
			PenaltyFitness[i] = PenaltyFitness[i - 1];
			Fitness[i] = Fitness[i - 1];
			P = Population[i];
			Population[i] = Population[i - 1];
			Population[i - 1] = P;
		}
		PenaltyFitness[i] = Penalty;
		Fitness[i] = Cost;
		P = Population[i];
		N = Alg->FirstNode;
		i = 1;
		do
			P[i++] = N->Id;
		while ((N = N->Suc) != Alg->FirstNode);
		P[0] = P[Alg->Dimension];
		PopulationSize++;
	}

	/*
	 * The ApplyCrossover function applies a specified crossover operator to two
	 * individuals.
	 */

	void ApplyCrossover(int i, int j, LKHAlg *Alg)
	{
		int *Pi, *Pj, k;

		Pi = Population[i];
		Pj = Population[j];
		for (k = 1; k <= Alg->Dimension; k++) {
			Alg->NodeSet[Pi[k - 1]].Suc = &Alg->NodeSet[Pi[k]];
			Alg->NodeSet[Pj[k - 1]].Next = &Alg->NodeSet[Pj[k]];
		}
		if (Alg->TraceLevel >= 1)
			Alg->printff("Crossover(%d,%d)\n", i + 1, j + 1);
		/* Apply the crossover operator */
		Crossover(Alg);
		if (Alg->ProblemType == LKH::SOP)
			Alg->SOP_RepairTour();
	}

#define Free(s) { free(s); s = 0; }

	/*
	 * The FreePopulation function frees the memory space allocated to the
	 * population.
	 */

	void FreePopulation()
	{
		if (Population) {
			int i;
			for (i = 0; i < MaxPopulationSize; i++)
				Free(Population[i]);
			Free(Population);
			Free(PenaltyFitness);
			Free(Fitness);
		}
		PopulationSize = 0;
	}

	/*
	 * The HasFitness function returns 1 if the population contains an
	 * individual with fitness equal to a given tour cost; otherwise 0.
	 *
	 * Since the population is sorted in fitness order the test may be
	 * made by binary search.
	 */

	int HasFitness(GainType Penalty, GainType Cost)
	{
		int Low = 0, High = PopulationSize - 1;
		while (Low < High) {
			int Mid = (Low + High) / 2;
			if (LargerFitness(Penalty, Cost, Mid))
				Low = Mid + 1;
			else
				High = Mid;
		}
		return High >= 0 && PenaltyFitness[High] == Penalty &&
			Fitness[High] == Cost;
	}

	/*
	 * Random01 is an auxiliary function for computing a random double number
	 * in the range [0;1).
	 */

	static double Random01(LKHAlg *Alg)
	{
		return ((double)Alg->Random()) / std::numeric_limits<int>::max();
	}

	/*
	 * The LinearSelection function is used to select an individual with
	 * random linear bias towards the best members of the population.
	 * The parameter Bias is a number between 1.0 and 2.0.
	 *
	 * See
	 *     Darrell Whitley,
	 *     The GENITOR algorithm and selection pressure:
	 *     Why rank-based allocation of reproductive trials is best.
	 *     Proceedings of the Third International Conference on Genetic Algorithms,
	 *     1989.
	 */

	int LinearSelection(int Size, double Bias, LKHAlg *Alg)
	{
		return (int)(Size *
			(Bias -
				sqrt((Bias * Bias - 4 * (Bias - 1) * Random01(Alg)))) /
			2 / (Bias - 1));
	}

	/*
	 * The MergeTourWithIndividual function attempts to find a short tour by
	 * merging the current tour with a specified inddividual of the population.
	 * The merging algorithm is the iterative partial transcription algrithm
	 * described by Mobius, Freisleben, Merz and Schreiber.
	 */

	GainType MergeTourWithIndividual(int i, LKHAlg *Alg)
	{
		int *Pi, k;

		assert(i >= 0 && i < PopulationSize);
		Pi = Population[i];
		for (k = 1; k <= Alg->Dimension; k++)
			Alg->NodeSet[Pi[k - 1]].Next = &Alg->NodeSet[Pi[k]];
		return (Alg->*(Alg->MergeWithTour))();
	}

	/*
	 * The PrintPopulation function prints the cost and gap to optimum for
	 * each individual of the population.
	 */

	void PrintPopulation(LKHAlg *Alg)
	{
		int i;
		Alg->printff("Population:\n");
		for (i = 0; i < PopulationSize; i++) {
			Alg->printff("%3d: ", i + 1);
			if (Alg->Penalty)
				Alg->printff(GainFormat "_" GainFormat,
					PenaltyFitness[i], Fitness[i]);
			else
				Alg->printff(GainFormat, Fitness[i]);
			if (Alg->Optimum != std::numeric_limits<GainType>::min() && Alg->Optimum != 0) {
				if (Alg->ProblemType != LKH::CCVRP && Alg->ProblemType != LKH::TRP &&
					Alg->ProblemType != LKH::MLP &&
					Alg->MTSPObjective != LKH::MINMAX &&
					Alg->MTSPObjective != LKH::MINMAX_SIZE)
					Alg->printff(", Gap = %0.4f%%",
						100.0 * (Fitness[i] - Alg->Optimum) / Alg->Optimum);
				else
					Alg->printff(", Gap = %0.4f%%",
						100.0 * (PenaltyFitness[i] - Alg->Optimum) / Alg->Optimum);
			}
			Alg->printff("\n");
		}
	}

	/*
	 * The ReplaceIndividualWithTour function replaces a given individual in
	 * the population by an indidual that represents the current tour.
	 * The population is kept sorted in increasing fitness order.
	 */

	void ReplaceIndividualWithTour(int i, GainType Penalty, GainType Cost, LKHAlg *Alg)
	{
		int j, *P;
		LKHAlg::Node *N;

		assert(i >= 0 && i < PopulationSize);
		PenaltyFitness[i] = Penalty;
		Fitness[i] = Cost;
		P = Population[i];
		N = Alg->FirstNode;
		for (j = 1; j <= Alg->Dimension; j++) {
			P[j] = N->Id;
			N = N->Suc;
		}
		P[0] = P[Alg->Dimension];
		while (i >= 1 && SmallerFitness(Penalty, Cost, i - 1)) {
			PenaltyFitness[i] = PenaltyFitness[i - 1];
			Fitness[i] = Fitness[i - 1];
			Population[i] = Population[i - 1];
			i--;
		}
		PenaltyFitness[i] = Cost;
		Fitness[i] = Cost;
		Population[i] = P;
		while (i < PopulationSize - 1 && LargerFitness(Penalty, Cost, i + 1)) {
			PenaltyFitness[i] = PenaltyFitness[i + 1];
			Population[i] = Population[i + 1];
			i++;
		}
		PenaltyFitness[i] = Penalty;
		Fitness[i] = Cost;
		Population[i] = P;
	}

	/*
	 * The DistanceToIndividual returns the number of different edges between
	 * the tour (given by OldSuc) and individual i.
	 */

	static int DistanceToIndividual(int i, LKHAlg *Alg)
	{
		int Count = 0, j, *P = Population[i];
		LKHAlg::Node *N;

		for (j = 0; j < Alg->Dimension; j++) {
			N = &Alg->NodeSet[P[j]];
			(N->Next = &Alg->NodeSet[P[j + 1]])->Prev = N;
		}
		N = Alg->FirstNode;
		do
			if (N->OldSuc != N->Next && N->OldSuc != N->Prev)
				Count++;
		while ((N = N->OldSuc) != Alg->FirstNode);
		return Count;
	}

	/*
	 * The ReplacementIndividual function returns the individual to be
	 * replaced with the current tour. The function implements the
	 * replacement strategy (CD/RW) proposed in
	 *
	 *      M. Lozano, F. Herrera, and J. R. Cano,
	 *      Replacement strategies to preserve useful diversity in
	 *      steady-state genetic algorithms.
	 *      Information Sciences 178 (2008) 4421–4433.
	 */

	int ReplacementIndividual(GainType Penalty, GainType Cost, LKHAlg *Alg)
	{
		int i, j, d, *P;
		int MinDist = std::numeric_limits<int>::max(), CMin = PopulationSize - 1;
		LKHAlg::Node *N = Alg->FirstNode;
		while ((N = N->OldSuc = N->Suc) != Alg->FirstNode);
		for (i = PopulationSize - 1;
			i >= 0 && SmallerFitness(Penalty, Cost, i); i--) {
			if ((d = DistanceToIndividual(i, Alg)) < MinDist) {
				CMin = i;
				MinDist = d;
			}
		}
		if (CMin == PopulationSize - 1)
			return CMin;
		P = Population[CMin];
		for (j = 0; j < Alg->Dimension; j++)
			Alg->NodeSet[P[j]].OldSuc = &Alg->NodeSet[P[j + 1]];
		for (i = 0; i < PopulationSize; i++)
			if (i != CMin && (d = DistanceToIndividual(i, Alg)) <= MinDist)
				return PopulationSize - 1;
		return CMin;
	}
}