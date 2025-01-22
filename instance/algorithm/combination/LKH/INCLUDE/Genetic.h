#ifndef _GENETIC_H
#define _GENETIC_H
#include "LKH.h"
/*
 * This header specifies the interface for the genetic algorithm part of LKH.
 */
namespace LKH {

	//extern thread_local unique_ptr<int> MaxPopulationSize; /* The maximum size of the population */
	//extern thread_local unique_ptr<int> PopulationSize;    /* The current size of the population */
	//extern thread_local unique_ptr<CrossoverFunction> Crossover;
	//extern thread_local unique_ptr<int*> Population;      /* Array of individuals (solution tours) */
	//extern thread_local unique_ptr<GainType> Fitness;     /* The fitness (tour cost) of each individual */
	//extern thread_local unique_ptr<GainType*> PenaltyFitness;     /* The fitness (tour cost) of each individual */

#ifdef __cplusplus
	extern "C" {
#endif
#define SmallerFitness(Penalty, Cost, i, info)\
    (((Penalty) < info.PenaltyFitness[i]) ||\
     ((Penalty) == info.PenaltyFitness[i] && (Cost) < info.Fitness[i]))

#define LargerFitness(Penalty, Cost, i, info)\
    (((Penalty) > info.PenaltyFitness[i]) ||\
     ((Penalty) == info.PenaltyFitness[i] && (Cost) > info.Fitness[i]))



		void AddToPopulation(GainType Penalty, GainType Cost, LKH::LKHAlg *Alg);
		void ApplyCrossover(int i, int j, LKH::LKHAlg *Alg);
		void FreePopulation(LKH::LKHAlg* Alg);
		int HasFitness(GainType Penalty, GainType Cost, LKH::LKHAlg* Alg);
		int LinearSelection(int Size, double Bias, LKH::LKHAlg *Alg);
		GainType MergeTourWithIndividual(int i, LKH::LKHAlg *Alg);
		void PrintPopulation(LKH::LKHAlg *Alg);
		void ReplaceIndividualWithTour(int i, GainType Penalty, GainType Cost, LKH::LKHAlg *Alg);
		int ReplacementIndividual(GainType Penalty, GainType Cost, LKH::LKHAlg *Alg);

		void ERXT(LKH::LKHAlg *Alg);
#ifdef __cplusplus
	}
#endif
}
#endif
