#include "./INCLUDE/LKH.h"
namespace LKH {
	static thread_local int TrialsMin, TrialsMax, TrialSum, Successes;
	static thread_local GainType CostMin, CostMax, CostSum;
	static thread_local GainType PenaltyMin, PenaltyMax, PenaltySum;
	static thread_local double TimeMin, TimeMax, TimeSum;

	void LKHAlg::InitializeStatistics()
	{
		TrialSum = Successes = 0;
		CostSum = 0;
		TimeSum = 0.0;
		TrialsMin = std::numeric_limits<int>::max();
		TrialsMax = 0;
		TimeMin = std::numeric_limits<double>::max();
		TimeMax = 0;
		CostMin = std::numeric_limits<GainType>::max();
		CostMax = std::numeric_limits<GainType>::min();
		PenaltySum = 0;
		PenaltyMin = std::numeric_limits<GainType>::max();
		PenaltyMax = std::numeric_limits<GainType>::min();
	}

	void LKHAlg::UpdateStatistics(GainType Cost, double Time)
	{
		if (Trial < TrialsMin)
			TrialsMin = Trial;
		if (Trial > TrialsMax)
			TrialsMax = Trial;
		TrialSum += Trial;
		if (Cost < CostMin)
			CostMin = Cost;
		if (Cost > CostMax)
			CostMax = Cost;
		CostSum += Cost;
		if (ProblemType != CCVRP && ProblemType != TRP &&
			ProblemType != MLP &&
			MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE) {
			if (CurrentPenalty == 0 && Cost <= Optimum)
				Successes++;
		}
		else if (CurrentPenalty <= Optimum)
			Successes++;
		if (CurrentPenalty < PenaltyMin)
			PenaltyMin = CurrentPenalty;
		if (CurrentPenalty > PenaltyMax)
			PenaltyMax = CurrentPenalty;
		PenaltySum += CurrentPenalty;
		if (Time < TimeMin)
			TimeMin = Time;
		if (Time > TimeMax)
			TimeMax = Time;
		TimeSum += Time;
	}

	void LKHAlg::PrintStatistics()
	{
		int _Runs = Runs, _TrialsMin = TrialsMin;
		double _TimeMin = TimeMin;
		GainType _Optimum = Optimum;

		printff("Successes/Runs = %d/%d \n", Successes, Runs);
		if (_Runs == 0)
			_Runs = 1;
		if (_TrialsMin > TrialsMax)
			_TrialsMin = 0;
		if (_TimeMin > TimeMax)
			_TimeMin = 0;
		if (ProblemType != CCVRP && ProblemType != TRP &&
			ProblemType != MLP &&
			MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE &&
			CostMin <= CostMax && CostMin != std::numeric_limits<GainType>::max()) {
			printff("Cost.min = " GainFormat ", Cost.avg = %0.2f, "
				"Cost.max = " GainFormat "\n",
				CostMin, (double)CostSum / _Runs, CostMax);
			if (_Optimum == std::numeric_limits<GainType>::min())
				_Optimum = BestCost;
			if (_Optimum != 0)
				printff("Gap.min = %0.4f%%, Gap.avg = %0.4f%%, "
					"Gap.max = %0.4f%%\n",
					100.0 * (CostMin - _Optimum) / _Optimum,
					100.0 * ((double)CostSum / _Runs -
						_Optimum) / _Optimum,
					100.0 * (CostMax - _Optimum) / _Optimum);
			if (Penalty && PenaltyMin != std::numeric_limits<GainType>::max())
				printff("Penalty.min = " GainFormat ", Penalty.avg = %0.2f, "
					"Penalty.max = " GainFormat "\n",
					PenaltyMin, (double)PenaltySum / _Runs, PenaltyMax);
		}
		else if (Penalty && PenaltyMin != std::numeric_limits<GainType>::max()) {
			printff("Penalty.min = " GainFormat ", Penalty.avg = %0.2f, "
				"Penalty.max = " GainFormat "\n",
				PenaltyMin, (double)PenaltySum / _Runs, PenaltyMax);
			if (_Optimum == std::numeric_limits<GainType>::min())
				_Optimum = BestPenalty;
			if (_Optimum != 0)
				printff("Gap.min = %0.4f%%, Gap.avg = %0.4f%%, "
					"Gap.max = %0.4f%%\n",
					100.0 * (PenaltyMin - _Optimum) / _Optimum,
					100.0 * ((double)PenaltySum / _Runs -
						_Optimum) / _Optimum,
					100.0 * (PenaltyMax - _Optimum) / _Optimum);
		}
		printff("Trials.min = %d, Trials.avg = %0.1f, Trials.max = %d\n",
			_TrialsMin, 1.0 * TrialSum / _Runs, TrialsMax);
		printff
		("Time.min = %0.2f sec., Time.avg = %0.2f sec., Time.max = %0.2f sec.\n",
			fabs(_TimeMin), fabs(TimeSum) / _Runs, fabs(TimeMax));
	}
}