#include "./INCLUDE/LKH.h"
namespace LKH {


	void LKHAlg::InitializeStatistics()
	{
		m_statistic_info.TrialSum = m_statistic_info.Successes = 0;
		m_statistic_info.CostSum = 0;
		m_statistic_info.TimeSum = 0.0;
		m_statistic_info.TrialsMin = std::numeric_limits<int>::max();
		m_statistic_info.TrialsMax = 0;
		m_statistic_info.TimeMin = std::numeric_limits<double>::max();
		m_statistic_info.TimeMax = 0;
		m_statistic_info.CostMin = std::numeric_limits<GainType>::max();
		m_statistic_info.CostMax = std::numeric_limits<GainType>::min();
		m_statistic_info.PenaltySum = 0;
		m_statistic_info.PenaltyMin = std::numeric_limits<GainType>::max();
		m_statistic_info.PenaltyMax = std::numeric_limits<GainType>::min();
	}

	void LKHAlg::UpdateStatistics(GainType Cost, double Time)
	{
		if (Trial < m_statistic_info.TrialsMin)
			m_statistic_info.TrialsMin = Trial;
		if (Trial > m_statistic_info.TrialsMax)
			m_statistic_info.TrialsMax = Trial;
		m_statistic_info.TrialSum += Trial;
		if (Cost < m_statistic_info.CostMin)
			m_statistic_info.CostMin = Cost;
		if (Cost > m_statistic_info.CostMax)
			m_statistic_info.CostMax = Cost;
		m_statistic_info.CostSum += Cost;
		if (ProblemType != CCVRP && ProblemType != TRP &&
			ProblemType != MLP &&
			MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE) {
			if (CurrentPenalty == 0 && Cost <= Optimum)
				m_statistic_info.Successes++;
		}
		else if (CurrentPenalty <= Optimum)
			m_statistic_info.Successes++;
		if (CurrentPenalty < m_statistic_info.PenaltyMin)
			m_statistic_info.PenaltyMin = CurrentPenalty;
		if (CurrentPenalty > m_statistic_info.PenaltyMax)
			m_statistic_info.PenaltyMax = CurrentPenalty;
		m_statistic_info.PenaltySum += CurrentPenalty;
		if (Time < m_statistic_info.TimeMin)
			m_statistic_info.TimeMin = Time;
		if (Time > m_statistic_info.TimeMax)
			m_statistic_info.TimeMax = Time;
		m_statistic_info.TimeSum += Time;
	}

	void LKHAlg::PrintStatistics()
	{
		int _Runs = Runs, _TrialsMin = m_statistic_info.TrialsMin;
		double _TimeMin = m_statistic_info.TimeMin;
		GainType _Optimum = Optimum;

		printff("m_statistic_info.Successes/Runs = %d/%d \n", m_statistic_info.Successes, Runs);
		if (_Runs == 0)
			_Runs = 1;
		if (_TrialsMin > m_statistic_info.TrialsMax)
			_TrialsMin = 0;
		if (_TimeMin > m_statistic_info.TimeMax)
			_TimeMin = 0;
		if (ProblemType != CCVRP && ProblemType != TRP &&
			ProblemType != MLP &&
			MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE &&
			m_statistic_info.CostMin <= m_statistic_info.CostMax && m_statistic_info.CostMin != std::numeric_limits<GainType>::max()) {
			printff("Cost.min = " GainFormat ", Cost.avg = %0.2f, "
				"Cost.max = " GainFormat "\n",
				m_statistic_info.CostMin, (double)m_statistic_info.CostSum / _Runs, m_statistic_info.CostMax);
			if (_Optimum == std::numeric_limits<GainType>::min())
				_Optimum = BestCost;
			if (_Optimum != 0)
				printff("Gap.min = %0.4f%%, Gap.avg = %0.4f%%, "
					"Gap.max = %0.4f%%\n",
					100.0 * (m_statistic_info.CostMin - _Optimum) / _Optimum,
					100.0 * ((double)m_statistic_info.CostSum / _Runs -
						_Optimum) / _Optimum,
					100.0 * (m_statistic_info.CostMax - _Optimum) / _Optimum);
			if (Penalty && m_statistic_info.PenaltyMin != std::numeric_limits<GainType>::max())
				printff("Penalty.min = " GainFormat ", Penalty.avg = %0.2f, "
					"Penalty.max = " GainFormat "\n",
					m_statistic_info.PenaltyMin, (double)m_statistic_info.PenaltySum / _Runs, m_statistic_info.PenaltyMax);
		}
		else if (Penalty && m_statistic_info.PenaltyMin != std::numeric_limits<GainType>::max()) {
			printff("Penalty.min = " GainFormat ", Penalty.avg = %0.2f, "
				"Penalty.max = " GainFormat "\n",
				m_statistic_info.PenaltyMin, (double)m_statistic_info.PenaltySum / _Runs, m_statistic_info.PenaltyMax);
			if (_Optimum == std::numeric_limits<GainType>::min())
				_Optimum = BestPenalty;
			if (_Optimum != 0)
				printff("Gap.min = %0.4f%%, Gap.avg = %0.4f%%, "
					"Gap.max = %0.4f%%\n",
					100.0 * (m_statistic_info.PenaltyMin - _Optimum) / _Optimum,
					100.0 * ((double)m_statistic_info.PenaltySum / _Runs -
						_Optimum) / _Optimum,
					100.0 * (m_statistic_info.PenaltyMax - _Optimum) / _Optimum);
		}
		printff("Trials.min = %d, Trials.avg = %0.1f, Trials.max = %d\n",
			_TrialsMin, 1.0 * m_statistic_info.TrialSum / _Runs, m_statistic_info.TrialsMax);
		printff
		("Time.min = %0.2f sec., Time.avg = %0.2f sec., Time.max = %0.2f sec.\n",
			fabs(_TimeMin), fabs(m_statistic_info.TimeSum) / _Runs, fabs(m_statistic_info.TimeMax));
	}
}