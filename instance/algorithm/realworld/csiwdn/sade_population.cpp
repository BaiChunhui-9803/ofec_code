#include "sade_population.h"
#include "../../../problem/realworld/csiwdn/csiwdn.h"
#include "amp_cc/amp_cc.h"

namespace ofec {
	SaDEPopulation::SaDEPopulation(size_t no, Environment *env, size_t dim) : 
		Population(no,env, dim){
		auto pro = env->problem();
		m_F.resize(no, std::vector<Real>(CAST_CSIWDN(pro)->numberSource()));
		m_cnt_success.resize(CAST_CSIWDN(pro)->numberSource());
		m_cnt_fail.resize(CAST_CSIWDN(pro)->numberSource());
		m_CRsuc.resize(CAST_CSIWDN(pro)->numberSource());
		mvv_CR.resize(no, std::vector<std::vector<Real>>(CAST_CSIWDN(pro)->numberSource(), std::vector<Real>(m_num_strategy)));
		m_CRm.resize(CAST_CSIWDN(pro)->numberSource(), std::vector<Real>(m_num_strategy, 0.5));
		m_pro_strategy.resize(CAST_CSIWDN(pro)->numberSource(), std::vector<Real>(m_num_strategy, 1. / m_num_strategy));
		m_strategy_selection.resize(no, std::vector<int>(CAST_CSIWDN(pro)->numberSource()));
		m_probability.resize(2);
		m_ST_data_obj.resize(CAST_CSIWDN(pro)->maxStartTimeSize(), std::make_pair(0.0, 0));
		m_duration_data_obj.resize(CAST_CSIWDN(pro)->maxDurationSize(), std::make_pair(0.0, 0));
	
		m_probability[0].resize(CAST_CSIWDN(pro)->maxStartTimeSize());
		m_probability[1].resize(CAST_CSIWDN(pro)->maxDurationSize());

		for (size_t j = 0; j < CAST_CSIWDN(pro)->numberSource(); ++j) {
			for (auto i = 0; i < m_pro_strategy.size(); ++i) {
				if (i > 0) m_pro_strategy[j][i] += m_pro_strategy[j][i - 1];
			}
		}
	}

	void SaDEPopulation::updateF(Environment *env, Random *rnd, const std::pair<int, int>& source_index) {
		int z = source_index.second, q = source_index.first;
		for (auto i = 0; i < size(); i++) {
			for (int j = q; j <= z; ++j) {
				m_F[i][j] = rnd->uniform.nextNonStd(0.5, 0.09);
			}
		}
	}

	void SaDEPopulation::updateCR(Environment *env, Random *rnd, const std::pair<int, int>& source_index) {
		int z = source_index.second, q = source_index.first;
		//if (m_iteration >= m_LP) {
		for (int j = q; j <= z; ++j) {
			if (m_CRsuc[j].size() > m_LP) {
				for (int k = 0; k < m_num_strategy; ++k) {
					std::vector<std::vector<Real>> t;
					for (auto it = m_CRsuc[j].begin(); it != m_CRsuc[j].end(); ++it) {
						for (auto i : it->at(k)) t[j].push_back(i);
					}
					if (t.size() > 0) {
						nth_element(t[j].begin(), t[j].begin() + t[j].size() / 2, t[j].end()); //中位数放中间
						m_CRm[j][k] = t[j][(t[j].size()) / 2];  //选中位数
					}
				}
			}
		}
		for (auto i = 0; i < size(); i++) {
			for (int j = q; j <= z; ++j) {
				for (int k = 0; k < m_num_strategy; ++k) {
					do {
						mvv_CR[i][j][k] = rnd->uniform.nextNonStd(m_CRm[j][k], 0.01);
					} while (mvv_CR[i][j][k] < 0 || mvv_CR[i][j][k]>1);
				}
			}
		}

	}

	void SaDEPopulation::updateProStrategy(Environment *env, const std::pair<int, int>& source_index) {

		auto pro = env->problem();
		std::vector<std::vector < std::list<Real>>> curmem(CAST_CSIWDN(pro)->numberSource(), std::vector<std::list<Real>>(m_num_strategy));
		std::vector<std::vector<int>> curSuc(CAST_CSIWDN(pro)->numberSource(), std::vector <int>(m_num_strategy)), curFail(CAST_CSIWDN(pro)->numberSource(), std::vector <int>(m_num_strategy));
		int z = source_index.second, q = source_index.first;
		for (auto i = 0; i < size(); i++) {
			for (int j = q; j <= z; ++j) {
				if (m_individuals[i]->isImproved()) {
					curmem[j][m_strategy_selection[i][j]].push_back(mvv_CR[i][j][m_strategy_selection[i][j]]);
					curSuc[j][m_strategy_selection[i][j]]++;
				}
				else {
					curFail[j][m_strategy_selection[i][j]]++;
				}
			}
		}

		for (int j = q; j <= z; ++j) {
			m_cnt_success[j].push_back(move(curSuc[j]));
			m_CRsuc[j].push_back(move(curmem[j]));
			m_cnt_fail[j].push_back(move(curFail[j]));

			if (m_CRsuc[j].size() > m_LP) {
				//if (m_iteration >= m_LP) {
				m_cnt_success[j].pop_front();
				m_CRsuc[j].pop_front();
				m_cnt_fail[j].pop_front();

				//update probability for all stategies
				for (int k = 0; k < m_num_strategy; ++k) {
					m_pro_strategy[j][k] = 0;
					std::vector<int> fail(CAST_CSIWDN(pro)->numberSource(), 0);
					for (auto &f : m_cnt_success[j]) m_pro_strategy[j][k] += f[k];
					for (auto &f : m_cnt_fail[j]) fail[j] += f[k];

					m_pro_strategy[j][k] = m_pro_strategy[j][k] / (m_pro_strategy[j][k] + fail[j]) + m_epsilon;
					if (k > 0) m_pro_strategy[j][k] += m_pro_strategy[j][k - 1];
				}
			}
		}
	}

	void SaDEPopulation::mutate(Environment *env, Random *rnd, const std::pair<int, int>& source_index) {
		int z = source_index.second, q = source_index.first;
		for (size_t i = 0; i < size(); ++i) {
			for (int j = q; j <= z; ++j) {
				Real p = rnd->uniform.next() * m_pro_strategy[j][m_num_strategy - 1];
				m_strategy_selection[i][j] = lower_bound(m_pro_strategy[j].begin(), m_pro_strategy[j].end(), p) - m_pro_strategy[j].begin();
				setMutationStrategy(DEMutationSratgy(m_strategy_selection[i][j]));
				updateBest(env);
				mutate_(env, rnd, i,j, m_F[i][j], m_probability);
			}
		}
	}

	void SaDEPopulation::setMutationStrategy(DEMutationSratgy rS) {
		m_mutation_strategy = rS;
	}

	void SaDEPopulation::recombine(Environment *env, Random *rnd, const std::pair<int, int>& source_index) {
		int z = source_index.second, q = source_index.first;
		for (size_t i = 0; i < size(); ++i) {
			for (size_t j = q; j <= z; ++j) {
				m_individuals[i]->recombine(env, rnd, j, mvv_CR[i][j][m_strategy_selection[i][j]]);
			}
		}
	}
	
	void SaDEPopulation::updateProbability(Environment *env) {

		auto pro = env->problem();
		/// update probability of starting time
		int z = 0;
		for (auto &i : m_individuals) {
			while ((z + 1) < CAST_CSIWDN(pro)->numberSource() && CAST_CSIWDN(pro)->phase() >= ((i->variable().startTime(z + 1) / CAST_CSIWDN(pro)->intervalTimeStep()))) {
				z++;
			}
			int idx = (i->variable().startTime(z) - CAST_CSIWDN(pro)->minStartTime()) / CAST_CSIWDN(pro)->patternStep();
			m_ST_data_obj[idx].first += i->objective()[0];
			++(m_ST_data_obj[idx].second);
		}
		std::vector<Real> mean_ST(m_ST_data_obj.size(), 0.0);
		size_t count_ST = 0;
		for (size_t i = 0; i < CAST_CSIWDN(pro)->maxStartTimeSize(); ++i) {
			if (m_ST_data_obj[i].second != 0) {
				mean_ST[i] = m_ST_data_obj[i].first / m_ST_data_obj[i].second;
				++count_ST;
			}
		}

		Real max = 0;
		for (auto &i : mean_ST)
			if (max < i) max = i;
		std::vector<Real> inverse;
		int size_ST = mean_ST.size();
		for (auto &i : mean_ST) {
			if (i == 0 || i == max)
				inverse.push_back(max / size_ST);
			else
				inverse.push_back(max - i);
		}
		Real sum = 0;
		for (auto &i : inverse)
			sum += i;
		m_probability[0].clear();
		for (auto &i : inverse)
			m_probability[0].push_back(i / sum);

		/// update probability of duration
		z = 0;
		for (auto &i : m_individuals) {
			while ((z + 1) < CAST_CSIWDN(pro)->numberSource() && CAST_CSIWDN(pro)->phase() >= ((i->variable().startTime(z + 1) / CAST_CSIWDN(pro)->intervalTimeStep()))) {
				z++;
			}
			int idx = (i->variable().duration(z) - CAST_CSIWDN(pro)->minDuration()) / CAST_CSIWDN(pro)->patternStep();
			m_duration_data_obj[idx].first += i->objective()[0];
			++(m_duration_data_obj[idx].second);
		}
		std::vector<Real> mean_duration(m_duration_data_obj.size(), 0.0);
		size_t count_duration = 0;
		for (size_t i = 0; i < CAST_CSIWDN(pro)->maxDurationSize(); ++i) {
			if (m_duration_data_obj[i].second != 0) {
				mean_duration[i] = m_duration_data_obj[i].first / m_duration_data_obj[i].second;
				++count_duration;
			}
		}

	    max = 0;
		for (auto &i : mean_duration)
			if (max < i) max = i;
		inverse.clear();
		int size_duration = mean_duration.size();
		for (auto &i : mean_duration) {
			if (i == 0 || i == max)
				inverse.push_back(max / size_duration);
			else
				inverse.push_back(max - i);
		}
		sum = 0;
		for (auto &i : inverse)
			sum += i;
		m_probability[1].clear();
		for (auto &i : inverse)
			m_probability[1].push_back(i / sum);
	}

	void SaDEPopulation::mutate_(Environment *env, Random *rnd, const int idx, const int jdx, Real F, const std::vector<std::vector<Real>> &prob) {
		std::vector<int> ridx;
		switch (m_mutation_strategy) {
		case rand_1:
			select(rnd,idx, 3, ridx);
			this->m_individuals[idx]->mutateSecondPart(env, F, jdx, prob, m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get());
			break;
		case best_1:
			select(rnd, idx, 2, ridx);
			this->m_individuals[idx]->mutateSecondPart(env, F, jdx, prob, m_best.front().get(), m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get());
			break;
		case target_to_best_1:
			select(rnd, idx, 2, ridx);
			this->m_individuals[idx]->mutateSecondPart(env, F, jdx, prob, m_individuals[idx].get(), m_best.front().get(), m_individuals[idx].get(), m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get());
			break;
		case best_2:
			select(rnd, idx, 4, ridx);
			this->m_individuals[idx]->mutateSecondPart(env, F, jdx, prob, m_best.front().get(), m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get(), m_individuals[ridx[3]].get());
			break;
		case rand_2:
			select(rnd, idx, 5, ridx);
			this->m_individuals[idx]->mutateSecondPart(env, F, jdx, prob, m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get(), m_individuals[ridx[3]].get(), m_individuals[ridx[4]].get());
			break;
		case rand_to_best_1:
			select(rnd, idx, 3, ridx);
			this->m_individuals[idx]->mutateSecondPart(env, F, jdx, prob, m_individuals[ridx[0]].get(), m_best.front().get(), m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get());
			break;
		case target_to_rand_1:
			select(rnd,idx, 3, ridx);
			this->m_individuals[idx]->mutateSecondPart(env,F, jdx, prob, m_individuals[idx].get(), m_individuals[ridx[0]].get(), m_individuals[idx].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get());
			break;
		}
	}

	void SaDEPopulation::select(Random *rnd,int base, int number, std::vector<int>& result) {
		std::vector<int> candidate;
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			if (base != i) candidate.push_back(i);
		}
		if (result.size() != number)	result.resize(number);
		for (size_t i = 0; i < number; ++i) {
			size_t idx = rnd->uniform.nextNonStd<size_t>(0, candidate.size() - i);
			result[i] = candidate[idx];
			if (idx != candidate.size() - (i + 1)) candidate[idx] = candidate[candidate.size() - (i + 1)];
		}
	}

	void SaDEPopulation::fillSolution(VarCSIWDN& indi,Environment *env, const std::pair<int, int>& source_index) {
		for (auto &i : m_individuals)
			i->coverFirstPart(indi, env, source_index);
	}

	void SaDEPopulation::select(Environment *env, bool is_stable, const std::pair<int, int>& source_index) {
		for (auto &i : m_individuals)
			i->select(env, is_stable, source_index);

	}

	bool SaDEPopulation::isFeasiblePopulation(Environment *env,const Real tar) {
		auto pro = env->problem();
		size_t phase = CAST_CSIWDN(pro)->phase();
		size_t interval = CAST_CSIWDN(pro)->interval();
		size_t num = phase*interval;
		Real temp = 0;

		for (size_t i = 0; i < num; ++i) {
			for (int j = 0; j < CAST_CSIWDN(pro)->numSensor(); ++j) {
				temp += pow(CAST_CSIWDN(pro)->observationConcentration()[j][i], 2);
			}
		}

		Real benchmark = tar * sqrt(temp / (CAST_CSIWDN(pro)->numSensor()*num));

		size_t count_feasible = 0, count_infeasible = 0;
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			if (m_individuals[i]->objective()[0] <= benchmark) ++count_feasible;
			else ++count_infeasible;
		}

		return count_feasible >= count_infeasible;
	}
}