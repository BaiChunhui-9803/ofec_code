#include "gl_adaptor_cont.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../template/framework/gl/gl_pop.h"
#include "../../../../../../utility/functional.h"

namespace ofec {
	AdaptorContGL::AdaptorContGL(Real alpha, size_t num_dim, size_t size_pop) :
		AdaptorGL<Solution<>>(alpha),
		m_proba(num_dim, std::vector<Infor>(m_num)),
		m_acc_proba(num_dim, std::vector<Real>(m_num)),
		m_num(100), m_limit(num_dim)
	{
		for (auto &i : m_limit) i.as.resize(size_pop);
		m_pos.resize(size_pop);
		for (auto &i : m_pos) i.resize(num_dim);
	}

	void AdaptorContGL::updateProbability(Environment *env,
		PopGL<Solution<>> &pop,
		const std::vector<Real> &weight,
		const std::vector<int> *index)
	{
		int popsize = pop.size();
		for (auto &j : m_proba) for (auto &i : j) {
			i.val = 0;
			i.idx.clear();
		}
		for (size_t j = 0; j < m_limit.size(); ++j) {
			m_limit[j].boundary.first = CAST_CONOP(env->problem())->range(j).second;
			m_limit[j].boundary.second = CAST_CONOP(env->problem())->range(j).first;
		}
		if (index == nullptr) {
			for (int j = 0; j < m_proba.size(); ++j) {
				for (int i = 0; i < popsize; ++i) {
					if (m_limit[j].boundary.first > pop[i].variable()[j])m_limit[j].boundary.first = pop[i].variable()[j];
					if (m_limit[j].boundary.second < pop[i].variable()[j])m_limit[j].boundary.second = pop[i].variable()[j];
				}
			}
		}
		else {
			for (int j = 0; j < m_proba.size(); ++j) {
				for (auto i : (*index)) {
					if (m_limit[j].boundary.first > pop[i].variable()[j])m_limit[j].boundary.first = pop[i].variable()[j];
					if (m_limit[j].boundary.second < pop[i].variable()[j])m_limit[j].boundary.second = pop[i].variable()[j];
				}
			}
		}
		if (index == nullptr) {
			for (int i = 0; i < m_proba.size(); ++i) {
				m_limit[i].range = m_limit[i].boundary.second - m_limit[i].boundary.first;
				m_limit[i].step = m_limit[i].range / m_num;
				if (m_limit[i].step > 0) {
					for (int j = 0; j < popsize; ++j) {
						if (!pop.isIndActive(j)) continue;
						int reg_idx = (pop[j].variable()[i] - m_limit[i].boundary.first) / m_limit[i].step;
						if (reg_idx == m_num) reg_idx--;
						m_proba[i][reg_idx].val += weight[j];
						m_proba[i][reg_idx].idx.push_back(j);

						m_pos[j][i] = reg_idx;
					}
				}
				else {
					for (int j = 0; j < popsize; ++j) {
						if (!pop.isIndActive(j)) continue;
						m_proba[i][0].val += weight[j];
						m_proba[i][0].idx.push_back(j);

						m_pos[j][i] = 0;
					}
				}
			}
		}
		else {
			for (int i = 0; i < m_proba.size(); ++i) {
				m_limit[i].range = m_limit[i].boundary.second - m_limit[i].boundary.first;
				m_limit[i].step = m_limit[i].range / m_num;
				if (m_limit[i].step > 0) {
					for (auto j : (*index)) {
						if (!pop.isIndActive(j)) continue;
						int reg_idx = (pop[j].variable()[i] - m_limit[i].boundary.first) / m_limit[i].step;
						if (reg_idx == m_num) reg_idx--;
						m_proba[i][reg_idx].val += weight[j];
						m_proba[i][reg_idx].idx.push_back(j);

						m_pos[j][i] = reg_idx;
					}
				}
				else {
					for (auto j : (*index)) {
						if (!pop.isIndActive(j)) continue;
						m_proba[i][0].val += weight[j];
						m_proba[i][0].idx.push_back(j);

						m_pos[j][i] = 0;
					}
				}
			}
		}
		accumlateProbability();
	}

	void AdaptorContGL::accumlateProbability() {
		for (auto i = 0; i < m_acc_proba.size(); ++i) {
			for (auto j = 0; j < m_acc_proba[i].size(); ++j) {
				m_acc_proba[i][j] = m_proba[i][j].val;
				if (j > 0) m_acc_proba[i][j] += m_acc_proba[i][j - 1];
			}
		}
	}

	void AdaptorContGL::createSolution(Environment *env, Random *rnd,
		PopGL<Solution<>> &pop,
		std::vector<Solution<>> &offspring)
	{
		Real p;
		for (int i = 0; i < pop.size(); i++) {
			p = rnd->uniform.next();
			if (p <= m_alpha)
				localSearch(env, rnd, i, pop, offspring);
			else
				globalSearch(env, rnd, i, pop, offspring);
		}
	}

	void AdaptorContGL::localSearch(Environment *env, Random *rnd, size_t i,
		PopGL<Solution<>> &pop,
		std::vector<Solution<>> &offspring)
	{
		Real center = 0, var = 0, x;
		const auto &range = CAST_CONOP(env->problem())->domain();
		for (int j = 0; j < offspring[i].variable().size(); j++) {
			center = pop[i].variable()[j];
			var = m_limit[j].as[i];
			if (center < range[j].limit.first) center = range[j].limit.first;
			else if (center > range[j].limit.second) center = range[j].limit.second;
			x = center + var * rnd->cauchy.next();
			if (x < range[j].limit.first) {
				x = (center + range[j].limit.first) / 2;
			}
			else if (x > range[j].limit.second) {
				x = (center + range[j].limit.second) / 2;
			}
			offspring[i].variable()[j] = x;
		}
	}

	void AdaptorContGL::globalSearch(Environment *env, Random *rnd, size_t i,
		PopGL<Solution<>> &pop,
		std::vector<Solution<>> &offspring)
	{
		Real center, var, x, p;
		const auto &range = CAST_CONOP(env->problem())->domain();
		for (int j = 0; j < offspring[i].variable().size(); j++) {
			center = 0, var = 0;
			std::vector<Real> env(m_num), acc_pro(m_num);
			for (int k = 0; k < m_num; ++k) {
				env[k] = m_proba[j][k].val * exp(-fabs(m_pos[i][j] - k) / m_num);
				acc_pro[k] = env[k];
				if (k > 0) acc_pro[k] += acc_pro[k - 1];
			}
			p = rnd->uniform.next() * acc_pro[m_num - 1];
			int idex = lower_bound(acc_pro.begin(), acc_pro.end(), p) - acc_pro.begin();
			if (m_proba[j][idex].idx.size() > 1) {
				for (auto k = 0; k < m_proba[j][idex].idx.size(); ++k) {
					center += pop[m_proba[j][idex].idx[k]].variable()[j];
				}
				center /= m_proba[j][idex].idx.size();
				for (auto k = 0; k < m_proba[j][idex].idx.size(); ++k) {
					var += (pop[m_proba[j][idex].idx[k]].variable()[j] - center) * (pop[m_proba[j][idex].idx[k]].variable()[j] - center);
				}
				var = sqrt(var / m_proba[j][idex].idx.size());
			}
			else {
				center = pop[i].variable()[j];
				var = m_limit[j].as[i];
			}
			if (center < range[j].limit.first) center = range[j].limit.first;
			else if (center > range[j].limit.second) center = range[j].limit.second;
			x = center + var * rnd->cauchy.next();
			if (x < range[j].limit.first) {
				x = (center + range[j].limit.first) / 2;
			}
			else if (x > range[j].limit.second) {
				x = (center + range[j].limit.second) / 2;
			}
			offspring[i].variable()[j] = x;
		}
	}

	int AdaptorContGL::updateSolution(Environment *env,
		PopGL<Solution<>> &pop,
		std::vector<Solution<>> &offspring, int &num_improve)
	{
		int rf = kNormalEval;
		int size_pop = pop.size();
		num_improve = 0;
		for (int i = 0; i < size_pop; i++) {
			if (!pop.isIndActive(i)) continue;
			rf = offspring[i].evaluate(env);
			if (dominate(offspring[i], pop[i], env->problem()->optimizeMode())) {
				pop[i] = offspring[i];
				pop.setImproved(i, true);
				num_improve++;
			}
			else {
				pop.setImproved(i, false);
			}
			if (rf & kTerminate) break;
		}
		return rf;
	}

	void AdaptorContGL::updateStep(Environment *env, PopGL<Solution<>> &pop) {
		int popsize = pop.size();
		int numDim = m_proba.size();
		std::vector<Real *> loc(popsize);
		std::vector<int> idx(popsize);
		for (int d = 0; d < numDim; ++d) {
			for (int i = 0; i < popsize; i++)
			{
				loc[i] = &pop[i].variable()[d];
			}
			mergeSort(loc, popsize, idx);
			Real lb = CAST_CONOP(env->problem())->range(d).first, ub = CAST_CONOP(env->problem())->range(d).second;

			for (int i = 0; i < popsize; i++)
			{
				if (i == 0) { //lower boundary
					m_limit[d].as[idx[i]] = *loc[idx[i + 1]] - *loc[idx[i]];
					//m_limit[d].as[idx[i]] = (*loc[idx[i + 1]] - lb) / 2;
				}
				else if (i == popsize - 1) { //upper boundary
					m_limit[d].as[idx[i]] = *loc[idx[i]] - *loc[idx[i - 1]];
					//m_limit[d].as[idx[i]] = (ub - *loc[idx[i - 1]]) / 2;
				}
				else {
					Real d1 = *loc[idx[i + 1]] - *loc[idx[i]], d2 = *loc[idx[i]] - *loc[idx[i - 1]];
					m_limit[d].as[idx[i]] = (d1 + d2) / 2;//d1 < d2 ? d1 : d2;
				}
			}
		}
	}
}