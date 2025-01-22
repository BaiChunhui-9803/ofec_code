#include "hts_cde.h"
#include "../../../../../../utility/linear_algebra/vector.h"
#include <numeric>

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void SC_CDE::addInputParameters() {
		m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
	}

	void SC_CDE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_seeds.clear();
		m_history_points.clear();
	}

	void SC_CDE::run_(Environment *env) {
		m_pop.reset(new PopCrowdingDE(m_pop_size, env));
		m_pop->initialize(env, m_random.get());
		m_pop->evaluate(env);
		for (size_t i = 0; i < m_pop_size; i++) {
			m_history_points.emplace_back(new Solution<>(m_pop->at(i)));
		}
		while (!terminating()) {
			speciation(env);
#ifdef OFEC_DATUM_SEEDS_H
			g_seeds.sols.clear();
			for (auto &seed : m_seeds) {
				g_seeds.sols.push_back(&seed);
			}
			datumUpdated(env, g_seeds);
#endif // OFEC_DATUM_SEEDS_H
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(2);
			for (size_t i = 0; i < m_pop_size; ++i) {
				g_multi_pop.pops[0].push_back(&m_pop->at(i));
			}
			for (auto &seed : m_seeds) {
				g_multi_pop.pops[1].push_back(&seed);
			}
			datumUpdated(env, g_multi_pop);
#endif
			m_pop->evolve(env, m_random.get());
			for (size_t i = 0; i < m_pop_size; i++) {
				m_history_points.emplace_back(new Solution<>(m_pop->at(i).trial()));
			}
			conserveSeeds(env);
		}
	}

	void SC_CDE::speciation(Environment *env) {
		m_seeds.clear();
		std::vector<size_t> order(m_pop->size());
		std::iota(order.begin(), order.end(), 0);
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMinimize) {
			std::sort(order.begin(), order.end(),
				[this](size_t i, size_t j) {
					return m_pop->at(i).objective(0) < m_pop->at(j).objective(0);
				}
			);
		}
		else {
			std::sort(order.begin(), order.end(),
				[this](size_t i, size_t j) {
					return m_pop->at(i).objective(0) > m_pop->at(j).objective(0);
				}
			);
		}
		for (size_t i : order) {
			bool found = false;
			for (auto &seed : m_seeds) {
				if (sameSpecies(seed, m_pop->at(i), m_history_points, env)) {
					found = true;
					break;
				}
			}
			if (!found) {
				m_seeds.push_back(m_pop->at(i));
			}
		}
	}

	void SC_CDE::conserveSeeds(Environment *env) {
		std::set<size_t> unmarked;
		for (size_t i = 0; i < m_pop_size; i++) {
			unmarked.insert(i);
		}
		for (auto& seed : m_seeds) {
			int p = -1;
			for (size_t i : unmarked) {
				if (sameSpecies(m_pop->at(i), seed, m_history_points, env)) {
					if (p == -1 || dominate(m_pop->at(p), m_pop->at(i), env->problem()->optimizeMode())) {
						p = i;
					}
				}
			}
			if (p != -1) {
				if (dominate(seed, m_pop->at(p), env->problem()->optimizeMode())) {
					dynamic_cast<Solution<>&>(m_pop->at(p)) = seed;
				}
			}
			else {
				for (size_t i : unmarked) {
					if (p == -1 || dominate(m_pop->at(p), m_pop->at(i), env->problem()->optimizeMode())) {
						p = i;
					}
				}
				dynamic_cast<Solution<>&>(m_pop->at(p)) = seed;
			}
			unmarked.erase(p);
		}
	}
	
	bool HTS_CDE::sameSpecies(const Solution<> &a, const Solution<> &b, 
		const std::list<std::shared_ptr<const Solution<>>> &H, Environment *env
) {
		return test(a, b, H, env);
	}

	bool HTS_CDE::test(const Solution<> &s1, const Solution<> &s2, const std::list<std::shared_ptr<const Solution<>>> &H, Environment *env) {
		std::list<std::shared_ptr<const Solution<>>> H_prime;
		std::vector<std::pair<Real, Real>> R = cubeRegion(s1, s2);
		for (auto &h : H) {
			if (isPointInOpenCubeRegion(*h, R)) {
				H_prime.push_back(h);
			}
		}
		if (H_prime.empty()) {
			return true;
		}
		else {
			auto m = mediumPoint(s1, s2, H_prime);
			if (dominate(s1, *m, env->problem()->optimizeMode()) && dominate(s2, *m, env->problem()->optimizeMode())) {
				return false;
			}
			else {
				return test(s1, *m, H_prime, env) && test(*m, s2, H_prime, env);
			}
		}
	}

	std::vector<std::pair<Real, Real>> HTS_CDE::cubeRegion(const Solution<> &s1, const Solution<> &s2) {
		std::vector<std::pair<Real, Real>> R(s1.variable().size());
		for (size_t j = 0; j < s1.variable().size(); ++j) {
			if (s1.variable()[j] < s2.variable()[j]) {
				R[j].first = s1.variable()[j];
				R[j].second = s2.variable()[j];
			}
			else {
				R[j].first = s2.variable()[j];
				R[j].second = s1.variable()[j];
			}
		}
		return R;
	}

	bool HTS_CDE::isPointInOpenCubeRegion(const Solution<> &p, const std::vector<std::pair<Real, Real>> &R) {
		for (size_t j = 0; j < R.size(); j++) {
			if (p.variable()[j] <= R[j].first || p.variable()[j] >= R[j].second) {
				return false;
			}
		}
		return true;
	}

	std::shared_ptr<const Solution<>> HTS_CDE::mediumPoint(const Solution<> &s1, const Solution<> &s2, 
		const std::list<std::shared_ptr<const Solution<>>> &H
	) {
		std::shared_ptr<const Solution<>> m;
		Real min_dis = std::numeric_limits<Real>::max();
		Vector x_a(s1.variable().vector()), x_b(s2.variable().vector());
		for (auto &h : H) {
			Vector x_m(h->variable().vector());
			Vector am = x_m - x_a, ab = x_b - x_a, bm = x_m - x_b, ba = x_a - x_b;
			Real norm_sa = pow(am.norm(), 2) * ab.norm() / (2 * am * ab);
			Real norm_tb = pow(bm.norm(), 2) * ba.norm() / (2 * bm * ba);
			Real dis = std::max(norm_sa, norm_tb);
			if (min_dis > dis) {
				min_dis = dis;
				m = h;
			}
		}
		return m;
	}
}
