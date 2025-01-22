#include "ande_pop.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	PopANDE::PopANDE(size_t size_pop, Environment *env, Real lambda, size_t Mits, size_t Cits) :
		PopulationDE<>(size_pop, env) , m_apc(lambda, Mits, Cits) {}

	void PopANDE::clustering(Environment *env) {
		m_apc.updateData(*this);
		m_apc.clustering(env->problem());
	}

	int PopANDE::evolve(Environment *env, Random *rnd) {
		std::vector<size_t> ridx;
		for (auto cluster : m_apc.clusters()) {
			if (cluster.size() < 4) 
				continue;
			for (size_t i : cluster) {
				selectInCandidates(3, cluster, ridx, rnd);
				m_individuals[i]->mutate(
					m_scaling_factor, 
					m_individuals[ridx[0]].get(), 
					m_individuals[ridx[1]].get(), 
					m_individuals[ridx[2]].get(), 
					env
				);
				recombine(i, rnd, env);
				m_individuals[i]->trial().evaluate(env);
				size_t idx_nearest = cluster[0];
				Real min_dis = m_individuals[i]->trial().variableDistance(*m_individuals[idx_nearest], env);
				Real temp_dis;
				for (size_t j = 1; j < cluster.size(); ++j) {
					temp_dis = m_individuals[i]->trial().variableDistance(*m_individuals[cluster[j]], env);
					if (min_dis > temp_dis) {
						min_dis = temp_dis;
						idx_nearest = cluster[j];
					}
				}
				if (dominate(m_individuals[i]->trial(), *m_individuals[idx_nearest], env->problem()->optimizeMode())) {
					m_individuals[idx_nearest]->solution() = m_individuals[i]->trial();
				}
			}
			CPA(cluster, env);
		}
		TLLS(m_apc.clusters(), env, rnd);
		return kNormalEval;
	}

	void PopANDE::CPA(const std::vector<size_t> &cluster, Environment *env) {
		if (cluster.size() < 4) return;
		/* Find lbest */
		size_t lbest = cluster[0];
		for (size_t i = 0; i < cluster.size(); i++) {
			if (dominate(*m_individuals[cluster[i]] , *m_individuals[lbest], env->problem()->optimizeMode()))
				lbest = cluster[i];
		}
		/* Find its neighbors */
		std::vector<Real> dis_to_lbest(cluster.size());
		std::vector<int> sequence;
		std::vector<int> neighbors;
		for (size_t i = 0; i < cluster.size(); i++) {
			dis_to_lbest[i] = m_individuals[lbest]->variableDistance(*m_individuals[cluster[i]], env);
		}
		mergeSort(dis_to_lbest, dis_to_lbest.size(), sequence);
		for (size_t i = 1; i < sequence.size(); ++i) {
			neighbors.emplace_back(sequence[i]);
			if (neighbors.size() >= 5)
				break;
		}
		/* Determine the contour value f */
		Real f, f_lbest = m_individuals[lbest]->objective(0), f_i;
		if (env->problem()->optimizeMode(0) == OptimizeMode::kMaximize)
			f = f_lbest + 0.2 * abs(f_lbest) + 0.1;
		else
			f = f_lbest - 0.2 * abs(f_lbest) - 0.1;
		/* Calculate the interpolated points */
		size_t D = env->problem()->numberVariables();
		const auto &x_lbest = m_individuals[lbest]->variable();
		std::vector<std::vector<Real>> inters(neighbors.size(), std::vector<Real>(D));
		for (size_t i = 0; i < inters.size(); i++) {
			const auto &x_i = m_individuals[neighbors[i]]->variable();
			f_i = m_individuals[neighbors[i]]->objective(0);
			for (size_t j = 0; j < D; j++) {
				inters[i][j] = x_lbest[j] + (f - f_lbest) / (f_i - f_lbest) * (x_i[j] - x_lbest[j]);
			}
		}
		/* Estimate the potential optima */
		std::vector<Real> optima(D, 0);
		for (size_t j = 0; j < D; j++) {
			for (const auto &inter : inters) {
				optima[j] += inter[j];
			}
			optima[j] /= inters.size();
		}
		/* Compete the potential optima with the lbest */
		Solution<> temp_sol(1, 0, D);
		temp_sol.variable().vector() = optima;
		env->problem()->validateSolution(temp_sol, Validation::kSetToBound);
		temp_sol.evaluate(env);
		if (dominate(temp_sol , *m_individuals[lbest], env->problem()->optimizeMode())) {
			//std::cout << "CPA works" << std::endl;
			m_individuals[lbest]->solution() = temp_sol;
		}
	}

	void PopANDE::TLLS(const std::vector<std::vector<size_t>> &clusters, Environment *env, Random *rnd) {
		/* Generate the sample standard deviation sigma */
		size_t D = env->problem()->numberVariables();
		Real sigma = pow(10, (-1 - (10 / (Real)D + 3) * env->algorithm()->evaluations() / m_MaxFEs));
		/* Calculate the niche-level search probability P */
		size_t n = clusters.size();
		std::vector<Real> P(n);
		std::vector<Real> f_niche_seed(n);
		for (size_t i = 0; i < n; i++) {
			size_t lbest = clusters[i][0];
			for (size_t k = 1; k < clusters[i].size(); k++) {
				if (dominate(*m_individuals[clusters[i][k]], *m_individuals[lbest], env->problem()->optimizeMode())) {
					lbest = clusters[i][k];
					f_niche_seed[i] = m_individuals[lbest]->objective(0);
				}
			}
		}
		std::vector<int> sequence;
		bool ascending = env->problem()->optimizeMode(0) == OptimizeMode::kMinimize;
		ascending = !ascending; // Add a exclamatory mark for sorting from worse to better
		mergeSort(f_niche_seed, n, sequence, ascending);
		std::vector<size_t> r(n);
		for (size_t i = 0; i < n; i++)
			r[sequence[i]] = i + 1;
		for (size_t i = 0; i < n; i++)
			P[i] = (Real)r[i] / n;
		/*  */
		Real rand;
		for (size_t i = 0; i < n; i++) {
			rand = rnd->uniform.next();
			if (rand < P[i]) {
				/*  Calculate the Solution-level local search probability */
				size_t n_i = clusters[i].size();
				std::vector<Real> P_i(n_i);
				std::vector<Real> f_inds_in_niche(n_i);
				for (size_t k = 0; k < n_i; k++)
					f_inds_in_niche[k] = m_individuals[clusters[i][k]]->objective(0);
				mergeSort(f_inds_in_niche, n_i, sequence, ascending);
				r.resize(n_i);
				for (size_t k = 0; k < n_i; k++)
					r[sequence[k]] = i + 1;
				for (size_t k = 0; k < n_i; k++)
					P_i[k] = (Real)r[k] / n_i;
				/* */
				for (size_t k = 0; k < n_i; k++) {
					rand = rnd->uniform.next();
					if (rand < P_i[k]) {
						const auto &x_ik = m_individuals[clusters[i][k]]->variable();
						Solution<> s1(1, 0, D), s2(1, 0, D);
						for (size_t j = 0; j < D; j++) {
							s1.variable()[j] = rnd->normal.nextNonStd(x_ik[j], sigma);
							s2.variable()[j] = rnd->normal.nextNonStd(x_ik[j], sigma);
						}
						env->problem()->validateSolution(s1, Validation::kSetToBound);
						env->problem()->validateSolution(s2, Validation::kSetToBound);
						s1.evaluate(env);
						s2.evaluate(env);
						s1 = dominate(s1, s2, env->problem()->optimizeMode()) ? s1 : s2;
						if (dominate(s1, *m_individuals[clusters[i][k]], env->problem()->optimizeMode())) {
							//std::cout << "TLLS works" << std::endl;
							m_individuals[clusters[i][k]]->solution() = s1;
						}
					}
				}

			}
		}
	}
}