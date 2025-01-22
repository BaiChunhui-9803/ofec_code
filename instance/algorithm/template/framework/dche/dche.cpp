#include "dche.h"
#include "../../../../../core/environment/environment.h"
#include "../../../../problem/continuous/free_peaks/free_peaks.h"
#include "../../../../../utility/clustering/nbc.h"
#include "../../../continuous/single_objective/multi_modal/hill_vallea/hill_vallea.h"
#include "../../../../../utility/functional.h"
#include <mlpack/core.hpp>
#include <mlpack/methods/kde/kde.hpp>
#include <algorithm>

namespace ofec {
	void DCHE::addInputParameters() {
		m_input_parameters.add("omniscient hills", new Bool(m_omniscient_hills, false));
		m_input_parameters.add("omniscient seeds", new Bool(m_omniscient_seeds, false));
		m_input_parameters.add("heuristic split", new Bool(m_heuristic_split, true));
	}

	void DCHE::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		m_seeds.clear();
		m_sp_tree.reset(new SPTree);
		m_sp_tree->setInitBox(CAST_CONOP(env->problem())->boundary());
		m_sp_tree->inputRatioData({ 1.0 });
		m_sp_tree->buildIndex();
		m_subspaces.clear();
		if (m_omniscient_hills) {
			if (CAST_FPs(env->problem()) == nullptr) {
				throw Exception("Omniscient hills are only available when solving the Free Peaks.");
			}
			else {
				std::map<size_t, size_t> add_id;
				m_sp_tree->addSubtree(0, *CAST_FPs(env->problem())->subspaceTree().tree, add_id);
				m_subspaces.resize(add_id.size());
				for (auto &pair : add_id) {
					m_subspaces[pair.second].reset(new Subspace(pair.second));
				}
			}
		}
		else {
			m_subspaces.emplace_back(new Subspace(0));
		}
		m_hills.clear();
		m_hills.emplace_back(new Hill);
		for (auto &subspace : m_subspaces) {
			m_hills.back()->addSubspace(subspace.get(), m_sp_tree.get(), true);

		}
		m_seeds.clear();
	}

	void DCHE::updateSeeds(Environment *env) {
		std::list<const Solution<> *> candidates;
		identifyCandidates(candidates, env);
		for (auto c : candidates) {
			size_t id_ssp = m_sp_tree->getRegionIdx(c->variable().vector());
			auto hill = m_subspaces.at(id_ssp)->assignedHills().front();
			bool flag_distinct = true;
			for (auto &s : m_seeds) {
				if (isSolutionInHill(s, hill)) {
					if (s.variableDistance(*c, env) < 1e-12 || HillVallEA::test(s, *c, 2, env)) {
						if (dominate(*c, s, env->problem()->optimizeMode())) {
							s = *c;
						}
						flag_distinct = false;
						break;
					}
				}
			}
			if (flag_distinct) {
				m_seeds.push_back(*c);
			}
		}
	}

	void DCHE::updateHills(Environment *env) {
		const size_t num_samples = env->problem()->numberVariables() * 100;
		NBC nbc;
		std::map<Hill *, std::list<size_t>> seeds_each_hill;
		for (size_t i = 0; i < m_seeds.size(); i++) {
			size_t id_subspace = m_sp_tree->getRegionIdx(m_seeds[i].variable().vector());
			for (auto hill : m_subspaces[id_subspace]->assignedHills()) {
				seeds_each_hill[hill].push_back(i);
			}
		}
		std::list<Hill *> new_hills;
		for (auto &hill : m_hills) {
			if (seeds_each_hill.count(hill.get()) == 0) {
				continue;
			}
			if (seeds_each_hill[hill.get()].size() > 1) {
				std::vector<std::unique_ptr<Solution<>>> samples;
				Solution<> *best_seed = nullptr;
				for (auto &id_seed : seeds_each_hill[hill.get()]) {
					if (best_seed == nullptr || dominate(m_seeds[id_seed], *best_seed, env->problem()->optimizeMode())) {
						best_seed = &m_seeds[id_seed];
					}
				}
				if (m_omniscient_hills) {
					std::map<size_t, std::vector<size_t>> smpls_each_subspace;
					for (size_t i = 0; i < num_samples; i++) {
						samples.emplace_back(dynamic_cast<Solution<>*>(env->problem()->createSolution()));
						do {
							randomSolutionInHill(*samples.back(), hill.get());
							samples.back()->evaluate(env, false);	// `effective = false` means cost no evaluations
						} while (dominate(*samples.back(), *best_seed, env->problem()->optimizeMode()));
						size_t id_subspace = m_sp_tree->getRegionIdx(samples.back()->variable().vector());
						smpls_each_subspace[id_subspace].push_back(i);
					}
					for (auto subspace : hill->subspaces()) {
						if (smpls_each_subspace.count(subspace->ID()) == 0) {
							auto &box = m_sp_tree->getBox(subspace->ID());
							auto new_sample = dynamic_cast<Solution<>*>(env->problem()->createSolution());
							for (size_t j = 0; j < box.size(); j++) {
								new_sample->variable()[j] = m_random->uniform.nextNonStd(box[j].first, box[j].second);
							}
							new_sample->evaluate(env, false);
							smpls_each_subspace[subspace->ID()].push_back(samples.size());
							samples.emplace_back(new_sample);
						}
					}
					std::vector<size_t> idx_centers;
					std::map<size_t, size_t> seed_subspace;	// index of seed among samples in subspace
					for (size_t i : seeds_each_hill[hill.get()]) {
						idx_centers.push_back(samples.size());
						size_t id_subspace = m_sp_tree->getRegionIdx(m_seeds[i].variable().vector());
						smpls_each_subspace[id_subspace].push_back(samples.size());
						seed_subspace[id_subspace] = samples.size();
						samples.emplace_back(new Solution<>(m_seeds[i]));
					}
					nbc.setData(samples);
					nbc.updateNbDistByKDTree(env);
					nbc.cutEdgesInGraph(idx_centers);
					nbc.updateClusters();
					auto &clusters = nbc.clusters();
					std::vector<size_t> cluster_each_sample(samples.size());
					for (size_t k = 0; k < clusters.size(); ++k) {
						for (size_t i : clusters[k]) {
							cluster_each_sample[i] = k;
						}
					}
					auto subspaces = hill->subspaces();
					hill->clear(m_sp_tree.get(), true);
					std::vector<Hill *> sub_hills = { hill.get() };
					for (size_t k = 1; k < clusters.size(); ++k) {
						new_hills.push_back(new Hill);
						sub_hills.push_back(new_hills.back());
					}
					for (auto subspace : subspaces) {
						if (seed_subspace.count(subspace->ID())) {
							size_t id_hill = cluster_each_sample[seed_subspace[subspace->ID()]];
							sub_hills[id_hill]->addSubspace(subspace, m_sp_tree.get(), true);
						}
						else {
							std::vector<size_t> num_each_cluster(clusters.size(), 0);
							for (size_t i : smpls_each_subspace[subspace->ID()]) {
								num_each_cluster[cluster_each_sample[i]]++;
							}
							size_t id_max_num = 0, max_num = num_each_cluster[0];
							for (size_t k = 1; k < clusters.size(); ++k) {
								if (max_num < num_each_cluster[k]) {
									id_max_num = k;
									max_num = num_each_cluster[k];
								}
							}
							sub_hills[id_max_num]->addSubspace(subspace, m_sp_tree.get(), true);
						}
					}
				}
				else {
					for (size_t i = 0; i < num_samples; i++) {
						samples.emplace_back(dynamic_cast<Solution<>*>(env->problem()->createSolution()));
						do {
							randomSolutionInHill(*samples.back(), hill.get());
							samples.back()->evaluate(env);
						} while (dominate(*samples.back(), *best_seed, env->problem()->optimizeMode()));
					}
					std::vector<size_t> idx_centers;
					for (size_t i : seeds_each_hill[hill.get()]) {
						idx_centers.push_back(samples.size());
						samples.emplace_back(new Solution<>(m_seeds[i]));
					}
					NBC nbc;
					nbc.setData(samples);
					nbc.updateNbDistByDistMat(env);
					nbc.cutEdgesInGraph(idx_centers);
					nbc.updateClusters();
					std::vector<size_t> sol_aff_clu(samples.size());
					for (size_t k = 0; k < nbc.clusters().size(); ++k) {
						for (size_t i : nbc.clusters()[k]) {
							sol_aff_clu[i] = k;
						}
					}
					std::map<size_t, std::set<size_t>> sols_in_ssp;
					for (size_t i = 0; i < samples.size(); ++i) {
						auto id_ssp = m_sp_tree->getRegionIdx(samples[i]->variable().vector());
						sols_in_ssp[id_ssp].insert(i);
					}
					std::map<size_t, int> ssp_aff_clu;
					auto ssps_aff_clu_uknown = hill->subspaces();
					while (!ssps_aff_clu_uknown.empty()) {
						auto subspace = ssps_aff_clu_uknown.front();
						if (sols_in_ssp[subspace->ID()].empty()) {
							ssp_aff_clu[subspace->ID()] = -1;
							ssps_aff_clu_uknown.remove(subspace);
						}
						else {
							while (true) {
								size_t num_total = 0;
								std::map<size_t, size_t> num_clu;
								for (size_t id_sol : sols_in_ssp[subspace->ID()]) {
									num_total++;
									if (num_clu.count(sol_aff_clu[id_sol])) {
										num_clu[sol_aff_clu[id_sol]]++;
									}
									else {
										num_clu[sol_aff_clu[id_sol]] = 1;
									}
								}
								int max_num_clu = -1; size_t max_num = 0;
								for (auto &p : num_clu) {
									if (max_num_clu == -1 || max_num < p.second) {
										max_num_clu = p.first;
										max_num = p.second;
									}
								}
								if (max_num == num_total) {
									ssp_aff_clu[subspace->ID()] = max_num_clu;
									ssps_aff_clu_uknown.remove(subspace);
									break;
								}
								else {
									Subspace *new_subspace = nullptr;
									if (m_heuristic_split) {
										std::vector<Solution<>*> sols_max_num_clu, sols_others;
										for (size_t id_sol : sols_in_ssp[subspace->ID()]) {
											if (sol_aff_clu[id_sol] == max_num_clu) {
												sols_max_num_clu.push_back(samples[id_sol].get());
											}
											else {
												sols_others.push_back(samples[id_sol].get());
											}
										}
										auto result = dimensionWithSmallestIntersection(sols_max_num_clu, sols_others, env->problem()->numberVariables());
										new_subspace = subspace->bisect(m_sp_tree.get(), &std::get<0>(result), &std::get<1>(result));
									}
									else {
										new_subspace = subspace->bisect(m_sp_tree.get());
									}
									ssps_aff_clu_uknown.push_back(new_subspace);
									m_subspaces.emplace_back(new_subspace);
									for (auto it = sols_in_ssp[subspace->ID()].begin(); it != sols_in_ssp[subspace->ID()].end();) {
										if (m_sp_tree->getRegionIdx(samples[*it]->variable().vector()) != subspace->ID()) {
											sols_in_ssp[new_subspace->ID()].insert(*it);
											it = sols_in_ssp[subspace->ID()].erase(it);
										}
										else {
											it++;
										}
									}
									if (sols_in_ssp[subspace->ID()].empty()) {
										ssp_aff_clu[subspace->ID()] = -1;
										ssps_aff_clu_uknown.remove(subspace);
										break;
									}
								}
							}
						}
					}
					std::map<size_t, Hill *> clu_to_hill;
					clu_to_hill[0] = hill.get();
					for (size_t k = 1; k < nbc.clusters().size(); ++k) {
						auto new_hill = new Hill;
						clu_to_hill[k] = new_hill;
						new_hills.push_back(new_hill);
					}
					std::list<Subspace *> subspaces_to_group;
					auto hill_kept_subspaces = hill->subspaces();
					for (Subspace *subspace : hill_kept_subspaces) {
						subspace->quitHill(hill.get(), m_sp_tree.get(), true);
						//subspace->clearAssignedHills(m_sp_tree.get(), true);
						subspaces_to_group.push_back(subspace);
					}
					std::set<size_t> empty_subspace_ids;
					for (Subspace *subspace : subspaces_to_group) {
						if (ssp_aff_clu[subspace->ID()] > -1) {
							subspace->joinHill(clu_to_hill[ssp_aff_clu[subspace->ID()]], m_sp_tree.get(), true);
						}
						else {
							empty_subspace_ids.insert(subspace->ID());
						}
					}
					while (!empty_subspace_ids.empty()) {
						std::map<size_t, std::list<Hill *>> border;
						for (size_t id_subspace : empty_subspace_ids) {
							for (Subspace *adj_subspace : m_subspaces[id_subspace]->adjacentSubspaces()) {
								for (Hill *hill : adj_subspace->assignedHills()) {
									if (std::find(border[id_subspace].begin(), border[id_subspace].end(), hill) == border[id_subspace].end()) {
										border[id_subspace].push_back(hill);
									}
								}
							}
						}
						for (auto &p : border) {
							for (Hill *hill : p.second) {
								m_subspaces[p.first]->joinHill(hill, m_sp_tree.get(), true);
							}
							empty_subspace_ids.erase(p.first);
						}
					}
				}
			}
		}
		for (auto hill : new_hills) {
			if (hill->subspaces().empty()) {
				throw Exception("Empty hill region.");
			}
			m_hills.emplace_back(hill);
		}
	}

	std::tuple<int, Real> DCHE::dimensionWithSmallestIntersection(const std::vector<Solution<>*> &sols1, const std::vector<Solution<>*> &sols2, size_t num_dims) {
		using namespace mlpack;
		arma::mat class1(num_dims, sols1.size()), class2(num_dims, sols2.size());
		for (size_t i = 0; i < sols1.size(); ++i) {
			class1.col(i) = arma::vec(sols1[i]->variable().vector().data(), num_dims, false);
		}
		for (size_t i = 0; i < sols2.size(); ++i) {
			class2.col(i) = arma::vec(sols2[i]->variable().vector().data(), num_dims, false);
		}
		std::vector<Real> intersections(num_dims), pivot(num_dims);
		for (size_t j = 0; j < num_dims; ++j) {
			double stddev = arma::stddev(arma::join_rows(class1.row(j), class2.row(j)));
			if (stddev == 0) {
				intersections[j] = 1.0;
				pivot[j] = arma::mean(arma::join_rows(class1.row(j), class2.row(j)));
				continue;
			}
			double bandwith = 1.06 * stddev * pow(sols1.size() + sols2.size(), -0.2);
			KDE<GaussianKernel, EuclideanDistance, arma::mat> kde1(0.05, 0.0, GaussianKernel(bandwith)), kde2(0.05, 0.0, GaussianKernel(bandwith));
			kde1.Train(class1.row(j));
			kde2.Train(class2.row(j));
			Real min_val = std::min(class1.row(j).min(), class2.row(j).min());
			Real max_val = std::max(class1.row(j).max(), class2.row(j).max());
			arma::mat query = arma::linspace(min_val, max_val, sols1.size() + sols2.size()).t();
			arma::vec estimations1, estimations2;
			kde1.Evaluate(query, estimations1);
			kde2.Evaluate(query, estimations2);
			arma::vec pdf1 = estimations1 / arma::accu(estimations1);
			arma::vec pdf2 = estimations2 / arma::accu(estimations2);
			intersections[j] = arma::accu(arma::min(pdf1, pdf2));
			arma::vec diff = pdf1 - pdf2;
			int index = -1;
			for (size_t i = 0; i < diff.n_elem - 1; ++i) {
				if (diff(i) * diff(i + 1) <= 0) {
					index = i;
					break;
				}
			}
			if (index == -1) {
				std::cout << "query:" << std::endl;
				std::cout << query << std::endl;
				std::cout << "estimations1:" << std::endl;
				std::cout << estimations1.t() << std::endl;
				std::cout << "estimations2:" << std::endl;
				std::cout << estimations2.t() << std::endl;
				std::cout << "pdf1:" << std::endl;
				std::cout << pdf1.t() << std::endl;
				std::cout << "pdf2:" << std::endl;
				std::cout << pdf2.t() << std::endl;
			}
			pivot[j] = (query(0, index) + query(0, index + 1)) / 2;
		}
		int min_dim = 0;
		Real min_intersection = intersections[0];
		for (size_t j = 1; j < num_dims; ++j) {
			if (intersections[j] < min_intersection) {
				min_intersection = intersections[j];
				min_dim = j;
			}
		}
		return { min_dim, pivot[min_dim] };
	}

	void DCHE::randomVariablesInHill(VariableVector<> &var, const Hill *hill) const {
		size_t id_ssp = hill->rouletteWheelSelection(m_random.get(), m_sp_tree.get())->ID();
		auto &box = m_sp_tree->getBox(id_ssp);
		for (size_t j = 0; j < box.size(); j++) {
			var[j] = m_random->uniform.nextNonStd(box[j].first, box[j].second);
		}
	}

	void DCHE::randomSolutionInHill(Solution<> &sol, const Hill *hill) const {
		randomVariablesInHill(sol.variable(), hill);
	}

	bool DCHE::isVariablesInHill(const VariableVector<> &var, const Hill *hill) const {
		size_t id_ssp = m_sp_tree->getRegionIdx(var.vector());
		auto &hills = m_subspaces.at(id_ssp)->assignedHills();
		return std::find(hills.begin(), hills.end(), hill) != hills.end();
	}

	bool DCHE::isSolutionInHill(const Solution<> &sol, const Hill *hill) const {
		return isVariablesInHill(sol.variable(), hill);
	}

	void DCHE::Hill::addSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace) {
		if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) != m_subspaces.end())
			return;
		m_subspaces.push_back(subspace);
		m_volume += sp_tree->getBoxVolume(subspace->ID());
		if (update_subspace) {
			subspace->joinHill(this, sp_tree);
		}
	}

	void DCHE::Hill::removeSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace) {
		if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) == m_subspaces.end())
			return;
		m_subspaces.remove(subspace);
		m_volume -= sp_tree->getBoxVolume(subspace->ID());
		if (update_subspace) {
			subspace->quitHill(this, sp_tree);
		}
	}

	void DCHE::Hill::merge(Hill *hill, const SPTree *sp_tree, bool update_subspace) {
		if (hill == this)
			return;
		for (auto subspace : hill->m_subspaces) {
			addSubspace(subspace, sp_tree, update_subspace);
		}
		hill->clear(sp_tree, update_subspace);
	}

	void DCHE::Hill::clear(const SPTree *sp_tree, bool update_subspace) {
		if (update_subspace) {
			for (auto subspace : m_subspaces) {
				subspace->quitHill(this, sp_tree);
			}
		}
		m_subspaces.clear();
		m_volume = 0;
	}

	DCHE::Subspace *DCHE::Hill::rouletteWheelSelection(Random *rnd, const SPTree *sp_tree) const {
		if (m_subspaces.empty())
			return nullptr;
		Real rand_pos = m_volume * rnd->uniform.next();
		Real accum = 0;
		for (auto subspace : m_subspaces) {
			accum += sp_tree->getBoxVolume(subspace->ID());
			if (rand_pos <= accum)
				return subspace;
		}
		return nullptr;
	}

	void DCHE::Subspace::joinHill(Hill *hill, const SPTree *sp_tree, bool update_hill) {
		if (std::find(m_assigned_hills.begin(), m_assigned_hills.end(), hill) != m_assigned_hills.end())
			return;
		m_assigned_hills.push_back(hill);
		if (update_hill) {
			hill->addSubspace(this, sp_tree);
		}
	}

	void DCHE::Subspace::quitHill(Hill *hill, const SPTree *sp_tree, bool update_hill) {
		if (std::find(m_assigned_hills.begin(), m_assigned_hills.end(), hill) == m_assigned_hills.end())
			return;
		m_assigned_hills.remove(hill);
		if (update_hill) {
			hill->removeSubspace(this, sp_tree);
		}
	}

	void DCHE::Subspace::clearAssignedHills(const SPTree *sp_tree, bool update_hill) {
		if (update_hill) {
			for (Hill *hill : m_assigned_hills)
				hill->removeSubspace(this, sp_tree);
		}
		m_assigned_hills.clear();
	}

	DCHE::Subspace* DCHE::Subspace::bisect(SPTree *sp_tree, int *dim, Real *pivot) {
		size_t id_new_subspace = sp_tree->splitRegion(m_id, dim, pivot);
		auto new_subspace = new Subspace(id_new_subspace);
		new_subspace->m_assigned_hills = m_assigned_hills;
		new_subspace->m_adjacent_subspaces = m_adjacent_subspaces;
		for (Hill *hill : m_assigned_hills) {
			hill->m_subspaces.push_back(new_subspace);
		}
		for (auto it = m_adjacent_subspaces.begin(); it != m_adjacent_subspaces.end();) {
			if (!sp_tree->checkAdjacency((*it)->m_id, m_id)) {
				(*it)->m_adjacent_subspaces.remove(this);
				it = this->m_adjacent_subspaces.erase(it);
			}
			else
				it++;
		}
		m_adjacent_subspaces.push_back(new_subspace);
		for (auto it = new_subspace->m_adjacent_subspaces.begin(); it != new_subspace->m_adjacent_subspaces.end();) {
			if (sp_tree->checkAdjacency((*it)->m_id, new_subspace->m_id)) {
				(*it)->m_adjacent_subspaces.push_back(new_subspace);
				it++;
			}
			else
				it = new_subspace->m_adjacent_subspaces.erase(it);
		}
		new_subspace->m_adjacent_subspaces.push_back(this);
		return new_subspace;
	}
}
