#include "../../../utility/catch.hpp"
#include "../../../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../../../instance/algorithm/continuous/single_objective/multi_modal/de_nrand_1/de_nrand_1_pop.h"

using namespace ofec;

using SPTree = nanoflann::KDTreeSpace<Real>;

class Subspace;
class Hill {
	friend class Subspace;
protected:
	std::list<Subspace*> m_subspaces;
	Real m_volume = 0;
public:
	void addSubspace(Subspace* subspace, const SPTree* sp_tree, bool update_subspace = false);
	void removeSubspace(Subspace* subspace, const SPTree* sp_tree, bool update_subspace = false);
	void merge(Hill* hill, const SPTree* sp_tree, bool update_subspace = false);
	void clear(const SPTree* sp_tree, bool update_subspace = false);
	Subspace* rouletteWheelSelection(Random* rnd, const SPTree* sp_tree) const;
	const std::list<Subspace*>& subspaces() const { return m_subspaces; }
	Real volume() const { return m_volume; }
};

class Subspace {
protected:
	std::list<Hill*> m_assigned_hills;
	std::list<Subspace*> m_adjacent_subspaces;
public:
	const size_t m_id;
	Subspace(size_t id) : m_id(id) {}
	void joinHill(Hill* hill, const SPTree* sp_tree, bool update_hill = false);
	void quitHill(Hill* hill, const SPTree* sp_tree, bool update_hill = false);
	Subspace* halve(SPTree* sp_tree);
	void clearAssignedHills(const SPTree* sp_tree, bool update_hill = false);
	const std::list<Hill*>& assignedHills() const { return m_assigned_hills; }
	const std::list<Subspace*>& adjacentSubspaces() const { return m_adjacent_subspaces; }
};

void Hill::addSubspace(Subspace* subspace, const SPTree* sp_tree, bool update_subspace) {
	if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) != m_subspaces.end())
		return;
	m_subspaces.push_back(subspace);
	m_volume += sp_tree->getBoxVolume(subspace->m_id);
	if (update_subspace) {
		subspace->joinHill(this, sp_tree);
	}
}

void Hill::removeSubspace(Subspace* subspace, const SPTree* sp_tree, bool update_subspace) {
	if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) == m_subspaces.end())
		return;
	m_subspaces.remove(subspace);
	m_volume -= sp_tree->getBoxVolume(subspace->m_id);
	if (update_subspace) {
		subspace->quitHill(this, sp_tree);
	}
}

void Hill::merge(Hill* hill, const SPTree* sp_tree, bool update_subspace) {
	for (auto subspace : hill->m_subspaces) {
		addSubspace(subspace, sp_tree, update_subspace);
	}
	hill->clear(sp_tree, update_subspace);
}

void Hill::clear(const SPTree* sp_tree, bool update_subspace) {
	if (update_subspace) {
		for (auto subspace : m_subspaces) {
			subspace->quitHill(this, sp_tree);
		}
	}
	m_subspaces.clear();
	m_volume = 0;
}

Subspace* Hill::rouletteWheelSelection(Random* rnd, const SPTree* sp_tree) const {
	if (m_subspaces.empty())
		return nullptr;
	Real rand_pos = m_volume * rnd->uniform.next();
	Real accum = 0;
	for (auto subspace : m_subspaces) {
		accum += sp_tree->getBoxVolume(subspace->m_id);
		if (rand_pos <= accum)
			return subspace;
	}
	return nullptr;
}

void Subspace::joinHill(Hill* hill, const SPTree* sp_tree, bool update_hill) {
	if (std::find(m_assigned_hills.begin(), m_assigned_hills.end(), hill) != m_assigned_hills.end())
		return;
	m_assigned_hills.push_back(hill);
	if (update_hill) {
		hill->addSubspace(this, sp_tree);
	}
}

void Subspace::quitHill(Hill* hill, const SPTree* sp_tree, bool update_hill) {
	if (std::find(m_assigned_hills.begin(), m_assigned_hills.end(), hill) == m_assigned_hills.end())
		return;
	m_assigned_hills.remove(hill);
	if (update_hill) {
		hill->removeSubspace(this, sp_tree);
	}
}

Subspace* Subspace::halve(SPTree* sp_tree) {
	size_t id_new_subspace = sp_tree->splitRegion(m_id);
	auto new_subspace = new Subspace(id_new_subspace);
	new_subspace->m_assigned_hills = m_assigned_hills;
	new_subspace->m_adjacent_subspaces = m_adjacent_subspaces;
	for (Hill* hill : m_assigned_hills) {
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

void Subspace::clearAssignedHills(const SPTree* sp_tree, bool update_hill) {
	if (update_hill) {
		for (Hill* hill : m_assigned_hills)
			hill->removeSubspace(this, sp_tree);
	}
	m_assigned_hills.clear();
}

void runDENrand1onFP(size_t num_vars, size_t pop_size) {
	ParameterMap param;
	param["problem name"] = std::string("free_peaks");
	param["number of variables"] = int(num_vars);

	std::shared_ptr<const ofec::ParameterMap> cparam(new ofec::ParameterMap(param));
	auto problem = Problem::generateByFactory(cparam, 0.5);
	problem->initialize();

	std::vector<Solution<>> sols;

	Random rand(0.5);
	PopDE_nrand_1 pop(pop_size, problem.get());
	pop.initialize(problem.get(), &rand);
	auto optima = CAST_CONOP(problem.get())->optima();
	for (size_t i = 0; i < optima->numberSolutions(); ++i)
		pop[i].variable() = optima->solution(i).variable();
	pop.evaluate(problem.get(), nullptr);
	for (size_t i = 0; i < pop_size; ++i)
		sols.emplace_back(pop[i]);
	for (size_t iter = 0; iter < num_vars * 10; ++iter) {
		pop.evolve(problem.get(), nullptr, &rand);
		for (size_t i = 0; i < pop_size; ++i)
			sols.emplace_back(pop[i]);
	}

	auto fp_tree = dynamic_cast<FreePeaks*>(problem.get())->subspaceTree().tree.get();
	std::vector<std::vector<size_t>> clusters(fp_tree->size());
	std::vector<size_t> sol_aff_clu(sols.size());
	for (size_t i = 0; i < sols.size(); ++i) {
		auto id_peak = fp_tree->getRegionIdx(sols[i].variable().vect());
		clusters[id_peak].push_back(i);
		sol_aff_clu[i] = id_peak;
	}
	
	auto sp_tree = std::make_unique<nanoflann::KDTreeSpace<double>>();
	sp_tree->setInitBox(fp_tree->getInitBox());
	sp_tree->inputRatioData(std::vector<Real>(1, 1.0));
	sp_tree->buildIndex();

	std::map<size_t, std::set<size_t>> sols_in_ssp;
	for (size_t i = 0; i < sols.size(); ++i) {
		auto id_ssp = sp_tree->getRegionIdx(sols[i].variable().vect());
		sols_in_ssp[id_ssp].insert(i);
	}
	std::map<size_t, int> ssp_aff_clu;

	std::vector<std::unique_ptr<Subspace>> subspaces;
	std::list<std::unique_ptr<Hill>> hills;
	subspaces.emplace_back(new Subspace(0));
	hills.emplace_back(new Hill);
	hills.front()->addSubspace(subspaces.front().get(), sp_tree.get(), true);

	auto hill_kept = hills.front().get();

	auto ssps_aff_clu_uknown = hill_kept->subspaces();
	while (!ssps_aff_clu_uknown.empty()) {
		auto subspace = ssps_aff_clu_uknown.front();
		if (sols_in_ssp[subspace->m_id].empty()) {
			ssp_aff_clu[subspace->m_id] = -1;
			ssps_aff_clu_uknown.remove(subspace);
		}
		else {
			while (true) {
				size_t num_total = 0;
				std::map<size_t, size_t> num_clu;
				for (size_t id_sol : sols_in_ssp[subspace->m_id]) {
					num_total++;
					if (num_clu.count(sol_aff_clu[id_sol]))
						num_clu[sol_aff_clu[id_sol]]++;
					else
						num_clu[sol_aff_clu[id_sol]] = 1;
				}
				bool split = true;
				for (auto& p : num_clu) {
					if (p.second >= 1.0 * num_total) {
						split = false;
						ssp_aff_clu[subspace->m_id] = p.first;
						ssps_aff_clu_uknown.remove(subspace);
						break;
					}
				}
				if (split) {
					auto new_subspace = subspace->halve(sp_tree.get());
					ssps_aff_clu_uknown.push_back(new_subspace);
					subspaces.emplace_back(new_subspace);
					for (auto it = sols_in_ssp[subspace->m_id].begin(); it != sols_in_ssp[subspace->m_id].end();) {
						if (sp_tree->getRegionIdx(sols[*it].variable().vect()) != subspace->m_id) {
							sols_in_ssp[new_subspace->m_id].insert(*it);
							it = sols_in_ssp[subspace->m_id].erase(it);
						}
						else
							it++;
					}
					if (sols_in_ssp[subspace->m_id].empty()) {
						ssp_aff_clu[subspace->m_id] = -1;
						ssps_aff_clu_uknown.remove(subspace);
						break;
					}
				}
				else {
					break;
				}
			}
		}
	}
	//std::map<size_t, Hill*> clu_to_hill;
	//clu_to_hill[0] = hill_kept;
	//for (size_t k = 1; k < m_clusters.size(); ++k) {
	//	auto new_hill = new Hill;
	//	clu_to_hill[k] = new_hill;
	//	m_hills.emplace_back(new_hill);
	//}
	//std::list<Subspace*> subspaces_to_group;
	//auto hill_kept_subspaces = hill_kept->subspaces();
	//for (Subspace* subspace : hill_kept_subspaces) {
	//	subspace->clearAssignedHills(sp_tree.get(), true);
	//	subspaces_to_group.push_back(subspace);
	//}
	//std::set<size_t> empty_subspace_ids;
	//for (Subspace* subspace : subspaces_to_group) {
	//	if (ssp_aff_clu[subspace->m_id] > -1)
	//		subspace->joinHill(clu_to_hill[ssp_aff_clu[subspace->m_id]], sp_tree.get(), true);
	//	else
	//		empty_subspace_ids.insert(subspace->m_id);
	//}
	//while (!empty_subspace_ids.empty()) {
	//	std::map<size_t, std::list<Hill*>> border;
	//	for (size_t id_subspace : empty_subspace_ids) {
	//		for (Subspace* adj_subspace : subspaces[id_subspace]->adjacentSubspaces()) {
	//			for (Hill* hill : adj_subspace->assignedHills()) {
	//				if (std::find(border[id_subspace].begin(), border[id_subspace].end(), hill) == border[id_subspace].end())
	//					border[id_subspace].push_back(hill);
	//			}
	//		}
	//	}
	//	for (auto& p : border) {
	//		for (Hill* hill : p.second)
	//			subspaces[p.first]->joinHill(hill, sp_tree.get(), true);
	//		empty_subspace_ids.erase(p.first);
	//	}
	//}
}

TEST_CASE("Test impact of dimensinality on update of hills", "[EDHE][SP]") {

}