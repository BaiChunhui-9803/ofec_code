#include "hghec.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/environment/environment.h"
#include "../../../../../../utility/functional.h"
#include "../../../../../../utility/clustering/nbc.h"
#include "../../../../../../utility/nondominated_sorting/fast_sort.h"
#include <numeric>
#include <queue>

namespace ofec {
	void HGHEC::addInputParameters() {
		m_input_parameters.add("use acceleration mode", new Bool(m_use_acceleration_mode, false));
		m_input_parameters.add("threshold number of solutions", new RangedSizeT(m_num_sols_thold, 5, 1000, 100));
		m_input_parameters.add("interval number of solutions", new RangedSizeT(m_num_sols_intvl, 1, 1000, 100));
	}

	void HGHEC::initialize_(Environment *env) {
		HGHE::initialize_(env);
		m_sp_tree.reset();
		m_subspaces.clear();
		m_hills.clear();
		m_id_to_ptr.clear();
		m_ptr_to_id.clear();
		m_phi = 6.0;
		m_keep_subspace_valley = true;
		m_use_his_explore_only = true;
	}

	void HGHEC::initHills(Environment *env) {
		m_sp_tree.reset(new SPTree);
		m_sp_tree->setInitBox(CAST_CONOP(env->problem())->boundary());
		m_sp_tree->inputRatioData(std::vector<Real>(1, 1.0));
		m_sp_tree->buildIndex();
		m_subspaces.emplace_back(new Subspace(0));
		m_hills.emplace_back(new Hill);
		m_hills.back()->addSubspace(m_subspaces.back().get(), m_sp_tree.get(), true);
		m_id_to_ptr[0] = m_hills.back().get();
		m_ptr_to_id[m_hills.back().get()] = 0;
	}

	void HGHEC::updateHills(Environment *env) {
		if (m_his_sols.empty()) {
			return;
		}
		for (auto &hill : m_hills) {
			if (m_use_acceleration_mode && !isTimeToLearn(hill.get())) {
				continue;
			}
			std::vector<const SolutionType*> his_sols;
			getHisSolsInHill(his_sols, hill.get());
			std::vector<const SolutionType *> seed_sols;
			identifySeedSolutions(his_sols, seed_sols, env);
			subdivideSpace(seed_sols, env);
			if (m_use_acceleration_mode) {
				updateNumLearn(seed_sols.size(), hill.get());
			}
		}
		assignSubspaceFitness();
		groupSubspace(env);
	}

	void HGHEC::randomVarInHill(VariableVector<Real> &var, Hill *hill) {
		Subspace *subspace = hill->rouletteWheelSelection(m_random.get(), m_sp_tree.get());
		auto &box = m_sp_tree->getBox(subspace->id);
		for (size_t j = 0; j < box.size(); j++) {
			var[j] = m_random->uniform.nextNonStd(box[j].first, box[j].second);
		}
	}

	bool HGHEC::isVarInHill(const VariableVector<Real>& var, Hill *hill) const {
		size_t id_ssp = m_sp_tree->getRegionIdx(var.vector());
		auto &hills = m_subspaces.at(id_ssp)->setassignedHills();
		return hills.count(hill);
	}

	const std::list<HGHEC::Hill*>& HGHEC::theHillsLocated(const VariableVector<Real> &var) const {
		return m_subspaces[m_sp_tree->getRegionIdx(var.vector())]->assignedHills();
	}

	void HGHEC::randomSolInHill(SolutionType &sol, Hill *hill) {
		randomVarInHill(sol.variable(), hill);
	}

	bool HGHEC::isSolInHill(const SolutionType &sol, Hill *hill) const {
		return isVarInHill(sol.variable(), hill);
	}

	const std::list<HGHEC::Hill*>& HGHEC::theHillsLocated(const SolutionType &sol) const {
		return theHillsLocated(sol.variable());
	}

	const HGHEC::SolutionType* HGHEC::archiveSolution(const SolutionType &sol, TaskEval task, Environment *env) {
		auto p_sol = HGHE::archiveSolution(sol, task, env);
		updateInfoSSP(p_sol, task, env);
		return p_sol;
	}

	void HGHEC::updateInfoSSP(const SolutionType *p_sol, TaskEval task, Environment *env) {
		auto id_ssp = m_sp_tree->getRegionIdx(p_sol->variable().vector());
		if (!m_subspaces[id_ssp]->best_sol || dominate(*p_sol, 
			*m_subspaces[id_ssp]->best_sol, env->problem()->optimizeMode())
		) {
			m_subspaces[id_ssp]->best_sol = p_sol;
			m_subspaces[id_ssp]->best_in_exploit = task == TaskEval::kExploit;
		}
		m_subspaces[id_ssp]->his_sols.push_back(p_sol);
		if (task == TaskEval::kExplore) {
			m_subspaces[id_ssp]->his_explore.push_back(p_sol);
		}
		else {
			m_subspaces[id_ssp]->his_exploit.push_back(p_sol);
		}
	}

	void HGHEC::subdivideSpace(const std::vector<const SolutionType*> &seed_sols, Environment *env) {
		if (seed_sols.size() < 2) {
			return;
		}
		size_t num_vars = env->problem()->numberVariables();
		std::vector<std::vector<Real>> max_size(seed_sols.size(), std::vector<Real>(num_vars, -1));
		std::vector<std::vector<std::vector<Real>>> dist_each_dim(num_vars);
		for (size_t j = 0; j < num_vars; ++j) {
			dist_each_dim[j].resize(seed_sols.size(), std::vector<Real>(seed_sols.size()));
		}
		for (size_t j = 0; j < num_vars; ++j) {
			for (size_t i = 0; i < seed_sols.size(); ++i) {
				for (size_t k = i + 1; k < seed_sols.size(); ++k) {
					dist_each_dim[j][i][k] = dist_each_dim[j][k][i] =
						fabs(seed_sols[i]->variable()[j] - seed_sols[k]->variable()[j]);
				}
			}
		}
		for (size_t i = 0; i < seed_sols.size(); ++i) {
			for (size_t k = i + 1; k < seed_sols.size(); ++k) {
				size_t dim_max = 0;
				for (size_t j = 1; j < num_vars; ++j) {
					if (dist_each_dim[j][i][k] > dist_each_dim[dim_max][i][k]) {
						dim_max = j;
					}
				}
				if (max_size[i][dim_max] == -1 || max_size[i][dim_max] > dist_each_dim[dim_max][i][k]) {
					max_size[i][dim_max] = dist_each_dim[dim_max][i][k];
				}
				if (max_size[k][dim_max] == -1 || max_size[k][dim_max] > dist_each_dim[dim_max][i][k]) {
					max_size[k][dim_max] = dist_each_dim[dim_max][i][k];
				}
			}
		}
		for (size_t i = 0; i < max_size.size(); ++i) {
			for (size_t j = 0; j < num_vars; ++j) {
				if (max_size[i][j] != -1) {
					max_size[i][j] = max_size[i][j] / 2;
				}
			}
		}
		for (size_t i = 0; i < seed_sols.size(); ++i) {
			while (true) {
				bool proposition = true;
				auto id_ssp = m_sp_tree->getRegionIdx(seed_sols[i]->variable().vector());
				auto &box = m_sp_tree->getBox(id_ssp);
				for (size_t j = 0; j < num_vars; ++j) {
					if (max_size[i][j] != -1 && (box[j].second - box[j].first) >= max_size[i][j]) {
						proposition = false;
						break;
					}
				}
				if (!proposition) {
					halveSubspace(m_subspaces[id_ssp].get(), env);
				}
				else {
					break;
				}
			};
		}
	}

	void HGHEC::assignSubspaceFitness() {
		size_t num_missing = 0;
		for (auto &subspace : m_subspaces) {
			if (subspace->best_sol != nullptr) {
				subspace->obj = subspace->best_sol->objective(0);
			}
			else {
				num_missing++;
			}
		}
		if (num_missing > 0) {
			bilinearInterpolation();
		}
	}

	HGHEC::Subspace* HGHEC::halveSubspace(Subspace *subspace, Environment *env) {
		auto new_subspace = subspace->halve(m_sp_tree.get());
		m_subspaces.emplace_back(new_subspace);
		/* update historical solutions */
		auto his_exploit = subspace->his_exploit;
		auto his_explore = subspace->his_explore;
		subspace->best_sol = nullptr;
		subspace->his_sols.clear();
		subspace->his_exploit.clear();
		subspace->his_explore.clear();
		new_subspace->best_sol = nullptr;
		new_subspace->his_sols.clear();
		new_subspace->his_exploit.clear();
		new_subspace->his_explore.clear();
		for (auto p_sol : his_exploit) {
			updateInfoSSP(p_sol, TaskEval::kExploit, env);
		}
		for (auto p_sol : his_explore) {
			updateInfoSSP(p_sol, TaskEval::kExplore, env);
		}
		return new_subspace;
	}

	void HGHEC::bilinearInterpolation() {
		std::vector<bool> missed(m_subspaces.size(), true);
		std::vector<int> num_smpled_nbrs(m_subspaces.size(), 0);
		for (size_t i = 0; i < m_subspaces.size(); i++) {
			if (m_subspaces[i]->best_sol != nullptr) {
				m_subspaces[i]->obj = m_subspaces[i]->best_sol->objective(0);
				for (Subspace *adj_subspace : m_subspaces[i]->adjacentSubspaces()) {
					num_smpled_nbrs[adj_subspace->id]++;
				}
				missed[i] = false;
			}
		}
		std::priority_queue<std::pair<int, Subspace*>> heap;
		for (auto &subspace : m_subspaces) {
			if (missed[subspace->id]) {
				heap.emplace(-num_smpled_nbrs[subspace->id], subspace.get());
			}
		}
		std::vector<bool> interpolated(m_subspaces.size(), false);
		while (!heap.empty()) {
			Subspace *subspace = heap.top().second;
			if (interpolated[subspace->id]) {
				heap.pop();
			}
			else {
				subspace->obj = 0;
				int num = 0;
				for (Subspace *adj_subspace : subspace->adjacentSubspaces()) {
					if (adj_subspace->best_sol != nullptr || interpolated[adj_subspace->id]) {
						subspace->obj += adj_subspace->obj;
						num++;
					}
				}
				subspace->obj /= num;
				interpolated[subspace->id] = true;
				heap.pop();
				for (Subspace *adj_subspace : subspace->adjacentSubspaces()) {
					if (missed[adj_subspace->id] && !interpolated[adj_subspace->id]) {
						num_smpled_nbrs[adj_subspace->id]++;
						heap.emplace(-num_smpled_nbrs[adj_subspace->id], adj_subspace);
					}
				}
			}
		}
	}

	bool HGHEC::isTimeToLearn(Hill *hill) {
		int num_sols = 0;
		for (Subspace* subspace : hill->subspaces()) {
			if (m_use_his_explore_only) {
				num_sols += subspace->his_explore.size();
			}
			else {
				num_sols += subspace->his_sols.size();
			}
		}
		int mean_learn = 0;
		for (Subspace *subspace : hill->subspaces()) {
			mean_learn += subspace->num_learn;
		}
		mean_learn /= hill->subspaces().size();
		return num_sols > m_num_sols_thold + m_num_sols_intvl * (mean_learn - 1);
	}

	void HGHEC::getHisSolsInHill(std::vector<const SolutionType*> &his_sols, Hill *hill) {
		his_sols.clear();
		if (m_use_his_explore_only) {
			std::vector<const SolutionType*> best_sols, explr_sols;
			for (Subspace* subspace : hill->subspaces()) {
				auto best_sol = subspace->best_sol;
				if (best_sol) {
					if (subspace->best_in_exploit) {
						best_sols.push_back(best_sol);
					}
					auto &explr_sols_ssp = subspace->his_explore;
					explr_sols.insert(explr_sols.end(), explr_sols_ssp.begin(), explr_sols_ssp.end());
				}
			}
			his_sols.insert(his_sols.end(), best_sols.begin(), best_sols.end());
			his_sols.insert(his_sols.end(), explr_sols.begin(), explr_sols.end());
		}
		else {
			for (Subspace* subspace : hill->subspaces()) {
				his_sols.insert(his_sols.end(), subspace->his_sols.begin(), subspace->his_sols.end());
			}
		}
	}

	void HGHEC::updateNumLearn(size_t num_seed_sols, Hill *hill) {
		for (Subspace *subspace : hill->subspaces()) {
			subspace->num_learn++;
		}
	}

	void HGHEC::identifySeedSolutions(const std::vector<const SolutionType*> &his_sols,
		std::vector<const SolutionType *> &seed_sols, Environment *env
	) {
		NBC nbc(m_phi, NBC::kByKDTree, NBC::kByOutlier);
		nbc.setData(his_sols);
		nbc.clustering(env);
		seed_sols.clear();
		for (size_t id_sol : nbc.clusterCenters()) {
			seed_sols.push_back(his_sols[id_sol]);
		}
	}

	void HGHEC::groupSubspace(Environment *env) {
		std::multimap<Real, size_t> sort_fitness;
		std::vector<std::map<size_t, size_t>> min_num_step_to_cluster(m_subspaces.size());
		std::list<std::set<size_t>> clusters;
		for (size_t i = 0; i < m_subspaces.size(); ++i) {
			sort_fitness.emplace(m_subspaces[i]->obj, i);
		}
		if (env->problem()->optimizeMode()[0] == OptimizeMode::kMinimize) {
			for (auto it = sort_fitness.begin(); it != sort_fitness.end(); ++it) {
				if (min_num_step_to_cluster[it->second].empty()) {
					std::set<size_t> cluster;
					findInferiorNeighbors1(it->second, cluster, min_num_step_to_cluster, clusters.size(), 0);
					clusters.emplace_back(std::move(cluster));
				}
			}
		}
		else {
			for (auto it = sort_fitness.rbegin(); it != sort_fitness.rend(); ++it) {
				if (min_num_step_to_cluster[it->second].empty()) {
					std::set<size_t> cluster;
					findInferiorNeighbors2(it->second, cluster, min_num_step_to_cluster, clusters.size(), 0);
					clusters.emplace_back(std::move(cluster));
				}
			}
		}


		std::vector<Hill*> hills;
		for (size_t k = 0; k < clusters.size(); ++k)
			hills.emplace_back(new Hill);
		if (m_keep_subspace_valley) {
			for (size_t i = 0; i < m_subspaces.size(); ++i) {
				m_subspaces[i]->clearAssignedHills(m_sp_tree.get());
				for (auto &it : min_num_step_to_cluster[i]) {
					m_subspaces[i]->joinHill(hills[it.first], m_sp_tree.get(), true);
				}
			}
		}
		else {
			for (size_t i = 0; i < m_subspaces.size(); ++i) {
				m_subspaces[i]->clearAssignedHills(m_sp_tree.get());
				size_t min_num_step = m_subspaces.size();
				size_t hill_located = m_hills.size();
				for (auto &it : min_num_step_to_cluster[i]) {
					if (min_num_step >= it.second) {
						min_num_step = it.second;
						hill_located = it.first;
					}
				}
				m_subspaces[i]->joinHill(hills[hill_located], m_sp_tree.get(), true);
			}
		}
		m_hills.clear();
		for (size_t k = 0; k < hills.size(); ++k) {
			m_hills.emplace_back(hills[k]);
		}
		m_ptr_to_id.clear();
		m_id_to_ptr.clear();
		for (auto &hill : m_hills) {
			m_ptr_to_id[hill.get()] = m_ptr_to_id.size();
			m_id_to_ptr[m_ptr_to_id[hill.get()]] = hill.get();
		}
	}

	void HGHEC::findInferiorNeighbors1(size_t center, std::set<size_t> &cluster,
		std::vector<std::map<size_t, size_t>> &min_num_step_to_cluster, size_t id_cluster, size_t num_step)
	{
		cluster.emplace(center);
		if (min_num_step_to_cluster[center].count(id_cluster) == 0 || 
			min_num_step_to_cluster[center][id_cluster] > num_step) 
		{
			min_num_step_to_cluster[center][id_cluster] = num_step;
		}
		for (Subspace *adj_subspace : m_subspaces[center]->adjacentSubspaces()) {
			if (adj_subspace->obj >= m_subspaces[center]->obj && cluster.count(adj_subspace->id) == 0) {
				findInferiorNeighbors1(adj_subspace->id, cluster, min_num_step_to_cluster, id_cluster, num_step + 1);
			}
		}
	}

	void HGHEC::findInferiorNeighbors2(size_t center, std::set<size_t> &cluster,
		std::vector<std::map<size_t, size_t>> &min_num_step_to_cluster, size_t id_cluster, size_t num_step)
	{
		cluster.emplace(center);
		if (min_num_step_to_cluster[center].count(id_cluster) == 0 || 
			min_num_step_to_cluster[center][id_cluster] > num_step)
		{
			min_num_step_to_cluster[center][id_cluster] = num_step;
		}
		for (Subspace *adj_subspace : m_subspaces[center]->adjacentSubspaces()) {
			if (adj_subspace->obj <= m_subspaces[center]->obj && cluster.count(adj_subspace->id) == 0) {
				findInferiorNeighbors2(adj_subspace->id, cluster, min_num_step_to_cluster, id_cluster, num_step + 1);
			}
		}
	}

	void HGHEC::Hill::addSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace) {
		if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) != m_subspaces.end()) {
			return;
		}
		m_subspaces.push_back(subspace);
		m_volume += sp_tree->getBoxVolume(subspace->id);
		if (update_subspace) {
			subspace->joinHill(this, sp_tree);
		}
	}

	void HGHEC::Hill::removeSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace) {
		if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) == m_subspaces.end())
			return;
		m_subspaces.remove(subspace);
		m_volume -= sp_tree->getBoxVolume(subspace->id);
		if (update_subspace) {
			subspace->quitHill(this, sp_tree);
		}
	}

	void HGHEC::Hill::merge(Hill *hill, const SPTree *sp_tree, bool update_subspace) {
		for (auto subspace : hill->m_subspaces) {
			addSubspace(subspace, sp_tree, update_subspace);
		}
		hill->clear(sp_tree, update_subspace);
	}

	void HGHEC::Hill::clear(const SPTree *sp_tree, bool update_subspace) {
		if (update_subspace) {
			for (auto subspace : m_subspaces) {
				subspace->quitHill(this, sp_tree);
			}
		}
		m_subspaces.clear();
		m_volume = 0;
	}

	HGHEC::Subspace* HGHEC::Hill::rouletteWheelSelection(Random *rnd, const SPTree *sp_tree) const {
		if (m_subspaces.empty()) {
			return nullptr;
		}
		Real rand_pos = m_volume * rnd->uniform.next();
		Real accum = 0;
		for (auto subspace : m_subspaces) {
			accum += sp_tree->getBoxVolume(subspace->id);
			if (rand_pos <= accum) {
				return subspace;
			}
		}
		return nullptr;
	}

	void HGHEC::Subspace::joinHill(Hill *hill, const SPTree *sp_tree, bool update_hill) {
		if (m_set_assigned_hills.count(hill)) {
			return;
		}
		m_assigned_hills.push_back(hill);
		m_set_assigned_hills.insert(hill);
		if (update_hill) {
			hill->addSubspace(this, sp_tree);
		}
	}

	void HGHEC::Subspace::quitHill(Hill *hill, const SPTree *sp_tree, bool update_hill) {
		if (!m_set_assigned_hills.count(hill)) {
			return;
		}
		m_assigned_hills.remove(hill);
		m_set_assigned_hills.erase(hill);
		if (update_hill) {
			hill->removeSubspace(this, sp_tree);
		}
	}

	void HGHEC::Subspace::clearAssignedHills(const SPTree *sp_tree, bool update_hill) {
		if (update_hill) {
			for (Hill *hill : m_assigned_hills) {
				hill->removeSubspace(this, sp_tree);
			}
		}
		m_assigned_hills.clear();
		m_set_assigned_hills.clear();
	}

	HGHEC::Subspace* HGHEC::Subspace::halve(SPTree *sp_tree) {
		size_t id_new_subspace = sp_tree->splitRegion(id);
		auto new_subspace = new Subspace(id_new_subspace);
		new_subspace->m_assigned_hills = m_assigned_hills;
		new_subspace->m_set_assigned_hills = m_set_assigned_hills;
		new_subspace->m_adjacent_subspaces = m_adjacent_subspaces;
		for (Hill *hill : m_assigned_hills) {
			hill->m_subspaces.push_back(new_subspace);
		}
		for (auto it = m_adjacent_subspaces.begin(); it != m_adjacent_subspaces.end();) {
			if (!sp_tree->checkAdjacency((*it)->id, id)) {
				(*it)->m_adjacent_subspaces.remove(this);
				it = this->m_adjacent_subspaces.erase(it);
			}
			else {
				it++;
			}
		}
		m_adjacent_subspaces.push_back(new_subspace);
		for (auto it = new_subspace->m_adjacent_subspaces.begin(); it != new_subspace->m_adjacent_subspaces.end();) {
			if (sp_tree->checkAdjacency((*it)->id, new_subspace->id)) {
				(*it)->m_adjacent_subspaces.push_back(new_subspace);
				it++;
			}
			else {
				it = new_subspace->m_adjacent_subspaces.erase(it);
			}
		}
		new_subspace->m_adjacent_subspaces.push_back(this);
		return new_subspace;
	}
}