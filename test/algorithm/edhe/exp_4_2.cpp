#include "../../../utility/catch.hpp"
#include "../../../utility/clustering/nbc.h"
#include "../../../instance/problem/continuous/free_peaks/free_peaks.h"
#include "../../../core/global.h"
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>

using namespace ofec;

using SPTree = nanoflann::KDTreeSpace<Real>;

class Subspace;
class Hill {
	friend class Subspace;
protected:
	std::list<Subspace*> m_subspaces;
	Real m_volume = 0;
public:
	void addSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace = false);
	void removeSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace = false);
	void merge(Hill *hill, const SPTree *sp_tree, bool update_subspace = false);
	void clear(const SPTree *sp_tree, bool update_subspace = false);
	Subspace *rouletteWheelSelection(Random *rnd, const SPTree *sp_tree) const;
	const std::list<Subspace *> &subspaces() const { return m_subspaces; }
	Real volume() const { return m_volume; }
};

class Subspace {
protected:
	std::list<Hill *> m_assigned_hills;
	std::list<Subspace *> m_adjacent_subspaces;
public:
	const size_t m_id;
	Subspace(size_t id) : m_id(id) {}
	void joinHill(Hill *hill, const SPTree *sp_tree, bool update_hill = false);
	void quitHill(Hill *hill, const SPTree *sp_tree, bool update_hill = false);
	Subspace *halve(SPTree *sp_tree);
	void clearAssignedHills(const SPTree *sp_tree, bool update_hill = false);
	const std::list<Hill *> &assignedHills() const { return m_assigned_hills; }
	const std::list<Subspace *> &adjacentSubspaces() const { return m_adjacent_subspaces; }
};

void Hill::addSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace) {
	if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) != m_subspaces.end())
		return;
	m_subspaces.push_back(subspace);
	m_volume += sp_tree->getBoxVolume(subspace->m_id);
	if (update_subspace) {
		subspace->joinHill(this, sp_tree);
	}
}

void Hill::removeSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace) {
	if (std::find(m_subspaces.begin(), m_subspaces.end(), subspace) == m_subspaces.end())
		return;
	m_subspaces.remove(subspace);
	m_volume -= sp_tree->getBoxVolume(subspace->m_id);
	if (update_subspace) {
		subspace->quitHill(this, sp_tree);
	}
}

void Hill::merge(Hill *hill, const SPTree *sp_tree, bool update_subspace) {
	if (hill == this)
		return;
	for (auto subspace : hill->m_subspaces) {
		addSubspace(subspace, sp_tree, update_subspace);
	}
	hill->clear(sp_tree, update_subspace);
}

void Hill::clear(const SPTree *sp_tree, bool update_subspace) {
	if (update_subspace) {
		for (auto subspace : m_subspaces) {
			subspace->quitHill(this, sp_tree);
		}
	}
	m_subspaces.clear();
	m_volume = 0;
}

Subspace *Hill::rouletteWheelSelection(Random *rnd, const SPTree *sp_tree) const {
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

void Subspace::joinHill(Hill *hill, const SPTree *sp_tree, bool update_hill) {
	if (std::find(m_assigned_hills.begin(), m_assigned_hills.end(), hill) != m_assigned_hills.end())
		return;
	m_assigned_hills.push_back(hill);
	if (update_hill) {
		hill->addSubspace(this, sp_tree);
	}
}

void Subspace::quitHill(Hill *hill, const SPTree *sp_tree, bool update_hill) {
	if (std::find(m_assigned_hills.begin(), m_assigned_hills.end(), hill) == m_assigned_hills.end())
		return;
	m_assigned_hills.remove(hill);
	if (update_hill) {
		hill->removeSubspace(this, sp_tree);
	}
}

void Subspace::clearAssignedHills(const SPTree *sp_tree, bool update_hill) {
	if (update_hill) {
		for (Hill *hill : m_assigned_hills)
			hill->removeSubspace(this, sp_tree);
	}
	m_assigned_hills.clear();
}

Subspace *Subspace::halve(SPTree *sp_tree) {
	size_t id_new_subspace = sp_tree->splitRegion(m_id);
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
	size_t size1 = new_subspace->m_adjacent_subspaces.size();
	size_t size2 = this->m_adjacent_subspaces.size();
	return new_subspace;
}

Real worstAccuracy(
	const FreePeaks *FPs, Random *rnd,
	const std::vector<const SolutionBase*> &sols_to_group, 
	const std::vector<std::vector<size_t>> &clusters) 
{
	auto sp_tree = std::make_unique<SPTree>();
	sp_tree->setInitBox(FPs->boundary());
	sp_tree->inputRatioData({ 1.0 });
	sp_tree->buildIndex();
	std::list<std::unique_ptr<Hill>> hills;
	std::vector<std::unique_ptr<Subspace>> subspaces;
	subspaces.emplace_back(new Subspace(0));
	hills.emplace_back(new Hill);
	hills.front()->addSubspace(subspaces.front().get(), sp_tree.get(), true);
	auto cur_hill = hills.front().get();

	std::vector<size_t> sol_aff_clu(sols_to_group.size());
	for (size_t k = 0; k < clusters.size(); ++k) {
		for (size_t i : clusters[k]) {
			sol_aff_clu[i] = k;
		}
	}
	std::map<size_t, std::set<size_t>> sols_in_ssp;
	for (size_t i = 0; i < sols_to_group.size(); ++i) {
		auto id_ssp = sp_tree->getRegionIdx(dynamic_cast<const Solution<>*>(sols_to_group[i])->variable().vect());
		sols_in_ssp[id_ssp].insert(i);
	}
	std::map<size_t, int> ssp_aff_clu;
	auto ssps_aff_clu_uknown = cur_hill->subspaces();
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
				for (auto &p : num_clu) {
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
						if (sp_tree->getRegionIdx(dynamic_cast<const Solution<>*>(sols_to_group[*it])->variable().vect()) != subspace->m_id) {
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
	std::map<size_t, Hill *> clu_to_hill;
	clu_to_hill[0] = cur_hill;
	for (size_t k = 1; k < clusters.size(); ++k) {
		auto new_hill = new Hill;
		clu_to_hill[k] = new_hill;
		hills.emplace_back(new_hill);
	}
	std::list<Subspace *> subspaces_to_group;
	auto hill_kept_subspaces = cur_hill->subspaces();
	for (Subspace *subspace : hill_kept_subspaces) {
		subspace->clearAssignedHills(sp_tree.get(), true);
		subspaces_to_group.push_back(subspace);
	}
	std::set<size_t> empty_subspace_ids;
	for (Subspace *subspace : subspaces_to_group) {
		if (ssp_aff_clu[subspace->m_id] > -1)
			subspace->joinHill(clu_to_hill[ssp_aff_clu[subspace->m_id]], sp_tree.get(), true);
		else
			empty_subspace_ids.insert(subspace->m_id);
	}
	while (!empty_subspace_ids.empty()) {
		std::map<size_t, std::list<Hill *>> border;
		for (size_t id_subspace : empty_subspace_ids) {
			for (Subspace *adj_subspace : subspaces[id_subspace]->adjacentSubspaces()) {
				for (Hill *hill : adj_subspace->assignedHills()) {
					if (std::find(border[id_subspace].begin(), border[id_subspace].end(), hill) == border[id_subspace].end())
						border[id_subspace].push_back(hill);
				}
			}
		}
		for (auto &p : border) {
			for (Hill *hill : p.second)
				subspaces[p.first]->joinHill(hill, sp_tree.get(), true);
			empty_subspace_ids.erase(p.first);
		}
	}
	Real min_ratio = 1.0;
	std::vector<Real> tmp_var(FPs->numberVariables());
	for (size_t k = 0; k < FPs->numberOptima(); ++k) {
		auto hill = subspaces[sp_tree->getRegionIdx(FPs->optima()->variable(k).vect())]->assignedHills().front();
		Real ratio = 0;
		for (size_t i = 0; i < 1000; ++i) {
			size_t id_ssp = hill->rouletteWheelSelection(rnd, sp_tree.get())->m_id;
			auto &box = sp_tree->getBox(id_ssp);
			for (size_t j = 0; j < box.size(); j++)
				tmp_var[j] = rnd->uniform.nextNonStd(box[j].first, box[j].second);
			if (FPs->subspaceTree().tree->getRegionIdx(tmp_var) == k)
				ratio++;
		}
		ratio /= 1000;
		if (ratio < min_ratio)
			min_ratio = ratio;
	}
	return min_ratio;
}

Real worstAccuracy2(
	const FreePeaks *FPs,
	const std::vector<const SolutionBase *> &sols_to_group,
	const std::vector<std::vector<size_t>> &clusters)
{
	Real min_ratio = 1.0;
	for (size_t k = 0; k < clusters.size(); ++k) {
		Real ratio = 0;
		size_t id_peak = FPs->subspaceTree().tree->getRegionIdx(
			dynamic_cast<const VariableVector<>&>(sols_to_group[clusters[k][0]]->variableBase()).vect());
		for (size_t i : clusters[k]) {
			if (id_peak == FPs->subspaceTree().tree->getRegionIdx(
				dynamic_cast<const VariableVector<>&>(sols_to_group[i]->variableBase()).vect()))
				ratio++;
		}
		ratio /= clusters[k].size();
		if (ratio < min_ratio)
			min_ratio = ratio;
	}
	return min_ratio;
}

void numSamplesOnIdentifyBndBoAs(
	size_t times,
	size_t num_vars,
	size_t min, size_t max, size_t step,
	size_t num_runs)
{
	std::cout << "Times=" << times << ", NumVars=" << num_vars << std::endl;

	ParameterMap param;
	param["problem name"] = std::string("free_peaks");
	param["number of variables"] = int(num_vars);
	param["generation type"] = std::string("read_file");
	param["dataFile1"] = "sop/2_s5_1_" + std::to_string(times);

	auto cparam = std::make_shared<const ParameterMap>(param);
	auto pro = Problem::generateByFactory(cparam, 0.5);
	auto FPs = dynamic_cast<FreePeaks*>(pro.get());
	FPs->initialize();

	std::list<std::unique_ptr<SolutionBase>> opts;
	std::vector<const SolutionBase*> ptr_opts;
	for (size_t k = 0; k < FPs->numberOptima(); ++k) {
		auto new_sol = FPs->createSolution(FPs->optima()->variable(k));
		new_sol->objective() = FPs->optima()->objective(k);
		opts.emplace_back(new_sol);
		ptr_opts.push_back(opts.back().get());
	}
	std::vector<size_t> centers(opts.size());
	std::iota(centers.begin(), centers.end(), 0);
	std::vector<std::vector<Real>> records(num_runs);
	Random rnd(0.5);
	NBC nbc;
	for (size_t id_run = 0; id_run < num_runs; ++id_run) {
		std::cout << "\tRunID=" << id_run << std::endl;
		std::list<std::unique_ptr<SolutionBase>> sols_to_group;
		std::vector<const SolutionBase*> ptr_sols_to_group(ptr_opts);
		std::vector<size_t> centers_(centers);
		for (size_t num_samples = min; num_samples <= max; num_samples += step) {
			while (sols_to_group.size() < num_samples) {
				auto new_sol = FPs->createSolution();
				new_sol->initialize(FPs, &rnd);
				new_sol->evaluate(FPs, nullptr, false);
				sols_to_group.emplace_back(new_sol);
				ptr_sols_to_group.push_back(sols_to_group.back().get());
			}
			nbc.setData(ptr_sols_to_group);
			nbc.updateNbDistByKDTree(FPs);
			nbc.cutEdgesInGraph(centers_);
			nbc.updateClusters();
			records[id_run].push_back(worstAccuracy(FPs, &rnd, ptr_sols_to_group, nbc.clusters()));
		}
	}
	std::stringstream out;
	for (size_t i = 0; i < (max - min) / step; ++i)
		out << "#smps=" << min + i * step << ',';
	out << "#smps=" << min + (max - min) / step * step << '\n';
	for (size_t id_run = 0; id_run < num_runs; ++id_run) {
		for (size_t i = 0; i < records[id_run].size() - 1; ++i)
			out << records[id_run][i] << ',';
		out << records[id_run].back() << '\n';
	}
	std::stringstream file_name;
	//file_name << g_working_directory << "/result/edhe/exp_4_2/2/" << times << '_' << num_vars << ".csv";
	file_name << g_working_directory << "/result/edhe/exp_4_2/" << times << '_' << num_vars << ".csv";
	std::ofstream out_file(file_name.str());
	out_file << out.str();
	out_file.close();
}

TEST_CASE("Impact of # of samples on identifing bnds of BoAs", "[EDHE][sample][boundary]") {
	//numSamplesOnIdentifyBndBoAs(1, 2, 50, 100, 5, 100);
	//numSamplesOnIdentifyBndBoAs(2, 2, 50, 100, 5, 100);
	//numSamplesOnIdentifyBndBoAs(3, 2, 50, 100, 5, 100);
	//numSamplesOnIdentifyBndBoAs(4, 2, 50, 100, 5, 100);

	//numSamplesOnIdentifyBndBoAs(2, 1, 50, 100, 5, 100);
	//numSamplesOnIdentifyBndBoAs(2, 3, 50, 100, 5, 100);
	//numSamplesOnIdentifyBndBoAs(2, 4, 50, 100, 5, 100);

	numSamplesOnIdentifyBndBoAs(1, 1, 4, 40, 4, 100);
	numSamplesOnIdentifyBndBoAs(1, 2, 24, 240, 24, 100);
	numSamplesOnIdentifyBndBoAs(1, 3, 144, 1440, 144, 100);
	numSamplesOnIdentifyBndBoAs(1, 4, 864, 8640, 864, 100);

	numSamplesOnIdentifyBndBoAs(2, 1, 4, 40, 4, 100);
	numSamplesOnIdentifyBndBoAs(2, 2, 48, 480, 48, 100);
	numSamplesOnIdentifyBndBoAs(2, 3, 576, 5760, 576, 100);
	numSamplesOnIdentifyBndBoAs(2, 4, 6912, 69120, 6912, 100);

	numSamplesOnIdentifyBndBoAs(3, 1, 4, 40, 4, 100);
	numSamplesOnIdentifyBndBoAs(3, 2, 72, 720, 72, 100);
	numSamplesOnIdentifyBndBoAs(3, 3, 1296, 12960, 1296, 100);
	//numSamplesOnIdentifyBndBoAs(3, 4, 23328, 233280, 23328, 100);

	numSamplesOnIdentifyBndBoAs(4, 1, 4, 40, 4, 100);
	numSamplesOnIdentifyBndBoAs(4, 2, 96, 960, 96, 100);
	numSamplesOnIdentifyBndBoAs(4, 3, 2304, 23040, 2304, 100);
	//numSamplesOnIdentifyBndBoAs(4, 4, 55296, 552960, 55296, 100);
}
