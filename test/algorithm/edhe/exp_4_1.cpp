#include "../../../utility/catch.hpp"
#include "../../../instance/algorithm/continuous/single_objective/multi_modal/ring_pso/ring_pso.h"
#include "../../../utility/clustering/nbc.h"

using namespace ofec;

void collectDataOnFreePeaks(
	const std::string &pro_name,
	int num_vars,
	size_t pop_size,
	Real rnd_seed,
	size_t max_num_iters = 40) 
{
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	auto cparam = std::make_shared<const ParameterMap>(v);
	auto pro = Problem::generateByFactory(cparam, 0.5);
	pro->initialize();
	Random rnd(rnd_seed);
	RingSwarm pop(pop_size, RingSwarm::Topology::R2, pro.get());
	pop.weight() = 0.7298;
	pop.accelerator1() = 2.05;
	pop.accelerator2() = 2.05;
	pop.initialize(pro.get(), &rnd);
	pop.initVelocity(pro.get(), &rnd);
	pop.evaluate(pro.get(), nullptr);
	pop.initPbest(pro.get());
	pop.setNeighborhood(&rnd);

	std::vector<std::vector<Real>> nbds;
	std::vector<std::vector<Real>> vars;
	std::vector<Real> objs;
	NBC nbc;
	nbc.setData(pop);
	nbc.updateNbDistByKDTree(pro.get());
	nbds.push_back(nbc.nearestBetterDis2());
	for (size_t i = 0; i < pop_size; ++i) {
		vars.push_back(pop[i].variable().vect());
		objs.push_back(pop[i].objective(0));
	}
	for (size_t iter = 0; iter < max_num_iters; ++iter) {
		pop.evolve(pro.get(), nullptr, &rnd);
		nbc.updateNbDistByKDTree(pro.get());
		nbds.push_back(nbc.nearestBetterDis2());
		for (size_t i = 0; i < pop_size; ++i) {
			vars.push_back(pop[i].variable().vect());
			objs.push_back(pop[i].objective(0));
		}
	}

	std::stringstream sout;
	sout << "iteration NBD ";
	for (size_t j = 0; j < num_vars; ++j)
		sout << 'x' + std::to_string(j + 1) + ' ';
	sout << "f\n";
	for (size_t iter = 0; iter < nbds.size(); ++iter) {
		for (size_t i = 0; i < pop_size; ++i) {
			sout << iter << ' ' << nbds[iter][i] << ' ';
			for (size_t j = 0; j < num_vars; ++j)
				sout << vars[iter * pop_size + i][j] << ' ';
			sout << objs[iter * pop_size + i] << '\n';
		}
	}

	std::stringstream fname;
	fname << g_working_directory << "/result/label-nbd-outlier/exp-4-1/info/PN(";
	fname << pro_name << ")_ND(" << num_vars << ")_PS(" << pop_size << ")_SD(" << rnd_seed << ").dat";
	std::ofstream fout(fname.str());
	fout << sout.str();
}

void screenSeedSolutions(
	const std::string &pro_name,
	int num_vars,
	size_t pop_size,
	Real rnd_seed,
	size_t max_num_iters = 40)
{
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	auto cparam = std::make_shared<const ParameterMap>(v);
	auto pro = Problem::generateByFactory(cparam, 0.5);
	pro->initialize();


	std::vector<std::vector<Solution<>>> seeds_each_iter(max_num_iters + 1);

	std::stringstream fname;
	fname << g_working_directory << "/result/label-nbd-outlier/exp-4-1/visual/PN(";
	fname << pro_name << ")_ND(" << num_vars << ")_PS(" << pop_size << ")_SD(" << rnd_seed << ")_sol_outlier.dat";
	std::ifstream fin(fname.str());
	if (fin) {
		std::string line;
		size_t iter;
		std::getline(fin, line);
		while (std::getline(fin, line)) {
			std::istringstream iss(line);
			Solution<> sol(1, 0, num_vars);
			iss >> iter;
			for (size_t j = 0; j < num_vars; ++j)
				iss >> sol.variable()[j];
			iss >> sol.objective(0);
			seeds_each_iter[iter].push_back(std::move(sol));
		}
	}

	std::vector<Solution<>> seeds;

	for (size_t iter = 0; iter < max_num_iters + 1; ++iter) {
		if (seeds.size() >= seeds_each_iter[iter].size()) {
			for (auto &seed : seeds_each_iter[iter]) {
				Real min_dis = std::numeric_limits<Real>::max();
				size_t nearest = seeds.size();
				for (size_t i = 0; i < seeds.size(); ++i) {
					Real dis = seed.variableDistance(seeds[i], pro.get());
					if (dis < min_dis) {
						min_dis = dis;
						nearest = i;
					}
				}
				if (seed.dominate(seeds[nearest], pro.get()))
					seeds[nearest] = seed;
			}
		}
		else {
			std::vector<bool> mapped(seeds_each_iter[iter].size(), false);
			for (auto &seed : seeds) {
				Real min_dis = std::numeric_limits<Real>::max();
				size_t nearest = seeds_each_iter[iter].size();
				for (size_t i = 0; i < seeds_each_iter[iter].size(); ++i) {
					Real dis = seed.variableDistance(seeds_each_iter[iter][i], pro.get());
					if (dis < min_dis) {
						min_dis = dis;
						nearest = i;
					}
				}
				mapped[nearest] = true;
				if (seeds_each_iter[iter][nearest].dominate(seed, pro.get()))
					seed = seeds_each_iter[iter][nearest];
			}
			for (size_t i = 0; i < seeds_each_iter[iter].size(); ++i) {
				if (!mapped[i])
					seeds.push_back(seeds_each_iter[iter][i]);
			}
		}
	}

	std::stringstream sout;
	for (size_t j = 0; j < num_vars; ++j)
		sout << 'x' + std::to_string(j + 1) + ' ';
	sout << "f\n";
	for (size_t i = 0; i < seeds.size(); ++i) {
		sout << seeds[i].variable() << ' ' << seeds[i].objective(0) << '\n';
	}

	std::stringstream fname2;
	fname2 << g_working_directory << "/result/label-nbd-outlier/exp-4-1/visual/PN(";
	fname2 << pro_name << ")_ND(" << num_vars << ")_PS(" << pop_size << ")_SD(" << rnd_seed << ")_seed_sols.dat";
	std::ofstream fout(fname2.str());
	fout << sout.str();
}

void clusteringBySeedSols(const std::string &pro_name,
	int num_vars,
	size_t pop_size,
	Real rnd_seed,
	size_t max_num_iters = 40)
{
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	auto cparam = std::make_shared<const ParameterMap>(v);
	auto pro = Problem::generateByFactory(cparam, 0.5);
	pro->initialize();

	std::vector<Solution<>> sols;
	std::vector<size_t> centers;

	std::stringstream fname1;
	fname1 << g_working_directory << "/result/label-nbd-outlier/exp-4-1/visual/PN(";
	fname1 << pro_name << ")_ND(" << num_vars << ")_PS(" << pop_size << ")_SD(" << rnd_seed << ")_seed_sols.dat";
	std::ifstream fin1(fname1.str());
	if (fin1) {
		std::string line;
		std::getline(fin1, line);
		while (std::getline(fin1, line)) {
			std::istringstream iss(line);
			Solution<> sol(1, 0, num_vars);
			for (size_t j = 0; j < num_vars; ++j)
				iss >> sol.variable()[j];
			iss >> sol.objective(0);
			sols.push_back(std::move(sol));
			centers.push_back(centers.size());
		}
	}

	std::stringstream fname2;
	fname2 << g_working_directory << "/result/label-nbd-outlier/exp-4-1/info/PN(";
	fname2 << pro_name << ")_ND(" << num_vars << ")_PS(" << pop_size << ")_SD(" << rnd_seed << ").dat";
	std::ifstream fin2(fname2.str());
	if (fin2) {
		std::string line;
		size_t iter;
		Real nbd;
		std::getline(fin2, line);
		while (std::getline(fin2, line)) {
			std::istringstream iss(line);
			Solution<> sol(1, 0, num_vars);
			iss >> iter >> nbd;
			for (size_t j = 0; j < num_vars; ++j)
				iss >> sol.variable()[j];
			iss >> sol.objective(0);
			sols.push_back(std::move(sol));
		}
	}

	NBC nbc;
	nbc.setData(sols);
	nbc.updateNbDistByKDTree(pro.get());
	nbc.cutEdgesInGraph(centers);
	nbc.updateClusters();

	std::stringstream sout;
	for (size_t j = 0; j < num_vars; ++j)
		sout << 'x' + std::to_string(j + 1) + ' ';
	sout << "f k\n";
	for (size_t k = 0; k < nbc.clusters().size(); ++k) {
		for (size_t id_sol : nbc.clusters()[k]) {
			sout << sols[id_sol].variable() << ' ' << sols[id_sol].objective(0) << ' ' << k << '\n';
		}
	}

	std::stringstream fname3;
	fname3 << g_working_directory << "/result/label-nbd-outlier/exp-4-1/visual/PN(";
	fname3 << pro_name << ")_ND(" << num_vars << ")_PS(" << pop_size << ")_SD(" << rnd_seed << ")_clusters.dat";
	std::ofstream fout(fname3.str());
	fout << sout.str();
}

TEST_CASE("Collect sequence of NBDs of evolving population", "[EDHE][data][collect]") {
	collectDataOnFreePeaks("MMOP_CEC2013_F09", 2, 8, 0.5);
	collectDataOnFreePeaks("MMOP_CEC2013_F09", 2, 15, 0.5);
	collectDataOnFreePeaks("MMOP_CEC2013_F09", 2, 20, 0.5);
}

TEST_CASE("Screen seed solutions ", "[EDHE][seed][screen]") {
	screenSeedSolutions("MMOP_CEC2013_F09", 2, 8, 0.5);
	screenSeedSolutions("MMOP_CEC2013_F09", 2, 15, 0.5);
	screenSeedSolutions("MMOP_CEC2013_F09", 2, 20, 0.5);
}

TEST_CASE("Clustering by seed solutions ", "[EDHE][clustering]") {
	clusteringBySeedSols("MMOP_CEC2013_F09", 2, 8, 0.5);
	clusteringBySeedSols("MMOP_CEC2013_F09", 2, 15, 0.5);
	clusteringBySeedSols("MMOP_CEC2013_F09", 2, 20, 0.5);
}
