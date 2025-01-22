#include "../../../utility/catch.hpp"
#include "../../../run/include_problem.h"
#include "../../../core/global.h"
#include "../../../core/problem/solution.h"
#include <fstream>
#include <sstream>
#include "../../../utility/clustering/nbc.h"
#include "../../../utility/kd-tree/kdtree_space.h"
#include "../../../utility/nondominated_sorting/filter_sort.h"
#include <iomanip>
#include "../../../instance/algorithm/multi_modal/de_nrand_1/de_nrand_1_pop.h"
#include "../../../utility/kd-tree/data_adaptor.h"

using namespace ofec;

void testNBD(const std::string &pro_name, int num_vars, const size_t num_sols = 500, const size_t iter = 0) {
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	//generateSamples(pro, 30);

	PopulationDE pop(num_sols, pro);
	pop.initialize(pro, rnd);
	pop.evaluate(pro, -1);
	pop.crossoverRate() = 0.6;
	pop.scalingFactor() = 0.5;
	for (size_t i = 0; i < iter; ++i)
		pop.evolve(pro, -1, rnd);

	std::vector<const SolutionBase *> data;
	for (size_t i = 0; i < num_sols; ++i)
		data.push_back(&pop[i]);

	//{
	//	auto file_name = g_working_directory + "/result/sol_" + pro_name
	//		+ "_iter_" + std::to_string(iter)
	//		+ "_num_" + std::to_string(num_sols)
	//		+ ".dat";
	//	std::ofstream out_file(file_name);
	//	std::stringstream out;
	//	size_t id_task = 1;
	//	for (size_t i = 0; i < num_vars; ++i)
	//		out << std::setw(10) << ('x' + std::to_string(i + 1));
	//	out << std::setw(10) << 'y';
	//	out << '\n';
	//	out << std::fixed << std::setprecision(2);
	//	for (auto &c : data) {
	//		auto s = dynamic_cast<const Solution<>*>(c);
	//		for (size_t j = 0; j < num_vars; ++j)
	//			out << std::setw(10) << s->variable()[j];
	//		out << std::setw(10) << s->objective(0);
	//		out << '\n';
	//	}
	//	out << '\n';
	//	out_file << out.str();
	//	out_file.close();
	//}

	//NBC nbc1(6.0, NBC::UpdateNBD::kByDistMat, NBC::CalThrshld::kByOutlier);
	//nbc1.setData(data, pro);

	//{
	//	auto start = std::chrono::steady_clock::now();
	//	for (size_t i = 0; i < 100; ++i)
	//		nbc1.clustering();
	//	auto end = std::chrono::steady_clock::now();
	//	float diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	//	std::cout << "kByDistMat Time used: " << diff / 1000 << " s" << std::endl;
	//}

	//auto &nbd1 = nbc1.nearestBetterDis();
	//{
	//	std::stringstream out_ss;
	//	for (Real val : nbd1)
	//		out_ss << val << '\n';
	//	auto file_name = g_working_directory + "/result/nbd_" + pro_name
	//		+ "_iter_" + std::to_string(iter)
	//		+ "_num_" + std::to_string(num_sols)
	//		+ ".dat";
	//	std::ofstream out_fs(file_name);
	//	out_fs << out_ss.str();
	//	out_fs.close();
	//}



	NBC nbc2(6.0, NBC::UpdateNBD::kByKDTree, NBC::CalThrshld::kByOutlier);
	nbc2.setData(data, pro);

	{
		auto start = std::chrono::steady_clock::now();
		for (size_t i = 0; i < 100; ++i)
			nbc2.clustering();
		std::cout << "Mean ¡À Stddev: " << nbc2.mean() << " ¡À " << nbc2.stddev() << std::endl;
		auto end = std::chrono::steady_clock::now();
		float diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
		std::cout << "Size of input: " << num_sols << std::endl;
		std::cout << "kByKDTree Time used: " << diff / 1000 / 100 << " ms" << std::endl;
	}


	auto &nbd2 = nbc2.nearestBetterDis();
	{
		std::stringstream out_ss;
		for (Real val : nbd2)
			out_ss << val << '\n';
		auto file_name = g_working_directory + "/result/nbd_" + pro_name
			+ "_iter_" + std::to_string(iter)
			+ "_num_sols_" + std::to_string(num_sols)
			+ "_num_vars_" + std::to_string(num_vars)
			+ "_by_kdtree.dat";
		std::ofstream out_fs(file_name);
		out_fs << out_ss.str();
		out_fs.close();
	}

	DEL_RND(rnd);
	DEL_PRO(pro);
	DEL_PARAM(id_param);
}

void testExploitiveNBD(const std::string &pro_name, int num_vars, const size_t num_sols, const size_t iter) {
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	std::list<Solution<>> his_sols;
	PopulationDE pop(num_sols, pro);
	pop.initialize(pro, rnd);
	pop.evaluate(pro, -1);
	pop.crossoverRate() = 0.6;
	pop.scalingFactor() = 0.5;
	for (size_t i = 0; i < num_sols; ++i)
		his_sols.push_back(pop[i]);
	for (size_t t = 1; t < iter; ++t) {
		pop.evolve(pro, -1, rnd);
		for (size_t i = 0; i < num_sols; ++i)
			his_sols.push_back(pop[i].trial());
	}

	std::vector<const SolutionBase *> data;
	for (auto &s : his_sols)
		data.push_back(&s);

	NBC nbc(6.0, NBC::UpdateNBD::kByKDTree, NBC::CalThrshld::kByOutlier);
	nbc.setData(data, pro);
	nbc.clustering();
	auto &nbd = nbc.nearestBetterDis();
	{
		std::stringstream out_ss;
		for (Real val : nbd)
			out_ss << val << '\n';
		auto file_name = g_working_directory + "/result/exploit_his_nbd_" + pro_name
			+ "_iter_" + std::to_string(iter)
			+ "_num_" + std::to_string(num_sols)
			+ "_by_kdtree.dat";
		std::ofstream out_fs(file_name);
		out_fs << out_ss.str();
		out_fs.close();
	}
	{
		auto file_name = g_working_directory + "/result/exploit_his_nbd_" + pro_name
			+ "_iter_" + std::to_string(iter)
			+ "_num_" + std::to_string(num_sols)
			+ "_by_kdtree_mean_and_stddev.dat";
		std::ofstream out_fs(file_name);
		out_fs << std::setw(10) << "mean" << std::setw(10) << "stddev" << '\n';
		out_fs << std::fixed << std::setprecision(2);
		out_fs << std::setw(10) << nbc.mean() << std::setw(10) << nbc.stddev() << '\n';
		out_fs.close();
	}

	DEL_RND(rnd);
	DEL_PRO(pro);
	DEL_PARAM(id_param);
}

void testExplorativeNBD(const std::string &pro_name, int num_vars, const size_t num_sols, const size_t iter) {
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	std::list<Solution<>> his_sols;
	std::vector<const SolutionBase *> data;
	PopulationDE pop(num_sols, pro);
	for (size_t t = 0; t < iter; ++t) {
		pop.initialize(pro, rnd);
		pop.evaluate(pro, -1);
		for (size_t i = 0; i < num_sols; ++i) {
			his_sols.push_back(pop[i]);
			data.push_back(&his_sols.back());
		}
	}

	NBC nbc(6.0, NBC::UpdateNBD::kByKDTree, NBC::CalThrshld::kByOutlier);
	nbc.setData(data, pro);
	nbc.clustering();
	auto &nbd = nbc.nearestBetterDis();
	{
		std::stringstream out_ss;
		for (Real val : nbd)
			out_ss << val << '\n';
		auto file_name = g_working_directory + "/result/explore_his_nbd_" + pro_name
			+ "_iter_" + std::to_string(iter)
			+ "_num_" + std::to_string(num_sols)
			+ "_by_kdtree.dat";
		std::ofstream out_fs;
		out_fs.open(file_name);
		out_fs << out_ss.str();
		out_fs.close();
	}
	{
		auto file_name = g_working_directory + "/result/explore_his_nbd_" + pro_name
			+ "_iter_" + std::to_string(iter)
			+ "_num_" + std::to_string(num_sols)
			+ "_by_kdtree_mean_and_stddev.dat";
		std::ofstream out_fs(file_name);
		out_fs << std::setw(10) << "mean" << std::setw(10) << "stddev" << '\n';
		out_fs << std::fixed << std::setprecision(2);
		out_fs << std::setw(10) << nbc.mean() << std::setw(10) << nbc.stddev() << '\n';
		out_fs.close();
	}

	DEL_RND(rnd);
	DEL_PRO(pro);
	DEL_PARAM(id_param);
}

void testMixedNBD(const std::string &pro_name, int num_vars, const size_t num_sols, const size_t iter1, const size_t iter2) {
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	std::list<Solution<>> his_sols;
	std::vector<const SolutionBase *> data;
	PopulationDE pop(num_sols, pro);
	for (size_t t = 0; t < iter2; ++t) {
		pop.initialize(pro, rnd);
		pop.evaluate(pro, -1);
		for (size_t i = 0; i < num_sols; ++i) {
			his_sols.push_back(pop[i]);
			data.push_back(&his_sols.back());
		}
	}
	for (size_t t = 0; t < iter1; ++t) {
		pop.evolve(pro, -1, rnd);
		for (size_t i = 0; i < num_sols; ++i) {
			his_sols.push_back(pop[i].trial());
			data.push_back(&his_sols.back());
		}
	}

	NBC nbc(6.0, NBC::UpdateNBD::kByKDTree, NBC::CalThrshld::kByOutlier);
	nbc.setData(data, pro);
	nbc.clustering();
	auto &nbd = nbc.nearestBetterDis();
	{
		std::stringstream out_ss;
		for (Real val : nbd)
			out_ss << val << '\n';
		auto file_name = g_working_directory + "/result/mixed_his_nbd_" + pro_name
			+ "_iter1_" + std::to_string(iter1)
			+ "_iter2_" + std::to_string(iter2)
			+ "_num_" + std::to_string(num_sols)
			+ "_by_kdtree.dat";
		std::ofstream out_fs;
		out_fs.open(file_name);
		out_fs << out_ss.str();
		out_fs.close();
	}
	{
		auto file_name = g_working_directory + "/result/mixed_his_nbd_" + pro_name
			+ "_iter1_" + std::to_string(iter1)
			+ "_iter2_" + std::to_string(iter2)
			+ "_num_" + std::to_string(num_sols)
			+ "_by_kdtree_mean_and_stddev.dat";
		std::ofstream out_fs(file_name);
		out_fs << std::setw(10) << "mean" << std::setw(10) << "stddev" << '\n';
		out_fs << std::fixed << std::setprecision(2);
		out_fs << std::setw(10) << nbc.mean() << std::setw(10) << nbc.stddev() << '\n';
		out_fs.close();
	}

	DEL_RND(rnd);
	DEL_PRO(pro);
	DEL_PARAM(id_param);
}

void testTholdNBC(const std::string &pro_name, int num_vars) {
	NBC nbc;

	const size_t num_sols = 1000, num_runs = 100;

	std::list<Solution<>> sols;
	std::vector<const SolutionBase *> data;
	for (size_t i = 0; i < num_sols; ++i) {
		sols.emplace_back(1, 0, 2);
		data.push_back(&sols.back());
	}

	auto file_name = g_working_directory + "/result/boa/" + pro_name + "(" + std::to_string(num_vars) + "D).dat";
	std::ifstream in(file_name);
	size_t num_div, num_bsns, num_ssp;
	std::vector<Real> obj;
	std::vector<std::list<size_t>> aff_bsn;
	std::string line, str_id_bsn;
	std::stringstream ss;
	size_t id_ssp;
	while (std::getline(in, line)) {
		if (line.find("number of divisions for each variable:") != std::string::npos)
			num_div = atoi(line.substr(line.find(":") + 1).c_str());
		if (line.find("total number of subspaces:") != std::string::npos) {
			num_ssp = atoi(line.substr(line.find(":") + 1).c_str());
			obj.resize(num_ssp);
			aff_bsn.resize(num_ssp);
		}
		if (line.find("total number of basins of attraction:") != std::string::npos)
			num_bsns = atoi(line.substr(line.find(":") + 1).c_str());
		if (line.find("SSp obj BoA") != std::string::npos) {
			while (std::getline(in, line)) {
				std::istringstream ss(line);
				ss >> id_ssp; ss >> obj[id_ssp];
				while (std::getline(ss, str_id_bsn, ','))
					aff_bsn[id_ssp].push_back(atoi(str_id_bsn.c_str()));
			}
		}
	}

	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);
	Real gap_x1, gap_x2;
	auto &range_x1 = CAST_CONOP(pro)->range(0);
	auto &range_x2 = CAST_CONOP(pro)->range(1);
	gap_x1 = (range_x1.second - range_x1.first) / num_div;
	gap_x2 = (range_x2.second - range_x2.first) / num_div;
	size_t row, col;

	std::vector<Real> fixed_thold_num_centers_in_bsn(num_bsns, 0);
	std::vector<Real> random_thold_num_centers_in_bsn(num_bsns, 0);
	std::vector<Real> outlier_thold_num_centers_in_bsn(num_bsns, 0);

	for (size_t i = 0; i < num_runs; ++i) {
		std::cout << "Current number of running: " << i + 1 << std::endl;
		for (auto &s : sols) {
			s.initialize(pro, rnd);
			s.evaluate(pro, -1, false);
		}
		nbc.setData(data, pro);

		nbc.setCalThrshld(NBC::CalThrshld::kByMean);
		nbc.setPhi(2.0);
		nbc.clustering();
		for (size_t c : nbc.clusterCenters()) {
			auto sol_c = dynamic_cast<const Solution<>*>(data[c]);
			row = floor((sol_c->variable()[0] - range_x1.first) / gap_x1);
			col = floor((sol_c->variable()[1] - range_x2.first) / gap_x2);
			id_ssp = row * num_div + col;
			for (size_t b : aff_bsn[id_ssp])
				fixed_thold_num_centers_in_bsn[b]++;
		}

		nbc.setCalThrshld(NBC::CalThrshld::kByMean);
		nbc.setPhi(rnd->uniform.nextNonStd(1.25, 2.25));
		nbc.clustering();
		for (size_t c : nbc.clusterCenters()) {
			auto sol_c = dynamic_cast<const Solution<>*>(data[c]);
			row = floor((sol_c->variable()[0] - range_x1.first) / gap_x1);
			col = floor((sol_c->variable()[1] - range_x2.first) / gap_x2);
			id_ssp = row * num_div + col;
			for (size_t b : aff_bsn[id_ssp])
				random_thold_num_centers_in_bsn[b]++;
		}

		nbc.setCalThrshld(NBC::CalThrshld::kByOutlier);
		nbc.setPhi(6.0);
		nbc.clustering();
		for (size_t c : nbc.clusterCenters()) {
			auto sol_c = dynamic_cast<const Solution<>*>(data[c]);
			row = floor((sol_c->variable()[0] - range_x1.first) / gap_x1);
			col = floor((sol_c->variable()[1] - range_x2.first) / gap_x2);
			id_ssp = row * num_div + col;
			for (size_t b : aff_bsn[id_ssp])
				outlier_thold_num_centers_in_bsn[b]++;
		}
	}

	for (size_t b = 0; b < num_bsns; ++b) {
		fixed_thold_num_centers_in_bsn[b] /= num_runs;
		random_thold_num_centers_in_bsn[b] /= num_runs;
		outlier_thold_num_centers_in_bsn[b] /= num_runs;
	}

	std::ofstream out_fstream;
	size_t id_bsn;

	out_fstream.open(g_working_directory + "/result/boa/" + pro_name + "(" + std::to_string(num_vars) + "D)_fixed.dat");
	out_fstream << "bsn num_center\n";
	id_bsn = 0;
	for (auto num : fixed_thold_num_centers_in_bsn)
		out_fstream << id_bsn++ << " " << num << '\n';
	out_fstream.close();

	out_fstream.open(g_working_directory + "/result/boa/" + pro_name + "(" + std::to_string(num_vars) + "D)_random.dat");
	out_fstream << "bsn num_center\n";
	id_bsn = 0;
	for (auto num : random_thold_num_centers_in_bsn)
		out_fstream << id_bsn++ << " " << num << '\n';
	out_fstream.close();

	out_fstream.open(g_working_directory + "/result/boa/" + pro_name + "(" + std::to_string(num_vars) + "D)_outlier.dat");
	out_fstream << "bsn num_center\n";
	id_bsn = 0;
	for (auto num : outlier_thold_num_centers_in_bsn)
		out_fstream << id_bsn++ << " " << num << '\n';
	out_fstream.close();

	DEL_PARAM(id_param);
	DEL_PRO(pro);
	DEL_RND(rnd);
}

void observeND(size_t D, size_t n) {
	std::vector<std::vector<Real>> data(n, std::vector<Real>(D));
	Random *rnd = ADD_RND(0.5);
	const size_t num_runs = 1000;
	std::vector<Real> NDs;
	for (size_t k = 0; k < num_runs; ++k) {
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < D; j++) {
				data[i][j] = rnd->uniform.next();
			}
		}
		KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<Real>>, Real> mat_index(D, data);
		std::vector<size_t> ret_indexes(2);
		std::vector<Real> out_dists_sqr(2);
		nanoflann::KNNResultSet<Real> result_set(2);
		for (size_t i = 0; i < n; ++i) {
			result_set.init(&ret_indexes[0], &out_dists_sqr[0]);
			mat_index.index->findNeighbors(result_set, &data[i][0], nanoflann::SearchParams(10));
			NDs.push_back(sqrt(out_dists_sqr[1]));
		}
	}
	Real mean_ND, stddev_ND;
	calMeanAndStd(NDs, mean_ND, stddev_ND);
	std::cout << std::fixed << std::setprecision(5);
	std::cout << std::setw(3) << "D" << std::setw(5) << "n" << std::setw(10) << "mean" << std::setw(10) << "stddev" << std::endl;
	std::cout << std::setw(3) << D << std::setw(5) << n << std::setw(10) << mean_ND << std::setw(10) << stddev_ND << std::endl;

	DEL_RND(rnd);
}

void observeND1Norm(size_t D, size_t n) {
	std::vector<std::vector<Real>> data(n, std::vector<Real>(D));
	Random *rnd = ADD_RND(0.5);
	const size_t num_runs = 1000 * 2 / n;
	std::vector<Real> NDs;
	for (size_t k = 0; k < num_runs; ++k) {
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < D; j++) {
				data[i][j] = rnd->uniform.next();
			}
		}
		KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<Real>>, Real, -1, nanoflann::metric_L1> mat_index(D, data);
		std::vector<size_t> ret_indexes(2);
		std::vector<Real> out_dists(2);
		nanoflann::KNNResultSet<Real> result_set(2);
		for (size_t i = 0; i < n; ++i) {
			result_set.init(&ret_indexes[0], &out_dists[0]);
			mat_index.index->findNeighbors(result_set, &data[i][0], nanoflann::SearchParams(10));
			NDs.push_back(out_dists[1]);
		}
	}
	Real mean_ND, stddev_ND;
	calMeanAndStd(NDs, mean_ND, stddev_ND);
	std::cout << std::fixed << std::setprecision(5);
	std::cout << std::setw(3) << "D" << std::setw(5) << "n" << std::setw(10) << "mean" << std::setw(10) << "stddev" << std::endl;
	std::cout << std::setw(3) << D << std::setw(5) << n << std::setw(10) << mean_ND << std::setw(10) << stddev_ND << std::endl;

	DEL_RND(rnd);
}

void observeNDpNorm(size_t D, size_t n, Real p) {
	std::vector<std::vector<Real>> data(n, std::vector<Real>(D));
	Random *rnd = ADD_RND(0.5);
	const size_t num_runs = 1000 * 2 / n;
	std::vector<Real> NDs;
	std::vector<std::vector<Real>> mat_dis(n, std::vector<Real>(n, 0));
	for (size_t id_run = 0; id_run < num_runs; ++id_run) {
		for (size_t i = 0; i < n; ++i) {
			for (size_t j = 0; j < D; j++) {
				data[i][j] = rnd->uniform.next();
			}
		}
		for (size_t i = 0; i < n; ++i) {
			for (size_t k = i + 1; k < n; ++k) {
				mat_dis[i][k] = mat_dis[k][i] = pNormDistance(data[i].begin(), data[i].end(), data[k].begin(), p);
			}
		}
		for (size_t i = 0; i < n; ++i) {
			Real min_dis = pow(D, 1 / p);
			for (size_t k = 0; k < n; ++k) {
				if (i != k && min_dis > mat_dis[i][k])
					min_dis = mat_dis[i][k];
			}
			NDs.push_back(min_dis);
		}
	}
	Real mean_ND, stddev_ND;
	calMeanAndStd(NDs, mean_ND, stddev_ND);
	std::cout << std::fixed << std::setprecision(5);
	//std::cout << std::setw(3) << "D" << std::setw(5) << "n" << std::setw(10) << "mean" << std::setw(10) << "stddev" << std::endl;
	std::cout << std::setw(3) << D << std::setw(5) << n << std::setw(10) << mean_ND << std::setw(10) << stddev_ND << std::endl;

	DEL_RND(rnd);
}

void observeEachDimND(size_t D, size_t n) {
	std::vector<std::vector<Real>> data(n, std::vector<Real>(D));
	Random *rnd = ADD_RND(0.5);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < D; j++) {
			data[i][j] = rnd->uniform.next();
		}
	}
	std::cout << std::fixed << std::setprecision(5);
	std::cout << std::setw(3) << "D" << std::setw(5) << "n" << std::setw(3) << "j" << std::setw(10) << "mean" << std::setw(10) << "stddev" << std::endl;
	for (size_t j = 0; j < D; ++j) {
		std::vector<Real> NDs;
		std::vector<std::vector<Real>> data_j(n, std::vector<Real>(1));
		for (size_t i = 0; i < n; ++i)
			data_j[i][0] = data[i][j];
		KDTreeVectorOfVectorsAdaptor<std::vector<std::vector<Real>>, Real> mat_index(1, data_j);
		std::vector<size_t> ret_indexes(2);
		std::vector<Real> out_dists_sqr(2);
		nanoflann::KNNResultSet<Real> result_set(2);
		for (size_t i = 0; i < n; ++i) {
			result_set.init(&ret_indexes[0], &out_dists_sqr[0]);
			mat_index.index->findNeighbors(result_set, &data[i][0], nanoflann::SearchParams(10));
			NDs.push_back(sqrt(out_dists_sqr[1]));
		}
		Real mean_ND, stddev_ND;
		calMeanAndStd(NDs, mean_ND, stddev_ND);
		std::cout << std::setw(3) << D << std::setw(5) << n << std::setw(3) << j + 1 << std::setw(10) << mean_ND << std::setw(10) << stddev_ND << std::endl;
	}
	DEL_RND(rnd);
}

void observeNBD(size_t D, size_t n) {
	ParameterMap v;
	v["problem name"] = std::string("Classic_normalized_sphere");
	v["number of variables"] = (int)D;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	PopulationDE pop(n, pro);
	NBC nbc(6.0, NBC::UpdateNBD::kByKDTree, NBC::CalThrshld::kByOutlier);

	const size_t num_runs = 1000;
	std::vector<Real> NBDs;
	for (size_t k = 0; k < num_runs; ++k) {
		std::list<Solution<>> his_sols;
		std::vector<const SolutionBase *> data;
		pop.initialize(pro, rnd);
		pop.evaluate(pro, -1);
		for (size_t i = 0; i < n; ++i) {
			his_sols.push_back(pop[i]);
			data.push_back(&his_sols.back());
		}
		nbc.setData(data, pro);
		nbc.clustering();
		auto &nbd = nbc.nearestBetterDis();
		NBDs.insert(NBDs.end(), nbd.begin(), nbd.end());
	}
	Real mean_NBD, stddev_NBD;
	calMeanAndStd(NBDs, mean_NBD, stddev_NBD);

	std::cout << std::fixed << std::setprecision(5);
	std::cout << std::setw(3) << "D" << std::setw(5) << "n" << std::setw(10) << "mean" << std::setw(10) << "stddev" << std::endl;
	std::cout << std::setw(3) << D << std::setw(5) << n << std::setw(10) << mean_NBD << std::setw(10) << stddev_NBD << std::endl;

	DEL_RND(rnd);
	DEL_PRO(pro);
	DEL_PARAM(id_param);
}

void observeNBDPop(size_t D, size_t n, size_t max_iter) {
	ParameterMap v;
	v["problem name"] = std::string("Classic_normalized_sphere");
	v["number of variables"] = (int)D;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	PopDE_nrand_1 pop(n, pro);
	NBC nbc(6.0, NBC::UpdateNBD::kByKDTree, NBC::CalThrshld::kByOutlier);

	const size_t num_runs = 100;

	std::stringstream out_sstream;
	out_sstream << std::setw(3) << "D" << std::setw(5) << "iter" << std::setw(10) << "mean" << std::setw(10) << "stddev" << std::endl;

	std::vector<std::vector<Real>> NBDs(max_iter + 1);
	for (size_t k = 0; k < num_runs; ++k) {
		pop.initialize(pro, rnd);
		pop.evaluate(pro, -1);
		size_t t = 0;
		while (t <= max_iter) {
			std::list<Solution<>> his_sols;
			std::vector<const SolutionBase*> data;
			for (size_t i = 0; i < n; ++i) {
				his_sols.push_back(pop[i]);
				data.push_back(&his_sols.back());
			}
			nbc.setData(data, pro);
			nbc.clustering();
			auto &nbd = nbc.nearestBetterDis();
			NBDs[t].insert(NBDs[t].end(), nbd.begin(), nbd.end());
			if (t == max_iter)
				break;
			pop.evolve(pro, -1, rnd);
			t++;
		}
	}
	out_sstream << std::fixed << std::setprecision(5);
	for (size_t t = 0; t <= max_iter; ++t) {
		Real mean_NBD, stddev_NBD;
		calMeanAndStd(NBDs[t], mean_NBD, stddev_NBD);
		out_sstream << std::setw(3) << D << std::setw(5) << t << std::setw(10) << mean_NBD << std::setw(10) << stddev_NBD << std::endl;
	}

	std::ofstream out_fstream;
	out_fstream.open(g_working_directory + "/result/nbd_pop_dim_" + std::to_string(D) + ".dat");
	out_fstream << out_sstream.str();
	out_fstream.close();

	DEL_RND(rnd);
	DEL_PRO(pro);
	DEL_PARAM(id_param);
}

TEST_CASE("Threshold NBC", "[phi][NBC]") {
	testTholdNBC("Classic_sphere", 2);
	testTholdNBC("Classic_six_hump_camel_back", 2);
	testTholdNBC("Classic_Rastrigin_modified", 2);
	testTholdNBC("Classic_Himmenblau", 2);
}

TEST_CASE("Calculate NBD", "[NBD][NBC]") {
	//testNBD("Classic_sphere", 2);
	//testNBD("Classic_valleys", 2);

	//testNBD("Classic_valleys", 2, 100);
	//testNBD("Classic_valleys", 2, 200);
	//testNBD("Classic_valleys", 2, 300);
	//testNBD("Classic_valleys", 2, 400);
	//testNBD("Classic_valleys", 2, 500);
	//testNBD("Classic_valleys", 2, 600);
	//testNBD("Classic_valleys", 2, 700);
	//testNBD("Classic_valleys", 2, 800);
	//testNBD("Classic_valleys", 2, 900);
	//testNBD("Classic_valleys", 2, 1000);
	//testNBD("Classic_valleys", 2, 1100);
	//testNBD("Classic_valleys", 2, 1200);
	//testNBD("Classic_valleys", 2, 1300);
	//testNBD("Classic_valleys", 2, 1400);
	//testNBD("Classic_valleys", 2, 1500);
	//testNBD("Classic_valleys", 2, 1600);
	//testNBD("Classic_valleys", 2, 1700);
	//testNBD("Classic_valleys", 2, 1800);
	//testNBD("Classic_valleys", 2, 1900);
	//testNBD("Classic_valleys", 2, 2000);

	//testExploitiveNBD("Classic_valleys", 2, 10, 4);
	//testExploitiveNBD("Classic_valleys", 2, 10, 8);
	//testExploitiveNBD("Classic_valleys", 2, 10, 12);

	//testExplorativeNBD("Classic_valleys", 2, 10, 4);
	//testExplorativeNBD("Classic_valleys", 2, 10, 8);
	//testExplorativeNBD("Classic_valleys", 2, 10, 12);
	//

	//testMixedNBD("Classic_valleys", 2, 10, 3, 1);
	//testMixedNBD("Classic_valleys", 2, 10, 6, 2);
	//testMixedNBD("Classic_valleys", 2, 10, 9, 3);

	testNBD("Classic_normalized_sphere", 20, 100);
	testNBD("Classic_normalized_sphere", 20, 1000);
}

TEST_CASE("Observe ND", "[ND]") {
	observeND(1, 2);
	observeND(1, 5);
	observeND(1, 10);
	observeND(1, 20);
	observeND(1, 50);
	observeND(1, 100);
	observeND(1, 200);
	observeND(1, 500);
	observeND(1, 1000);
	std::cout << std::endl;

	observeND(2, 2);
	observeND(2, 5);
	observeND(2, 10);
	observeND(2, 20);
	observeND(2, 50);
	observeND(2, 100);
	observeND(2, 200);
	observeND(2, 500);
	observeND(2, 1000);
	std::cout << std::endl;

	observeND(3, 2);
	observeND(3, 5);
	observeND(3, 10);
	observeND(3, 20);
	observeND(3, 50);
	observeND(3, 100);
	observeND(3, 200);
	observeND(3, 500);
	observeND(3, 1000);
	std::cout << std::endl;

	observeND(5, 2);
	observeND(5, 5);
	observeND(5, 10);
	observeND(5, 20);
	observeND(5, 50);
	observeND(5, 100);
	observeND(5, 200);
	observeND(5, 500);
	observeND(5, 1000);
	std::cout << std::endl;

	observeND(10, 2);
	observeND(10, 5);
	observeND(10, 10);
	observeND(10, 20);
	observeND(10, 50);
	observeND(10, 100);
	observeND(10, 200);
	observeND(10, 500);
	observeND(10, 1000);
	std::cout << std::endl;

	observeND(20, 2);
	observeND(20, 5);
	observeND(20, 10);
	observeND(20, 20);
	observeND(20, 50);
	observeND(20, 100);
	observeND(20, 200);
	observeND(20, 500);
	observeND(20, 1000);
}

TEST_CASE("Observe ND 1-norm", "[ND]") {
	observeND1Norm(1, 2);
	observeND1Norm(1, 5);
	observeND1Norm(1, 10);
	observeND1Norm(1, 20);
	observeND1Norm(1, 50);
	observeND1Norm(1, 100);
	observeND1Norm(1, 200);
	observeND1Norm(1, 500);
	observeND1Norm(1, 1000);
	std::cout << std::endl;

	observeND1Norm(2, 2);
	observeND1Norm(2, 5);
	observeND1Norm(2, 10);
	observeND1Norm(2, 20);
	observeND1Norm(2, 50);
	observeND1Norm(2, 100);
	observeND1Norm(2, 200);
	observeND1Norm(2, 500);
	observeND1Norm(2, 1000);
	std::cout << std::endl;

	observeND1Norm(3, 2);
	observeND1Norm(3, 5);
	observeND1Norm(3, 10);
	observeND1Norm(3, 20);
	observeND1Norm(3, 50);
	observeND1Norm(3, 100);
	observeND1Norm(3, 200);
	observeND1Norm(3, 500);
	observeND1Norm(3, 1000);
	std::cout << std::endl;

	observeND1Norm(5, 2);
	observeND1Norm(5, 5);
	observeND1Norm(5, 10);
	observeND1Norm(5, 20);
	observeND1Norm(5, 50);
	observeND1Norm(5, 100);
	observeND1Norm(5, 200);
	observeND1Norm(5, 500);
	observeND1Norm(5, 1000);
	std::cout << std::endl;

	observeND1Norm(10, 2);
	observeND1Norm(10, 5);
	observeND1Norm(10, 10);
	observeND1Norm(10, 20);
	observeND1Norm(10, 50);
	observeND1Norm(10, 100);
	observeND1Norm(10, 200);
	observeND1Norm(10, 500);
	observeND1Norm(10, 1000);
	std::cout << std::endl;

	observeND1Norm(20, 2);
	observeND1Norm(20, 5);
	observeND1Norm(20, 10);
	observeND1Norm(20, 20);
	observeND1Norm(20, 50);
	observeND1Norm(20, 100);
	observeND1Norm(20, 200);
	observeND1Norm(20, 500);
	observeND1Norm(20, 1000);
	std::cout << std::endl;
}

TEST_CASE("Observe ND 3-norm", "[ND]") {
	observeNDpNorm(1, 2, 3);
	observeNDpNorm(1, 5, 3);
	observeNDpNorm(1, 10, 3);
	observeNDpNorm(1, 20, 3);
	observeNDpNorm(1, 50, 3);
	observeNDpNorm(1, 100, 3);
	observeNDpNorm(1, 200, 3);
	observeNDpNorm(1, 500, 3);
	observeNDpNorm(1, 1000, 3);
	std::cout << std::endl;

	observeNDpNorm(2, 2, 3);
	observeNDpNorm(2, 5, 3);
	observeNDpNorm(2, 10, 3);
	observeNDpNorm(2, 20, 3);
	observeNDpNorm(2, 50, 3);
	observeNDpNorm(2, 100, 3);
	observeNDpNorm(2, 200, 3);
	observeNDpNorm(2, 500, 3);
	observeNDpNorm(2, 1000, 3);
	std::cout << std::endl;

	observeNDpNorm(3, 2, 3);
	observeNDpNorm(3, 5, 3);
	observeNDpNorm(3, 10, 3);
	observeNDpNorm(3, 20, 3);
	observeNDpNorm(3, 50, 3);
	observeNDpNorm(3, 100, 3);
	observeNDpNorm(3, 200, 3);
	observeNDpNorm(3, 500, 3);
	observeNDpNorm(3, 1000, 3);
	std::cout << std::endl;

	observeNDpNorm(5, 2, 3);
	observeNDpNorm(5, 5, 3);
	observeNDpNorm(5, 10, 3);
	observeNDpNorm(5, 20, 3);
	observeNDpNorm(5, 50, 3);
	observeNDpNorm(5, 100, 3);
	observeNDpNorm(5, 200, 3);
	observeNDpNorm(5, 500, 3);
	observeNDpNorm(5, 1000, 3);
	std::cout << std::endl;

	observeNDpNorm(10, 2, 3);
	observeNDpNorm(10, 5, 3);
	observeNDpNorm(10, 10, 3);
	observeNDpNorm(10, 20, 3);
	observeNDpNorm(10, 50, 3);
	observeNDpNorm(10, 100, 3);
	observeNDpNorm(10, 200, 3);
	observeNDpNorm(10, 500, 3);
	observeNDpNorm(10, 1000, 3);
	std::cout << std::endl;

	observeNDpNorm(20, 2, 3);
	observeNDpNorm(20, 5, 3);
	observeNDpNorm(20, 10, 3);
	observeNDpNorm(20, 20, 3);
	observeNDpNorm(20, 50, 3);
	observeNDpNorm(20, 100, 3);
	observeNDpNorm(20, 200, 3);
	observeNDpNorm(20, 500, 3);
	observeNDpNorm(20, 1000, 3);
}

TEST_CASE("Observe ND 5-norm", "[ND]") {
	observeNDpNorm(1, 2, 5);
	observeNDpNorm(1, 5, 5);
	observeNDpNorm(1, 10, 5);
	observeNDpNorm(1, 20, 5);
	observeNDpNorm(1, 50, 5);
	observeNDpNorm(1, 100, 5);
	observeNDpNorm(1, 200, 5);
	observeNDpNorm(1, 500, 5);
	observeNDpNorm(1, 1000, 5);
	std::cout << std::endl;

	observeNDpNorm(2, 2, 5);
	observeNDpNorm(2, 5, 5);
	observeNDpNorm(2, 10, 5);
	observeNDpNorm(2, 20, 5);
	observeNDpNorm(2, 50, 5);
	observeNDpNorm(2, 100, 5);
	observeNDpNorm(2, 200, 5);
	observeNDpNorm(2, 500, 5);
	observeNDpNorm(2, 1000, 5);
	std::cout << std::endl;

	observeNDpNorm(3, 2, 5);
	observeNDpNorm(3, 5, 5);
	observeNDpNorm(3, 10, 5);
	observeNDpNorm(3, 20, 5);
	observeNDpNorm(3, 50, 5);
	observeNDpNorm(3, 100, 5);
	observeNDpNorm(3, 200, 5);
	observeNDpNorm(3, 500, 5);
	observeNDpNorm(3, 1000, 5);
	std::cout << std::endl;

	observeNDpNorm(5, 2, 5);
	observeNDpNorm(5, 5, 5);
	observeNDpNorm(5, 10, 5);
	observeNDpNorm(5, 20, 5);
	observeNDpNorm(5, 50, 5);
	observeNDpNorm(5, 100, 5);
	observeNDpNorm(5, 200, 5);
	observeNDpNorm(5, 500, 5);
	observeNDpNorm(5, 1000, 5);
	std::cout << std::endl;

	observeNDpNorm(10, 2, 5);
	observeNDpNorm(10, 5, 5);
	observeNDpNorm(10, 10, 5);
	observeNDpNorm(10, 20, 5);
	observeNDpNorm(10, 50, 5);
	observeNDpNorm(10, 100, 5);
	observeNDpNorm(10, 200, 5);
	observeNDpNorm(10, 500, 5);
	observeNDpNorm(10, 1000, 5);
	std::cout << std::endl;

	observeNDpNorm(20, 2, 5);
	observeNDpNorm(20, 5, 5);
	observeNDpNorm(20, 10, 5);
	observeNDpNorm(20, 20, 5);
	observeNDpNorm(20, 50, 5);
	observeNDpNorm(20, 100, 5);
	observeNDpNorm(20, 200, 5);
	observeNDpNorm(20, 500, 5);
	observeNDpNorm(20, 1000, 5);
}

TEST_CASE("Observe Each Dim ND", "[ND]") {
	observeEachDimND(20, 1000);
}

TEST_CASE("Observe NBD", "[NBD]") {
	observeNBD(1, 2);
	observeNBD(1, 5);
	observeNBD(1, 10);
	observeNBD(1, 20);
	observeNBD(1, 50);
	observeNBD(1, 100);
	observeNBD(1, 200);
	observeNBD(1, 500);
	observeNBD(1, 1000);
	std::cout << std::endl;

	observeNBD(2, 2);
	observeNBD(2, 5);
	observeNBD(2, 10);
	observeNBD(2, 20);
	observeNBD(2, 50);
	observeNBD(2, 100);
	observeNBD(2, 200);
	observeNBD(2, 500);
	observeNBD(2, 1000);
	std::cout << std::endl;

	observeNBD(3, 2);
	observeNBD(3, 5);
	observeNBD(3, 10);
	observeNBD(3, 20);
	observeNBD(3, 50);
	observeNBD(3, 100);
	observeNBD(3, 200);
	observeNBD(3, 500);
	observeNBD(3, 1000);
	std::cout << std::endl;

	observeNBD(5, 2);
	observeNBD(5, 5);
	observeNBD(5, 10);
	observeNBD(5, 20);
	observeNBD(5, 50);
	observeNBD(5, 100);
	observeNBD(5, 200);
	observeNBD(5, 500);
	observeNBD(5, 1000);
	std::cout << std::endl;

	observeNBD(10, 2);
	observeNBD(10, 5);
	observeNBD(10, 10);
	observeNBD(10, 20);
	observeNBD(10, 50);
	observeNBD(10, 100);
	observeNBD(10, 200);
	observeNBD(10, 500);
	observeNBD(10, 1000);
	std::cout << std::endl;

	observeNBD(20, 2);
	observeNBD(20, 5);
	observeNBD(20, 10);
	observeNBD(20, 20);
	observeNBD(20, 50);
	observeNBD(20, 100);
	observeNBD(20, 200);
	observeNBD(20, 500);
	observeNBD(20, 1000);
}

TEST_CASE("Observe NBD of population", "[NBD]") {
	//observeNBDPop(1, 20, 2);
	//observeNBDPop(2, 20, 4);
	//observeNBDPop(3, 20, 6);
	//observeNBDPop(5, 20, 10);
	//observeNBDPop(10, 20, 20);
	//observeNBDPop(20, 20, 40);
	observeNBDPop(50, 20, 100);
	observeNBDPop(100, 20, 200);
	observeNBDPop(200, 20, 400);
	//observeNBDPop(500, 20, 1000);
	//observeNBDPop(1000, 20, 2000);
	//observeNBDPop(2000, 20, 4000);
}