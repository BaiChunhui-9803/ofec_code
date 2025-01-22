#include "../../../utility/catch.hpp"
#include "../../../core/problem/encoding.h"
#include <set>
#include "../../../instance/algorithm/template/classic/de/population.h"
#include "../../../core/global.h"
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace ofec;

void testExploreOperator(const std::string &pro_name, const size_t num_sols) {
	int num_vars = 2;


	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;

	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	Solution<> tmp_sol(1, 0, num_vars);
	std::list<Solution<>> explore_sols;
	std::vector<const SolutionBase*> explore_data;
	while (explore_sols.size() < num_sols) {
		tmp_sol.initialize(pro, rnd);
		tmp_sol.evaluate(pro, -1, false);
		explore_sols.push_back(tmp_sol);
		explore_data.push_back(&explore_sols.back());
	}

	PopulationDE<> pop(10, pro);
	pop.crossoverRate() = 0.6;
	pop.scalingFactor() = 0.5;
	pop.mutationStrategy() = de::MutateStrategy::kRand1;
	pop.initialize(pro, rnd);
	pop.evaluate(pro, -1);
	std::list<Solution<>> exploit_sols;
	std::vector<const SolutionBase*> exploit_data;
	bool fill = false;
	while (!fill) {
		pop.evolve(pro, -1, rnd);
		for (size_t i = 0; i < pop.size(); ++i) {
			exploit_sols.push_back(pop[i].trial());
			exploit_data.push_back(&exploit_sols.back());
			if (exploit_sols.size() >= num_sols) {
				fill = true;
				break;
			}
		}
	}

	std::vector<const SolutionBase*> data;
	const size_t num_explore = num_sols * 0.5;
	const size_t num_exploit = num_sols - num_explore;
	rnd->uniform.shuffle(explore_data.begin(), explore_data.end(), num_explore);
	rnd->uniform.shuffle(exploit_data.begin(), exploit_data.end(), num_exploit);
	data.insert(data.end(), explore_data.begin(), explore_data.begin() + num_explore);
	data.insert(data.end(), exploit_data.begin(), exploit_data.begin() + num_exploit);

	{
		auto file_name = g_working_directory + "/result/explore_and_exploit_sol_" + pro_name
			+ "_num_" + std::to_string(num_sols)
			+ ".dat";
		std::ofstream out_file(file_name);
		std::stringstream out;
		size_t id_task = 1;
		for (size_t i = 0; i < num_vars; ++i)
			out << std::setw(10) << ('x' + std::to_string(i + 1));
		out << std::setw(10) << 'y';
		out << '\n';
		out << std::fixed << std::setprecision(2);
		for (auto &c : data) {
			auto s = dynamic_cast<const Solution<>*>(c);
			for (size_t j = 0; j < num_vars; ++j)
				out << std::setw(10) << s->variable()[j];
			out << std::setw(10) << s->objective(0);
			out << '\n';
		}
		out << '\n';
		out_file << out.str();
		out_file.close();
	}

	DEL_PARAM(id_param);
	DEL_PRO(pro);
	DEL_RND(rnd);
}

TEST_CASE("Explore operator", "[explore][operator]") {
	testExploreOperator("Classic_valleys", 20);
	testExploreOperator("Classic_valleys", 50);
	testExploreOperator("Classic_valleys", 100);
	testExploreOperator("Classic_valleys", 500);
}