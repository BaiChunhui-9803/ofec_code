#include "nbn_utility_continous.h"
#include <algorithm>
#include "../../../../core/problem/continuous/continuous.h"
#include "../../../../run/interface.h"
#include "../../../function/custom_function.h"
namespace ofec {
	namespace nbn {
		namespace continous {

			void fitlerSameSolutions(std::vector<ofec::SolutionBase*>& sols,
				std::vector<ofec::SolutionBase*>& filterSols, std::vector<int>& mapSolId, ofec::Environment* env) {
				auto pro=  env->problem();
			
				using SolutionType = Continuous::SolutionType;
					std::vector<int> sortedIds(sols.size());
					for (int idx(0); idx < sortedIds.size(); ++idx) {
						sortedIds[idx] = idx;
					}
					std::vector<SolutionType*> cursols(sols.size());
					for (int idx(0); idx < cursols.size(); ++idx) {
						cursols[idx] = dynamic_cast<SolutionType*>(sols[idx]);
					}
					std::sort(sortedIds.begin(), sortedIds.end(), [&](int a, int b) {
						if (sols[a]->fitness() == sols[b]->fitness()) {
							if (cursols[a]->variable().vector() == cursols[b]->variable().vector()) {
								return a < b;
							}
							else return cursols[a]->variable().vector() < cursols[b]->variable().vector();
						}
						else return sols[a]->fitness() < sols[b]->fitness();
						});

					filterSols.push_back(sols[sortedIds.front()]);
					mapSolId.resize(sols.size());
					//originIds.push_back(sortedIds.front());
					int mapsolid = 0;
					mapSolId[sortedIds.front()] = mapsolid;
					//mapSolId[solId] = otherId;
					for (int idx(1); idx < sortedIds.size(); ++idx) {
						bool isSame = pro->same(sols[sortedIds[idx - 1]]->variableBase(), sols[sortedIds[idx]]->variableBase());
						if (!isSame) {
							filterSols.push_back(sols[sortedIds[idx]]);
							//	originIds.push_back(sortedIds[idx]);
							mapsolid++;
						}
						mapSolId[sortedIds[idx]] = mapsolid;
					}
			}

			void initSolutionsRadius(ofec::SolutionBase& curSol, const std::vector<std::pair<Real, Real>>& initBoundary, ofec::Random* rnd) {
				auto& curConSol = dynamic_cast<ofec::Continuous::SolutionType&>(curSol);
				for (int idx(0); idx < initBoundary.size(); ++idx) {
					curConSol.variable()[idx] = rnd->uniform.nextNonStd<double>(initBoundary[idx].first, initBoundary[idx].second);
				}

			}

			void generateSolutionsConRadiusOneThread(
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				int from, int to,
				ofec::Environment* env,
				ofec::Random* rnd
			) {
				using namespace std;
				using namespace ofec;
				auto conPro = CAST_CONOP(env->problem());
				auto& curCsol = dynamic_cast<const ofec::Continuous::SolutionType&>(centerSol);

				auto boundary = conPro->boundary();
				auto initRange = conPro->boundary();
				
				for (int idx(0); idx < initRange.size(); ++idx) {

					double halfRadius = (initRange[idx].second - initRange[idx].first) * sampleRadius / 2.0;
					initRange[idx].first = curCsol.variable()[idx] - halfRadius;
					initRange[idx].second = curCsol.variable()[idx] + halfRadius;
				}

				for (int idx(0); idx < initRange.size(); ++idx) {
					if (initRange[idx].first < boundary[idx].first) {
						initRange[idx].first = boundary[idx].first;
					}
					if (initRange[idx].second > boundary[idx].second) {
						initRange[idx].second = boundary[idx].second;
					}

				}


				std::vector<int> data;
				for (int idx(from); idx < to; ++idx) {
					sols[idx].reset(env->problem()->createSolution(centerSol));
					initSolutionsRadius(*sols[idx], initRange, rnd);
				}
			}



			void generateSolutionsConRadius(
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				ofec::Environment* env,
				ofec::Random* rnd
			) {
				using namespace std;
				using namespace ofec;
			//	sols.resize(1e6);

				int num_task = std::thread::hardware_concurrency();
				std::cout << "generateSolutionsConRadius num task\t" << num_task << std::endl;
				std::vector<std::shared_ptr<Environment>> envs(num_task);
				std::vector<std::shared_ptr<Random>> rnds(num_task);
				for (auto& it : envs) {
					genenrateProblem(env->problem()->name(), env->problem()->archivedParameters(),it);
				}
				for (auto& it : rnds) {
					it.reset(new Random(rnd->uniform.next()));
				}
				std::vector<int> tasks;
				std::vector<std::thread> thrds;

				UTILITY::assignThreads(sols.size(), num_task, tasks);
				std::pair<int, int> from_to;
				for (size_t i = 0; i < num_task; ++i) {
					from_to.first = tasks[i];
					from_to.second = tasks[i + 1];

					thrds.push_back(std::thread(
						generateSolutionsConRadiusOneThread, std::cref(centerSol), sampleRadius, std::ref(sols),
						tasks[i], tasks[i + 1], envs[i].get(), rnds[i].get()));
				}
				for (auto& thrd : thrds)
					thrd.join();
			}

			void genenrateProblem(
				const std::string& proname,
				const ofec::ParameterMap& params,
				std::shared_ptr<ofec::Environment>& env) {
				using namespace ofec;


				env.reset(Environment::create());
				env->recordInputParameters();
				env->initialize(0.5);
				env->setProblem(ofec::Factory<ofec::Problem>::produce(proname));
				env->problem()->inputParameters().input(params);
				env->problem()->setName(proname);
				env->problem()->recordInputParameters();
				env->initializeProblem(0.5);
			}
		}
	}
}