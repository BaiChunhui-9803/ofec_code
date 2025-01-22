#include "cswidn_utility.h"
#include "../../../../instance/problem/realworld/csiwdn/csiwdn.h"
#include "../../../function/custom_function.h"
#include "../../../../run/interface.h"

namespace ofec {
	namespace nbn {
		namespace cswidn {


			// 函数用于将十进制数转换为M进制数
			void decimalToM(int decimal, int M, int numSource, std::vector<int>& data) {

				data.clear();
				// 持续除以M，直到decimal为0
				while (numSource--) {
					int remainder = decimal % M;
					data.push_back(remainder);
					//char symbol = (remainder < 10) ? ('0' + remainder) : ('A' + (remainder - 10));
					//result = symbol + result; // 将余数添加到结果的前面
					decimal /= M; // 更新十进制数
				}

			}

			void copyCSWIDN_Environment(std::shared_ptr<ofec::Environment>& envCopy, ofec::Environment* env) {
				using namespace std;
				using namespace ofec;
				ofec::ParameterMap params;

				env->problem()->inputParameters().output(params);

				envCopy.reset(Environment::create());
				envCopy->recordInputParameters();
				envCopy->initialize();
				envCopy->setProblem(ofec::Factory<ofec::Problem>::produce(env->problem()->name()));
				envCopy->problem()->inputParameters().input(params);
				envCopy->problem()->recordInputParameters();
				envCopy->initializeProblem(0.5);
				auto cswidnPro = CAST_CSIWDN(env->problem());
				auto cswidnProCopy = CAST_CSIWDN(envCopy->problem());

				//std::string filernd = UTILITY::getCurrentSystemTime() + "_" + UTILITY::createRandomString(50, m_rnd.get());
				//auto filernd = m_winName + "_" + std::to_string(total_file++);
				//cswidnPro->setFileRnd(filernd);

				cswidnProCopy->setDistanceType(cswidnPro->distanceType());
				cswidnProCopy->setPhase(cswidnPro->phase());
				cswidnProCopy->setSourceIdx(cswidnPro->sourceIdx());
				cswidnProCopy->setProcessName(cswidnPro->processName());



			}



			void generateSolRasiu(std::shared_ptr<ofec::SolutionBase>& initSol,
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				ofec::Environment* env,
				ofec::Random* rnd) {

				using namespace std;
				using namespace ofec;
				auto cswidnPro = CAST_CSIWDN(env->problem());
				auto& curCsol = dynamic_cast<const ofec::CSIWDN::solutionType&>(centerSol);

				double halfRadius = sampleRadius / 2.0 * (cswidnPro->maxMultiplier() - cswidnPro->minMultiplier());

				std::pair<double, double> proRange = { cswidnPro->minMultiplier(), cswidnPro->maxMultiplier() };

				std::vector<std::vector<std::pair<double, double>>> targetRange(cswidnPro->numSource());
				for (int idx(0); idx < targetRange.size(); ++idx) {
					targetRange[idx].resize(curCsol.variable().multiplier(idx).size());

					for (int idy(0); idy < targetRange[idx].size(); ++idy) {
						targetRange[idx][idy].first = curCsol.variable().multiplier(idx)[idy] - halfRadius;
						targetRange[idx][idy].second = curCsol.variable().multiplier(idx)[idy] + halfRadius;
					}
				}

				for (auto& it : targetRange) {
					for (auto& it2 : it) {
						if (it2.first < proRange.first) {
							double left = proRange.first - it2.first;
							it2.first = proRange.first;
							it2.second += left;
						}
						else if (it2.second > proRange.second) {
							double left = it2.second - proRange.second;
							it2.second = proRange.second;
							it2.first -= left;
						}
					}
				}



				//	initSol.reset(env->problem()->createSolution(centerSol));
				auto& cursol = dynamic_cast<ofec::CSIWDN::solutionType&>(*initSol);
				//	cursol.variable().flagLocation() = false;
					//cswidnPro->initSolutionMultiplier(cursol, targetRange, rnd);

				VarCSIWDN& var = cursol.variable();
				for (size_t z = 0; z < cswidnPro->numSource(); z++) {
					for (int j(0); j < var.multiplier(z).size(); ++j) {
						if (rnd->uniform.next() < 0.5) {
							var.multiplier(z)[j] = targetRange[z][j].first;
						}
						else var.multiplier(z)[j] = targetRange[z][j].second;
						//var.multiplier(z)[j] = rnd->uniform.nextNonStd<double>(range[z][j].first, range[z][j].second);
					}
					//for (auto& j : var.multiplier(z)) {
					//	j = rnd->uniform.nextNonStd<float>(m_min_multiplier, m_max_multiplier);
					//}
				}
			}



			void generateSolIntRasiu(std::shared_ptr<ofec::SolutionBase>& initSol,
				const ofec::SolutionBase& centerSol,
				int sampleRadius,
				ofec::Environment* env,
				ofec::Random* rnd) {

				using namespace std;
				using namespace ofec;
				auto cswidnPro = CAST_CSIWDN(env->problem());
				auto& curCsol = dynamic_cast<const ofec::CSIWDN::solutionType&>(centerSol);

				auto numSource = cswidnPro->numberSource();
				std::vector<int> shuffleIds(numSource);
				for (int idx(0); idx < shuffleIds.size(); ++idx) {
					shuffleIds[idx] = idx;
				}
				rnd->uniform.shuffle(shuffleIds.begin(), shuffleIds.end());
				shuffleIds.resize(sampleRadius);

				auto& cursol = dynamic_cast<ofec::CSIWDN::solutionType&>(*initSol);
				//cursol.variable().flagLocation() = false;


				for (auto z : shuffleIds) {
					int beforeNode = cursol.variable().index(z);
					cursol.variable().index(z) = rnd->uniform.nextNonStd<int>(1, cswidnPro->numberNode());
					if (cursol.variable().index(z) >= beforeNode) ++cursol.variable().index(z);
				}

			}








			void generateSolutionsIntsSub(
				const ofec::SolutionBase& centerSol,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, int from, int to, ofec::Environment* env, ofec::Random* rnd) {
				using namespace std;
				using namespace ofec;

				auto cswidnPro = CAST_CSIWDN(env->problem());

				std::vector<int> data;
				for (int idx(from); idx < to; ++idx) {
					sols[idx].reset(env->problem()->createSolution(centerSol));

					auto& cursol = dynamic_cast<ofec::CSIWDN::solutionType&>(*sols[idx]);
					cursol.variable().flagLocation() = false;

					decimalToM(idx, cswidnPro->numberNode(), cswidnPro->numSource(), data);
					for (auto& it : data) {
						++it;
					}
					for (size_t z = 0; z < cswidnPro->numSource(); z++) {
						cursol.variable().index(z) = data[z];
					}

					//	cursol.evaluate(env, false);
				}

			}



			void generateSolutionsConsRandomSub(
				const ofec::SolutionBase& centerSol,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, 
				int from, int to, ofec::Environment* env, ofec::Random* rnd) {
				using namespace std;
				using namespace ofec;

				auto cswidnPro = CAST_CSIWDN(env->problem());

				std::vector<int> data;
				for (int idx(from); idx < to; ++idx) {
					sols[idx].reset(env->problem()->createSolution(centerSol));
					auto& cursol = dynamic_cast<ofec::CSIWDN::solutionType&>(*sols[idx]);
					cursol.variable().flagLocation() = false;
					cswidnPro->initSolutionMultiplier(cursol, rnd);
					//			CAST_CSIWDN(pro)->initSolutionMultiplier(cursol, rnd);
				}

			}



			void generateSolutionsConsRadiusSub(
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				int from, int to, ofec::Environment* env, ofec::Random* rnd) {
				using namespace std;
				using namespace ofec;
				auto cswidnPro = CAST_CSIWDN(env->problem());
				auto& curCsol = dynamic_cast<const ofec::CSIWDN::solutionType&>(centerSol);

				auto& optSol = cswidnPro->optima()->solution(0);

				double halfRadius = sampleRadius / 2.0 * (cswidnPro->maxMultiplier() - cswidnPro->minMultiplier());

				std::pair<double, double> proRange = { cswidnPro->minMultiplier(), cswidnPro->maxMultiplier() };

				std::vector<std::vector<std::pair<double, double>>> targetRange(cswidnPro->numSource());
				for (int idx(0); idx < targetRange.size(); ++idx) {
					targetRange[idx].resize(optSol.variable().multiplier(idx).size());

					for (int idy(0); idy < targetRange[idx].size(); ++idy) {
						targetRange[idx][idy].first = optSol.variable().multiplier(idx)[idy] - halfRadius;
						targetRange[idx][idy].second = optSol.variable().multiplier(idx)[idy] + halfRadius;
					}
				}

				for (auto& it : targetRange) {
					for (auto& it2 : it) {
						if (it2.first < proRange.first) {
							//	double left = proRange.first - it2.first;
							it2.first = proRange.first;
							//	it2.second += left;
						}
						else if (it2.second > proRange.second) {
							//	double left = it2.second - proRange.second;
							it2.second = proRange.second;
							//	it2.first -= left;
						}
					}
				}



				std::vector<int> data;
				for (int idx(from); idx < to; ++idx) {
					sols[idx].reset(env->problem()->createSolution(centerSol));
					auto& cursol = dynamic_cast<ofec::CSIWDN::solutionType&>(*sols[idx]);
					cursol.variable().flagLocation() = false;
					cswidnPro->initSolutionMultiplier(cursol, targetRange, rnd);
					//			CAST_CSIWDN(pro)->initSolutionMultiplier(cursol, rnd);
				}

			}


			void generateSolutionsInts(
				const ofec::SolutionBase& centerSol,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, 
				ofec::Environment* env,
				ofec::Random *rnd
				) {
				using namespace std;
				using namespace ofec;


				auto cswidnPro = CAST_CSIWDN(env->problem());
				sols.resize(pow(cswidnPro->numberNode(), cswidnPro->numSource()));
				int num_task = std::thread::hardware_concurrency();
				// for test 
				//num_task = 1;
				std::cout << "generateSolutionsInts num task\t" << num_task << std::endl;
				std::vector<std::shared_ptr<Environment>> envs(num_task);
				std::vector<std::shared_ptr<Random>> rnds(200);
				for (auto& it : envs) {
					copyCSWIDN_Environment(it, env);
				//	auto curpro = CAST_CSIWDN(it->problem());
				//	curpro->setPhase(cswidnPro->phase());
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
						generateSolutionsIntsSub, std::cref(centerSol), std::ref(sols),
						tasks[i], tasks[i + 1], envs[i].get(), rnds[i].get()));
				}
				for (auto& thrd : thrds)
					thrd.join();


			}
			void generateSolutionsConRandom(
				const ofec::SolutionBase& centerSol,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				ofec::Environment* env,
				ofec::Random* rnd
			) {
				using namespace std;

				auto cswidnPro = CAST_CSIWDN(env->problem());

				int num_task = std::thread::hardware_concurrency();

				std::vector<std::shared_ptr<Environment>> envs(num_task);
				std::vector<std::shared_ptr<Random>> rnds(200);
				for (auto& it : envs) {
					copyCSWIDN_Environment(it,env);
					//auto curpro = CAST_CSIWDN(it->problem());
					//curpro->setPhase(cswidnPro->phase());
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
						generateSolutionsConsRandomSub, std::cref(centerSol), std::ref(sols),
						tasks[i], tasks[i + 1], envs[i].get(), rnds[i].get()));
				}
				for (auto& thrd : thrds)
					thrd.join();
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
				auto cswidnPro = CAST_CSIWDN(env->problem());
				sols.resize(1e6);



				int num_task = std::thread::hardware_concurrency();
				std::cout << "generateSolutionsConRadius num task\t" << num_task << std::endl;
				std::vector<std::shared_ptr<Environment>> envs(num_task);
				std::vector<std::shared_ptr<Random>> rnds(num_task);
				for (auto& it : envs) {
					copyCSWIDN_Environment(it,env);
					//auto curpro = CAST_CSIWDN(it->problem());
					//curpro->setPhase(cswidnPro->phase());
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
						generateSolutionsConsRadiusSub, std::cref(centerSol), sampleRadius, std::ref(sols),
						tasks[i], tasks[i + 1], envs[i].get(), rnds[i].get()));
				}
				for (auto& thrd : thrds)
					thrd.join();
			}


		}
	}
}