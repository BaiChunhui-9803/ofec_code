#include "sl_swarm.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	int SwarmSL::ms_updateFre = -1;
	double SwarmSL::ms_learnRatio = -1;
	float SwarmSL::ms_ratioLearnToGbest = -1;

	SwarmSL::SwarmSL(size_t size_pop, Environment *env) : Swarm(size_pop, env), mv_sel(size_pop) {}

	void SwarmSL::updateParameters(Random *rnd) {
		std::vector<int> rindex(m_individuals.size());
		std::iota(rindex.begin(), rindex.end(), 0);
		rnd->uniform.shuffle(rindex.begin(), rindex.end());
		if (ms_learnRatio < 0)
			for (int i = 0; i < m_individuals.size(); i++) {
				m_individuals[rindex[i]]->m_learnRatio = (1 - exp(-pow(i * 1.6 / m_individuals.size(), 4))) > 0.05 ? (1 - exp(-pow(i * 1.6 / m_individuals.size(), 4))) : 0.05;
			}
		else
			for (int i = 0; i < m_individuals.size(); i++) {
				m_individuals[i]->m_learnRatio = ms_learnRatio;
			}

		rnd->uniform.shuffle(rindex.begin(), rindex.end());

		if (ms_updateFre == -1)
			for (int i = 0; i < m_individuals.size(); i++) {
				m_individuals[rindex[i]]->m_updateFre = 10 * exp(-pow(i * 1.6 / m_individuals.size(), 4)) > 1 ? 10 * exp(-pow(i * 1.6 / m_individuals.size(), 4)) : 1;
			}
		else
			for (int i = 0; i < m_individuals.size(); i++) {
				m_individuals[i]->m_updateFre = ms_updateFre;
			}

		updateLearnToGbest(rnd);
	}

	void SwarmSL::setParameters(Random *rnd) {
		for (int i = 0; i < m_individuals.size(); i++) {
			if (ms_learnRatio < 0)
				m_individuals[i]->m_learnRatio = (1 - exp(-pow(i * 1.6 / m_individuals.size(), 4))) > 0.05 ? (1 - exp(-pow(i * 1.6 / m_individuals.size(), 4))) : 0.05;
			else
				m_individuals[i]->m_learnRatio = ms_learnRatio;

			if (ms_updateFre == -1)
				m_individuals[i]->m_updateFre = 10 * exp(-pow(i * 1.6 / m_individuals.size(), 4)) > 1 ? 10 * exp(-pow(i * 1.6 / m_individuals.size(), 4)) : 1;
			else
				m_individuals[i]->m_updateFre = ms_updateFre;
		}

		m_numLearnToGbest = m_individuals.size() * ms_ratioLearnToGbest;
		MRandToBest(m_numLearnToGbest, rnd);
		for (int i = 0; i < m_individuals.size(); i++) m_individuals[i]->setSelRatio();
	}

	void SwarmSL::resize(size_t size_pop, Environment* env) {
		Swarm::resize(size_pop, env);
		mv_sel.resize(size_pop);
	}

	void SwarmSL::updateLearnToGbest(Random *rnd) {
		int PrenumLearnToGbest = 0, CurnumLearnToGbest = 0;
		for (int i = 0; i < m_individuals.size(); i++)
			if (m_individuals[i]->type())
				PrenumLearnToGbest++;

		std::vector<bool> preLearnToGbestInfor(m_individuals.size());
		for (int i = 0; i < m_individuals.size(); i++) preLearnToGbestInfor[i] = m_individuals[i]->type();
		std::vector<int> preLearnToGbest;
		int k = 0;
		if (PrenumLearnToGbest > 0) {
			preLearnToGbest.resize(PrenumLearnToGbest);
			k = 0;
			for (int i = 0; i < m_individuals.size(); i++) {
				if (m_individuals[i]->type()) preLearnToGbest[k++] = i;
				if (k >= PrenumLearnToGbest) break;
			}
		}

		MRandToBest(m_numLearnToGbest, rnd);

		for (int i = 0; i < m_individuals.size(); i++) if (m_individuals[i]->type()) CurnumLearnToGbest++;
		std::vector<int> curLearnToGbest;

		if (CurnumLearnToGbest > 0) {
			curLearnToGbest.resize(CurnumLearnToGbest);
			k = 0;
			for (int i = 0; i < m_individuals.size(); i++) {
				if (m_individuals[i]->type()) curLearnToGbest[k++] = i;
				if (k >= CurnumLearnToGbest) break;
			}
		}


		for (int i = 0; i < PrenumLearnToGbest; i++) {
			if (!m_individuals[preLearnToGbest[i]]->type()) 		m_individuals[preLearnToGbest[i]]->learnToNonLearn();

		}
		for (int i = 0; i < CurnumLearnToGbest; i++) {
			if (!preLearnToGbestInfor[curLearnToGbest[i]])	m_individuals[curLearnToGbest[i]]->nonLearnToLearn();

		}

	}

	int SwarmSL::evolve(Environment *env, Random *rnd) {
		int rf = kNormalEval;
		if (m_individuals.size() <= 0)
			return kTerminate;

		Solution<>  t_p;
		int numDim = env->problem()->numberVariables();
		std::vector<double> avg_v (numDim);

		for (int j = 0; j < numDim; j++) {
			avg_v[j] = 0;
			for (int i = 0; i < m_individuals.size(); i++) {
				avg_v[j] += fabs(m_individuals[i]->velocity()[j]);
			}
			avg_v[j] /= m_individuals.size();
		}
		int k;
		std::vector<int> rindex(m_individuals.size());
		std::iota(rindex.begin(), rindex.end(), 0);
		rnd->uniform.shuffle(rindex.begin(), rindex.end());

		for (int j = 0; j < m_individuals.size(); j++) {
			int i = rindex[j];
			k = i;
			if (m_individuals[i]->m_itersUnimpr > m_individuals[i]->m_updateFre) {
				m_individuals[i]->updateSelectionRatioProg(rnd);
				m_individuals[i]->m_itersUnimpr = 0;
			}
			mv_sel[i] = m_individuals[i]->selectOperator(rnd);
			t_p = *m_individuals[i];
			int nearest;
			switch (mv_sel[i] + 1) {
			case 1:
				m_individuals[i]->nextVelocity(&m_individuals[i]->pbest(), m_weight, m_accelerator1, m_accelerator2, rnd);
				m_individuals[i]->move();
				m_individuals[i]->clampVelocity(env, rnd);
				rf = m_individuals[i]->evaluate(env);
				if (rf != kNormalEval) {
					avg_v.clear();
		//			delete[] avg_v;
		//			avg_v = 0;
					return rf;
				}
				break;
			case 2:
				m_individuals[i]->normalMutation(avg_v.data(), rnd, env);
				rf = m_individuals[i]->evaluate(env);
				if (rf != kNormalEval) {
					avg_v.clear();
					//			delete[] avg_v;
					//			avg_v = 0;
					return rf;
				}
				break;
			case 3:
				nearest = rnd->uniform.nextNonStd<int>(0, m_individuals.size());
				if (m_individuals.size() > 1) {
					while (nearest == i) {
						nearest = rnd->uniform.nextNonStd<int>(0, m_individuals.size());
					}
				}
				if (dominate(m_individuals[i]->pbest(), m_individuals[nearest]->pbest(), env->problem()->optimizeMode())) {
					i = nearest;
					t_p = *m_individuals[i];
					m_individuals[i]->nextVelocity(&m_individuals[k]->pbest(), m_weight, m_accelerator1, m_accelerator2, rnd);
					m_individuals[i]->move();
					m_individuals[i]->clampVelocity(env, rnd);
					rf = m_individuals[i]->evaluate(env);
					if (rf != kNormalEval) {
						avg_v.clear();
						//			delete[] avg_v;
						//			avg_v = 0;
						return rf;
					}
				}
				else {
					m_individuals[i]->nextVelocity(&m_individuals[nearest]->pbest(), m_weight, m_accelerator1, m_accelerator2, rnd);
					m_individuals[i]->move();
					m_individuals[i]->clampVelocity(env, rnd);
					rf = m_individuals[i]->evaluate(env);
					if (rf != kNormalEval) {
						avg_v.clear();
						//			delete[] avg_v;
						//			avg_v = 0;
						return rf;
					}
				}
				break;
			case 4:
				m_individuals[i]->nextVelocity(&this->m_best_individual->pbest(), m_weight, m_accelerator1, m_accelerator2, rnd);
				m_individuals[i]->move();
				m_individuals[i]->clampVelocity(env, rnd);
				rf = m_individuals[i]->evaluate(env);
				if (rf != kNormalEval) {
					avg_v.clear();
					//			delete[] avg_v;
					//			avg_v = 0;
					return rf;
				}
				break;
			default:
				avg_v.clear();
				//			delete[] avg_v;
				//			avg_v = 0;
				throw Exception("operator selection error @SwarmSL::evolve()");
				break;
			}

			if (dominate(*m_individuals[i], m_individuals[i]->pbest(), env->problem()->optimizeMode())) {
				m_individuals[i]->pbest() = *m_individuals[i];
				if (dominate(*m_individuals[i], this->m_best_individual->pbest(), env->problem()->optimizeMode())) {
					this->m_best_individual = m_individuals[i].get();
				}
			}

			m_individuals[i]->mv_prog[mv_sel[i]].m_numSelected++;

			if (dominate(*m_individuals[i], t_p, env->problem()->optimizeMode())) {
				m_individuals[i]->mv_prog[mv_sel[i]].m_numSuccess++;
				m_individuals[i]->mv_prog[mv_sel[i]].m_rewards += fabs(t_p.objective(0) - m_individuals[i]->objective(0));
			}
			else {
				m_individuals[i]->m_itersUnimpr++;
			}

			if (dominate(*m_individuals[i], t_p, env->problem()->optimizeMode())) {
				rf = updateBestByRatio(i, m_individuals[i]->m_learnRatio, env, rnd);
				if (rf != kNormalEval) {
					avg_v.clear();
					//			delete[] avg_v;
					//			avg_v = 0;
					return rf;
				}
			}
		}

		updateParameters(rnd);
		m_iteration++;

		avg_v.clear();
		//			delete[] avg_v;
		//			avg_v = 0;

		return rf;
	}


	int SwarmSL::getNumLearnTogbest() {

		int n = 0;
		for (int i = 0; i < m_individuals.size(); i++) if (m_individuals[i]->type()) n++;

		return n;
	}

	void SwarmSL::MRandToBest(int num, Random *rnd) {
		if (num <= 0) {
			for (int i = 0; i < m_individuals.size(); i++)	m_individuals[i]->setType(0);
			return;
		}
		if (num > m_individuals.size()) num = m_individuals.size();

		for (int i = 0; i < m_individuals.size(); i++)	m_individuals[i]->setType(0);
		std::vector<int> l;
		for (int i = 0; i < m_individuals.size(); i++) l.push_back(i);

		for (int j = 0; j < num; j++) {
			int k = rnd->uniform.nextNonStd<int>(0, l.size() - 1);
			m_individuals[l[k]]->setType(1);
			l.erase(l.begin() + k);
		}
		l.clear();
	}

	int SwarmSL::updateBestByRatio(const int p, double ratio, Environment *env, Random *rnd) {
		int rf = kNormalEval;
		Solution<> x;
		for (int j = 0; j < this->m_best_individual->pbest().variable().size(); j++) {
			if (ratio<1 && rnd->uniform.next() >ratio) continue;
			x = this->m_best_individual->pbest();
			x.variable()[j] = this->m_individuals[p]->variable()[j];
			rf = x.evaluate(env);

			if (dominate(x, this->m_best_individual->pbest(), env->problem()->optimizeMode()))   
				this->m_best_individual->pbest() = x;
			if (rf != kNormalEval) break;
		}
		return rf;
	}

	void SwarmSL::calculateNumLearning(const int sfes, Algorithm *alg) {
		if (m_individuals.size() == 0) return;
		if (ms_ratioLearnToGbest < 0) {
			int cur_evals = alg->evaluations();
			int max_evals = alg->maximumEvaluations();
			m_numLearnToGbest = (int)(m_individuals.size() * (1 - exp(-100 * pow((double)(cur_evals - sfes) / (max_evals - sfes), 3))));
		}
	}
}
