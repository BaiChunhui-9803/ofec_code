#include "sl_particle.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include <algorithm>
#include "../../../../../../core/environment/environment.h"

namespace ofec {
	const int ParticleSL::ms_numOperators = 4;

	ParticleSL::ParticleSL() : 
		Particle(), 
		m_itersUnimpr(0), 
		m_updateFre(0), 
		m_learnRatio(0),
		mv_prog(ms_numOperators),
		mv_monitor(ms_numOperators) {}

	ParticleSL::ParticleSL(size_t num_obj, size_t num_con, size_t size_var) :
		Particle(num_obj, num_con, size_var),
		m_itersUnimpr(0),
		m_updateFre(0),
		m_learnRatio(0),
		mv_prog(ms_numOperators),
		mv_monitor(ms_numOperators)	{}

	ParticleSL::ParticleSL(const ParticleSL &rhs) :
		Particle(rhs),
		m_itersUnimpr(rhs.m_itersUnimpr),
		m_updateFre(rhs.m_updateFre),
		m_learnRatio(rhs.m_learnRatio),
		mv_prog(ms_numOperators),
		mv_monitor(ms_numOperators)
	{
		std::copy(rhs.mv_prog.begin(), rhs.mv_prog.end(), mv_prog.begin());
		std::copy(rhs.mv_monitor.begin(), rhs.mv_monitor.end(), mv_monitor.begin());
	}

	ParticleSL::ParticleSL(ParticleSL &&rhs) noexcept :
		Particle(std::move(rhs)),
		m_itersUnimpr(rhs.m_itersUnimpr),
		m_updateFre(rhs.m_updateFre),
		m_learnRatio(rhs.m_learnRatio),
		mv_prog(ms_numOperators),
		mv_monitor(ms_numOperators)
	{
		std::copy(rhs.mv_prog.begin(), rhs.mv_prog.end(), mv_prog.begin());
		std::copy(rhs.mv_monitor.begin(), rhs.mv_monitor.end(), mv_monitor.begin());
	}


	ParticleSL& ParticleSL::operator=(const ParticleSL &other) {
		if (this != &other) {
			Particle::operator=(other);
			m_itersUnimpr = other.m_itersUnimpr;
			m_updateFre = other.m_updateFre;
			m_learnRatio = other.m_learnRatio;
			std::copy(other.mv_prog.begin(), other.mv_prog.end(), mv_prog.begin());
			std::copy(other.mv_monitor.begin(), other.mv_monitor.end(), mv_monitor.begin());
		}
		return *this;
	}

	ParticleSL& ParticleSL::operator=(ParticleSL &&other) noexcept {
		if (this != &other) {
			Particle::operator=(std::move(other));
			m_itersUnimpr = other.m_itersUnimpr;
			m_updateFre = other.m_updateFre;
			m_learnRatio = other.m_learnRatio;
			std::copy(other.mv_prog.begin(), other.mv_prog.end(), mv_prog.begin());
			std::copy(other.mv_monitor.begin(), other.mv_monitor.end(), mv_monitor.begin());
		}
		return *this;
	}

	void ParticleSL::setSelRatio() {
		for (int i = 0; i < ms_numOperators; i++) {
			if (m_type) {
				mv_prog[i].m_ratio = 1. / ms_numOperators;
			}
			else {
				if (i == ms_numOperators - 1) mv_prog[i].m_ratio = 0.0;
				else mv_prog[i].m_ratio = 1. / (ms_numOperators - 1);
			}

			mv_prog[i].m_minRatio = 0.001;
			mv_prog[i].m_numSelected = 0;
			mv_prog[i].m_numSuccess = 0;
			mv_prog[i].m_rewards = 0;
			mv_monitor[i] = mv_prog[i];
		}
	}

	void ParticleSL::nonLearnToLearn() {
		for (int m = 0; m < ms_numOperators; m++) {
			mv_prog[m].initialize(ms_numOperators);
		}
		double sum = 0;
		for (int m = 0; m < ms_numOperators - 1; m++) sum += mv_monitor[m].m_ratio;
		int m;
		for (m = 0; m < ms_numOperators - 1; m++) {
			mv_monitor[m].m_ratio = (double)(ms_numOperators - 1.) / ms_numOperators * (mv_monitor[m].m_ratio / sum);
		}
		mv_monitor[m].initialize(ms_numOperators);
	}

	void ParticleSL::learnToNonLearn() {
		double sum = 0;
		for (int j = 0; j < ms_numOperators - 1; j++) {
			sum += mv_prog[j].m_ratio;
		}
		for (int j = 0; j < ms_numOperators - 1; j++) {
			mv_prog[j].m_ratio = mv_prog[j].m_ratio / sum;
		}
		mv_prog[ms_numOperators - 1].m_ratio = 0;
		sum = 0;
		for (int j = 0; j < ms_numOperators - 1; j++) {
			sum += mv_monitor[j].m_ratio;
		}
		for (int j = 0; j < ms_numOperators - 1; j++) {
			mv_monitor[j].m_ratio = mv_monitor[j].m_ratio / sum;
		}
		mv_monitor[ms_numOperators - 1].m_ratio = 0;

	}

	void ParticleSL::updateSelectionRatioMonitor(Random *rnd) {
		if (m_type) {
			slpso::Progress::updateProgress(rnd, mv_monitor, ms_numOperators);
		}
		else {
			slpso::Progress::updateProgress(rnd, mv_monitor, ms_numOperators - 1);
		}

		for (int j = 0; j < ms_numOperators; j++) {
			mv_monitor[j].m_numSelected = 0;
			mv_monitor[j].m_numSuccess = 0;
			mv_monitor[j].m_rewards = 0;
		}
	}

	void ParticleSL::updateSelectionRatioProg(Random *rnd) {
		if (m_type) {
			slpso::Progress::updateProgress(rnd, mv_prog, ms_numOperators);
		}
		else {
			slpso::Progress::updateProgress(rnd, mv_prog, ms_numOperators - 1);
		}
		for (int j = 0; j < ms_numOperators; j++) {
			mv_prog[j].m_numSelected = 0;
			mv_prog[j].m_numSuccess = 0;
			mv_prog[j].m_rewards = 0;
		}
	}

	int ParticleSL::selectOperator(Random *rnd) {
		if (m_type) {
			return slpso::Progress::getAction(rnd, ms_numOperators, mv_prog);
		}
		else {
			return slpso::Progress::getAction(rnd, ms_numOperators - 1, mv_prog);
		}
	}

	void ParticleSL::nextVelocity(const Solution<> *lbest, Real w, Real c1, Real c2, Random *rnd) {
		for (size_t j = 0; j < variable().size(); j++) {
			m_vel[j] = w * m_vel[j] + c1 * rnd->uniform.next() * (lbest->variable()[j] - variable()[j]);
			if (m_vel[j] > m_vel_max[j].max) {
				m_vel[j] = m_vel_max[j].max;
			}
			else if (m_vel[j] < m_vel_max[j].min) {
				m_vel[j] = m_vel_max[j].min;
			}
		}
	}

	void ParticleSL::normalMutation(double *avg_v, Random *rnd, Environment *env) {
		for (size_t j = 0; j < variable().size(); j++) {
			variable()[j] += avg_v[j] * rnd->uniform.next();
			auto &range = CAST_CONOP(env->problem())->range(j);
			if (variable()[j] > range.second || variable()[j] < range.first) {
				variable()[j] = rnd->uniform.nextNonStd(range.first, range.second);
			}
		}
	}
}
