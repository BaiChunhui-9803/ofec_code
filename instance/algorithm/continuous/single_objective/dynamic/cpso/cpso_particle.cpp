#include "cpso_particle.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
    CPSOParticle::CPSOParticle(size_t num_obj, size_t num_con, size_t size_var):Particle(num_obj, num_con, size_var), m_vmax(size_var){}
    CPSOParticle::~CPSOParticle() = default;
    CPSOParticle::CPSOParticle(const CPSOParticle& rhs):Particle(rhs){
        m_vmax = rhs.m_vmax;
    }

    CPSOParticle& CPSOParticle::operator=(const CPSOParticle& rhs){
        if (this == &rhs) return *this;
        Particle::operator=(rhs);
        m_vmax = rhs.m_vmax;
        return *this;
    }

    void CPSOParticle::initializeVmax(Problem *pro, Random *rnd) {
        for (size_t i = 0; i < variable().size(); i++) {
            auto& range = CAST_CONOP(pro)->range(i);
            m_vmax[i].m_max = (range.second - range.first)/ 5;
            m_vmax[i].m_min = -m_vmax[i].m_max;
        }
    }

    void CPSOParticle::nextVelocity(const Solution<>* lbest, Real w, Real c1, Real c2, Random *rnd) {
        for (size_t j = 0; j < variable().size(); j++) {
            m_vel[j] = w * m_vel[j]
                + c1 * rnd->uniform.next() * (m_pbest.variable()[j] - variable()[j])
                + c2 * rnd->uniform.next() * (lbest->variable()[j] - variable()[j]);

            if (m_vel[j] > m_vmax[j].m_max) m_vel[j] = m_vmax[j].m_max;
            else if (m_vel[j] < m_vmax[j].m_min) m_vel[j] = m_vmax[j].m_min;
        }
    }

    void CPSOParticle::initVelocity(Problem *pro, Random *rnd) {
        for (size_t i = 0; i < variable().size(); i++) {
            auto& range = CAST_CONOP(pro)->range(i);
        	m_vel[i] = (range.second - range.first) * (-0.5 + rnd->uniform.next()) / 2;
        }
    }

    void CPSOParticle::clampVelocity(Problem *pro, Random *rnd){
        for (size_t j = 0; j < variable().size(); j++){
	        auto& range = CAST_CONOP(pro)->range(j);
            Real x = variable()[j] - m_vel[j];
            if(x == range.first || x == range.second){
                Real p = rnd->uniform.next();
                x = p * (m_pbest.variable()[j]) + (1 - p) * (x);
            }
            Real cur_x = x;

            x = (x)+m_vel[j];

            while ((x) > range.second || (x) < range.first) {
                x = (cur_x - range.first < range.second - cur_x) ? cur_x - (cur_x - range.first) * rnd->uniform.next() : cur_x + (range.second - cur_x) * rnd->uniform.next();
            }
            variable()[j] = x;
        }
    }
}
