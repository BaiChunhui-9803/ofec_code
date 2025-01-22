#include "mqso_particle.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
    //MQSOParticle::MQSOParticle() :
    //    MQSOParticle(CONTINUOUS_CAST->objective_size(), CONTINUOUS_CAST->num_constraints(), CONTINUOUS_CAST->variable_size()) {}

    MQSOParticle::MQSOParticle(size_t num_obj, size_t num_con, size_t size_var) :
        Particle(num_obj, num_con, size_var),
        mv_repulse(size_var) {}

    MQSOParticle::~MQSOParticle() {

    }
    MQSOParticle::MQSOParticle(const MQSOParticle& other):Particle(other)
    {
        mv_repulse.resize(other.variable().size());
        m_type = other.m_type;
        copy(other.mv_repulse.begin(), other.mv_repulse.end(), mv_repulse.begin());
    }

    void MQSOParticle::nextVelocity(const Solution<>* lbest, Real w, Real c1, Real c2, Random *rnd) {
        for (size_t j = 0; j < variable().size(); j++) {
            if (m_type == ParticleType::PARTICLE_NEUTRAL) {
                m_vel[j] = w * m_vel[j] + c1 * rnd->uniform.next() * (m_pbest.variable()[j] - variable()[j]) + c2 * rnd->uniform.next() * (lbest->variable()[j] - variable()[j]);
            }
            else if (m_type == ParticleType::PARTICLE_CHARGED) {
                m_vel[j] = w * m_vel[j] + c1 * rnd->uniform.next() * (m_pbest.variable()[j] - variable()[j]) + c2 * rnd->uniform.next() * (lbest->variable()[j] - variable()[j]) + mv_repulse[j];
            }
            if (m_vel[j] > m_vMax[j].m_max)	m_vel[j] = m_vMax[j].m_max;
            else if (m_vel[j] < m_vMax[j].m_min)		m_vel[j] = m_vMax[j].m_min;
        }
    }
    void MQSOParticle::move() {
        if (m_type == ParticleType::PARTICLE_NEUTRAL || m_type == ParticleType::PARTICLE_CHARGED) {
            for (size_t j = 0; j < variable().size(); j++) {
                variable()[j] += m_vel[j];
            }
        }
    }
    void MQSOParticle::quantumMove(const Solution<>* lbest, Real m_Rcloud, Random *rnd)
    {
        for (size_t j = 0; j < variable().size(); j++) {
            Real temp = rnd->uniform.nextNonStd<Real>(-1, 1);
            variable()[j] = lbest->variable()[j] + m_Rcloud * temp;
        }
    }

    void MQSOParticle::initializeVmax(Problem *pro) {
        for (int i = 0; i < CAST_CONOP(pro)->numberVariables(); i++) {
            double l, u; 
            l = CAST_CONOP(pro)->range(i).first;
            u = CAST_CONOP(pro)->range(i).second;
            // for static optimization problems

            m_vMax.resize(CAST_CONOP(pro)->numberVariables());
            if (pro->name().find("FUN_") != std::string::npos)
                m_vMax[i].m_max = (u - l) / 2;
            else // for dynamic problems
                m_vMax[i].m_max = (u - l) / 5;
            m_vMax[i].m_min = -m_vMax[i].m_max;

        }
    }
    void MQSOParticle::initVelocity(Problem *pro, Random *rnd) {
        for (size_t i = 0; i < variable().size(); i++) {
            auto& range = CAST_CONOP(pro)->range(i);
            m_vel[i] = (range.second - range.first) * (-0.5 + rnd->uniform.next()) / 2;
        }
    }
}