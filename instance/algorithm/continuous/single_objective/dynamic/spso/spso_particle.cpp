#include "spso_particle.h"
#include "../../../../../../core/problem/continuous/continuous.h"

namespace ofec {
    SPSOParticle::SPSOParticle(size_t num_obj, size_t num_con, size_t size_var):Particle(num_obj, num_con, size_var){}
    SPSOParticle::~SPSOParticle() {}
    SPSOParticle::SPSOParticle(const SPSOParticle& rhs):Particle(rhs){
        Particle::operator=(rhs);
        m_lbestidx = rhs.m_lbestidx;
        m_seedship = rhs.m_seedship;
        m_speciesship = rhs.m_speciesship;
    }

    SPSOParticle& SPSOParticle::operator=(const SPSOParticle& rhs){
        if (this == &rhs) return *this;
        Particle::operator=(rhs);
        m_lbestidx = rhs.m_lbestidx;
        m_seedship = rhs.m_seedship;
        m_speciesship = rhs.m_speciesship;
        return *this;
    }
}