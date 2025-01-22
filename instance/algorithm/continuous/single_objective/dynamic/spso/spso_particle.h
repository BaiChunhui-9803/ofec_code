#ifndef SPSO_PARTICLE_H
#define SPSO_PARTICLE_H

#include "../../../../template/classic/pso/particle.h"
#include <vector>

namespace ofec {
    class SPSOParticle : public Particle
    {
    public:
        SPSOParticle(size_t num_obj, size_t num_con, size_t size_var);
        ~SPSOParticle();

        SPSOParticle(const SPSOParticle& other);
        SPSOParticle& operator = (const SPSOParticle& other);

        bool isSeed() { return m_seedship; }
        void setSeed(bool seed) { m_seedship = seed; }
        bool isSpecies() { return m_speciesship; }
        void setSpecies(bool species) { m_speciesship = species; }
    private:
        int m_lbestidx;
        bool m_seedship = false;
        bool m_speciesship = false;
    };
}

#endif // SPSO_PARTICLE_H
