#ifndef MQSO_PARTICLE_H
#define MQSO_PARTICLE_H

#include "../../../../template/classic/pso/particle.h"
#include <vector>

namespace ofec {
    class MQSOParticle : public Particle
    {
        struct VMax {
            double m_min, m_max;
            VMax() :m_min(0), m_max(0) {}
        };
    public:
        MQSOParticle() = default;
        MQSOParticle(size_t num_obj, size_t num_con, size_t size_var);
        ~MQSOParticle();

        MQSOParticle(const MQSOParticle& other);

        enum class ParticleType { PARTICLE_NEUTRAL, PARTICLE_QUANTUM, PARTICLE_CHARGED };

        //void increaseDimension();
        //void decreaseDimension();
        void nextVelocity(const Solution<>* best, Real w, Real c1, Real c2, Random *rnd) override;
        void move() override;
        void quantumMove(const Solution<>* lbest, Real m_Rcloud, Random *rnd);
        void setType(ParticleType t) {
            m_type = t;
        };

        ParticleType& getType() {
            return m_type;
        }

        std::vector<double>& getRepulse() {
            return mv_repulse;
        }

        int gSign(Real val) {
            return (Real(0) < val) - (val < Real(0));
        }

        void initializeVmax(Problem *pro);
        virtual void initVelocity(Problem *pro, Random *rnd) override;

    private:
        ParticleType m_type;
        std::vector<double> mv_repulse;
        std::vector<VMax> m_vMax;

    };
}

#endif // MQSO_PARTICLE_H
