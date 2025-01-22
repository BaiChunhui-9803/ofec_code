//
// Created by Dell on 2022/11/17.
//

#ifndef OFEC_DFDL_PARTICLE_H
#define OFEC_DFDL_PARTICLE_H

#include "../../../../template/classic/pso/particle.h"
namespace ofec {
    class DFDLParticle : public Particle {
    public:
        DFDLParticle(size_t num_obj, size_t num_con, size_t size_var): Particle(num_obj, num_con, size_var) {};
        DFDLParticle(const Solution<> & s): Particle(s) {};
        void initialize(int id, Problem *pro, Random *rnd) override;
        void initVelocityMax(Problem *pro, Random *rnd) override;
        void initVelocity(Problem *pro, Random *rnd) override;
        void initNichingVelocity(Problem *pro, Random *rnd, Real radius);
        int cauchyMove(Problem *pro, Algorithm *alg, Random *rnd, double radius=-1);
        void brownianMove(Problem *pro, Random *rnd, double radius);
    };
}


#endif //OFEC_DFDL_PARTICLE_H
