#include "dfdl_swarm.h"

#include <memory>

namespace ofec {
    DFDLSwarm::DFDLSwarm(size_t pop, Problem *pro, PopulationType pop_type): Swarm<DFDLParticle>(pop, pro) {
        m_pop_type = pop_type;
        m_max_sub_size = 5;
        m_center = std::make_unique<DFDLParticle>(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
    }

    void DFDLSwarm::setNeighborhood(Random *rnd) {
        m_link.assign(this->m_individuals.size(), std::vector<bool>(this->m_individuals.size(), true));
    }

    void DFDLSwarm::checkOverCrowd(Problem *pro, Random *rnd) {
        if (size() <= m_max_sub_size) return;
        this->sortPbest(pro);
        for(unsigned int i(m_max_sub_size); i < size();i++){
            m_individuals.erase(m_individuals.begin() + i);
            i--;
        }
        setNeighborhood(rnd);
    }

    void DFDLSwarm::checkConverged(Problem *pro) {
        updateCurVelocity(pro);
        updateCurRadius(pro);
        if(m_curRadius < m_seed_radius) {
            setSeedFlag(true);
        }
        if(m_curRadius < m_converge_radius) {
            setConvergeFlag(true);
            m_stagnanted_count = 0;
        }
        else {
            setConvergeFlag(false);
        }
    }

    void DFDLSwarm::checkStagnanted(Real avg_radius) {
        if(m_convergeFlag) {
            m_stagnanted_count = 0;
            m_stagnantedFlag = false;
        }
        else {
            if (m_stagnanted_count > this->size() && this->m_curRadius > avg_radius * 10 || m_stagnanted_count > this->size() * 10)
                m_stagnantedFlag = true;
            else
                m_stagnantedFlag = false;
        }
    }

    void DFDLSwarm::sortPbest(Problem *pro) {
        DFDLParticle temp(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
        for (int i = 0; i < this->size() - 1; i++) {
            for (int j = this->size() - 1; j > 0; j--) {
                if (m_individuals[j]->pbest().dominate(m_individuals[j - 1]->pbest(), pro)) {
                    temp = (*m_individuals[j]);
                    (*m_individuals[j]) = (*m_individuals[j - 1]);
                    (*m_individuals[j - 1]) = temp;
                }
            }
        }
    }

    void DFDLSwarm::merge(DFDLSwarm &pop, Problem *pro, Random *rnd) {
        const size_t pre_size = this->size();
        this->m_individuals.resize(pre_size + pop.size());
        setNeighborhood(rnd);
        for (size_t i = 0; i < pop.size(); i++) {
            std::swap(this->m_individuals[pre_size + i], pop.m_individuals[i]);
        }
    }

    size_t DFDLSwarm::bestIndex(Problem *pro) {
        int l = 0;
        for (int i = 0; i < this->m_individuals.size(); ++i) {
            if (this->m_individuals[i]->pbest().dominate(this->m_individuals[l]->pbest(), pro)) {
                l = i;
            }
        }
        return l;
    }

    void DFDLSwarm::initializeRadius() {
        m_initRadius = m_curRadius;
    }

    void DFDLSwarm::initializeParameters(Problem *pro, Random *rnd) {
        m_accelerator1 = 1.496; // accelerators
        m_accelerator2 = 1.496;
        m_weight = 0.7298;	 // inertia weight
        initVelocityMax(pro, rnd);
        initPbest(pro);
        setNeighborhood(rnd);
        updateCurRadius(pro);
        initializeRadius();
    }

    void DFDLSwarm::initializeNichingPops(Problem *pro, Random *rnd) {
        m_accelerator1 = 1.496; // accelerators
        m_accelerator2 = 1.496;
        m_weight = 0.7298;	 // inertia weight
        initVelocityMax(pro, rnd);
        setNeighborhood(rnd);
        updateCurRadius(pro);
        initializeRadius();
        for (int j = 0; j < this->size(); j++) {
            this->m_individuals[j]->initVelocity(pro, rnd);
//            this->m_individuals[j]->initNichingVelocity(pro, rnd, m_initRadius);
        }
    }

    void DFDLSwarm::updateCurRadius(Problem *pro) {
        m_curRadius = 0;
        computeCenter(pro);
        if (this->size() < 2) return;
        for (int j = 0; j < this->size(); j++) m_curRadius += this->m_individuals[j]->pbest().variableDistance(*m_center, pro);
        m_curRadius = m_curRadius / this->size();
    }

    void DFDLSwarm::updateCurVelocity(Problem *pro) {
        m_curVelocity = 0;
        if (this->size() < 2) return;
        for (int j = 0; j < this->size(); j++)
            for (int k = 0; k < this->m_individuals[j]->velocity().size(); k++)
                m_curVelocity += abs(this->m_individuals[j]->velocity()[k]);
        m_curVelocity = m_curVelocity/this->size();
    }

    void DFDLSwarm::computeCenter(Problem *pro) {
        if (this->size() < 1) return;
        if (pro->hasTag(ProblemTag::kConOP)) {
            for (size_t i = 0; i < CAST_CONOP(pro)->numberVariables(); i++) {
                double x = 0.;
                for (int j = 0; j < this->size(); j++) {
                    auto chr = this->m_individuals[j]->pbest().variable();
                    x = chr[i] + x;
                }
                m_center->variable()[i] = x / this->size();
            }
            //m_center.evaluate(true, caller::Algorithm);
        }
    }

    int DFDLSwarm::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
        int rf;
        Real before_radius(m_curRadius);
        rf = Swarm<DFDLParticle>::evolve(pro, alg, rnd);
        updateCurRadius(pro);
        updateCurVelocity(pro);
        if(m_curRadius / (before_radius + 0.00001) >= 0.9 && !m_convergeFlag) m_stagnanted_count++;     // 舍弃了一个判定条件 m_curVelocity < 10e-1
        else m_stagnanted_count = 0;
//        std::cout << "------------------------------------- " << std::endl;
//        std::cout << "m_stagnanted_count: " << m_stagnanted_count << std::endl;
        return rf;
    }
    std::ostream & operator<<(std::ostream &out, DFDLSwarm &A) {
        out << "solution: " << std::endl;
        for (auto &ind : A.m_individuals) {
            for (auto &var: ind.variable()) {
                out << var << " ";
            }
            out << ind.objective(0) << " ";
            out << std::endl;
        }
        out << std::endl;
        out << "pbest: " << std::endl;
        for (auto &ind : A.m_individuals) {
            for (auto &var: ind->pbest().variable()) {
                out << var << " ";
            }
            out << ind->pbest().objective(0) << " ";
            out << std::endl;
        }
        return out;
    }
}