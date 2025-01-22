#ifndef OFEC_FAMF_POP_H
#define OFEC_FAMF_POP_H

#include "hibernated_zone.h"

namespace ofec {
	template<typename TPopulation>
	class FAMFPop : public TPopulation {
	public:
		using IndividualType = typename TPopulation::IndividualType;

	private:
		bool f_hiber = false;
		bool stop_flag = 1;

	protected:
		int m_countRadiusDecreace;
		bool m_stagnantFlag;

		double m_initialRadius;
		double m_curRadius;

		Real m_convergThreshold = 0.0001;
		Real m_convFactor = 0.005;

		IndividualType m_center;

	public:
		FAMFPop(int Popsize) :TPopulation(Popsize), m_countRadiusDecreace(0), m_stagnantFlag(false) {
			estimateInitialRadius();
		};
		//FAMFPop(int Popsize,Real c1, Real c2, Real w) :TPopulation(Popsize,c1,c2,w), m_countRadiusDecreace(0), m_stagnantFlag(false) {
		//	estimateInitialRadius();
		//};
		bool isStagnant(double avgRadius);//*
		void estimateInitialRadius();
		virtual ~FAMFPop() {}
		void computeCenter(Problem *pro);
		void updateCurRadius(Problem *pro);
		bool updateHiberFlag();//*
		bool isHibernating() { return f_hiber; }//*
		void merge_pop(FAMFPop& pop, int max_size);
		int best_idx(Problem *pro);
		//void deratingFitness(unique_ptr<TPopulation>& single){};//*
		bool isConverged() { return m_curRadius <= m_convergThreshold; }
		bool isStagnant();//*
		void setHiberFlag(bool flag);//*
		void wakeUp();
		double getInitialraidus() { return this->m_initialRadius; }
		double getCurrentraidus() { return this->m_curRadius; }
		IndividualType& get_center() { return m_center; }
		// r1 domiante r2
		Real update_radius(Real  r1, Real r2) {
			if (r1 > r2) {
				r1 -= (r1 - r2) * r2 / r1;
			}
			else {
				r1 += (r2 - r1) * r1 / r2;
			}

			return r1;
		}

		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override{
			Real before_radius(m_curRadius);
			int rf = TPopulation::evolve(pro,alg,rnd);
			if (rf != kNormalEval)return rf;
			updateCurRadius(pro);
			if (m_curRadius / (before_radius + 0.00001) < 0.9) m_countRadiusDecreace = 0;
			else m_countRadiusDecreace++;
			this->updateBest(pro);
			std::unique_ptr<IndividualType> new_indiv = std::make_unique<IndividualType>(pro->numberObjectives(), pro->numberConstraints(), CAST_CONOP(pro)->numberVariables());
			*new_indiv = *(this->m_individuals[best_idx(pro)]);
			if(this->m_stagnantFlag)
			{
				new_indiv->cauchyMove(pro, rnd);
				new_indiv->evaluate(pro, alg);
				*(this->m_individuals[best_idx(pro)]) = *new_indiv;
			}
			else {
				new_indiv->brwonianMove(pro, rnd, m_curRadius);
				new_indiv->evaluate(pro, alg);
				if (new_indiv->dominate(*(this->m_individuals[best_idx(pro)]), pro)) {
					*(this->m_individuals[best_idx(pro)]) = *new_indiv;
				}
			}
			this->updateBest(*new_indiv, pro);
			if (rf != kNormalEval) return rf;
			return rf;
		}

		void degrade_explored_areas(HibernatedZone& h) {
			for (auto& it : m_best) {
				h.derateFitness(*it);
			}
			for (auto& it : m_individuals) {
				h.derateFitness(*it);
				// for pso
				//h.derateFitness(*it);
			}
		}
	};


	template<typename TPopulation>
	bool FAMFPop<TPopulation>::isStagnant(double avgRadius)
	{
		if (m_countRadiusDecreace >= 1. * this->size() && this->m_curRadius >= avgRadius && this->m_curRadius > m_convFactor * CONTINUOUS_CAST->domain_size() || m_countRadiusDecreace >= 10 * this->size())
			m_stagnantFlag = true;
		else
			m_stagnantFlag = false;
		return m_stagnantFlag;
	}

	template<typename TPopulation>
	void FAMFPop<TPopulation>::estimateInitialRadius() {
		if (this->size() <= 1) {
			this->m_initialRadius = 0;
			return;
		}
		this->updateCurRadius();
		this->m_initialRadius = this->m_curRadius;
	}

	template<typename TPopulation>
	void FAMFPop<TPopulation>::merge_pop(FAMFPop& pop, int max_size){
		int better(0), worse(0);
		for (auto& it : m_best) {
			for (auto& it2 : pop.m_best) {

				auto r = objective_compare(it->objective(), it2->objective(), global::ms_global->m_problem->opt_mode());
				if (r == kDominant) {
					better++;
				}
				else if (r == kDominated) {
					worse++;
				}
			}
		}
		if (better >= worse) {
			m_initialRadius = update_radius(m_initialRadius, pop.m_initialRadius);
		}
		else {
			m_initialRadius = update_radius(pop.m_initialRadius, m_initialRadius);
		}
		for (auto& it : pop.best()) {
			updateBest(*it, pro);
		}
		size_t cur_size(m_individuals.size());
		m_individuals.resize(m_individuals.size() + pop.size());
		for (size_t i = cur_size; i < m_individuals.size(); ++i) {
			swap(m_individuals[i], pop.m_individuals[i - cur_size]);
		}
		if (m_individuals.size() > max_size) {
			sort();
			std::sort(m_individuals.begin(), m_individuals.end(), [](std::unique_ptr<IndividualType>& a, std::unique_ptr<IndividualType>& b) {
				return a->fitness() < b->fitness();
				});
			m_individuals.resize(max_size);
		}
	}

	template<typename TPopulation>
	inline int FAMFPop<TPopulation>::best_idx(Problem *pro){
		int l = 0;
		for (int i = 0; i < this->m_individuals.size(); ++i) {
			if (this->m_individuals[i].dominate(this->m_individuals[l], pro)) {
				l = i;
			}
		}
		return l;
	}

	//template<typename TPopulation>
	//void FAMFPop<TPopulation>::deratingFitness(unique_ptr<TPopulation>& single){

	//}

	template<typename TPopulation>
	inline void FAMFPop<TPopulation>::computeCenter(Problem *pro) {
		if (this->size() < 1) return;
		if (pro->hasTag(ProblemTag::kConOP)) {
			for (size_t i = 0; i < CAST_CONOP(pro)->numberVariables(); i++) {
				double x = 0.;
				for (int j = 0; j < this->size(); j++) {
					auto chr = this->m_individuals[j]->variable();
					x = chr[i] + x;
				}
				m_center.variable()[i] = x / this->size();
			}
			//m_center.evaluate(true, caller::Algorithm);
		}
	}

	template<typename TPopulation>
	inline bool FAMFPop<TPopulation>::isStagnant()
	{
		return false;
	}

	template<typename TPopulation>
	inline void FAMFPop<TPopulation>::setHiberFlag(bool flag){
		f_hiber = flag;
	}

	template<typename TPopulation>
	void FAMFPop<TPopulation>::wakeUp(){
		f_hiber = false;
	}

	template<typename TPopulation>
	void FAMFPop<TPopulation>::updateCurRadius(Problem *pro) {
		m_curRadius = 0;
		computeCenter(pro);
		if (this->size() < 2) return;
		for (int j = 0; j < this->size(); j++) m_curRadius += this->m_individuals[j]->variableDistance(m_center,pro);
		m_curRadius = m_curRadius / this->size();
	}

	template<typename TPopulation>
	bool FAMFPop<TPopulation>::updateHiberFlag()
	{
		return false;
	}
}

//namespace OFEC {
//
//	class FAMFPop_DE : public DE::population<OFEC::FAMFSolution_DE> {
//	public:
//		//using IndividualType = OFEC::FAMFSolution_DE;
//		FAMFPop_DE(size_t size_pop) : DE::population<FAMFSolution_DE>(size_pop) {
//			m_F = 0.5, m_CR = 0.6;
//			m_mutation_strategy = DE::mutation_strategy::best_2;
//			CONTINUOUS_CAST->set_eval_monitor_flag(true);
//		}
//	};
//
//	class FAMFPop_PSO : public OFEC::swarm<OFEC::FAMFSolution_PSO> {
//	public:
//		//using IndividualType = OFEC::FAMFSolution_DE;
//		FAMFPop_PSO(size_t size_pop) : OFEC::swarm<FAMFSolution_PSO>(size_pop) {
//			m_C1 = 1.4, m_C2 = 1.4;						//accelerators
//			m_weight = 1;								// inertia weight	
//			CONTINUOUS_CAST->set_eval_monitor_flag(true);
//		}
//		FAMFPop_PSO(size_t size_pop, Real c1, Real c2, Real w) : OFEC::swarm<FAMFSolution_PSO>(size_pop) {
//			m_C1 = c1, m_C2 = c2;						//accelerators
//			m_weight = w;								// inertia weight	
//			CONTINUOUS_CAST->set_eval_monitor_flag(true);
//		}
//		int evolve() override {
//			int rf;
//			if (m_iteration == 0)
//				update_best();
//			rf = swarm::evolve();
//			return rf;
//		}
//		void set_neighborhood() override {}
//		solution<>& neighborhood_best(int) override { return *(m_best[0]); }
//	};
//
//}


#endif