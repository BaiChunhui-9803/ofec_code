#ifndef MP_COM_POP_H
#define MP_COM_POP_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/global.h"

namespace ofec {
	template<typename TPopulation>
	class MP_COM_Pop : public TPopulation {
	public:
		using SolutionType = typename TPopulation::SolutionType;
	protected:
		std::list<std::unique_ptr<SolutionType>> m_discard_indis;
	public:
		virtual void addSol(int idx,const SolutionType& indi) {
			if (idx != -1)
				*this->m_individuals[idx] = indi;
		}
		virtual void addSol(const SolutionType& indi) {
			this->m_individuals.emplace_back(new SolutionType(indi));
		}
		virtual void setRadius(Real radius) {}
		std::list<std::unique_ptr<SolutionType>>& getDiscardIndis() {
			return m_discard_indis;
		}
	//	const SolutionType& getCenter()const = 0;
	/*	virtual void mergePop(Problem *pro, Algorithm *alg, MP_COM_Pop& pop) = 0;
		virtual Real distance(Problem *pro, Algorithm *alg, const MP_COM_Pop& pop) = 0;*/
		virtual void initPopAround(Problem *pro,Algorithm *alg,const SolutionType& center, Real radius, size_t* popSize=nullptr) = 0;
	protected:
		Real m_innerRadius = 1.0;
		Real m_outerRadius = 0.0;

	};
}


#endif