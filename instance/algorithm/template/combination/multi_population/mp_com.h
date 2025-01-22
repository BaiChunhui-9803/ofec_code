#ifndef MP_COM_H
#define MP_COM_H

#include "../../../../../core/algorithm/multi_population.h"
#include "../../../../../core/global.h"
#include <vector>

namespace ofec {
	template <typename TPopulation>
	class MP_COM : public MultiPopulation<TPopulation> {
		TPopulation m_mainPop;
		Real m_OuterRadius = 0.5;
		
	public:

		

	};
}


#endif