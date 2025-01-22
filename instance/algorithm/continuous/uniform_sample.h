/********* Begin Register Information **********
{
	"name": "uniform sample",
	"identifier": "UniformSample",
	"problem tags": [ "ConOP", "MOP", "SOP", "MMOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_UNIFORM_SAMPLE_H
#define OFEC_UNIFORM_SAMPLE_H

#include "../../../core/algorithm/algorithm.h"
#include "../../../core/problem/solution.h"

namespace ofec {
	class UniformSample : public Algorithm {
	public:
		virtual void record() override;

	protected:
		virtual void run_() override;

		std::list<Solution<>> m_samples;
	};
}

#endif // !OFEC_UNIFORM_SAMPLE_H
