#ifndef OFEC_FAMF_Solution_H
#define OFEC_FAMF_Solution_H

namespace ofec {
	template <typename TInd>
	class FAMFSolution : public TInd {
	protected:
		// whether the Solution fall into explored area
		bool m_explored = false;
		std::vector<Real> m_Real_obj;

	public:
		FAMFSolution(size_t size_obj, size_t size_con, size_t size_var) : TInd(size_obj, size_con, size_var) {}
		FAMFSolution(const FAMFSolution& rhs) : TInd(rhs){}
		void brwonianMove(Problem *pro, Random *rnd, double radius) {
			for (size_t j = 0; j < variable().size(); j++) {
				variable()[j] += rnd->uniform.nextNonStd<Real>(0, radius);
			}
			pro->validateSolution(*this, Validation::kSetToBound,rnd);
		}
		void cauchyMove(Problem *pro, Random *rnd, double radius = -1) {
			const auto& range(CAST_CONOP(pro)->boundary());
			for (size_t i = 0; i < variable().size(); i++) {
				if (radius < 0) {
					variable()[i] += rnd->uniform.nextNonStd<Real>(0, (range[i].second - range[i].first) / 2);
				}
				else {
					variable()[i] += rnd->uniform.nextNonStd<Real>(0, radius);
				}
			}
			pro->validateSolution(*this, Validation::kSetToBound, rnd);
		}

	};
}

#endif