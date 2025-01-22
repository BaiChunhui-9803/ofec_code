#ifndef FLA_RANDOM_WALK_H
#define FLA_RANDOM_WALK_H



#include "../../../core/problem/continuous/continuous.h"


namespace ofec {
	class FLA_RandomWalkBase {
		//
		//struct DirectionBase {
		//};

	protected:
		std::shared_ptr<Problem> m_pro;
		std::shared_ptr<Random> m_rand;
		
		double m_radius = 0;
		int m_stepSize= 0;

		// 
		bool m_direction_flag =false;

		std::vector<std::array<double,2>> m_disturb_radius;
		std::vector<int> m_change_direction;

		std::function<void(ofec::
			
			& sol, Problem *pro)> m_eval_fun;
		
	public:
		

		void setDirectionFlag(bool direction_flag) {
			m_direction_flag = direction_flag;
		}

		void initialize(const std::shared_ptr<Problem>& pro, double seed, double radius,
			int stepSize,
			const std::function<void(SolBase& sol, Problem *pro)>& eval_fun,
			bool flagDirection=false) {
			m_eval_fun = eval_fun;
			m_pro = pro;;
			m_rand.reset(new Random(seed));
			m_radius = radius;
			m_stepSize = stepSize;
			m_direction_flag = flagDirection;

			if (m_pro->hasTag(ProTag::kConOP)) {
				const auto& boundary = CAST_CONOP(m_pro.get())->boundary();
				m_disturb_radius.resize(boundary.size());
				for (int idx(0); idx < m_disturb_radius.size(); ++idx) {
					m_disturb_radius[idx].front() = 0;
					m_disturb_radius[idx].back() = (boundary[idx].second - boundary[idx].first) * radius;
					m_disturb_radius[idx].front() = -m_disturb_radius[idx].back() / 2.0;
					m_disturb_radius[idx].back()= m_disturb_radius[idx].back() / 2.0;

				}
				/*if (m_direction_flag) */{
					m_change_direction.resize(boundary.size());
					for (int idx(0); idx < m_change_direction.size(); ++idx) {
						m_change_direction[idx] = m_rand->uniform.nextNonStd<int>(0, 2);
						if (m_change_direction[idx] == 0) {
							m_change_direction[idx] = -1;
						}
					}
				}
			}
		}

		~FLA_RandomWalkBase() = default;

		virtual void randomWalk(
			const std::unique_ptr<ofec::SolBase>& fromSol,
			std::unique_ptr<ofec::SolBase>& toSol) {
			//using namespace ofec;
			if (m_pro->hasTag(ProTag::kConOP)) {
				
				auto& fromConSol = dynamic_cast<const ofec::Solution<>&>(*fromSol);
				auto& toConSol = dynamic_cast<ofec::Solution<>&>(*toSol);

				const auto& boundary = CAST_CONOP(m_pro.get())->boundary();
				if (m_direction_flag) {
					
					for (int idx(0); idx < boundary.size(); ++idx) {
						toConSol.variable()[idx] = fromConSol.variable()[idx] +
							m_rand->uniform.nextNonStd<double>(0, m_disturb_radius[idx].back()) * m_change_direction[idx];

						if (toConSol.variable()[idx] > boundary[idx].second) {
							toConSol.variable()[idx] = 2 * boundary[idx].second - toConSol.variable()[idx];
							m_change_direction[idx] = -m_change_direction[idx];

						}
						else if (toConSol.variable()[idx] < boundary[idx].first) {
							toConSol.variable()[idx] = 2 * boundary[idx].first - toConSol.variable()[idx];
							m_change_direction[idx] = -m_change_direction[idx];
						}
					}
				}
				else {
					std::pair<double, double> change_range;
					for (int idx(0); idx < boundary.size(); ++idx) {
						change_range.first = std::max<double>(m_disturb_radius[idx].front(), boundary[idx].first- fromConSol.variable()[idx] );
						change_range.second = std::min<double>(m_disturb_radius[idx].back(), boundary[idx].second - fromConSol.variable()[idx]);
						toConSol.variable()[idx] = fromConSol.variable()[idx] +
							m_rand->uniform.nextNonStd<double>(m_disturb_radius[idx].front(), m_disturb_radius[idx].back());
						
					}
				}
			}
			else {
				// to do
			}
		
		};

		virtual void getRandomWalkLink(
			std::vector<std::unique_ptr<ofec::SolBase>>& sols) {
			//using namespace ofec;
			std::unique_ptr<SolBase> curSol(m_pro->createSolution());
			curSol->initialize(m_pro.get(),m_rand.get());
			curSol->evaluate(m_pro.get());
			sols.clear();
			sols.push_back(std::move(curSol));
			for (int idx(0); idx < m_stepSize; ++idx) {
				curSol.reset(m_pro->createSolution(*sols.back()));
				randomWalk(sols.back(), curSol);
				curSol->evaluate(m_pro.get());
				m_eval_fun(*curSol, m_pro.get());
				sols.push_back(std::move(curSol));
			}
		};

	};
}


#endif //  RANDOM_WALK_H
