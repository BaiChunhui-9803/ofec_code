#ifndef OFEC_BOUNDARY_MAP_H
#define OFEC_BOUNDARY_MAP_H

#include "../../core/problem/problem.h"
#include "../../core/problem/continuous/continuous.h"
#include "../../utility/functional.h"


namespace ofec {

#define CAST_PROBASE(pro) dynamic_cast<ProblemBoundaryBase*>(pro)

	class ProblemBoundaryBase : virtual public Continuous {
	//protected: ProblemBoundaryBase() {
	//	addInputParameters();
	//}private:


	//OFEC_CONCRETE_INSTANCE(ProblemBoundaryBase)
	protected:



		virtual void updateDomain(Environment* env) {
			m_origin_domain = m_domain;
			setDomain(-100, 100);
			m_nor_domain = m_domain;
			auto optimaData = dynamic_cast<Optima<>*>(m_optima.get());;

			for (size_t idx(0); idx < optimaData->numberSolutions(); ++idx) {
				mapFromB1ToB2(optimaData->solution(idx).variable().data(), m_origin_domain, m_nor_domain);
			}
		}
		//virtual void evaluateObjective(Real* vars, std::vector<Real>& objs) {
		//	mapFromB1ToB2(vars, m_domain, m_origin_domain);
		//	//	TPro::evaluateObjective(vars, objs);
		//}
		//virtual void evaluateObjectiveAndConstraint(Real* vars, std::vector<Real>& objs, std::vector<Real>& cons) {
		//	mapFromB1ToB2(vars, m_domain, m_origin_domain);
		//}

		void mapFromB1ToB2(Real* vars, const Domain<Real>& b1, const Domain<Real>& b2) const{
			for (int idx(0); idx < m_number_variables; ++idx) {
				vars[idx] = mapRealInside(vars[idx], b1[idx].limit.first, b1[idx].limit.second, b2[idx].limit.first, b2[idx].limit.second);
			}
		}

		bool isTabuInside(const VariableBase& vars) {
			for (int idx(0); idx < m_tabuSolutions.size(); ++idx) {
				double dis = variableDistance(vars, m_tabuSolutions[idx]->variableBase());
				if (dis <= m_radius[idx])return true;
			}
			return false;
		}

	//	void addInputParameters();

		void inputFile();

	public:
		void setTabuInfo(const std::vector<ofec::SolutionBase*>& sols, std::vector<double>& radius) {
			m_tabuSolutions.clear();
			for (auto& it : sols) {
				m_tabuSolutions.emplace_back(Problem::createSolution(*it));
			}
			m_radius = radius;
		}



	protected:
		Domain<Real> m_origin_domain;		// search domain
		Domain<Real> m_nor_domain;		// search domain

	protected:
		bool m_tabuSetup = false;
		
		std::vector<double> m_radius;
		std::vector<std::shared_ptr<ofec::SolutionBase>> m_tabuSolutions;
	};
	


	template<typename TPro>
	class ProblemBoundary : public TPro,public ProblemBoundaryBase {

		OFEC_CONCRETE_INSTANCE(ProblemBoundary)

	//public: 
	//	static ProblemBoundary* create() {
	//		return new ProblemBoundary();
	//	}
	//
	//protected: 
	//	ProblemBoundary () {
	//	addInputParameters();
	//	};
	protected:


	protected:
		void addInputParameters() {}
		virtual void initializeAfter_(Environment* env) {
			TPro::initializeAfter_(env);
			updateDomain(env);

			//if (m_tabuSetup) {
			//	inputFile();
			//}

		}
		virtual void evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons)const override {
			const VariableVector<Real>& x = dynamic_cast<const VariableType&>(vars);


		/*	if (isTabuInside(vars)) {

				std::cout << "tabu inside" << std::endl;
				objs = m_default_objective;
				cons = m_default_contrait;
			}
			else */
			
			{


				std::vector<Real> x_(x.begin(), x.end()); //for parallel running
				mapFromB1ToB2(x_.data(), ProblemBoundaryBase::m_nor_domain, ProblemBoundaryBase::m_origin_domain);


				if (boundaryViolated(vars)) {
					objs = m_default_objective;
					cons = m_default_contrait;
				}
				else {
					if (cons.empty()) {
						evaluateObjective(x_.data(), objs);
					}
					else {
						evaluateObjectiveAndConstraint(x_.data(), objs, cons);
					}
				}
			}
		}
		

		
	};




}

#endif // !OFEC_BOUNDARY_MAP_H