#ifndef MP_POP_GL_COM_H
#define MP_POP_GL_COM_H
#include "../gl_pop_seq.h"
#include "../../../../../utility/functional.h"
#include "../../../template/combination/multi_population/mp_com_pop.h"
#include "../../../template/combination/sequence/sequence_algorithm.h"
namespace ofec {
	template<typename TAdaptor>
	class PopGLMpCom : public MP_COM_Pop<PopGLSeq<TAdaptor>> {
	public:
		using AdaptorType = typename TAdaptor;
		using OpType = typename TAdaptor::OpType;
		using SolutionType = typename TAdaptor::SolutionType;
		using InterpreterType = typename TAdaptor::InterpreterType;

	protected:
		int m_maxIndiSize = 50;
	public:
		void setMaxIndiSize(int size) {
			m_maxIndiSize = size;
		}
		void setInnerRadius(Real radius) {
			PopGLSeq<TAdaptor>::setRadius(radius);
			m_innerRadius = radius;
		}
		Real getInnerRadius() {
			return m_innerRadius;
		}
		virtual bool betterThan(Problem *pro,Algorithm *alg,const PopGLMpCom& pop) {
			auto& adaptor(dynamic_cast<AdaptorType&>(*m_adaptor));
			auto& op = GET_ASeq(alg).getOp<OpType>();
			return op.better(pro,
				dynamic_cast<AdaptorType&>(*m_adaptor).getCenter(),
				dynamic_cast<const AdaptorType&>(*pop.m_adaptor).getCenter());
		}

		virtual void addIndisAround(Problem *pro, Algorithm *alg) ;
		virtual void initPopAround(Problem *pro, Algorithm *alg, 
			const SolutionType& center, Real radius, size_t* popSize = nullptr) override;
		virtual void mergePop(Problem *pro, Algorithm *alg, PopGLMpCom& pop) ;
		virtual void mergePop(Problem *pro, Algorithm *alg, std::vector<std::unique_ptr<SolutionType>>& pop);
		virtual Real distance(Problem *pro, Algorithm *alg, const PopGLMpCom& pop);
		virtual Real avgRadius(Problem *pro);
	};

	template<typename TAdaptor>
	void PopGLMpCom<TAdaptor>::addIndisAround(Problem *pro, Algorithm *alg) {
		int originSize = this->size();
		resize(m_maxIndiSize, pro );
		//m_individuals.resize(m_maxIndiSize);
		auto& adaptor(dynamic_cast<AdaptorType&>(*m_adaptor));
		auto& center(adaptor.getCenter());
		auto& op = GET_ASeq(alg).getOp<OpType>();
		auto& interperter = GET_ASeq(alg).getInterpreter<InterpreterType>();
		Real radius = getRadius();
		auto rnd = alg.idRandom();
		for (int id(originSize); id < m_maxIndiSize; ++id) {
			//m_individuals[id].reset(new SolutionType());
			m_individuals[id]->reset();
			m_individuals[id]->initialize(id, pro, rnd);
			op.learn_from_other(*m_individuals[id], rnd, pro, interperter, center, m_innerRadius);
		}
		for (int id(0); id < m_maxIndiSize; ++id) {
			m_individuals[id]->setId(id);
		}
		for (auto& it : m_individuals) {
			interperter.stepFinal(pro, *it);
			op.evaluate(pro, alg, rnd, *it);
		}
		initializeMemory(pro, alg);
		m_offspring.clear();
		for (int i = 0; i < this->size(); i++) {
			m_offspring.push_back(*this->m_individuals[i]);
		}
	}

	template<typename TAdaptor>
	void PopGLMpCom<TAdaptor>::initPopAround(Problem *pro, Algorithm *alg, const SolutionType& center, Real radius, size_t* popSize){
		
		if (popSize != nullptr) {
			assign(*popSize, pro);
			resize(*popSize, pro);
		}
		for (int id(0); id < m_individuals.size(); ++id) {
			m_individuals[id]->setId(id);
		}
		m_his.clear();
		setRadius(radius);
		auto& op =GET_ASeq(alg).getOp<OpType>();
		auto& interperter = GET_ASeq(alg).getInterpreter<InterpreterType>();
		*m_individuals.front() = center;
		auto rnd = alg->idRandom();
		for (int id(1); id < m_individuals.size(); ++id) {
			m_individuals[id]->reset();
			m_individuals[id]->initialize(id, pro, rnd);
			op.learn_from_other(*m_individuals[id], rnd, pro, interperter, center,radius);
		}
		for (auto& it : m_individuals) {
			interperter.stepFinal(pro, *it);
			op.evaluate(pro, alg, rnd, *it);
		}
		initializeMemory(pro, alg);
		m_offspring.clear();
		for (int i = 0; i < this->size(); i++) {
			m_offspring.push_back(*this->m_individuals[i]);
		}

	}


	template<typename TAdaptor>
	void PopGLMpCom<TAdaptor>::mergePop(Problem *pro, Algorithm *alg, PopGLMpCom& pop) {
	//	auto& op = GET_ASeq(alg).getOp<OpType>();
	//	auto& interperter = GET_ASeq(alg).getInterpreter<InterpreterType>();

	//	int curbestId = calIdx<std::unique_ptr<SolutionType>>(m_individuals, [&](const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b) {
	//		return op.better(pro, *a, *b);
	//	});
	//	int otherbestId = calIdx<std::unique_ptr<SolutionType>>(pop.m_individuals, [&](const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b) {
	//		return op.better(pro, *a, *b);
	//	});
	//	auto& adaptor(dynamic_cast<AdaptorType&>(*m_adaptor));
	//	adaptor.setCenter(*m_individuals.front());
	//	if (op.better(pro, *m_individuals[curbestId], *pop.m_individuals[otherbestId])) {
	//		adaptor.setCenter(*m_individuals[curbestId]);
	//	}
	//	else {
	//		adaptor.setCenter(*pop.m_individuals[otherbestId]);
	//	}
	//	auto& center(adaptor.getCenter());
	//	int originSize(m_individuals.size());
	//	std::vector<std::unique_ptr<SolutionType>> newIndis;
	//	for (auto& it : pop.m_individuals) {
	//		if (center.variableDistance(*it, pro) <= m_innerRadius) {
	//			newIndis.emplace_back(std::move(it));
	//		}
	//	}
	//	for (auto& it : m_individuals) {
	//		if (center.variableDistance(*it, pro) <= m_innerRadius) {
	//			newIndis.emplace_back(std::move(it));
	//		}
	//	}
	//	std::sort(newIndis.begin(), newIndis.end(), [&](
	//		const std::unique_ptr<SolutionType>& a,
	//		const std::unique_ptr<SolutionType>& b) {
	//		return op.better(pro, *a, *b);
	//	});

	//	if(newIndis.size()> m_maxIndiSize)
	//	newIndis.resize(m_maxIndiSize);
	//	swap(m_individuals, newIndis);
	//	for (int idx(0); idx < m_individuals.size(); ++idx) {
	//		m_individuals[idx]->setId(idx);
	////		m_individuals[idx]->setActive(true);
	//	}
		mergePop(pro, alg, pop.m_individuals);
		pop.clear();
		pop.setActive(false);

	}
	template<typename TAdaptor>
	void PopGLMpCom<TAdaptor>::mergePop(Problem *pro, Algorithm *alg, std::vector<std::unique_ptr<SolutionType>>& pop) {
		if (pop.empty())return;
		auto& op = GET_ASeq(alg).getOp<OpType>();
		auto& interperter = GET_ASeq(alg).getInterpreter<InterpreterType>();

		int curbestId = calIdx<std::unique_ptr<SolutionType>>(m_individuals, [&](const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b) {
			return op.better(pro, *a, *b);
		});
		int otherbestId = calIdx<std::unique_ptr<SolutionType>>(pop, [&](const std::unique_ptr<SolutionType>& a, const std::unique_ptr<SolutionType>& b) {
			return op.better(pro, *a, *b);
		});
		auto& adaptor(dynamic_cast<AdaptorType&>(*m_adaptor));
		if (curbestId==-1||!op.better(pro, *m_individuals[curbestId], *pop[otherbestId])) {
			//adaptor.setCenter(*pop[otherbestId]);
			adaptor.getCenter() = *pop[otherbestId];

		}
		else {
			adaptor.getCenter() = *pop[curbestId];
			//adaptor.setCenter(*m_individuals[curbestId]);
		}
		auto& center(adaptor.getCenter());
		int originSize(m_individuals.size());
		std::vector<std::unique_ptr<SolutionType>> newIndis;
		
		for (auto& it : pop) {
			if (center.variableDistance(*it, pro) <= m_innerRadius) {
				bool isUnique(true);
				for (auto& it2 : newIndis) {
					if (it->same(*it2, pro)) {
						isUnique = false;
						break;
					}
				}
				if(isUnique)newIndis.emplace_back(std::move(it));
			}
		}
		for (auto& it : m_individuals) {
			if (center.variableDistance(*it, pro) <= m_innerRadius) {
				bool isUnique(true);
				for (auto& it2 : newIndis) {
					if (it->same(*it2, pro)) {
						isUnique = false;
						break;
					}
				}
				if (isUnique)newIndis.emplace_back(std::move(it));
			}
		}
		std::sort(newIndis.begin(), newIndis.end(), [&](
			const std::unique_ptr<SolutionType>& a,
			const std::unique_ptr<SolutionType>& b) {
			return op.better(pro, *a, *b);
		});
		if(newIndis.size()>= m_maxIndiSize)
		newIndis.resize(m_maxIndiSize);
		swap(m_individuals, newIndis);
		for (int idx(0); idx < m_individuals.size(); ++idx) {
			m_individuals[idx]->setId(idx);
			//		m_individuals[idx]->setActive(true);
		}
		resize(m_individuals.size(), pro);

		initializeMemory(pro, alg);
		m_offspring.clear();
		for (int i = 0; i < this->size(); i++) {
			m_offspring.push_back(*this->m_individuals[i]);
		}

	}



	template<typename TAdaptor>
	Real PopGLMpCom<TAdaptor>::distance(Problem *pro, Algorithm *alg, const PopGLMpCom& pop) {
		return dynamic_cast<AdaptorType&>(*m_adaptor).getCenter().variableDistance(dynamic_cast<AdaptorType&>(*m_adaptor).getCenter(),pro);
	//	return 0;
	}

	template<typename TAdaptor>
	Real PopGLMpCom<TAdaptor>::avgRadius(Problem *pro) {
		Real avgDis(0);
		auto& center(dynamic_cast<AdaptorType&>(*m_adaptor).getCenter());
		for (auto& it : m_individuals) {
			avgDis += center.variableDistance(*it, pro);
		}
		return avgDis /= m_individuals.size();
	}
}
#endif