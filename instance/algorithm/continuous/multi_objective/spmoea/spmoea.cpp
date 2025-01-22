#include "spmoea.h"
#include"../../../../../utility/clustering/dbscan.h"

namespace ofec {

	void Subspace::updateSubspaceInfo(Population<Solution<>>& new_pop, Problem* pro) {
		//������ʷ�⼰ǰ������
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < new_pop.size(); ++i) {
			temp_pop.append(new_pop[i]);
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			objs.emplace_back(&temp_pop[i].objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			temp_pop[i].setFitness(rank[i]);
		}
		Population<Solution<>> front_pop;
		std::vector<size_t> add_front_inx;
		std::vector<size_t> add_behind_inx;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			if (temp_pop[i].fitness() == 0) {
				front_pop.append(temp_pop[i]);
				add_front_inx.push_back(i);
			}
			else {
				add_behind_inx.push_back(i);
			}
		}
		//�ӿռ�ǰ�ؼ��������£��ӿռ��������
		if (getSubspaceFrontInx().empty()) {
			for (size_t i = 0; i < front_pop.size(); ++i) {
				getSubspaceFrontInx().push_back(getSubspaceHisSols().size());
				getSubspaceHisSols().emplace_back(new Solution<>(front_pop[i]));
			}
		}
		else {
			auto temp_front_inx = getSubspaceFrontInx();
			std::vector<std::shared_ptr<Solution<>>> temp_front_sols;
			for (size_t j = 0; j < getSubspaceFrontInx().size(); ++j) {
				temp_front_sols.emplace_back(getSubspaceHisSols()[getSubspaceFrontInx()[j]]);
			}
			getSubspaceFrontInx().clear();
			std::vector<size_t> add_his_inx;

			//�ȿ�ÿ�����Ƿ񱻵�ǰǰ��֧��
			std::vector<size_t> temp_his_sols(temp_front_sols.size(), 1);
			std::vector<size_t> temp_pop_sols(add_front_inx.size(), 0);
			for (size_t j = 0; j < add_front_inx.size(); ++j) {
				for (size_t i = 0; i < temp_front_sols.size(); ++i) {
					if (temp_his_sols[i] == 1) {
						Dominance dominanceship = temp_pop[add_front_inx[j]].compare(*temp_front_sols[i], pro->optimizeMode());
						if (dominanceship == Dominance::kNonDominated) {
							if (i == temp_his_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else if (dominanceship == Dominance::kDominant) {
							temp_his_sols[i] = 0;
							if (i == temp_his_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else if (dominanceship == Dominance::kEqual) {
							//temp_his_sols[i] = 0;
							if (i == temp_his_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else {
							add_his_inx.push_back(j);
							break;
						}
					}
					else {
						if (i == temp_his_sols.size() - 1) {
							temp_pop_sols[j] = 1;
						}
					}
				}
			}

			//�����ӿռ�ǰ�ؽ�����
			size_t count = 0;
			for (size_t i = 0; i < temp_pop_sols.size(); ++i) {
				if (temp_pop_sols[i] == 1) {
					count++;
					getSubspaceFrontInx().push_back(getSubspaceHisSols().size());
					getSubspaceHisSols().emplace_back(new Solution<>(temp_pop[add_front_inx[i]]));
					//getMO_HLC().getSubspaceInfo(idx).m_history_inds.emplace_back(new Solution<>(temp_pop[add_front_inx[i]]));
				}
			}
			for (size_t i = 0; i < temp_his_sols.size(); ++i) {
				if (temp_his_sols[i] == 1) {
					getSubspaceFrontInx().push_back(temp_front_inx[i]);
				}
				else {
					getSubspaceHisSols()[temp_front_inx[i]]->setType(-1);
				}
			}
			//���뱻��̭�Ľ�����ʷ��
			for (size_t i = 0; i < add_his_inx.size(); ++i) {
				count++;
				getSubspaceHisSols().emplace_back(new Solution<>(temp_pop[add_front_inx[add_his_inx[i]]]));
			}
			if (count != add_front_inx.size()) {
				size_t a = 1;
			}
		}
	}

	void SPMOEA::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		if (m_pop_size < 4)
			throw MyExcept("Population size must be over 4 for reproduction.");

		m_num_region_var = v.get<int>("number of subspaces");
		m_num_region_obj = v.get<int>("number of obj subspaces");

		m_mo_hlc.reset(new MO_HLC(CAST_CONOP(m_problem.get())->numberVariables(), CAST_CONOP(m_problem.get())->numberObjectives()));
		m_multi_pop.reset(new MultiPopulation<SPMOEA_pop>());
		m_exploit_box.reset(new Subspace(CAST_CONOP(m_problem.get())->boundary()));

		initiVarSpace(m_problem.get());
		
		auto pro_name = CAST_CONOP(m_problem.get())->name();
		if (pro_name == "MOP_TEST3") {
			auto& v = *m_problem->getParam();
			size_t peak_num= v.get<int>("numPeak");
			m_visit_count = std::vector<size_t>(peak_num);
		}

		std::vector<Real> max_ref_point;
		for (size_t i = 0; i < m_problem->numberObjectives(); ++i) {
			Real temp = 0.;
			for (size_t j = 0; j < m_problem->optimaBase()->numberObjectives(); ++j) {
				if (m_problem->optimaBase()->objective(j)[i] > temp) {
					temp = m_problem->optimaBase()->objective(j)[i];
		        }
	        }
			max_ref_point.push_back(temp);
        }
		setMaxRefPoint(max_ref_point);
	}

	void SPMOEA::run_() {
		initPop(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
		updateBuffer();
#endif
		while (!terminating()) {
			evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void SPMOEA::record() {
		//std::vector<Real> entry;
		//entry.push_back(m_evaluations);
		//entry.push_back(m_R1);
		//entry.push_back(m_R2);
		//entry.push_back(m_R3);
		////Real IGD = m_problem->optima().invertGenDist(*m_pop);
		//entry.push_back(m_IGD.back());
		//dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);

		std::vector<Real> entry;
		entry.push_back(m_evaluations);

		//// ����HV
		//std::vector<std::vector<Real>> ref_optima;
		//for (size_t i = 0; i < 1001; ++i) {
		//	std::vector<Real> temp;
		//	temp.push_back(0.001*i);
		//	temp.push_back(1./(1+std::pow(std::tan(OFEC_PI/2*0.001*i),3)));
		//	ref_optima.emplace_back(temp);
		//}

		//std::vector<std::vector<Real>> test_points;
		//for (size_t i = 0; i < 21; ++i) {
		//	std::vector<Real> temp;
		//	temp.push_back(0.05 * i+0.6);
		//	temp.push_back(0.6+1. / (1 + std::pow(std::tan(OFEC_PI / 2 * 0.05 * i), 3)));
		//	test_points.emplace_back(temp);
		//}
		////����HV
		//Real temp_hv = hypervolumeVector(test_points, ref_optima, m_problem.get(), m_random.get());
		//Real temp_hv2 = hypervolumeVector(ref_optima, ref_optima, m_problem.get(), m_random.get());
		
		//HV
		//�ȵõ�PF�������е����ֵ��Ϊref_point
		//auto& maxRefPoint = getMaxRefPoint();
		////����HV
		//std::vector<std::vector<Real>> test_points;
		//for (size_t i = 0; i < 10; ++i) {
		//	std::vector<Real> temp;
		//	Real xx = (i + 1.) / 11;
		//	temp.push_back(xx);
		//	temp.push_back(1-std::pow(xx,3));
		//	test_points.emplace_back(temp);
		//}
		//Real temp_hv= hypervolumeVector(test_points, maxRefPoint, m_problem.get(), m_random.get());
		
		/*Real HV = hypervolumeVector(temp_objs, maxRefPoint, m_problem.get(), m_random.get());
		entry.push_back(HV);*/
		size_t num_his_front_sols = getHisFrontSols().size();
		entry.push_back(num_his_front_sols);
#ifdef OFEC_DEMO
		std::cout << "�ۻ�ǰ�ؽ����" << num_his_front_sols << std::endl;
		//std::cout << m_evaluations<<"     "<< "IGD:  " << IGD << std::endl;
#else
		Population<Solution<>> temp_pop;
		//std::vector<std::vector<Real>> temp_objs;
		std::vector<std::vector<Real>> temp_decs;
		for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
			temp_pop.append(*m_history_front_sols[i]);
			temp_decs.emplace_back(m_history_front_sols[i]->variable().vect());
		}
		//IGD
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			pop_objs.push_back(temp_pop[i].objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);
		//Real IGD = m_problem->optimaBase()->invertGenDist(temp_pop);
		entry.push_back(temp_IGD);
		//IGDX
		if (m_problem->isOptimumSolutionGiven()) {
			Real distance = 0;
			size_t num_opts = m_problem->optimaBase()->numberSolutions();
			for (size_t i = 0; i < num_opts; ++i) {
				auto& opt_decs = dynamic_cast<const Optima<>&>(*m_problem->optimaBase()).solution(i).variable();
				Real min_d = std::numeric_limits<Real>::max();
				for (int j = 0; j < temp_decs.size(); ++j) {
					Real d = euclideanDistance(opt_decs.begin(), opt_decs.end(), temp_decs[j].begin());
					if (d < min_d) {
						min_d = d;
					}
				}
				distance += min_d;
			}
			Real IGDX = distance / num_opts;
			entry.push_back(IGDX);
		}
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
#endif

	}

#ifdef OFEC_DEMO
	void SPMOEA::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(1);
			for (size_t i = 0; i < m_multi_pop->size(); ++i)
				for (size_t j = 0; j < m_multi_pop->at(i).size(); ++j) {
					m_solution[0].push_back(&m_multi_pop->at(i)[j]);
				}
			/*for (size_t i = 0; i < m_pop->getOffspring().size(); ++i)
				m_solution[1].push_back(&m_pop->getOffspring()[i]);*/
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA::initPop(Problem *pro, Algorithm *alg, Random *rnd) {
		//����ÿ���ӿռ��������ɸ��岢����
		size_t number_objectives = CAST_CONOP(pro)->numberObjectives();
		if (number_objectives == 2) {
			setArchiveNum(100);
		}
		else if (number_objectives == 3) {
			setArchiveNum(300);
		}
		else {
			setArchiveNum(500);
		}
		size_t num = m_mo_hlc->numSubspace();
		size_t pop_size = num;
		if (m_pop_size > num) {
			pop_size = m_pop_size;
		}
		SPMOEA_pop temp_pop(pop_size, pro);
		//temp_pop.initialize(pro,rnd);
		size_t count = 0;
		for (size_t i = 0; i <num ; ++i) {
			//Solution<> ind(CAST_CONOP(pro)->numberObjectives(),0);
			auto bound = m_mo_hlc->subspaceTree().getBox(i);
			for (size_t k = 0; k < bound.size(); ++k) {
				temp_pop[count].variable().vect()[k]=bound[k].first + 1./2 * (bound[k].second - bound[k].first);
			}
			count++;
			//temp_pop.append(ind);
		}
		while (count < pop_size) {
			size_t inx = std::floor(num * rnd->uniform.next());
			//auto space_inx = m_mo_hlc->getCluster(i)[inx];
			//Solution<> ind(CAST_CONOP(pro)->numberObjectives(),0);
			auto bound = m_mo_hlc->subspaceTree().getBox(inx);
			for (size_t k = 0; k < bound.size(); ++k) {
				temp_pop[count].variable().vect()[k]=bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
			}
			count++;
			//temp_pop.append(ind);
		}
		temp_pop.evaluate(pro, alg);
		NDSort(temp_pop);
		//�����ӿռ���Ϣ
		updateVarSpaceInfo(temp_pop, pro, rnd);
		//������ʷǰ��
		updateHistoryInfo(temp_pop, pro);
		m_multi_pop->append(temp_pop);
		//��ʼ���Ӵ�
		for (size_t i = 0; i < m_multi_pop->size(); ++i) {
			for (size_t j = 0;j< m_multi_pop->at(i).size(); ++j) {
				m_multi_pop->at(i).getOffspring()[j] = m_multi_pop->at(i)[j];
				m_multi_pop->at(i).getOffspring()[m_multi_pop->at(i).size() + j] = m_multi_pop->at(i)[j];
			}
			////�����Ⱥ
			//for (size_t j = 0; j < m_multi_pop->at(i).size(); ++j) {
			//	std::cout << m_multi_pop->at(i)[j].variable()[0] << "      " << m_multi_pop->at(i)[j].variable()[1] << std::endl;
			//}
			//std::cout << "******************" << std::endl;
		}
		updateObjRange(temp_pop,pro);

		updateObjSpace();
	}


	SPMOEA_pop::SPMOEA_pop(size_t size_pop, Problem *pro) : Population<Solution<>>(size_pop,pro,CAST_CONOP(pro)->numberVariables()),\
		m_offspring(2 * size_pop, pro, CAST_CONOP(pro)->numberVariables())
	{
	}

	void SPMOEA_pop::initialize(Problem *pro, Random *rnd) {
		/*************************/
		/*       ��Ⱥ��ʼ��      */
		/*************************/
		Population<Solution<>>::initialize(pro, rnd);
	}

	void SPMOEA::initiVarSpace(Problem *pro) {
		/*************************************/
		/*      MO_HLC�����ռ仮�ֳ�ʼ��     */
		/*************************************/
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		std::vector<std::pair<Real, Real>> m_var_boundary;
		for (size_t i = 0; i < var_num; ++i) {
			m_var_boundary.emplace_back(CAST_CONOP(pro)->range(i));
		}
		m_mo_hlc->initialVarSpace(m_var_boundary, m_num_region_var);
		m_pop_var_range = m_var_boundary;
		
		m_var_space_volume = 1.;
		for (size_t i = 0; i <var_num; ++i) {
			Real dim_volume = m_var_boundary[i].second- m_var_boundary[i].first;
			m_var_space_volume *= dim_volume;
		}
	}

	int SPMOEA::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		//visualization options
		/*auto& v = *m_param;
		m_add_neighbor = v.get<bool>(("cluster add neighbors"));
		m_sample_in_basin = v.get<bool>(("sample in basin"));
		m_subspace_select = v.get<bool>(("select by subspace"));
		m_evolve_by_potential = v.get<bool>(("evolve by potential"));*/
		bool m_traditional = true;
		if (m_traditional) {
			//generate solutions by operators, such as GA, DE
			generateOffspring(pro, rnd);
		}
		else {
			//generating solutions by learning methods
			bool m_prediction = false;
			if (!m_prediction) {
				//multi-population method and operators
				//�ӿռ�Ԥ�����������ӿռ�ǰ�ؽ⣬�ӿռ��ڴ������������࣬����Ⱥ����

			}
			else {
				//multi-popution and prediction
				//��������Ԥ�⣬����Ԥ�����
				//����Ⱥ�ݻ���ʷ�켣
			}
		}

		int tag = EvaluationTag::kNormalEval;
		for (size_t i = 0; i < m_multi_pop->size(); ++i) {
			SPMOEA_pop offspring_pop(m_multi_pop->at(i).size(), pro);
			for (size_t j = 0; j < m_multi_pop->at(i).size(); j++) {
				//�����Ӵ�
				tag = m_multi_pop->at(i).getOffspring()[j].evaluate(pro, alg);
				if (tag != EvaluationTag::kNormalEval)
					break;
				offspring_pop[j] = m_multi_pop->at(i).getOffspring()[j];
			}
			//�����ӿռ���Ϣ
			NDSort(offspring_pop);
			updateVarSpaceInfo(offspring_pop, pro, rnd);
			//������ʷ��Ϣ
			updateHistoryInfo(offspring_pop, pro);
			updateObjRange(offspring_pop, pro);
		}

		//����Ⱥ�ڲ�������̭ѡ��
		for (size_t i = 0; i < m_multi_pop->size(); ++i) {
			for (size_t j = 0; j < m_multi_pop->at(i).size(); ++j) {
				m_multi_pop->at(i).getOffspring()[m_multi_pop->at(i).size() + j] = m_multi_pop->at(i)[j];
			}
			//m_multi_pop->at(i).envirSelection(pro, rnd, 0);
			m_multi_pop->at(i).envirSelection(m_multi_pop->at(i), m_multi_pop->at(i).getOffspring(),pro);
			////�����Ⱥ
			//for (size_t j = 0; j < m_multi_pop->at(i).size(); ++j) {
			//	std::cout << m_multi_pop->at(i)[j].variable()[0] << "      " << m_multi_pop->at(i)[j].variable()[1] << std::endl;
			//}
			//std::cout << "************************" << std::endl;
			//m_multi_pop->at(i).iteration()++;
		}

		////��������Ⱥ����λ�ã�����ж�Ϊ��������Ե����ǰ���������򻮷ֳ�������ǡ�exploit���ӿռ�
		//for (size_t i = 0; i < m_multi_pop->size(); ++i) {
		//	if (m_multi_pop->at(i).getPopState() == "active") {
		//		if (subspaceSeperable(i)) {
		//			//����i���ӿռ����»���
		//			splitSubspace(i);
		//			//���µ��ӿռ���������Ⱥ

		//		}
		//	}
		//}

		//���ݸ�������Ⱥ���ӿռ��������������ⲿ�浵,
		//ʹ����ʷ���з�֧������archive
		updateArchive(m_archive_num, pro);
		updateObjSpace();

		recordMetrics(pro, alg);

		return tag;
	}

	void SPMOEA::generatePop(const std::vector<std::vector<size_t>>& space_clusters, Problem *pro, Algorithm *alg, Random *rnd) {
		for (size_t i = 0; i < space_clusters.size(); ++i) {
			size_t num = space_clusters[i].size();
			size_t pop_size = getPopsize();
			/*if (m_pop_size > num) {
				pop_size = m_pop_size;
			}*/
			SPMOEA_pop temp_pop(pop_size, pro);
			size_t count = 0;
			std::vector<std::vector<Real>> sols;
			/*for (size_t j = 0; j < num; ++j) {
				auto bound = m_mo_hlc->subspaceTree().getBox(space_clusters[i][j]);
				for (size_t k = 0; k < bound.size(); ++k) {
					temp_pop[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
				}
				count++;
			}*/
			while (count < pop_size) {
				size_t inx = std::floor(num * rnd->uniform.next());
				auto bound = getMO_HLC().subspaceTree().getBox(space_clusters[i][inx]);
				for (size_t k = 0; k < bound.size(); ++k) {
					temp_pop[count].variable().vect()[k] = bound[k].first + rnd->uniform.next() * (bound[k].second - bound[k].first);
				}
				sols.emplace_back(temp_pop[count].variable().vect());
				count++;
			}
			temp_pop.evaluate(pro, alg);
			NDSort(temp_pop);
			std::vector<Real> middle_point;
			for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
				middle_point.push_back((CAST_CONOP(pro)->range(j).first + CAST_CONOP(pro)->range(j).second) / 2);
			}
			visitCount(sols, middle_point);
			////�����ӿռ���Ϣ
			updateVarSpaceInfo(temp_pop, pro, rnd);
			//������ʷǰ��
			updateHistoryInfo(temp_pop, pro);
			//��ʼ���Ӵ�
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring()[j] = temp_pop[j];
				temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			}
			/*temp_pop.setRate(m_cr, m_mr);
			temp_pop.setEta(m_ceta, m_meta);*/
			getPop().append(temp_pop);
		}
	}

	void SPMOEA::generateOffspring(Problem *pro,Random *rnd) {
		for (size_t i = 0; i < m_multi_pop->size(); ++i) {
			auto sub_bound = m_mo_hlc->subspaceTree().getBox(i);
			for (size_t j = 0; j < m_multi_pop->at(i).size(); j += 2) {
				auto off = sampleByGA(m_multi_pop->at(i), sub_bound,2,pro,rnd);
				//std::vector<size_t> p(2);
				//do {//���ظ�����
				//	do {
				//		p[0] = m_multi_pop->at(i).tournamentSelection(pro, rnd);
				//		p[1] = m_multi_pop->at(i).tournamentSelection(pro, rnd);
				//	} while (p[1] == p[0]);
				//	m_multi_pop->at(i).crossover(p[0], p[1], m_multi_pop->at(i).getOffspring()[j], m_multi_pop->at(i).getOffspring()[j + 1], pro, rnd);
				//	m_multi_pop->at(i).mutate(m_multi_pop->at(i).getOffspring()[j], pro, rnd);
				//	m_multi_pop->at(i).mutate(m_multi_pop->at(i).getOffspring()[j + 1], pro, rnd);
				//	//�Ӵ�Խ�紦��
				//	repairSol(m_multi_pop->at(i).getOffspring()[j].variable().vect(), m_pop_var_range, rnd);
				//	repairSol(m_multi_pop->at(i).getOffspring()[j + 1].variable().vect(), m_pop_var_range, rnd);
				//} while (m_multi_pop->at(i).getOffspring()[j].variable().vect() == m_multi_pop->at(i)[p[0]].variable().vect() || \
				//	m_multi_pop->at(i).getOffspring()[j].variable().vect() == m_multi_pop->at(i)[p[1]].variable().vect() || \
				//	m_multi_pop->at(i).getOffspring()[j + 1].variable().vect() == m_multi_pop->at(i)[p[0]].variable().vect() || \
				//	m_multi_pop->at(i).getOffspring()[j + 1].variable().vect() == m_multi_pop->at(i)[p[1]].variable().vect());

				////�Ӵ�Խ�紦��
				///*m_multi_pop->at(i).crossover(p[0], p[1], m_multi_pop->at(i).getOffspring()[j], m_multi_pop->at(i).getOffspring()[j + 1], pro, rnd);
				//m_multi_pop->at(i).mutate(m_multi_pop->at(i).getOffspring()[j], pro, rnd);
				//m_multi_pop->at(i).mutate(m_multi_pop->at(i).getOffspring()[j + 1], pro, rnd);
				//repairSol(m_multi_pop->at(i).getOffspring()[j].variable().vect(), m_pop_var_range, rnd);
				//repairSol(m_multi_pop->at(i).getOffspring()[j + 1].variable().vect(), m_pop_var_range, rnd);*/
				//m_multi_pop->at(i).getOffspring()[j].setCounter(0);
				//m_multi_pop->at(i).getOffspring()[j + 1].setCounter(0);
			}
			////�����Ⱥ
			//for (size_t j = 0; j < m_multi_pop->at(i).size(); ++j) {
			//	std::cout << m_multi_pop->at(i).getOffspring()[j].variable()[0] << "      " << m_multi_pop->at(i).getOffspring()[j].variable()[1] << std::endl;
			//}
			//std::cout << "*******************" << std::endl;
		}
	}

	std::vector<Real> SPMOEA::vectorOperator(std::vector<Real>& sol1, std::vector<Real>& sol2, Real offset, Problem* pro, Algorithm* alg, Random* rnd) {
		//�м����
		std::vector<Real> off;
		off.resize(sol1.size());
		std::vector<Real> delta;
		for (size_t i = 0; i < sol1.size(); ++i) {
			delta.push_back(sol2[i] - sol1[i]);
		}
		auto ratio = rnd->uniform.next();
		for (size_t i = 0; i < sol1.size(); ++i) {	
			off[i] = sol1[i] + ratio * delta[i];
		}
		//���Ŷ����Ŷ�������ratio���
		Real raodong = ratio > 0.5 ? 1 - ratio : ratio;
		for (size_t i = 0; i < off.size(); ++i) {
			auto rand = rnd->uniform.next();
			off[i] += ((2*rand-1)*(raodong/2)*delta[i]);
		}
		
		//������
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < off.size(); ++i) {
			if (off[i] < bound[i].first) {
				off[i] = 2 * bound[i].first - off[i];
			}
			else if (off[i] > bound[i].second) {
				off[i] = 2 * bound[i].second - off[i];
			}
		}

		size_t M = CAST_CONOP(pro)->numberObjectives();
		size_t c = CAST_CONOP(pro)->numberConstraints();
		Solution<> sol(M, c, bound.size());
		Solution<> ind2(sol);
		ind2.variable().vect() = off;
		ind2.evaluate(pro, alg);
		Solution<> ind1(ind2);
		ind1.variable().vect() = sol1;
		std::vector<std::shared_ptr<Solution<>>> temp_pair;
		temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
		temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
		m_interactive_sol_pair.emplace_back(temp_pair);

		return off;
	}

	std::vector<Real> SPMOEA::vectorPredict(std::vector<Real>& sol1, std::vector<Real>& sol2, std::vector<std::pair<Real, Real>>& bound, Random *rnd) {
		//����Ԥ��
		std::vector<Real> off;
		off.resize(sol1.size());
		std::vector<Real> delta;
		for (size_t i = 0; i < sol1.size(); ++i) {
			delta.push_back(sol2[i] - sol1[i]);
		}
		auto ratio = rnd->uniform.next();
		for (size_t i = 0; i < sol1.size(); ++i) {
			off[i] = sol1[i] + (1. + ratio) * delta[i];
		}
		//������
		for (size_t i = 0; i < off.size(); ++i) {
			if (off[i] < bound[i].first) {
				off[i] = 2 * bound[i].first - off[i];
			}
			else if (off[i] > bound[i].second) {
				off[i] = 2 * bound[i].second - off[i];
			}
		}
		return off;
	}

	//GA���Ӳ���
	std::vector<std::vector<Real>> SPMOEA::sampleByGA(SPMOEA_pop& pop, std::vector<std::pair<Real, Real>>& bound, size_t off_num, Problem *pro, Random *rnd) {
		std::vector<std::vector<Real>> all_off;
		//auto search_bound = CAST_CONOP(pro)->boundary();
		//����SBX��Ⱥ
		PopSBX<> temp_pop(pop.size(), pro);
		for (size_t i = 0; i < pop.size(); ++i) {
			temp_pop[i].variable().vect() = pop[i].variable().vect();
			temp_pop[i].objective() = pop[i].objective();
		}
		temp_pop.setRate(0.9,1./CAST_CONOP(pro)->numberVariables());//���ý�����ʺͱ������
		temp_pop.setEta(1, 1);//����SBX����
		for (size_t j = 0; j < off_num; ++j) {
			std::vector<size_t> p(2);
			std::vector<bool> flag(2, true);
			do {//���ظ�����
				do {
					p[0] = temp_pop.tournamentSelection(pro, rnd);
					p[1] = temp_pop.tournamentSelection(pro, rnd);
				} while (p[1] == p[0]);
				temp_pop.crossover(p[0], p[1], pop.getOffspring()[0], pop.getOffspring()[1], pro, rnd);
				temp_pop.mutate(pop.getOffspring()[0], pro, rnd);
				temp_pop.mutate(pop.getOffspring()[1], pro, rnd);
				//�Ӵ�Խ�紦��
				/*repairSol(pop.getOffspring()[0].variable().vect(), bound, rnd);
				repairSol(pop.getOffspring()[1].variable().vect(), bound, rnd);*/
				flag[0] = flag[1] = true;
				for (size_t k = 0; k < pop.getOffspring()[0].variable().vect().size(); ++k) {
					if (pop.getOffspring()[0].variable().vect()[k] < bound[k].first || pop.getOffspring()[0].variable().vect()[k] > bound[k].second) {
						flag[0] = false;
						break;
					}
				}
				for (size_t k = 0; k < pop.getOffspring()[1].variable().vect().size(); ++k) {
					if (pop.getOffspring()[1].variable().vect()[k] < bound[k].first || pop.getOffspring()[1].variable().vect()[k] > bound[k].second) {
						flag[1] = false;
						break;
					}
				}

			} while (pop.getOffspring()[0].variable().vect() == temp_pop[p[0]].variable().vect() || \
				pop.getOffspring()[0].variable().vect() == temp_pop[p[1]].variable().vect() || \
				pop.getOffspring()[1].variable().vect() == temp_pop[p[0]].variable().vect() || \
				pop.getOffspring()[1].variable().vect() == temp_pop[p[1]].variable().vect()||(!flag[0]&&!flag[1]));

			//ѡ����bound�ڵĽ⣬������bound�ڣ��������淽ʽѡ��
			//ÿ��ȷ��һ���⣬��2���Ӵ���ѡ������Ÿ�����ģ�����������廥��֧�䣬ѡ���븸����С��������
			if (flag[0] && (!flag[1])) {
				all_off.emplace_back(pop.getOffspring()[0].variable().vect());
			}
			else if (flag[1] && (!flag[0])) {
				all_off.emplace_back(pop.getOffspring()[1].variable().vect());
			}
			else {
				bool random_select = true;
				if (random_select) {
					Real rnd_num = rnd->uniform.next();
					if (rnd_num > 0.5) {
						all_off.emplace_back(pop.getOffspring()[1].variable().vect());
					}
					else {
						all_off.emplace_back(pop.getOffspring()[0].variable().vect());
					}
				}
				else {
					Dominance p_ship = objectiveCompare(temp_pop[p[0]].objective(), temp_pop[p[1]].objective(), CAST_CONOP(pro)->optimizeMode());
					auto point1 = pop.getOffspring()[0].variable().vect();
					auto point2 = pop.getOffspring()[1].variable().vect();
					if (p_ship == Dominance::kDominant) {
						//ѡ��p[0]����
						Real s1 = euclideanDistance(point1.begin(), point1.end(), temp_pop[p[0]].variable().vect().begin());
						Real s2 = euclideanDistance(point2.begin(), point2.end(), temp_pop[p[0]].variable().vect().begin());
						if (s1 > s2) {
							all_off.emplace_back(pop.getOffspring()[1].variable().vect());
						}
						else {
							all_off.emplace_back(pop.getOffspring()[0].variable().vect());
						}
					}
					else if (p_ship == Dominance::kDominated) {
						//ѡ��p[1]����
						Real s1 = euclideanDistance(point1.begin(), point1.end(), temp_pop[p[1]].variable().vect().begin());
						Real s2 = euclideanDistance(point2.begin(), point2.end(), temp_pop[p[1]].variable().vect().begin());
						if (s1 > s2) {
							all_off.emplace_back(pop.getOffspring()[1].variable().vect());
						}
						else {
							all_off.emplace_back(pop.getOffspring()[0].variable().vect());
						}
					}
					else if (p_ship == Dominance::kNonDominated) {
						//ѡ�븸������С�����С��
						Real s11 = euclideanDistance(point1.begin(), point1.end(), temp_pop[p[0]].variable().vect().begin());
						Real s12 = euclideanDistance(point1.begin(), point1.end(), temp_pop[p[1]].variable().vect().begin());
						Real s1 = s11 < s12 ? s11 : s12;
						Real s21 = euclideanDistance(point2.begin(), point2.end(), temp_pop[p[0]].variable().vect().begin());
						Real s22 = euclideanDistance(point2.begin(), point2.end(), temp_pop[p[1]].variable().vect().begin());
						Real s2 = s21 < s22 ? s21 : s22;
						if (s1 > s2) {
							all_off.emplace_back(pop.getOffspring()[1].variable().vect());
						}
						else {
							all_off.emplace_back(pop.getOffspring()[0].variable().vect());
						}
					}
				}
			}
		}
		return all_off;
	}

	//DE���ӻ�����Ⱥ����
	std::vector<std::vector<Real>> SPMOEA::sampleByDE(SPMOEA_pop& pop, std::vector<std::pair<Real, Real>>& bound, size_t off_num, size_t k,Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> all_off;
		//auto search_bound = CAST_CONOP(pro)->boundary();
		//����DE��Ⱥ
		PopMODE temp_pop(pop.size(), pro);
		for (size_t i = 0; i < pop.size(); ++i) {
			temp_pop[i].variable().vect() = pop[i].variable().vect();
			temp_pop[i].objective() = pop[i].objective();
		}
		temp_pop.setParamMODE(1./CAST_CONOP(pro)->numberVariables(), 20);//���ñ�����ʺͲ���:0.1,20
		temp_pop.setParameter(1, 0.5);//���ý����ʺ���������:0.6,0.5
		for (size_t j = 0; j < off_num; ++j) {
			std::vector<size_t> select_parent(3);
			//����һ��DE����
			IndDE indDE(pop[0]);
			do {//���ظ�����
				select_parent[0] = temp_pop.tournamentSelection(pro, rnd,k);
				while (1) {
					select_parent[1] = temp_pop.tournamentSelection(pro, rnd,k);
					if (select_parent[1] != select_parent[0])
						break;
				}
				while (1) {
					select_parent[2] = temp_pop.tournamentSelection(pro, rnd,k);
					if (select_parent[2] != select_parent[0] && select_parent[2] != select_parent[1])
						break;
				}
				temp_pop.crossMutate(select_parent, indDE, pro, rnd);
				//Խ�紦��
				repairSol(indDE.variable().vect(), bound, rnd);
			} while (indDE.variable().vect() == temp_pop[select_parent[0]].variable().vect() || \
				indDE.variable().vect() == temp_pop[select_parent[1]].variable().vect() || \
				indDE.variable().vect() == temp_pop[select_parent[2]].variable().vect());
			
			if (ifSame(indDE.variable().vect(), temp_pop[select_parent[0]].variable().vect())) {
				size_t a = 1;
			}
			//pop.getOffspring()[j].setCounter(0);
			Solution<> ind1(temp_pop[select_parent[0]]);
			Solution<> ind2(temp_pop[select_parent[1]]);
			Solution<> ind3(temp_pop[select_parent[2]]);
			indDE.evaluate(pro,alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
			temp_pair.emplace_back(std::make_shared<Solution<>>(indDE));
			m_interactive_sol_pair.emplace_back(temp_pair);

			all_off.emplace_back(indDE.variable().vect());
		}
		return all_off;
	}

	//DE���ӻ��ڸ������
	std::vector<std::vector<Real>> SPMOEA::sampleByDE(SPMOEA_pop& pop, size_t base_inx, std::vector<std::pair<Real, Real>>& bound, size_t off_num, size_t k, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> all_off;
		//auto search_bound = CAST_CONOP(pro)->boundary();
		//����DE��Ⱥ
		PopMODE temp_pop(pop.size(), pro);
		for (size_t i = 0; i < pop.size(); ++i) {
			temp_pop[i].variable().vect() = pop[i].variable().vect();
			temp_pop[i].objective() = pop[i].objective();
		}
		temp_pop.setParamMODE(1./CAST_CONOP(pro)->numberVariables(), 20);//���ñ�����ʺͲ���
		temp_pop.setParameter(1, 0.5);//���ý����ʺ���������
		for (size_t j = 0; j < off_num; ++j) {
			std::vector<size_t> select_parent(3);
			//����һ��DE����
			IndDE indDE(pop[0]);
			select_parent[0] = base_inx;
			do {//���ظ�����
				while (1) {
					select_parent[1] = temp_pop.tournamentSelection(pro, rnd, k);
					if (select_parent[1] != select_parent[0])
						break;
				}
				while (1) {
					select_parent[2] = temp_pop.tournamentSelection(pro, rnd, k);
					if (select_parent[2] != select_parent[0] && select_parent[2] != select_parent[1])
						break;
				}
				temp_pop.crossMutate(select_parent, indDE, pro, rnd);
				//Խ�紦��
				repairSol(indDE.variable().vect(), bound, rnd);
			} while (indDE.variable().vect() == temp_pop[select_parent[0]].variable().vect() || \
				indDE.variable().vect() == temp_pop[select_parent[1]].variable().vect() || \
				indDE.variable().vect() == temp_pop[select_parent[2]].variable().vect());

			if (ifSame(indDE.variable().vect(), temp_pop[select_parent[0]].variable().vect())) {
				size_t a = 1;
			}
			//pop.getOffspring()[j].setCounter(0);
			Solution<> ind1(temp_pop[select_parent[0]]);
			Solution<> ind2(temp_pop[select_parent[1]]);
			Solution<> ind3(temp_pop[select_parent[2]]);
			indDE.evaluate(pro, alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
			temp_pair.emplace_back(std::make_shared<Solution<>>(indDE));
			m_interactive_sol_pair.emplace_back(temp_pair);

			all_off.emplace_back(indDE.variable().vect());
		}
		return all_off;
	}

	//�����ӿռ��ڸ���Ľ���
	std::vector<std::vector<Real>> SPMOEA::sampleBetweenSpaces(std::vector<size_t>& front_link_space1, std::vector<size_t>& front_link_space2, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		std::vector<std::vector<size_t>> pair_spaces;
		pair_spaces.emplace_back(front_link_space1);
		pair_spaces.emplace_back(front_link_space2);
		auto select_pair_spaces = findCloseSpaces(pair_spaces, rnd);
		//�ֱ���ÿһ������ͨ���Ͻ�����С�����ӿռ�Ľ���
		size_t inx1 = std::get<0>(select_pair_spaces[0]);
		size_t inx2 = std::get<1>(select_pair_spaces[0]);
		int method = 1;
		if (method == 1) {//���������ӿռ����������ʽ����
			//�ӿռ�ǰ�ظ��彻��
			auto& front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
			auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_in_subspace;
			size_t s1 = std::floor(front_sols1.size() * rnd->uniform.next());
			size_t s2 = std::floor(front_sols2.size() * rnd->uniform.next());
			auto p1 = front_sols1[s1]->variable().vect();
			auto p2 = front_sols2[s2]->variable().vect();

			for (size_t k = 0; k < sample_num; ++k) {
				auto off = vectorOperator(p1, p2, 0, pro, alg, rnd);
				out_off.emplace_back(off);
			}
		}
		else if (method == 2) {//�������ӿռ����ӷ�ʽ����
			auto& front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
			auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_in_subspace;
			//����Щ�������һ����Ⱥ�������½�
			SPMOEA_pop temp_pop(0, pro);
			for (size_t j = 0; j < front_sols1.size(); ++j) {
				temp_pop.append(*front_sols1[j]);
			}
			for (size_t j = 0; j < front_sols2.size(); ++j) {
				temp_pop.append(*front_sols2[j]);
			}
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring().append(temp_pop[j]);
			}
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				temp_pop.getOffspring().append(temp_pop[j]);
			}
			size_t kk = 2;
			std::vector<std::vector<Real>> off;
			if (temp_pop.size() < 3) {
				off = sampleByGA(temp_pop, bound, 5, pro, rnd);
			}
			else {
				off = sampleByDE(temp_pop, bound, 5, 2, pro,alg, rnd);
			}
			out_off = off;
		}
		return out_off;
	}

	//���������ӿռ��ڲ�������
	std::vector<std::vector<Real>> SPMOEA::sampleInSpace(size_t parent_space, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg,Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		for (size_t i = 0; i < sample_num; ++i) {
			auto& front_sols = getMO_HLC().getSubspaceInfo(parent_space).m_subspace_front_sol;
			auto& his_sols = getMO_HLC().getSubspaceInfo(parent_space).m_history_inds;
			if (his_sols.size() < 2) {
				//�������ֲ��Ŷ�
				auto off = sampleInRange(*front_sols[0], 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {
				//��һ���ӿռ�ǰ�ؽ�Ϊ���㣬��һ��֧���Ϊ��������������½�
				size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				//�����Ͻ��ĸ���
				auto base_sol = his_sols[inx1]->variable().vect();
				Solution<> ind2(*his_sols[inx1]);
				std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
				Real min_span = INT16_MAX;
				auto box = getMO_HLC().subspaceTree().getBox(parent_space);
				for (size_t j = 0; j < box.size(); ++j) {
					if (min_span > (box[j].second - box[j].first)) {
						min_span = box[j].second - box[j].first;
					}
				}
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto p2 = his_sols[j]->variable().vect();
					auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
					if (j ==inx1) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				std::vector<size_t> candidate_inx;
				for (size_t j = 0; j < temp_dist.size(); ++j) {
					if (temp_dist[j] < min_span) {
						candidate_inx.push_back(j);
					}
				}
				//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
				size_t inx2 = 0;
				if (candidate_inx.size() > 0) {
					inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
				}
				else {
					do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					} while (inx1 == inx2);
				}
				//�Ƚ��������֧���ϵ
				std::vector<Real> delta_v;
				std::vector<Real> sol2= his_sols[inx2]->variable().vect();
				Solution<> ind1(*his_sols[inx2]);
				auto dominance_ship = objectiveCompare(his_sols[inx1]->objective(), his_sols[inx2]->objective(), CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(2 * (sol2[j] - base_sol[j]));
					}
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real coff = 1.;
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						coff *= -1;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(coff * (base_sol[j] - sol2[j]));
					}
				}
				else if (dominance_ship == Dominance::kEqual) {
					size_t a = 1;
					//�������ֲ��Ŷ�
					auto sol = base_sol;
					auto box = getMO_HLC().subspaceTree().getBox(parent_space);
					//��ÿһά������ߵ�ֵ
					std::vector<Real> min_v;
					for (size_t j = 0; j < box.size(); ++j) {
						auto v1 = std::fabs(sol[j] - box[j].first);
						auto v2 = std::fabs(sol[j] - box[j].second);
						min_v.push_back(v1 > v2 ? v2 : v1);
					}
					std::vector<Real> new_sol(box.size(), 0.);
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < box.size(); ++j) {
						dim_inx.push_back(j);
					}
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - new_sol[j]);
					}
				}
				std::vector<Real> new_sol(delta_v.size());
				Real ratio =rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}

				Solution<> ind3(ind2);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				out_off.emplace_back(new_sol);
			}
		}
		return out_off;
	}

	//����ĳ�������������Խ��������
	std::vector<std::vector<Real>> SPMOEA::sampleBySolution(Solution<>& sol, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			if (his_sols.size() < 2) {
				auto off = sampleInRange(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {
				auto base_sol = sol.variable().vect();
				Solution<> ind2(sol);
				//�����Ͻ��ĸ���
				std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
				Real min_span = INT16_MAX;
				for (size_t j = 0; j < box.size(); ++j) {
					if (min_span > (box[j].second - box[j].first)) {
						min_span = box[j].second - box[j].first;
					}
				}
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto p2 = his_sols[j]->variable().vect();
					auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
					if (dist <= 0.) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				std::vector<size_t> candidate_inx;
				for (size_t j = 0; j < temp_dist.size(); ++j) {
					if (temp_dist[j] < min_span) {
						candidate_inx.push_back(j);
					}
				}
				//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
				size_t inx2 = 0;
				std::vector<Real> sol2;
				if (candidate_inx.size() > 0) {
					inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
					sol2 = his_sols[inx2]->variable().vect();
				}
				else {
					do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						sol2 = his_sols[inx2]->variable().vect();
					} while (ifSame(base_sol, sol2));
				}
				std::vector<Real> delta_v;
				//�Ƚ��������֧���ϵ
				Solution<> ind1(ind2);
				ind1.variable() = his_sols[inx2]->variable();
				ind1.objective() = his_sols[inx2]->objective();

				auto dominance_ship = objectiveCompare(sol.objective(), his_sols[inx2]->objective(), CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(2 * (sol2[j] - base_sol[j]));
					}
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real coff = 1.;
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						coff *= -1;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(coff * (base_sol[j] - sol2[j]));
					}
				}
				else if (dominance_ship == Dominance::kEqual) {
					size_t a = 1;
					//�������ֲ��Ŷ�
					auto sol = base_sol;
					//��ÿһά������ߵ�ֵ
					std::vector<Real> min_v;
					for (size_t j = 0; j < box.size(); ++j) {
						auto v1 = std::fabs(sol[j] - box[j].first);
						auto v2 = std::fabs(sol[j] - box[j].second);
						min_v.push_back(v1 > v2 ? v2 : v1);
					}
					std::vector<Real> new_sol(box.size(), 0.);
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < box.size(); ++j) {
						dim_inx.push_back(j);
					}
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - new_sol[j]);
					}
				}
				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}
				/*for (size_t j = 0; j < new_sol.size(); ++j) {
					auto temp = new_sol[j];
					new_sol[j] += ((2 * rand - 1) * (raodong / 2) * delta_v[j]);
					if (new_sol[j] < bound[j].first || new_sol[j] > bound[j].second) {
						new_sol[j] = temp;
					}
				}*/
				
				Solution<> ind3(ind2);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				out_off.emplace_back(new_sol);
			}
		}
		return out_off;
	}

	//����ĳ���������������������ӿռ�ǰ�ؽ⿿������չ���ƽ�
	std::vector<std::vector<Real>> SPMOEA::sampleInSolution(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			if (his_sols.size() < 3) {
				auto off=sampleInRange(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//���ӿռ��ǰ�ؽ��ƽ�
				auto base_sol = sol.variable().vect();
				Solution<> ind2(sol);
				auto& front_sols = getMO_HLC().getSubspaceInfo(space).m_subspace_front_sol;
				//�ȿ��ý��ǲ����ӿռ�ǰ�ؽ�
				bool front_flag = false;
				for (size_t j = 0; j < front_sols.size(); ++j) {
					auto temp_sol = front_sols[j]->variable().vect();
					if (ifSame(base_sol,temp_sol)) {
						front_flag = true;
						break;
					}
				}
				Solution<> ind1(ind2);
				std::vector<Real> delta_v;
				if (front_flag) {//����ѡ��������
					//�����Ͻ��ĸ���
					std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
					Real min_span = INT16_MAX;
					for (size_t j = 0; j < box.size(); ++j) {
						if (min_span > (box[j].second - box[j].first)) {
							min_span = box[j].second - box[j].first;
						}
					}
					for (size_t j = 0; j < his_sols.size(); ++j) {
						auto p2 = his_sols[j]->variable().vect();
						auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
						if (dist <= 0.) {
							temp_dist.push_back(INT16_MAX);
						}
						else {
							temp_dist.push_back(dist);
						}
					}
					std::vector<size_t> candidate_inx;
					for (size_t j = 0; j < temp_dist.size(); ++j) {
						if (temp_dist[j] < min_span) {
							candidate_inx.push_back(j);
						}
					}
					//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
					size_t inx2 = 0;
					std::vector<Real> sol2;
					if (candidate_inx.size() > 0) {
						inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
						sol2 = his_sols[inx2]->variable().vect();
					}
					else {
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
							sol2 = his_sols[inx2]->variable().vect();
						} while (ifSame(base_sol, sol2));
					}
					//�Ƚ��������֧���ϵ
					
					ind1.variable() = his_sols[inx2]->variable();
					ind1.objective() = his_sols[inx2]->objective();

					auto dominance_ship = objectiveCompare(sol.objective(), his_sols[inx2]->objective(), CAST_CONOP(pro)->optimizeMode());
					if (dominance_ship == Dominance::kDominant) {
						for (size_t j = 0; j < base_sol.size(); ++j) {
							delta_v.push_back(base_sol[j] - sol2[j]);
						}
					}
					else if (dominance_ship == Dominance::kDominated) {
						for (size_t j = 0; j < base_sol.size(); ++j) {
							delta_v.push_back(2 * (sol2[j] - base_sol[j]));
						}
					}
					else if (dominance_ship == Dominance::kNonDominated) {
						Real coff = 1.;
						Real rand = rnd->uniform.next();
						if (rand > 0.5) {
							coff *= -1;
						}
						for (size_t j = 0; j < base_sol.size(); ++j) {
							delta_v.push_back(coff * (base_sol[j] - sol2[j]));
						}
					}
					else if (dominance_ship == Dominance::kEqual) {
						size_t a = 1;
						//�������ֲ��Ŷ�
						auto sol = base_sol;
						//��ÿһά������ߵ�ֵ
						std::vector<Real> min_v;
						for (size_t j = 0; j < box.size(); ++j) {
							auto v1 = std::fabs(sol[j] - box[j].first);
							auto v2 = std::fabs(sol[j] - box[j].second);
							min_v.push_back(v1 > v2 ? v2 : v1);
						}
						std::vector<Real> new_sol(box.size(), 0.);
						std::vector<size_t> dim_inx;
						for (size_t j = 0; j < box.size(); ++j) {
							dim_inx.push_back(j);
						}
						while (dim_inx.size() > 0) {
							std::vector<size_t> temp;
							for (size_t j = 0; j < dim_inx.size(); ++j) {
								new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
								if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
									temp.push_back(dim_inx[j]);
								}
							}
							dim_inx = temp;
						}
						for (size_t j = 0; j < base_sol.size(); ++j) {
							delta_v.push_back(base_sol[j] - new_sol[j]);
						}
					}
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					size_t inx2=0;
					for (size_t j = 0; j < front_sols.size(); ++j) {
						auto objs = front_sols[j]->objective();
						auto dominance_ship = objectiveCompare(sol.objective(), objs, CAST_CONOP(pro)->optimizeMode());
						if (dominance_ship == Dominance::kDominated) {
							inx2 = j;
							break;
						}
					}
					auto &sol2 = front_sols[inx2]->variable().vect();
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(sol2[j]- base_sol[j]);
					}
					ind1.variable() = front_sols[inx2]->variable();
					ind1.objective() = front_sols[inx2]->objective();
				}
				
				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}

				/*for (size_t j = 0; j < new_sol.size(); ++j) {
					auto temp = new_sol[j];
					new_sol[j] += ((2 * rand - 1) * (raodong / 2) * delta_v[j]);
					if (new_sol[j] < bound[j].first || new_sol[j] > bound[j].second) {
						new_sol[j] = temp;
					}
				}*/

				Solution<> ind3(ind2);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				out_off.emplace_back(new_sol);
			}
		}
		return out_off;
	}

	//�����ӿռ�������������������ӿռ�ǰ�ؽ⿿������չ���ƽ�����������
	std::vector<std::vector<Real>> SPMOEA::sampleByImprove(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//���ӿռ��ǰ�ؽ��ƽ�
				auto base_sol = sol.variable().vect();
				//�ȿ��ý��ǲ����ӿռ�ǰ�ؽ�
				size_t sol_inx = 0;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto temp_sol = his_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						sol_inx = j;
						break;
					}
				}
				//�ӿռ���ʷ������ࣺǰ�غͷ�ǰ��
				std::vector<size_t> subspace_front_inx;
				std::vector<size_t> subspace_behind_inx;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					if (his_sols[j]->getType() == 0) {//type==0, is front sol
						subspace_front_inx.push_back(j);
					}
					else {
						subspace_behind_inx.push_back(j);
					}
				}
				Real rand = rnd->uniform.next();
				PopDE<> temp_pop;
				std::vector<size_t> candidates;
				candidates.push_back(sol_inx);
				IndDE ind1(sol);
				ind1.variable() = his_sols[sol_inx]->variable();
				ind1.objective() = his_sols[sol_inx]->objective();
				temp_pop.append(ind1);
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					//��һ������ǰ�ƻ���չ
					if (rand < 0.7) {//��չ����ǰ�Ž⽻��
						//һ��ĸ���ѡ���ӿռ�����չ��һ��ĸ���ѡ�������rank�ӿռ���չ
						std::vector<size_t> equal_space;
						auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
						for (auto ii : neighs) {
							if (getMO_HLC().getSubspaceInfo(ii).m_best_rank == getMO_HLC().getSubspaceInfo(space).m_best_rank && getMO_HLC().getSubspaceInfo(space).m_best_rank!=INT16_MAX) {
								equal_space.push_back(ii);
							}
						}
						Real rand2 = rnd->uniform.next();
						if (rand2 < 0.5) {
							//������չ
							if (equal_space.empty()) {
								size_t sele_inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								candidates.push_back(sele_inx1);
								IndDE ind2(sol);
								ind2.variable() = his_sols[sele_inx1]->variable();
								ind2.objective() = his_sols[sele_inx1]->objective();
								temp_pop.append(ind2);
								size_t sele_inx2;
								while (1) {
									sele_inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
									if (sele_inx2 != candidates[1]) {
										candidates.push_back(sele_inx2);
										break;
									}
								}
								IndDE ind3(sol);
								ind3.variable() = his_sols[sele_inx2]->variable();
								ind3.objective() = his_sols[sele_inx2]->objective();
								temp_pop.append(ind3);
							}
							else {
								size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
								IndDE ind2(sol);
								ind2.variable() = getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]->variable();
								ind2.objective() = getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]->objective();
								temp_pop.append(ind2);
								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								IndDE ind3(sol);
								ind3.variable() = his_sols[subspace_front_inx[sele_inx2]]->variable();
								ind3.objective() = his_sols[subspace_front_inx[sele_inx2]]->objective();
								temp_pop.append(ind3);
							}
						}
						else {
							//�ӿռ��ڽ���
							if (subspace_front_inx.size() < 3) {
								size_t sele_inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								candidates.push_back(sele_inx1);
								IndDE ind2(sol);
								ind2.variable() = his_sols[sele_inx1]->variable();
								ind2.objective() = his_sols[sele_inx1]->objective();
								temp_pop.append(ind2);
								size_t sele_inx2;
								while (1) {
									sele_inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
									if (sele_inx2 != candidates[1]) {
										candidates.push_back(sele_inx2);
										break;
									}
								}
								IndDE ind3(sol);
								ind3.variable() = his_sols[sele_inx2]->variable();
								ind3.objective() = his_sols[sele_inx2]->objective();
								temp_pop.append(ind3);
							}
							else {
								size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								candidates.push_back(subspace_front_inx[sele_inx1]);
								IndDE ind2(sol);
								ind2.variable() = his_sols[subspace_front_inx[sele_inx1]]->variable();
								ind2.objective() = his_sols[subspace_front_inx[sele_inx1]]->objective();
								temp_pop.append(ind2);
								size_t sele_inx2;
								while (1) {
									sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
									if (subspace_front_inx[sele_inx2] != candidates[1]) {
										candidates.push_back(subspace_front_inx[sele_inx2]);
										break;
									}
								}
								IndDE ind3(sol);
								ind3.variable() = his_sols[subspace_front_inx[sele_inx2]]->variable();
								ind3.objective() = his_sols[subspace_front_inx[sele_inx2]]->objective();
								temp_pop.append(ind3);
							}
							
						}
					}
					else {//ǰ�ƣ�����Ž⽻��������������ǰ��֧����Ľ⽻��
						//�������ӿռ�ǰ�ؽ���
						std::vector<size_t> better_space;
						auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
						for (auto ii : neighs) {
							if (getMO_HLC().getSubspaceInfo(ii).m_best_rank < getMO_HLC().getSubspaceInfo(space).m_best_rank) {
								better_space.push_back(ii);
							}
						}
						if (better_space.empty()) {//�ֲ������ӿռ�
							size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							candidates.push_back(subspace_front_inx[sele_inx1]);
							IndDE ind2(sol);
							ind2.variable() = his_sols[subspace_front_inx[sele_inx1]]->variable();
							ind2.objective() = his_sols[subspace_front_inx[sele_inx1]]->objective();
							temp_pop.append(ind2);
							size_t sele_inx2;
							while (1) {
								sele_inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								if (sele_inx2 != candidates[1]) {
									candidates.push_back(sele_inx2);
									break;
								}
							}
							IndDE ind3(sol);
							ind3.variable() = his_sols[sele_inx2]->variable();
							ind3.objective() = his_sols[sele_inx2]->objective();
							temp_pop.append(ind3);
						}
						else {//ѡ������һ�������ӿռ��ǰ�ؽ�
							size_t sele_space = better_space[(size_t)std::floor(better_space.size()*rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							IndDE ind2(sol);
							ind2.variable() =  getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]->variable();
							ind2.objective() = getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]->objective();
							temp_pop.append(ind2);
							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							candidates.push_back(subspace_front_inx[sele_inx2]);
							IndDE ind3(sol);
							ind3.variable() = his_sols[subspace_front_inx[sele_inx2]]->variable();
							ind3.objective() = his_sols[subspace_front_inx[sele_inx2]]->objective();
							temp_pop.append(ind3);
						}
					}
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					//��һ������ǰ�ƻ���չ
					if (rand < 1.1) {//ǰ�ƣ���ǰ�Ž⽻��
						size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						candidates.push_back(subspace_front_inx[sele_inx1]);
						IndDE ind2(sol);
						ind2.variable() = his_sols[subspace_front_inx[sele_inx1]]->variable();
						ind2.objective() = his_sols[subspace_front_inx[sele_inx1]]->objective();
						temp_pop.append(ind2);
						size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						IndDE ind3(sol);
						ind3.variable() = his_sols[subspace_behind_inx[sele_inx2]]->variable();
						ind3.objective() = his_sols[subspace_behind_inx[sele_inx2]]->objective();
						temp_pop.append(ind3);
						candidates.push_back(subspace_behind_inx[sele_inx2]);
					}
					//else {//��չ������Ž⽻��
					//	if (subspace_behind_inx.size() < 3) {
					//		//�ӿռ��ڽ���
					//		size_t sele_inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					//		candidates.push_back(sele_inx1);
					//		IndDE ind2(sol);
					//		ind2.variable() = his_sols[sele_inx1]->variable();
					//		ind2.objective() = his_sols[sele_inx1]->objective();
					//		temp_pop.append(ind2);
					//		size_t sele_inx2;
					//		while (1) {
					//			sele_inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					//			if (sele_inx2 != candidates[1]) {
					//				candidates.push_back(sele_inx2);
					//				break;
					//			}	
					//		}
					//		IndDE ind3(sol);
					//		ind3.variable() = his_sols[sele_inx2]->variable();
					//		ind3.objective() = his_sols[sele_inx2]->objective();
					//		temp_pop.append(ind3);
					//	}
					//	else {
					//		size_t sele_inx1 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
					//		candidates.push_back(subspace_behind_inx[sele_inx1]);
					//		IndDE ind2(sol);
					//		ind2.variable() = his_sols[subspace_behind_inx[sele_inx1]]->variable();
					//		ind2.objective() = his_sols[subspace_behind_inx[sele_inx1]]->objective();
					//		temp_pop.append(ind2);
					//		size_t sele_inx2;
					//		while (1) {
					//			sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
					//			if (subspace_behind_inx[sele_inx2] != candidates[1]) {
					//				candidates.push_back(subspace_behind_inx[sele_inx2]);
					//				break;
					//			}	
					//		}
					//		IndDE ind3(sol);
					//		ind3.variable() = his_sols[subspace_behind_inx[sele_inx2]]->variable();
					//		ind3.objective() = his_sols[subspace_behind_inx[sele_inx2]]->objective();
					//		temp_pop.append(ind3);
					//	}
					//}
				}
				temp_pop[0].mutate(1.2*rnd->uniform.next(), &temp_pop[0], &temp_pop[1], &temp_pop[2], pro);
				temp_pop.recombine(0, rnd, pro);
				temp_pop[0].trial().evaluate(pro, this);

				Solution<> ind4(sol);
				ind4.variable() = temp_pop[0].trial().variable();
				ind4.objective() = temp_pop[0].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[0]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[1]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[2]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				m_interactive_sol_pair.emplace_back(temp_pair);
				if (ifSame(temp_pop[0].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[0].trial().variable().vect());
			}
		}
		return out_off;
	}

	//�����ӿռ�������������������ӿռ�ǰ�ؽ⿿������չ���ƽ�����������
	std::vector<std::vector<Real>> SPMOEA::samplePushorExt(Solution<>& sol, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		std::vector<size_t> equal_space;
		auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank == getMO_HLC().getSubspaceInfo(space).m_best_rank && getMO_HLC().getSubspaceInfo(space).m_best_rank != INT16_MAX) {
				equal_space.push_back(ii);
			}
		}
		std::vector<size_t> better_space;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank < getMO_HLC().getSubspaceInfo(space).m_best_rank) {
				better_space.push_back(ii);
			}
		}
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//���������ӿռ�λ�ø���ѡ���ƽ�����չ
				auto base_sol = sol.variable().vect();
				//�ȿ��ý��ǲ����ӿռ�ǰ�ؽ�
				size_t sol_inx = 0;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto temp_sol = his_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						sol_inx = j;
						break;
					}
				}
				//�ӿռ���ʷ������ࣺǰ�غͷ�ǰ��
				std::vector<size_t> subspace_front_inx;
				std::vector<size_t> subspace_behind_inx;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					if (his_sols[j]->getType() == 0) {//type==0, is front sol
						subspace_front_inx.push_back(j);
					}
					else {
						subspace_behind_inx.push_back(j);
					}
				}
				Real rand = rnd->uniform.next();
				PopDE<> temp_pop;
				std::vector<size_t> candidates;
				candidates.push_back(sol_inx);
				Real extend_prob = 0.3;
				std::vector<size_t> candidate_inx;
				temp_pop.append(IndDE(*his_sols[sol_inx]));
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					//��һ������ǰ�ƻ���չ
					if (rand < extend_prob) {//��չ����ǰ�Ž⽻��
						//���ӿռ�ǰ��������ͬ��rank�ӿռ�ǰ�غϳ���Ⱥ
						for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
							if (sol_inx != subspace_front_inx[j]) {
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[j]]));
							}
						}
						size_t nei_space_inx;
						if (equal_space.size()>0) {
							nei_space_inx = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
							auto& front_sols2 = getMO_HLC().getSubspaceInfo(nei_space_inx).m_subspace_front_sol;
							for (size_t j = 0; j < front_sols2.size(); ++j) {
								temp_pop.append(IndDE(*front_sols2[j]));
							}
						}
						if (temp_pop.size() < 3) {//�ڲ�����
							for (size_t j = 0; j < subspace_behind_inx.size(); ++j) {
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[j]]));
							}
						}
						candidate_inx.push_back(0);
						candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
						size_t inx = 0;
						do {
							inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
						} while (inx==candidate_inx[1]);
						candidate_inx.push_back(inx);
					}
					else {//ǰ�ƣ�����Ž⽻��������������ǰ��֧����Ľ⽻��
						//�������ӿռ�ǰ�ؽ���
						if (better_space.empty()) {//�ֲ������ӿռ䣬�ڲ�����
							if (subspace_behind_inx.empty()) {
								for (size_t j = 0; j < his_sols.size(); ++j) {
									if (sol_inx != j) {
										temp_pop.append(IndDE(*his_sols[j]));
									}
								}
								candidate_inx.push_back(0);
								candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
								size_t inx = 0;
								do {
									inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
								} while (inx == candidate_inx[1]);
								candidate_inx.push_back(inx);
							}
							else {
								candidate_inx = {0,1,2};
								size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));
								
								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
							}
						}
						else {//ѡ������һ�������ӿռ��ǰ�ؽ�
							candidate_inx = {0,1,2};
							size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));
						
							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
						}
					}
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					//��һ������ǰ�ƻ���չ
					if (rand < 1-extend_prob) {//ǰ�ƣ���ǰ�Ž⽻��
						candidate_inx = {0,1,2};
						size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));
						size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
					}
					else {
						//�ӿռ��ڽ���
						for (size_t j = 0; j < his_sols.size(); ++j) {
							if (sol_inx != j) {
								temp_pop.append(IndDE(*his_sols[j]));
							}
						}
						candidate_inx.push_back(0);
						candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
						size_t inx = 0;
						do {
							inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
						} while (inx == candidate_inx[1]);
						candidate_inx.push_back(inx);
					}
				}
				temp_pop[candidate_inx[0]].mutate(temp_pop.scalingFactor(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);
				temp_pop[candidate_inx[0]].trial().evaluate(pro, this);

				Solution<> ind4(sol);
				ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				m_interactive_sol_pair.emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
		}
		return out_off;
	}

	//�����ӿռ�������������������ӿռ�ǰ�ؽ��ƽ�����������
	std::vector<std::vector<Real>> SPMOEA::samplePush(Solution<>& sol, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		std::vector<size_t> equal_space;
		auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank == getMO_HLC().getSubspaceInfo(space).m_best_rank && getMO_HLC().getSubspaceInfo(space).m_best_rank != INT16_MAX) {
				equal_space.push_back(ii);
			}
		}
		std::vector<size_t> better_space;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank < getMO_HLC().getSubspaceInfo(space).m_best_rank) {
				better_space.push_back(ii);
			}
		}
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//���������ӿռ�λ�ø���ѡ���ƽ�����
				auto base_sol = sol.variable().vect();
				//�ȿ��ý��ǲ����ӿռ�ǰ�ؽ�
				size_t sol_inx = 0;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto temp_sol = his_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						sol_inx = j;
						break;
					}
				}
				//�ӿռ���ʷ������ࣺǰ�غͷ�ǰ��
				std::vector<size_t> subspace_front_inx;
				std::vector<size_t> subspace_behind_inx;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					if (his_sols[j]->getType() == 0) {//type==0, is front sol
						subspace_front_inx.push_back(j);
					}
					else {
						subspace_behind_inx.push_back(j);
					}
				}
				Real rand = rnd->uniform.next();
				PopDE<> temp_pop;
				Real external_prob = 0.5;
				std::vector<size_t> candidate_inx;
				temp_pop.append(IndDE(*his_sols[sol_inx]));
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					if (rand < external_prob) {//��һ���ĸ�������������ӿռ佻��
						if (better_space.empty()) {//�ֲ������ӿռ䣬�ڲ�����
							if (subspace_behind_inx.empty()) {
								//�ӿռ��ڽ���
								candidate_inx = { 0,1,2 };
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {
								candidate_inx = { 0,1,2 };
								size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
							}
						}
						else {//ѡ������һ�������ӿռ��ǰ�ؽ�
							candidate_inx = { 0,1,2 };
							size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
						}
					}
					else {//ǰ�ƣ�����Ž⽻��
						if (subspace_behind_inx.empty()) {
							if (better_space.empty()) {//�ֲ������ӿռ䣬�ڲ�����
								//�ӿռ��ڽ���
								candidate_inx = {0,1,2};
								size_t inx1=(size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {//ѡ������һ�������ӿռ��ǰ�ؽ�
								candidate_inx = { 0,1,2 };
								size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
							}
						}
						else {
							candidate_inx = { 0,1,2 };
							size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));

							size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
						}
					}
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					//��һ������ǰ�ƻ���չ
					if (rand < external_prob) {//ǰ�ƣ���ǰ�Ž⽻��
						candidate_inx = { 0,1,2 };
						size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));
						size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
					}
					else {//��������õ��ӿռ佻��
						if (better_space.empty()) {//�ֲ������ӿռ䣬�ڲ�����
							candidate_inx = { 0,1,2 };
							size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));

							size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
						}
						else {//ѡ������һ�������ӿռ��ǰ�ؽ�
							candidate_inx = { 0,1,2 };
							size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
						}
					}
				}
				temp_pop[candidate_inx[0]].mutate(rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);
				temp_pop[candidate_inx[0]].trial().evaluate(pro, this);

				Solution<> ind4(sol);
				ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				m_interactive_sol_pair.emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
		}
		return out_off;
	}

	//���ӿռ��ڲ���������ʷ�����ӿռ�ǰ�ؽ�Ľ�����չ���ƽ�����������
	std::vector<std::vector<Real>> SPMOEA::sampleInSubspace(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space_inx);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space_inx).m_history_inds;
			if (his_sols.size() < 5) {
				auto off = sampleInRange(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
					if (ifSame(sol.variable().vect(), ii)) {
						size_t a = 1;
					}
				}
			}
			else {//���ӿռ��ǰ�ؽ��ƽ�
				auto base_sol = sol.variable().vect();
				Solution<> ind1(sol);
				auto& front_sols = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
				//�ȿ��ý��ǲ����ӿռ�ǰ�ؽ�
				bool front_flag = false;
				for (size_t j = 0; j < front_sols.size(); ++j) {
					auto temp_sol = front_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						front_flag = true;
						break;
					}
				}
				Solution<> ind2(ind1);
				Solution<> ind3(ind1);
				std::vector<Real> delta_v;
				if (front_flag) {
					//Ѱ��֧�����������Ҳ��������ҷ�֧������
					size_t inx2 = 0;
					std::vector<Real> sol2;
					//std::vector<size_t> dominant_inx;
					//std::vector<size_t> non_dominant_inx;
					//for (size_t j = 0; j < his_sols.size(); ++j) {
					//	auto dominance_ship = objectiveCompare(sol.objective(), his_sols[j]->objective(), CAST_CONOP(pro)->optimizeMode());
					//	if (dominance_ship == Dominance::kDominant) {
					//		dominant_inx.push_back(j);
					//	}
					//}
					//if (dominant_inx.size() > 0) {
					//	inx2 = dominant_inx[(size_t)std::floor(dominant_inx.size() * rnd->uniform.next())];
					//}
					//else {
					//	//��ʷ�������ѡ
					//	do {
					//		inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					//	} while (ifSame(base_sol,his_sols[inx2]->variable().vect()));
					//}
					//��ʷ�������ѡ
					//�����Ͻ��ĸ���
					std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
					Real min_span = INT16_MAX;
					for (size_t j = 0; j < box.size(); ++j) {
						if (min_span > (box[j].second - box[j].first)) {
							min_span = box[j].second - box[j].first;
						}
					}
					for (size_t j = 0; j < his_sols.size(); ++j) {
						auto p2 = his_sols[j]->variable().vect();
						auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
						if (dist <= 0.) {
							temp_dist.push_back(INT16_MAX);
						}
						else {
							temp_dist.push_back(dist);
						}
					}
					std::vector<size_t> candidate_inx;
					for (size_t j = 0; j < temp_dist.size(); ++j) {
						if (temp_dist[j] < min_span) {
							candidate_inx.push_back(j);
						}
					}
					//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
					if (candidate_inx.size() > 0) {
						inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
						sol2 = his_sols[inx2]->variable().vect();
					}
					else {
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
							sol2 = his_sols[inx2]->variable().vect();
						} while (ifSame(base_sol, sol2));
					}
					/*do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					} while (ifSame(base_sol, his_sols[inx2]->variable().vect()));
					sol2 = his_sols[inx2]->variable().vect();*/
					for (size_t j = 0; j < sol2.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
					ind2.variable() = his_sols[inx2]->variable();
					ind2.objective() = his_sols[inx2]->objective();
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					/*std::vector<size_t> dominant_inx;
					for (size_t j = 0; j < front_sols.size(); ++j) {
						auto objs = front_sols[j]->objective();
						auto dominance_ship = objectiveCompare(sol.objective(), objs, CAST_CONOP(pro)->optimizeMode());
						if (dominance_ship == Dominance::kDominated) {
							dominant_inx.push_back(j);
						}
					}
					size_t inx2 = dominant_inx[(size_t)std::floor(dominant_inx.size() * rnd->uniform.next())];*/
					//ѡ��һ����ʷ��ǰ�ؽ⣬
					std::vector<size_t> his_inx;
					for (size_t j = 0; j < his_sols.size(); ++j) {
						his_inx.push_back(j);
					}
					//����˳��
					rnd->uniform.shuffle(his_inx.begin(), his_inx.end());
					//ѡ��һ����ǰ�ؽ�
					size_t select_his_inx=0;
					for (size_t j = 0; j < his_inx.size(); ++j) {
						auto& temp_sol1 = his_sols[his_inx[j]]->variable().vect();
						bool flag = true;
						for (size_t k = 0; k < front_sols.size(); ++k) {
							auto &temp_sol2 = front_sols[k]->variable().vect();
							if (ifSame(temp_sol1, temp_sol2)) {
								flag = false;
								break;
							}
						}
						if (flag) {
							select_his_inx = j;
							break;
						}
					}
					auto &sol1= his_sols[his_inx[select_his_inx]]->variable().vect();
					base_sol = sol1;
					size_t inx2 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
					auto& sol2 = front_sols[inx2]->variable().vect();
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(sol2[j] - base_sol[j]);
					}
					ind2.variable() = his_sols[his_inx[select_his_inx]]->variable();
					ind2.objective() = his_sols[his_inx[select_his_inx]]->objective();
					ind3.variable() = front_sols[inx2]->variable();
					ind3.objective() = front_sols[inx2]->objective();
				}

				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				if (!front_flag) {
					for (size_t j = 0; j < delta_v.size(); ++j) {
						delta_v[j] *= 1.5;
					}
				}
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				auto mutate_sol = new_sol;
				for (size_t j = 0; j < mutate_sol.size(); ++j) {
					//ÿһά����ʼ��ľ�����Сֵ�����Ŷ���Χ
					Real dist1 = std::fabs(new_sol[j] - base_sol[j]);
					Real dist2 = std::fabs(new_sol[j] - base_sol[j] - delta_v[j]);
					auto amplitude = dist1 < dist2 ? dist1 : dist2;
					auto rand = rnd->uniform.next();
					mutate_sol[j] += ((2 * rand - 1) * amplitude);
				}
				//��Խ�紦��
				for (size_t j = 0; j < mutate_sol.size(); ++j) {
					if (mutate_sol[j] < bound[j].first) {
						mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
					}
					else if (mutate_sol[j] > bound[j].second) {
						mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
					}
				}
				////��new_solΪԲ�ģ���base_sol�ľ���Ϊ�뾶���Ʊ��췶Χ
				//Real radius1 = 0.;
				//auto& end_sol = ind3.variable().vect();
				//for (size_t j = 0; j < new_sol.size(); ++j) {
				//	radius1 += std::pow((new_sol[j] - end_sol[j]), 2);
				//}
				//radius1 = std::sqrt(radius1);
				//auto mutate_sol = new_sol;
				//Real temp_dist = 0.;
				//Real vectorcos = 0.;
				//do {
				//	for (size_t j = 0; j < new_sol.size(); ++j) {
				//		Real randn = rnd->uniform.next();
				//		mutate_sol[j] = new_sol[j] + (2 * randn - 1) * radius1;
				//	}
				//	vectorcos = vectorAngle(new_sol,end_sol,mutate_sol);
				//	temp_dist = euclideanDistance(new_sol.begin(), new_sol.end(), mutate_sol.begin());
				//} while (temp_dist > radius1||vectorcos<0);
				////��Խ�紦��
				//for (size_t j = 0; j < mutate_sol.size(); ++j) {
				//	if (mutate_sol[j] < bound[j].first) {
				//		mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
				//	}
				//	else if (mutate_sol[j] > bound[j].second) {
				//		mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
				//	}
				//}
				Solution<> ind4(ind2);
				ind4.variable().vect() = mutate_sol;
				ind4.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				m_interactive_sol_pair.emplace_back(temp_pair);

				out_off.emplace_back(mutate_sol);
				if (ifSame(sol.variable().vect(), mutate_sol)) {
					size_t a = 1;
				}
			}
		}
		return out_off;
	}

	//�����ӿռ�������Ը���ǰ�ƻ���չ
	std::vector<std::vector<Real>> SPMOEA::sampleInSubspace(size_t space1, size_t ind_inx, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto bound = CAST_CONOP(pro)->boundary();
		//�ҳ���ǰ�������ӿռ�
		std::vector<size_t> nei1;
		nei1.push_back(space1);
		auto neigh = getMO_HLC().getSubspaceInfo(space1).m_sub_neighbors;
		for (auto jj : neigh) {
			if (std::find(nei1.begin(), nei1.end(), space1) != nei1.end()) {
				nei1.push_back(jj);
			}
		}
		auto& his_sols = getMO_HLC().getSubspaceInfo(space1).m_history_inds;
		auto& front_sols = getMO_HLC().getSubspaceInfo(space1).m_subspace_front_sol;
		auto& sol = getPop()[0][ind_inx].variable().vect();
		auto& box = getMO_HLC().subspaceTree().getBox(space1);
		//�Ը���ǰ�ƻ���չ
		Real rand = rnd->uniform.next();
		if (rand > 0.5) {//ǰ��
			//Ѱ��֧�䷽��
			if (his_sols.size() < 5) {
				auto off1 = sampleInRange(getPop()[0][ind_inx], 1,box,pro, alg, rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//Ѱ��֧�����������Ҳ��������ҷ�֧������
				std::map<size_t,std::vector<size_t>> dominant_inx;
				std::map<size_t, std::vector<std::vector<Real>>> dominant_directions;
				std::vector<size_t> non_dominant_inx;
				std::vector<std::vector<Real>> dominance_vectors;
				for (size_t p = 0; p < front_sols.size(); ++p) {
					dominance_vectors.clear();
					for (size_t j = 0; j < his_sols.size(); ++j) {
						auto dominance_ship = objectiveCompare(front_sols[p]->objective(), his_sols[j]->objective(), CAST_CONOP(pro)->optimizeMode());
						if (dominance_ship == Dominance::kDominant) {
							dominant_inx[p].push_back(j);
							std::vector<Real> temp;
							for (size_t k = 0; k < sol.size(); ++k) {
								temp.push_back(front_sols[p]->variable()[k] - his_sols[j]->variable()[k]);
							}
							dominance_vectors.emplace_back(temp);
						}
					}
					dominant_directions.insert(std::make_pair<>(p,dominance_vectors));
				}
				//ѡһ��ǰ�ؽ�
				size_t inx1 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
				if (dominant_directions[inx1].empty()) {
					SPMOEA_pop temp_pop(0, pro);
					for (size_t k = 0; k < front_sols.size(); ++k) {
						temp_pop.append(*front_sols[k]);
					}
					for (size_t k = 0; k < 2; ++k) {
						temp_pop.append(*his_sols[k]);
					}
					size_t kk = 2;
					auto off = sampleByDE(temp_pop, bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						out_off.emplace_back(off[k]);
						if (ifSame(sol, off[k])) {
							size_t a = 1;
						}
					}
				}
				else {
					//����ǰ������
					Solution<> ind1(getPop()[0][ind_inx]);
					auto evolve_direction = dominant_directions[inx1];
					Solution<> ind2(*his_sols[dominant_inx[inx1][0]]);
					Solution<> ind3(*front_sols[inx1]);
					std::vector<Real> sol1 = ind2.variable().vect();
					std::vector<Real> sol2 = ind3.variable().vect();
					std::vector<Real> new_sol(sol1.size());
					std::vector<Real> delta_v(sol1.size(), 0.);
					for (size_t j = 0; j < delta_v.size(); ++j) {
						for (size_t k = 0; k < evolve_direction.size(); ++k) {
							delta_v[j] += evolve_direction[k][j];
						}
						delta_v[j] = delta_v[j]/ evolve_direction.size();
					}
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < delta_v.size(); ++j) {
						dim_inx.push_back(j);
					}
					Real ratio = rnd->uniform.next();
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = sol2[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
						if (dim_inx.size() > 0) {
							ratio /= 1.2;
						}
					}
					auto mutate_sol = new_sol;
					for (size_t j = 0; j < mutate_sol.size(); ++j) {
						//ÿһά����ʼ��ľ�����Сֵ�����Ŷ���Χ
						Real dist1 = std::fabs(new_sol[j] - sol2[j]);
						Real dist2 = std::fabs(new_sol[j] - sol2[j]-delta_v[j]);
						auto amplitude = dist1 < dist2 ? dist1 : dist2;
						auto rand = rnd->uniform.next();
						mutate_sol[j] += ((2 * rand - 1) * amplitude);
					}
					//��Խ�紦��
					for (size_t j = 0; j < mutate_sol.size(); ++j) {
						if (mutate_sol[j] < bound[j].first) {
							mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
						}
						else if (mutate_sol[j] > bound[j].second) {
							mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
						}
					}
					if (ifSame(sol, mutate_sol)) {
						size_t a = 1;
					}

					Solution<> ind4(ind2);
					ind4.variable().vect() = mutate_sol;
					ind4.evaluate(pro, alg);
					std::vector<std::shared_ptr<Solution<>>> temp_pair;
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
					m_interactive_sol_pair.emplace_back(temp_pair);

					out_off.emplace_back(mutate_sol);
				}
			}
		}
		else {//��չ
			//������ǰ�ظ��幹������Ⱥ
			SPMOEA_pop temp_pop(0, pro);
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_subspace_front_sol;
				for (size_t k = 0; k < front_sol.size(); ++k) {
					if (!ifSame(sol, front_sol[k]->variable().vect())) {
						temp_pop.append(*front_sol[k]);
					}
				}
			}
			temp_pop.append(getPop()[0][ind_inx]);
			for (size_t i = 0; i < sample_num; ++i) {
				if (temp_pop.size() < 5) {
					//���岻���������ռ����
					auto off1 = sampleInRange(temp_pop[temp_pop.size() - 1], 1,bound, pro, alg, rnd);
					for (auto jj : off1) {
						out_off.emplace_back(jj);
						if (ifSame(sol, jj)) {
							size_t a = 1;
						}
					}
				}
				else {
					size_t kk = 2;
					auto off = sampleByDE(temp_pop, bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						out_off.emplace_back(off[k]);
						if (ifSame(sol, off[k])) {
							size_t a = 1;
						}
					}
				}
			}
		}
		return out_off;
	}

	//����ǰ���ӿռ�������Ը���ǰ�ƻ���չ
	std::vector<std::vector<Real>> SPMOEA::sampleInFrontSubspace(size_t space1, size_t ind_inx, std::vector<size_t>& front_link_spaces, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto bound = CAST_CONOP(pro)->boundary();
		//�ҳ���ǰ�������ӿռ�
		std::vector<size_t> nei1;
		nei1.push_back(space1);
		auto neigh = getMO_HLC().getSubspaceInfo(space1).m_sub_neighbors;
		for (auto jj : neigh) {
			if (std::find(front_link_spaces.begin(), front_link_spaces.end(), space1) != front_link_spaces.end()) {
				nei1.push_back(jj);
			}
		}
		auto& his_sols = getMO_HLC().getSubspaceInfo(space1).m_history_inds;
		auto& front_sols = getMO_HLC().getSubspaceInfo(space1).m_front_sol_in_subspace;
		auto& sol = getPop()[0][ind_inx].variable().vect();
		auto& box = getMO_HLC().subspaceTree().getBox(space1);
		//�Ը���ǰ�ƻ���չ
		Real rand = rnd->uniform.next();
		if (rand > 0.5) {//ǰ��
			//Ѱ��֧�䷽��
			if (his_sols.size() < 5) {
				auto off1 = sampleInRange(getPop()[0][ind_inx], 1, box, pro, alg, rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//Ѱ��֧�����������Ҳ��������ҷ�֧������
				std::map<size_t, std::vector<size_t>> dominant_inx;
				std::map<size_t, std::vector<std::vector<Real>>> dominant_directions;
				std::vector<size_t> non_dominant_inx;
				std::vector<std::vector<Real>> dominance_vectors;
				for (size_t p = 0; p < front_sols.size(); ++p) {
					dominance_vectors.clear();
					for (size_t j = 0; j < his_sols.size(); ++j) {
						auto dominance_ship = objectiveCompare(front_sols[p]->objective(), his_sols[j]->objective(), CAST_CONOP(pro)->optimizeMode());
						if (dominance_ship == Dominance::kDominant) {
							dominant_inx[p].push_back(j);
							std::vector<Real> temp;
							for (size_t k = 0; k < sol.size(); ++k) {
								temp.push_back(front_sols[p]->variable()[k] - his_sols[j]->variable()[k]);
							}
							dominance_vectors.emplace_back(temp);
						}
					}
					dominant_directions.insert(std::make_pair<>(p, dominance_vectors));
				}
				//ѡһ��ǰ�ؽ�
				size_t inx1 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
				if (dominant_directions[inx1].empty()) {
					SPMOEA_pop temp_pop(0, pro);
					for (size_t k = 0; k < front_sols.size(); ++k) {
						temp_pop.append(*front_sols[k]);
					}
					for (size_t k = 0; k < 2; ++k) {
						temp_pop.append(*his_sols[k]);
					}
					size_t kk = 2;
					auto off = sampleByDE(temp_pop, bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						out_off.emplace_back(off[k]);
						if (ifSame(sol, off[k])) {
							size_t a = 1;
						}
					}
				}
				else {
					//����ǰ������
					Solution<> ind1(getPop()[0][ind_inx]);
					auto evolve_direction = dominant_directions[inx1];
					Solution<> ind2(*his_sols[dominant_inx[inx1][0]]);
					Solution<> ind3(*front_sols[inx1]);
					std::vector<Real> sol1 = ind2.variable().vect();
					std::vector<Real> sol2 = ind3.variable().vect();
					std::vector<Real> new_sol(sol1.size());
					std::vector<Real> delta_v(sol1.size(), 0.);
					for (size_t j = 0; j < delta_v.size(); ++j) {
						for (size_t k = 0; k < evolve_direction.size(); ++k) {
							delta_v[j] += evolve_direction[k][j];
						}
						delta_v[j] = delta_v[j] / evolve_direction.size();
					}
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < delta_v.size(); ++j) {
						dim_inx.push_back(j);
					}
					Real ratio = rnd->uniform.next();
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = sol2[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
						if (dim_inx.size() > 0) {
							ratio /= 1.2;
						}
					}
					auto mutate_sol = new_sol;
					for (size_t j = 0; j < mutate_sol.size(); ++j) {
						//ÿһά����ʼ��ľ�����Сֵ�����Ŷ���Χ
						Real dist1 = std::fabs(new_sol[j] - sol2[j]);
						Real dist2 = std::fabs(new_sol[j] - sol2[j] - delta_v[j]);
						auto amplitude = dist1 < dist2 ? dist1 : dist2;
						auto rand = rnd->uniform.next();
						mutate_sol[j] += ((2 * rand - 1) * amplitude);
					}
					//��Խ�紦��
					for (size_t j = 0; j < mutate_sol.size(); ++j) {
						if (mutate_sol[j] < bound[j].first) {
							mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
						}
						else if (mutate_sol[j] > bound[j].second) {
							mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
						}
					}
					if (ifSame(sol, mutate_sol)) {
						size_t a = 1;
					}

					Solution<> ind4(ind2);
					ind4.variable().vect() = mutate_sol;
					ind4.evaluate(pro, alg);
					std::vector<std::shared_ptr<Solution<>>> temp_pair;
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
					m_interactive_sol_pair.emplace_back(temp_pair);

					out_off.emplace_back(mutate_sol);
				}
			}
		}
		else {//��չ
			SPMOEA_pop temp_pop(0, pro);
			Real randn = rnd->uniform.next();
			if (randn > 0.5) {
				temp_pop.append(getPop()[0][ind_inx]);
				for (size_t j = 0; j < nei1.size(); ++j) {
					auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_subspace_front_sol;
					for (size_t k = 0; k < front_sol.size(); ++k) {
						if (!ifSame(sol, front_sol[k]->variable().vect())) {
							temp_pop.append(*front_sol[k]);
						}
					}
				}
			}
			else {
				//������ͨ�ӿռ���ԵĲ�ؾ������
				auto manifold_dist = calSpaceManifoldDist(front_link_spaces);//�Ż�����
				//�ҳ���Ե�����Զ���ӿռ䣨��չ��
				size_t inx1, inx2;
				std::vector<size_t> temp_max;
				for (size_t j = 0; j < manifold_dist.size(); ++j) {
					temp_max.push_back(*std::max_element(manifold_dist[j].begin(), manifold_dist[j].end()));
				}
				auto index1 = std::distance(temp_max.begin(), std::max_element(temp_max.begin(), temp_max.end()));
				auto index2 = std::distance(manifold_dist[index1].begin(), std::max_element(manifold_dist[index1].begin(), manifold_dist[index1].end()));
				inx1 = front_link_spaces[index1];//��ʵ�ӿռ�����
				inx2 = front_link_spaces[index2];
				Real rand = rnd->uniform.next();
				size_t select_space = rand > 0.5 ? inx1 : inx2;
				//�ҳ���ǰ�������ӿռ�
				std::vector<size_t> nei;
				nei.push_back(select_space);
				auto neigh = getMO_HLC().getSubspaceInfo(select_space).m_sub_neighbors;
				for (auto jj : neigh) {
					nei.push_back(jj);
				}
				//������ǰ�ظ��幹������Ⱥ
				for (size_t j = 0; j < nei.size(); ++j) {
					auto& front_sol = getMO_HLC().getSubspaceInfo(nei[j]).m_subspace_front_sol;
					for (size_t k = 0; k < front_sol.size(); ++k) {
						if (!ifSame(sol, front_sol[k]->variable().vect())) {
							temp_pop.append(*front_sol[k]);
						}
					}
				}
			}
			for (size_t i = 0; i < sample_num; ++i) {
				if (temp_pop.size() < 5) {
					//���岻���������ռ����
					auto off1 = sampleInRange(temp_pop[0], 1, pro, alg, rnd);
					for (auto jj : off1) {
						out_off.emplace_back(jj);
						if (ifSame(sol, jj)) {
							size_t a = 1;
						}
					}
				}
				else {
					size_t kk = 2;
					auto off = sampleByDE(temp_pop,0, bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						out_off.emplace_back(off[k]);
						if (ifSame(sol, off[k])) {
							size_t a = 1;
						}
					}
				}
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA::sampleInSubspace(size_t space_inx, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto& bound = getMO_HLC().subspaceTree().getBox(space_inx);
		size_t M = CAST_CONOP(pro)->numberObjectives();
		size_t c = CAST_CONOP(pro)->numberConstraints();
		for (size_t i = 0; i < sample_num; ++i) {
			//�ӿռ����������
			std::vector<Real> new_sol(bound.size(), 0.);
			for (size_t j = 0; j < bound.size(); ++j) {
				new_sol[j] = bound[j].first + rnd->uniform.next() * (bound[j].second-bound[j].first);
			}
			Solution<> sol(M,c,bound.size());
			Solution<> ind2(sol);
			ind2.variable().vect() = new_sol;
			ind2.evaluate(pro, alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			//temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			m_interactive_sol_pair.emplace_back(temp_pair);
			out_off.emplace_back(new_sol);
		}
		return out_off;

	}

	std::vector<std::vector<Real>> SPMOEA::sampleBySubspace(Solution<>& sol, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& bound = getMO_HLC().subspaceTree().getBox(space_inx);
		for (size_t i = 0; i < sample_num; ++i) {
			//�ӿռ����������
			Solution<> ind1(sol);
			ind1.variable() = sol.variable();
			ind1.objective() = sol.objective();
			std::vector<Real> new_sol(bound.size(), 0.);
			for (size_t j = 0; j < bound.size(); ++j) {
				new_sol[j] = bound[j].first + rnd->uniform.next() * (bound[j].second - bound[j].first);
			}
			Solution<> ind2(sol);
			ind2.variable().vect() = new_sol;
			ind2.evaluate(pro, alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			m_interactive_sol_pair.emplace_back(temp_pair);
			out_off.emplace_back(new_sol);
		}
		return out_off;

	}

	//����ĳ���������ӿռ������Ŷ����Ҳ����������߽�
	std::vector<std::vector<Real>> SPMOEA::sampleInRange(Solution<>& sol, size_t sample_num,const std::vector<std::pair<Real,Real>>&bound, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		//auto& box = getMO_HLC().subspaceTree().getBox(space);
		//auto& bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			Solution<> ind1(sol);
			ind1.variable() = sol.variable();
			ind1.objective() = sol.objective();
			//�������ֲ��Ŷ�
			//��ÿһά������ߵ�ֵ
			std::vector<Real> min_v;
			for (size_t j = 0; j < bound.size(); ++j) {
				auto v1 = std::fabs(sol.variable()[j] - bound[j].first);
				auto v2 = std::fabs(sol.variable()[j] - bound[j].second);
				min_v.push_back(v1 > v2 ? v2 : v1);
			}
			std::vector<Real> new_sol(bound.size(), 0.);
			std::vector<size_t> dim_inx;
			for (size_t j = 0; j < bound.size(); ++j) {
				dim_inx.push_back(j);
			}
			while (dim_inx.size() > 0) {
				std::vector<size_t> temp;
				for (size_t j = 0; j < dim_inx.size(); ++j) {
					new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol.variable()[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
					if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
						temp.push_back(dim_inx[j]);
					}
				}
				dim_inx = temp;
			}
			Solution<> ind2(ind1);
			ind2.variable().vect() = new_sol;
			ind2.evaluate(pro, alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			m_interactive_sol_pair.emplace_back(temp_pair);
			out_off.emplace_back(new_sol);
		}
		return out_off;
	}

	//����ĳ���������ӿռ������Ŷ����Ҳ����������߽�
	std::vector<std::vector<Real>> SPMOEA::sampleInRange(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			Solution<> ind1(sol);
			ind1.variable() = sol.variable();
			ind1.objective() = sol.objective();
			//�������ֲ��Ŷ�
			//��ÿһά������ߵ�ֵ
			std::vector<Real> min_v;
			for (size_t j = 0; j < box.size(); ++j) {
				auto v1 = std::fabs(sol.variable()[j] - box[j].first);
				auto v2 = std::fabs(sol.variable()[j] - box[j].second);
				min_v.push_back(v1 > v2 ? v2 : v1);
			}
			std::vector<Real> new_sol(box.size(), 0.);
			
			do {
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < box.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol.variable()[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
				}
			} while (ifSame(new_sol, sol.variable().vect()));
			if (ifSame(new_sol, sol.variable().vect())) {
				size_t a = 1;
			}
			Solution<> ind2(ind1);
			ind2.variable().vect() = new_sol;
			ind2.evaluate(pro, alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			m_interactive_sol_pair.emplace_back(temp_pair);
			out_off.emplace_back(new_sol);
		}
		return out_off;
	}

	//��ͨǰ���ӿռ����
	std::vector<std::vector<Real>> SPMOEA::sampleInLinkSpaces(std::vector<size_t>& front_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		int method = 1;
		std::vector<Real> front_sol_density;
		std::vector<size_t> front_sol_num;
		std::vector<Real> front_volume;
		for (size_t i = 0; i < front_spaces.size(); ++i) {
			auto sol_num=getMO_HLC().getSubspaceInfo(front_spaces[i]).m_subspace_front_sol.size();
			auto volume = getMO_HLC().subspaceTree().getBoxVolume(front_spaces[i]);
			front_volume.push_back(volume);
			front_sol_num.push_back(sol_num);
			front_sol_density.push_back(sol_num/volume);
		}
		for (size_t i = 0; i < sample_num; ++i) {
			//ѡ��һ���ӿռ�
			size_t sele_inx = 0;
			if (method == 1) {//���ѡ��
				sele_inx = (size_t)std::floor(front_spaces.size() * rnd->uniform.next());
			}
			else if(method==2) {//�����ӿռ��ڵ�ǰ����ѡ��
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_spaces[inx]);
			}
			//�ҳ���ǰ�������ӿռ�
			std::vector<size_t> nei1;
			auto neigh = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_sub_neighbors;
			for (auto jj : neigh) {
				if (std::find(front_spaces.begin(), front_spaces.end(), jj) != front_spaces.end()) {
					nei1.push_back(jj);
				}
			}
			std::vector<std::vector<Real>> off;

			auto& front_sols = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_subspace_front_sol;
			auto& his_sols = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_history_inds;
			if (his_sols.size() < 2) {
				auto off1 = sampleInSpace(front_spaces[sele_inx], bound, 1, pro, alg,rnd);
				for (auto jj : off1) {
					off.emplace_back(jj);
				}
			}
			else {
				//��һ���ӿռ�ǰ�ؽ�Ϊ���㣬��һ��֧���Ϊ��������������½�
				size_t inx1 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
				//�����Ͻ��ĸ���
				auto base_sol = front_sols[inx1]->variable().vect();
				Solution<> ind2(*front_sols[inx1]);
				std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
				Real min_span = INT16_MAX;
				auto box = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
				for (size_t j = 0; j < box.size(); ++j) {
					if (min_span > (box[j].second - box[j].first)) {
						min_span = box[j].second - box[j].first;
					}
				}
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto p2 = his_sols[j]->variable().vect();
					auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
					if (dist <= 0.) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				std::vector<size_t> candidate_inx;
				for (size_t j = 0; j < temp_dist.size(); ++j) {
					if (temp_dist[j] < min_span) {
						candidate_inx.push_back(j);
					}
				}
				//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
				size_t inx2 = 0;
				std::vector<Real> sol2;
				if (candidate_inx.size() > 0) {
					inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
					sol2 = his_sols[inx2]->variable().vect();
				}
				else {
					do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						sol2 = his_sols[inx2]->variable().vect();
					} while (ifSame(base_sol,sol2));
				}
				std::vector<Real> delta_v;
				Solution<> ind1(*his_sols[inx2]);
				//�Ƚ��������֧���ϵ
				auto dominance_ship = objectiveCompare(front_sols[inx1]->objective(), his_sols[inx2]->objective(), CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(2 * (sol2[j] - base_sol[j]));
					}
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real coff = 1.;
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						coff *= -1;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(coff * (base_sol[j] - sol2[j]));
					}
				}
				else if (dominance_ship == Dominance::kEqual) {
					size_t a = 1;
					//�������ֲ��Ŷ�
					auto sol = base_sol;
					auto box = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
					//��ÿһά������ߵ�ֵ
					std::vector<Real> min_v;
					for (size_t j = 0; j < box.size(); ++j) {
						auto v1 = std::fabs(sol[j] - box[j].first);
						auto v2 = std::fabs(sol[j] - box[j].second);
						min_v.push_back(v1 > v2 ? v2 : v1);
					}
					std::vector<Real> new_sol(box.size(), 0.);
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < box.size(); ++j) {
						dim_inx.push_back(j);
					}
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j]- new_sol[j]);
					}
				}
			    std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}

				Solution<> ind3(ind1);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro,alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				off.emplace_back(new_sol);
			}

			//auto& front_sols1 = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_front_sol_sub;
			//if (front_sols1.size() < 2) {
			//	//�������ֲ��Ŷ�
			//	auto sol = front_sols1[0]->variable().vect();
			//	auto box = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
			//	//��ÿһά������ߵ�ֵ
			//	std::vector<Real> min_v;
			//	for (size_t j = 0; j < box.size(); ++j) {
			//		auto v1 = std::fabs(sol[j] - box[j].first);
			//		auto v2 = std::fabs(sol[j] - box[j].second);
			//		min_v.push_back(v1 > v2 ? v2 : v1);
			//	}
			//	std::vector<Real> new_sol;
			//	bool flag = false;
			//	do {//ʹ�����ӿռ���
			//		new_sol.clear();
			//      flag = false;
			//		for (size_t j = 0; j < box.size(); ++j) {
			//			new_sol.push_back(rnd->normal.nextNonStd(sol[j], std::pow(min_v[j] / 6,2)));
			//		}
			//		for (size_t j = 0; j < new_sol.size(); ++j) {
			//			if (new_sol[j] < box[j].first || new_sol[j] > box[j].second) {
			//				flag = true;
			//				break;
			//			}
			//		}
			//	} while (flag);
			//	off.emplace_back(new_sol);
			//}
			//else {//GA�����������ӿռ佻��
			//	if (nei1.size() == 0) {//��ѡ�ӿռ��ڽ���
			//		SPMOEA_pop temp_pop(0, pro);
			//		for (size_t k = 0; k < front_sols1.size(); ++k) {
			//			temp_pop.append(*front_sols1[k]);
			//		}
			//		for (size_t j = 0; j < temp_pop.size(); ++j) {
			//			temp_pop.getOffspring().append(temp_pop[j]);
			//		}
			//		for (size_t j = 0; j < temp_pop.size(); ++j) {
			//			temp_pop.getOffspring().append(temp_pop[j]);
			//		}
			//		size_t kk = 2;
			//		auto search_bound = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
			//		off = sampleByGA(temp_pop, search_bound, 1,pro, rnd);
			//	}
			//	else {//��ѡ�ӿռ��������ӿռ佻��
			//		//���ѡ��һ���ӿռ�
			//		size_t inx1 = (size_t)std::floor(front_sols1.size() * rnd->uniform.next());
			//		size_t inx2= (size_t)std::floor(nei1.size() * rnd->uniform.next());
			//		auto& front_sols2 = getMO_HLC().getSubspaceInfo(nei1[inx2]).m_front_sol_sub;
			//		size_t inx3 = (size_t)std::floor(front_sols2.size() * rnd->uniform.next());
			//		SPMOEA_pop temp_pop(0, pro);
			//		temp_pop.append(*front_sols1[inx1]);
			//		temp_pop.append(*front_sols2[inx3]);
			//		for (size_t j = 0; j < temp_pop.size(); ++j) {
			//			temp_pop.getOffspring().append(temp_pop[j]);
			//		}
			//		for (size_t j = 0; j < temp_pop.size(); ++j) {
			//			temp_pop.getOffspring().append(temp_pop[j]);
			//		}
			//		size_t kk = 2;
			//		auto search_bound = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
			//		off = sampleByGA(temp_pop, search_bound, 1, pro, rnd);
			//	}
			//}

			//nei1.push_back(front_spaces[sele_inx]);
			////����Щ�������һ����Ⱥ�������½�
			//SPMOEA_pop temp_pop(0, pro);
			//for (size_t j = 0; j < nei1.size(); ++j) {
			//	auto& front_sols = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_sub;
			//	for (size_t k = 0; k < front_sols.size(); ++k) {
			//		temp_pop.append(*front_sols[k]);
			//	}
			//}
			//for (size_t j = 0; j < temp_pop.size(); ++j) {
			//	temp_pop.getOffspring().append(temp_pop[j]);
			//}
			//for (size_t j = 0; j < temp_pop.size(); ++j) {
			//	temp_pop.getOffspring().append(temp_pop[j]);
			//}
			//std::vector<std::vector<Real>> off;
			//if (temp_pop.size() < 2) {
			//	//�������ֲ��Ŷ�
			//	auto sol = temp_pop[0].variable().vect();
			//	auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			//	auto box = getMO_HLC().subspaceTree().getBox(space_inx);
			//	//��ÿһά������ߵ�ֵ
			//	std::vector<Real> min_v;
			//	for (size_t j = 0; j < box.size(); ++j) {
			//		auto v1 = std::fabs(sol[j] - box[j].first);
			//		auto v2 = std::fabs(sol[j] - box[j].second);
			//		min_v.push_back(v1 > v2 ? v2 : v1);
			//	}
			//	std::vector<Real> new_sol;
			//	bool flag = false;
			//	do {//ʹ�����ӿռ���
			//		new_sol.clear();
			//      flag = false;
			//		for (size_t j = 0; j < box.size(); ++j) {
			//			new_sol.push_back(rnd->normal.nextNonStd(sol[j], std::pow(min_v[j] / 6,2)));
			//		}
			//		for (size_t j = 0; j < new_sol.size(); ++j) {
			//			if (new_sol[j] < box[j].first|| new_sol[j] > box[j].second) {
			//				flag = true;
			//				break;
			//			}
			//		}
			//	} while (flag);
			//	off.emplace_back(new_sol);
			//}
			//else {
			//	size_t kk = 2;
			//	auto search_bound = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
			//	off = sampleByDE(temp_pop,search_bound, 1,2, pro, rnd);
			//}

			for (auto jj : off) {
				out_off.emplace_back(jj);
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA::sampleInSpace(size_t parent_space, size_t sample_num, Problem *pro, Algorithm *alg,Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto &box = getMO_HLC().subspaceTree().getBox(parent_space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& front_sols = getMO_HLC().getSubspaceInfo(parent_space).m_subspace_front_sol;
			auto& his_sols = getMO_HLC().getSubspaceInfo(parent_space).m_history_inds;
			if (his_sols.size() < 3) {
				//�������ֲ��Ŷ�
				auto off = sampleInRange(*front_sols[0], 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {
				//��һ���ӿռ�ǰ�ؽ�Ϊ���㣬��һ��֧���Ϊ��������������½�
				size_t inx1 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
				//�����Ͻ��ĸ���
				auto base_sol = front_sols[inx1]->variable().vect();
				Solution<> ind2(*front_sols[inx1]);
				std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
				Real min_span = INT16_MAX;
				for (size_t j = 0; j < box.size(); ++j) {
					if (min_span > (box[j].second - box[j].first)) {
						min_span = box[j].second - box[j].first;
					}
				}
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto p2 = his_sols[j]->variable().vect();
					auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
					if (dist <= 0.) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				std::vector<size_t> candidate_inx;
				for (size_t j = 0; j < temp_dist.size(); ++j) {
					if (temp_dist[j] < min_span) {
						candidate_inx.push_back(j);
					}
				}
				//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
				size_t inx2 = 0;
				std::vector<Real> sol2;
				if (candidate_inx.size() > 0) {
					inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
					sol2 = his_sols[inx2]->variable().vect();
				}
				else {
					do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						sol2 = his_sols[inx2]->variable().vect();
					} while (ifSame(base_sol, sol2));
				}
				std::vector<Real> delta_v;
				Solution<> ind1(*his_sols[inx2]);
				//�Ƚ��������֧���ϵ
				auto dominance_ship = objectiveCompare(front_sols[inx1]->objective(), his_sols[inx2]->objective(), CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(2 * (sol2[j] - base_sol[j]));
					}
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real coff = 1.;
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						coff *= -1;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(coff * (base_sol[j] - sol2[j]));
					}
				}
				else if (dominance_ship == Dominance::kEqual) {
					//�������ֲ��Ŷ�
					auto sol = base_sol;
					//��ÿһά������ߵ�ֵ
					std::vector<Real> min_v;
					for (size_t j = 0; j < box.size(); ++j) {
						auto v1 = std::fabs(sol[j] - box[j].first);
						auto v2 = std::fabs(sol[j] - box[j].second);
						min_v.push_back(v1 > v2 ? v2 : v1);
					}
					std::vector<Real> new_sol(box.size(), 0.);
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < box.size(); ++j) {
						dim_inx.push_back(j);
					}
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - new_sol[j]);
					}
				}
				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}
				/*for (size_t j = 0; j < new_sol.size(); ++j) {
					auto temp = new_sol[j];
					new_sol[j] += ((2 * rand - 1) * (raodong / 2) * delta_v[j]);
					if (new_sol[j] < bound[j].first || new_sol[j] > bound[j].second) {
						new_sol[j] = temp;
					}
				}*/

				Solution<> ind3(ind1);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				out_off.emplace_back(new_sol);
			}
		}
		return out_off;
	}

	//ǰ���ӿռ����
	std::vector<std::vector<Real>> SPMOEA::sampleInFrontSpaces(std::vector<size_t>& front_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg,Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		int method = 1;
		std::vector<Real> front_sol_density;
		std::vector<size_t> front_sol_num;
		std::vector<Real> front_volume;
		for (size_t i = 0; i < front_spaces.size(); ++i) {
			auto sol_num = getMO_HLC().getSubspaceInfo(front_spaces[i]).m_subspace_front_sol.size();
			auto volume = getMO_HLC().subspaceTree().getBoxVolume(front_spaces[i]);
			front_volume.push_back(volume);
			front_sol_num.push_back(sol_num);
			front_sol_density.push_back(sol_num / volume);
		}
		for (size_t i = 0; i < sample_num; ++i) {
			//ѡ��һ���ӿռ�
			size_t sele_inx = 0;
			if (method == 1) {//���ѡ��
				sele_inx = (size_t)std::floor(front_spaces.size() * rnd->uniform.next());
			}
			else if (method == 2) {//�����ӿռ��ڵ�ǰ����ѡ��
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_spaces[inx]);
			}
			//�ҳ���ǰ�������ӿռ�
			std::vector<size_t> nei1;
			auto neigh = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_sub_neighbors;
			for (auto jj : neigh) {
				if (std::find(front_spaces.begin(), front_spaces.end(), jj) != front_spaces.end()) {
					nei1.push_back(jj);
				}
			}
			std::vector<std::vector<Real>> off;

			auto& front_sols = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_subspace_front_sol;
			auto& his_sols = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_history_inds;
			if (his_sols.size() < 2) {
				//�������ֲ��Ŷ�
				auto off1 = sampleInRange(*his_sols[0], 1, pro, alg, rnd);
				for (auto ii : off1) {
					off.emplace_back(ii);
				}
			}
			else {
				//��һ���ӿռ�ǰ�ؽ�Ϊ���㣬��һ��֧���Ϊ��������������½�
				size_t inx1 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
				//�����Ͻ��ĸ���
				auto base_sol = front_sols[inx1]->variable().vect();
				Solution<> ind2(*front_sols[inx1]);
				std::vector<Real> temp_dist;//����С���ӿռ���С��ȵĵ�
				Real min_span = INT16_MAX;
				auto box = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
				for (size_t j = 0; j < box.size(); ++j) {
					if (min_span > (box[j].second - box[j].first)) {
						min_span = box[j].second - box[j].first;
					}
				}
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto p2 = his_sols[j]->variable().vect();
					auto dist = euclideanDistance(base_sol.begin(), base_sol.end(), p2.begin());
					if (dist <= 0.) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				std::vector<size_t> candidate_inx;
				for (size_t j = 0; j < temp_dist.size(); ++j) {
					if (temp_dist[j] < min_span) {
						candidate_inx.push_back(j);
					}
				}
				//���ѡȡһ����Χ�ڵĵ�Ƚ�֧���ϵ
				size_t inx2 = 0;
				std::vector<Real> sol2;
				if (candidate_inx.size() > 0) {
					inx2 = candidate_inx[(size_t)std::floor(candidate_inx.size() * rnd->uniform.next())];
					sol2 = his_sols[inx2]->variable().vect();
				}
				else {
					do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						sol2 = his_sols[inx2]->variable().vect();
					} while (ifSame(base_sol, sol2));
				}
				std::vector<Real> delta_v;
				Solution<> ind1(*his_sols[inx2]);
				//�Ƚ��������֧���ϵ
				auto dominance_ship = objectiveCompare(front_sols[inx1]->objective(), his_sols[inx2]->objective(), CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(2 * (sol2[j] - base_sol[j]));
					}
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real coff = 1.;
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						coff *= -1;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(coff * (base_sol[j] - sol2[j]));
					}
				}
				else if (dominance_ship == Dominance::kEqual) {
					size_t a = 1;
					//�������ֲ��Ŷ�
					auto sol = base_sol;
					auto box = getMO_HLC().subspaceTree().getBox(front_spaces[sele_inx]);
					//��ÿһά������ߵ�ֵ
					std::vector<Real> min_v;
					for (size_t j = 0; j < box.size(); ++j) {
						auto v1 = std::fabs(sol[j] - box[j].first);
						auto v2 = std::fabs(sol[j] - box[j].second);
						min_v.push_back(v1 > v2 ? v2 : v1);
					}
					std::vector<Real> new_sol(box.size(), 0.);
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < box.size(); ++j) {
						dim_inx.push_back(j);
					}
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = rnd->normal.nextNonStd(sol[dim_inx[j]], std::pow(min_v[dim_inx[j]], 2));
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j] - new_sol[j]);
					}
				}
				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}

				Solution<> ind3(ind2);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				off.emplace_back(new_sol);
			}
			for (auto jj : off) {
				out_off.emplace_back(jj);
			}
		}
		return out_off;
	}

	//��ͨǰ���ӿռ��Ե��չ����
	std::vector<std::vector<Real>> SPMOEA::sampleByExtend(size_t space1, size_t space2, size_t off_num, Problem *pro, Random *rnd) {
		//���ӿռ��ڵ�ǰ�ؽ������Ԥ�ⷽ��,space1Ϊ��ʼ�ӿռ�
		std::vector<std::vector<Real>> all_off;
		auto& front_sols1 = getMO_HLC().getSubspaceInfo(space1).m_subspace_front_sol;//����Ƿ����
		auto& front_sols2 = getMO_HLC().getSubspaceInfo(space2).m_subspace_front_sol;//����Ƿ����
		std::vector<Real> center1;
		std::vector<Real> center2;
		size_t dim = CAST_CONOP(pro)->numberVariables();
		for (size_t i = 0; i < dim; ++i) {
			Real sum = 0.;
			for (size_t j = 0; j < front_sols1.size(); ++j) {
				sum += (front_sols1[j]->variable().vect()[i]);
			}
			center1.push_back(sum / (Real)front_sols1.size());
		}
		for (size_t i = 0; i < dim; ++i) {
			Real sum = 0.;
			for (size_t j = 0; j < front_sols2.size(); ++j) {
				sum += (front_sols2[j]->variable().vect()[i]);
			}
			center2.push_back(sum / (Real)front_sols2.size());
		}
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < off_num; ++i) {
			auto sol = vectorPredict(center1, center2, bound, rnd);
			//Ԥ��������ӿռ�
			auto sp_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			//��Ԥ����ӿռ���Ԥ���Ϊ���Ĳ���
			std::vector<Real> off;
			for (size_t k = 0; k < sol.size(); ++k) {
				auto nd_rand = rnd->normal.nextNonStd(sol[k],std::abs(center1[k]-center2[k])/3);
				if (nd_rand < bound[k].first) {
					nd_rand = 2 * bound[k].first - nd_rand;
				}
				else if (nd_rand > bound[k].second) {
					nd_rand = 2 * bound[k].second - nd_rand;
				}
				off.push_back(nd_rand);
			}
			all_off.emplace_back(off);
		}
		return all_off;
	}

	//�ӿռ估�������ӿռ��������
	std::vector<std::vector<Real>> SPMOEA::sampleByNeighSpace(size_t space1, size_t off_num, Problem *pro, Algorithm* alg, Random *rnd) {
		//����ѡ�ӿռ��ڵ������������ӿռ�������Ϊ�������н���������һ���������Ӵ�
		std::vector<std::vector<Real>> all_off;
		//auto& front_sols1 = getMO_HLC().getSubspaceInfo(space1).m_subspace_front_sol;//����Ƿ����
		std::vector<Real> p1;
		auto& box1 = getMO_HLC().subspaceTree().getBox(space1);
		for (size_t i = 0; i < box1.size(); ++i) {
			p1.push_back(box1[i].first+(box1[i].second-box1[i].first)*rnd->uniform.next());
		}
		//�������ӿռ������ѡ��һ���ӿռ�
		auto& neigh = getMO_HLC().getSubspaceInfo(space1).m_sub_neighbors;
		size_t inx = (size_t)std::floor(neigh.size()*rnd->uniform.next());
		size_t count = 0;
		std::vector<std::pair<Real, Real>> box2;
		for (auto ii : neigh) {
			if (count == inx) {
				box2 = getMO_HLC().subspaceTree().getBox(ii);
				break;
			}
			count++;
		}
		std::vector<Real> p2;
		for (size_t i = 0; i < box2.size(); ++i) {
			p2.push_back(box2[i].first + (box2[i].second - box2[i].first) * rnd->uniform.next());
		}
		//����λ�ý�������һ���������Ӵ�
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < off_num; ++i) {
			auto off=vectorOperator(p1,p2,0.1,pro,alg,rnd);
			all_off.emplace_back(off);
		}

		return all_off;
	}

	//�ӿռ估�������ӿռ��ǰ�ؽ������ʷ�⽻��
	std::vector<std::vector<Real>> SPMOEA::sampleByNeighSpaceFront(size_t space1,size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t off_num, Problem *pro,Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//�ҳ��������ӿռ�
		std::vector<size_t> nei1;
		nei1.push_back(space1);
		auto neigh = getMO_HLC().getSubspaceInfo(space1).m_sub_neighbors;
		for (auto jj : neigh) {
			nei1.push_back(jj);
		}
		//������ǰ�ظ��幹������Ⱥ
		auto& sol = getPop()[0][ind_inx].variable().vect();
		SPMOEA_pop temp_pop(0, pro);
		for (size_t j = 0; j < nei1.size(); ++j) {
			if (getMO_HLC().getSubspaceInfo(nei1[j]).m_history_inds.size()>0) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_subspace_front_sol;
				auto inx = (size_t)std::floor(front_sol.size() * rnd->uniform.next());
				if(nei1[j]==space1) {
					if (front_sol.size() > 1) {
						while (ifSame(sol, front_sol[inx]->variable().vect())) {
							inx = (size_t)std::floor(front_sol.size() * rnd->uniform.next());
						}
					}
				}
				temp_pop.append(*front_sol[inx]);
			}
		}
		
		temp_pop.append(getPop()[0][ind_inx]);
		for (size_t i = 0; i < off_num; ++i) {
			if (temp_pop.size() < 3) {
				auto off1 = sampleInRange(temp_pop[temp_pop.size()-1], 1, pro, alg, rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				size_t kk = 2;
				auto off = sampleByDE(temp_pop, temp_pop.size() - 1,bound, 1, kk, pro, alg, rnd);
				for (size_t k = 0; k < off.size(); ++k) {
					out_off.emplace_back(off[k]);
					if (ifSame(sol, off[k])) {
						size_t a = 1;
					}
				}
			}
		}

		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA::sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		int method = 1;
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		size_t con_num = CAST_CONOP(pro)->numberConstraints();
		std::vector<Real> front_sol_density;
		std::vector<size_t> front_sol_num;
		std::vector<Real> front_volume;
		for (size_t i = 0; i < front_link_spaces.size(); ++i) {
			auto sol_num = getMO_HLC().getSubspaceInfo(front_link_spaces[i]).m_history_inds.size();
			auto volume = getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[i]);
			front_volume.push_back(volume);
			front_sol_num.push_back(sol_num);
			front_sol_density.push_back(sol_num / volume);
		}
		for (size_t i = 0; i < sample_num; ++i) {
			//ѡ��һ���ӿռ�
			size_t sele_inx = 0;
			if (method == 1) {//���ѡ��
				sele_inx = (size_t)std::floor(front_link_spaces.size() * rnd->uniform.next());
			}
			else if (method == 2) {//�����ӿռ��ڵ�ǰ����ѡ��
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[inx]);
			}
			//�ҳ���ǰ�������ӿռ�
			std::vector<size_t> nei1;
			auto neigh = getMO_HLC().getSubspaceInfo(front_link_spaces[sele_inx]).m_sub_neighbors;
			for (auto jj : neigh) {
				if (std::find(front_link_spaces.begin(), front_link_spaces.end(), jj) != front_link_spaces.end()) {
					nei1.push_back(jj);
				}
			}
			nei1.push_back(front_link_spaces[sele_inx]);

			std::vector<std::vector<Real>> off;
			if (nei1.size()<2) {//�ӿռ��ڲ�����
				auto off1 = sampleInSpace(front_link_spaces[sele_inx],1,pro,alg,rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
				}
			}
			else {
				auto& front_sols1 = getMO_HLC().getSubspaceInfo(front_link_spaces[sele_inx]).m_front_sol_in_subspace;
				size_t inx1= (size_t)std::floor(front_sols1.size() * rnd->uniform.next());
				auto &base_sol = front_sols1[inx1]->variable().vect();
				Solution<> ind2(*front_sols1[inx1]);
				std::vector<std::vector<Real>> all_sols;
				std::vector<std::vector<Real>> all_objs;
				for (size_t j = 0; j < nei1.size(); ++j) {
					auto &front_sol2= getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
					for (size_t k = 0; k < front_sol2.size(); ++k) {
						auto& sol = front_sol2[k]->variable().vect();
						all_sols.emplace_back(sol);
						auto& obj = front_sol2[k]->objective();
						all_objs.emplace_back(obj);
					}
				}
				//ѡȡһ�����ظ���
				std::vector<Real> sol2;
				size_t inx2;
				do {
					inx2 = (size_t)std::floor(all_sols.size() * rnd->uniform.next());
					sol2 = all_sols[inx2];
				} while (ifSame(base_sol,sol2));
				//�Ƚ�����֧���ϵ
				Solution<> ind1(ind2);
				ind1.variable().vect() = all_sols[inx2];
				ind1.objective() = all_objs[inx2];
				std::vector<Real> delta_v;
				auto dominance_ship = objectiveCompare(front_sols1[inx1]->objective(),all_objs[inx2],CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(base_sol[j]-sol2[j]);
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(2*(sol2[j]- base_sol[j]));
					}
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real coff = 1.;
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						coff *= -1;
					}
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(coff*(base_sol[j] - sol2[j]));
					}
				}
				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				Real raodong = ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				auto rand = rnd->uniform.next();
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						rand /= 1.2;
					}
				}
				/*for (size_t j = 0; j < new_sol.size(); ++j) {
					auto temp = new_sol[j];
					new_sol[j] += ((2 * rand - 1) * (raodong / 2) * delta_sol[j]);
					if (new_sol[j] < bound[j].first || new_sol[j] > bound[j].second) {
						new_sol[j] = temp;
					}
				}*/

				Solution<> ind3(ind2);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);
				
				out_off.emplace_back(new_sol);
			}
		}
		return out_off;
	}

	//ǰ���ӿռ���������ѡ���彻��
	std::vector<std::vector<Real>> SPMOEA::sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces,std::vector<Real>& sol, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
		//�ҳ���ǰ�������ӿռ�
		std::vector<size_t> nei1;
		auto neigh = getMO_HLC().getSubspaceInfo(space_inx).m_sub_neighbors;
		for (auto jj : neigh) {
			if (std::find(front_link_spaces.begin(), front_link_spaces.end(), jj) != front_link_spaces.end()) {
				nei1.push_back(jj);
			}
		}
		nei1.push_back(space_inx);
		//������ǰ�ظ��幹������Ⱥ
		size_t pop_size = 0;
		for (size_t j = 0; j < nei1.size(); ++j) {
			auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_subspace_front_sol;
			pop_size += front_sol.size();
		}
		SPMOEA_pop temp_pop(pop_size, pro);
		size_t count = 0;
		for (size_t j = 0; j < nei1.size(); ++j) {
			auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_subspace_front_sol;
			for (size_t k = 0; k < front_sol.size(); ++k) {
				temp_pop[count].variable() = front_sol[k]->variable();
				temp_pop[count].objective() = front_sol[k]->objective();
				count++;
			}
		}
		bool in_pop = false;
		size_t sol_inx = 0;
		for (size_t j = 0; j < pop_size; ++j) {
			if (ifSame(sol, temp_pop[j].variable().vect())) {
				in_pop = true;
				sol_inx = j;
				break;
			}
		}
		if (!in_pop) {
			temp_pop.append(temp_pop[0]);
			temp_pop.back().variable() = sol;
			sol_inx = pop_size;
		}
		if (pop_size < 5) {
			auto off1 = sampleInRange(temp_pop[sol_inx], sample_num, pro, alg, rnd);
			for (auto jj : off1) {
				out_off.emplace_back(jj);
				if (ifSame(sol,jj)) {
					size_t a = 1;
				}
			}
		}
		else {
			size_t kk = 2;
			for (size_t i = 0; i < sample_num; ++i) {
				auto off = sampleByDE(temp_pop, sol_inx, bound, 1, kk, pro, alg, rnd);
				for (size_t k = 0; k < off.size(); ++k) {
					out_off.emplace_back(off[k]);
					if (ifSame(sol, off[k])) {
						size_t a = 1;
					}
				}
			}
		}
		return out_off;
	}

	//�ӿռ��ǰ���ӿռ����򽻻��γ�һ���⵱���Ӵ�
	std::vector<std::vector<Real>> SPMOEA::sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		int method = 1;
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		size_t con_num = CAST_CONOP(pro)->numberConstraints();
		std::vector<Real> front_sol_density;
		std::vector<size_t> front_sol_num;
		std::vector<Real> front_volume;
		for (size_t i = 0; i < front_link_spaces.size(); ++i) {
			auto sol_num = getMO_HLC().getSubspaceInfo(front_link_spaces[i]).m_history_inds.size();
			auto volume = getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[i]);
			front_volume.push_back(volume);
			front_sol_num.push_back(sol_num);
			front_sol_density.push_back(sol_num / volume);
		}
		auto& sol = getPop()[0][ind_inx].variable().vect();
		for (size_t i = 0; i < sample_num; ++i) {
			//ѡ��һ���ӿռ�
			size_t sele_inx = 0;
			if (method == 1) {//���ѡ��
				sele_inx = (size_t)std::floor(front_link_spaces.size() * rnd->uniform.next());
			}
			else if (method == 2) {//�����ӿռ��ڵ�ǰ����ѡ��
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				//front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[inx]);
			}
			//�ҳ���ǰ�������ӿռ�
			std::vector<size_t> nei1;
			/*nei1.push_back(front_link_spaces[sele_inx]);
			auto neigh = getMO_HLC().getSubspaceInfo(front_link_spaces[sele_inx]).m_sub_neighbors;
			for (auto jj : neigh) {
				if (std::find(front_link_spaces.begin(), front_link_spaces.end(), jj) != front_link_spaces.end()) {
					nei1.push_back(jj);
				}
			}*/
			for (auto jj : front_link_spaces) {
				nei1.push_back(jj);
			}
			//������ǰ�ظ��幹������Ⱥ
			size_t pop_size = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				pop_size += front_sol.size();
			}
			PopDE<> temp_pop(pop_size, pro);
			size_t count = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				for (size_t k = 0; k < front_sol.size(); ++k) {
					temp_pop[count].variable() = front_sol[k]->variable();
					temp_pop[count].objective() = front_sol[k]->objective();
					count++;
				}
			}
			if (pop_size < 5) {
				size_t inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				auto off1 = sampleInRange(temp_pop[inx], sample_num, pro, alg, rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				std::vector<size_t> candidate_inx;
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				size_t inx = 0;
				do {
					inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				} while (inx == candidate_inx[1]);
				candidate_inx.push_back(inx);
				temp_pop[candidate_inx[0]].mutate(1.2*rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);
				temp_pop[candidate_inx[0]].trial().evaluate(pro, this);
				Solution<> ind4(temp_pop[candidate_inx[0]].trial());
				/*ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();*/
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				m_interactive_sol_pair.emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
				
				/*size_t kk = 2;
				auto off = sampleByDE(temp_pop, bound, 1, kk, pro, alg, rnd);
				for (size_t k = 0; k < off.size(); ++k) {
					out_off.emplace_back(off[k]);
					if (ifSame(sol, off[k])) {
						size_t a = 1;
					}
				}*/
			}
		}
		return out_off;
	}

	//ǰ���ӿռ���ͨ������
	std::vector<std::vector<Real>> SPMOEA::sampleInFrontNeighSpace2(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		int method = 1;
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		size_t obj_num = CAST_CONOP(pro)->numberObjectives();
		size_t con_num = CAST_CONOP(pro)->numberConstraints();
		std::vector<Real> front_sol_density;
		std::vector<size_t> front_sol_num;
		std::vector<Real> front_volume;
		for (size_t i = 0; i < front_link_spaces.size(); ++i) {
			auto sol_num = getMO_HLC().getSubspaceInfo(front_link_spaces[i]).m_history_inds.size();
			auto volume = getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[i]);
			front_volume.push_back(volume);
			front_sol_num.push_back(sol_num);
			front_sol_density.push_back(sol_num / volume);
		}
		auto& sol = getPop()[0][ind_inx].variable().vect();
		for (size_t i = 0; i < sample_num; ++i) {
			//ѡ��һ���ӿռ�
			size_t sele_inx = 0;
			if (method == 1) {//���ѡ��
				sele_inx = (size_t)std::floor(front_link_spaces.size() * rnd->uniform.next());
			}
			else if (method == 2) {//�����ӿռ��ڵ�ǰ����ѡ��
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				//front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[inx]);
			}
			//�ҳ���ǰ�������ӿռ�
			std::vector<size_t> nei1;
			nei1.push_back(front_link_spaces[sele_inx]);
			auto neigh = getMO_HLC().getSubspaceInfo(front_link_spaces[sele_inx]).m_sub_neighbors;
			for (auto jj : neigh) {
				if (std::find(front_link_spaces.begin(), front_link_spaces.end(), jj) != front_link_spaces.end()) {
					nei1.push_back(jj);
				}
			}
			//������ǰ�ظ��幹������Ⱥ
			size_t pop_size = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				pop_size += front_sol.size();
			}
			SPMOEA_pop temp_pop(pop_size, pro);
			size_t count = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				for (size_t k = 0; k < front_sol.size(); ++k) {
					temp_pop[count].variable() = front_sol[k]->variable();
					temp_pop[count].objective() = front_sol[k]->objective();
					count++;
				}
			}
			if (pop_size < 5) {
				size_t inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				auto off1 = sampleInRange(temp_pop[inx], sample_num, pro, alg, rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				size_t kk = 2;
				auto off = sampleByDE(temp_pop, bound, 1, kk, pro, alg, rnd);
				for (size_t k = 0; k < off.size(); ++k) {
					out_off.emplace_back(off[k]);
					if (ifSame(sol, off[k])) {
						size_t a = 1;
					}
				}
			}
		}
		return out_off;
	}

	//ǰ���ӿռ估�������ӿռ佻��
	std::vector<std::vector<Real>> SPMOEA::sampleByFrontNeighSpace(size_t space1, std::vector<size_t>& front_link_spaces, std::vector<std::vector<size_t>>& manifold_dist, size_t sample_num, Problem *pro, Algorithm* alg, Random *rnd) {
		//����ѡ�ӿռ��ڵ������������ӿռ�������Ϊ�������н���������һ���������Ӵ�
		std::vector<std::vector<Real>> all_off;
		auto search_bound = CAST_CONOP(pro)->boundary();
		size_t index1 = front_link_spaces[space1];
		//�ҳ������ؾ���Ϊ1���ӿռ�
	    //��ȡ��Щ�ӿռ��ǰ�ؽ�
		std::vector<std::vector<Real>> all_sols;
		std::vector<size_t> nei1;
		for (size_t k = 0; k < manifold_dist[space1].size(); ++k) {
			if (manifold_dist[space1][k] == 1) {
				for (size_t p = 0; p < getMO_HLC().getSubspaceInfo(front_link_spaces[k]).m_front_sol_in_subspace.size(); ++p) {
					all_sols.emplace_back(getMO_HLC().getSubspaceInfo(front_link_spaces[k]).m_front_sol_in_subspace[p]->variable().vect());
				}
				nei1.push_back(front_link_spaces[k]);
			}
		}
		auto& pp = getMO_HLC().getSubspaceInfo(index1).m_front_sol_in_subspace;
		for (size_t i = 0; i < sample_num; ++i) {
			//���ѡ��ѡ���ӿռ��е�һ��ǰ�ؽ���Ϊ����
			auto p1 = (size_t)std::floor(pp.size() * rnd->uniform.next());
			//���ѡ�������һ����Ϊ��������������������ʽ�����Ӵ�
			auto p2 = (size_t)std::floor(all_sols.size() * rnd->uniform.next());
			auto off = vectorOperator(pp[p1]->variable().vect(), all_sols[p2],0.1,pro,alg, rnd);
			all_off.emplace_back(off);
		}

		return all_off;
	}

	//ǰ���ӿռ估�������ӿռ佻��
	std::vector<std::vector<Real>> SPMOEA::sampleByFrontNeighSpace(size_t space1,size_t ind_inx, std::vector<size_t>& front_link_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		//�ҳ���ǰ�������ӿռ�
		std::vector<size_t> nei1;
		nei1.push_back(space1);
		auto neigh = getMO_HLC().getSubspaceInfo(space1).m_sub_neighbors;
		for (auto jj : neigh) {
			if (std::find(front_link_spaces.begin(), front_link_spaces.end(), space1) != front_link_spaces.end()) {
				nei1.push_back(jj);
			}
		}
		auto& his_sols = getMO_HLC().getSubspaceInfo(space1).m_history_inds;
		auto& front_sols = getMO_HLC().getSubspaceInfo(space1).m_front_sol_in_subspace;
		auto& sol = getPop()[0][ind_inx].variable().vect();
		//�Ը���ǰ�ƻ���չ
		Real rand = rnd->uniform.next();
		if (rand > 0.5) {//ǰ��
			if (his_sols.size() < 5) {
				auto off1 = sampleInRange(getPop()[0][ind_inx], 1, pro, alg, rnd);
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//Ѱ��֧�����������Ҳ��������ҷ�֧������
				size_t inx1 = (size_t)std::floor(front_sols.size() * rnd->uniform.next());
				std::vector<size_t> dominant_inx;
				std::vector<size_t> non_dominant_inx;
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto dominance_ship = objectiveCompare(front_sols[inx1]->objective(), his_sols[j]->objective(), CAST_CONOP(pro)->optimizeMode());
					if (dominance_ship == Dominance::kDominant) {
						dominant_inx.push_back(j);
					}
				}
				if (dominant_inx.size() < 1) {
					SPMOEA_pop temp_pop(0, pro);
					for (size_t k = 0; k < his_sols.size(); ++k) {
						temp_pop.append(*his_sols[k]);
					}
					temp_pop.append(getPop()[0][ind_inx]);
					size_t kk = 2;
					auto off = sampleByDE(temp_pop, temp_pop.size() - 1, bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						out_off.emplace_back(off[k]);
						if (ifSame(sol, off[k])) {
							size_t a = 1;
						}
					}
				}
				else {
					//ѡ�Ͻ��ĸ��彻��
					Solution<> ind1(getPop()[0][ind_inx]);
					size_t inx2 = dominant_inx[(size_t)std::floor(dominant_inx.size() * rnd->uniform.next())];
					Solution<> ind2(*his_sols[inx2]);
					Solution<> ind3(*front_sols[inx1]);
					std::vector<Real> sol1 = ind2.variable().vect();
					std::vector<Real> sol2 = ind3.variable().vect();
					std::vector<Real> new_sol(sol1.size());
					std::vector<Real> delta_v(sol1.size());
					for (size_t j = 0; j < delta_v.size(); ++j) {
						delta_v[j] = 1.5*(sol2[j] - sol1[j]);
					}
					std::vector<size_t> dim_inx;
					for (size_t j = 0; j < delta_v.size(); ++j) {
						dim_inx.push_back(j);
					}
					Real ratio = rnd->uniform.next();
					while (dim_inx.size() > 0) {
						std::vector<size_t> temp;
						for (size_t j = 0; j < dim_inx.size(); ++j) {
							new_sol[dim_inx[j]] = sol1[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
							if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
								temp.push_back(dim_inx[j]);
							}
						}
						dim_inx = temp;
						if (dim_inx.size() > 0) {
							ratio /= 1.2;
						}
					}
					auto mutate_sol = new_sol;
					for (size_t j = 0; j < mutate_sol.size(); ++j) {
						//ÿһά����ʼ��ľ�����Сֵ�����Ŷ���Χ
						Real dist1 = std::fabs(new_sol[j] - sol1[j]);
						Real dist2 = std::fabs(new_sol[j] - sol2[j]);
						auto amplitude = dist1 < dist2 ? dist1 : dist2;
						auto rand = rnd->uniform.next();
						mutate_sol[j] += ((2 * rand - 1) * amplitude);
					}
					//��Խ�紦��
					for (size_t j = 0; j < mutate_sol.size(); ++j) {
						if (mutate_sol[j] < bound[j].first) {
							mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
						}
						else if (mutate_sol[j] > bound[j].second) {
							mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
						}
					}
					if (ifSame(sol, mutate_sol)) {
						size_t a = 1;
					}

					Solution<> ind4(ind2);
					ind4.variable().vect() = mutate_sol;
					ind4.evaluate(pro, alg);
					std::vector<std::shared_ptr<Solution<>>> temp_pair;
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
					m_interactive_sol_pair.emplace_back(temp_pair);

					out_off.emplace_back(mutate_sol);
				}
			}
		}
		else {//��չ
			//������ǰ�ظ��幹������Ⱥ
			SPMOEA_pop temp_pop(0, pro);
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				for (size_t k = 0; k < front_sol.size(); ++k) {
					temp_pop.append(*front_sol[k]);
				}
			}

			temp_pop.append(getPop()[0][ind_inx]);
			for (size_t i = 0; i < sample_num; ++i) {
				if (temp_pop.size() < 3) {
					auto off1 = sampleInRange(temp_pop[temp_pop.size() - 1], 1, pro, alg, rnd);
					for (auto jj : off1) {
						out_off.emplace_back(jj);
						if (ifSame(sol, jj)) {
							size_t a = 1;
						}
					}
				}
				else {
					size_t kk = 2;
					auto off = sampleByDE(temp_pop, temp_pop.size() - 1, bound, 1, kk, pro, alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						out_off.emplace_back(off[k]);
						if (ifSame(sol, off[k])) {
							size_t a = 1;
						}
					}
				}
			}
		}
		return out_off;
	}

	//��ͨ�ӿռ���չ������������
	std::vector<std::vector<Real>> SPMOEA::sampleByLinkSpaces(std::vector<size_t>& front_spaces, std::vector<std::vector<size_t>>& manifold_dist, size_t sample_num,Problem *pro, Algorithm* alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		////����ͨ�ռ��ڲ���
		//for (size_t i = 0; i < sample_num; ++i) {
		//	//ѡ��һ���ӿռ�
		//	size_t sele_inx = (size_t)std::floor(front_spaces.size() * rnd->uniform.next());
		//	//�ҳ��������ξ�����С�������ӿռ�
		//	std::vector<size_t> nei1;
		//	auto neigh = getMO_HLC().getSubspaceInfo(front_spaces[sele_inx]).m_sub_neighbors;
		//	for (auto jj:neigh) {
		//		if (std::find(front_spaces.begin(),front_spaces.end(),jj)!=front_spaces.end()) {
		//			nei1.push_back(jj);
		//		}
		//	}
		//	nei1.push_back(front_spaces[sele_inx]);
		//	//����Щ�������һ����Ⱥ�������½�
		//	SPMOEA_pop temp_pop(0,pro);
		//	for (size_t j = 0; j < nei1.size(); ++j) {
		//		auto& front_sols = getMO_HLC().getSubspaceInfo(nei1[j]).m_subspace_front_sol;
		//		temp_pop.append(*front_sols[j]);
		//	}
		//	std::vector<std::vector<Real>> off;
		//	if (temp_pop.size() < 2) {
		//		//�������ֲ��Ŷ�

		//	}
		//	else if(temp_pop.size() < 3) {
		//		off = sampleByGA(temp_pop, 1, pro, rnd);
		//	}
		//	else {
		//		off = sampleByDE(temp_pop, 1, pro, rnd);
		//	}
		//	for (auto jj : off) {
		//		out_off.emplace_back(jj);
		//	}
		//}

		//������ͨ�ӿռ���ԵĲ�ؾ������
		//�ҳ���Ե�����Զ���ӿռ䣨��չ��
		size_t inx1, inx2;
		std::vector<size_t> temp_max;
		for (size_t j = 0; j < manifold_dist.size(); ++j) {
			temp_max.push_back(*std::max_element(manifold_dist[j].begin(), manifold_dist[j].end()));
		}
		auto index1 = std::distance(temp_max.begin(), std::max_element(temp_max.begin(), temp_max.end()));
		auto index2 = std::distance(manifold_dist[index1].begin(), std::max_element(manifold_dist[index1].begin(), manifold_dist[index1].end()));
		inx1 = front_spaces[index1];
		inx2 = front_spaces[index2];
		
		////�������ӿռ佻��
		//auto off1=sampleByNeighSpace(inx1,sample_num/2,pro,rnd);
		//auto off2 = sampleByNeighSpace(inx2, sample_num/2, pro, rnd);
		//for (auto ii : off1) {
		//	out_off.emplace_back(ii);
		//}
		//for (auto jj : off2) {
		//	out_off.emplace_back(jj);
		//}

		//��ԵԤ��
		if (index1 == index2) {//�������ӿռ佻��
			//��ȡ�����ӿռ��ǰ�ؽ�
			auto off1=sampleByNeighSpace(inx1,sample_num/2,pro,alg,rnd);
			auto off2 = sampleByNeighSpace(inx2, sample_num / 2, pro,alg, rnd);
			for (auto ii : off1) {
				out_off.emplace_back(ii);
			}
			for (auto jj : off2) {
				out_off.emplace_back(jj);
			}
		}
		else {
			//���������
			//�ҳ����ξ�����С�������ӿռ�
			std::vector<size_t> nei1;
			for (size_t j = 0; j < manifold_dist[index1].size(); ++j) {
				if (manifold_dist[index1][j] == 1) {
					nei1.push_back(front_spaces[j]);
				}
			}
			std::vector<size_t> nei2;
			for (size_t j = 0; j < manifold_dist[index2].size(); ++j) {
				if (manifold_dist[index2][j] == 1) {
					nei2.push_back(front_spaces[j]);
				}
			}
			//ѡ��һ���ӿռ�
			size_t sele1 = (size_t)std::floor(nei1.size() * rnd->uniform.next());
			size_t sele2 = (size_t)std::floor(nei2.size() * rnd->uniform.next());
			auto off1 = sampleByExtend(nei1[sele1], inx1, sample_num / 2, pro, rnd);
			auto off2 = sampleByExtend(nei2[sele2], inx2, sample_num / 2, pro, rnd);
			for (auto jj : off1) {
				out_off.emplace_back(jj);
			}
			for (auto jj : off2) {
				out_off.emplace_back(jj);
			}
		}
		return out_off;
	}

	void SPMOEA::manifoldSpread(std::vector<size_t>& front_link_spaces,Problem *pro,Random *rnd) {
		//������ͨ�ӿռ���ԵĲ�ؾ������
		auto manifold_dist = calSpaceManifoldDist(front_link_spaces);//�Ż�����
		//�ҳ���Ե�����Զ���ӿռ䣨��չ��
		size_t inx1, inx2;
		std::vector<size_t> temp_max;
		for (size_t j = 0; j < manifold_dist.size(); ++j) {
			temp_max.push_back(*std::max_element(manifold_dist[j].begin(), manifold_dist[j].end()));
		}
		auto index1 = std::distance(temp_max.begin(), std::max_element(temp_max.begin(), temp_max.end()));
		auto index2 = std::distance(manifold_dist[index1].begin(), std::max_element(manifold_dist[index1].begin(), manifold_dist[index1].end()));
		inx1 = front_link_spaces[index1];//��ʵ�ӿռ�����
		inx2 = front_link_spaces[index2];
		//���㷽��Ͳ���������2������
		size_t layer_num = 2;
		std::vector<std::vector<std::vector<size_t>>> all_direction_layer_nei;
		std::vector<std::vector<size_t>> layer_nei1;//����������
		std::vector<std::vector<size_t>> layer_nei2;
		std::vector<size_t> first_nei2;
		std::vector<size_t> second_nei1;
		std::vector<size_t> second_nei2;
		for (size_t j = 0; j < layer_num; ++j) {
			std::vector<size_t> temp;
			for (size_t i = 0; i < manifold_dist[index1].size(); ++i) {
				if (manifold_dist[index1][i] == j+1) {
					temp.push_back(i);
				}
				layer_nei1.emplace_back(temp);
			}
			temp.clear();
			for (size_t i = 0; i < manifold_dist[index2].size(); ++i) {
				if (manifold_dist[index2][i] == j + 1) {
					temp.push_back(i);
				}
				layer_nei2.emplace_back(temp);
			}
			all_direction_layer_nei.emplace_back(layer_nei1);
			all_direction_layer_nei.emplace_back(layer_nei2);
		}
		//ĩ���ӿռ������
		std::vector<std::vector<Real>> space_center;
		//ĩ���ӿռ�����
		std::vector<size_t> space_inx = { inx1,inx2 };
		size_t dim = CAST_CONOP(pro)->numberVariables();
		for (size_t i = 0; i < space_inx.size(); ++i) {
			std::vector<Real> center1;
			auto& front_sols = getMO_HLC().getSubspaceInfo(space_inx[i]).m_front_sol_in_subspace;
			for (size_t k = 0; k < dim; ++k) {
				Real sum = 0.;
				for (size_t j = 0; j < front_sols.size(); ++j) {
					sum += (front_sols[j]->variable().vect()[k]);
				}
				center1.push_back(sum / (Real)front_sols.size());
			}
			space_center.emplace_back(center1);
		}
		//���򰴱������
		std::vector<std::vector<std::vector<Real>>> all_delta;
		for (size_t p = 0; p < all_direction_layer_nei.size(); ++p) {
			std::vector<std::vector<Real>> delta1;//��ͬ��ķ���
			for (size_t i = 0; i < all_direction_layer_nei[p].size(); ++i) {
				Real ratio = 1. / std::pow((1. + i), 2);
				if (all_direction_layer_nei[p][i].size() > 0) {
					for (size_t j = 0; j < all_direction_layer_nei[p][i].size(); ++j) {
						//�����ӿռ��ڵ�ǰ�ؽ������
						auto& front_sols2 = getMO_HLC().getSubspaceInfo(front_link_spaces[all_direction_layer_nei[p][i][j]]).m_front_sol_in_subspace;
						std::vector<Real> center0;
						for (size_t q = 0; q < dim; ++q) {
							Real sum = 0.;
							for (size_t y = 0; y < front_sols2.size(); ++y) {
								sum += (front_sols2[y]->variable().vect()[q]);
							}
							center0.push_back(sum / (Real)front_sols2.size());
						}
						std::vector<Real> temp_direction;
						for (size_t q = 0; q < dim; ++q) {
							temp_direction.push_back(ratio * (space_center[p][q] - center0[q]));
						}
						delta1.emplace_back(temp_direction);
					}
				}
			}
			all_delta.emplace_back(delta1);
		}
		
		//����
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t p = 0; p < all_delta.size(); ++p) {
			auto& new_sol1 = space_center[p];
			auto& end_v1 = space_center[p];
			Real rate = rnd->uniform.next();
			for (size_t i = 0; i < all_delta[p].size(); ++i) {
				for (size_t j = 0; j < dim; ++j) {
					new_sol1[j] += (rate * all_delta[p][i][j]);
					end_v1[j] += all_delta[p][i][j];
					if (new_sol1[j] < bound[j].first) {
						new_sol1[j] = 2 * bound[j].first - new_sol1[j];
					}
					else if (new_sol1[j] > bound[j].second) {
						new_sol1[j] = 2 * bound[j].second - new_sol1[j];
					}
				}
			}
			//�ֲ�����
			auto mutate_sol = new_sol1;
			for (size_t j = 0; j < mutate_sol.size(); ++j) {
				//ÿһά����ʼ��ľ�����Сֵ�����Ŷ���Χ
				Real dist1 = std::fabs(new_sol1[j] - space_center[p][j]);
				Real dist2 = std::fabs(new_sol1[j] - end_v1[j]);
				auto amplitude = dist1 < dist2 ? dist1 : dist2;
				auto rand = rnd->uniform.next();
				mutate_sol[j] += ((2 * rand - 1) * amplitude);
			}
			//��Խ�紦��
			for (size_t j = 0; j < mutate_sol.size(); ++j) {
				if (mutate_sol[j] < bound[j].first) {
					mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
				}
				else if (mutate_sol[j] > bound[j].second) {
					mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
				}
			}

		}

	}

	//������ͨ�ӿռ���չ������
	std::vector<std::vector<Real>> SPMOEA::sampleByManifoldSpaces(std::vector<size_t>& front_link_spaces, std::vector<std::vector<size_t>>& manifold_dist, size_t sample_num, Problem *pro, Algorithm* alg,Random *rnd) {
		std::vector<std::vector<Real>> all_off;
		/*�м����*/
		//��ͨ�ӿռ��ڲ�����1�����ѡ��һ���ӿռ䣬���������ӿռ佻����2�����ݽ���Ŀ��ռ�ķֲ�����
		//�Ȼ�ȡ��ͨ�ӿռ��ڵ�ǰ�ؽ���Ŀ��ռ�ķֲ���Ȼ�����Ŀ��ռ�ǰ�ؾ���֮���ϡ��Ȳ���
		auto search_bound = CAST_CONOP(pro)->boundary();
		int method = 2;
		if (front_link_spaces.size() == 1) {//��ͨ��ֻ��һ��ǰ���ӿռ�
			//�����ӿռ����
			size_t index1 = front_link_spaces[0];
			//�ҳ������ӿռ�
			all_off = sampleByNeighSpace(index1, sample_num, pro, alg,rnd);

			//auto neighs = getMO_HLC().getSubspaceInfo(index1).m_sub_neighbors;
			////��ǰ�ӿռ�ǰ�ؽ�����������ӿռ�ǰ�ؽ⽻��
			//size_t pop_num = getMO_HLC().getSubspaceInfo(index1).m_subspace_front_sol.size();
			//for (auto ii : neighs) {
			//	pop_num += getMO_HLC().getSubspaceInfo(ii).m_subspace_front_sol.size();
			//}
			//SPMOEA_pop temp_pop(pop_num, pro);
			//size_t count = 0;
			//for (auto ii : neighs) {
			//	auto& front_sols = getMO_HLC().getSubspaceInfo(ii).m_subspace_front_sol;
			//	for (size_t j = 0; j < front_sols.size(); ++j) {
			//		temp_pop[count] = *front_sols[j];
			//		count++;
			//	}
			//}
			//auto& front_sols = getMO_HLC().getSubspaceInfo(index1).m_subspace_front_sol;
			//for (size_t j = 0; j < front_sols.size(); ++j) {
			//	temp_pop[count] = *front_sols[j];
			//	count++;
			//}
			//all_off = sampleByGA(temp_pop, sample_num, pro, rnd);
		}
		else {
			if (method == 1) {
				//���ѡ��ǰ���ӿռ�
				for (size_t j = 0; j < sample_num; ++j) {
					size_t sele_space_inx = (size_t)std::floor(front_link_spaces.size() * rnd->uniform.next());
					auto off = sampleByFrontNeighSpace(sele_space_inx, front_link_spaces, manifold_dist, 1,pro,alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
			}
			else if (method == 2) {
				//����ǰ���ӿռ�Ĳ����ܶ�ѡ�񽻻��ӿռ�
				std::vector<size_t> front_space_fre;
				for (size_t k = 0; k < front_link_spaces.size(); ++k) {
					front_space_fre.push_back(getMO_HLC().getSubspaceInfo(front_link_spaces[k]).m_sub_freq);
				}
				size_t select_num = 0;
				std::vector<size_t> sele_inx;
				while (select_num < sample_num) {
					size_t sele_space_inx = std::distance(front_space_fre.begin(), std::min_element(front_space_fre.begin(), front_space_fre.end()));
					sele_inx.push_back(sele_space_inx);
					front_space_fre[sele_space_inx]++;
					select_num++;
				}
				for (size_t j = 0; j < sele_inx.size(); ++j) {
					auto off = sampleByFrontNeighSpace(sele_inx[j], front_link_spaces, manifold_dist, 1, pro,alg, rnd);
					for (size_t k = 0; k < off.size(); ++k) {
						all_off.emplace_back(off[k]);
					}
				}
			}
			else {
				//����Ŀ��ռ�ֲ�ѡ��ǰ���ӿռ�
				Population<Solution<>> temp_pop;
				for (size_t j = 0; j < front_link_spaces.size(); ++j) {
					auto& front_sols = getMO_HLC().getSubspaceInfo(front_link_spaces[j]).m_subspace_front_sol;
					for (size_t k = 0; k < front_sols.size(); ++k) {
						temp_pop.append(*front_sols[k]);
					}
				}
				SPMOEA::NDSort(temp_pop);
				Population<Solution<>> front_pop;
				std::vector<std::vector<Real>> obj_sols;
				for (size_t j = 0; j < temp_pop.size(); ++j) {
					if (temp_pop[j].fitness() == 0) {
						front_pop.append(temp_pop[j]);
						obj_sols.emplace_back(temp_pop[j].objective());
					}
				}
				//��ǰ�ص���о���

				//�ҳ����������

				//��������
			}
			
		}
		
		return all_off;
	}

	std::vector<std::vector<Real>> SPMOEA::sampleByNeighGradient(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			if (his_sols.size() < 3) {
				auto off = sampleInRange(sol, 1, pro, alg, rnd);
				for (auto ii : off) {
					out_off.emplace_back(ii);
					if (ifSame(sol.variable().vect(), ii)) {
						size_t a = 1;
					}
				}
			}
			else {//���ӿռ��ǰ�ؽ��ƽ�
				auto base_sol = sol.variable().vect();
				Solution<> ind2(sol);
				auto& front_sols = getMO_HLC().getSubspaceInfo(space).m_subspace_front_sol;
				//Ѱ���ӿռ��ڵ��ݶȷ���
				
				//�ȿ��ý��ǲ����ӿռ�ǰ�ؽ�
				bool front_flag = false;
				for (size_t j = 0; j < front_sols.size(); ++j) {
					auto temp_sol = front_sols[j]->variable().vect();
					if (ifSame(base_sol, temp_sol)) {
						front_flag = true;
						break;
					}
				}
				Solution<> ind1(ind2);
				std::vector<Real> delta_v;
				if (front_flag) {
					//Ѱ��֧�����������Ҳ��������ҷ�֧������
					size_t inx2 = 0;
					std::vector<Real> sol2;
					//std::vector<size_t> dominant_inx;
					//std::vector<size_t> non_dominant_inx;
					//for (size_t j = 0; j < his_sols.size(); ++j) {
					//	auto dominance_ship = objectiveCompare(sol.objective(), his_sols[j]->objective(), CAST_CONOP(pro)->optimizeMode());
					//	if (dominance_ship == Dominance::kDominant) {
					//		dominant_inx.push_back(j);
					//	}
					//}
					//if (dominant_inx.size() > 0) {
					//	inx2 = dominant_inx[(size_t)std::floor(dominant_inx.size() * rnd->uniform.next())];
					//}
					//else {
					//	//��ʷ�������ѡ
					//	do {
					//		inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					//	} while (ifSame(base_sol,his_sols[inx2]->variable().vect()));
					//}
					//��ʷ�������ѡ
					do {
						inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					} while (ifSame(base_sol, his_sols[inx2]->variable().vect()));
					sol2 = his_sols[inx2]->variable().vect();
					for (size_t j = 0; j < sol2.size(); ++j) {
						delta_v.push_back(base_sol[j] - sol2[j]);
					}
					ind1.variable() = his_sols[inx2]->variable();
					ind1.objective() = his_sols[inx2]->objective();
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					std::vector<size_t> dominant_inx;
					for (size_t j = 0; j < front_sols.size(); ++j) {
						auto objs = front_sols[j]->objective();
						auto dominance_ship = objectiveCompare(sol.objective(), objs, CAST_CONOP(pro)->optimizeMode());
						if (dominance_ship == Dominance::kDominated) {
							dominant_inx.push_back(j);
						}
					}
					size_t inx2 = dominant_inx[(size_t)std::floor(dominant_inx.size() * rnd->uniform.next())];
					auto& sol2 = front_sols[inx2]->variable().vect();
					for (size_t j = 0; j < base_sol.size(); ++j) {
						delta_v.push_back(sol2[j] - base_sol[j]);
					}
					ind1.variable() = front_sols[inx2]->variable();
					ind1.objective() = front_sols[inx2]->objective();
				}

				std::vector<Real> new_sol(delta_v.size());
				Real ratio = rnd->uniform.next();
				if (!front_flag) {
					for (size_t j = 0; j < delta_v.size(); ++j) {
						delta_v[j] *= 2;
					}
				}
				std::vector<size_t> dim_inx;
				for (size_t j = 0; j < delta_v.size(); ++j) {
					dim_inx.push_back(j);
				}
				while (dim_inx.size() > 0) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < dim_inx.size(); ++j) {
						new_sol[dim_inx[j]] = base_sol[dim_inx[j]] + ratio * delta_v[dim_inx[j]];
						if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
							temp.push_back(dim_inx[j]);
						}
					}
					dim_inx = temp;
					if (dim_inx.size() > 0) {
						ratio /= 1.2;
					}
				}
				//Real raodong=ratio > 0.5 ? 1 - ratio : ratio;//���Ŷ����Ŷ�������ratio���
				//auto rand = rnd->uniform.next();
				//for (size_t j = 0; j < delta_v.size(); ++j) {
				//	dim_inx.push_back(j);
				//}
				//while (dim_inx.size() > 0) {
				//	std::vector<size_t> temp;
				//	for (size_t j = 0; j < dim_inx.size(); ++j) {
				//		new_sol[dim_inx[j]] += ((2 * rand - 1) * (raodong / 2) * delta_v[dim_inx[j]]);
				//		if (new_sol[dim_inx[j]] < bound[dim_inx[j]].first || new_sol[dim_inx[j]] > bound[dim_inx[j]].second) {
				//			temp.push_back(dim_inx[j]);
				//		}
				//	}
				//	dim_inx = temp;
				//	if (dim_inx.size() > 0) {
				//		rand /= 1.2;
				//	}
				//}
				//��new_solΪԲ�ģ���base_sol�ľ���Ϊ�뾶���Ʊ��췶Χ
				Real radius = 0.;
				for (size_t j = 0; j < new_sol.size(); ++j) {
					radius += std::pow((new_sol[j] - base_sol[j]), 2);
				}
				radius = std::sqrt(radius);
				auto mutate_sol = new_sol;
				Real temp_dist = 0.;
				do {
					for (size_t j = 0; j < new_sol.size(); ++j) {
						Real randn = rnd->uniform.next();
						mutate_sol[j] = new_sol[j] + (2 * randn - 1) * radius / 2;
					}
					temp_dist = euclideanDistance(new_sol.begin(), new_sol.end(), mutate_sol.begin());
				} while (temp_dist > radius);
				//��Խ�紦��
				for (size_t j = 0; j < mutate_sol.size(); ++j) {
					if (mutate_sol[j] < bound[j].first) {
						mutate_sol[j] = 2 * bound[j].first - mutate_sol[j];
					}
					else if (mutate_sol[j] > bound[j].second) {
						mutate_sol[j] = 2 * bound[j].second - mutate_sol[j];
					}
				}
				Solution<> ind3(ind2);
				ind3.variable().vect() = new_sol;
				ind3.evaluate(pro, alg);
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
				m_interactive_sol_pair.emplace_back(temp_pair);

				out_off.emplace_back(mutate_sol);
				if (ifSame(sol.variable().vect(), mutate_sol)) {
					size_t a = 1;
				}
			}
		}
		return out_off;
	}

	void SPMOEA::updateInteractiveSols() {
		//��������ݻ��켣
		std::vector<std::pair<std::shared_ptr<Solution<>>, std::shared_ptr<Solution<>>>> temp_ind;
		auto& front_clusters = getFrontRegionLinkSpace();
		auto front_spaces = getFrontSpace();
		for (size_t i = 0; i < getPop().size(); ++i) {
			//for (size_t k = 0; k < getPop()[i].size(); ++k) {
			//	temp_ind.emplace_back(std::make_pair(std::make_shared<Solution<>>(getPop()[i][k]), std::make_shared<Solution<>>(getPop()[i].getOffspring()[k])));
			//}
			////�����Ӵ�����
			//for (size_t k = 0; k < getPop()[i].size(); ++k) {
			//	getPop()[i].getOffspring()[getPop()[i].size() + k] = getPop()[i][k];
			//}
			//SPMOEA::NDSort(getPop()[i].getOffspring());
			//for (size_t k = 0; k < getPop()[i].size(); ++k) {
			//	getPop()[i][k].setFitness(getPop()[i].getOffspring()[getPop()[i].size() + k].fitness());
			//}
			//for (size_t k = 0; k < temp_ind.size(); ++k) {
			//	temp_ind[k].first->setFitness(getPop()[i][k].fitness());
			//	temp_ind[k].second->setFitness(getPop()[i].getOffspring()[k].fitness());
			//}
			auto& interactive_sols = getInteractiveSols();
			for (size_t k = 0; k < interactive_sols.size(); ++k) {
				temp_ind.emplace_back(std::make_pair<>(interactive_sols[k].front(), interactive_sols[k].back()));
			}
			//�켣������
			for (size_t k = 0; k < getPop()[i].size(); ++k) {
				getPop()[i].getOffspring()[getPop()[i].size() + k] = *(interactive_sols[k].front());
			}
			SPMOEA::NDSort(getPop()[i].getOffspring());
			/*for (size_t k = 0; k < getPop()[i].size(); ++k) {
				getPop()[i][k].setFitness(getPop()[i].getOffspring()[getPop()[i].size() + k].fitness());
			}*/
			for (size_t k = 0; k < temp_ind.size(); ++k) {
				temp_ind[k].first->setFitness(getPop()[i].getOffspring()[k+temp_ind.size()].fitness());
				temp_ind[k].second->setFitness(getPop()[i].getOffspring()[k].fitness());
			}
			getGenEvolveLocus().emplace_back(temp_ind);
			//�����ռ��ķֲ�
			std::vector<size_t> parents_in_space;
			std::vector<size_t> offs_in_space;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto sol = getPop()[i][j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				parents_in_space.push_back(space_inx);
			}
			getParentSpaces() = parents_in_space;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto sol = getPop()[i].getOffspring()[j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				offs_in_space.push_back(space_inx);
			}
			getOffspringSpaces() = offs_in_space;
			//�������ڵ��ӿռ�
			std::vector<size_t> ind_space;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto sol = getPop()[i][j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				if (std::find(ind_space.begin(), ind_space.end(), space_inx) == ind_space.end()) {
					ind_space.push_back(space_inx);
				}
			}
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto sol = getPop()[i].getOffspring()[j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				if (std::find(ind_space.begin(), ind_space.end(), space_inx) == ind_space.end()) {
					ind_space.push_back(space_inx);
				}
			}
			//
			std::vector<std::vector<size_t>> plot_spaces;
			for (size_t j = 0; j < front_clusters.size(); ++j) {
				std::vector<size_t> temp;
				for (size_t k = 0; k < ind_space.size(); ++k) {
					if (std::find(front_clusters[j].begin(), front_clusters[j].end(), ind_space[k]) != front_clusters[j].end()) {
						temp.push_back(ind_space[k]);
					}
				}
				plot_spaces.emplace_back(temp);
				//plot_spaces.emplace_back(front_clusters[j]);
			}
			for (size_t k = 0; k < ind_space.size(); ++k) {
				if (std::find(front_spaces.begin(), front_spaces.end(), ind_space[k]) == front_spaces.end()) {
					std::vector<size_t> temp;
					temp.push_back(ind_space[k]);
					plot_spaces.emplace_back(temp);
				}
			}
			getIndSpaces() = plot_spaces;
		}
	}

	void SPMOEA::setPopInitRange(SPMOEA_pop& pop, const std::vector<std::pair<Real, Real>>& pre_bound, const std::vector<std::pair<Real, Real>>& new_bound,Problem *pro) {
		//������Ⱥ��ʼ������λ��
		for (size_t j = 0; j < pro->numberVariables(); ++j) {
			for (size_t k = 0; k < pop.size(); ++k) {
				auto temp = pop[k].variable().vect()[j];
				pop[k].variable().vect()[j] = (temp - pre_bound[j].first) / (pre_bound[j].second - pre_bound[j].first) * (new_bound[j].second - new_bound[j].first) + new_bound[j].first;
			}
		}
	}

	void SPMOEA::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//������Ⱥ����ʷǰ�ؽ��ڱ߽��ڵı������������Դ
		size_t max_pop_size = getPopsize();
		//updateFrontRegionLinkSpace();
		if (m_stage_last_time > 0) {
			auto front_link_spaces = getFrontRegionLinkSpace();
			size_t M = CAST_CONOP(pro)->numberObjectives();
			size_t total_size = 0;
			for (size_t i = 0; i < front_link_spaces.size(); ++i) {
				total_size += (1 * front_link_spaces[i].size());
			}
			assign_pop_resource.push_back(max_pop_size);
		}
		else {
			assign_pop_resource.push_back(max_pop_size);
		}
		updatePopResource(assign_pop_resource);
	}

	void SPMOEA::updateRegionLinkSpace(size_t inx) {
		//getRegionLinkSpace().clear();
		auto cluster = getMO_HLC().getCluster(inx);
		//���ҵ�ǰ���ӿռ�
		std::vector<size_t> front_space;
		for (size_t i = 0; i < cluster[0].size(); ++i) {
			if (getMO_HLC().getSubspaceInfo(cluster[0][i]).m_best_rank == 0) {
				front_space.push_back(cluster[0][i]);
			}
		}
		//�ٽ������ӿռ����
		if (front_space.empty()) {
			//��������ǰ����ͨ����
			for (size_t i = 0; i < cluster[0].size(); ++i) {
				if (getMO_HLC().getSubspaceInfo(cluster[0][i]).m_sub_best_rank == 0) {
					front_space.push_back(cluster[0][i]);
				}
			}
		}
		auto link_space = clusterRegionFrontSpace(front_space);
		getRegionLinkSpace()[inx] = link_space;
	}

	void SPMOEA::updateFrontRegionLinkSpace(Problem* pro, Random* rnd) {
		//���ҵ�ǰ���ӿռ�
		std::vector<size_t> front_space = getFrontSpace();
		auto link_space = clusterRegionFrontSpace(front_space);
		//auto link_space = linearClusterFrontSpace(front_space,pro);//�����ӿռ����
		getFrontRegionLinkSpace() = link_space;
	}

	void SPMOEA::updateFrontSpace() {
		//ʹ����ʷ��֧������ǰ���ӿռ�
		m_front_space.clear();
		std::vector<size_t> front_space;
		auto& his_front = getHisFrontSols();
		for (size_t i = 0; i < his_front.size(); ++i) {
			auto temp_sol = his_front[i]->variable().vect();
			auto inx = getMO_HLC().subspaceTree().getRegionIdx(temp_sol);
			if (std::find(front_space.begin(), front_space.end(), inx) == front_space.end()) {
				front_space.push_back(inx);
			}
		}
		m_front_space = front_space;
		/*for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).m_best_rank = INT16_MAX;
		}*/

		//����ǰ�����ӿռ�Ľ�
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).m_front_sol_in_subspace.clear();
		}
		/*std::vector<std::vector<std::pair<Real, Real>>> front_boxes;
		for (size_t i = 0; i < m_front_space.size(); ++i) {
			front_boxes.emplace_back(getMO_HLC().subspaceTree().getBox(m_front_space[i]));
		}*/
		for (size_t i = 0; i < his_front.size(); ++i) {
			auto var = his_front[i]->variable().vect();
			auto idx = getMO_HLC().subspaceTree().getRegionIdx(var);
			getMO_HLC().getSubspaceInfo(idx).m_front_sol_in_subspace.emplace_back(his_front[i]);
			getMO_HLC().getSubspaceInfo(idx).m_best_rank = 0;
		}
	}

	std::vector<std::tuple<size_t, size_t, Real>> SPMOEA::findCloseSpaces(std::vector<std::vector<size_t>>& clusters, Random *rnd) {
		std::vector<std::tuple<size_t, size_t, Real>> candidate;
		std::vector<std::vector<std::tuple<size_t, size_t, Real>>> all_pairs;
		for (size_t i = 0; i < clusters.size(); ++i) {
			//��ͨ���ӿռ�����λ��
			std::vector<std::vector<Real>> center_pos1;
			for (size_t j = 0; j < clusters[i].size(); ++j) {
				auto box = getMO_HLC().subspaceTree().getBox(clusters[i][j]);
				std::vector<Real> center;
				for (size_t k = 0; k < box.size(); ++k) {
					center.push_back((box[k].first + box[k].second) / 2);
				}
				center_pos1.emplace_back(center);
			}
			std::vector<std::tuple<size_t, size_t, Real>> temp_candidate;
			for (size_t j = 0; j < clusters.size(); ++j) {
				//������ͨ�ӿռ�������ӿռ��
				if (j < i) {
					temp_candidate.emplace_back(std::make_tuple(0, 0, 0.));
				}
				else if (j == i) {
					temp_candidate.emplace_back(std::make_tuple(0, 0, (Real)INT16_MAX));
				}
				else {
					std::vector<std::vector<Real>> center_pos2;
					for (size_t p = 0; p < clusters[j].size(); ++p) {
						auto box = getMO_HLC().subspaceTree().getBox(clusters[j][p]);
						std::vector<Real> center;
						for (size_t k = 0; k < box.size(); ++k) {
							center.push_back((box[k].first + box[k].second) / 2);
						}
						center_pos2.emplace_back(center);
					}
					std::vector<std::vector<Real>> dist;
					size_t inx1;
					size_t inx2;
					std::vector<Real> min_dist;
					for (size_t k = 0; k < center_pos1.size(); ++k) {
						std::vector<Real> temp_dist;
						for (size_t p = 0; p < center_pos2.size(); ++p) {
							temp_dist.push_back(euclideanDistance(center_pos1[k].begin(), center_pos1[k].end(), center_pos2[p].begin()));
						}
						min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
						dist.emplace_back(temp_dist);
					}
					inx1 = std::distance(min_dist.begin(), std::min_element(min_dist.begin(), min_dist.end()));
					Real min_vv = *std::min_element(dist[inx1].begin(), dist[inx1].end());
					std::vector<size_t> same_inx;
					for (size_t k = 0; k < dist[inx1].size(); ++k) {
						if (dist[inx1][k] == min_vv) {
							same_inx.push_back(k);
						}
					}
					auto select_inx = (size_t)std::floor(same_inx.size() * rnd->uniform.next());
					temp_candidate.emplace_back(std::make_tuple(clusters[i][inx1], clusters[j][same_inx[select_inx]], dist[inx1][same_inx[select_inx]]));
				}
			}
			all_pairs.emplace_back(temp_candidate);
		}
		for (size_t i = 0; i < all_pairs.size(); ++i) {
			for (size_t j = 0; j < i; ++j) {
				all_pairs[i][j] = all_pairs[j][i];
			}
		}
		for (size_t i = 0; i < all_pairs.size(); ++i) {
			size_t min_inx = 0;
			Real min_d = (Real)INT16_MAX;
			for (size_t j = 0; j < all_pairs[i].size(); ++j) {
				if (std::get<2>(all_pairs[i][j]) < min_d) {
					min_d = std::get<2>(all_pairs[i][j]);
					min_inx = j;
				}
			}
			candidate.emplace_back(all_pairs[i][min_inx]);
		}
		return candidate;
	}

	std::vector<std::vector<size_t>> SPMOEA::calSpaceManifoldDist(std::vector<size_t>& spaces) {
		std::vector<std::vector<size_t>> manifold_dist;
		//for (size_t i = 0; i < spaces.size(); ++i) {
		//	std::vector<size_t> temp_dist;
		//	for (size_t j = 0; j < spaces.size(); ++j) {
		//		if (j < i) {
		//			temp_dist.push_back(INT16_MAX);
		//		}
		//		else if (j == i) {
		//			temp_dist.push_back(0);
		//		}
		//		else {
		//			auto nei = getMO_HLC().getSubspaceInfo(spaces[i]).m_sub_neighbors;
		//			std::vector<size_t> neigh;
		//			std::vector<size_t> has_index;
		//			has_index.push_back(spaces[i]);
		//			for (auto ii : nei) {
		//				neigh.push_back(ii);
		//				has_index.push_back(ii);
		//			}
		//			size_t count_layer = 1;
		//			while (std::find(neigh.begin(), neigh.end(), spaces[j]) == neigh.end()) {
		//				count_layer++;
		//				//����neigh
		//				std::vector<size_t> temp_neigh;
		//				for (size_t k = 0; k < neigh.size(); ++k) {
		//					auto temp_nei = getMO_HLC().getSubspaceInfo(neigh[k]).m_sub_neighbors;
		//					for (auto jj : temp_nei) {
		//						if (std::find(spaces.begin(), spaces.end(), jj) != spaces.end()) {
		//							if (std::find(has_index.begin(), has_index.end(), jj) == has_index.end()) {
		//								has_index.push_back(jj);
		//								temp_neigh.push_back(jj);
		//							}
		//						}
		//					}
		//				}
		//				neigh = temp_neigh;
		//			}
		//			temp_dist.push_back(count_layer);
		//		}
		//	}
		//	manifold_dist.emplace_back(temp_dist);
		//}

		//for (size_t i = 0; i < spaces.size(); ++i) {
		//	for (size_t j = 0; j < i; ++j) {
		//		manifold_dist[i][j] = manifold_dist[j][i];
		//	}
		//}

		////�ȵõ�ÿһ���ӿռ���ڽӱ�
		//size_t total_count = spaces.size();
		//for (size_t i = 0; i < spaces.size(); ++i) {
		//	std::vector<size_t> temp_dist(spaces.size());
		//	temp_dist[i] = 0;
		//	auto nei = getMO_HLC().getSubspaceInfo(spaces[i]).m_sub_neighbors;
		//	for (auto ii : nei) {
		//		if (std::find(spaces.begin(), spaces.end(), ii) != spaces.end()) {
		//			size_t inx = std::distance(spaces.begin(), std::find(spaces.begin(), spaces.end(), ii));
		//			temp_dist[i] = 1;
		//			total_count++;
		//		}
		//	}
		//	manifold_dist.emplace_back(temp_dist);
		//}
		////�ٸ����ڽӱ�������ξ������
		//while (1) {
		//	for (size_t i = 0; i < manifold_dist.size(); ++i) {
		//		for (size_t j = 0; j < manifold_dist[i].size(); ++j) {
		//			if (manifold_dist[i][j] == 0 && i!=j) {
		//				//�ҳ���j����1���У�
		//				std::vector<size_t> inx;
		//				//size_t count = 1;
		//				for (size_t k = i + 1; k < manifold_dist.size(); ++k) {
		//					if (manifold_dist[k][j] == 1) {
		//						inx.push_back(k);
		//					}
		//				}
		//				while (!inx.empty()) {
		//					std::vector<size_t> temp_inx;
		//					for (size_t k = 0; k < inx.size(); ++k) {
		//						if (manifold_dist[i][k] > 0) {
		//							 auto temp= manifold_dist[i][inx[k]] + manifold_dist[inx[k]][j];
		//							temp_inx.push_back(temp);
		//						}
		//					}
		//					if (temp_inx.size() > 0) {
		//						manifold_dist[i][j] = *std::min_element(temp_inx.begin(), temp_inx.end());
		//						total_count++;
		//						break;
		//					}
		//				}
		//			}
		//		}
		//	}
		//	//����Ƿ����
		//	if (total_count == spaces.size() * spaces.size()) {
		//		break;
		//	}
		//}
		
		
		//�ȵõ�ÿһ���ӿռ���ڽӱ�
		for (size_t i = 0; i < spaces.size(); ++i) {
			std::vector<size_t> temp_dist(spaces.size(),10e9);
			temp_dist[i] = 0;
			auto nei = getMO_HLC().getSubspaceInfo(spaces[i]).m_sub_neighbors;
			for (auto ii : nei) {
				if (std::find(spaces.begin(), spaces.end(), ii) != spaces.end()) {
					size_t inx = std::distance(spaces.begin(), std::find(spaces.begin(), spaces.end(), ii));
					temp_dist[inx] = 1;
				}
			}
			manifold_dist.emplace_back(temp_dist);
		}
		//�ٸ����ڽӱ�������ξ������
		size_t dist_count = 2;
		while (1) {
			for (size_t i = 0; i < manifold_dist.size(); ++i) {
				for (size_t j = 0; j < manifold_dist[i].size(); ++j) {
					if (manifold_dist[i][j] >= 10e9 && i != j) {
						//����ȡ��i��Ԫ��
						std::vector<size_t> row_inx;
						for (size_t k = 0; k < manifold_dist[i].size(); ++k) {
							if (k != i && k != j) {
								row_inx.push_back(manifold_dist[i][k]);
							}
						}
						//����ȡ��j��Ԫ��
						std::vector<size_t> col_inx;
						for (size_t k = 0; k < manifold_dist[i].size(); ++k) {
							if (k != i && k != j) {
								col_inx.push_back(manifold_dist[k][j]);
							}
						}
						//�õ�����Ͳ�ȡ��Сֵ
						std::vector<size_t> sum;
						for (size_t k = 0; k < col_inx.size(); ++k) {
							sum.push_back(row_inx[k] + col_inx[k]);
						}
						auto min_dist= *std::min_element(sum.begin(), sum.end());
						if (min_dist == dist_count) {
							manifold_dist[i][j] = min_dist;
							manifold_dist[j][i] = min_dist;
						}
						////�õ�����Ͳ�ȡ��Сֵ
						//std::vector<size_t> sum;
						//for (size_t k = 0; k < col_inx.size(); ++k) {
						//	sum.push_back(row_inx[k] + col_inx[k]);
						//}
						//auto min_dist = *std::min_element(sum.begin(), sum.end());
						//manifold_dist[i][j] = min_dist;
						//manifold_dist[j][i] = min_dist;
					}
				}
			}
			dist_count++;
			//ȡ�������ֵ
			std::vector<size_t> max_v;
			for (size_t i = 0; i < manifold_dist.size(); ++i) {
				max_v.push_back(*std::max_element(manifold_dist[i].begin(), manifold_dist[i].end()));
			}
			auto vv = *std::max_element(max_v.begin(), max_v.end());
			if (vv < 10e9) {
				break;
			}
		}

		return manifold_dist;
	}

	bool SPMOEA::ifFrontChanged(std::vector<size_t>& pre_front_spaces) {
		if (pre_front_spaces.size() != m_front_space.size()) {
			setFrontLastGens(1);
			return true;
		}
		else {
			for (size_t i = 0; i < m_front_space.size(); ++i) {
				if (std::find(pre_front_spaces.begin(), pre_front_spaces.end(), m_front_space[i]) == pre_front_spaces.end()) {
					setFrontLastGens(1);
					return true;
				}
			}
			m_front_last_gens++;
			return false;
		}
	}


	/*����ѡ��*/
	//���Ǿ��߿ռ��Ŀ��ռ���ϡ���ѡ���
	void SPMOEA::sparseSelection(Problem *pro) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<size_t>> layer_inds_inx;//�ֲ㴢���������
			std::vector<size_t> flag(getPop()[i].getOffspring().size(), 0);
			int temp_rank = 0;
			while (std::find(flag.begin(), flag.end(), 0) != flag.end()) {
				std::vector<size_t> temp_inx;
				for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
					if (getPop()[i].getOffspring()[j].fitness() == temp_rank) {
						temp_inx.push_back(j);
						flag[j] = 1;
					}
				}
				temp_rank++;
				layer_inds_inx.emplace_back(temp_inx);
			}
			//��һ��ֵ
			int space = 1;//1�����߿ռ䣻2��Ŀ��ռ�
			std::vector<std::vector<Real>> points;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				if (space == 1) {
					points.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
				}
				else if (space == 2) {
					points.emplace_back(getPop()[i].getOffspring()[j].objective());
				}
			}
			dataNormalize(points);
			//������ѡ�񣬲��Ҵ�ѡ������ѡ���������ռ�ľ��벻��̫��
			std::vector<size_t> select_inx;
			std::vector<std::vector<size_t>> layer_inds_status;//�������ѡ��״̬
			for (size_t j = 0; j < layer_inds_inx.size(); ++j) {
				std::vector<size_t> temp(layer_inds_inx[j].size(),0);
				layer_inds_status.emplace_back(temp);
			}
			//���е�ռ��ڵľ���
			std::vector<std::vector<Real>> p_dist;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto p1 = points[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = points[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				p_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					if (k < j) {
						p_dist[j][k] = p_dist[k][j];
					}
					else if (k == j) {
						p_dist[j][k] = INT16_MAX;
					}
				}
			}
			//�ҳ���ǰǰ���ӿռ��е����ά��ֵ����Ϊ�����ռ��ӵ���뾶
			Real niche_radius = INT16_MAX;
			//size_t num_var = CAST_CONOP(pro)->numberVariables();
			/*for (size_t j = 0; j < getFrontSpace().size(); ++j) {
				auto box = getMO_HLC().subspaceTree().getBox(getFrontSpace()[j]);
				for (size_t k = 0; k < num_var; ++k) {
					Real dim_span = box[k].second - box[k].first;
					if (niche_radius > dim_span) {
						niche_radius = dim_span;
					}
				}
			}*/
			//������Ⱥ��С�����ƽ��ֵ
			Real mean_min_dist = 0.;
			for (size_t j = 0; j < p_dist.size(); ++j) {
				auto min_d = *std::min_element(p_dist[j].begin(), p_dist[j].end());
				mean_min_dist += min_d;
			}
			mean_min_dist /= p_dist.size();
			//ѭ��ѡ��
			niche_radius = mean_min_dist;
			niche_radius *= 1.1;
			while (select_inx.size() < getPop()[i].size()) {
				niche_radius /= 1.1;
				for (size_t j = 0; j < layer_inds_inx.size(); ++j) {
					for (size_t k = 0; k < layer_inds_inx[j].size(); ++k) {
						if (select_inx.empty()) {
							select_inx.push_back(layer_inds_inx[j][k]);
							layer_inds_status[j][k] = 1;
						}
						else if (layer_inds_status[j][k] == 0) {
							//��������ѡ��ľ���
							bool sele_flag = true;
							for (size_t p = 0; p < select_inx.size(); ++p) {
								if (p_dist[layer_inds_inx[j][k]][select_inx[p]] < niche_radius) {
									sele_flag = false;
									break;
								}
							}
							if (sele_flag) {
								select_inx.push_back(layer_inds_inx[j][k]);
								layer_inds_status[j][k] = 1;
							}
						}
						if (select_inx.size() >= getPop()[i].size()) {
							break;
						}
					}
					if (select_inx.size() >= getPop()[i].size()) {
						break;
					}
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//�������Ӵ��Ƚ�ѡ���õģ���֧��ķ�һ�����ϡ��ָ��ѡ��
	void SPMOEA::localSelection(Problem *pro) {
		//auto& interactive_sol = getInteractiveSols();
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			Population<Solution<>> residual_pop;
			std::vector<size_t> sele_inx;
			std::vector<size_t> residual_inx;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp,off,CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					sele_inx.push_back(getPop()[i].size()+j);
				}
				else if (dominance_ship == Dominance::kDominated) {
					sele_inx.push_back(j);
				}
				else if (dominance_ship==Dominance::kNonDominated) {
					residual_inx.push_back(j);
					residual_inx.push_back(getPop()[i].size() + j);
					residual_pop.append(getPop()[i][j]);
					residual_pop.append(getPop()[i].getOffspring()[j]);
				}
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			SPMOEA::NDSort(residual_pop);
			std::vector<std::vector<size_t>> residual_layer_inds_inx;//�ֲ㴢���������
			std::vector<size_t> flag(residual_pop.size(), 0);
			int temp_rank = 0;
			while (std::find(flag.begin(), flag.end(), 0) != flag.end()) {
				std::vector<size_t> temp_inx;
				for (size_t j = 0; j < residual_pop.size(); ++j) {
					if (residual_pop[j].fitness() == temp_rank) {
						temp_inx.push_back(j);
						flag[j] = 1;
					}
				}
				temp_rank++;
				residual_layer_inds_inx.emplace_back(temp_inx);
			}
			//��һ�����߱���ֵ
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������ѡ�񣬲��Ҵ�ѡ������ѡ���������ռ�ľ��벻��̫��
			std::vector<std::vector<size_t>> layer_inds_status;//�������ѡ��״̬
			for (size_t j = 0; j < residual_layer_inds_inx.size(); ++j) {
				std::vector<size_t> temp(residual_layer_inds_inx[j].size(), 0);
				layer_inds_status.emplace_back(temp);
			}
			//ʣ��������е��������ռ��ڵľ���
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < residual_inx.size(); ++j) {
				auto p1 = all_vars[residual_inx[j]];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					auto p2 = all_vars[k];
					auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
					if (dist == 0.) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			Real niche_radius = INT16_MAX;
			size_t num_var = CAST_CONOP(pro)->numberVariables();
			//������Ⱥ��С�����ƽ��ֵ
			Real mean_min_dist = 0.;
			for (size_t j = 0; j < residual_inx.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_min_dist += min_d;
			}
			mean_min_dist /= var_dist.size();
			//ѭ��ѡ������ѡ��Զ���ҿ�ǰ�ĸ���
			niche_radius = mean_min_dist;
			niche_radius *= 1.1;
			while (sele_inx.size() < getPop()[i].size()) {
				niche_radius /= 1.1;
				for (size_t j = 0; j < residual_layer_inds_inx.size(); ++j) {
					for (size_t k = 0; k < residual_layer_inds_inx[j].size(); ++k) {
						if (sele_inx.empty()) {
							sele_inx.push_back(residual_inx[residual_layer_inds_inx[j][k]]);
							layer_inds_status[j][k] = 1;
						}
						else if (layer_inds_status[j][k] == 0) {
							//��������ѡ��ľ���
							bool sele_flag = true;
							for (size_t p = 0; p < sele_inx.size(); ++p) {
								if (var_dist[residual_layer_inds_inx[j][k]][sele_inx[p]] < niche_radius) {
									sele_flag = false;
									break;
								}
							}
							if (sele_flag) {
								sele_inx.push_back(residual_inx[residual_layer_inds_inx[j][k]]);
								layer_inds_status[j][k] = 1;
							}
						}
						if (sele_inx.size() >= getPop()[i].size()) {
							break;
						}
					}
					if (sele_inx.size() >= getPop()[i].size()) {
						break;
					}
				}
			}
			//ʹ��ѡ�����½���¸����ݻ��켣

			for (size_t j = 0; j < sele_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
			}
		}
	}

	//�������Ӵ��ȣ���֧��ķ�һ�𣬱Ƚ�����ѡ���ھ��ߺ�Ŀ��ռ�ľ��룬�ٲ���֧����ѡ��
	void SPMOEA::localCrowdSelection(Problem *pro,Random *rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			Population<Solution<>> residual_pop;
			std::vector<size_t> select_inx;
			std::vector<size_t> residual_inx;
			//���ѡ�����
			std::vector<size_t> no_select_pair;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					/*auto temp_dist = var_dist[getPop()[i].size() + j];
					temp_dist.erase(temp_dist.begin()+j);
					auto min_dist = *std::min_element(temp_dist.begin(), temp_dist.end());
					if (min_dist > mean_var_min_dist) {
						select_inx.push_back(getPop()[i].size() + j);
					}
					else {
						select_inx.push_back(j);
					}*/
					select_inx.push_back(getPop()[i].size() + j);
				}
				else if (dominance_ship == Dominance::kDominated) {
					/*auto temp_dist = var_dist[j];
					temp_dist.erase(temp_dist.begin() + getPop()[i].size()+ j);
					auto min_dist = *std::min_element(temp_dist.begin(), temp_dist.end());
					if (min_dist > mean_var_min_dist) {
						select_inx.push_back(j);
					}
					else {
						select_inx.push_back(getPop()[i].size() + j);
					}*/
					select_inx.push_back(j);
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					no_select_pair.push_back(j);
					/*residual_inx.push_back(getPop()[i].size() + j);
					residual_pop.append(getPop()[i][j]);
					residual_pop.append(getPop()[i].getOffspring()[j]);*/
				}
				
			}
		    //ʣ�������ѡ�����С�����������ռ����������֧�䣬���ѡ�����rank�����
			while (no_select_pair.size() > 0) {
				size_t inx = no_select_pair.back();
				std::vector<Real> parent_two_min_dist;
				Real min_var_dist = INT16_MAX;
				Real min_obj_dist = INT16_MAX;
				if (select_inx.empty()) {
					min_var_dist = *std::min_element(var_dist[getPop()[i].size() + inx].begin(), var_dist[getPop()[i].size() + inx].end());
					min_obj_dist = *std::min_element(obj_dist[getPop()[i].size() + inx].begin(), obj_dist[getPop()[i].size() + inx].end());
				}
				else {
					for (size_t j = 0; j < select_inx.size(); ++j) {
						if (min_var_dist > var_dist[getPop()[i].size() + inx][select_inx[j]]) {
							min_var_dist = var_dist[getPop()[i].size() + inx][select_inx[j]];
						}
						if (min_obj_dist > obj_dist[getPop()[i].size() + inx][select_inx[j]]) {
							min_obj_dist = obj_dist[getPop()[i].size() + inx][select_inx[j]];
						}
					}
				}
				parent_two_min_dist.push_back(1. / min_var_dist);
				parent_two_min_dist.push_back(1. / min_obj_dist);
				std::vector<Real> off_two_min_dist;
				min_var_dist = INT16_MAX;
				min_obj_dist = INT16_MAX;
				if (select_inx.empty()) {
					min_var_dist = *std::min_element(var_dist[inx].begin(), var_dist[inx].end());
					min_obj_dist = *std::min_element(obj_dist[inx].begin(), obj_dist[inx].end());
				}
				else {
					for (size_t j = 0; j < select_inx.size(); ++j) {
						if (min_var_dist > var_dist[inx][select_inx[j]]) {
							min_var_dist = var_dist[inx][select_inx[j]];
						}
						if (min_obj_dist > obj_dist[inx][select_inx[j]]) {
							min_obj_dist = obj_dist[inx][select_inx[j]];
						}
					}
				}
				off_two_min_dist.push_back(1. / min_var_dist);
				off_two_min_dist.push_back(1. / min_obj_dist);
				Dominance dominance_ship = objectiveCompare(parent_two_min_dist, off_two_min_dist, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					select_inx.push_back(getPop()[i].size() + inx);
				}
				else if (dominance_ship == Dominance::kDominated) {
					select_inx.push_back(inx);
				}
				else {
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						select_inx.push_back(getPop()[i].size() + inx);
					}
					else {
						select_inx.push_back(inx);
					}
				}
				no_select_pair.pop_back();
			}

			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//���ڸ���ֲ�ѡ��
	void SPMOEA::localCrowdSelection2(Problem *pro, Random *rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			Population<Solution<>> residual_pop;
			std::vector<size_t> select_inx;
			std::vector<size_t> residual_inx;
			//���ѡ�����
			auto front_spaces = getFrontSpace();
			std::vector<size_t> no_select_pair;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(all_vars[getPop()[i].size() + j]);
					bool flag = false;
					//�ȿ������Ƿ�Ϊǰ���ӿռ�
					if (std::find(front_spaces.begin(), front_spaces.end(), space_inx) != front_spaces.end()) {
						flag = true;
					}
					if (flag) {
						//�ٿ��Ƿ�Ϊ��ʷǰ�ؽ�
						bool if_front = false;
						auto& sub_front_sols = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
						//�������Ƿ�Ϊ�ӿռ�ǰ�ظ���
						for (size_t k = 0; k < sub_front_sols.size(); ++k) {
							auto& sols = sub_front_sols[k]->variable().vect();
							if (ifSame(all_vars[getPop()[i].size() + j], sols)) {
								if_front = true;
								break;
							}
						}
						if (if_front) {
							select_inx.push_back(j);
						}
						else {
							select_inx.push_back(getPop()[i].size() + j);
						}
					}
					else {
						select_inx.push_back(getPop()[i].size() + j);
						//��������ѡ��ľ���
					}
				}
				else if (dominance_ship == Dominance::kDominated) {
					/*auto temp_dist = var_dist[j];
					temp_dist.erase(temp_dist.begin() + getPop()[i].size() + j);
					auto min_dist = *std::min_element(temp_dist.begin(), temp_dist.end());
					if (min_dist > mean_var_min_dist) {
						select_inx.push_back(j);
					}
					else {
						select_inx.push_back(getPop()[i].size() + j);
					}*/
					//��������ѡ�����
					if (select_inx.empty()) {
						select_inx.push_back(j);
					}
					else {
						bool flag = false;
						for (size_t k = 0; k < select_inx.size(); ++k) {
							if (obj_dist[j][select_inx[k]]<mean_obj_min_dist) {
								flag=true;
								break;
							}
						}
						if (flag) {
							select_inx.push_back(getPop()[i].size()+j);
						}
						else {
							select_inx.push_back(j);
						}
					}
					//select_inx.push_back(j);
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					//�ȱȽ�rank
					size_t p_rank = getPop()[i].getOffspring()[getPop()[i].size() + j].fitness();
					size_t o_rank = getPop()[i].getOffspring()[j].fitness();
					if (p_rank < o_rank) {
						select_inx.push_back(getPop()[i].size() + j);
					}
					else if (p_rank > o_rank) {
						select_inx.push_back(j);
					}
					else {
						//��������ѡ�����
						Real p_min_var_dist = INT16_MAX;
						Real p_min_obj_dist = INT16_MAX;
						Real o_min_var_dist = INT16_MAX;
						Real o_min_obj_dist = INT16_MAX;
						if (select_inx.empty()) {
							Real rand = rnd->uniform.next();
							if (rand > 0.5) {
								select_inx.push_back(getPop()[i].size() + j);
							}
							else {
								select_inx.push_back(j);
							}
						}
						else {
							for (size_t k = 0; k < select_inx.size(); ++k) {
								if (p_min_var_dist > var_dist[getPop()[i].size() + j][select_inx[k]]) {
									p_min_var_dist = var_dist[getPop()[i].size() + j][select_inx[k]];
								}
								if (p_min_obj_dist > obj_dist[getPop()[i].size() + j][select_inx[k]]) {
									p_min_obj_dist = obj_dist[getPop()[i].size() + j][select_inx[k]];
								}
							}
							for (size_t k = 0; k < select_inx.size(); ++k) {
								if (o_min_var_dist > var_dist[j][select_inx[k]]) {
									o_min_var_dist = var_dist[j][select_inx[k]];
								}
								if (o_min_obj_dist > obj_dist[j][select_inx[k]]) {
									o_min_obj_dist = obj_dist[j][select_inx[k]];
								}
							}
							std::vector<Real> parent_two_min_dist, off_two_min_dist;
							parent_two_min_dist.push_back(1. / p_min_var_dist);
							parent_two_min_dist.push_back(1. / p_min_obj_dist);
							off_two_min_dist.push_back(1. / o_min_var_dist);
							off_two_min_dist.push_back(1. / o_min_obj_dist);
							Dominance dominance_ship = objectiveCompare(parent_two_min_dist, off_two_min_dist, CAST_CONOP(pro)->optimizeMode());
							if (dominance_ship == Dominance::kDominant) {
								select_inx.push_back(getPop()[i].size() + j);
							}
							else if (dominance_ship == Dominance::kDominated) {
								select_inx.push_back(j);
							}
							else {
								Real rand = rnd->uniform.next();
								if (rand > 0.5) {
									select_inx.push_back(getPop()[i].size() + j);
								}
								else {
									select_inx.push_back(j);
								}
							}
						}
					}
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//����֧����ѡ�񣬷�֧�����ԱȽϺ�����ѡ����Ŀ��������ռ�ľ���
	void SPMOEA::localCrowdSelection3(Problem *pro, Random *rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			Population<Solution<>> residual_pop;
			std::vector<size_t> select_inx;
			std::vector<size_t> residual_inx;
			//���ѡ�����
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					select_inx.push_back(getPop()[i].size() + j);
				}
				else if (dominance_ship == Dominance::kDominated) {
					select_inx.push_back(j);
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					//�ȱȽ�rank
					size_t p_rank = getPop()[i].getOffspring()[getPop()[i].size() + j].fitness();
					size_t o_rank = getPop()[i].getOffspring()[j].fitness();
					/*if (p_rank < o_rank) {
						select_inx.push_back(getPop()[i].size() + j);
					}
					else if (p_rank > o_rank) {
						select_inx.push_back(j);
					}
					else {
						
					}*/
					//��������ѡ�����
					Real p_min_var_dist = INT16_MAX;
					Real p_min_obj_dist = INT16_MAX;
					Real o_min_var_dist = INT16_MAX;
					Real o_min_obj_dist = INT16_MAX;
					if (select_inx.empty()) {
						Real rand = rnd->uniform.next();
						if (rand > 0.5) {
							select_inx.push_back(getPop()[i].size() + j);
						}
						else {
							select_inx.push_back(j);
						}
					}
					else {
						for (size_t k = 0; k < select_inx.size(); ++k) {
							if (p_min_var_dist > var_dist[getPop()[i].size() + j][select_inx[k]]) {
								p_min_var_dist = var_dist[getPop()[i].size() + j][select_inx[k]];
							}
							if (p_min_obj_dist > obj_dist[getPop()[i].size() + j][select_inx[k]]) {
								p_min_obj_dist = obj_dist[getPop()[i].size() + j][select_inx[k]];
							}
						}
						for (size_t k = 0; k < select_inx.size(); ++k) {
							if (o_min_var_dist > var_dist[j][select_inx[k]]) {
								o_min_var_dist = var_dist[j][select_inx[k]];
							}
							if (o_min_obj_dist > obj_dist[j][select_inx[k]]) {
								o_min_obj_dist = obj_dist[j][select_inx[k]];
							}
						}
						std::vector<Real> parent_two_min_dist, off_two_min_dist;
						parent_two_min_dist.push_back(1. / p_min_var_dist);
						parent_two_min_dist.push_back(1. / p_min_obj_dist);
						off_two_min_dist.push_back(1. / o_min_var_dist);
						off_two_min_dist.push_back(1. / o_min_obj_dist);
						Dominance dominance_ship = objectiveCompare(parent_two_min_dist, off_two_min_dist, CAST_CONOP(pro)->optimizeMode());
						if (dominance_ship == Dominance::kDominant) {
							select_inx.push_back(getPop()[i].size() + j);
						}
						else if (dominance_ship == Dominance::kDominated) {
							select_inx.push_back(j);
						}
						else {
							Real rand = rnd->uniform.next();
							if (rand > 0.5) {
								select_inx.push_back(getPop()[i].size() + j);
							}
							else {
								select_inx.push_back(j);
							}
						}
					}
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//�ۺ��������أ��������Ӵ���֧���ϵ���Ƿ�λ���ӿռ�ǰ�أ��ӿռ�����ܶȣ���ͨ�ӿռ��ڸ���ı���
	void SPMOEA::localCrowdSelection4(Problem *pro, Random *rnd) {
		auto& front_clusters = getFrontRegionLinkSpace();
		std::vector<size_t> reserve_front(front_clusters.size(), 0);
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			Population<Solution<>> residual_pop;
			std::vector<size_t> select_inx;
			std::vector<size_t> residual_inx;
			auto front_spaces = getFrontSpace();
			//���ѡ�����
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					select_inx.push_back(getPop()[i].size() + j);
				}
				else if (dominance_ship == Dominance::kDominated) {
					select_inx.push_back(j);
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					auto p_space = getMO_HLC().subspaceTree().getRegionIdx(getPop()[i][j].variable().vect());
					auto& p_sol = getPop()[i][j].variable().vect();
					auto& o_sol = getPop()[i].getOffspring()[j].variable().vect();
					auto o_space = getMO_HLC().subspaceTree().getRegionIdx(getPop()[i].getOffspring()[j].variable().vect());
					auto& p_front_sols = getMO_HLC().getSubspaceInfo(p_space).m_subspace_front_sol;
					auto& o_front_sols = getMO_HLC().getSubspaceInfo(o_space).m_subspace_front_sol;
					if (std::find(front_spaces.begin(), front_spaces.end(), p_space) != front_spaces.end()) {
						//�����ĸ���ͨ��
						size_t p_cluster_inx;
						for (size_t k = 0; k < front_clusters.size(); ++k) {
							if (std::find(front_clusters[k].begin(), front_clusters[k].end(), p_space) != front_clusters[k].end()) {
								p_cluster_inx = k;
								break;
							}
						}
						if (std::find(front_spaces.begin(), front_spaces.end(), o_space) != front_spaces.end()) {
							size_t o_cluster_inx;
							for (size_t k = 0; k < front_clusters.size(); ++k) {
								if (std::find(front_clusters[k].begin(), front_clusters[k].end(), o_space) != front_clusters[k].end()) {
									o_cluster_inx = k;
									break;
								}
							}
							if ((p_cluster_inx == o_cluster_inx)|| (reserve_front[p_cluster_inx] == 1 && reserve_front[o_cluster_inx] == 1)) {
								//�����ӿռ��ǰ�ؽ�
								bool p_status = false;
								for (size_t k = 0; k < p_front_sols.size(); ++k) {
									if (ifSame(p_sol, p_front_sols[k]->variable().vect())) {
										p_status = true;
										break;
									}
								}
								bool o_status = false;
								for (size_t k = 0; k < o_front_sols.size(); ++k) {
									if (ifSame(o_sol, o_front_sols[k]->variable().vect())) {
										o_status = true;
										break;
									}
								}
								if (p_status && (!o_status)) {
									select_inx.push_back(getPop()[i].size() + j);
								}
								else if (o_status && (!p_status)) {
									select_inx.push_back(j);
								}
								else {
									/*Real randn = rnd->uniform.next();
									if (randn > 0.5) {
										select_inx.push_back(j);
									}
									else {
										select_inx.push_back(getPop()[i].size() + j);
									}*/
									//�Ƚϲ����ܶ�
									Real p_density = (Real)getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(p_space);
									Real o_density = (Real)getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(o_space);
									if (p_density < o_density) {
										select_inx.push_back(getPop()[i].size() + j);
									}
									else {
										select_inx.push_back(j);
									}
								}
								reserve_front[p_cluster_inx] = 1;
							}
							else {
								if (reserve_front[p_cluster_inx] == 0 && reserve_front[o_cluster_inx] == 1) {
									select_inx.push_back(getPop()[i].size() + j);
									reserve_front[p_cluster_inx] = 1;
								}
								else if (reserve_front[p_cluster_inx] == 1 && reserve_front[o_cluster_inx] == 0) {
									select_inx.push_back(j);
									reserve_front[o_cluster_inx] = 1;
								}
								else if (reserve_front[p_cluster_inx] == 0 && reserve_front[o_cluster_inx] == 0) {
									if (front_clusters[p_cluster_inx].size() < front_clusters[o_cluster_inx].size()) {
										select_inx.push_back(getPop()[i].size() + j);
										reserve_front[p_cluster_inx] = 1;
									}
									else {
										select_inx.push_back(j);
										reserve_front[o_cluster_inx] = 1;
									}
								}
							}
						}
						else {
							//�Ƚϲ����ܶ�
							Real p_density = (Real)getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size()/ getMO_HLC().subspaceTree().getBoxVolume(p_space);
							Real o_density = (Real)getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(o_space);
							if (p_density < o_density) {
								select_inx.push_back(getPop()[i].size() + j);
								reserve_front[p_cluster_inx] = 1;
							}
							else {
								select_inx.push_back(j);
							}
						}
					}
					else {
						if (std::find(front_spaces.begin(), front_spaces.end(), o_space) != front_spaces.end()) {
							size_t o_cluster_inx;
							for (size_t k = 0; k < front_clusters.size(); ++k) {
								if (std::find(front_clusters[k].begin(), front_clusters[k].end(), o_space) != front_clusters[k].end()) {
									o_cluster_inx = k;
									break;
								}
							}
							if (reserve_front[o_cluster_inx] == 1) {
								select_inx.push_back(getPop()[i].size() + j);
							}
							else {
								select_inx.push_back(j);
								reserve_front[o_cluster_inx] = 1;
							}
						}
						else {
							if (p_space == o_space) {
								select_inx.push_back(j);
							}
							else {
								//�Ƚϲ����ܶ�
								Real p_density = (Real)getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(p_space);
								Real o_density = (Real)getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(o_space);
								if (p_density < o_density) {
									select_inx.push_back(getPop()[i].size() + j);
								}
								else {
									select_inx.push_back(j);
								}
							}
							////�����ӿռ��ǰ�ؽ�
							//bool p_status = false;
							//for (size_t k = 0; k < p_front_sols.size(); ++k) {
							//	if (ifSame(p_sol, p_front_sols[k]->variable().vect())) {
							//		p_status = true;
							//		break;
							//	}
							//}
							//bool o_status = false;
							//for (size_t k = 0; k < o_front_sols.size(); ++k) {
							//	if (ifSame(o_sol, o_front_sols[k]->variable().vect())) {
							//		o_status = true;
							//		break;
							//	}
							//}
							//if (p_status && (!o_status)) {
							//	select_inx.push_back(getPop()[i].size() + j);
							//}
							//else if (o_status && (!p_status)) {
							//	select_inx.push_back(j);
							//}
							//else {
							//	//�Ƚϲ����ܶ�
							//	Real p_density = (Real)getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(p_space);
							//	Real o_density = (Real)getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / getMO_HLC().subspaceTree().getBoxVolume(o_space);
							//	if (p_density < o_density) {
							//		select_inx.push_back(getPop()[i].size() + j);
							//	}
							//	else {
							//		select_inx.push_back(j);
							//	}
							//}
						}
					}
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//��ͬ����ͨ�ӿռ��б������壬�ۺ��������أ��������Ӵ���֧���ϵ���Ƿ�λ���ӿռ�ǰ�أ��ӿռ�����ܶȣ���ͨ�ӿռ��ڸ���ı���
	void SPMOEA::localCrowdSelection5(size_t select_num, Problem *pro, Random *rnd) {
		auto& front_clusters = getFrontRegionLinkSpace();
		std::vector<size_t> reserve_front(front_clusters.size(), 0);
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<size_t>> layer_inds_inx;//�ֲ㴢���������
			std::vector<size_t> flag(getPop()[i].getOffspring().size(), 0);
			std::vector<size_t> select_flag(getPop()[i].getOffspring().size(), 0);
			int temp_rank = 0;
			while (std::find(flag.begin(), flag.end(), 0) != flag.end()) {
				std::vector<size_t> temp_inx;
				for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
					if (getPop()[i].getOffspring()[j].fitness() == temp_rank) {
						temp_inx.push_back(j);
						flag[j] = 1;
					}
				}
				temp_rank++;
				layer_inds_inx.emplace_back(temp_inx);
			}
			//�õ��������ӿռ��ӳ��
			std::vector<size_t> space_index;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto& sol = getPop()[i].getOffspring()[j].variable().vect();
				auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				space_index.push_back(inx);
			}
			std::map<size_t, std::vector<size_t>> space_ind;
			std::map<size_t, std::vector<size_t>> space_ind_select_flag;
			for (size_t j = 0; j < space_index.size(); ++j) {
				if (space_ind[space_index[j]].empty()) {
					std::vector<size_t> temp;
					temp.push_back(j);
					space_ind.insert(std::make_pair<>(space_index[j],temp));
				}
				space_ind[space_index[j]].push_back(j);
			}
			for (auto jj : space_ind) {
				std::vector<size_t> temp(jj.second.size(), 0);
				space_ind_select_flag.insert(std::make_pair<>(jj.first,temp));
			}
			//��һ���ռ����
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			Population<Solution<>> residual_pop;
			std::vector<size_t> select_inx;
			std::vector<size_t> residual_inx;
			auto front_spaces = getFrontSpace();
			//ѡ�����
			for (size_t j = 0; j < front_clusters.size(); ++j) {
				for (auto jj : space_ind) {
					if (std::find(front_clusters[j].begin(), front_clusters[j].end(), jj.first) != front_clusters[j].end()) {
						//��ѡһ���ӿռ俿ǰ��
						std::vector<int> ranks;
						std::vector<std::vector<Real>*> objs;
						for (size_t k = 0; k < jj.second.size(); ++k) {
							objs.emplace_back(&all_objs[jj.second[k]]);
						}
						ofec::nd_sort::fastSort<Real>(objs, ranks, CAST_CONOP(pro)->optimizeMode());
						std::vector<size_t> first_ind;
						for (size_t k = 0; k < ranks.size(); ++k) {
							if (ranks[k] == 0) {
								first_ind.push_back(k);
							}
						}
						size_t sele = (size_t)std::floor(first_ind.size()*rnd->uniform.next());
						select_inx.push_back(jj.second[sele]);
						select_flag[jj.second[sele]] = 1;
					}
				}
			}
			if (space_ind.size() <= select_num) {

			}
			else {

			}
			while (select_inx.size() < select_num) {
				for (size_t j = 0; j < layer_inds_inx.size(); ++j) {
					for (size_t k = 0; k < layer_inds_inx[j].size(); ++k) {
						auto ind_space = space_index[layer_inds_inx[j][k]];
					}
					
				}
			}

			
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//�ۺϲ�ͬ�������¸���ı���״̬������ۺ���������ѡ�����ذ������ӿռ�ֲ�ѡ������Ŀ��ռ�ѡ������
	//���߿ռ�ӵ����ѡ�����������Ӵ��ݻ�ѡ����
	void SPMOEA::ensembleSelection(size_t select_num, Problem *pro, Random *rnd) {
		auto& front_clusters = getFrontRegionLinkSpace();
		std::vector<size_t> reserve_front(front_clusters.size(), 0);
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////����archive�з��ظ���
			//std::vector<size_t> add_archive_inx;
			//for (size_t j = 0; j < getArchiveSols().size(); ++j) {
			//	auto ind1 = getArchiveSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1, ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_archive_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_archive_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getArchiveSols()[add_archive_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<size_t>> layer_inds_inx;//�ֲ㴢���������
			std::vector<size_t> flag(getPop()[i].getOffspring().size(), 0);
			int temp_rank = 0;
			while (std::find(flag.begin(), flag.end(), 0) != flag.end()) {
				std::vector<size_t> temp_inx;
				for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
					if (getPop()[i].getOffspring()[j].fitness() == temp_rank) {
						temp_inx.push_back(j);
						flag[j] = 1;
					}
				}
				temp_rank++;
				layer_inds_inx.emplace_back(temp_inx);
			}

			//�õ��������ӿռ��ӳ��
			std::vector<size_t> var_space_index;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto& sol = getPop()[i].getOffspring()[j].variable().vect();
				auto inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				var_space_index.push_back(inx);
			}
			std::map<size_t, std::vector<size_t>> var_space_indi;
			std::map<size_t, std::vector<size_t>> space_ind_select_flag;
			for (size_t j = 0; j < var_space_index.size(); ++j) {
				if (var_space_indi[var_space_index[j]].empty()) {
					std::vector<size_t> temp;
					temp.push_back(j);
					var_space_indi.insert(std::make_pair<>(var_space_index[j], temp));
				}
				var_space_indi[var_space_index[j]].push_back(j);
			}
			for (auto jj : var_space_indi) {
				std::vector<size_t> temp(jj.second.size(), 0);
				space_ind_select_flag.insert(std::make_pair<>(jj.first, temp));
			}
			//��һ���ռ����
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			std::vector<size_t> select_inx;
			/*********************************************************************
			                                 �����ش��
			*********************************************************************/
			//�����ӿռ�������
			std::vector<Real> var_subspace_select_score(getPop()[i].getOffspring().size(),0.);
			for (auto jj : var_space_indi) {
				std::vector<std::vector<Real>> temp_sol;
				for (size_t j = 0; j < jj.second.size(); ++j) {
					temp_sol.emplace_back(all_objs[jj.second[j]]);
				}
				std::vector<int> ranks;
				std::vector<std::vector<Real>*> objs;
				for (size_t j = 0; j < temp_sol.size(); ++j) {
					objs.emplace_back(&temp_sol[j]);
				}
				int layer_num=ofec::nd_sort::fastSort<Real>(objs, ranks, CAST_CONOP(pro)->optimizeMode());
				for (size_t j = 0; j < ranks.size(); ++j) {
					//var_subspace_select_score[jj.second[j]] = (1 - 0.3 * ranks[j]) > 0 ? (1 - 0.3 * ranks[j]) : 0.;
					var_subspace_select_score[jj.second[j]] = (1. - ranks[j] / 5) > 0 ? (1. - ranks[j] / 5) : 0;
				}
			}

			//Ŀ��ռ�������
			std::vector<Real> obj_space_select_score(getPop()[i].getOffspring().size(), 0.);
			////Ŀ���ӿռ�������
			//size_t M = CAST_CONOP(pro)->numberObjectives();
			//size_t num_obj_spaces = numObjRegion();
			//std::vector<Real> subspace_ratios(num_obj_spaces, 1. / num_obj_spaces);
			//std::vector<std::pair<Real, Real>> obj_bound;
			//for (size_t i = 0; i < M; ++i) {
			//	obj_bound.emplace_back(std::make_pair<>(0., 1.));
			//}
			//nanoflann::KDTreeSpace<Real> obj_tree(subspace_ratios, obj_bound);// ����Ŀ���ӿռ�
			//obj_tree.buildIndex();
			//std::vector<size_t> obj_space_index;// ȷ���������ڵ��ӿռ�
			//for (size_t j = 0; j < all_objs.size(); ++j) {
			//	obj_space_index.emplace_back(obj_tree.getRegionIdx(all_objs[j]));
			//}
			//std::map<size_t, std::vector<size_t>> obj_space_indi;//�õ��ӿռ�͸����ӳ��
			//std::vector<std::vector<size_t>> obj_subspace_layer;//���ݸ���λ�û����ӿռ�㼶
			//std::vector<size_t> deal_flag(all_objs.size(), 0);
			//while (std::find(deal_flag.begin(), deal_flag.end(), 0) != deal_flag.end()) {
			//	//��ʣ������ǰ�ţ�ȷ���ӿؼ��㼶
			//	std::vector<size_t> residual_ind_inx;
			//	for (size_t i = 0; i < deal_flag.size(); ++i) {
			//		if (deal_flag[i] == 0) {
			//			residual_ind_inx.push_back(i);
			//		}
			//	}
			//	std::vector<std::vector<Real>> new_ind;
			//	for (size_t i = 0; i < residual_ind_inx.size(); ++i) {
			//		new_ind.emplace_back(all_objs[residual_ind_inx[i]]);
			//	}
			//	//ʣ�����ǰ������
			//	auto f_inx = getNondominatedSetIndex(new_ind, CAST_CONOP(pro)->optimizeMode());
			//	std::vector<size_t> new_front_ind_inx;
			//	for (size_t i = 0; i < f_inx.size(); ++i) {
			//		new_front_ind_inx.push_back(residual_ind_inx[f_inx[i]]);
			//	}
			//	//ʣ������ǰ�������ӿռ�
			//	std::vector<size_t> new_front_space_inx;
			//	for (size_t i = 0; i < new_front_ind_inx.size(); ++i) {
			//		size_t inx = obj_space_index[new_front_ind_inx[i]];
			//		if (std::find(new_front_space_inx.begin(), new_front_space_inx.end(), inx) == new_front_space_inx.end()) {
			//			new_front_space_inx.emplace_back(inx);
			//		}
			//	}
			//	//�µ��ӿռ��
			//	for (size_t i = 0; i < new_front_space_inx.size(); ++i) {
			//		std::vector<size_t> temp_space;
			//		for (size_t j = 0; j < obj_space_index.size(); ++j) {
			//			if (deal_flag[j] == 0) {
			//				if (obj_space_index[j] == new_front_space_inx[i]) {
			//					temp_space.push_back(j);
			//					deal_flag[j] = 1;
			//				}
			//			}
			//		}
			//		obj_space_indi.insert(std::make_pair<>(new_front_space_inx[i], temp_space));
			//	}
			//	obj_subspace_layer.emplace_back(new_front_space_inx);
			//}
			//for (auto jj : obj_space_indi) {
			//	std::vector<std::vector<Real>> temp_obj;
			//	for (size_t j = 0; j < jj.second.size(); ++j) {
			//		temp_obj.emplace_back(all_objs[jj.second[j]]);
			//	}
			//	std::vector<int> ranks;
			//	std::vector<std::vector<Real>*> objs;
			//	for (size_t j = 0; j < temp_obj.size(); ++j) {
			//		objs.emplace_back(&temp_obj[j]);
			//	}
			//	ofec::nd_sort::fastSort<Real>(objs, ranks, CAST_CONOP(pro)->optimizeMode());
			//	size_t space_layer=0;
			//	for (size_t j = 0; j < obj_subspace_layer.size(); ++j) {
			//		if (std::find(obj_subspace_layer[j].begin(), obj_subspace_layer[j].end(), jj.first) != obj_subspace_layer[j].end()) {
			//			space_layer = j;
			//		}
			//	}
			//	for (size_t j = 0; j < ranks.size(); ++j) {
			//		obj_space_select_score[jj.second[j]] = ((1-0.5*space_layer) - 0.2 * ranks[j]) > 0 ? ((1 - 0.5 * space_layer) - 0.2 * ranks[j]) : 0.;
			//	}
			//}
			//ֱ�Ӱ�����ֵ���
			size_t layer_index = 0;
			size_t count = 0;
			for (size_t j = 0; j < layer_inds_inx.size(); ++j) {
				if (count < select_num) {
					count += layer_inds_inx[j].size();
				}
				else {
					layer_index = j;
					break;
				}
			}
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto rank = getPop()[i].getOffspring()[j].fitness();
				obj_space_select_score[j] = (1.- rank/5.)>0? (1. - rank / 5.):0;
				//if (rank < layer_index) {
				//	obj_space_select_score[j] = 1.;
				//}
				//else {
				//	//obj_space_select_score[j] = (1 - 0.2 * rank) > 0 ? (1 - 0.2 * rank) : 0.;
				//	obj_space_select_score[j] = (1. - (rank-layer_index)/5.) > 0 ? (1. - (rank - layer_index) / 5.) : 0;
				//}
			}

			//���߿ռ�ӵ���ȴ�֣������ӿռ��ܶ�
			std::vector<Real> var_space_crowdist_select_score1(getPop()[i].getOffspring().size(), 0.);
			std::vector<Real> var_min_dist;//ƽ��ӵ������
			size_t kk = 2;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				Real temp = 0.;
				for (size_t k = 0; k < kk; ++k) {
					auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
					temp += min_d;
					auto inx = std::distance(var_dist[j].begin(), std::min_element(var_dist[j].begin(), var_dist[j].end()));
					var_dist[j][inx] = INT16_MAX;
				}
				var_min_dist.push_back(temp/kk);
			}
			Real min_min_var_dist= *std::min_element(var_min_dist.begin(), var_min_dist.end());
			Real max_min_var_dist = *std::max_element(var_min_dist.begin(), var_min_dist.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				var_space_crowdist_select_score1[j] = (var_min_dist[j]- min_min_var_dist) / (max_min_var_dist- min_min_var_dist);
			}

			/*std::vector<Real> var_min_dist_sort;
			for (size_t j = 0; j < var_min_dist.size(); ++j) {
				size_t count = 0;
				for (size_t k = 0; k < var_min_dist.size(); ++k) {
					if (k != j) {
						if (var_min_dist[k] > var_min_dist[j]) {
							count++;
						}

					}
				}
				var_min_dist_sort.push_back(count);
			}
			Real max_min_var_dist = *std::max_element(var_min_dist_sort.begin(), var_min_dist_sort.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				var_space_crowdist_select_score1[j] = 1 - var_min_dist_sort[j]/ max_min_var_dist;
			}*/

			std::vector<Real> var_space_crowdist_select_score2(getPop()[i].getOffspring().size(), 0.);
			std::vector<Real> space_density;//�ӿռ��ܶ�
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto& sol = getPop()[i].getOffspring()[j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				auto volume = getMO_HLC().subspaceTree().getBoxVolume(space_inx);
				space_density.push_back(getMO_HLC().getSubspaceInfo(space_inx).m_history_inds.size()/volume);
			}
			std::vector<Real> var_density_sort;
			for (size_t j = 0; j <space_density.size(); ++j) {
				size_t count = 0;
				for (size_t k = 0; k < space_density.size(); ++k) {
					if (k != j) {
						if (space_density[k] < space_density[j]) {
							count++;
						}

					}
				}
				var_density_sort.push_back(count);
			}
			Real max_density_sort = *std::max_element(var_density_sort.begin(), var_density_sort.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				var_space_crowdist_select_score2[j] = 1 - var_density_sort[j] / max_density_sort;
			}

			//Ŀ��ռ�ӵ���ȴ��
			std::vector<Real> obj_space_crowdist_select_score(getPop()[i].getOffspring().size(), 0.);
			std::vector<Real> obj_min_dist;
			size_t oo = 2;
			for (size_t j = 0; j < obj_dist.size(); ++j) {
				Real temp = 0.;
				for (size_t k = 0; k < oo; ++k) {
					auto min_d = *std::min_element(obj_dist[j].begin(), obj_dist[j].end());
					temp += min_d;
					auto inx = std::distance(obj_dist[j].begin(), std::min_element(obj_dist[j].begin(), obj_dist[j].end()));
					obj_dist[j][inx] = INT16_MAX;
				}
				obj_min_dist.push_back(temp / oo);
			}
			Real min_min_obj_dist = *std::min_element(obj_min_dist.begin(), obj_min_dist.end());
			Real max_min_obj_dist = *std::max_element(obj_min_dist.begin(), obj_min_dist.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				obj_space_crowdist_select_score[j] = (obj_min_dist[j] - min_min_obj_dist) / (max_min_obj_dist - min_min_obj_dist);
			}

			/*std::vector<Real> obj_min_dist_sort;
			for (size_t j = 0; j < obj_min_dist.size(); ++j) {
				size_t count = 0;
				for (size_t k = 0; k < obj_min_dist.size(); ++k) {
					if (k != j) {
						if (obj_min_dist[k] > obj_min_dist[j]) {
							count++;
						}

					}
				}
				obj_min_dist_sort.push_back(count);
			}
			Real max_min_obj_dist = *std::max_element(obj_min_dist_sort.begin(), obj_min_dist_sort.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				obj_space_crowdist_select_score[j] = 1 - obj_min_dist_sort[j] / max_min_obj_dist;
			}*/

			//����ѡ����
			std::vector<Real> evolve_select_score(getPop()[i].getOffspring().size(), 0.);
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					evolve_select_score[j]=0.;
					evolve_select_score[getPop()[i].size()+j]=1.;
				}
				else if (dominance_ship == Dominance::kDominated) {
					evolve_select_score[j] =1.;
					evolve_select_score[getPop()[i].size() + j] =0.;
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					evolve_select_score[j] =0.5;
					evolve_select_score[getPop()[i].size() + j] =0.5;
				}
			}
			for (size_t j = 2*getPop()[i].size(); j < getPop()[i].getOffspring().size(); ++j) {
				evolve_select_score[j] = 1.;
			}

			//�ۺϸ�����ѡ��
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				//var_subspace_select_score[j] = 0.;
				//obj_space_select_score[j] = 0.;
				//var_space_crowdist_select_score[j] = 0.;
				obj_space_crowdist_select_score[j] = 0.;
				evolve_select_score[j] = 0.;
			}	

			std::vector<Real> all_scores(getPop()[i].getOffspring().size(), 0.);
			for (size_t j = 0; j < all_scores.size(); ++j) {
				all_scores[j] = 0.25*var_subspace_select_score[j]+0.25*obj_space_select_score[j] + 0.25 * var_space_crowdist_select_score1[j]+ 0.*var_space_crowdist_select_score2[j]+ 0.25*obj_space_crowdist_select_score[j] + 0.1*evolve_select_score[j];
			}
			//����ֵ�ɸߵ���ѡ��
			std::vector<size_t> select_flag(all_scores.size(), 0);
			for (size_t j = 0; j < all_scores.size(); ++j) {
				if (select_inx.size() >= select_num) {
					break;
				}
				else {
					Real max_v = 0.;
					size_t inx=0;
					for (size_t k = 0; k < all_scores.size(); ++k) {
						if (select_flag[k] == 0) {
							if (all_scores[k] > max_v) {
								max_v = all_scores[k];
								inx = k;
							}
						}
					}
					select_inx.push_back(inx);
					select_flag[inx] = 1;
				}
			}

			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//��������ѡ�񣬷�֧��ʱ���ѡ���Ӵ�����ʱ�����ݸ����Ƿ�Ϊǰ�ؾ����Ӵ������ĸ���
	void SPMOEA::SolutionSelection(Problem* pro, Random* rnd) {
		//�������ѡ����Ҫ��������ѡ��ķֲ�
		//�ۺϵ�ǰ��Ⱥ���Ӵ���Ⱥ����ʷ��֧�����linkspace�ķֲ�����ѡ��
		//�ӿռ������ѡһ����Ȼ������ѡ����linkspace�ӿռ���ѡ����Զ��һ����ѡ��
		auto& his_front_sols = getHisFrontSols();
		//�����н�ֲ����ӿռ��Ƿ�Ϊǰ���ӿռ������䣬���ڸ����ӿռ����ѡ��һ����ǰ��ǰ�Ľ⣬Ȼ��ѡ���ǰ���ӿռ��нϺõĽ⣬
		// ����������α���ǰ���ӿռ��linkspace,ѡ������ѡ����Զ��һ��
		
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			Population<Solution<>> residual_pop;
			std::vector<size_t> sele_inx;
			std::vector<size_t> residual_inx;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				bool flag1 = false;//���������Ƿ�Ϊǰ�ظ���
				auto& sol1 = getPop()[i][j].variable().vect();
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					auto sol2 = his_front_sols[k]->variable().vect();
					if (ifSame(sol1, sol2)) {
						flag1 = true;
						break;
					}
				}
				bool flag2 = false;//�Ӵ������Ƿ�Ϊǰ�ظ���
				auto& sol2 = getPop()[i][j].variable().vect();
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					auto sol3 = his_front_sols[k]->variable().vect();
					if (ifSame(sol2, sol3)) {
						flag2 = true;
						break;
					}
				}
				
				if (dominance_ship == Dominance::kDominant) {//����֧���Ӵ�
					//���ݸ����Ƿ�Ϊǰ�ؽ�����滻����
					//Real rand = rnd->uniform.next();
					/*if (flag1) {
						if (rand < 0.1) {
							sele_inx.push_back(j);
						}
						else {
							sele_inx.push_back(getPop()[i].size() + j);
						}
					}
					else {
						if (rand < 0.2) {
							sele_inx.push_back(j);
						}
						else {
							sele_inx.push_back(getPop()[i].size() + j);
						}
					}*/
					sele_inx.push_back(getPop()[i].size() + j);
				}
				else if (dominance_ship == Dominance::kDominated) {//�Ӵ�֧�丸��
					/*if (flag2) {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					else {
						sele_inx.push_back(j);
					}*/
					sele_inx.push_back(j);
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					auto pp_sol = getPop()[i][j].variable().vect();
					size_t p_space = getMO_HLC().subspaceTree().getRegionIdx(pp_sol);
					auto off_sol = getPop()[i].getOffspring()[j].variable().vect();
					size_t o_space = getMO_HLC().subspaceTree().getRegionIdx(off_sol);
					Real p_span = 0.;
					Real o_span = 0.;
					auto p_box = getMO_HLC().subspaceTree().getBox(p_space);
					for (size_t j = 0; j < p_box.size(); ++j) {
						p_span += (p_box[j].second - p_box[j].first);
					}
					p_span /= p_box.size();
					auto o_box = getMO_HLC().subspaceTree().getBox(o_space);
					for (size_t j = 0; j < o_box.size(); ++j) {
						o_span += (o_box[j].second - o_box[j].first);
					}
					o_span /= o_box.size();

					if (flag1 && flag2) {
						//�Ƚ�ǰ�ؽ�ĸ���
						Real density1 = getMO_HLC().getSubspaceInfo(p_space).m_front_sol_in_subspace.size() / p_span;
						Real density2 = getMO_HLC().getSubspaceInfo(o_space).m_front_sol_in_subspace.size() / o_span;
						if (density1 < density2) {
							sele_inx.push_back(getPop()[i].size() + j);
						}
						else {
							sele_inx.push_back(j);
						}
					}
					else if (flag1 && (!flag2)) {
						sele_inx.push_back(j);
					}
					else if (flag2 && (!flag1)) {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					else {
						/*Real density1 = getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size() / p_span;
						Real density2 = getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / o_span;
						if (density1 < density2) {
							sele_inx.push_back(getPop()[i].size() + j);
						}
						else {
							sele_inx.push_back(j);
						}*/
						size_t rank1 = getMO_HLC().getSubspaceInfo(p_space).m_best_rank;
						size_t rank2 = getMO_HLC().getSubspaceInfo(o_space).m_best_rank;
						if (rank1 < rank2) {
							Real rand = rnd->uniform.next();
							if (rand < 0.2) {
								sele_inx.push_back(j);
							}
							else {
								sele_inx.push_back(getPop()[i].size() + j);
							}
							//sele_inx.push_back(getPop()[i].size() + j);
						}
						else if (rank1 > rank2) {
							Real rand = rnd->uniform.next();
							if (rand > 0.2) {
								sele_inx.push_back(j);
							}
							else {
								sele_inx.push_back(getPop()[i].size() + j);
							}
							//sele_inx.push_back(j);
						}
						else {
							Real rand = rnd->uniform.next();
							if (rand > 0.5) {
								sele_inx.push_back(j);
							}
							else {
								sele_inx.push_back(getPop()[i].size() + j);
							}
						}
					}
				}
			}
			//ʹ��ѡ�����½���¸����ݻ��켣
			for (size_t j = 0; j < sele_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
			}
		}
	}

	//���ڸ����K����ѡ��
	void SPMOEA::knnSelection(Problem* pro, Random* rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}

		}
	}

	//��������ָ���ѡ�񣬱Ƚϸ������Ӵ���HV����ֵѡ��
	void SPMOEA::indicatorSelection(Problem *pro, Random *rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			//��һ�������ľ���ֵ��Ŀ��ֵ
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������Ⱥ�����ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < all_vars.size(); ++j) {
				auto p1 = all_vars[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_vars[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_vars.size(); ++j) {
				for (size_t k = 0; k < all_vars.size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
					}
					else if (k == j) {
						var_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_var_min_dist = 0.;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_var_min_dist += min_d;
			}
			mean_var_min_dist /= var_dist.size();
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> all_objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(all_objs);
			//������ȺĿ��ռ���С�����ƽ��ֵ
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto p1 = all_objs[j];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k <= j) {
						temp_dist.push_back(0.);
					}
					else {
						auto p2 = all_objs[k];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp_dist.push_back(dist);
					}
				}
				obj_dist.emplace_back(temp_dist);
			}
			for (size_t j = 0; j < all_objs.size(); ++j) {
				for (size_t k = 0; k < all_objs.size(); ++k) {
					if (k < j) {
						obj_dist[j][k] = obj_dist[k][j];
					}
					else if (k == j) {
						obj_dist[j][k] = INT16_MAX;
					}
				}
			}
			Real mean_obj_min_dist = 0.;
			for (size_t j = 0; j < all_objs.size(); ++j) {
				auto min_d = *std::min_element(all_objs[j].begin(), all_objs[j].end());
				mean_obj_min_dist += min_d;
			}
			mean_obj_min_dist /= all_objs.size();
			Population<Solution<>> residual_pop;
			std::vector<size_t> select_inx;
			std::vector<size_t> residual_inx;
			//���ѡ�����
			std::vector<Real> reference_point(CAST_CONOP(pro)->numberObjectives(), 1.01);
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				//�Ƚϸ������Ӵ���hvֵ���׶�
				Real p_hv = 1., o_hv = 1.;
				for (size_t k = 0; k < reference_point.size(); ++k) {
					p_hv *= (reference_point[k]-all_objs[getPop()[i].size()+j][k]);
					o_hv *= (reference_point[k] - all_objs[j][k]);
				}
				if (p_hv > o_hv) {
					select_inx.push_back(getPop()[i].size() + j);
				}
				else {
					select_inx.push_back(j);
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//�����Ӵ�Ŀ���ӿռ�ѡ��
	void SPMOEA::subspaceSelection(Problem *pro, Random *rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			SPMOEA::NDSort(getPop()[i].getOffspring());
			Population<Solution<>> temp_pop;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				temp_pop.append(getPop()[i].getOffspring()[j]);
			}
			auto select_inx = selectFromSpace(temp_pop, getPop()[i].size(), pro,rnd);
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//���ڸ���ѡ�񣬷�֧��ĵ�ѡ������ѡ��������ĵ�
	void SPMOEA::localIndSelection(Problem *pro,Random *rnd) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			SPMOEA::NDSort(getPop()[i].getOffspring());
			Population<Solution<>> residual_pop;
			std::vector<size_t> sele_inx;
			std::vector<size_t> residual_inx;

			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				auto& pp = getPop()[i][j].objective();
				auto& off = getPop()[i].getOffspring()[j].objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					sele_inx.push_back(getPop()[i].size() + j);
				}
				else if (dominance_ship == Dominance::kDominated) {
					sele_inx.push_back(j);
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					//��������������ѡ��Ĺ�ϵ
					if (getPop()[i].getOffspring()[j].fitness() < getPop()[i].getOffspring()[getPop()[i].size() + j].fitness()) {
						sele_inx.push_back(j);
					}
					else if (getPop()[i].getOffspring()[j].fitness() > getPop()[i].getOffspring()[getPop()[i].size() + j].fitness()) {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					else {
						if (sele_inx.empty()) {
							Real rand = rnd->uniform.next();
							if (rand > 0.5) {
								sele_inx.push_back(getPop()[i].size() + j);
							}
							else {
								sele_inx.push_back(j);
							}
						}
						else {

						}
					}
					
					residual_inx.push_back(j);
					residual_inx.push_back(getPop()[i].size() + j);
					residual_pop.append(getPop()[i][j]);
					residual_pop.append(getPop()[i].getOffspring()[j]);
				}
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			SPMOEA::NDSort(residual_pop);
			std::vector<std::vector<size_t>> residual_layer_inds_inx;//�ֲ㴢���������
			std::vector<size_t> flag(residual_pop.size(), 0);
			int temp_rank = 0;
			while (std::find(flag.begin(), flag.end(), 0) != flag.end()) {
				std::vector<size_t> temp_inx;
				for (size_t j = 0; j < residual_pop.size(); ++j) {
					if (residual_pop[j].fitness() == temp_rank) {
						temp_inx.push_back(j);
						flag[j] = 1;
					}
				}
				temp_rank++;
				residual_layer_inds_inx.emplace_back(temp_inx);
			}
			//��һ�����߱���ֵ
			std::vector<std::vector<Real>> all_vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(all_vars);
			//������ѡ�񣬲��Ҵ�ѡ������ѡ���������ռ�ľ��벻��̫��
			std::vector<std::vector<size_t>> layer_inds_status;//�������ѡ��״̬
			for (size_t j = 0; j < residual_layer_inds_inx.size(); ++j) {
				std::vector<size_t> temp(residual_layer_inds_inx[j].size(), 0);
				layer_inds_status.emplace_back(temp);
			}
			//ʣ��������е��������ռ��ڵľ���
			std::vector<std::vector<Real>> var_dist;
			for (size_t j = 0; j < residual_inx.size(); ++j) {
				auto p1 = all_vars[residual_inx[j]];
				std::vector<Real> temp_dist;
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					auto p2 = all_vars[k];
					auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
					if (dist == 0.) {
						temp_dist.push_back(INT16_MAX);
					}
					else {
						temp_dist.push_back(dist);
					}
				}
				var_dist.emplace_back(temp_dist);
			}
			Real niche_radius = INT16_MAX;
			size_t num_var = CAST_CONOP(pro)->numberVariables();
			//������Ⱥ��С�����ƽ��ֵ
			Real mean_min_dist = 0.;
			for (size_t j = 0; j < residual_inx.size(); ++j) {
				auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
				mean_min_dist += min_d;
			}
			mean_min_dist /= var_dist.size();
			//ѭ��ѡ������ѡ��Զ���ҿ�ǰ�ĸ���
			niche_radius = mean_min_dist;
			niche_radius *= 1.1;
			while (sele_inx.size() < getPop()[i].size()) {
				niche_radius /= 1.1;
				for (size_t j = 0; j < residual_layer_inds_inx.size(); ++j) {
					for (size_t k = 0; k < residual_layer_inds_inx[j].size(); ++k) {
						if (sele_inx.empty()) {
							sele_inx.push_back(residual_inx[residual_layer_inds_inx[j][k]]);
							layer_inds_status[j][k] = 1;
						}
						else if (layer_inds_status[j][k] == 0) {
							//��������ѡ��ľ���
							bool sele_flag = true;
							for (size_t p = 0; p < sele_inx.size(); ++p) {
								if (var_dist[residual_layer_inds_inx[j][k]][sele_inx[p]] < niche_radius) {
									sele_flag = false;
									break;
								}
							}
							if (sele_flag) {
								sele_inx.push_back(residual_inx[residual_layer_inds_inx[j][k]]);
								layer_inds_status[j][k] = 1;
							}
						}
						if (sele_inx.size() >= getPop()[i].size()) {
							break;
						}
					}
					if (sele_inx.size() >= getPop()[i].size()) {
						break;
					}
				}
			}
			//ʹ��ѡ�����½���¸����ݻ��켣

			for (size_t j = 0; j < sele_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
			}
		}
	}

	//Ŀ��ռ�;��߿ռ�ӵ������ѡ��
	void SPMOEA::doubleCrowdSelection(Problem *pro) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ط��ظ���
			//std::vector<size_t> add_his_inx;
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	auto ind1 = getHisFrontSols()[j]->variable().vect();
			//	bool sele_flag = true;
			//	for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
			//		auto ind2 = getPop()[i].getOffspring()[k].variable().vect();
			//		if (ifSame(ind1,ind2)) {
			//			sele_flag = false;
			//			break;
			//		}
			//	}
			//	if (sele_flag) {
			//		add_his_inx.push_back(j);
			//	}
			//}
			//for (size_t j = 0; j < add_his_inx.size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[add_his_inx[j]]);
			//}
			//Ŀ��ռ��κ������ռ�����ѡ��
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<size_t>> layer_inds_inx;//�ֲ㴢���������
			std::vector<size_t> flag(getPop()[i].getOffspring().size(), 0);
			int temp_rank = 0;
			while (std::find(flag.begin(), flag.end(), 0) != flag.end()) {
				std::vector<size_t> temp_inx;
				for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
					if (getPop()[i].getOffspring()[j].fitness() == temp_rank) {
						temp_inx.push_back(j);
						flag[j] = 1;
					}
				}
				temp_rank++;
				layer_inds_inx.emplace_back(temp_inx);
			}
			//��һ��Ŀ��ֵ
			std::vector<std::vector<Real>> objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				objs.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			dataNormalize(objs);
			//��һ�����߱���ֵ
			std::vector<std::vector<Real>> vars;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				objs.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			dataNormalize(vars);
			//������ѡ�񣬲��Ҵ�ѡ������ѡ���������ռ�ľ��벻��̫��
			std::vector<size_t> select_inx;
			std::vector<std::vector<size_t>> layer_inds_status;//�������ѡ��״̬
			for (size_t j = 0; j < layer_inds_inx.size(); ++j) {
				std::vector<size_t> temp(layer_inds_inx[j].size(), 0);
				layer_inds_status.emplace_back(temp);
			}
			//�ҳ���ǰǰ���ӿռ��е����ά��ֵ����Ϊ�����ռ��ӵ���뾶
			Real niche_radius = INT16_MAX;
			size_t num_var = CAST_CONOP(pro)->numberVariables();
			for (size_t j = 0; j < getFrontSpace().size(); ++j) {
				auto box = getMO_HLC().subspaceTree().getBox(getFrontSpace()[j]);
				for (size_t k = 0; k < num_var; ++k) {
					Real dim_span = box[k].second - box[k].first;
					if (niche_radius > dim_span) {
						niche_radius = dim_span;
					}
				}
			}
			//���е�֮��ľ���
			std::vector<std::vector<Real>> var_dist;
			std::vector<std::vector<Real>> obj_dist;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto p1 = getPop()[i].getOffspring()[j].variable().vect();
				auto pp1 = getPop()[i].getOffspring()[j].objective();
				std::vector<Real> temp_dist1;
				std::vector<Real> temp_dist2;
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					if (k <= j) {
						temp_dist1.push_back(0.);
						temp_dist2.push_back(0.);
					}
					else {
						auto p2 = getPop()[i].getOffspring()[k].variable().vect();
						auto pp2 = getPop()[i].getOffspring()[k].objective();
						auto dist1 = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						auto dist2 = euclideanDistance(pp1.begin(), pp1.end(), pp2.begin());
						temp_dist1.push_back(dist1);
						temp_dist2.push_back(dist2);
					}
				}
				var_dist.emplace_back(temp_dist1);
				obj_dist.emplace_back(temp_dist2);
			}
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					if (k < j) {
						var_dist[j][k] = var_dist[k][j];
						obj_dist[j][k] = obj_dist[k][j];
					}
				}
			}
			//ѭ��ѡ��
			niche_radius *= 1.5;
			while (select_inx.size() < getPop()[i].size()) {
				niche_radius /= 1.5;
				for (size_t j = 0; j < layer_inds_inx.size(); ++j) {
					for (size_t k = 0; k < layer_inds_inx[j].size(); ++k) {
						if (select_inx.empty()) {
							select_inx.push_back(layer_inds_inx[j][k]);
							layer_inds_status[j][k] = 1;
						}
						else if (layer_inds_status[j][k] == 0) {
							//��������ѡ��ľ���
							bool sele_flag = true;
							for (size_t p = 0; p < select_inx.size(); ++p) {
								if (var_dist[layer_inds_inx[j][k]][select_inx[p]] < niche_radius) {
									sele_flag = false;
									break;
								}
							}
							if (sele_flag) {
								select_inx.push_back(layer_inds_inx[j][k]);
								layer_inds_status[j][k] = 1;
							}
						}
						if (select_inx.size() >= getPop()[i].size()) {
							break;
						}
					}
					if (select_inx.size() >= getPop()[i].size()) {
						break;
					}
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	//��Ŀ��ѡ��֧��������͹�һ����ӵ��ֵ
	void SPMOEA::multiObjSelection(Problem *pro) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			////������ʷǰ�ؽ�
			//for (size_t j = 0; j < getHisFrontSols().size(); ++j) {
			//	getPop()[i].getOffspring().append(*getHisFrontSols()[j]);
			//}
			//���ö�Ŀ������ķ�ʽѡ��: rankֵ�͹�һ�������ռ��
			SPMOEA::NDSort(getPop()[i].getOffspring());
			//��������dominance��
			std::vector<size_t> dominance_num(getPop()[i].getOffspring().size(), 0);
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
					auto dominaceship = objectiveCompare(getPop()[i].getOffspring()[j].objective(), getPop()[i].getOffspring()[k].objective(), CAST_CONOP(pro)->optimizeMode());
					if (dominaceship == Dominance::kDominant) {
						dominance_num[j]++;
					}
				}
			}
			//���������Ŀ��ռ��niche countֵ
			std::vector<std::vector<Real>> objs;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				objs.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
			}
			//��һ��
			dataNormalize(objs);
			std::vector<Real> min_mean_dist(objs.size(), 0.);
			size_t niche_count = 10;
			for (size_t j = 0; j < objs.size(); ++j) {
				std::vector<Real> all_dist;
				for (size_t k = 0; k < objs.size(); ++k) {
					auto dist = euclideanDistance(objs[j].begin(), objs[j].end(), objs[k].begin());
					all_dist.push_back(dist);
				}
				std::sort(all_dist.begin(), all_dist.end());//ȡ��С�ļ�������ֵ
				Real sum = 0.;
				for (size_t k = 0; k < niche_count; ++k) {
					sum += all_dist[k];
				}
				min_mean_dist[j] = sum / niche_count;
			}
			//����2��Ŀ��
			std::vector<std::vector<Real>> new_objs(min_mean_dist.size());
			for (size_t j = 0; j < min_mean_dist.size(); ++j) {
				//new_objs[j].push_back(1./(dominance_num[j]+1));
				new_objs[j].push_back((Real)getPop()[i].getOffspring()[j].fitness());
				new_objs[j].push_back(1. / min_mean_dist[j]);
			}
			//����ָ����������
			std::vector<int> rank(new_objs.size());
			nd_sort::fastSort(new_objs, rank, CAST_CONOP(pro)->optimizeMode());
			//���¸���
			std::vector<std::vector<Real>> first_pop;
			std::vector<size_t> first_inx;
			for (size_t j = 0; j < rank.size(); ++j) {
				if (rank[j] == 0) {
					first_pop.emplace_back(new_objs[j]);
					first_inx.push_back(j);
				}
			}
			if (first_pop.size() > getPop()[i].size()) {
				auto select_inx = selectMaxMinFromFront(first_pop, getPop()[i].size());
				for (size_t j = 0; j < select_inx.size(); ++j) {
					getPop()[i][j] = getPop()[i].getOffspring()[first_inx[select_inx[j]]];
				}
			}
			else {
				size_t count = 0;
				for (size_t j = 0; j < first_pop.size(); ++j) {
					getPop()[i][j] = getPop()[i].getOffspring()[first_inx[j]];
					count++;
				}
				int temp_rank = 1;
				while (count < getPop()[i].size()) {
					for (size_t j = 0; j < rank.size(); ++j) {
						if (rank[j] == temp_rank) {
							getPop()[i][count] = getPop()[i].getOffspring()[j];
							count++;
						}
						if (count >= getPop()[i].size()) {
							break;
						}
					}
					if (count < getPop()[i].size()) {
						temp_rank++;
					}
				}
			}
		}
	}

	std::vector<size_t> SPMOEA::neighLocalSelection(size_t neighbor_num,Problem* pro) {
		//1���ֲ�ѡ�������Ⱥ
		std::vector<size_t> candidates;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].getOffspring().size() - getPop()[i].size() + j] = getPop()[i][j];
			}
			NDSort(getPop()[i].getOffspring());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto nn = getPop()[i].getOffspring().nearestNeighbour(j, pro, neighbor_num);
				std::vector<size_t> ind_inx;
				for (auto& n : nn)
					ind_inx.push_back(n.second);//�����������Ⱥ�е�����
				ind_inx.push_back(j);
				//k���ڱȽ�
				std::vector<std::vector<Real>*> objs;
				for (size_t k = 0; k < ind_inx.size(); ++k) {
					objs.emplace_back(&getPop()[i].getOffspring()[ind_inx[k]].objective());
				}
				std::vector<int> rank;
				ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
				for (size_t k = 0; k < ind_inx.size(); ++k) {
					if (rank[k] == 0) {
						if (std::find(candidates.begin(), candidates.end(), ind_inx[k]) == candidates.end()) {
							candidates.push_back(ind_inx[k]);
						}
					}
				}
			}

			std::vector<std::vector<Real>> candidate_vars;
			for (size_t j = 0; j < candidates.size(); ++j) {
				candidate_vars.emplace_back(getPop()[i].getOffspring()[candidates[j]].variable().vect());
			}
			std::vector<std::vector<Real>> residual_vars;
			std::vector<size_t> residuals;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				if (std::find(candidates.begin(), candidates.end(), j) == candidates.end()) {
					residuals.push_back(j);
					residual_vars.emplace_back(getPop()[i].getOffspring()[j].variable().vect());
				}
			}
			while (candidates.size() > m_pop_size) {
				//����ȥ����ӵ���ĵ�
				auto crowd_pair = crowdedPointPair(candidate_vars);
				//�Ƚ�ӵ�����Ŀ��ֵ����̭���
				size_t delete_inx = 0;
				auto obj1 = getPop()[i].getOffspring()[candidates[crowd_pair.first]].objective();
				auto obj2 = getPop()[i].getOffspring()[candidates[crowd_pair.second]].objective();
				auto opt_mode = CAST_CONOP(pro)->optimizeMode();
				auto ship = objectiveCompare(obj1,obj2,opt_mode);
				if (ship==Dominance::kDominant) {
					delete_inx = crowd_pair.second;
				}
				else if (ship == Dominance::kDominated) {
					delete_inx = crowd_pair.first;
				}
				else if (ship==Dominance::kNonDominated) {
					size_t rank1 = getPop()[i].getOffspring()[candidates[crowd_pair.first]].fitness();
					size_t rank2 = getPop()[i].getOffspring()[candidates[crowd_pair.second]].fitness();
					if (rank1>rank2) {
						delete_inx = crowd_pair.first;
					}
					else if(rank1 < rank2) {
						delete_inx = crowd_pair.second;
					}
					else {
						//�Ƚϸ���k���ڵľ���ֵ
						auto nn1 = getPop()[i].getOffspring().nearestNeighbour(candidates[crowd_pair.first], pro, neighbor_num);
						Real neigh_dist1=0.;
						for (auto& n : nn1)
							neigh_dist1+=n.first;//�������
						auto nn2 = getPop()[i].getOffspring().nearestNeighbour(candidates[crowd_pair.second], pro, neighbor_num);
						Real neigh_dist2=0.;
						for (auto& n : nn2)
							neigh_dist2 += n.first;//�������
						if (neigh_dist1 < neigh_dist2) {
							delete_inx = crowd_pair.first;
						}
						else {
							delete_inx = crowd_pair.second;
						}
					}
				}
				candidate_vars.erase(candidate_vars.begin() + delete_inx);
				candidates.erase(candidates.begin() + delete_inx);
			}
			while (candidates.size() < m_pop_size) {
				//���μ�������ѡ����Զ�ĵ�
				auto inx = selectFarthestPoint(candidate_vars, residual_vars);
				candidate_vars.emplace_back(residual_vars[inx]);
				residual_vars.erase(residual_vars.begin() + inx);
				candidates.push_back(residuals[inx]);
				residuals.erase(residuals.begin() + inx);
			}
			/*for (size_t j = 0; j < m_pop_size; ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[candidates[j]];
			}*/
			return candidates;
		}
	}

	std::vector<size_t> SPMOEA::getNonlapNeighSpaces(std::vector<size_t>& spaces,std::vector<size_t> &cur_inx, std::vector<size_t> &dealt_spaces) {
		std::vector<size_t> all_neighs;
		for (size_t i = 0; i < cur_inx.size(); ++i) {
			auto neighs = getMO_HLC().getSubspaceInfo(spaces[cur_inx[i]]).m_sub_neighbors;
			for (auto ii : neighs) {
				all_neighs.push_back(ii);
			}
		}
		std::vector<size_t> out_neighs;
		for (size_t i = 0; i < all_neighs.size(); ++i) {
			if (std::find(dealt_spaces.begin(), dealt_spaces.end(), all_neighs[i]) == dealt_spaces.end()) {
				out_neighs.push_back(all_neighs[i]);
			}
		}
		return all_neighs;
	}

	std::map<size_t, Real> SPMOEA::frontSpaceLinkProbability(size_t inx1) {
		//�����ӿռ��е�ǰ�ؽ�ķֲ�
		auto front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
		//�ҳ���ǰ�������ӿռ�
		std::vector<size_t> nei;
		auto neigh = getMO_HLC().getSubspaceInfo(inx1).m_sub_neighbors;
		auto front_spaces = getFrontSpace();
		for (auto jj : neigh) {
			if (std::find(front_spaces.begin(), front_spaces.end(), jj) != front_spaces.end()) {
				nei.push_back(jj);
			}
		}
		std::vector<Real> all_min_dist;//�ӿռ�֮�����̾���
		Real min_dist = (Real)INT16_MAX;//�ӿռ�֮�����̾����е���С����
		for (size_t i = 0; i < nei.size(); ++i) {
			auto front_sols2 = getMO_HLC().getSubspaceInfo(nei[i]).m_front_sol_in_subspace;
			//���������ӿռ�ǰ���е���̾���ֵ
			Real min_v = INT16_MAX;
			size_t idx1 = 0, idx2 = 0;
			for (size_t i = 0; i < front_sols1.size(); ++i) {
				for (size_t j = 0; j < front_sols2.size(); ++j) {
					auto dist = euclideanDistance(front_sols1[i]->variable().vect().begin(), front_sols1[i]->variable().vect().end(), front_sols2[j]->variable().vect().begin());
					min_v = min_v < dist ? min_v : dist;
					idx1 = i;
					idx2 = j;
				}
			}
			if (min_v < min_dist) {
				min_dist = min_v;
			}
			all_min_dist.push_back(min_v);
		}
		std::vector<Real> normal_dist;
		Real total = 0.;
		for (size_t i = 0; i < all_min_dist.size(); ++i) {
			auto temp = std::pow(min_dist / all_min_dist[i], 1.);
			normal_dist.push_back(temp);
			total += temp;
		}
		std::map<size_t, Real> out_put;
		for (size_t i = 0; i < all_min_dist.size(); ++i) {
			out_put.insert(std::make_pair<>(nei[i],normal_dist[i]/total));
		}

		return out_put;
	}

	std::map<size_t, Real> SPMOEA::spaceLinkProbability(size_t inx1,std::vector<size_t> neis) {
		//�����ӿռ��е�ǰ�ؽ�ķֲ�
		auto front_spaces = getFrontSpace();
		std::vector<std::vector<Real>> vars;
		if (std::find(front_spaces.begin(), front_spaces.end(), inx1) == front_spaces.end()) {
			auto &front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
			for (size_t i= 0; i < front_sols1.size(); ++i) {
				vars.emplace_back(front_sols1[i]->variable().vect());
			}
		}
		else {
			auto &front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
			for (size_t i = 0; i < front_sols1.size(); ++i) {
				vars.emplace_back(front_sols1[i]->variable().vect());
			}
		}
		
		std::vector<Real> all_min_dist;//�ӿռ�֮�����̾���
		Real min_dist = (Real)INT16_MAX;//�ӿռ�֮�����̾����е���С����
		for (size_t i = 0; i < neis.size(); ++i) {
			std::vector<std::vector<Real>> vars2;
			if (std::find(front_spaces.begin(), front_spaces.end(), neis[i]) == front_spaces.end()) {
				auto &front_sols2 = getMO_HLC().getSubspaceInfo(neis[i]).m_subspace_front_sol;
				for (size_t j = 0; j < front_sols2.size(); ++j) {
					vars2.emplace_back(front_sols2[j]->variable().vect());
				}
			}
			else {
				auto& front_sols2 = getMO_HLC().getSubspaceInfo(neis[i]).m_front_sol_in_subspace;
				for (size_t j = 0; j < front_sols2.size(); ++j) {
					vars2.emplace_back(front_sols2[j]->variable().vect());
				}
			}
			
			//���������ӿռ�ǰ���е���̾���ֵ
			Real min_v = INT16_MAX;
			size_t idx1 = 0, idx2 = 0;
			for (size_t i = 0; i < vars.size(); ++i) {
				for (size_t j = 0; j < vars2.size(); ++j) {
					auto dist = euclideanDistance(vars[i].begin(), vars[i].end(), vars2[j].begin());
					min_v = min_v < dist ? min_v : dist;
					idx1 = i;
					idx2 = j;
				}
			}
			if (min_v < min_dist) {
				min_dist = min_v;
			}
			all_min_dist.push_back(min_v);
		}
		std::vector<Real> normal_dist;
		Real total = 0.;
		for (size_t i = 0; i < all_min_dist.size(); ++i) {
			auto temp = std::pow(min_dist / all_min_dist[i], 1.);
			normal_dist.push_back(temp);
			total += temp;
		}
		std::map<size_t, Real> out_put;
		for (size_t i = 0; i < all_min_dist.size(); ++i) {
			out_put.insert(std::make_pair<>(neis[i], normal_dist[i] / total));
		}

		return out_put;
	}

	//���boolֵ����������ʣ�inx2��inx1�����ӿռ�ĸ���
	Real SPMOEA::ifContinuousSpace(size_t inx1, size_t inx2) {
		Real probability = 0.;
		//�����ӿռ��е�ǰ�ؽ�ķֲ�
		auto front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
		auto front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_in_subspace;
		//�����ӿռ��֮�����С����
		//���ӿռ�߽��ƽ������
		auto bound1 = getMO_HLC().subspaceTree().getBox(inx1);
		auto bound2 = getMO_HLC().subspaceTree().getBox(inx2);
		//���������ӿռ�ǰ���е���̾���ֵ
		Real min_v = INT16_MAX;
		size_t idx1 = 0,idx2=0;
		for (size_t i = 0; i < front_sols1.size(); ++i) {
			for (size_t j = 0; j < front_sols2.size(); ++j) {
				auto dist = euclideanDistance(front_sols1[i]->variable().vect().begin(), front_sols1[i]->variable().vect().end(), front_sols2[j]->variable().vect().begin());
				min_v = min_v < dist ? min_v : dist;
				idx1 = i;
				idx2 = j;
			}
		}
		//���ӿռ��ڵ�ƽ����С����
		Real mean_dist1 = 0.;
		Real mean_dist2 = 0;
		std::vector<Real> min_dist1;
		std::vector<Real> min_dist2;
		for (size_t i = 0; i < front_sols1.size(); ++i) {
			std::vector<Real> temp;
			auto& sol1 = front_sols1[i]->variable().vect();
			for (size_t j = 0; j < front_sols1.size(); ++j) {
				if (i==j) {
					temp.push_back((Real)INT16_MAX);
				}
				else{
					auto& sol2 = front_sols1[j]->variable().vect();
					temp.push_back(euclideanDistance(sol1.begin(),sol1.end(),sol2.begin()));
				}
			}
			auto dist = *std::min_element(temp.begin(), temp.end());
			mean_dist1 += dist;
			min_dist1.push_back(dist);
		}
		for (size_t i = 0; i < front_sols2.size(); ++i) {
			std::vector<Real> temp;
			auto& sol1 = front_sols2[i]->variable().vect();
			for (size_t j = 0; j < front_sols2.size(); ++j) {
				if (i == j) {
					temp.push_back((Real)INT16_MAX);
				}
				else {
					auto& sol2 = front_sols2[j]->variable().vect();
					temp.push_back(euclideanDistance(sol1.begin(), sol1.end(), sol2.begin()));
				}
			}
			auto dist = *std::min_element(temp.begin(), temp.end());
			mean_dist2 += dist;
			min_dist2.push_back(dist);
		}
		mean_dist1 /= front_sols1.size();
		mean_dist2 /= front_sols2.size();
		Real v = mean_dist1 / min_v;
		probability = 1. / (1 + std::exp(-1 * 5 * (v - 1. / 2)));
		size_t count1 = 0;
		for (size_t i = 0; i < min_dist1.size(); ++i) {
			if (min_dist1[i] > min_v) {
				count1++;
			}
		}
		size_t count2 = 0;
		for (size_t i = 0; i < min_dist2.size(); ++i) {
			if (min_dist2[i] > min_v) {
				count2++;
			}
		}
		Real ratio1= (Real)count1 / min_dist1.size();
		Real ratio2 = (Real)count2 / min_dist2.size();
		probability = ratio1 > ratio2 ? ratio1 : ratio2;

		return probability;
	}

	Real SPMOEA::ifSingleSegSpace(size_t inx1) {
		//�ӿռ��ڵ���������ֲ�������ж��ӿռ����Ƿ���ڶ��
		Real probability = 0.;
		//���ĵ��λ��
		auto& bound = getMO_HLC().subspaceTree().getBox(inx1);
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
		if (front_sols.size() < 3) {
			probability = 1.;
		}
		else {
			std::vector<Real> center_point;
			for (size_t i = 0; i < bound.size(); ++i) {
				Real s = 0.;
				for (size_t j = 0; j < front_sols.size(); ++j) {
					s += front_sols[j]->variable()[i];
				}
				center_point.push_back(s / front_sols.size());
			}
			//���������ĵ�����ľ���
			Real center_min_dist = (Real)INT16_MAX;
			for (size_t i = 0; i < front_sols.size(); ++i) {
				auto& sol = front_sols[i]->variable().vect();
				Real temp = euclideanDistance(center_point.begin(), center_point.end(), sol.begin());
				center_min_dist = center_min_dist < temp ? center_min_dist : temp;
			}
			//ǰ�ص���С����
			Real mean_min_dist = 0;
			std::vector<Real> min_dist;
			for (size_t i = 0; i < front_sols.size(); ++i) {
				std::vector<Real> temp;
				auto& sol1 = front_sols[i]->variable().vect();
				for (size_t j = 0; j < front_sols.size(); ++j) {
					if (i == j) {
						temp.push_back((Real)INT16_MAX);
					}
					else {
						auto& sol2 = front_sols[j]->variable().vect();
						temp.push_back(euclideanDistance(sol1.begin(), sol1.end(), sol2.begin()));
					}
				}
				auto dist = *std::min_element(temp.begin(), temp.end());
				mean_min_dist += dist;
				min_dist.push_back(dist);
			}
			mean_min_dist /= front_sols.size();
			Real v = mean_min_dist / center_min_dist;
			probability = 1. / (1 + std::exp(-1 * 2 * (v - 1)));
		}
		
		return probability;
	}

	std::pair<bool, std::pair<size_t, Real>> SPMOEA::ifMultiSegSpace(size_t inx1) {
		//����ӳ���ӿռ��ڵĵ㣬ͨ���жϼ����ռά�ȿ�ȵı������ж��ӿռ����Ƿ���ڶ��
		std::pair<bool, std::pair<size_t, Real>> out_put;
		auto& bound = getMO_HLC().subspaceTree().getBox(inx1);
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
		if (front_sols.size() < 10) {
			out_put.first = false;
			out_put.second = std::make_pair<>(0,0.);
		}
		else {
			//ӳ������ά�ж��Ƿ�ɷ�
			std::vector<Real> dim_interval;
			std::vector<Real> dim_span;
			std::vector<std::pair<size_t, size_t>> dim_inds;
			for (size_t i = 0; i < bound.size(); ++i) {
				std::vector<std::vector<Real>> map_value;
				for (size_t j = 0; j < front_sols.size(); ++j) {
					std::vector<Real> temp;
					temp.push_back(front_sols[j]->variable()[i]);
					temp.push_back(j);
					map_value.emplace_back(temp);
				}
				//�����ļ��ֵ
				//std::sort(map_value.begin(),map_value.end());
				sortrows(map_value,0);
				Real interval_v = 0.;
				size_t idx1, idx2;
				for (size_t j = 1; j < map_value.size(); ++j) {
					Real temp = map_value[j][0] - map_value[j - 1][0];
					if (temp > interval_v) {
						interval_v = temp;
						idx1 = (size_t)map_value[j-1][1];
						idx2 = (size_t)map_value[j][1];
					}
				}
				//Real ratio = interval_v / (map_value.back()[0] - map_value.front()[0]);
				Real ratio = interval_v;
				dim_inds.emplace_back(std::make_pair<>(idx1,idx2));
				dim_interval.emplace_back(ratio);
				dim_span.push_back(map_value.back()[0] - map_value.front()[0]);
			}
			Real max_ratio = *std::max_element(dim_interval.begin(),dim_interval.end());
			size_t dim_inx = std::distance(dim_interval.begin(), std::max_element(dim_interval.begin(), dim_interval.end()));
			bool flag = false;
			if (max_ratio / (dim_span[dim_inx] - dim_span[dim_inx]) > 0.4) {
				flag = true;
				std::pair<size_t, Real> split_pos;
				split_pos.first = dim_inx;
				Real pos1 = front_sols[dim_inds[dim_inx].first]->variable()[dim_inx];
				Real pos2 = front_sols[dim_inds[dim_inx].second]->variable()[dim_inx];
				split_pos.second = (pos1 + pos2) / 2 - bound[dim_inx].first;
				out_put.first = flag;
				out_put.second = split_pos;
			}
			else {
				out_put.first = false;
				out_put.second = std::make_pair<>(0, 0.);
			}
		}

		return out_put;
	}

	bool SPMOEA::ifLinearLinkSpace(size_t inx1, size_t inx2) {
		//�ϲ������ӿռ��ǰ�ؽ��ж����Ի��̶�
		//obtain front sols in inx subspace
		std::vector<std::vector<Real>> all_vars;
		auto& front_sols1 = getMO_HLC().getSubspaceInfo(inx1).m_front_sol_in_subspace;
		for (size_t i = 0; i < front_sols1.size(); ++i) {
			all_vars.emplace_back(front_sols1[i]->variable().vect());
		}
		auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_in_subspace;
		for (size_t i = 0; i < front_sols2.size(); ++i) {
			all_vars.emplace_back(front_sols2[i]->variable().vect());
		}
		if (all_vars.size() < 3 || front_sols1.size() == 1 || front_sols2.size() == 1) {
			return true;
		}
		else {
			std::vector<std::vector<Real>> temp_vars = all_vars;
			auto& bound1 = getMO_HLC().subspaceTree().getBox(inx1);
			auto& bound2 = getMO_HLC().subspaceTree().getBox(inx2);
			auto bound = bound1;
			for (size_t i = 0; i < bound.size(); ++i) {
				bound[i].first = bound1[i].first < bound2[i].first ? bound1[i].first : bound2[i].first;
				bound[i].second = bound1[i].second > bound2[i].second ? bound1[i].second : bound2[i].second;
			}
			dataNormalizeInBound(temp_vars, bound);
			//��һ�����ݺ���������������ϣ���������ֵ�ֽ���ϵ������
			Eigen::MatrixXd Apoints = Eigen::MatrixXd::Random(temp_vars.size(), bound.size());
			Eigen::MatrixXd bconstant = Eigen::MatrixXd::Ones(temp_vars.size(), 1);
			Eigen::MatrixXd Xk;
			for (size_t i = 0; i < temp_vars.size(); ++i) {
				for (size_t j = 0; j < bound.size(); ++j) {
					Apoints(i, j) = temp_vars[i][j];
				}
			}
			Xk = Apoints.colPivHouseholderQr().solve(bconstant);
			//�õõ���������Ϸ������
			std::vector<Real> cal_error(temp_vars.size(), 0.);
			for (size_t i = 0; i < temp_vars.size(); ++i) {
				auto temp = Apoints.row(i) * Xk;
				Real error = std::fabs(temp.value() - 1);
				if (error > 0.2) {//���ٷֱ�
					return false;
				}
			}
		}
		return true;
	}

	bool SPMOEA::ifLinearLinkSpace(size_t inx,std::vector<size_t> clustered_inx1, size_t inx2,Problem* pro, Random* rnd) {
		//�ۺ��Ѿ��������ӿռ��ǰ�ؽ���������ӿռ��ǰ�ؽ��ж����Ի��̶�
		// �����ж��ͳ�ƽ���ж�
		// �����ж����ýǶ�ƫ���ƽ���ж��������Է��������
		std::vector<std::vector<Real>> all_vars;
		size_t var_dim = CAST_CONOP(pro)->numberVariables();
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_front_sol_in_subspace;
		for (size_t j = 0; j < front_sols.size(); ++j) {
			all_vars.emplace_back(front_sols[j]->variable().vect());
		}
		auto& neighs = getMO_HLC().getSubspaceInfo(inx).m_sub_neighbors;
		std::vector<size_t> work_inx;
		for (size_t i = 0; i < clustered_inx1.size(); ++i) {
			//ѡ��inx�ӿռ�������ӿռ��ڵ�ǰ�ؽ����
			if (std::find(neighs.begin(), neighs.end(), clustered_inx1[i]) != neighs.end()) {
				work_inx.push_back(clustered_inx1[i]);
			}
		}
		if (work_inx.size() > 0) {
			size_t select_inx = (size_t)std::floor(work_inx.size()*rnd->uniform.next());
			auto& front_sols1 = getMO_HLC().getSubspaceInfo(work_inx[select_inx]).m_front_sol_in_subspace;
			for (size_t j = 0; j < front_sols1.size(); ++j) {
				all_vars.emplace_back(front_sols1[j]->variable().vect());
			}
		}
		//�����ӿռ�����
		auto& front_sols2 = getMO_HLC().getSubspaceInfo(inx2).m_front_sol_in_subspace;
		std::vector<std::vector<Real>> test_vars;
		for (size_t i = 0; i < front_sols2.size(); ++i) {
			test_vars.emplace_back(front_sols2[i]->variable().vect());
		}
		bool flag = true;

		if (all_vars.size()+test_vars.size() <= var_dim) {
			return true;
		}
		else {
			std::vector<std::vector<Real>> temp_vars = all_vars;
			auto bound = getMO_HLC().subspaceTree().getBox(inx2);
			//boundȡ������½�
			for (size_t i = 0; i < bound.size(); ++i) {
				bound[i].first = CAST_CONOP(pro)->boundary()[i].second;
				bound[i].second = CAST_CONOP(pro)->boundary()[i].first;
				for (size_t j = 0; j < temp_vars.size(); ++j) {
					if (bound[i].first > temp_vars[j][i]) {
						bound[i].first = temp_vars[j][i];
					}
					if (bound[i].second < temp_vars[j][i]) {
						bound[i].second = temp_vars[j][i];
					}
				}
				for (size_t j = 0; j < test_vars.size(); ++j) {
					if (bound[i].first > test_vars[j][i]) {
						bound[i].first = test_vars[j][i];
					}
					if (bound[i].second < test_vars[j][i]) {
						bound[i].second = test_vars[j][i];
					}
				}
			}
			dataNormalizeInBound(temp_vars, bound);
			//svd�����ֽ⣬��ȡ����������
			std::vector<Real> vector_data;
			for (size_t i = 0; i < temp_vars.size(); ++i) {
				for (size_t j = 0; j < var_dim; ++j) {
					vector_data.push_back(temp_vars[i][j]);
				}
			}
			Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix_data(vector_data.data(), temp_vars.size(), var_dim);
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix_data, Eigen::ComputeThinU | Eigen::ComputeFullV);
			Eigen::MatrixXd singular_values = svd.singularValues();
			Eigen::MatrixXd Vvectors = svd.matrixV();

			//�Ƚ�����ֵ�Ĵ�С��ȷ�����Գ̶Ⱥ�ά��
			std::vector<Real> svd_value;
			for (size_t i = 0; i < singular_values.size(); ++i) {
				svd_value.push_back(singular_values(i));
			}
			std::vector<std::vector<Real>> eigenvectors;
			for (size_t i = 0; i < svd_value.size(); ++i) {
				if (svd_value[i] > 10e-3) {
					std::vector<Real> row_data;
					for (size_t j = 0; j < var_dim; ++j) {
						row_data.push_back(Vvectors(i, j));
					}
					eigenvectors.emplace_back(row_data);//��������������ڹ�һ�����ԭ��
				}
			}
			//���Ե�����������ļн�
			std::vector<Real> center_p;
			for (size_t i = 0; i < var_dim; ++i) {
				Real s = 0.;
				for (size_t j = 0; j < temp_vars.size(); ++j) {
					s += temp_vars[j][i];
				}
				center_p.push_back(s / temp_vars.size());
				for (size_t j = 0; j < eigenvectors.size(); ++j) {
					eigenvectors[j][i] += (s / temp_vars.size());//��������ƽ��
				}
			}
			//ͳ�Ʋ��Ե�����������֮��ĽǶ�
			std::vector<Real> angles;
			std::vector<Real> vertical_dist;
			for (size_t j = 0; j < eigenvectors.size(); ++j) {
				for (size_t i = 0; i < test_vars.size(); ++i) {
					Real ang = vectorAngle(center_p, eigenvectors[j], test_vars[i]);
					//Real l = euclideanDistance(test_vars[i].begin(), test_vars[i].end(), corner.begin());
					//vertical_dist.push_back(l * std::sin(ang));
					angles.push_back(ang / OFEC_PI * 180);
				}
				size_t count = 0;
				for (size_t i = 0; i < angles.size(); ++i) {
					if (angles[i] > 20 || angles[i] < 160) {//�н�
						count++;
					}
				}
				if (count > angles.size() / 2) {
					flag = false;
				}
			}

			////�Ƚ��������жϣ��ٽ��г�ƽ���ж�
		 //   //ͳ�������нǣ�������Զ�����㣬����ʹ��������������
			//Real max_dist = 0;
			//std::pair<size_t, size_t> max_inx;
			//for (size_t i = 0; i < all_vars.size(); ++i) {
			//	for (size_t j = i + 1; j < all_vars.size(); ++j) {
			//		Real temp_dist = euclideanDistance(all_vars[i].begin(), all_vars[i].end(), all_vars[j].begin());
			//		if (temp_dist > max_dist) {
			//			max_dist = temp_dist;
			//			max_inx.first = i;
			//			max_inx.second = j;
			//		}
			//	}
			//}
			////����ǵ�λ��
			//auto corner = all_vars[max_inx.first];
			///*for (size_t i = 0; i < var_dim; ++i) {
			//	corner[i] = 2 * corner[i] - all_vars[max_inx.second][i];
			//}*/
			////ͳ�Ʋ��Ե�������֮��ľ���
			//std::vector<Real> angles;
			//std::vector<Real> vertical_dist;
			//for (size_t i = 0; i < test_vars.size(); ++i) {
			//	Real ang = vectorAngle(corner, all_vars[max_inx.second], test_vars[i]);
			//	Real l = euclideanDistance(test_vars[i].begin(), test_vars[i].end(), corner.begin());
			//	vertical_dist.push_back(l * std::sin(ang));
			//	angles.push_back(ang / OFEC_PI * 180);
			//}
			//size_t count = 0;
			//for (size_t i = 0; i < vertical_dist.size(); ++i) {
			//	if (vertical_dist[i] > max_dist*0.5) {//���ٷֱ�
			//		count++;
			//	}
			//}
			//if (count > angles.size() / 2) {
			//	flag = false;
			//}
			//if (!flag) {
			//	
			//	

			//	////��һ�����ݺ���������������ϣ���������ֵ�ֽ���ϵ������
			//	//Eigen::MatrixXd Apoints = Eigen::MatrixXd::Random(temp_vars.size(), bound.size());
			//	//Eigen::MatrixXd bconstant = Eigen::MatrixXd::Ones(temp_vars.size(), 1);
			//	//Eigen::MatrixXd Xk1, Xk2;
			//	//for (size_t i = 0; i < temp_vars.size(); ++i) {
			//	//	for (size_t j = 0; j < bound.size(); ++j) {
			//	//		Apoints(i, j) = temp_vars[i][j];
			//	//	}
			//	//}
			//	//Xk1 = Apoints.colPivHouseholderQr().solve(bconstant);
			//	////Xk2 = Apoints.jacobiSvd(Eigen::ComputeThinV | Eigen::ComputeThinU).solve(bconstant);

			//	//std::vector<Real> cal_error(temp_vars.size(), 0.);
			//	//Real mean_error = 0.;
			//	//for (size_t i = 0; i < temp_vars.size(); ++i) {
			//	//	auto temp = Apoints.row(i) * Xk1;
			//	//	Real error = std::fabs(temp.value() - 1);
			//	//	cal_error[i] = error;
			//	//	mean_error += error;
			//	//}
			//	//mean_error /= temp_vars.size();

			//	////// 1���õõ���������Ϸ������
			//	////size_t count = 0;
			//	////for (size_t i = 0; i < cal_error.size(); ++i) {
			//	////	if (cal_error[i] > 0.2) {//���ٷֱ�
			//	////		count++;
			//	////	}
			//	////}
			//	////if (count > cal_error.size() / 5) {
			//	////	return false;
			//	////}

			//	//// 2��ʹ�ô������ӿռ��������
			//	//std::vector<Real> test_error(test_vars.size(), 0.);
			//	//Eigen::MatrixXd test_Apoints = Eigen::MatrixXd::Random(test_vars.size(), bound.size());
			//	//for (size_t i = 0; i < test_vars.size(); ++i) {
			//	//	for (size_t j = 0; j < bound.size(); ++j) {
			//	//		test_Apoints(i, j) = test_vars[i][j];
			//	//	}
			//	//}
			//	//for (size_t i = 0; i < test_vars.size(); ++i) {
			//	//	auto temp = test_Apoints.row(i) * Xk1;
			//	//	Real error = std::fabs(temp.value() - 1);
			//	//	test_error[i] = error;
			//	//}
			//	//size_t count = 0;
			//	//for (size_t i = 0; i < test_error.size(); ++i) {
			//	//	if (test_error[i] > mean_error) {//���ٷֱ�
			//	//		count++;
			//	//	}
			//	//}
			//	//if (count > test_error.size() / 2) {
			//	//	return false;
			//	//}
			//	////// 3������ֵ�ֽ⣬��������ֵ�ж����Գ̶�
			// //   //std::vector<Real> vector_data;
			// //   //for (size_t i = 0; i < temp_vars.size(); ++i) {
			// //   //	for (size_t j = 0; j < var_dim; ++j) {
			// //   //		vector_data.push_back(temp_vars[i][j]);
			// //   //	}
			// //   //}
			// //   //Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> matrix_data(vector_data.data(), temp_vars.size(), var_dim);
			// //   //Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix_data, Eigen::ComputeThinV | Eigen::ComputeThinU);
			// //   //Eigen::MatrixXd singular_values = svd.singularValues();
			// //   ////�Ƚ�����ֵ�Ĵ�С��ȷ�����Գ̶Ⱥ�ά��
			// //   //std::vector<Real> qr_value;
			// //   //for (size_t i = 0; i < singular_values.size(); ++i) {
			// //   //	qr_value.push_back(singular_values(i));
			// //   //}
			// //   ////ͨ�����һ������ֵ�Ĵ�С��ǰһ������ֵ�Ĵ�С�ж�
			// //   ////������Ч������ֵ������С��ά��
			// //   //for (size_t i = qr_value.size() - 1; i > 1; --i) {
			// //   //	Real qr1 = qr_value[i];
			// //   //	Real qr2 = qr_value[i - 1];
			// //   //	if (qr1 / qr2 < 0.1 || qr1 < 10e-5) {
			// //   //		return true;
			// //   //	}
			// //   //}
			//}
		}
		return flag;
	}

	//���������ռ�ǰ���ӿռ�ľ��ࣺ�������+������Էֲ�
	std::vector<std::vector<size_t>> SPMOEA::linearClusterFrontSpace(std::vector<size_t>& frontspace,Problem* pro, Random* rnd) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//���
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			////����ѡ���ӿռ��ڸ��������ٵ��ӿռ俪ʼ��
			//std::vector<size_t> num_front_sols;
			//for (size_t i = 0; i < select_flag.size(); ++i) {
			//	if (select_flag[i] == 1) {
			//		num_front_sols.push_back(INT16_MAX);
			//	}
			//	else {
			//		num_front_sols.push_back(getMO_HLC().getSubspaceInfo(frontspace[i]).m_front_sol_in_subspace.size());
			//	}
			//}
			//auto inx = std::distance(num_front_sols.begin(), std::min_element(num_front_sols.begin(), num_front_sols.end()));
			//head_cluster.push_back(frontspace[inx]);//�ȼ���һ���ӿռ�
			//select_flag[inx] = 1;
			
			//�����ӿռ�žۣ��������ӿռ䵥������
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					head_cluster.push_back(frontspace[i]);//�ȼ���һ���ӿռ�
					select_flag[i] = 1;
					break;
				}
			}
			for (size_t i = 0; i < select_flag.size(); ++i) {
				count += select_flag[i];
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
			if (getMO_HLC().getSubspaceInfo(head_cluster[0]).m_linear_flag==0) {
				clustered.emplace_back(head_cluster);
			}
			else {
				auto temp_cluster = head_cluster;
				while (count < select_flag.size()) {
					std::vector<size_t> temp;
					for (size_t j = 0; j < temp_cluster.size(); ++j) {
						size_t inx = temp_cluster[j];
						std::list<size_t> neighbors;
						getMO_HLC().subspaceTree().findNeighbor(inx, neighbors);
						for (size_t k = 0; k < frontspace.size(); ++k) {
							if (select_flag[k] == 0) {
								if (std::find(neighbors.begin(), neighbors.end(), frontspace[k]) != neighbors.end()) {
									//�ж��Ƿ����������ӿռ����Թ�ϵ
									if (ifLinearLinkSpace(inx, head_cluster, frontspace[k], pro, rnd)) {
										head_cluster.push_back(frontspace[k]);
										temp.push_back(frontspace[k]);
										select_flag[k] = 1;
										count++;
									}
								}
							}
						}
					}
					if (temp.empty()) {//û���µ��������
						break;
					}
					else {
						temp_cluster = temp;
					}
				}
				clustered.emplace_back(head_cluster);
			}
		}
		return clustered;
	}

	//���������ӿռ�ľ��ࣺ��������࣬���ǿ���ǰ�ؽ�ķֲ���
	std::vector<std::vector<size_t>> SPMOEA::clusterRegionFrontSpace(std::vector<size_t>& frontspace) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//���
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					begin_space = i;
					head_cluster.push_back(frontspace[begin_space]);//�ȼ���һ���ӿռ�
					select_flag[i] = 1;
					break;
				}
			}
			for (size_t i = 0; i < select_flag.size(); ++i) {
				count += select_flag[i];
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
			auto temp_cluster = head_cluster;
			while (count < select_flag.size()) {
				std::vector<size_t> temp;
				for (size_t j = 0; j < temp_cluster.size(); ++j) {
					size_t inx = temp_cluster[j];
					std::list<size_t> neighbors;
					getMO_HLC().subspaceTree().findNeighbor(inx, neighbors);
					for (size_t k = 0; k < frontspace.size(); ++k) {
						if (select_flag[k] == 0) {
							if (std::find(neighbors.begin(), neighbors.end(), frontspace[k]) != neighbors.end()) {
								head_cluster.push_back(frontspace[k]);
								temp.push_back(frontspace[k]);
								select_flag[k] = 1;
								count++;
							}
						}
					}
				}
				if (temp.empty()) {//û���µ��������
					//clustered.emplace_back(head_cluster);
					break;
				}
				else {
					temp_cluster = temp;
				}
			}
			clustered.emplace_back(head_cluster);
			/*if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}*/
		}
		//getMO_HLC().setFrontClusters(clustered);
		return clustered;
	}

	void SPMOEA::updateVarSpaceInfo(Population<Solution<>>& pop, Problem *pro, Random *rnd) {
		m_mo_hlc->updateSubspaceInfo(pop, pro, rnd);
		//�����ӿռ���������ֵ��ȡ�ӿռ�Ĵ�����廹��ǰ�ظ��壿
		std::vector<int> pre_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			pre_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
		std::vector<std::shared_ptr<Solution<>>> temp_pop;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			for (size_t j = 0; j < m_mo_hlc->getSubspaceInfo(i).m_represent_sol.size(); ++j) {
				temp_pop.emplace_back(m_mo_hlc->getSubspaceInfo(i).m_represent_sol[j]);
			}
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			objs.emplace_back(&temp_pop[i]->objective());
		}
		std::vector<int> rank;
		int num_layer=ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		getMO_HLC().setLayerNum(num_layer);
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			temp_pop[i]->setFitness(rank[i]);
		}
		//�����ӿռ��������ֵ
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			m_mo_hlc->getSubspaceInfo(i).m_best_rank = INT16_MAX;
		}
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			auto temp_sol = temp_pop[i]->variable().vect();
			auto inx = m_mo_hlc->subspaceTree().getRegionIdx(temp_sol);
			if (temp_pop[i]->fitness() < m_mo_hlc->getSubspaceInfo(inx).m_best_rank) {
				m_mo_hlc->getSubspaceInfo(inx).m_best_rank = temp_pop[i]->fitness();
			}
		}
		//����û�и�����ӿռ䣬��������õ�����ֵ
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank == INT16_MAX) {
				auto neigh_space = m_mo_hlc->getSubspaceInfo(i).m_sub_neighbors;
				std::vector<int> neigh_ranks;
				auto iter = neigh_space.begin();
				for (size_t j = 0; j < neigh_space.size(); ++j) {
					if (m_mo_hlc->getSubspaceInfo(*iter).m_best_rank<INT16_MAX) {
						neigh_ranks.push_back(m_mo_hlc->getSubspaceInfo(*iter).m_best_rank);
					}
					++iter;
				}
				if (neigh_ranks.empty()) {
					m_mo_hlc->getSubspaceInfo(i).m_best_rank = num_layer;
				}
				else if (neigh_ranks.size() == 1) {
					m_mo_hlc->getSubspaceInfo(i).m_best_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end()) + 1;
				}
				else {
					int min_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end());
					int max_rank = *std::max_element(neigh_ranks.begin(), neigh_ranks.end());
					m_mo_hlc->getSubspaceInfo(i).m_best_rank = std::ceil(min_rank+max_rank);
				}
			}
		}
		std::vector<int> past_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			past_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
	}

	void SPMOEA::updateVarSpaceRank(Problem* pro, Random* rnd) {
		//�����ӿռ���������ֵ��ȡ�ӿռ�Ĵ�����壬�ų�ǰ�ؿռ�
		std::vector<int> pre_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			pre_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
		auto front_space = getFrontSpace();
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (std::find(front_space.begin(),front_space.end(),i)==front_space.end()) {
				m_mo_hlc->getSubspaceInfo(i).m_best_rank = INT16_MAX;
			}
		}
		std::vector<std::shared_ptr<Solution<>>> temp_pop;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank != 0) {
				for (size_t j = 0; j < m_mo_hlc->getSubspaceInfo(i).m_represent_sol.size(); ++j) {
					temp_pop.emplace_back(m_mo_hlc->getSubspaceInfo(i).m_represent_sol[j]);
				}
			}
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			objs.emplace_back(&temp_pop[i]->objective());
		}
		std::vector<int> rank;
		int num_layer = ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		getMO_HLC().setLayerNum(num_layer);
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			temp_pop[i]->setFitness(rank[i]);
		}
		//�����ӿռ��������ֵ
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			auto temp_sol = temp_pop[i]->variable().vect();
			auto inx = m_mo_hlc->subspaceTree().getRegionIdx(temp_sol);
			if ((temp_pop[i]->fitness()+1) < m_mo_hlc->getSubspaceInfo(inx).m_best_rank) {
				m_mo_hlc->getSubspaceInfo(inx).m_best_rank = temp_pop[i]->fitness()+1;
			}
		}
		//����û�и�����ӿռ䣬��������õ�����ֵ
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank == INT16_MAX) {
				/*auto neigh_space = m_mo_hlc->getSubspaceInfo(i).m_sub_neighbors;
				std::vector<int> neigh_ranks;
				auto iter = neigh_space.begin();
				for (size_t j = 0; j < neigh_space.size(); ++j) {
					if (m_mo_hlc->getSubspaceInfo(*iter).m_best_rank < INT16_MAX) {
						neigh_ranks.push_back(m_mo_hlc->getSubspaceInfo(*iter).m_best_rank);
					}
					++iter;
				}
				if (neigh_ranks.empty()) {
					m_mo_hlc->getSubspaceInfo(i).m_best_rank = num_layer;
				}
				else if (neigh_ranks.size() == 1) {
					m_mo_hlc->getSubspaceInfo(i).m_best_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end()) + 1;
				}
				else {
					int min_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end());
					int max_rank = *std::max_element(neigh_ranks.begin(), neigh_ranks.end());
					m_mo_hlc->getSubspaceInfo(i).m_best_rank = std::ceil(min_rank + max_rank);
				}*/
				m_mo_hlc->getSubspaceInfo(i).m_best_rank = num_layer+1;
			}
		}
		std::vector<int> past_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			past_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
	}

	void SPMOEA::updateVarSpaceDominanceRank(Problem* pro, Random* rnd) {
		//�����ӿռ���������ֵ��ȡ�ӿռ�Ĵ�����壬�ų�ǰ�ؿռ�
		std::vector<int> pre_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			pre_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank != 0) {
				m_mo_hlc->getSubspaceInfo(i).m_best_rank = INT16_MAX;
			}
		}
		std::vector<std::shared_ptr<Solution<>>> temp_pop;
		std::vector<size_t> all_count;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank != 0) {
				size_t count = 1;
				auto obj1 = m_mo_hlc->getSubspaceInfo(i).m_represent_sol[0]->objective();
				for (size_t j = 0; j < m_mo_hlc->numSubspace(); ++j) {
					if (j != i) {
						auto obj2= m_mo_hlc->getSubspaceInfo(j).m_represent_sol[0]->objective();
						auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
						if (ship == Dominance::kDominated) {
							count++;
						}
					}
				}
				m_mo_hlc->getSubspaceInfo(i).m_best_rank = count;
				all_count.push_back(count);
			}
		}
		size_t max_count = *std::max_element(all_count.begin(), all_count.end());
		//����û�и�����ӿռ䣬��������õ�����ֵ
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			if (m_mo_hlc->getSubspaceInfo(i).m_best_rank == INT16_MAX) {
				m_mo_hlc->getSubspaceInfo(i).m_best_rank = max_count + 1;
			}
		}
		std::vector<int> past_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			past_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
	}

	void SPMOEA::updateSubspaceFrontSol(Population<Solution<>>& pop, Problem *pro, Random *rnd) {
		std::map<size_t, std::vector<size_t>> pop2space;
		for (size_t i = 0; i < pop.size(); ++i) {
			auto var = pop[i].variable().vect();
			auto idx = getMO_HLC().subspaceTree().getRegionIdx(var);
			if (pop2space[idx].empty()) {
				std::vector<size_t> temp_inx;
				pop2space.insert(std::make_pair(idx, temp_inx));
			}
			pop2space[idx].emplace_back(i);
		}
		//�����ӿռ�ǰ�ؽ�ʹ����
		for (auto& sub : pop2space) {
			size_t idx = sub.first;
			Population<Solution<>> temp_pop;
			for (size_t i = 0; i < sub.second.size(); ++i) {
				temp_pop.append(pop[sub.second[i]]);
			}
			std::vector<std::vector<Real>*> objs;
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				objs.emplace_back(&temp_pop[i].objective());
			}
			std::vector<int> rank;
			ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				temp_pop[i].setFitness(rank[i]);
			}
			Population<Solution<>> front_pop;
			std::vector<size_t> add_front_inx;
			std::vector<size_t> add_behind_inx;
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				if (temp_pop[i].fitness() == 0) {
					front_pop.append(temp_pop[i]);
					add_front_inx.push_back(i);
				}
				else {
					add_behind_inx.push_back(i);
				}
				getMO_HLC().getSubspaceInfo(idx).m_sub_freq++;
			}
			//�ӿռ�ǰ�ؼ��������£��ӿռ��������
			if (getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol.empty()) {
				for (size_t i = 0; i < front_pop.size(); ++i) {
					//getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol.emplace_back(new Solution<>(front_pop[i]));
					getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx.emplace_back(getMO_HLC().getSubspaceInfo(idx).m_history_inds.size());
					getMO_HLC().getSubspaceInfo(idx).m_history_inds.emplace_back(new Solution<>(front_pop[i]));
				}
			}
			else {
				auto temp_front_sols = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol;
				auto temp_front_inx =  getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx;
				getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol.clear();
				getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx.clear();
				std::vector<size_t> add_his_inx;
				/*std::vector<size_t> final_add_front_inx;
				std::vector<size_t> temp_add_inx;
				for (size_t i = 0; i < temp_front_sols.size(); ++i) {
					temp_add_inx.push_back(i);
				}
				std::vector<size_t> final_add;
				
				for (size_t i = 0; i < add_front_inx.size(); ++i) {
					bool flag = false;
					auto& obj2 = temp_pop[add_front_inx[i]].objective();
					for (size_t j = 0; j < temp_front_sols.size(); ++j) {
						auto& obj1 = temp_front_sols[j]->objective();
						auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
						if (ship == Dominance::kDominant) {
							flag = true;
							break;
						}
					}
					if (!flag) {
						final_add = temp_add_inx;
						temp_add_inx.clear();
						final_add_front_inx.push_back(add_front_inx[i]);
						for (size_t j = 0; j < final_add.size(); ++j) {
							if (!temp_pop[add_front_inx[i]].dominate(temp_front_sols[final_add[j]]->solut(), pro)) {
								temp_add_inx.push_back(final_add[j]);
							}
						}
					}
					else {
						add_his_inx.push_back(add_front_inx[i]);
					}
				}*/

				//�ȿ�ÿ�����Ƿ񱻵�ǰǰ��֧��
				std::vector<size_t> temp_his_sols(temp_front_sols.size(), 1);
				std::vector<size_t> temp_pop_sols(add_front_inx.size(), 0);
				for (size_t j = 0; j < add_front_inx.size(); ++j) {
					for (size_t i = 0; i < temp_front_sols.size(); ++i) {
						if (temp_his_sols[i] == 1) {
							Dominance dominanceship = temp_pop[add_front_inx[j]].compare(*temp_front_sols[i], pro->optimizeMode());
							if (dominanceship == Dominance::kNonDominated) {
								if (i == temp_his_sols.size() - 1) {
									temp_pop_sols[j] = 1;
								}
							}
							else if (dominanceship == Dominance::kDominant) {
								temp_his_sols[i] = 0;
								if (i == temp_his_sols.size() - 1) {
									temp_pop_sols[j] = 1;
								}
							}
							else if (dominanceship==Dominance::kEqual) {
								//temp_his_sols[i] = 0;
								if (i == temp_his_sols.size() - 1) {
									temp_pop_sols[j] = 1;
								}
							}
							else {
								add_his_inx.push_back(j);
								break;
							}
						}
						else {
							if (i == temp_his_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
					}
				}

				//�����ӿռ�ǰ�ؽ�����
				size_t count = 0;
				for (size_t i = 0; i < temp_pop_sols.size(); ++i) {
					if (temp_pop_sols[i] == 1) {
						count++;
						getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx.push_back(getMO_HLC().getSubspaceInfo(idx).m_history_inds.size());
						getMO_HLC().getSubspaceInfo(idx).m_history_inds.emplace_back(new Solution<>(temp_pop[add_front_inx[i]]));
					}
				}
				for (size_t i = 0; i < temp_his_sols.size(); ++i) {
					if (temp_his_sols[i] == 1) {
						getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx.push_back(temp_front_inx[i]);
					}
					else {
						getMO_HLC().getSubspaceInfo(idx).m_history_inds[temp_front_inx[i]]->setType(-1);
					}
				}
				//���뱻��̭�Ľ�����ʷ��
				for (size_t i = 0; i < add_his_inx.size(); ++i) {
					count++;
					getMO_HLC().getSubspaceInfo(idx).m_history_inds.emplace_back(new Solution<>(temp_pop[add_front_inx[add_his_inx[i]]]));
				}
				if (count != add_front_inx.size()) {
					size_t a = 1;
				}
			}
			for (size_t i = 0; i < add_behind_inx.size(); ++i) {
				getMO_HLC().getSubspaceInfo(idx).m_history_inds.emplace_back(new Solution<>(temp_pop[add_behind_inx[i]]));
			}
			//�����ӿռ�ǰ�ؽ�
			for (size_t i = 0; i < getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx.size(); ++i) {
				getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol.emplace_back(getMO_HLC().getSubspaceInfo(idx).m_history_inds[getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx[i]]);
				getMO_HLC().getSubspaceInfo(idx).m_history_inds[getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx[i]]->setType(0);
			}

			///*�����ӿռ�ǰ���Ƿ���ȷ*/
			//auto front_inx = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx;
			//for (size_t i = 0; i < getMO_HLC().getSubspaceInfo(idx).m_history_inds.size(); ++i) {
			//	if (std::find(front_inx.begin(), front_inx.end(), i) == front_inx.end()) {
			//		auto p1 = getMO_HLC().getSubspaceInfo(idx).m_history_inds[i]->objective();
			//		for (size_t j = 0; j < front_inx.size(); ++j) {
			//			auto p2 = getMO_HLC().getSubspaceInfo(idx).m_history_inds[front_inx[j]]->objective();
			//			Dominance ship = objectiveCompare(p2, p1, pro->optimizeMode());
			//			if (ship == Dominance::kDominant||ship==Dominance::kEqual) {
			//				break;
			//			}
			//			if (j == front_inx.size() - 1) {
			//				size_t a = 1;
			//			}
			//		}
			//	}
			//}

			//�����ӿռ����⣬���ѡ����⣬�������������ӿռ��С���
			size_t num_var = CAST_CONOP(pro)->numberVariables();
			//select_num = (num_var + 1) * std::floor(std::pow(subspaceTree().getBoxVolume(idx) / min_volume, 1. / 3));
			size_t total_num = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol.size();
			size_t select_num = total_num;
			select_num = (size_t)std::ceil(getMO_HLC().getSubspaceInfo(idx).m_represent_num * std::pow(num_var, 2. / 4));
			//select_num = 1;
			getMO_HLC().getSubspaceInfo(idx).m_represent_sol.clear();
			getMO_HLC().getSubspaceInfo(idx).m_subspace_represent_inx.clear();
			if (total_num <= select_num) {
				getMO_HLC().getSubspaceInfo(idx).m_subspace_represent_inx = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx;
			}
			else {
				bool random_select = true;
				if (random_select) {
					for (size_t j = 0; j < select_num; ++j) {
						size_t inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol.size() * rnd->uniform.next());
						getMO_HLC().getSubspaceInfo(idx).m_subspace_represent_inx.push_back(getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx[inx]);
					}
				}
				else {
					std::vector<size_t> selected_inx;
					//�ȼ���ĳһĿ���µ���ֵ
					Real min_v = INT16_MAX;
					size_t add_inx = 0;
					for (size_t i = 0; i < total_num; ++i) {
						if (getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol[i]->objective()[0] < min_v) {
							min_v = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol[i]->objective()[0];
							add_inx = i;
						}
					}
					selected_inx.push_back(add_inx);
					//�ȵõ�δѡ������
					std::vector<size_t> no_select_inx;//δѡ������
					for (size_t i = 0; i < total_num; ++i) {
						if (std::find(selected_inx.begin(), selected_inx.end(), i) == selected_inx.end()) {
							no_select_inx.push_back(i);
						}
					}
					while (selected_inx.size() < select_num) {//ѡ���뵱ǰ��ѡ����Զ�ĵ�
						//����δѡ������ѡ�����С��������ֵ
						std::vector<Real> min_dist;
						for (size_t i = 0; i < no_select_inx.size(); ++i) {
							std::vector<Real> temp_dist;
							for (size_t j = 0; j < selected_inx.size(); ++j) {
								auto p1 = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol[no_select_inx[i]]->objective();
								auto p2 = getMO_HLC().getSubspaceInfo(idx).m_subspace_front_sol[selected_inx[j]]->objective();
								temp_dist.push_back(euclideanDistance(p1.begin(), p1.end(), p2.begin()));
							}
							Real min_d = *std::min_element(temp_dist.begin(), temp_dist.end());
							min_dist.push_back(min_d);
						}
						auto s_inx = std::distance(min_dist.begin(), std::max_element(min_dist.begin(), min_dist.end()));
						selected_inx.push_back(no_select_inx[s_inx]);
						no_select_inx.erase(no_select_inx.begin() + s_inx);
					}
					for (size_t i = 0; i < selected_inx.size(); ++i) {
						getMO_HLC().getSubspaceInfo(idx).m_subspace_represent_inx.push_back(getMO_HLC().getSubspaceInfo(idx).m_subspace_front_inx[selected_inx[i]]);
					}
				}
			}
			//�����ӿռ�����
			for (size_t i = 0; i < getMO_HLC().getSubspaceInfo(idx).m_subspace_represent_inx.size(); ++i) {
				getMO_HLC().getSubspaceInfo(idx).m_represent_sol.emplace_back(getMO_HLC().getSubspaceInfo(idx).m_history_inds[getMO_HLC().getSubspaceInfo(idx).m_subspace_represent_inx[i]]);
			}
			
			getMO_HLC().getSubspaceInfo(idx).add_flag = "yes";
		}
	}

	void SPMOEA::updateHisVarSpaceInfo(std::vector<std::shared_ptr<Solution<>>>& pop, Problem *pro, Random *rnd) {
		m_mo_hlc->updateSubSpaceInfo(pop, pro, rnd);
		//�����ӿռ���������ֵ��ȡ�ӿռ�Ĵ�����廹��ǰ�ظ��壿
		std::vector<int> pre_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			pre_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
		std::vector<std::shared_ptr<Solution<>>> temp_pop;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			for (size_t j = 0; j < m_mo_hlc->getSubspaceInfo(i).m_represent_sol.size(); ++j) {
				temp_pop.emplace_back(m_mo_hlc->getSubspaceInfo(i).m_represent_sol[j]);
			}
		}
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			objs.emplace_back(&temp_pop[i]->objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			temp_pop[i]->setFitness(rank[i]);
		}
		//�����ӿռ��������ֵ
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			m_mo_hlc->getSubspaceInfo(i).m_best_rank = INT16_MAX;
		}
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			auto temp_sol = temp_pop[i]->variable().vect();
			auto inx = m_mo_hlc->subspaceTree().getRegionIdx(temp_sol);
			if (temp_pop[i]->fitness() < m_mo_hlc->getSubspaceInfo(inx).m_best_rank) {
				m_mo_hlc->getSubspaceInfo(inx).m_best_rank = temp_pop[i]->fitness();
			}
		}
		std::vector<int> past_rank;
		for (size_t i = 0; i < m_mo_hlc->numSubspace(); ++i) {
			past_rank.push_back(m_mo_hlc->getSubspaceInfo(i).m_best_rank);
		}
	}

	std::vector<size_t> SPMOEA::getFrontSubspace() {
		std::vector<size_t> front_subspace;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank == 0) {
				front_subspace.push_back(i);
			}
		}
		return front_subspace;
	}

	std::vector<size_t> SPMOEA::getBehindSubspace() {
		std::vector<size_t> behind_subspace;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			if (getMO_HLC().getSubspaceInfo(i).m_best_rank > 0) {
				behind_subspace.push_back(i);
			}
		}
		return behind_subspace;
	}

	/*void SPMOEA::updateVarSpaceInfo(Solution<>& sol) {
		m_mo_hlc->updateSubspaceInfo(sol);
	}*/

	void SPMOEA::updateObjRange(Population<Solution<>> &pop, Problem *pro) {
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		////���µ�ǰ��ȺĿ��ֵ��Χ
		//if (m_pop_obj_range.empty()) {
		//	for (int i = 0; i < num_obj; ++i) {
		//		m_pop_obj_range.push_back(std::make_pair<Real, Real>(INT16_MAX, -1 * INT16_MAX));
		//	}
		//}
		//for (int j = 0; j < num_obj; ++j) {
		//	for (int k = 0; k < pop.size(); ++k) {
		//		if (m_pop_obj_range[j].first > pop[k].objective()[j]) {
		//			m_pop_obj_range[j].first = pop[k].objective()[j];
		//		}
		//		if (m_pop_obj_range[j].second < pop[k].objective()[j]) {
		//			m_pop_obj_range[j].second = pop[k].objective()[j];
		//		}
		//	}
		//}
		//������ʷǰ��Ŀ��ֵ��Χ
		if (m_front_obj_range.empty()) {
			for (int i = 0; i < num_obj; ++i) {
				m_front_obj_range.push_back(std::make_pair<Real, Real>(INT16_MAX, INT16_MAX));
			}
		}
		for (int j = 0; j < pop.size(); ++j) {
			for (int k = 0; k < num_obj; ++k) {
				if (m_front_obj_range[k].first > pop[j].objective()[k]) {
					for (int p = 0; p < num_obj; ++p) {
						if (p == k) {
							m_front_obj_range[p].first = pop[j].objective()[p];
						}
						else {
							m_front_obj_range[p].second = pop[j].objective()[p];
						}
					}
				}
			}
		}
		//������ʷ��Ŀ��ֵ��Χ
		if (m_his_obj_range.empty()) {
			for (int i = 0; i < num_obj; ++i) {
				m_his_obj_range.push_back(std::make_pair<Real, Real>(INT16_MAX, -1 * INT16_MAX));
			}
		}
		for (int j = 0; j < num_obj; ++j) {
			for (int k = 0; k < pop.size(); ++k) {
				if (m_his_obj_range[j].first > pop[k].objective()[j]) {
					m_his_obj_range[j].first = pop[k].objective()[j];
				}
				if (m_his_obj_range[j].second < pop[k].objective()[j]) {
					m_his_obj_range[j].second = pop[k].objective()[j];
				}
			}
		}
		//CAST_CONOP(pro)->setObjRange(m_his_obj_range);
	}

	void SPMOEA::updatePopDistribute(size_t inx, Problem *pro) {
		getPop()[inx].updatePopDistribute(pro);
	}

	void SPMOEA_pop::updatePopDistribute(Problem *pro) {
		//ʹ�ø���������ĵ�������ķֲ�������ʹ�������ռ��Ŀ��ռ�����ĵ�λ��
		/*auto bound = getSearchRange();
		std::vector<Real> ref_point;
		for (size_t i = 0; i < bound.size(); ++i) {
			ref_point.push_back(bound[i].first);
		}*/
		std::vector<Real> var_center;//����ڲο��������
		for (size_t i = 0; i < CAST_CONOP(pro)->numberVariables(); ++i) {
			Real dim_sum = 0;
			for (size_t j = 0; j < this->size(); ++j) {
				dim_sum += this->m_individuals[j]->variable().vect()[i];
			}
			var_center.push_back(dim_sum / this->size());
		}
		std::vector<Real> obj_center;//����ȺĿ��ռ������0���������
		for (size_t i = 0; i < CAST_CONOP(pro)->numberObjectives(); ++i) {
			Real dim_sum = 0;
			for (size_t j = 0; j < this->size(); ++j) {
				dim_sum += this->m_individuals[j]->objective()[i];
			}
			obj_center.push_back(dim_sum / this->size());
		}
		std::vector<std::vector<Real>> pop_distribute;
		pop_distribute.emplace_back(var_center);
		pop_distribute.emplace_back(obj_center);
		getPopDistribute().push_back(pop_distribute);
	}

	void SPMOEA_pop::updatePopHisInfo(Population<Solution<>> &new_pop,Problem *pro) {
		//ʹ���Ӵ�������Ⱥ����ʷǰ�أ�ÿ����ǰ��
		NDSort(new_pop,pro);
		Population<Solution<>> front_pop;
		for (size_t i = 0; i < new_pop.size(); ++i) {
			if (new_pop[i].fitness() == 0) {
				front_pop.append(new_pop[i]);
			}
		}
		//�ӿռ�ǰ�ظ��£��ӿռ��������
		if (m_pop_history_front_sols.empty()) {
			for (size_t i = 0; i < front_pop.size(); ++i) {
				m_pop_history_front_sols.emplace_back(std::make_shared<Solution<>>(front_pop[i]));
			}
		}
		else {
			auto temp_front_sols = m_pop_history_front_sols;
			m_pop_history_front_sols.clear();
			std::vector<size_t> temp_add_inx;
			for (size_t i = 0; i < temp_front_sols.size(); ++i) {
				temp_add_inx.push_back(i);
			}
			std::vector<size_t> final_add;
			for (size_t i = 0; i < front_pop.size(); ++i) {
				bool flag = false;
				for (size_t j = 0; j < temp_front_sols.size(); ++j) {
					if (temp_front_sols[j]->dominate(front_pop[i], pro)) {
						flag = true;
						break;
					}
				}
				if (!flag) {
					final_add = temp_add_inx;
					temp_add_inx.clear();
					m_pop_history_front_sols.emplace_back(std::make_shared<Solution<>>(front_pop[i]));
					for (size_t j = 0; j < final_add.size(); ++j) {
						if (!front_pop[i].dominate(*temp_front_sols[final_add[j]], pro)) {
							temp_add_inx.push_back(final_add[j]);
						}
					}
				}
			}
			for (size_t i = 0; i < temp_add_inx.size(); ++i) {
				m_pop_history_front_sols.emplace_back(temp_front_sols[temp_add_inx[i]]);
			}
		}
		//ʹ���Ӵ�������Ⱥ��ÿ����ǰ��
		if (m_pop_gen_front_sols.empty()) {
			m_pop_gen_front_sols.emplace_back(std::make_shared<Population<Solution<>>>(front_pop));
		}
		else {
			Population<Solution<>> combined_pop;
			for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				combined_pop.append(*this->m_individuals[i]);
			}
			for (size_t i = 0; i < new_pop.size(); ++i) {
				combined_pop.append(new_pop[i]);
			}
			NDSort(combined_pop,pro);
			Population<Solution<>> fro_pop;
			for (size_t i = 0; i < combined_pop.size(); ++i) {
				if (combined_pop[i].fitness() == 0) {
					fro_pop.append(combined_pop[i]);
				}
			}
			m_pop_gen_front_sols.emplace_back(std::make_shared<Population<Solution<>>>(fro_pop));
		}
	}

	void SPMOEA::updatePopdist(size_t inx) {
		getPop()[inx].addPopdist();
	}

	void SPMOEA::updateNewPop(Population<Solution<>>& new_pop) {
		m_new_pops.clear();
		m_new_pops = new_pop;
	}

	void SPMOEA_pop::addPopdist() {
		auto& v = getPopdist();
		//�ȸ��·ֲ�״̬���ټ���ֲ����죬������ŷʽ�����ʾ
		auto dd = getPopDistribute().back();
		std::vector<Real> origin0(dd[0].size(), 0.);
		std::vector<Real> origin1(dd[1].size(), 0.);
		auto var_dist = euclideanDistance(dd[0].begin(), dd[0].end(), origin0.begin());
		auto obj_dist = euclideanDistance(dd[1].begin(), dd[1].end(), origin1.begin());
		if (empty(v)) {
			std::vector<Real> var;
			var.push_back(var_dist);
			v.emplace_back(var);
			std::vector<Real> obj;
			obj.push_back(obj_dist);
			v.emplace_back(obj);
		}
		else {
			v[0].push_back(var_dist);
			v[1].push_back(obj_dist);
		}
	}

	bool SPMOEA_pop::popConverged() {
		//��������Ⱥ����ʷ�ֲ�
		auto var_dist = getPopdist()[0];
		auto obj_dist = getPopdist()[1];
		if (var_dist.size() <= m_stag_gen) {
			return false;
		}
		else {//���ݽ�ռ��Ŀ��ռ�ǰ�صķֲ��ж�����״̬
			std::vector<Real> test_v;
			for (size_t i = var_dist.size()-m_stag_gen; i < var_dist.size(); ++i) {
				test_v.push_back(var_dist[i]);
			}
			Real ee0 = *std::max_element(test_v.begin(), test_v.end()) - *std::min_element(test_v.begin(), test_v.end());
			std::vector<Real> test_o;
			for (size_t i = obj_dist.size() - m_stag_gen; i < obj_dist.size(); ++i) {
				test_o.push_back(obj_dist[i]);
			}
			Real ee1 = *std::max_element(test_o.begin(), test_o.end()) - *std::min_element(test_o.begin(), test_o.end());
			if (ee0 < 10e-3 && ee1 < 10e-3) {
				return true;
			}
			else {
				return false;
			}
		}
	}

	bool SPMOEA::popConverged(size_t inx) {
		//��������Ⱥ����ʷ�ֲ�
		return getPop()[inx].popConverged();
	}

	void SPMOEA::updateObjSpace() {
		auto bound = getFrontObjRange();
		getMO_HLC().initialObjSpace(bound,m_num_region_obj);
	}

	bool SPMOEA::updateCluster() {
		return false;
	}

	void SPMOEA::clusterSubspace() {
		m_mo_hlc->rankClustering();
	}

	void SPMOEA::updateHistoryInfo(Population<Solution<>>& new_pop, Problem *pro) {
		//������ʷ��
		for (size_t i = 0; i < new_pop.size(); ++i) {
			m_historical_sols.emplace_back(std::make_shared<Solution<>>(new_pop[i]));
		}
		//������ʷǰ��
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < new_pop.size(); ++i) {
			if (new_pop[i].fitness() == 0) {
				temp_pop.append(new_pop[i]);
			}
		}
		if (m_history_front_sols.empty()) {
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				m_history_front_sols.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
			}
		}
		else {
			//�ȿ�ÿ�����Ƿ񱻵�ǰǰ��֧��
			std::vector<size_t> temp_his_sols(m_history_front_sols.size(), 1);
			std::vector<size_t> temp_pop_sols(temp_pop.size(), 0);
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
					if (temp_his_sols[i] == 1) {
						Dominance dominanceship = temp_pop[j].compare(*m_history_front_sols[i], CAST_CONOP(pro)->optimizeMode());
						if (dominanceship == Dominance::kNonDominated) {
							if (i == m_history_front_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else if (dominanceship == Dominance::kDominant) {
							temp_his_sols[i] = 0;
							if (i == m_history_front_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else {
							break;
						}
					}
					else {
						if (i == m_history_front_sols.size() - 1) {
							temp_pop_sols[j] = 1;
						}
					}
				}
			}
			auto temp = m_history_front_sols;
			m_history_front_sols.clear();
			for (size_t i = 0; i < temp.size(); ++i) {
				if (temp_his_sols[i] == 1) {
					m_history_front_sols.emplace_back(temp[i]);
				}
			}
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				if (temp_pop_sols[i] == 1) {
					m_history_front_sols.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
				}
			}
		}
		//����ÿ��ǰ��
		m_gen_front_sols.emplace_back(std::make_shared<Population<Solution<>>>(temp_pop));
	}

	void SPMOEA::updateArchive(size_t num, Problem *pro) {
		//����ʷǰ����ѡ��һ�������ĸ�����Ϊarchive
		//����ѡ��ʽ��1��ӵ��ѡ��2���ο�����ѡ��3���ӿռ���ѡ��
		if (m_history_front_sols.size() <= num) {
			m_archive.clear();
			for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
				m_archive.emplace_back(m_history_front_sols[i]);
			}
		}
		else if (m_history_front_sols.size() <= 10*num) {
			std::vector<size_t> selected_inx;
			//�ȼ���ĳһĿ���µ���ֵ
			Real min_v = INT16_MAX;
			size_t add_inx = 0;
			for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
				if (m_history_front_sols[i]->objective()[0] < min_v) {
					min_v = m_history_front_sols[i]->objective()[0];
					add_inx = i;
				}
			}
			selected_inx.push_back(add_inx);
			//�ȵõ�δѡ������
			std::vector<size_t> no_select_inx;//δѡ������
			for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
				if (std::find(selected_inx.begin(), selected_inx.end(), i) == selected_inx.end()) {
					no_select_inx.push_back(i);
				}
			}
			while (selected_inx.size() < num) {//ѡ���뵱ǰ��ѡ����Զ�ĵ�
				//����δѡ������ѡ�����С��������ֵ
				std::vector<Real> min_dist;
				for (size_t i = 0; i < no_select_inx.size(); ++i) {
					std::vector<Real> temp_dist;
					for (size_t j = 0; j < selected_inx.size(); ++j) {
						auto p1 = m_history_front_sols[no_select_inx[i]]->objective();
						auto p2 = m_history_front_sols[selected_inx[j]]->objective();
						temp_dist.push_back(euclideanDistance(p1.begin(), p1.end(), p2.begin()));
					}
					Real min_d = *std::min_element(temp_dist.begin(), temp_dist.end());
					min_dist.push_back(min_d);
				}
				auto s_inx = std::distance(min_dist.begin(), std::max_element(min_dist.begin(), min_dist.end()));
				selected_inx.push_back(no_select_inx[s_inx]);
				no_select_inx.erase(no_select_inx.begin() + s_inx);
			}
			m_archive.clear();
			for (size_t i = 0; i < selected_inx.size(); ++i) {
				m_archive.emplace_back(m_history_front_sols[selected_inx[i]]);
			}
		}
		else {
			//Ŀ���ӿռ�ѡ��
			std::vector<size_t> select_index = selectIndiFromFront(m_history_front_sols, num, pro);
			m_archive.clear();
			for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
				if (std::find(select_index.begin(), select_index.end(), i) != select_index.end()) {
					m_archive.emplace_back(m_history_front_sols[i]);
				}
			}
		}
		//std::cout << "�ۻ�ǰ�ؽ����" << m_history_front_sols.size() << std::endl;
	}

	int SPMOEA::subspaceSeperable(size_t inx,Problem *pro) {
		int flag = 0;
		//����ӿռ����Ƿ���ж���ֲ��ṹ,���ӿռ���ʷ���������
		auto& front_ind = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		auto& his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;

		
		if (his_ind.size() > 40) {
			//�����ܶȾ���:DBSCAN
		    //���е��ƽ����С����ֵ
			std::vector<Real> min_dist;
			for (size_t i = 0; i < his_ind.size(); ++i) {
				std::vector<Real> temp_dist;
				for (size_t j = 0; j < his_ind.size(); ++j) {
					if (j != i) {
						auto dist = euclideanDistance(his_ind[i]->objective().begin(), his_ind[i]->objective().end(), his_ind[j]->objective().begin());
						temp_dist.push_back(dist);
					}
				}
				min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
			}
			Real sum_dist = 0.;
			for (size_t i = 0; i < min_dist.size(); ++i) {
				sum_dist += min_dist[i];
			}
			Real mean_dist = sum_dist / his_ind.size();
			//��ƽ��������Ϊ�ܶȾ���ľ���ֵ
			size_t minPts = 5;//ÿ��������С�ĸ�����
			Real epsilon = mean_dist;
			std::vector<std::vector<Real>> cluster_sol;
			for (size_t i = 0; i < his_ind.size(); ++i) {
				cluster_sol.emplace_back(his_ind[i]->variable().vect());
			}
			DBSCAN dscluster(minPts, epsilon, cluster_sol);
			dscluster.run();
			std::vector<int> cluster_id;//ÿ������������,��ʵ���δ�����
			for (size_t i = 0; i < dscluster.m_points.size(); ++i) {
				cluster_id.push_back(dscluster.m_points[i]->clusterID);
			}
			std::vector<int> cluster_num;//���е���
			cluster_num.push_back(cluster_id[0]);
			for (size_t i = 1; i < cluster_id.size(); ++i) {
				if (std::find(cluster_num.begin(), cluster_num.end(), cluster_id[i]) == cluster_num.end()) {
					cluster_num.push_back(cluster_id[i]);
				}
			}
			//����ͬ����ĵ����һ��
			std::vector<std::vector<size_t>> ind2cluster;
			for (size_t i = 0; i < cluster_num.size(); ++i) {
				std::vector<size_t> temp_cluster;
				for (size_t j = 0; j < cluster_id.size(); ++j) {
					if (cluster_id[j] == cluster_num[i]) {
						temp_cluster.push_back(j);
					}
				}
				ind2cluster.emplace_back(temp_cluster);
			}
			//���ݷ��࣬�γ��ж�
			std::vector<size_t> effective_cluster;
			for (size_t i = 0; i < cluster_num.size(); ++i) {
				if (cluster_num[i] > 0) {
					effective_cluster.push_back(cluster_num[i]);
				}
			}
			if (effective_cluster.size() > 1) {
				flag = 2;
			}
		}
		

		////�������
		//std::vector<std::vector<Real>*> objs;
		//for (size_t i = 0; i < his_ind.size(); ++i) {
		//	objs.emplace_back(&his_ind[i]->objective());
		//}
		//std::vector<int> rank;
		//int layer_num = ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		//for (size_t i = 0; i < his_ind.size(); ++i) {
		//	his_ind[i]->setFitness(rank[i]);
		//}
		////����������
		//size_t count = 0;
		//std::vector<std::vector<size_t>> layers_inx;
		//while (count < layer_num) {
		//	std::vector<size_t> temp_layer;
		//	for (size_t i = 0; i < his_ind.size(); ++i) {
		//		if (rank[i] == count) {
		//			temp_layer.push_back(i);
		//		}
		//	}
		//	layers_inx.emplace_back(temp_layer);
		//	count++;
		//}
		////�ӵ�2�㿪ʼ�жϣ�ǰһ��֧�䵱ǰ�����С����
		//for (size_t i = 1; i < std::floor(layers_inx.size() / 2); ++i) {
		//	std::vector<size_t> temp1 = layers_inx[i - 1];
		//	std::vector<size_t> temp2 = layers_inx[i];
		//	std::vector<Real> ind_dominant_dist;
		//	for (size_t j = 0; j < temp2.size(); ++j) {
		//		std::vector<size_t> dominant_ind;
		//		for (size_t k = 0; k < temp1.size(); ++k) {
		//			if (his_ind[temp1[k]]->dominate(*his_ind[temp2[j]], pro)) {
		//				dominant_ind.push_back(temp1[k]);
		//			}
		//		}
		//		//�����ռ�ľ���
		//		std::vector<Real> temp_dist;
		//		for (size_t k = 0; k < dominant_ind.size(); ++k) {
		//			temp_dist.push_back(euclideanDistance(his_ind[temp2[j]]->variable().begin(), his_ind[temp2[j]]->variable().end(), his_ind[dominant_ind[k]]->variable().begin()));
		//		}
		//		ind_dominant_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
		//	}
		//	//����֧���ľ���ֵ

		//}

		return flag;
	}

	void SPMOEA::divideSubspace(size_t inx,size_t num) {
		getMO_HLC().spaceDivide(inx, num);
	}

	void SPMOEA::splitSubspace(size_t inx, int dim, Real pos) {
		getMO_HLC().spaceSplit(inx,dim,pos);
	}

	void SPMOEA_pop::envirSelection(Population<Solution<>>& pop, Population<Solution<>>& compound_pop,Problem* pro) {
		//survivorSelection(pop, compound_pop);
		envirSelectionByNSGAII(pop, compound_pop,pro->optimizeMode());
	}

	void SPMOEA_pop::envirSelection(Problem *pro, Random *rnd, bool b) {
		if (b) {//��������ʱ�Ļ���ѡ������ӿռ�
			//����Ŀ���ӿռ����ѡ��ѡ������Ϊ��ǰ���ӿռ�ľ����Ժ�֧��ռ��������
			//�ȶԸ������Ӵ����з�֧������
			/*for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				m_offspring[i + this->m_individuals.size()] = *this->m_individuals[i];
			}*/
			//nondominatedSorting(m_offspring);
			popNondominatedSorting(m_offspring,pro->optimizeMode());
			/*auto select_pop = selectIndi(m_offspring, this->m_individuals.size(),pro);
			for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				*this->m_individuals[i] = m_offspring[select_pop[i]];
				this->m_individuals[i]->setCounter(this->m_individuals[i]->surviveAge() + 1);
			}*/
		}
		else {
			//survivorSelection(*this, m_offspring);
			envirSelectionByNSGAII(*this, m_offspring, pro->optimizeMode());
			for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				this->m_individuals[i]->setTimeEvaluate(this->m_individuals[i]->timeEvaluate() + 1);
			}
		}
	}

	void SPMOEA_pop::envirSelection(std::vector<std::pair<Real, Real>>& bound, Problem *pro, Random *rnd) {

	}

	void SPMOEA::repairSol(std::vector<Real>& sol, std::vector<std::pair<Real, Real>>& sub_bound, Random *rnd) {
		/*size_t var_num = sub_bound.size();
		std::vector<std::pair<Real, Real>> m_var_boundary;
		for (size_t i = 0; i < var_num; ++i) {
			sub_bound.emplace_back(CAST_CONOP(pro)->range(i));
		}*/
		for (size_t i = 0; i < sol.size(); ++i) {
			if (sol[i] < sub_bound[i].first || (!isnormal(sol[i]))) {
				sol[i] = sub_bound[i].first + 0.5 * rnd->uniform.next() * (sub_bound[i].second - sub_bound[i].first);
			}
			if (sol[i] > sub_bound[i].second || (!isnormal(sol[i]))) {
				sol[i] = sub_bound[i].second - 0.5 * rnd->uniform.next() * (sub_bound[i].second - sub_bound[i].first);
			}

			/*if (sol[i] < m_pop_var_range[i].first || (!isnormal(sol[i]))) {
				sol[i] = m_pop_var_range[i].first + 0.5 * rnd->uniform.next() * (m_pop_var_range[i].second - m_pop_var_range[i].first);
			}
			if (sol[i] > m_pop_var_range[i].second || (!isnormal(sol[i]))) {
				sol[i] = m_pop_var_range[i].second - 0.5 * rnd->uniform.next() * (m_pop_var_range[i].second - m_pop_var_range[i].first);
			}*/
		}
	}

	void SPMOEA::repairSol(std::vector<Real>& sol, std::vector<size_t> space_inx, Random *rnd) {
		bool flag = false;
		for (size_t i = 0; i < space_inx.size(); ++i) {
			auto bound = getMO_HLC().subspaceTree().getBox(space_inx[i]);
			if (normalSol(sol, bound)) {
				flag = true;
				break;
			}
		}
		if (!flag) {
			size_t inx = (size_t)std::floor(rnd->uniform.next() * space_inx.size());
			auto bound = getMO_HLC().subspaceTree().getBox(space_inx[inx]);
			for (size_t i = 0; i < sol.size(); ++i) {
				sol[i] = bound[i].first + rnd->uniform.next() * (bound[i].second - bound[i].first);
			}
		}
	}

	bool SPMOEA::normalSol(std::vector<Real>& sol, std::vector<std::pair<Real, Real>>& sub_bound) {
		bool flag = true;
		for (size_t i = 0; i < sol.size(); ++i) {
			if (sol[i] < sub_bound[i].first || sol[i] > sub_bound[i].second||(!isnormal(sol[i]))) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	Real SPMOEA::spaceCoverRatio(size_t inx, size_t subspace_num) {
		auto bound= getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> subspace_ratio(subspace_num, 1. / subspace_num);
		KDTree space_tree(subspace_ratio,bound);
		space_tree.buildIndex();
		std::vector<int> space_count(subspace_num, 0);
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		for (size_t i = 0; i < his_ind.size(); ++i) {
			auto pos = his_ind[i]->variable().vect();
			int ind_inx = space_tree.getRegionIdx(pos);
			space_count[ind_inx]++;
		}
		size_t total_num = 0;
		for (size_t i = 0; i < space_count.size(); ++i) {
			if (space_count[i] > 0) {
				total_num += 1;
			}
		}
		return (Real)total_num / subspace_num;
	}

	void SPMOEA::updateEE(size_t explore_num, size_t exploit_num) {
		if (m_explore_exploit_ratio.empty()) {
			std::vector<Real> explore_nums;
			explore_nums.push_back((Real)explore_num/(explore_num+exploit_num));
			m_explore_exploit_ratio.emplace_back(explore_nums);
			std::vector<Real> exploit_nums;
			exploit_nums.push_back((Real)exploit_num / (explore_num + exploit_num));
			m_explore_exploit_ratio.emplace_back(exploit_nums);
		}
		else {
			m_explore_exploit_ratio[0].push_back((Real)explore_num / (explore_num + exploit_num));
			m_explore_exploit_ratio[1].push_back((Real)exploit_num / (explore_num + exploit_num));
		}
	}

	//��ǰ�ظ�����ѡ��
	std::vector<size_t> SPMOEA::selectIndiFromFront(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd) {
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		auto ind_space_att = spaceAttach(pop, m_num_region_obj, pro, true);
		std::vector<size_t> selected_index;//�Ѿ�ѡ����������
		std::vector<size_t> front_space_index;//ǰ���ӿռ�����
		std::map<size_t, std::vector<size_t>> space_select_inx;//�ӿռ���ѡ��������
		std::map<size_t, std::vector<size_t>> space_no_select_inx;//�ӿռ�δѡ��������
		std::map<size_t, std::vector<size_t>> nei_space_select_inx;//�����ӿռ���ѡ������ʵ����
		std::map<size_t, Matrix> subspace_dist_matrix;//�ӿռ�δѡ���������ӿռ���ѡ��ľ������
		//size_t front_space_num = 0;//ǰ���ӿռ��ڵĸ�������
		for (const auto& ind : ind_space_att[0]) {
			front_space_index.emplace_back(ind.first);
			//front_space_num += ind.second.size();
		}
		//��ѡ��ǰ���ӿռ�ĳһ��Ŀ��ļ�ֵ��
		size_t obj_inx = 0;
		for (auto& sub : ind_space_att[0]) {
			int select_inx = -1;
			if (sub.second.size() > 1) {
				std::vector<Real> sub_obj;
				for (size_t i = 0; i < sub.second.size(); ++i) {
					sub_obj.push_back(pop[sub.second[i]].objective()[obj_inx]);
				}
				select_inx = std::distance(sub_obj.begin(), std::min_element(sub_obj.begin(), sub_obj.end()));
			}
			else {
				select_inx = 0;
			}
			std::vector<size_t> selected_inx;
			selected_inx.push_back(select_inx);
			selected_index.push_back(sub.second[select_inx]);
			std::vector<size_t> no_select_inx;
			for (size_t i = 0; i < sub.second.size(); ++i) {
				if (i != select_inx) {
					no_select_inx.push_back(i);
				}
			}
			space_select_inx.insert(std::make_pair(sub.first, selected_inx));
			space_no_select_inx.insert(std::make_pair(sub.first, no_select_inx));
		}
		//�ٳ�ʼ���ӿռ��������������ѡ��
		for (auto& sub : ind_space_att[0]) {
			std::list<size_t> neighbors;
			std::vector<size_t> front_neighbors;//ǰ������
			front_neighbors.push_back(sub.first);
			m_mo_hlc->getObjspaceTree().findNeighbor(sub.first, neighbors);
			for (size_t i = 0; i < front_space_index.size(); ++i) {
				if (std::find(neighbors.begin(), neighbors.end(), front_space_index[i]) != neighbors.end()) {
					front_neighbors.push_back(front_space_index[i]);
				}
			}
			std::vector<size_t> temp_neigh_select_inx;//������ѡ�������ʵ����
			for (size_t i = 0; i < front_neighbors.size(); ++i) {
				auto neigh_select_inx = space_select_inx[front_neighbors[i]];
				auto neigh_ind_inx = ind_space_att[0][front_neighbors[i]];
				for (size_t j = 0; j < neigh_select_inx.size(); ++j) {
					temp_neigh_select_inx.push_back(neigh_ind_inx[neigh_select_inx[j]]);
				}
			}
			nei_space_select_inx.insert(std::make_pair(sub.first, temp_neigh_select_inx));
			auto subspace_no_select_inx = space_no_select_inx[sub.first];//�ӿռ�δѡ���������
			Matrix dist_matrix(subspace_no_select_inx.size(), temp_neigh_select_inx.size());
			for (size_t i = 0; i < subspace_no_select_inx.size(); ++i) {
				for (size_t j = 0; j < temp_neigh_select_inx.size(); ++j) {
					auto p1 = pop[sub.second[subspace_no_select_inx[i]]].objective();
					auto p2 = pop[temp_neigh_select_inx[j]].objective();
					dist_matrix[i][j] = euclideanDistance(p1.begin(), p1.end(), p2.begin());
				}
			}
			subspace_dist_matrix.insert(std::make_pair(sub.first, dist_matrix));
		}
		//Ȼ�����ÿ���ӿռ���δѡ������������ѡ����ľ���ͳ��ֵѡ��Ҫѡ���ӿռ�
		size_t subspace_inx = 0;
		while (selected_index.size() < select_num) {
			//����ÿ���ӿռ���δѡ������������ѡ����ľ�����С��������ֵ������Ҹ�����ӿռ�
			std::vector<std::tuple<size_t, Real, int>> compare_info;//�ֱ��ǣ��ӿռ��������ӿռ���С��������ֵ����Ӧ��δѡ���������
			for (auto& sub : ind_space_att[0]) {
				Real row_min_max;
				int row_min_max_inx;
				if (subspace_dist_matrix[sub.first].data().size() == 0) {
					row_min_max = -1;
					row_min_max_inx = -1;
				}
				else if (subspace_dist_matrix[sub.first].data().size() != 0 && subspace_dist_matrix[sub.first].data()[0].size() == 0){
					row_min_max = -1;
					row_min_max_inx = -1;
				}
				else {
					std::vector<Real> row_min_dist;
					for (size_t i = 0; i < subspace_dist_matrix[sub.first].data().size(); i++) {
						auto temp_dist = subspace_dist_matrix[sub.first][i].vect();
						row_min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
					}
					row_min_max = *std::max_element(row_min_dist.begin(), row_min_dist.end());
					row_min_max_inx = std::distance(row_min_dist.begin(), std::max_element(row_min_dist.begin(), row_min_dist.end()));
				}
				compare_info.emplace_back(std::tuple<size_t, Real, int>(sub.first, row_min_max, row_min_max_inx));
			}
			//��������ӿռ���ѡ������
			for (size_t i = 0; i < compare_info.size() - 1; ++i) {
				for (size_t j = 0; j < compare_info.size() - 1 - i; ++j) {
					if (std::get<1>(compare_info[j]) < std::get<1>(compare_info[j + 1])) {
						auto temp = compare_info[j];
						compare_info[j] = compare_info[j + 1];
						compare_info[j + 1] = temp;
					}
				}
			}
			size_t count = 0;
			for (size_t i = 0; i < compare_info.size(); ++i) {
				if (!space_no_select_inx[std::get<0>(compare_info[i])].empty()) {
					subspace_inx = std::get<0>(compare_info[i]);
					count = i;
					break;
				}
			}
			//��ѡ����ӿռ���ѡ��������ĸ���
			auto ind_inx = space_no_select_inx[subspace_inx][std::get<2>(compare_info[count])];
			auto real_select_ind_inx = ind_space_att[0][subspace_inx][ind_inx];
			selected_index.push_back(real_select_ind_inx);
			//����ÿ���ӿռ����ѡ����Ϣ
			std::list<size_t> neighbors;
			std::vector<size_t> front_neighbors;//ǰ������
			m_mo_hlc->getObjspaceTree().findNeighbor(subspace_inx, neighbors);
			for (size_t i = 0; i < front_space_index.size(); ++i) {
				if (std::find(neighbors.begin(), neighbors.end(), front_space_index[i]) != neighbors.end()) {
					front_neighbors.push_back(front_space_index[i]);
				}
			}
			//���������ӿռ���Ϣ
			for (size_t i = 0; i < front_neighbors.size(); ++i) {
				//���һ��
				auto no_select_ind_inx = space_no_select_inx[front_neighbors[i]];//���ӿռ�δѡ������
				std::vector<size_t> real_ind_inx;
				for (size_t j = 0; j < no_select_ind_inx.size(); ++j) {
					real_ind_inx.push_back(ind_space_att[0][front_neighbors[i]][no_select_ind_inx[j]]);
				}
				auto row_num = subspace_dist_matrix[front_neighbors[i]].data().size();
				for (size_t j = 0; j < row_num; ++j) {
					auto p1 = pop[real_ind_inx[j]].objective();
					auto p2 = pop[real_select_ind_inx].objective();
					auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
					subspace_dist_matrix[front_neighbors[i]][j].pushBack(dist);
				}
				//����������ѡ��
				nei_space_select_inx[front_neighbors[i]].push_back(real_select_ind_inx);
			}
			//������ѡ�ӿռ���Ϣ
			//ȥ����ѡ��������
			auto it1 = subspace_dist_matrix[subspace_inx].data().begin();
			subspace_dist_matrix[subspace_inx].data().erase(it1 + std::get<2>(compare_info[count]));
			//�����ӿռ���ѡ��
			space_select_inx[subspace_inx].push_back(ind_inx);
			//�����ӿռ�δѡ��
			auto it2 = space_no_select_inx[subspace_inx].begin();
			space_no_select_inx[subspace_inx].erase(it2 + std::get<2>(compare_info[count]));
			//����������ѡ��
			nei_space_select_inx[subspace_inx].push_back(real_select_ind_inx);
			//�������ѡ�����
			std::vector<size_t> no_select_real_inx;
			for (size_t j = 0; j < space_no_select_inx[subspace_inx].size(); ++j) {
				no_select_real_inx.push_back(ind_space_att[0][subspace_inx][space_no_select_inx[subspace_inx][j]]);
			}
			for (size_t j = 0; j < no_select_real_inx.size(); ++j) {
				auto p1 = pop[no_select_real_inx[j]].objective();
				auto p2 = pop[real_select_ind_inx].objective();
				auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
				subspace_dist_matrix[subspace_inx][j].pushBack(dist);
			}
		}
		return selected_index;
	}

	//��ǰ�ظ�����ѡ��
	std::vector<size_t> SPMOEA::selectIndiFromFront(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t select_num, Problem *pro) {
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		auto ind_space_att = spaceAttach(pop,m_num_region_obj,pro,true);//ÿ���ӿռ�����Щ����
		std::vector<size_t> selected_index;//�Ѿ�ѡ����������
		std::vector<size_t> front_space_index;//ǰ���ӿռ�����
		std::map<size_t, std::vector<size_t>> space_select_inx;//�ӿռ���ѡ��������
		std::map<size_t, std::vector<size_t>> space_no_select_inx;//�ӿռ�δѡ��������
		std::map<size_t, std::vector<size_t>> nei_space_select_inx;//�����ӿռ���ѡ������ʵ����
		std::map<size_t, Matrix> subspace_dist_matrix;//�ӿռ�δѡ���������ӿռ���ѡ��ľ������
		//size_t front_space_num = 0;//ǰ���ӿռ��ڵĸ�������
		for (const auto& ind : ind_space_att[0]) {
			front_space_index.emplace_back(ind.first);
			//front_space_num += ind.second.size();
		}
		//��ѡ��ǰ���ӿռ�ĳһ��Ŀ��ļ�ֵ��
		size_t obj_inx = 0;
		for (auto& sub : ind_space_att[0]) {
			int select_inx = -1;
			if (sub.second.size() > 1) {
				std::vector<Real> sub_obj;
				for (size_t i = 0; i < sub.second.size(); ++i) {
					sub_obj.push_back(pop[sub.second[i]]->objective()[obj_inx]);
				}
				select_inx = std::distance(sub_obj.begin(), std::min_element(sub_obj.begin(), sub_obj.end()));
			}
			else {
				select_inx = 0;
			}
			std::vector<size_t> selected_inx;
			selected_inx.push_back(select_inx);
			selected_index.push_back(sub.second[select_inx]);
			std::vector<size_t> no_select_inx;
			for (size_t i = 0; i < sub.second.size(); ++i) {
				if (i != select_inx) {
					no_select_inx.push_back(i);
				}
			}
			space_select_inx.insert(std::make_pair(sub.first, selected_inx));
			space_no_select_inx.insert(std::make_pair(sub.first, no_select_inx));
		}
		//�ٳ�ʼ���ӿռ��������������ѡ��
		for (auto& sub : ind_space_att[0]) {
			std::list<size_t> neighbors;
			std::vector<size_t> front_neighbors;//ǰ������
			m_mo_hlc->getObjspaceTree().findNeighbor(sub.first, neighbors);
			front_neighbors.push_back(sub.first);
			for (size_t i = 0; i < front_space_index.size(); ++i) {
				if (std::find(neighbors.begin(), neighbors.end(), front_space_index[i]) != neighbors.end()) {
					front_neighbors.push_back(front_space_index[i]);
				}
			}
			std::vector<size_t> temp_neigh_select_inx;//������ѡ�������ʵ����
			for (size_t i = 0; i < front_neighbors.size(); ++i) {
				auto neigh_select_inx = space_select_inx[front_neighbors[i]];
				auto neigh_ind_inx = ind_space_att[0][front_neighbors[i]];
				for (size_t j = 0; j < neigh_select_inx.size(); ++j) {
					temp_neigh_select_inx.push_back(neigh_ind_inx[neigh_select_inx[j]]);
				}
			}
			nei_space_select_inx.insert(std::make_pair(sub.first, temp_neigh_select_inx));
			auto subspace_no_select_inx = space_no_select_inx[sub.first];//�ӿռ�δѡ���������
			Matrix dist_matrix(subspace_no_select_inx.size(), temp_neigh_select_inx.size());
			for (size_t i = 0; i < subspace_no_select_inx.size(); ++i) {
				for (size_t j = 0; j < temp_neigh_select_inx.size(); ++j) {
					auto p1 = pop[sub.second[subspace_no_select_inx[i]]]->objective();
					auto p2 = pop[temp_neigh_select_inx[j]]->objective();
					dist_matrix[i][j] = euclideanDistance(p1.begin(), p1.end(), p2.begin());
				}
			}
			subspace_dist_matrix.insert(std::make_pair(sub.first, dist_matrix));
		}
		//Ȼ�����ÿ���ӿռ���δѡ������������ѡ����ľ���ͳ��ֵѡ��Ҫѡ���ӿռ�
		size_t subspace_inx = 0;
		while (selected_index.size() < select_num) {
			//����ÿ���ӿռ���δѡ������������ѡ����ľ�����С��������ֵ������Ҹ�����ӿռ�
			std::vector<std::tuple<size_t, Real, int>> compare_info;//�ֱ��ǣ��ӿռ��������ӿռ���С��������ֵ����Ӧ��δѡ���������
			for (auto& sub : ind_space_att[0]) {
				Real row_min_max;
				int row_min_max_inx;
				if (m_evaluations > 21000) {
					size_t a = 1;
				}
				if (subspace_dist_matrix[sub.first].data().size() == 0 || (subspace_dist_matrix[sub.first].data().size() != 0 && subspace_dist_matrix[sub.first].data()[0].size() == 0)) {
					row_min_max = -1;
					row_min_max_inx = -1;
				}
				else {
					std::vector<Real> row_min_dist;
					for (size_t i = 0; i < subspace_dist_matrix[sub.first].data().size(); i++) {
						auto temp_dist = subspace_dist_matrix[sub.first][i].vect();
						row_min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
					}
					row_min_max = *std::max_element(row_min_dist.begin(), row_min_dist.end());
					row_min_max_inx = std::distance(row_min_dist.begin(), std::max_element(row_min_dist.begin(), row_min_dist.end()));
				}
				compare_info.emplace_back(std::tuple<size_t, Real, int>(sub.first, row_min_max, row_min_max_inx));
			}
			//��������ӿռ���ѡ������
			for (size_t i = 0; i < compare_info.size() - 1; ++i) {
				for (size_t j = 0; j < compare_info.size() - 1 - i; ++j) {
					if (std::get<1>(compare_info[j]) < std::get<1>(compare_info[j + 1])) {
						auto temp = compare_info[j];
						compare_info[j] = compare_info[j + 1];
						compare_info[j + 1] = temp;
					}
				}
			}
			size_t count = 0;
			for (size_t i = 0; i < compare_info.size(); ++i) {
				if (!space_no_select_inx[std::get<0>(compare_info[i])].empty()) {
					subspace_inx = std::get<0>(compare_info[i]);
					count = i;
					break;
				}
			}
			//��ѡ����ӿռ���ѡ��������ĸ���
			auto ind_inx = space_no_select_inx[subspace_inx][std::get<2>(compare_info[count])];
			auto real_select_ind_inx = ind_space_att[0][subspace_inx][ind_inx];
			selected_index.push_back(real_select_ind_inx);
			//����ÿ���ӿռ����ѡ����Ϣ
			std::list<size_t> neighbors;
			std::vector<size_t> front_neighbors;//ǰ������
			m_mo_hlc->getObjspaceTree().findNeighbor(subspace_inx, neighbors);
			for (size_t i = 0; i < front_space_index.size(); ++i) {
				if (std::find(neighbors.begin(), neighbors.end(), front_space_index[i]) != neighbors.end()) {
					front_neighbors.push_back(front_space_index[i]);
				}
			}
			//���������ӿռ���Ϣ
			for (size_t i = 0; i < front_neighbors.size(); ++i) {
				//���һ��
				auto no_select_ind_inx = space_no_select_inx[front_neighbors[i]];//���ӿռ�δѡ������
				std::vector<size_t> real_ind_inx;
				for (size_t j = 0; j < no_select_ind_inx.size(); ++j) {
					real_ind_inx.push_back(ind_space_att[0][front_neighbors[i]][no_select_ind_inx[j]]);
				}
				auto row_num = subspace_dist_matrix[front_neighbors[i]].data().size();
				for (size_t j = 0; j < row_num; ++j) {
					auto p1 = pop[real_ind_inx[j]]->objective();
					auto p2 = pop[real_select_ind_inx]->objective();
					auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
					subspace_dist_matrix[front_neighbors[i]][j].pushBack(dist);
				}
				//����������ѡ��
				nei_space_select_inx[front_neighbors[i]].push_back(real_select_ind_inx);
			}
			//������ѡ�ӿռ���Ϣ
			//ȥ����ѡ��������
			auto it1 = subspace_dist_matrix[subspace_inx].data().begin();
			subspace_dist_matrix[subspace_inx].data().erase(it1 + std::get<2>(compare_info[count]));
			//�����ӿռ���ѡ��
			space_select_inx[subspace_inx].push_back(ind_inx);
			//�����ӿռ�δѡ��
			auto it2 = space_no_select_inx[subspace_inx].begin();
			space_no_select_inx[subspace_inx].erase(it2 + std::get<2>(compare_info[count]));
			//����������ѡ��
			nei_space_select_inx[subspace_inx].push_back(real_select_ind_inx);
			//�������ѡ�����
			std::vector<size_t> no_select_real_inx;
			for (size_t j = 0; j < space_no_select_inx[subspace_inx].size(); ++j) {
				no_select_real_inx.push_back(ind_space_att[0][subspace_inx][space_no_select_inx[subspace_inx][j]]);
			}
			for (size_t j = 0; j < no_select_real_inx.size(); ++j) {
				auto p1 = pop[no_select_real_inx[j]]->objective();
				auto p2 = pop[real_select_ind_inx]->objective();
				auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
				subspace_dist_matrix[subspace_inx][j].pushBack(dist);
			}
		}
		return selected_index;
	}

	//��С��󷽷���ǰ�ص���ѡ��
	std::vector<size_t> SPMOEA::selectMaxMinFromFront(const Population<Solution<>>& pop, size_t select_num) {
		std::vector<size_t> select_index;//�Ѿ�ѡ����������
		std::vector<size_t> remain_index;
		//��ѡ��ĳ��Ŀ������С��һ����
		size_t obj_inx = 0;
		size_t selec_inx = 0;
		Real min_v = INT16_MAX;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (min_v < pop[i].objective()[obj_inx]) {
				min_v = pop[i].objective()[obj_inx];
				selec_inx = i;
			}
		}
		select_index.push_back(selec_inx);
		for (size_t i = 0; i < pop.size(); ++i) {
			if (i != selec_inx) {
				remain_index.push_back(i);
			}
		}
		//��ѡ������ѡ�����������Զ�ĵ�
		std::vector<std::vector<Real>> remain2select;
		std::vector<Real> dist;
		for (size_t i = 0; i < remain_index.size(); ++i) {
			auto p1 = pop[select_index.back()].objective();
			auto p2 = pop[remain_index[i]].objective();
			auto temp_dist = euclideanDistance(p1.begin(),p1.end(),p2.begin());
			dist.push_back(temp_dist);
		}
		remain2select.emplace_back(dist);
		while (select_index.size() < select_num) {
			//�ҳ�δѡ�㵽��ѡ�����С��������ֵ
			std::vector<std::vector<Real>> dist2select;
			for (size_t i = 0; i < remain2select[0].size(); ++i) {
				std::vector<Real> temp_dist;
				for (size_t j = 0; j < select_index.size(); ++j) {
					temp_dist.push_back(remain2select[j][i]);
				}
				dist2select.emplace_back(temp_dist);
			}
			//�ҳ�ÿ��δѡ�����С����ֵ
			std::vector<Real> min_dist;
			for (size_t i = 0; i < dist2select.size(); ++i) {
				auto temp = *std::min_element(dist2select[i].begin(), dist2select[i].end());
				min_dist.push_back(temp);
			}
			//�ҳ���С��������ֵ
			auto sele_inx = std::distance(min_dist.begin(),std::max_element(min_dist.begin(),min_dist.end()));
			//������ѡ���δѡ��
			select_index.push_back(remain_index[sele_inx]);
			remain_index.erase(remain_index.begin()+sele_inx);
			//���δѡ�㵽��ѡ��ľ��룬ɾ��ѡ����ľ���
			for (size_t i = 0; i < remain2select.size(); ++i) {
				remain2select[i].erase(remain2select[i].begin() + sele_inx);
			}
			std::vector<Real> add_dist;
			for (size_t i = 0; i < remain_index.size(); ++i) {
				auto pp1 = pop[select_index.back()].objective();
				auto pp2 = pop[remain_index[i]].objective();
				auto temp_dist = euclideanDistance(pp1.begin(), pp1.end(), pp2.begin());
				add_dist.push_back(temp_dist);
			}
			remain2select.emplace_back(add_dist);
		}

		return select_index;
	}

	std::vector<size_t> SPMOEA::selectMaxMinFromFront(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t select_num) {
		std::vector<size_t> select_index;//�Ѿ�ѡ����������
		std::vector<size_t> remain_index;
		//��ѡ��ĳ��Ŀ������С��һ����
		size_t obj_inx = 0;
		size_t selec_inx = 0;
		Real min_v = INT16_MAX;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (min_v < pop[i]->objective()[obj_inx]) {
				min_v = pop[i]->objective()[obj_inx];
				selec_inx = i;
			}
		}
		select_index.push_back(selec_inx);
		for (size_t i = 0; i < pop.size(); ++i) {
			if (i != selec_inx) {
				remain_index.push_back(i);
			}
		}
		//��ѡ������ѡ�����������Զ�ĵ�
		std::vector<std::vector<Real>> remain2select;
		std::vector<Real> dist;
		for (size_t i = 0; i < remain_index.size(); ++i) {
			auto p1 = pop[select_index.back()]->objective();
			auto p2 = pop[remain_index[i]]->objective();
			auto temp_dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
			dist.push_back(temp_dist);
		}
		remain2select.emplace_back(dist);
		while (select_index.size() < select_num) {
			//�ҳ�δѡ�㵽��ѡ�����С��������ֵ
			std::vector<std::vector<Real>> dist2select;
			for (size_t i = 0; i < remain2select[0].size(); ++i) {
				std::vector<Real> temp_dist;
				for (size_t j = 0; j < select_index.size(); ++j) {
					temp_dist.push_back(remain2select[j][i]);
				}
				dist2select.emplace_back(temp_dist);
			}
			//�ҳ�ÿ��δѡ�����С����ֵ
			std::vector<Real> min_dist;
			for (size_t i = 0; i < dist2select.size(); ++i) {
				auto temp = *std::min_element(dist2select[i].begin(), dist2select[i].end());
				min_dist.push_back(temp);
			}
			//�ҳ���С��������ֵ
			auto sele_inx = std::distance(min_dist.begin(),std::max_element(min_dist.begin(), min_dist.end()));
			//������ѡ���δѡ��
			select_index.push_back(remain_index[sele_inx]);
			remain_index.erase(remain_index.begin() + sele_inx);
			//���δѡ�㵽��ѡ��ľ��룬ɾ��ѡ����ľ���
			for (size_t i = 0; i < remain2select.size(); ++i) {
				remain2select[i].erase(remain2select[i].begin() + sele_inx);
			}
			std::vector<Real> add_dist;
			for (size_t i = 0; i < remain_index.size(); ++i) {
				auto pp1 = pop[select_index.back()]->objective();
				auto pp2 = pop[remain_index[i]]->objective();
				auto temp_dist = euclideanDistance(pp1.begin(), pp1.end(), pp2.begin());
				add_dist.push_back(temp_dist);
			}
			remain2select.emplace_back(add_dist);
		}

		return select_index;
	}

	
	std::vector<size_t> SPMOEA::selectMaxMinFromFront(const std::vector <std::vector<Real>>& data, size_t select_num) {
		//�����ݹ�һ��
		auto normal_data = data;
		dataNormalize(normal_data);
		std::vector<size_t> select_index;//�Ѿ�ѡ����������
		std::vector<size_t> remain_index;
		//��ѡ��ĳ��Ŀ������С��һ����
		size_t obj_inx = 0;
		size_t selec_inx = 0;
		Real min_v = std::numeric_limits<Real>::max();
		for (size_t i = 0; i < normal_data.size(); ++i) {
			if (min_v < normal_data[i][obj_inx]) {
				min_v = normal_data[i][obj_inx];
				selec_inx = i;
			}
		}
		select_index.push_back(selec_inx);
		for (size_t i = 0; i < normal_data.size(); ++i) {
			if (i != selec_inx) {
				remain_index.push_back(i);
			}
		}
		//��ѡ������ѡ�����������Զ�ĵ�
		std::vector<std::vector<Real>> remain2select;
		std::vector<Real> dist;
		for (size_t i = 0; i < remain_index.size(); ++i) {
			auto p1 = normal_data[select_index.back()];
			auto p2 = normal_data[remain_index[i]];
			auto temp_dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
			dist.push_back(temp_dist);
		}
		remain2select.emplace_back(dist);
		while (select_index.size() < select_num) {
			//�ҳ�δѡ�㵽��ѡ�����С��������ֵ
			std::vector<std::vector<Real>> dist2select;
			for (size_t i = 0; i < remain2select[0].size(); ++i) {
				std::vector<Real> temp_dist;
				for (size_t j = 0; j < select_index.size(); ++j) {
					temp_dist.push_back(remain2select[j][i]);
				}
				dist2select.emplace_back(temp_dist);
			}
			//�ҳ�ÿ��δѡ�����С����ֵ
			std::vector<Real> min_dist;
			for (size_t i = 0; i < dist2select.size(); ++i) {
				auto temp = *std::min_element(dist2select[i].begin(), dist2select[i].end());
				min_dist.push_back(temp);
			}
			//�ҳ���С��������ֵ
			auto sele_inx = std::distance(min_dist.begin(), std::max_element(min_dist.begin(), min_dist.end()));
			//������ѡ���δѡ��
			select_index.push_back(remain_index[sele_inx]);
			remain_index.erase(remain_index.begin() + sele_inx);
			//���δѡ�㵽��ѡ��ľ��룬ɾ��ѡ����ľ���
			for (size_t i = 0; i < remain2select.size(); ++i) {
				remain2select[i].erase(remain2select[i].begin() + sele_inx);
			}
			std::vector<Real> add_dist;
			for (size_t i = 0; i < remain_index.size(); ++i) {
				auto pp1 = normal_data[select_index.back()];
				auto pp2 = normal_data[remain_index[i]];
				auto temp_dist = euclideanDistance(pp1.begin(), pp1.end(), pp2.begin());
				add_dist.push_back(temp_dist);
			}
			remain2select.emplace_back(add_dist);
		}

		return select_index;
	}

	std::vector<size_t> SPMOEA::selectMaxMinFromFront(const std::vector<std::vector<Real>>& data,std::vector<size_t> ranks, size_t select_num) {
		//��������ranks
		std::vector<std::vector<size_t>> all_ranks;
		auto max_rank = *std::max_element(ranks.begin(),ranks.end());
		all_ranks.resize(max_rank+1);
		for (size_t i = 0; i < ranks.size(); ++i) {
			auto r = ranks[i];
			all_ranks[r].push_back(i);
		}
		//���ڼ��㳬��
		size_t count = 0;
		size_t layer = 0;
		for (size_t i = 0; i < all_ranks.size(); ++i) {
			if (all_ranks[i].size() + count >= select_num) {
				layer = i;
				break;
			}
			else {
				count += all_ranks[i].size();
			}
		}
		//�ӵ�layer��ѡ
		std::vector<size_t> select_index;//�Ѿ�ѡ����������
		for (size_t i = 0; i <= layer; ++i) {
			auto inx = all_ranks[i];
			if (i < layer) {
				for (auto jj : inx) {
					select_index.push_back(jj);
				}
			}
			else {
				//��ǰ�����max-min����ѡ��ʣ�����
				std::vector<std::vector<Real>> layer_data;
				for (size_t j = 0; j < inx.size(); ++j) {
					layer_data.emplace_back(data[inx[j]]);
				}
				auto sele_inx = selectMaxMinFromFront(layer_data, select_num - select_index.size());
				for (size_t j = 0; j < sele_inx.size(); ++j) {
					select_index.push_back(inx[sele_inx[j]]);
				}
			}
		}

		return select_index;
	}

	std::vector<size_t> SPMOEA::selectByNorm(const std::vector<std::vector<Real>>& data, size_t select_num, Problem* pro, Random* rnd) {
		//���е��һ��
		auto normal_data = data;
		dataNormalize(normal_data);
		std::vector<Real> origin(normal_data[0].size(),0.);
		Real r = 0;
		Real step = 0.05;
		std::vector<size_t> select_inx;
		while (select_inx.size() < select_num) {
			//���ڷ�Χ�ڵĸ���
			size_t num = 0;
			for (size_t i = 0; i < normal_data.size(); ++i) {
				Real temp = euclideanDistance(normal_data[i].begin(), normal_data[i].end(), origin.begin());
				if (temp<r + step && temp>=r) {
					num++;
					select_inx.push_back(i);
				}
			}
			r += step;
		}
		std::vector<std::vector<Real>> candidate_objs;
		for (size_t i = 0; i < select_inx.size(); ++i) {
			candidate_objs.emplace_back(normal_data[select_inx[i]]);
		}
		while (select_inx.size() > select_num) {
			//����ɾ����ӵ����
			auto crowd_pair = crowdedPointPair(candidate_objs);
			size_t delete_inx = 0;
			auto obj1 = candidate_objs[crowd_pair.first];
			auto obj2 = candidate_objs[crowd_pair.second];
			auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
			if (ship == Dominance::kDominant) {
				delete_inx = crowd_pair.second;
			}
			else if (ship == Dominance::kDominated) {
				delete_inx = crowd_pair.first;
			}
			else if (ship == Dominance::kNonDominated) {
				Real rand = rnd->uniform.next();
				if (rand < 0.5) {
					delete_inx = crowd_pair.first;
				}
				else {
					delete_inx = crowd_pair.second;
				}
			}
			candidate_objs.erase(candidate_objs.begin() + delete_inx);
			select_inx.erase(select_inx.begin()+delete_inx);
		}
		
		return select_inx;
	}

	//Ŀ���ӿռ�ѡ��
	std::vector<size_t> SPMOEA::selectIndiFromSpace(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd) {
		//��ͳ��ǰ�Ÿ��������
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		auto pop_att = popAttach(pop,m_num_region_obj,pro);
		std::vector<size_t> select_index;//�Ѿ�ѡ����������
		std::vector<size_t> behind_space_index;//�����ӿռ�����
		std::vector<size_t> front_space_index;//ǰ���ӿռ�����
		size_t front_space_num = 0;//ǰ���ӿռ��ڵĸ�������
		std::vector<size_t> first_ind_index;//�ȵõ���һ�Ÿ��������
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].fitness() == 0)
				first_ind_index.emplace_back(i);
		}
		for (const auto& ind : pop_att[0][0]) {
			front_space_index.emplace_back(ind.first);
			front_space_num += ind.second.size();
		}
		for (const auto& ind : pop_att[0][1]) {
			behind_space_index.emplace_back(ind.first);
		}
		std::map<size_t, std::vector<size_t>> map_first_indi;//�õ�ǰ�Ÿ���ֲ�����Щ�ӿռ�
		for (const auto& sp : pop_att[0][0]) {
			std::pair<size_t, std::vector<size_t>> temp;
			temp.first = sp.first;
			for (const auto& in : sp.second) {
				if (pop[in].fitness() == 0) {
					temp.second.push_back(in);
				}
			}
			map_first_indi.insert(temp);
		}
		std::map<size_t, std::vector<size_t>> map_selected_indi;//ǰ���ӿռ��Ѿ�ѡ�Ľ⣬˳������
		std::map<size_t, std::vector<size_t>> map_no_selected_indi;//ǰ���ӿռ��Ѿ�ѡ�Ľ⣬˳������
		//����ǰ���ӿռ��ڲ��ĸ����������и���ѡ��ѡ��ʱ���Ǹ�������
		if (first_ind_index.size() > select_num) {//ǰ�ظ����Ƿ����
			Population<Solution<>> temp_pop(first_ind_index.size(),pro);
			for (size_t i = 0; i < first_ind_index.size(); ++i) {
				temp_pop[i] = pop[first_ind_index[i]];
			}
			std::vector<size_t> temp_select;
			if (first_ind_index.size() > 2 * select_num) {
				temp_select = selectIndiFromFront(temp_pop, select_num, pro, rnd);
			}
			else {
				temp_select = selectMaxMinFromFront(temp_pop, select_num);
			}
			for (size_t i = 0; i < select_num; ++i) {
				select_index.push_back(temp_select[i]);
			}
		}
		else if(front_space_num>select_num){//ǰ���ӿռ��ڵĸ������Ƿ����
			bool flag = 1;
			if (flag) {//Ŀ��ռ�ѡ��
				//��ѡǰ�ظ���
				select_index = first_ind_index;
				//������ѡ��ǰ���ӿռ��е���������
				auto front_subspace_ind = pop_att[0][0];
				for (const auto& i : front_subspace_ind) {
					std::pair<size_t, std::vector<size_t>> temp;
					temp.first = i.first;
					std::vector<size_t> f_ind_idx;
					for (size_t j = 0; j < i.second.size(); ++j) {
						if (find(first_ind_index.begin(), first_ind_index.end(), i.second[j]) != first_ind_index.end()) {
							f_ind_idx.push_back(i.second[j]);
						}
					}
					temp.second = f_ind_idx;
					map_selected_indi.insert(temp);//ǰ���ӿռ��Ѿ�ѡ�Ľ�
				}
				while (select_index.size() < select_num) {
					for (const auto& i : front_subspace_ind) {
						if (map_selected_indi[i.first].size() < i.second.size()) {
							//��ʣ�������ѡ��һ��
							for (size_t j = 0; j < i.second.size(); ++j) {
								if (find(map_selected_indi[i.first].begin(), map_selected_indi[i.first].end(), i.second[j]) == map_selected_indi[i.first].end()) {
									select_index.push_back(i.second[j]);
									map_selected_indi[i.first].push_back(i.second[j]);
									break;
								}
							}
						}
						if (select_index.size() >= select_num) {
							break;
						}
					}
				}
			}
			else {//��ռ�ѡ�񣬸��ݽ�ռ�ֲ�ѡ�񣬻��߸����ӿռ��Ǳ��ֵѡ��
				//��ǰ���ӿռ�ĵ�ֳ������֣�ǰ�ؽ�ͷ�ǰ�ؽ�

				//�����н�ӳ�䵽�����ռ䣬�ҳ���ǰ�ؽ���Զ�ĸ���
			}
		}
		else {
			//�ȼ���ǰ���ӿռ���������
			for (const auto& ind : pop_att[0][0]) {
				for (const auto& p : ind.second) {
					select_index.push_back(p);
				}
			}
			//���ں����ӿռ�����ѡ��
			std::vector<size_t> space_select_num(pop_att[0][1].size(), 0);//ÿ�������ӿռ���ѡ��ĸ���
			std::vector<std::vector<size_t>> behind_spaces;
			for (const auto& sp : pop_att[0][1]) {
				behind_spaces.push_back(sp.second);
			}
			while (select_index.size() < select_num) {
				for (size_t i = 0; i < behind_spaces.size(); ++i) {
					if (space_select_num[i] < behind_spaces[i].size()) {
						size_t select_idx = std::floor(behind_spaces[i].size() * rnd->uniform.next());
						if (std::find(select_index.begin(), select_index.end(), behind_spaces[i][select_idx]) == select_index.end()) {
							select_index.push_back(behind_spaces[i][select_idx]);
							space_select_num[i]++;
						}
					}
					if (select_index.size() >= select_num) {
						break;
					}
				}
			}
		}
		return select_index;
	}

	//Ŀ���ӿռ�ѡ��
	std::vector<size_t> SPMOEA::selectFromSpace(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd) {
		// �ȹ�һ������
		std::vector<std::vector<Real>> pop_obj;
		std::vector<std::vector<Real>> pop_var;
		for (size_t i = 0; i < pop.size(); ++i) {
			pop_obj.emplace_back(pop[i].objective());
			pop_var.emplace_back(pop[i].variable().vect());
		}
		std::vector<size_t> rank;
		for (size_t i = 0; i < pop.size(); ++i) {
			rank.push_back(pop[i].fitness());
		}
		//dataNormalize(pop_obj);
		//// �ٹ����ӿռ�
		//size_t M = CAST_CONOP(pro)->numberObjectives();
		//size_t num_obj_spaces = numObjRegion();
		//std::vector<Real> subspace_ratios(num_obj_spaces, 1. / num_obj_spaces);
		//std::vector<std::pair<Real, Real>> obj_bound;
		//for (size_t i = 0; i < M; ++i) {
		//	obj_bound.emplace_back(std::make_pair<>(0., 1.));
		//}
		//nanoflann::KDTreeSpace<Real> obj_tree(subspace_ratios, obj_bound);
		//obj_tree.buildIndex();
		//// ȷ��ǰ���ӿռ估ÿ���ӿռ��ڵĸ���
		//std::vector<size_t> front_inx;//ǰ�Ÿ�������
		//for (size_t i = 0; i < pop.size(); ++i) {
		//	if (pop[i].fitness() == 0) {
		//		front_inx.push_back(i);
		//	}
		//}
		////Ȼ��õ�ÿ���������ĸ��ӿռ�
		//std::vector<size_t> space_index;
		//for (size_t i = 0; i < pop.size(); ++i) {
		//	space_index.emplace_back(obj_tree.getRegionIdx(pop_obj[i]));
		//}
		////�õ��ӿռ�͸����ӳ��
		//std::map<size_t, std::vector<size_t>> space_indi;
		////���ݸ���λ�û����ӿռ�㼶
		//std::vector<std::vector<size_t>> subspace_layer;
		//std::vector<size_t> deal_flag(pop.size(), 0);
		///*for (size_t i = 0; i < space_index.size(); ++i) {
		//	space_indi[space_index[i]].push_back(i);
		//}*/
		//while (std::find(deal_flag.begin(), deal_flag.end(), 0) != deal_flag.end()) {
		//	//��ʣ������ǰ�ţ�ȷ���ӿؼ��㼶
		//	std::vector<size_t> residual_ind_inx;
		//	for (size_t i = 0; i < deal_flag.size(); ++i) {
		//		if (deal_flag[i] == 0) {
		//			residual_ind_inx.push_back(i);
		//		}
		//	}
		//	std::vector<std::vector<Real>> new_ind;
		//	for (size_t i = 0; i < residual_ind_inx.size(); ++i) {
		//		new_ind.emplace_back(pop_obj[residual_ind_inx[i]]);
		//	}
		//	//ʣ�����ǰ������
		//	auto f_inx = getNondominatedSetIndex(new_ind, CAST_CONOP(pro)->optimizeMode());
		//	std::vector<size_t> new_front_ind_inx;
		//	for (size_t i = 0; i < f_inx.size(); ++i) {
		//		new_front_ind_inx.push_back(residual_ind_inx[f_inx[i]]);
		//	}
		//	//ʣ������ǰ�������ӿռ�
		//	std::vector<size_t> new_front_space_inx;
		//	for (size_t i = 0; i < new_front_ind_inx.size(); ++i) {
		//		size_t inx = space_index[new_front_ind_inx[i]];
		//		if (std::find(new_front_space_inx.begin(), new_front_space_inx.end(), inx) == new_front_space_inx.end()) {
		//			new_front_space_inx.emplace_back(inx);
		//		}
		//	}
		//	//�µ��ӿռ��
		//	for (size_t i = 0; i < new_front_space_inx.size(); ++i) {
		//		std::vector<size_t> temp_space;
		//		for (size_t j = 0; j < space_index.size(); ++j) {
		//			if (deal_flag[j] == 0) {
		//				if (space_index[j] == new_front_space_inx[i]) {
		//					temp_space.push_back(j);
		//					deal_flag[j] = 1;
		//				}
		//			}
		//		}
		//		space_indi.insert(std::make_pair<>(new_front_space_inx[i], temp_space));
		//	}
		//	subspace_layer.emplace_back(new_front_space_inx);
		//}
		////��ǰ���ӿռ俪ʼ��ÿ���ӿռ䱣��

		// �������ӿռ������ѡ��
		std::vector<size_t> select_index = selectFromSpace(pop_obj, rank, select_num, pro, rnd);

		return select_index;
	}

	std::vector<size_t> SPMOEA::selectFromSpace(const std::vector<std::vector<Real>>& pop, std::vector<size_t> rank, size_t select_num, Problem *pro, Random *rnd) {
		// �ȹ�һ������
		auto pop_obj = pop;
		dataNormalize(pop_obj);
		// �ٹ����ӿռ�
		size_t M = CAST_CONOP(pro)->numberObjectives();
		size_t num_obj_spaces = numObjRegion();
		std::vector<Real> subspace_ratios(num_obj_spaces, 1. / num_obj_spaces);
		std::vector<std::pair<Real, Real>> obj_bound;
		for (size_t i = 0; i < M; ++i) {
			obj_bound.emplace_back(std::make_pair<>(0.,1.));
		}
		nanoflann::KDTreeSpace<Real> obj_tree(subspace_ratios,obj_bound);
		obj_tree.buildIndex();
		// ȷ��ǰ���ӿռ估ÿ���ӿռ��ڵĸ���
		std::vector<size_t> front_inx;//ǰ�Ÿ�������
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				front_inx.push_back(i);
			}
		}
		//Ȼ��õ�ÿ���������ĸ��ӿռ�
		std::vector<size_t> space_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			space_index.emplace_back(obj_tree.getRegionIdx(pop_obj[i]));
		}
		//�õ��ӿռ�͸����ӳ��
		std::map<size_t, std::vector<size_t>> space_indi;
		//���ݸ���λ�û����ӿռ�㼶
		std::vector<std::vector<size_t>> subspace_layer;
		std::vector<size_t> deal_flag(pop.size(), 0);
		/*for (size_t i = 0; i < space_index.size(); ++i) {
			space_indi[space_index[i]].push_back(i);
		}*/
		while (std::find(deal_flag.begin(), deal_flag.end(), 0) != deal_flag.end()) {
			//��ʣ������ǰ�ţ�ȷ���ӿؼ��㼶
			std::vector<size_t> residual_ind_inx;
			for (size_t i = 0; i < deal_flag.size(); ++i) {
				if (deal_flag[i] == 0) {
					residual_ind_inx.push_back(i);
				}
			}
			std::vector<std::vector<Real>> new_ind;
			for (size_t i = 0; i < residual_ind_inx.size(); ++i) {
				new_ind.emplace_back(pop_obj[residual_ind_inx[i]]);
			}
			//ʣ�����ǰ������
			auto f_inx = getNondominatedSetIndex(new_ind, CAST_CONOP(pro)->optimizeMode());
			std::vector<size_t> new_front_ind_inx;
			for (size_t i = 0; i < f_inx.size(); ++i) {
				new_front_ind_inx.push_back(residual_ind_inx[f_inx[i]]);
			}
			//ʣ������ǰ�������ӿռ�
			std::vector<size_t> new_front_space_inx;
			for (size_t i = 0; i < new_front_ind_inx.size(); ++i) {
				size_t inx = space_index[new_front_ind_inx[i]];
				if (std::find(new_front_space_inx.begin(), new_front_space_inx.end(), inx) == new_front_space_inx.end()) {
					new_front_space_inx.emplace_back(inx);
				}
			}
			//�µ��ӿռ��
			for (size_t i = 0; i < new_front_space_inx.size(); ++i) {
				std::vector<size_t> temp_space;
				for (size_t j = 0; j < space_index.size(); ++j) {
					if (deal_flag[j] == 0) {
						if (space_index[j] == new_front_space_inx[i]) {
							temp_space.push_back(j);
							deal_flag[j] = 1;
						}
					}
				}
				space_indi.insert(std::make_pair<>(new_front_space_inx[i], temp_space));
			}
			subspace_layer.emplace_back(new_front_space_inx);
		}
		//�ҳ�ÿ���ӿռ�����������������ĵ���ѡ��
		std::map<size_t,std::vector<std::vector<Real>>> space_dist_matrix;
		std::map<size_t, std::vector<size_t>> space_select_inx;//�ӿռ��������
		std::map<size_t, std::vector<size_t>> space_no_select_inx;//�ӿռ�δѡ��������
		std::map<size_t, Real> space_max_min_dist;//�ӿռ�����ѡ����С��������ֵ
		std::map<size_t, size_t> space_max_min_dist_inx;
		for (auto ss : space_indi) {
			std::vector<std::vector<Real>> temp_matrix;
			auto ind_inx = ss.second;
			std::vector<Real> min_dist;
			for (size_t i = 0; i < ind_inx.size(); ++i) {
				std::vector<Real> temp;
				auto p1 = pop_obj[ind_inx[i]];
				for (size_t j = 0; j < ind_inx.size(); ++j) {
					if (j == i) {
						temp.push_back(INT16_MAX);
					}
					else {
						auto p2 = pop_obj[ind_inx[j]];
						auto dist = euclideanDistance(p1.begin(), p1.end(), p2.begin());
						temp.push_back(dist);
					}
				}
				auto min_v = *std::min_element(temp.begin(),temp.end());
				min_dist.push_back(min_v);
				temp_matrix.emplace_back(temp);
			}
			auto max_v = *std::max_element(min_dist.begin(), min_dist.end());
			auto max_min_inx = std::distance(min_dist.begin(), std::max_element(min_dist.begin(),min_dist.end()));
			std::vector<size_t> temp_inx1;
			//temp_inx1.push_back(max_min_inx);
			std::vector<size_t> temp_inx2;
			for (size_t i = 0; i < ss.second.size(); ++i) {
				temp_inx2.push_back(i);
			}
			space_max_min_dist.insert(std::make_pair<>(ss.first,max_v));
			space_max_min_dist_inx.insert(std::make_pair<>(ss.first, max_min_inx));
			space_select_inx.insert(std::make_pair<>(ss.first, temp_inx1));
			space_no_select_inx.insert(std::make_pair<>(ss.first, temp_inx2));
			space_dist_matrix.insert(std::make_pair<>(ss.first,temp_matrix));
		}
		//Ȼ�����αȽ�ÿ���ӿռ�������ѡ����������ӿռ䣬������ѡ��
		Real ratio = 0.2;
		size_t total_size = 0;
		while (total_size < select_num) {
			//�Ƚϸ����ӿռ�����ѡ�����С��������ֵ
			std::vector<Real> space_max_dist;
			for (size_t i = 0; i < subspace_layer.size(); ++i) {
				for (size_t j = 0; j < subspace_layer[i].size(); ++j) {
					size_t space_inx = subspace_layer[i][j];
					auto &temp_matrix = space_dist_matrix[space_inx];
					auto &no_sele = space_no_select_inx[space_inx];
					auto &sele= space_select_inx[space_inx];
					if (no_sele.size() == 0) {
						space_max_dist.push_back(0.);
						space_max_min_dist[space_inx] = 0.;
						space_max_min_dist_inx[space_inx] = 0;
					}
					else if (sele.size() == 0) {//�����ӿռ���С�������ֵ
						space_max_dist.push_back(space_max_min_dist[space_inx]/(1+i*ratio));
					}
					else {
						std::vector<std::vector<Real>> temp_mat;
						for (size_t p = 0; p < no_sele.size(); ++p) {
							std::vector<Real> tt1;
							for (size_t q = 0; q < sele.size(); ++q) {
								tt1.push_back(temp_matrix[no_sele[p]][sele[q]]);
							}
							temp_mat.emplace_back(tt1);
						}
						std::vector<Real> min_dist;
						for (size_t k = 0; k < temp_mat.size(); ++k) {
							auto min_v = *std::min_element(temp_mat[k].begin(), temp_mat[k].end());
							min_dist.push_back(min_v);
						}
						auto max_min= *std::max_element(min_dist.begin(), min_dist.end());
						auto max_inx = std::distance(min_dist.begin(), std::max_element(min_dist.begin(), min_dist.end()));
						space_max_dist.push_back(max_min/(1+i*ratio));
						space_max_min_dist[space_inx] = max_min/(1+i*ratio);
						space_max_min_dist_inx[space_inx] = max_inx;
					}
				}
			}
			Real max_v = 0.;
			size_t max_space_inx;
			for (auto ss : space_max_min_dist) {
				if (ss.second > max_v) {
					max_space_inx = ss.first;
					max_v = ss.second;
				}
			}
			//����ѡ��״̬
			size_t index = space_max_min_dist_inx[max_space_inx];
			space_select_inx[max_space_inx].push_back(space_no_select_inx[max_space_inx][index]);
			space_no_select_inx[max_space_inx].erase(space_no_select_inx[max_space_inx].begin()+index);
			total_size++;
		}
		// �������ӿռ������ѡ��
		std::vector<size_t> select_index;
		for (auto ss : space_select_inx) {
			for (size_t i = 0; i < ss.second.size(); ++i) {
				select_index.push_back(space_indi[ss.first][ss.second[i]]);
			}
		}
		return select_index;
	}

	std::vector<size_t> SPMOEA::selectIndiFromSpace(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t select_num, Problem *pro,Random *rnd) {
		//��ͳ��ǰ�Ÿ��������
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		auto pop_att = popAttach(pop,m_num_region_obj,pro);
		std::vector<size_t> select_index;//�Ѿ�ѡ����������
		std::vector<size_t> behind_space_index;//�����ӿռ�����
		std::vector<size_t> front_space_index;//ǰ���ӿռ�����
		size_t front_space_num = 0;//ǰ���ӿռ��ڵĸ�������
		std::vector<size_t> first_ind_index;//�ȵõ���һ�Ÿ��������
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i]->fitness() == 0)
				first_ind_index.emplace_back(i);
		}
		for (const auto& ind : pop_att[0][0]) {
			front_space_index.emplace_back(ind.first);
			front_space_num += ind.second.size();
		}
		for (const auto& ind : pop_att[0][1]) {
			behind_space_index.emplace_back(ind.first);
		}
		std::map<size_t, std::vector<size_t>> map_first_indi;//�õ�ǰ�Ÿ���ֲ�����Щ�ӿռ�
		for (const auto& sp : pop_att[0][0]) {
			std::pair<size_t, std::vector<size_t>> temp;
			temp.first = sp.first;
			for (const auto& in : sp.second) {
				if (pop[in]->fitness() == 0) {
					temp.second.push_back(in);
				}
			}
			map_first_indi.insert(temp);
		}
		std::map<size_t, std::vector<size_t>> map_selected_indi;//ǰ���ӿռ��Ѿ�ѡ�Ľ⣬˳������
		std::map<size_t, std::vector<size_t>> map_no_selected_indi;//ǰ���ӿռ��Ѿ�ѡ�Ľ⣬˳������
		//����ǰ���ӿռ��ڲ��ĸ����������и���ѡ��ѡ��ʱ���Ǹ�������
		if (first_ind_index.size() > select_num) {//ǰ�ظ����Ƿ����
			select_index = selectIndiFromFront(pop, select_num, pro);
		}
		else if (front_space_num > select_num) {//ǰ���ӿռ��ڵĸ������Ƿ����
			bool flag = 1;
			if (flag) {//Ŀ��ռ�ѡ��
				//��ѡǰ�ظ���
				select_index = first_ind_index;
				//������ѡ��ǰ���ӿռ��е���������
				auto front_subspace_ind = pop_att[0][0];
				for (const auto& i : front_subspace_ind) {
					std::pair<size_t, std::vector<size_t>> temp;
					temp.first = i.first;
					std::vector<size_t> f_ind_idx;
					for (size_t j = 0; j < i.second.size(); ++j) {
						if (find(first_ind_index.begin(), first_ind_index.end(), i.second[j]) != first_ind_index.end()) {
							f_ind_idx.push_back(i.second[j]);
						}
					}
					temp.second = f_ind_idx;
					map_selected_indi.insert(temp);//ǰ���ӿռ��Ѿ�ѡ�Ľ�
				}
				while (select_index.size() < select_num) {
					for (const auto& i : front_subspace_ind) {
						if (map_selected_indi[i.first].size() < i.second.size()) {
							//��ʣ�������ѡ��һ��
							for (size_t j = 0; j < i.second.size(); ++j) {
								if (find(map_selected_indi[i.first].begin(), map_selected_indi[i.first].end(), i.second[j]) == map_selected_indi[i.first].end()) {
									select_index.push_back(i.second[j]);
									map_selected_indi[i.first].push_back(i.second[j]);
									break;
								}
							}
						}
						if (select_index.size() >= select_num) {
							break;
						}
					}
				}
			}
			else {//��ռ�ѡ�񣬸��ݽ�ռ�ֲ�ѡ�񣬻��߸����ӿռ��Ǳ��ֵѡ��
				//��ǰ���ӿռ�ĵ�ֳ������֣�ǰ�ؽ�ͷ�ǰ�ؽ�

				//�����н�ӳ�䵽�����ռ䣬�ҳ���ǰ�ؽ���Զ�ĸ���
			}
		}
		else {
			//�ȼ���ǰ���ӿռ���������
			for (const auto& ind : pop_att[0][0]) {
				for (const auto& p : ind.second) {
					select_index.push_back(p);
				}
			}
			//���ں����ӿռ�����ѡ��
			std::vector<size_t> space_select_num(pop_att[0][1].size(), 0);//ÿ�������ӿռ���ѡ��ĸ���
			std::vector<std::vector<size_t>> behind_spaces;
			for (const auto& sp : pop_att[0][1]) {
				behind_spaces.push_back(sp.second);
			}
			while (select_index.size() < select_num) {
				for (size_t i = 0; i < behind_spaces.size(); ++i) {
					if (space_select_num[i] < behind_spaces[i].size()) {
						size_t select_idx = std::floor(behind_spaces[i].size() * rnd->uniform.next());
						if (std::find(select_index.begin(), select_index.end(), behind_spaces[i][select_idx]) == select_index.end()) {
							select_index.push_back(behind_spaces[i][select_idx]);
							space_select_num[i]++;
						}
					}
					if (select_index.size() >= select_num) {
						break;
					}
				}
			}
		}
		return select_index;
	}

	space_attach SPMOEA::spaceAttach(const Population<Solution<>>& pop,size_t num_spaces,Problem *pro, bool b) {
		space_attach space_all_indi;
		//�ȵõ�pop�ķ�Χ������pop��Χ���ֿռ���
		size_t dims;
		if (b) {
			dims = CAST_CONOP(pro)->numberObjectives();
		}
		else {
			dims = CAST_CONOP(pro)->numberVariables();
		}
		std::vector<std::pair<Real, Real>> bound_range(dims);
		std::vector<std::pair<Real, Real>> m_boundary;
		if (b) {
			for (int i = 0; i < dims; ++i) {
				bound_range[i].first = 1.0e14;
				for (int j = 0; j < pop.size(); ++j) {
					if (pop[j].objective()[i] < bound_range[i].first) {
						bound_range[i].first = pop[j].objective()[i];
					}
				}
			}
			for (int i = 0; i < dims; ++i) {
				bound_range[i].second = -1 * 1.0e14;
				for (int j = 0; j < pop.size(); ++j) {
					if (pop[j].objective()[i] > bound_range[i].second) {
						bound_range[i].second = pop[j].objective()[i];
					}
				}
			}
			if (m_normalize) {
				for (size_t i = 0; i < dims; ++i) {
					m_boundary.push_back(std::make_pair<>(0., 1.));
				}
			}
			else {
				m_boundary = bound_range;
			}
		}
		else {
			for (size_t i = 0; i < dims; ++i) {
				m_boundary.emplace_back(CAST_CONOP(pro)->domain().range(i).limit);
			}
		}
		//���ݱ߽緶Χ�ͷ��������ӿռ�
		if (b) {//Ŀ��ռ�
			m_mo_hlc->updateObjTree();
			m_mo_hlc->initialObjSpace(m_boundary, num_spaces);
		}
		else {
			m_mo_hlc->updateVarSelectionTree();
			m_mo_hlc->initialVarSelectionSpace(m_boundary, num_spaces);
		}
		//�õ���һ�Ÿ��������
		std::vector<size_t> first_ind_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].fitness() == 0)
				first_ind_index.emplace_back(i);
		}
		//Ȼ��õ�ÿ���������ĸ��ӿռ�
		std::vector<size_t> space_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			std::vector<Real> temp;
			if (b) {
				temp = pop[i].objective();
				if (m_normalize) {
					for (size_t j = 0; j < temp.size(); ++j) {
						temp[j] = (temp[j] - bound_range[j].first) / (bound_range[j].second - bound_range[j].first);
					}
				}
				space_index.emplace_back(m_mo_hlc->getObjspaceTree().getRegionIdx(temp));
			}
			else {
				temp = pop[i].variable().vect();
				space_index.emplace_back(m_mo_hlc->getVarSelectionTree().getRegionIdx(temp));
			}
		}
		//�õ�ǰ���ӿռ�����
		std::vector<size_t> nondominate_space;
		for (size_t i = 0; i < first_ind_index.size(); ++i) {
			size_t temp1 = space_index[first_ind_index[i]];
			if (nondominate_space.empty())
				nondominate_space.emplace_back(temp1);
			else if (std::find(nondominate_space.begin(), nondominate_space.end(), temp1) == nondominate_space.end()) {
				nondominate_space.emplace_back(temp1);
			}
		}
		//ǰ���ӿռ��ڵĸ�������
		size_t front_indi_num = 0;
		//�õ�ǰ���ӿռ��ڰ��������и���
		std::map<size_t, std::vector<size_t>> front_indi;//keyΪ�ӿռ�������valueΪ��������
		std::vector<size_t> redi_space_index = space_index;
		for (size_t i = 0; i < nondominate_space.size(); ++i) {
			std::pair<size_t, std::vector<size_t>> temp;
			temp.first = nondominate_space[i];
			for (size_t j = 0; j < space_index.size(); ++j) {
				if (space_index[j] == temp.first) {
					temp.second.emplace_back(j);
					++front_indi_num;
				}
			}
			front_indi.insert(temp);
		}
		space_all_indi.emplace_back(front_indi);

		//��ǰ���ӿռ�����
		std::vector<size_t> dominated_space;
		//�õ���ǰ���ӿռ��ڰ����ĸ���
		//std::vector<std::pair<size_t, std::vector<size_t>>> behind_indi;
		std::map<size_t, std::vector<size_t>> behind_indi;
		if (front_indi_num < pop.size()) {
			for (size_t i = 0; i < space_index.size(); ++i) {
				size_t temp = space_index[i];
				if (std::find(nondominate_space.begin(), nondominate_space.end(), temp) == nondominate_space.end()) {
					if (dominated_space.empty())
						dominated_space.emplace_back(temp);
					else if (std::find(dominated_space.begin(), dominated_space.end(), temp) == dominated_space.end()) {
						dominated_space.emplace_back(temp);
					}
				}
			}
			for (size_t i = 0; i < dominated_space.size(); ++i) {
				std::pair<size_t, std::vector<size_t>> temp;
				temp.first = dominated_space[i];
				for (size_t j = 0; j < space_index.size(); ++j) {
					if (space_index[j] == temp.first)
						temp.second.emplace_back(j);
				}
				behind_indi.insert(temp);
			}
		}
		space_all_indi.emplace_back(behind_indi);
		return space_all_indi;
	}

	space_attach SPMOEA::spaceAttach(const std::vector<std::shared_ptr<Solution<>>>& pop,size_t num_spaces, Problem *pro,bool b) {
		space_attach space_all_indi;
		//�ȵõ�pop�ķ�Χ������pop��Χ���ֿռ���
		size_t dims;
		if (b) {
			dims = CAST_CONOP(pro)->numberObjectives();
		}
		else {
			dims = CAST_CONOP(pro)->numberVariables();
		}
		std::vector<std::pair<Real, Real>> bound_range(dims);
		std::vector<std::pair<Real, Real>> m_boundary;
		if (b) {
			for (int i = 0; i < dims; ++i) {
				bound_range[i].first = 1.0e14;
				for (int j = 0; j < pop.size(); ++j) {
					if (pop[j]->objective()[i] < bound_range[i].first) {
						bound_range[i].first = pop[j]->objective()[i];
					}
				}
			}
			for (int i = 0; i < dims; ++i) {
				bound_range[i].second = -1 * 1.0e14;
				for (int j = 0; j < pop.size(); ++j) {
					if (pop[j]->objective()[i] > bound_range[i].second) {
						bound_range[i].second = pop[j]->objective()[i];
					}
				}
			}
			if (m_normalize) {
				for (size_t i = 0; i < dims; ++i) {
					m_boundary.push_back(std::make_pair<>(0., 1.));
				}
			}
			else {
				m_boundary = bound_range;
			}
		}
		else {
			for (size_t i = 0; i < dims; ++i) {
				m_boundary.emplace_back(CAST_CONOP(pro)->domain().range(i).limit);
			}
		}
		//���ݱ߽緶Χ�ͷ��������ӿռ�
		if (b) {//Ŀ��ռ�
			m_mo_hlc->updateObjTree();
			m_mo_hlc->initialObjSpace(m_boundary, num_spaces);
		}
		else {
			m_mo_hlc->updateVarSelectionTree();
			m_mo_hlc->initialVarSelectionSpace(m_boundary, num_spaces);
		}
		//�õ���һ�Ÿ��������
		std::vector<size_t> first_ind_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i]->fitness() == 0)
				first_ind_index.emplace_back(i);
		}
		//Ȼ��õ�ÿ���������ĸ��ӿռ�
		std::vector<size_t> space_index;
		for (size_t i = 0; i < pop.size(); ++i) {
			std::vector<Real> temp;
			if (b) {
				temp = pop[i]->objective();
				if (m_normalize) {
					for (size_t j = 0; j < temp.size(); ++j) {
						temp[j] = (temp[j] - bound_range[j].first) / (bound_range[j].second - bound_range[j].first);
					}
				}
				space_index.emplace_back(m_mo_hlc->getObjspaceTree().getRegionIdx(temp));
			}
			else {
				temp = pop[i]->variable().vect();
				space_index.emplace_back(m_mo_hlc->getVarSelectionTree().getRegionIdx(temp));
			}
		}
		//�õ�ǰ���ӿռ�����
		std::vector<size_t> nondominate_space;
		for (size_t i = 0; i < first_ind_index.size(); ++i) {
			size_t temp1 = space_index[first_ind_index[i]];
			if (nondominate_space.empty())
				nondominate_space.emplace_back(temp1);
			else if (std::find(nondominate_space.begin(), nondominate_space.end(), temp1) == nondominate_space.end()) {
				nondominate_space.emplace_back(temp1);
			}
		}
		//ǰ���ӿռ��ڵĸ�������
		size_t front_indi_num = 0;
		//�õ�ǰ���ӿռ��ڰ��������и���
		std::map<size_t, std::vector<size_t>> front_indi;//keyΪ�ӿռ�������valueΪ��������
		std::vector<size_t> redi_space_index = space_index;
		for (size_t i = 0; i < nondominate_space.size(); ++i) {
			std::pair<size_t, std::vector<size_t>> temp;
			temp.first = nondominate_space[i];
			for (size_t j = 0; j < space_index.size(); ++j) {
				if (space_index[j] == temp.first) {
					temp.second.emplace_back(j);
					++front_indi_num;
				}
			}
			front_indi.insert(temp);
		}
		space_all_indi.emplace_back(front_indi);

		//��ǰ���ӿռ�����
		std::vector<size_t> dominated_space;
		//�õ���ǰ���ӿռ��ڰ����ĸ���
		//std::vector<std::pair<size_t, std::vector<size_t>>> behind_indi;
		std::map<size_t, std::vector<size_t>> behind_indi;
		if (front_indi_num < pop.size()) {
			for (size_t i = 0; i < space_index.size(); ++i) {
				size_t temp = space_index[i];
				if (std::find(nondominate_space.begin(), nondominate_space.end(), temp) == nondominate_space.end()) {
					if (dominated_space.empty())
						dominated_space.emplace_back(temp);
					else if (std::find(dominated_space.begin(), dominated_space.end(), temp) == dominated_space.end()) {
						dominated_space.emplace_back(temp);
					}
				}
			}
			for (size_t i = 0; i < dominated_space.size(); ++i) {
				std::pair<size_t, std::vector<size_t>> temp;
				temp.first = dominated_space[i];
				for (size_t j = 0; j < space_index.size(); ++j) {
					if (space_index[j] == temp.first)
						temp.second.emplace_back(j);
				}
				behind_indi.insert(temp);
			}
		}
		space_all_indi.emplace_back(behind_indi);
		return space_all_indi;
	}

	pop_attach SPMOEA::popAttach(const Population<Solution<>>& pop,size_t num_spaces,Problem *pro) {
		pop_attach all_indi;
		all_indi.emplace_back(spaceAttach(pop, num_spaces, pro, true));
		all_indi.emplace_back(spaceAttach(pop, num_spaces, pro, false));
		return all_indi;
	}

	pop_attach SPMOEA::popAttach(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t num_spaces, Problem *pro) {
		pop_attach all_indi;
		all_indi.emplace_back(spaceAttach(pop, num_spaces, pro, true));
		all_indi.emplace_back(spaceAttach(pop, num_spaces, pro, false));
		return all_indi;
	}

	int SPMOEA_pop::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		return 0;
	}

	void SPMOEA_pop::NDSort(Population<Solution<>>& pop, Problem* pro) {
		//nondominatedSorting(pop);
		popNondominatedSorting(pop,pro->optimizeMode());
	}

	//�ȹ�һ��Ŀ��ֵ
	void SPMOEA_pop::evalDense(Population<Solution<>>& pop,Population<Solution<>>& pop_combined, Problem* pro) {
		int pops = 0;  //indicate parent population size be 0
		int size = pop_combined.size();
		int rank = 0;
		while (true) {
			int count = 0;
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
					count++;
			int size2 = pops + count;
			if (size2 > pop.size()) {
				break;
			}
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
				{
					pop[pops] = pop_combined[i];
					++pops;
				}
			rank++;
			if (pops >= pop.size()) break;
		}
		if (pops < pop.size()) {
			std::vector<int> list;
			// save the Solutions in the overflowed front
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
					list.push_back(i);
			int s2 = list.size();
			std::vector<Real> density(s2);
			std::vector<Real> obj(s2);
			std::vector<int> idx(s2);
			std::vector<int> idd(s2);
			for (size_t i = 0; i < s2; i++) {
				idx[i] = i;
				density[i] = 0;
			}
			for (size_t j = 0; j < pro->numberObjectives(); j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop_combined[list[i]].objective()[j];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = pop.size() - pops;
			mergeSort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				pop[pops] = pop_combined[list[idx[i]]];
				++pops;
			}
			density.clear();
			idx.clear();
			list.clear();
		}
	}

	std::vector<size_t> SPMOEA_pop::denseSelect(Population<Solution<>>& pop, Population<Solution<>>& pop_combined, Problem* pro) {
		std::vector<size_t> sele_inx;
		int pops = 0;  //indicate parent population size be 0
		int size = pop_combined.size();
		int rank = 0;
		while (true) {
			int count = 0;
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
					count++;
			int size2 = pops + count;
			if (size2 > pop.size()) {
				break;
			}
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
				{
					pop[pops] = pop_combined[i];
					++pops;
					sele_inx.push_back(i);
				}
			rank++;
			if (pops >= pop.size()) break;
		}
		if (pops < pop.size()) {
			std::vector<int> list;
			// save the Solutions in the overflowed front
			for (size_t i = 0; i < size; i++)
				if (pop_combined[i].fitness() == rank)
					list.push_back(i);
			int s2 = list.size();
			std::vector<Real> density(s2);
			std::vector<Real> obj(s2);
			std::vector<int> idx(s2);
			std::vector<int> idd(s2);
			for (size_t i = 0; i < s2; i++) {
				idx[i] = i;
				density[i] = 0;
			}
			for (size_t j = 0; j < pro->numberObjectives(); j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop_combined[list[i]].objective()[j];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = pop.size() - pops;
			mergeSort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				pop[pops] = pop_combined[list[idx[i]]];
				++pops;
				sele_inx.push_back(list[idx[i]]);
			}
			density.clear();
			idx.clear();
			list.clear();
		}
		return sele_inx;
	}

	void SPMOEA::NDSort(Population<Solution<>>& pop) {
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i) {
			objs.emplace_back(&pop[i].objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < pop.size(); ++i) {
			pop[i].setFitness(rank[i]);
		}
	}

	size_t SPMOEA::layerNDSort(Population<Solution<>>& pop) {
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i) {
			objs.emplace_back(&pop[i].objective());
		}
		std::vector<int> rank;
		size_t layer_num=ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < pop.size(); ++i) {
			pop[i].setFitness(rank[i]);
		}
		return layer_num;
	}

	void SPMOEA::visitCount(const std::vector<std::vector<Real>>& sols, std::vector<Real>& middle_point) {
		/*for (size_t i = 0; i < sols.size(); ++i) {
			size_t parseBinary = 0;
			for (size_t j = 0; j < sols[i].size(); ++j) {
				if (sols[i][j] > middle_point[j]) {
					parseBinary += std::pow(2, sols[i].size() - 1 - j);
				}
			}
			m_visit_count[parseBinary]++;
		}*/
	}

	void SPMOEA::recordMetrics(Problem *pro, Algorithm *alg) {
		/************************************/
		/*            ����ָ�����          */
		/************************************/
		//ʹ��archive��������ָ��
		Population<Solution<>> temp_pop;
		/*for (size_t i = 0; i < m_archive.size(); ++i) {
			temp_pop.append(*m_archive[i]);
		}*/
		for (size_t i = 0; i < getPop()[0].size(); ++i) {
			temp_pop.append(getPop()[0][i]);
		}
		/*for (size_t i = 0; i < m_history_front_sols.size(); ++i) {
			temp_pop.append(*m_history_front_sols[i]);
		}*/
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < temp_pop.size(); ++i) {
			pop_objs.push_back(temp_pop[i].objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);
		m_IGD.push_back(temp_IGD);
		//m_IGD.push_back(CAST_CONOP(pro)->optima()->invertGenDist(temp_pop));
		std::cout << alg->evaluations() << "  " << m_IGD.back() << std::endl;
		if (m_IGD.size() > 1) {
			auto delta_igd = std::fabs(m_IGD.back() - m_IGD[m_IGD.size() - 2]);
			auto ratio=delta_igd / m_IGD[m_IGD.size() - 2];
			std::cout << "delta IGD ratio:     " << ratio << std::endl;
		}
		//record();//store metrics data
		/*for (size_t i = 0; i < m_visit_count.size(); ++i) {
			if (m_visit_count[i] > 0) {
				std::cout << (Real)m_visit_count[i]/ alg->evaluations() << "    ";
			}
		}*/
		std::cout <<"��ʷ������� "<< m_history_front_sols.size() << std::endl;
		std::cout << std::endl;
	}

	//�������ݵĵ����ڶ���Ϊ1�����һ������ʵֵ
	std::vector<Real> SPMOEA::linearRegression(const std::vector<std::vector<Real>>& input_data, Real learn_ratio, int iterate_times) {
		std::vector<std::vector<Real>> training_set = input_data;
		std::vector<Real> weigh(input_data[0].size() - 1, 1.0);
		//�ݶ��½�
		gradient_descent(training_set, input_data[0].size() - 1, weigh, learn_ratio, iterate_times);//�ڶ�������Ϊ������
		return weigh;
	}

	// �ݶ��½�
	void SPMOEA::gradient_descent(const std::vector<std::vector<Real>>& training_set, int feature_num, std::vector<Real>& w, Real learn_ratio, int iterator_time) {
		int sample_num = training_set.size();
		size_t count = 0;
		//Real loss = 1.;
		while (count<iterator_time) {
			//����ǰJ��ֵ
			Real temp = Theta(training_set, feature_num, w);

			std::vector<Real> del_theta(feature_num, 0.);
			for (int i = 0; i < feature_num; i++) {
				for (int j = 0; j < sample_num; j++) {
					del_theta[i] += ((predict(w, training_set[j], feature_num) - training_set[j][feature_num]) * training_set[j][i]);
				}
			}
			//w[i]�ĸ��±�������е�del_theta��������˲ſ��ԣ���Ȼ���µĻ�Ӱ��û���µ�
			//���������ڴ����ڱ�ʾ���������forѭ�����ܺ�����ĺϲ���
			for (int i = 0; i < feature_num; i++)
			{
				w[i] -= learn_ratio * del_theta[i] / (Real)sample_num;
			}
			//������J��ֵ
			Real temp1 = Theta(training_set, feature_num, w);
			//loss = temp1;
			//���ε���J��ֵ�仯С��0.001��ѭ����ֹ
			if (std::fabs(temp1) < 0.000001) {
				break;
			}
			/*if (std::fabs(temp1-temp) < 0.0001) {
				break;
			}*/
			count++;
			//std::cout <<"iterate_time: "<< count<<"         " << "J_Theta = " << Theta(training_set, feature_num, w) << std::endl;
		}
	}

	// ��ʧ����
	Real SPMOEA::Theta(const std::vector<std::vector<Real>>& training_set, int featue_num, const std::vector<Real>& w) {
		Real sum = 0;
		for (int i = 0; i < training_set.size(); i++) {
			sum += (training_set[i][featue_num] - predict(w, training_set[i], featue_num)) * (training_set[i][featue_num] - predict(w, training_set[i], featue_num));
		}
		return sum / (2 * training_set.size());
	}

	Real SPMOEA::predict(const std::vector<Real>& weigh, const std::vector<Real>& data, int feature_num) {
		Real sum = 0;
		for (int i = 0; i < feature_num; i++) {
			sum += weigh[i] * data[i];
		}
		return sum;
	}

	std::pair<int, Real> SPMOEA::findSplitDimAndPos(int inx, Problem* pro) {
		//��������һά�ϻ����ӿռ�ǰ�ؽ�ʹ�û��ֺ��ӿռ��ķֲ��������Ի�
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_front_sol_in_subspace;
		auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> front_ind_ratio;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			Real middle_v = (space_bound[j].first + space_bound[j].second) / 2;
			size_t upper_count = 0, lower_count = 0;
			for (size_t k = 0; k < front_sols.size(); ++k) {
				if (front_sols[k]->variable()[j] >= middle_v) {
					upper_count++;
				}
				else {
					lower_count++;
				}
			}
			size_t max_v = upper_count > lower_count ? upper_count : lower_count;
			size_t min_v = upper_count < lower_count ? upper_count : lower_count;
			if (max_v == 0) {
				max_v += 1;
			}
			front_ind_ratio.push_back((Real)min_v / max_v);
		}
		int dim = std::distance(front_ind_ratio.begin(), std::min_element(front_ind_ratio.begin(), front_ind_ratio.end()));

		//ʹ�ÿ��ȷ���ָ��ά��
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			dim_span.push_back(space_bound[j].second - space_bound[j].first);
		}
		dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));

		Real pos = space_bound[dim].first + (space_bound[dim].second - space_bound[dim].first) / 2;

		return std::make_pair<>(dim, pos);
	}
}
