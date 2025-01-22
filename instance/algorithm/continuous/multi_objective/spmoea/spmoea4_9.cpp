#include "spmoea4_9.h"
#include "../../../../../utility/clustering/dbscan.h"
#include "../../../../../utility/clustering/kmeans.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA4_9::initialize_() {
		SPMOEA::initialize_();
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
		auto& v = *m_param;

		m_greedy_rate = v.get<Real>("greedy rate");
		m_disturbance_rate = v.get<Real>("disturbance rate");
		//m_neigh_interactive_rate = v.get<Real>("neigh interactive rate");
		m_subspace_interactive_rate = (1 - m_disturbance_rate) / 2;
		m_neigh_interactive_rate = 1 - m_disturbance_rate - m_subspace_interactive_rate;
		//m_switch_period = m_num_push + m_num_extend;

		//m_neighs = 5;
	}

	void SPMOEA4_9::run_() {
		initPop(m_problem.get(), this, m_random.get());
		//#ifdef OFEC_DEMO
		//		updateBuffer();
		//#endif
		while (!terminating()) {
			evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void SPMOEA4_9::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA4_9::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(getPop().size());
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_solution[i].push_back(&getPop()[i][j]);
				}
			}
			m_off_solution.clear();
			m_off_solution.resize(getPop().size());
			for (size_t i = 0; i < getPop().size(); ++i) {
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					m_off_solution[i].push_back(&getPop()[i].getOffspring()[j]);
				}
			}
			////��������Ľ�
			//m_off_solution.resize(m_pop_index_clusters.size());
			//for (size_t i = 0; i < m_pop_index_clusters.size(); ++i) {
			//	for (size_t j = 0; j < m_pop_index_clusters[i].size(); ++j) {
			//		if (m_pop_index_clusters[i][j] < getPop()[0].size()) {
			//			m_off_solution[i].push_back(&getPop()[0].getOffspring()[m_pop_index_clusters[i][j]].phenotype());
			//		}
			//		else {
			//			m_off_solution[i].push_back(&getPop()[0][m_pop_index_clusters[i][j]-getPop()[0].size()].phenotype());
			//		}
			//	}
			//}
			m_his_solution.clear();
			auto& his_sols = getHisFrontSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_his_solution.push_back(his_sols[i].get());
			}

			/*for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				auto& subspace_front_sols = getMO_HLC().getSubspaceInfo(i).m_subspace_front_sol;
				for (size_t j = 0; j < subspace_front_sols.size(); ++j) {
					m_solution.back().push_back(&subspace_front_sols[j]->phenotype());
				}
			}*/

			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA4_9::initiVarSpace(Problem* pro) {
		/*************************************/
		/*      MO_HLC�����ռ仮�ֳ�ʼ��     */
		/*************************************/
		size_t var_num = CAST_CONOP(pro)->numberVariables();
		std::vector<std::pair<Real, Real>> m_var_boundary;
		for (size_t i = 0; i < var_num; ++i) {
			m_var_boundary.emplace_back(CAST_CONOP(pro)->range(i));
		}
		getMO_HLC().initialVarSpace(m_var_boundary, numVarRegion());
		m_pop_var_range = m_var_boundary;
		m_num_e_e.resize(2);

	}

	void SPMOEA4_9::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//�ڸ����ӿռ�����һ������Ⱥ
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
		size_t pop_num = getPopsize();
		//size_t space_num = getMO_HLC().numSubspace();//��ʼ���ӿռ�����
		Population<Solution<>> new_pop;
		//ȫ����Ⱥ
		SPMOEA_pop temp_pop(pop_num, pro);
		temp_pop.initialize(pro, rnd);
		std::vector<std::vector<Real>> sols;
		auto bound = CAST_CONOP(pro)->boundary();
		auto num_var = CAST_CONOP(pro)->numberVariables();
		for (size_t j = 0; j < pop_num; ++j) {
			temp_pop[j].evaluate(pro, alg);
			sols.emplace_back(temp_pop[j].variable().vect());
		}
		SPMOEA::NDSort(temp_pop);
		//��ʼ���Ӵ�
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			temp_pop.getOffspring()[j] = temp_pop[j];
			temp_pop.getOffspring()[temp_pop.size() + j] = temp_pop[j];
			new_pop.append(temp_pop[j]);
		}

		temp_pop.setPopState("active");
		//��ʼ������Ⱥ����ʷ��ÿ��ǰ�ؽ�
		Population<Solution<>> front_pop;
		for (size_t j = 0; j < temp_pop.size(); ++j) {
			if (temp_pop[j].fitness() == 0) {
				front_pop.append(temp_pop[j]);
				temp_pop.getPopHisFrontSols().emplace_back(std::make_shared<Solution<>>(temp_pop[j]));
			}
		}
		temp_pop.getPopGenFrontSols().emplace_back(std::make_shared<Population<Solution<>>>(front_pop));
		getPop().append(temp_pop);
		getPop().back().setSearchBox(bound);

		//������ʷ��Ϣ
		SPMOEA::NDSort(new_pop);
		updateHistoryInfo(new_pop, pro);
		SPMOEA::updateObjRange(new_pop, pro);
		SPMOEA::updateNewPop(new_pop);
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		//�����ӿռ���Ϣ
		for (size_t i = 0; i < getPop().size(); ++i) {
			Population<Solution<>> temp_off_pop;
			for (size_t j = 0; j < getPop()[i].getOffspring().size() - getPop()[i].size(); ++j) {
				temp_off_pop.append(getPop()[i].getOffspring()[j]);
			}
			SPMOEA::updateSubspaceFrontSol(temp_off_pop, pro, rnd);
		}
		updateFrontSpace();
		updateVarSpaceRank(pro, rnd);
		updateEE(pop_num, 0);
		SPMOEA::recordMetrics(pro, alg);

		updateFrontRegionLinkSpace(pro, rnd);
		//�����ӿռ仮�־���
		clusterSubspace();
		setAddNumSpace(0);
		Real total_volume = getVarSpaceVolume();
		Real front_ratio = 0.;
		size_t front_count = 0;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < getFrontSpace().size(); ++i) {
			auto v = getMO_HLC().subspaceTree().getBoxVolume(getFrontSpace()[i]);
			front_ratio += (v / total_volume);
		}
		setFrontSpaceRatio(front_ratio);
		//m_divide_granularity = getSplitGranularity();
		m_divide_granularity = 15 - 5. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
		//������ʷ�ݻ��켣
		for (size_t j = 0; j < getPop()[0].size(); ++j) {
			std::vector<std::shared_ptr<Solution<>>> temp_ind;
			temp_ind.emplace_back(std::make_shared<Solution<>>(getPop()[0][j]));
			getHisEvolveLocus().emplace_back(temp_ind);
		}
		//���ӱ���
		std::vector<Real> operator_ratio(2, 1.);
		operator_ratio[1] = 0.;
		getOperatorRatio().emplace_back(operator_ratio);
		//�������ڵ��ӿռ�
		std::vector<size_t> ind_space;
		for (size_t j = 0; j < getPop()[0].size(); ++j) {
			auto sol = getPop()[0][j].variable().vect();
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			if (std::find(ind_space.begin(), ind_space.end(), space_inx) == ind_space.end()) {
				ind_space.push_back(space_inx);
			}
		}
		std::vector<std::vector<size_t>> plot_spaces;
		auto front_clusters = getFrontRegionLinkSpace();
		auto front_spaces = getFrontSpace();
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

		//��ʵPS���ӿռ�֮��Ĺ�ϵ
		//���ŵ�һάֱ�ӹ����
		size_t num = 10000;
		//auto bound = CAST_CONOP(pro)->boundary();
		Real delta = (bound[0].second - 0.001 - bound[0].first) / num;
		std::vector<size_t> real_front;
		std::vector<size_t> multiseg_front;
		size_t temp_space = INT16_MAX;
		for (size_t j = 0; j < num - 1; ++j) {
			std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), bound[0].first + 0.0001 + (j + 1) * delta);
			auto sol = CAST_CONOP(pro)->createVar(s);
			//�ж�sol���ĸ��ӿռ�
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			if (space_inx != temp_space) {
				if (std::find(real_front.begin(), real_front.end(), space_inx) == real_front.end()) {
					real_front.push_back(space_inx);
				}
				else {
					if (space_inx != real_front.back()) {
						if (std::find(multiseg_front.begin(), multiseg_front.end(), space_inx) == multiseg_front.end()) {
							multiseg_front.push_back(space_inx);
						}
					}
				}
				m_real_all_front_subspaces.push_back(space_inx);
			}
			temp_space = space_inx;
		}
		m_real_front_subspaces = real_front;
		m_multi_segment_subspaces = multiseg_front;

		m_operator_results;
		for (size_t i = 0; i < 8; ++i) {
			std::vector<size_t> temp = { 0,0,0 };
			m_operator_results.emplace_back(temp);
		}

		//SPMOEA::record();
	}

	int SPMOEA4_9::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		/********************************************************************************
							   ��������Ⱥǰ���ӿռ���������������Դ
		********************************************************************************/
		//�����������ͺ�����ͨ���򳤶ȷ��������Դ
		getInteractiveSols().clear();
		std::vector<size_t> assign_pop_resource;
		//size_t switch_period = 3;
		PopResourceAssign(assign_pop_resource, m_switch_period, pro);
		/********************************************************************************
											 ����Ⱥ�ݻ�
		********************************************************************************/
		//1.�Ӵ���̽���뿪���ı�����2.�Ӵ���̽���뿪���ķ�ʽ
		generateOffspring(pro, alg, rnd, assign_pop_resource);

		Population<Solution<>> offspring_pop;
		int tag = EvaluationTag::kNormalEval;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); j++) {
				offspring_pop.append(getPop()[i].getOffspring()[j]);
			}
		}
		SPMOEA::NDSort(offspring_pop);
		SPMOEA::updateHistoryInfo(offspring_pop, pro);
		updateSubspaceFrontSol(offspring_pop, pro, rnd);//�Ӵ������ӿռ�ǰ�ؽ�ʹ����
		updateFrontSpace();//ʹ������Ⱥ�Ӵ������ӿռ���Ϣ
		auto pre_front_space = getFrontSpace();
		/**********************************************************************************
					����Ƿ�ϸ���ӿռ�:���������Ⱥ���������ռ��ǰ���ӿռ�ϸ��
		***********************************************************************************/
		//�����Ӵ��ɹ�������Ӧ����
		auto pre_num_spaces = getMO_HLC().numSubspace();
		//adaptiveSplitSpace(pro, rnd);
		Real total_volume = getVarSpaceVolume();
		size_t partition_period = (size_t)std::floor(4 * std::sqrt(pro->numberVariables()));
		partition_period = std::ceil(m_maximum_evalutions / getPopsize() / (pro->numberVariables() * std::log(m_divide_granularity) / std::log(2)));
		//partition_period = 10000;
		//�ӿռ�ϸ�ַ�ʽ��1�����Ȼ��֣�2����λ��֣�3�������Գ̶Ȼ���
		// 1�����Ȼ���
		if ((m_divide_iteration + 1) % partition_period == 0) {
			auto divide_flag = spaceSubdivision(pro, rnd);
			if (divide_flag) {
				updateFrontSpace();
				//�����ӿռ�����
				for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
					getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
				}
			}
			updateVarSpaceRank(pro, rnd);
		}

		//�����ӿռ�������
		auto cur_num_spaces = getMO_HLC().numSubspace();
		//if (cur_num_spaces > pre_num_spaces) {
		//	updateFrontSpace();
		//	//�����ӿռ�����
		//	for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
		//		getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
		//	}
		//	//������Ⱥ�ֲ����Ӵ��ֲ�

		//}

		auto cur_front_spaces = getFrontSpace();
		/*for (size_t i = 0; i < cur_front_spaces.size(); ++i) {
			updateLinkSubspace(cur_front_spaces[i], cur_front_spaces, rnd);
		}*/
		setAddNumSpace(cur_num_spaces - pre_num_spaces);
		/**********************************************************************************
							   �ۺϸ�������Ⱥ��Ϣ�������ӿռ���Ϣ
		**********************************************************************************/
		/*if ((m_divide_iteration + 1) % 1 == 0) {
			updateVarSpaceRank(pro, rnd);
		}*/
		//updateVarSpaceDominanceRank(pro, rnd);
		updateFrontRegionLinkSpace(pro, rnd);
		setFrontLastGens(1);
		clusterSubspace();
		SPMOEA::updateNewPop(offspring_pop);
		SPMOEA::updateObjRange(offspring_pop, pro);
		//ʹ����ʷ���з�֧������archive
		//SPMOEA::updateArchive(archiveNum(), pro);
		//updateObjSpace();

		//testCoverage(pro);

		Real front_ratio = 0.;
		size_t front_count = 0;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < getFrontSpace().size(); ++i) {
			auto v = getMO_HLC().subspaceTree().getBoxVolume(getFrontSpace()[i]);
			front_ratio += (v / total_volume);
			front_count++;
		}
		setFrontSpaceRatio(front_ratio);
		/********************************************************************************************
										 ��Ⱥ���£�����Ⱥ�ڲ�������̭ѡ��
		*********************************************************************************************/
		updateInteractiveSols();
		auto his_front_sols = getHisFrontSols();
		if ((m_divide_iteration + 1) % partition_period == 0 && his_front_sols.size() > 2 * getPop()[0].size()) {
			for (size_t i = 0; i < getPop().size(); ++i) {
				Population<Solution<>> temp_pop;
				for (size_t j = 0; j < getPop()[i].size(); ++j) {
					temp_pop.append(getPop()[i][j]);
					temp_pop.append(getPop()[i].getOffspring()[j]);
				}
				auto his_front_sols = getHisFrontSols();
				for (size_t j = 0; j < his_front_sols.size(); ++j) {
					temp_pop.append(*his_front_sols[j]);
				}
				std::vector<std::vector<Real>> temp_data;
				std::vector<size_t> temp_ranks;
				SPMOEA::NDSort(temp_pop);
				for (size_t j = 0; j < temp_pop.size(); ++j) {
					temp_data.emplace_back(temp_pop[j].objective());
					temp_ranks.push_back(temp_pop[j].fitness());
				}
				auto sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, getPop()[i].size());
				for (size_t j = 0; j < sele_inx.size(); ++j) {
					getPop()[i][j] = temp_pop[sele_inx[j]];
				}
			}
		}
		else {
			subspaceClusterSelection(getPop()[0].size(), pro, rnd);
			//objspaceClusterSelection(getPop()[0].size(), pro, rnd, false);
			//clusterSelection(getPop()[0].size(), pro, rnd,2);
			//selectNewPop();
			//ensembleSelect(getPop()[0].size(), pro, rnd);
			//doubleClusterSelection(getPop()[0].size(), pro, rnd);
		}

		//if (m_local_selection) {
		//	//�ھֲ��������ӿռ�ѡ��
		//	for (size_t i = 0; i < getPop().size(); ++i) {
		//		std::vector<std::vector<std::pair<Real, Real>>> space_bounds=m_exploit_space_bounds;
		//		/*for (size_t j = 0; j < m_expliot_pos.size(); ++j) {
		//			auto space = getMO_HLC().subspaceTree().getRegionIdx(m_expliot_pos[j]);
		//			if (std::find(space_inx.begin(), space_inx.end(), space) == space_inx.end()) {
		//				space_inx.push_back(space);
		//			}
		//		}*/
		//		Population<Solution<>> temp_pop;
		//		size_t select_num = getPop()[i].size();
		//		std::vector<size_t> pop_inx;
		//		for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//			temp_pop.append(getPop()[i][j]);
		//			temp_pop.append(getPop()[i].getOffspring()[j]);
		//		}
		//		std::vector<std::vector<Real>> temp_data;
		//		std::vector<size_t> temp_ranks;
		//		SPMOEA::NDSort(temp_pop);
		//		for (size_t j = 0; j < temp_pop.size(); ++j) {
		//			temp_data.emplace_back(temp_pop[j].objective());
		//			temp_ranks.push_back(temp_pop[j].fitness());
		//		}

		//		std::vector<size_t> sele_inx;
		//		sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, select_num);
		//		for (size_t j = 0; j < sele_inx.size(); ++j) {
		//			getPop()[i][j] = temp_pop[sele_inx[j]];
		//		}

		//		////ѡ�������Щ�ӿռ��и��õ�һ�������Ľ�
		//		//Population<Solution<>> temp_pop;
		//		//size_t select_num = 0;
		//		//std::vector<size_t> pop_inx;
		//		//for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//		//	auto& sol = getPop()[i][j].variable().vect();
		//		//	for (size_t k = 0; k < space_bounds.size(); ++k) {
		//		//		if (indInRange(space_bounds[k],sol)) {
		//		//			select_num++;
		//		//			pop_inx.push_back(j);
		//		//			temp_pop.append(getPop()[i][j]);
		//		//			break;
		//		//		}
		//		//	}
		//		//}
		//		//for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//		//	auto& sol = getPop()[i].getOffspring()[j].variable().vect();
		//		//	for (size_t k = 0; k < space_bounds.size(); ++k) {
		//		//		if (indInRange(space_bounds[k], sol)) {
		//		//			temp_pop.append(getPop()[i].getOffspring()[j]);
		//		//			break;
		//		//		}
		//		//	}
		//		//}
		//		//size_t temp_size = temp_pop.size();
		//		//auto his_front_sols = getHisFrontSols();
		//		//for (size_t j = 0; j < his_front_sols.size(); ++j) {
		//		//	auto& sol = his_front_sols[j]->variable().vect();
		//		//	for (size_t k = 0; k < space_bounds.size(); ++k) {
		//		//		if (indInRange(space_bounds[k], sol)) {
		//		//			temp_pop.append(*his_front_sols[j]);
		//		//			break;
		//		//		}
		//		//	}
		//		//}
		//		//std::vector<std::vector<Real>> temp_data;
		//		//std::vector<size_t> temp_ranks;
		//		//SPMOEA::NDSort(temp_pop);
		//		//for (size_t j = 0; j < temp_pop.size(); ++j) {
		//		//	temp_data.emplace_back(temp_pop[j].objective());
		//		//	temp_ranks.push_back(temp_pop[j].fitness());
		//		//}

		//		//std::vector<size_t> sele_inx;
		//		//if (select_num > 0) {
		//		//	sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, select_num);
		//		//	for (size_t j = 0; j < sele_inx.size(); ++j) {
		//		//		getPop()[i][pop_inx[j]] = temp_pop[sele_inx[j]];
		//		//	}
		//		//}
		//		//else {
		//		//	std::vector<size_t> first_inx;
		//		//	size_t min_rank = *std::min_element(temp_ranks.begin(), temp_ranks.end());
		//		//	for (size_t j = 0; j < temp_ranks.size(); ++j) {
		//		//		if (temp_ranks[j] == min_rank) {
		//		//			first_inx.push_back(j);
		//		//		}
		//		//	}
		//		//	std::vector<size_t> selected_inx;
		//		//	for (size_t j = 0; j < first_inx.size(); ++j) {
		//		//		size_t se_inx = (size_t)std::floor(getPop()[i].size() * rnd->uniform.next());
		//		//		while (std::find(selected_inx.begin(), selected_inx.end(), se_inx) != selected_inx.end()) {
		//		//			se_inx = (size_t)std::floor(getPop()[i].size() * rnd->uniform.next());
		//		//		}
		//		//		getPop()[i][se_inx] = temp_pop[first_inx[j]];
		//		//		selected_inx.push_back(se_inx);
		//		//	}
		//		//}
		//	}
		//}
		//else {
		//	auto his_front_sols = getHisFrontSols();
		//	if ((m_divide_iteration + 1) % partition_period == 0 && his_front_sols.size() > 2 * getPop()[0].size()) {
		//		for (size_t i = 0; i < getPop().size(); ++i) {
		//			Population<Solution<>> temp_pop;
		//			for (size_t j = 0; j < getPop()[i].size(); ++j) {
		//				temp_pop.append(getPop()[i][j]);
		//				temp_pop.append(getPop()[i].getOffspring()[j]);
		//			}
		//			auto his_front_sols = getHisFrontSols();
		//			for (size_t j = 0; j < his_front_sols.size(); ++j) {
		//				temp_pop.append(*his_front_sols[j]);
		//			}
		//			std::vector<std::vector<Real>> temp_data;
		//			std::vector<size_t> temp_ranks;
		//			SPMOEA::NDSort(temp_pop);
		//			for (size_t j = 0; j < temp_pop.size(); ++j) {
		//				temp_data.emplace_back(temp_pop[j].objective());
		//				temp_ranks.push_back(temp_pop[j].fitness());
		//			}
		//			auto sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, getPop()[i].size());
		//			for (size_t j = 0; j < sele_inx.size(); ++j) {
		//				getPop()[i][j] = temp_pop[sele_inx[j]];
		//			}
		//		}
		//	}
		//	else {
		//		subspaceClusterSelection(getPop()[0].size(), pro, rnd);
		//		//objspaceClusterSelection(getPop()[0].size(), pro, rnd, false);
		//		//clusterSelection(getPop()[0].size(), pro, rnd,2);
		//		//selectNewPop();
		//		//ensembleSelect(getPop()[0].size(), pro, rnd);
		//		//doubleClusterSelection(getPop()[0].size(), pro, rnd);
		//	}
		//}

		/*********************************************************************************************
												 ��¼������Ϣ
		**********************************************************************************************/
		SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA4_9::selectNewPop() {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			SPMOEA::NDSort(getPop()[i].getOffspring());
			std::vector<std::vector<Real>> all_data;
			std::vector<size_t> ranks;
			size_t first_num = 0;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_data.emplace_back(getPop()[i].getOffspring()[j].objective());
				ranks.push_back(getPop()[i].getOffspring()[j].fitness());
				if (getPop()[i].getOffspring()[j].fitness() == 0) {
					first_num++;
				}
			}
			std::cout << "the number of front sols: " << first_num << std::endl;
			std::vector<size_t> sele_inx;
			sele_inx = selectMaxMinFromFront(all_data, ranks, getPop()[i].size());
			for (size_t j = 0; j < sele_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
			}
		}
	}

	//�ۺϲ�ͬ�������¸���ı���״̬������ۺ���������ѡ�����ذ������ӿռ�ֲ�ѡ������Ŀ��ռ�ѡ������
	//���߿ռ�ӵ����ѡ�����������Ӵ��ݻ�ѡ����
	void SPMOEA4_9::ensembleSelect(size_t select_num, Problem* pro, Random* rnd) {
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

			std::vector<size_t> select_inx;
			/*********************************************************************
											 �����ش��
			*********************************************************************/
			//�����ӿռ�������
			std::vector<Real> var_subspace_select_score(getPop()[i].getOffspring().size(), 0.);
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
				int layer_num = ofec::nd_sort::fastSort<Real>(objs, ranks, CAST_CONOP(pro)->optimizeMode());
				for (size_t j = 0; j < ranks.size(); ++j) {
					//var_subspace_select_score[jj.second[j]] = (1 - 0.3 * ranks[j]) > 0 ? (1 - 0.3 * ranks[j]) : 0.;
					//var_subspace_select_score[jj.second[j]] = (1. - ranks[j] / 5) > 0 ? (1. - ranks[j] / 5) : 0;
					var_subspace_select_score[jj.second[j]] = 1. - ranks[j] / (Real)layer_num;
				}
			}

			//Ŀ��ռ�������
			std::vector<Real> obj_space_select_score(getPop()[i].getOffspring().size(), 0.);
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
				//obj_space_select_score[j] = (1. - rank / 5.) > 0 ? (1. - rank / 5.) : 0;
				if (rank < layer_index) {
					obj_space_select_score[j] = 1.;
				}
				else {
					//obj_space_select_score[j] = (1 - 0.2 * rank) > 0 ? (1 - 0.2 * rank) : 0.;
					obj_space_select_score[j] = (1. - (rank - layer_index) / 5.) > 0 ? (1. - (rank - layer_index) / 5.) : 0;
				}
			}

			//���߿ռ�ӵ���ȴ�֣������ӿռ��ܶ�
			std::vector<Real> var_space_crowdist_select_score1(getPop()[i].getOffspring().size(), 0.);
			std::vector<Real> var_min_dist;//ƽ��ӵ������
			size_t kk = m_neighs;
			for (size_t j = 0; j < var_dist.size(); ++j) {
				Real temp = 0.;
				for (size_t k = 0; k < kk; ++k) {
					auto min_d = *std::min_element(var_dist[j].begin(), var_dist[j].end());
					temp += min_d;
					auto inx = std::distance(var_dist[j].begin(), std::min_element(var_dist[j].begin(), var_dist[j].end()));
					var_dist[j][inx] = INT16_MAX;
				}
				var_min_dist.push_back(temp / kk);
			}
			Real min_min_var_dist = *std::min_element(var_min_dist.begin(), var_min_dist.end());
			Real max_min_var_dist = *std::max_element(var_min_dist.begin(), var_min_dist.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				var_space_crowdist_select_score1[j] = (var_min_dist[j] - min_min_var_dist) / (max_min_var_dist - min_min_var_dist);
			}

			std::vector<Real> var_space_crowdist_select_score2(getPop()[i].getOffspring().size(), 0.);
			std::vector<Real> space_density;//�ӿռ��ܶ�
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto& sol = getPop()[i].getOffspring()[j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				auto volume = getMO_HLC().subspaceTree().getBoxVolume(space_inx);
				space_density.push_back(getMO_HLC().getSubspaceInfo(space_inx).m_history_inds.size() / volume);
			}
			std::vector<size_t> var_density_sort;
			for (size_t j = 0; j < space_density.size(); ++j) {
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
			auto max_density_sort = *std::max_element(var_density_sort.begin(), var_density_sort.end());
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				if (max_density_sort < 1) {
					var_space_crowdist_select_score2[j] = 0.;
				}
				else {
					var_space_crowdist_select_score2[j] = 1. - (Real)var_density_sort[j] / (Real)max_density_sort;
				}

			}

			//Ŀ��ռ�ӵ���ȴ��
			std::vector<Real> obj_space_crowdist_select_score(getPop()[i].getOffspring().size(), 0.);
			std::vector<Real> obj_min_dist;
			size_t oo = m_neighs;
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

			//�ۺϸ�����ѡ��
			std::vector<Real> all_scores(getPop()[i].getOffspring().size(), 0.);
			for (size_t j = 0; j < all_scores.size(); ++j) {
				all_scores[j] = 0.2 * var_subspace_select_score[j] + 0.2 * obj_space_select_score[j] + 0.2 * var_space_crowdist_select_score1[j] + 0.2 * var_space_crowdist_select_score2[j] + 0.2 * obj_space_crowdist_select_score[j];
			}
			//����ֵ�ɸߵ���ѡ��
			std::vector<size_t> select_flag(all_scores.size(), 0);
			for (size_t j = 0; j < all_scores.size(); ++j) {
				if (select_inx.size() >= select_num) {
					break;
				}
				else {
					Real max_v = 0.;
					size_t inx = 0;
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

	void SPMOEA4_9::selectInHisRange(size_t select_num, Problem* pro, Random* rnd) {
		auto front_obj_range = getFrontObjRange();
		for (size_t i = 0; i < front_obj_range.size(); ++i) {
			auto span = front_obj_range[i].second - front_obj_range[i].first;
			front_obj_range[i].second += span * 0.2;
		}
		auto front_spaces = getFrontSpace();
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			std::vector<size_t> in_bound_inx;
			std::vector<size_t> out_bound_inx;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto& obj = getPop()[i].getOffspring()[j].objective();
				if (indInRange(front_obj_range, obj)) {
					in_bound_inx.push_back(j);
				}
				else {
					out_bound_inx.push_back(j);
				}
			}
			std::vector<size_t> select_inx;
			//���ڲ������㹻ʱ�����ڲ�ѡ
			if (in_bound_inx.size() > getPop()[i].size()) {

			}
			else {
				//�ڲ�ȫѡ���ٴ��ⲿѡ��ʣ���
				for (size_t j = 0; j < in_bound_inx.size(); ++j) {
					select_inx.push_back(in_bound_inx[j]);
				}
				//�ⲿ����ѡ��
				size_t sele_num = getPop()[i].size() - select_inx.size();
				Population<Solution<>> temp_pop;
				for (size_t j = 0; j < out_bound_inx.size(); ++j) {
					temp_pop.append(getPop()[i].getOffspring()[out_bound_inx[j]]);
				}
				SPMOEA::NDSort(temp_pop);
				std::vector<std::vector<Real>> temp_data;
				std::vector<size_t> temp_ranks;
				for (size_t j = 0; j < temp_pop.size(); ++j) {
					temp_data.emplace_back(temp_pop[j].objective());
					temp_ranks.push_back(temp_pop[j].fitness());
				}
				std::vector<size_t> sele_inx;
				sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, sele_num);
				for (size_t j = 0; j < sele_inx.size(); ++j) {
					select_inx.push_back(out_bound_inx[sele_inx[j]]);
				}
			}
			for (size_t j = 0; j < select_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_inx[j]];
			}
		}
	}

	void SPMOEA4_9::subspaceClusterSelection(size_t select_num, Problem* pro, Random* rnd) {
		auto front_spaces = getFrontSpace();
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			std::map<size_t, std::vector<size_t>> pop2space;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto var = getPop()[i].getOffspring()[j].variable().vect();
				auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
				if (pop2space[space_idx].empty()) {
					std::vector<size_t> temp_inx;
					pop2space.insert(std::make_pair(space_idx, temp_inx));
				}
				pop2space[space_idx].emplace_back(j);
			}

			////������Щ�ӿռ�ǰ�ؽ���������
			//std::vector<size_t> spaces;
			//for (auto jj : pop2space) {
			//	spaces.push_back(jj.first);
			//}
			//std::vector<std::vector<Real>> space_min_dist;
			//for (size_t j = 0; j < spaces.size(); ++j) {
			//	std::vector<std::vector<Real>> vars;
			//	if (std::find(front_spaces.begin(), front_spaces.end(), spaces[j]) == front_spaces.end()) {
			//		auto& front_sols1 = getMO_HLC().getSubspaceInfo(spaces[j]).m_subspace_front_sol;
			//		for (size_t k = 0; k < front_sols1.size(); ++k) {
			//			vars.emplace_back(front_sols1[k]->variable().vect());
			//		}
			//	}
			//	else {
			//		auto& front_sols1 = getMO_HLC().getSubspaceInfo(spaces[j]).m_front_sol_in_subspace;
			//		for (size_t k = 0; k < front_sols1.size(); ++k) {
			//			vars.emplace_back(front_sols1[k]->variable().vect());
			//		}
			//	}
			//	/*auto& front_sols1 = getMO_HLC().getSubspaceInfo(spaces[j]).m_subspace_front_sol;
			//	for (size_t k = 0; k < front_sols1.size(); ++k) {
			//		vars.emplace_back(front_sols1[k]->variable().vect());
			//	}*/
			//	std::vector<Real> temp_min_dist;//�ӿռ�֮�����̾���
			//	//temp_min_dist.resize(spaces.size());
			//	Real min_dist = std::numeric_limits<Real>::max();//�ӿռ�֮�����̾����е���С����
			//	for (size_t k = 0; k < spaces.size(); ++k) {
			//		if (k <= j) {
			//			temp_min_dist.push_back(0.);
			//		}
			//		else {
			//			std::vector<std::vector<Real>> vars2;
			//			if (std::find(front_spaces.begin(), front_spaces.end(), spaces[k]) == front_spaces.end()) {
			//				auto& front_sols2 = getMO_HLC().getSubspaceInfo(spaces[k]).m_subspace_front_sol;
			//				for (size_t p = 0; p < front_sols2.size(); ++p) {
			//					vars2.emplace_back(front_sols2[p]->variable().vect());
			//				}
			//			}
			//			else {
			//				auto& front_sols2 = getMO_HLC().getSubspaceInfo(spaces[k]).m_front_sol_in_subspace;
			//				for (size_t p = 0; p < front_sols2.size(); ++p) {
			//					vars2.emplace_back(front_sols2[p]->variable().vect());
			//				}
			//			}
			//			/*auto& front_sols2 = getMO_HLC().getSubspaceInfo(spaces[k]).m_subspace_front_sol;
			//			for (size_t p = 0; p < front_sols2.size(); ++p) {
			//				vars2.emplace_back(front_sols2[p]->variable().vect());
			//			}*/
			//			//���������ӿռ�ǰ���е���̾���ֵ
			//			//��һ��
			//			dataNormalizeInBound(vars,bound);
			//			dataNormalizeInBound(vars2, bound);
			//			Real min_v = std::numeric_limits<Real>::max();
			//			size_t idx1 = 0, idx2 = 0;
			//			for (size_t p = 0; p < vars.size(); ++p) {
			//				for (size_t q = 0; q < vars2.size(); ++q) {
			//					auto dist = euclideanDistance(vars[p].begin(), vars[p].end(), vars2[q].begin());
			//					min_v = min_v < dist ? min_v : dist;
			//					idx1 = p;
			//					idx2 = q;
			//				}
			//			}
			//			if (min_v < min_dist) {
			//				min_dist = min_v;
			//			}
			//			temp_min_dist.push_back(min_v);
			//		}
			//		
			//	}
			//	space_min_dist.emplace_back(temp_min_dist);
			//}
			//for (size_t j = 0; j < space_min_dist.size(); ++j) {
			//	for (size_t k = 0; k <= j; ++k) {
			//		if (k == j) {
			//			space_min_dist[j][k] = std::numeric_limits<Real>::max();
			//		}
			//		else {
			//			space_min_dist[j][k] = space_min_dist[k][j];
			//		}
			//	}
			//}

			//������Щ�ӿռ��ڽ�֮����������
			std::vector<size_t> spaces;
			for (auto jj : pop2space) {
				spaces.push_back(jj.first);
			}
			std::vector<std::vector<Real>> space_min_dist;
			for (size_t j = 0; j < spaces.size(); ++j) {
				std::vector<std::vector<Real>> vars;
				auto ind_inx = pop2space[spaces[j]];
				for (size_t k = 0; k < ind_inx.size(); ++k) {
					vars.emplace_back(getPop()[i].getOffspring()[ind_inx[k]].variable().vect());
				}
				std::vector<Real> temp_min_dist;//�ӿռ��ڽ�֮�����̾���
				//temp_min_dist.resize(spaces.size());
				Real min_dist = std::numeric_limits<Real>::max();//�ӿռ�֮�����̾����е���С����
				for (size_t k = 0; k < spaces.size(); ++k) {
					if (k <= j) {
						temp_min_dist.push_back(0.);
					}
					else {
						std::vector<std::vector<Real>> vars2;
						auto ind_inx2 = pop2space[spaces[k]];
						for (size_t p = 0; p < ind_inx2.size(); ++p) {
							vars2.emplace_back(getPop()[i].getOffspring()[ind_inx2[p]].variable().vect());
						}
						//���������ӿռ��ڽ����̾���ֵ
						//��һ��
						dataNormalizeInBound(vars, bound);
						dataNormalizeInBound(vars2, bound);
						Real min_v = std::numeric_limits<Real>::max();
						size_t idx1 = 0, idx2 = 0;
						for (size_t p = 0; p < vars.size(); ++p) {
							for (size_t q = 0; q < vars2.size(); ++q) {
								auto dist = euclideanDistance(vars[p].begin(), vars[p].end(), vars2[q].begin());
								min_v = min_v < dist ? min_v : dist;
								idx1 = p;
								idx2 = q;
							}
						}
						if (min_v < min_dist) {
							min_dist = min_v;
						}
						temp_min_dist.push_back(min_v);
					}
				}
				space_min_dist.emplace_back(temp_min_dist);
			}
			for (size_t j = 0; j < space_min_dist.size(); ++j) {
				for (size_t k = 0; k <= j; ++k) {
					if (k == j) {
						space_min_dist[j][k] = std::numeric_limits<Real>::max();
					}
					else {
						space_min_dist[j][k] = space_min_dist[k][j];
					}
				}
			}

			//�����ӿռ�֮�����С�����γɵıջ�����
			SpaceDirectedGraph space_graph;
			for (size_t j = 0; j < spaces.size(); ++j) {
				//��ӽڵ�ͱ�
				size_t next_space = std::distance(space_min_dist[j].begin(), std::min_element(space_min_dist[j].begin(), space_min_dist[j].end()));
				space_graph.addNode(spaces[j]);
				//��ӱ�
				Real dist;
				if (j == next_space) {
					dist = 0.;
				}
				else {
					dist = space_min_dist[j][next_space];
				}
				space_graph.addEdge(spaces[j], spaces[next_space], dist);
			}

			std::vector<std::vector<size_t>> path;
			for (auto ss : space_graph.getGraphNode()) {
				//�õ��Դ��ӿռ�Ľڵ�·��
				std::vector<size_t> temp_path;
				temp_path.push_back(ss.first);
				size_t next_node_inx = ss.second->getBehindNeighbors().front()->getSpaceInx();
				while (std::find(temp_path.begin(), temp_path.end(), next_node_inx) == temp_path.end()) {
					temp_path.push_back(next_node_inx);
					next_node_inx = space_graph.getGraphNode()[next_node_inx]->getBehindNeighbors().front()->getSpaceInx();
				}
				path.emplace_back(temp_path);
			}

			std::vector<std::vector<size_t>> common_flag(path.size());
			for (size_t j = 0; j < path.size(); ++j) {
				common_flag[j].resize(path.size());
				std::set<size_t> set1;
				for (size_t k = 0; k < path[j].size(); ++k) {
					set1.insert(path[j][k]);
				}
				for (size_t k = j + 1; k < path.size(); ++k) {
					std::set<size_t> set2 = set1;
					for (size_t p = 0; p < path[k].size(); ++p) {
						set2.insert(path[k][p]);
					}
					if (set2.size() < set1.size() + path[k].size()) {
						common_flag[j][k] = 1;
					}
				}
			}
			std::vector<std::vector<size_t>> space_clusters;
			std::vector<size_t> merge_flag(path.size(), 0);
			while (std::find(merge_flag.begin(), merge_flag.end(), 0) != merge_flag.end()) {
				size_t begin_index;
				for (size_t j = 0; j < merge_flag.size(); ++j) {
					if (merge_flag[j] == 0) {
						begin_index = j;
						break;
					}
				}
				std::vector<size_t> temp_cluster = path[begin_index];
				merge_flag[begin_index] = 1;
				std::vector<size_t> merge_inx;
				for (size_t j = 0; j < common_flag[begin_index].size(); ++j) {
					if (common_flag[begin_index][j] == 1) {
						merge_inx.push_back(j);
					}
				}
				while (merge_inx.size() > 0) {
					for (size_t j = 0; j < merge_inx.size(); ++j) {
						if (merge_flag[merge_inx[j]] == 0) {
							for (size_t k = 0; k < path[merge_inx[j]].size(); ++k) {
								if (std::find(temp_cluster.begin(), temp_cluster.end(), path[merge_inx[j]][k]) == temp_cluster.end()) {
									temp_cluster.push_back(path[merge_inx[j]][k]);
								}
							}
							merge_flag[merge_inx[j]] = 1;
						}
					}
					std::vector<size_t> temp_merge_inx;
					for (size_t j = 0; j < merge_inx.size(); ++j) {
						for (size_t k = 0; k < common_flag[merge_inx[j]].size(); ++k) {
							if (merge_flag[k] == 0) {
								if (common_flag[merge_inx[j]][k] == 1) {
									temp_merge_inx.push_back(k);
								}
							}

						}
					}
					merge_inx = temp_merge_inx;
				}
				space_clusters.emplace_back(temp_cluster);
			}

			//�Ӿ�����ӿռ���ѡ�����
			//ÿһ��ѡȡǰ�أ�Ȼ���������ѡ����ʷ��֧��⸨��ѡ��
			auto his_front_sols = getHisFrontSols();
			auto front_obj_range = getFrontObjRange();
			for (size_t i = 0; i < front_obj_range.size(); ++i) {
				auto span = front_obj_range[i].second - front_obj_range[i].first;
				front_obj_range[i].second += span * 0.2;
			}
			std::vector<size_t> cluster_ind_num;
			std::vector<size_t> survive_num;
			for (size_t j = 0; j < space_clusters.size(); ++j) {
				size_t num = 0;
				for (size_t k = 0; k < space_clusters[j].size(); ++k) {
					auto ind_inx = pop2space[space_clusters[j][k]];
					num += ind_inx.size();
				}
				cluster_ind_num.push_back(num);
			}
			std::cout << "the number clusters in var space: " << space_clusters.size() << std::endl;
			std::vector<size_t> select_index;
			std::vector<std::vector<size_t>> pop_inx_clusters;//�����ĸ�������
			std::vector<std::vector<size_t>> sele_in_clusters;
			std::vector<std::map<size_t, std::vector<size_t>>> cluster_layer_index;//����ֲ�����
			for (size_t j = 0; j < space_clusters.size(); ++j) {
				Population<Solution<>> temp_pop;
				std::vector<size_t> pop_inx;
				for (size_t k = 0; k < space_clusters[j].size(); ++k) {
					auto ind_inx = pop2space[space_clusters[j][k]];
					for (size_t p = 0; p < ind_inx.size(); ++p) {
						temp_pop.append(getPop()[i].getOffspring()[ind_inx[p]]);
						pop_inx.push_back(ind_inx[p]);
					}
				}
				pop_inx_clusters.emplace_back(pop_inx);
				std::map<size_t, std::vector<size_t>> layer_index;//�ֲ�λ������
				//��ʷǰ�ظ�������
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					temp_pop.append(*his_front_sols[k]);
				}
				auto layer_num = layerNDSort(temp_pop);
				std::vector<size_t> pop_rank;
				for (size_t k = 0; k < pop_inx.size(); ++k) {
					pop_rank.push_back(temp_pop[k].fitness());
				}
				for (size_t k = 0; k < pop_inx.size(); ++k) {
					auto& obj = temp_pop[k].objective();
					if (!indInRange(front_obj_range, obj)) {
						pop_rank[k] += layer_num;
					}
				}
				size_t max_rank = *std::max_element(pop_rank.begin(), pop_rank.end());
				for (size_t k = 0; k < 1 + max_rank; ++k) {
					std::vector<size_t> t_inx;
					for (size_t p = 0; p < pop_rank.size(); ++p) {
						if (pop_rank[p] == k) {
							t_inx.push_back(p);
						}
					}
					if (!t_inx.empty()) {
						layer_index.insert(std::make_pair<>(k, t_inx));
					}
				}
				cluster_layer_index.emplace_back(layer_index);
			}
			m_pop_index_clusters.clear();
			m_pop_index_clusters = pop_inx_clusters;
			//������ѡ����������ӵ������̭�����
			size_t select_rank = 0;
			while (select_index.size() < getPop()[i].size()) {
				//�����������ѡ��
				//�ȵõ���ǰ�������
				std::vector<size_t> candidate;
				for (size_t j = 0; j < cluster_layer_index.size(); ++j) {
					for (size_t k = 0; k < cluster_layer_index[j][select_rank].size(); ++k) {
						candidate.push_back(pop_inx_clusters[j][cluster_layer_index[j][select_rank][k]]);
					}
				}
				if (select_index.size() + candidate.size() <= getPop()[i].size()) {
					//��ǰ��ȫѡ
					for (size_t j = 0; j < candidate.size(); ++j) {
						select_index.push_back(candidate[j]);
					}
					select_rank++;
				}
				else {
					//�ӵ�ǰ��ѡ��һ������
					size_t num = getPop()[i].size() - select_index.size();
					std::vector<std::vector<Real>> temp_data;
					std::vector<size_t> temp_ranks(candidate.size(), 0);
					for (size_t j = 0; j < candidate.size(); ++j) {
						temp_data.emplace_back(getPop()[i].getOffspring()[candidate[j]].objective());
					}
					std::vector<size_t> sele_inx;
					sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, num);
					for (size_t j = 0; j < sele_inx.size(); ++j) {
						select_index.push_back(candidate[sele_inx[j]]);
					}
				}
			}
			for (size_t j = 0; j < select_index.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_index[j]];
			}
		}
	}

	void SPMOEA4_9::objspaceClusterSelection(size_t select_num, Problem* pro, Random* rnd, bool add_his) {
		auto front_spaces = getFrontSpace();
		size_t num_obj = pro->numberObjectives();
		auto his_front_sols = getHisFrontSols();
		std::vector<std::vector<Real>> all_data;
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				all_data.emplace_back(getPop()[i].getOffspring()[j].objective());
			}
			if (add_his) {//������ʷ��֧���
				for (size_t j = 0; j < his_front_sols.size(); ++j) {
					all_data.emplace_back(his_front_sols[j]->objective());
				}
			}

			auto normal_data = all_data;
			dataNormalize(normal_data);
			std::vector<Real> min_dist;
			std::vector<std::vector<Real>> all_dist;
			Real max_dist = 0;
			for (size_t i = 0; i < normal_data.size(); ++i) {
				std::vector<Real> temp_dist;
				for (size_t j = 0; j < normal_data.size(); ++j) {
					if (j != i) {
						auto dist = euclideanDistance(normal_data[i].begin(), normal_data[i].end(), normal_data[j].begin());
						temp_dist.push_back(dist);
					}
					else {
						temp_dist.push_back(0.);
					}
				}
				all_dist.emplace_back(temp_dist);
				min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
				auto temp = *std::max_element(temp_dist.begin(), temp_dist.end());
				if (temp > max_dist) {
					max_dist = temp;
				}
			}
			Real sum_dist = 0.;
			for (size_t j = 0; j < min_dist.size(); ++j) {
				sum_dist += min_dist[j];
			}
			//dbscan����
			//�������ξ�����max_dist
			Real mean_dist = max_dist / normal_data.size();
			//��ƽ��������Ϊ�ܶȾ���ľ���ֵ
			size_t minPts = 3;//ÿ��������С�ĸ�����
			Real epsilon = mean_dist;

			DBSCAN dscluster(minPts, epsilon, normal_data);
			dscluster.run();
			std::vector<int> cluster_id;//ÿ������������,��ʵ���δ�����
			for (size_t j = 0; j < dscluster.m_points.size(); ++j) {
				cluster_id.push_back(dscluster.m_points[j]->clusterID);
			}
			std::vector<int> cluster_num;//���е���
			cluster_num.push_back(cluster_id[0]);
			for (size_t j = 0; j < cluster_id.size(); ++j) {
				if (std::find(cluster_num.begin(), cluster_num.end(), cluster_id[j]) == cluster_num.end()) {
					cluster_num.push_back(cluster_id[j]);
				}
			}
			std::vector<std::vector<size_t>> ind2cluster;//�����������
			std::vector<size_t> noise_inx;
			//����ͬ����ĵ����һ��
			for (size_t j = 0; j < cluster_num.size(); ++j) {
				if (cluster_num[j] > 0) {
					std::vector<size_t> temp_cluster;
					for (size_t k = 0; k < cluster_id.size(); ++k) {
						if (cluster_id[k] == cluster_num[j]) {
							temp_cluster.push_back(k);
						}
					}
					ind2cluster.emplace_back(temp_cluster);
				}
			}
			for (size_t k = 0; k < cluster_id.size(); ++k) {
				if (cluster_id[k] < 0) {
					noise_inx.push_back(k);
				}
			}
			ind2cluster.emplace_back(noise_inx);
			//�Ӿ�����ӿռ���ѡ�����
			//ÿһ��ѡȡǰ�أ�Ȼ���������ѡ����ʷ��֧��⸨��ѡ��
			auto front_obj_range = getFrontObjRange();
			for (size_t j = 0; j < front_obj_range.size(); ++j) {
				auto span = front_obj_range[j].second - front_obj_range[j].first;
				front_obj_range[j].second += span * 0.2;
			}
			std::vector<size_t> cluster_ind_num;
			std::vector<size_t> survive_num;
			for (size_t j = 0; j < ind2cluster.size(); ++j) {
				cluster_ind_num.push_back(ind2cluster[j].size());
			}
			std::cout << "the number clusters in var space: " << ind2cluster.size() << std::endl;
			std::vector<size_t> select_index;
			std::vector<std::vector<size_t>> pop_inx_clusters;//�����ĸ�������
			std::vector<std::map<size_t, std::vector<size_t>>> cluster_layer_index;//����ֲ�����
			for (size_t j = 0; j < ind2cluster.size(); ++j) {
				Population<Solution<>> temp_pop;
				for (size_t k = 0; k < ind2cluster[j].size(); ++k) {
					temp_pop.append(getPop()[i].getOffspring()[ind2cluster[j][k]]);
				}
				pop_inx_clusters.emplace_back(ind2cluster[j]);
				std::map<size_t, std::vector<size_t>> layer_index;//�ֲ�λ������
				//��ʷǰ�ظ�������
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					temp_pop.append(*his_front_sols[k]);
				}
				auto layer_num = layerNDSort(temp_pop);
				std::vector<size_t> pop_rank;
				for (size_t k = 0; k < ind2cluster[j].size(); ++k) {
					pop_rank.push_back(temp_pop[k].fitness());
				}
				for (size_t k = 0; k < ind2cluster[j].size(); ++k) {
					auto& obj = temp_pop[k].objective();
					if (!indInRange(front_obj_range, obj)) {
						pop_rank[k] += layer_num;
					}
				}
				size_t max_rank = *std::max_element(pop_rank.begin(), pop_rank.end());
				for (size_t k = 0; k < 1 + max_rank; ++k) {
					std::vector<size_t> t_inx;
					for (size_t p = 0; p < pop_rank.size(); ++p) {
						if (pop_rank[p] == k) {
							t_inx.push_back(p);
						}
					}
					if (!t_inx.empty()) {
						layer_index.insert(std::make_pair<>(k, t_inx));
					}
				}
				cluster_layer_index.emplace_back(layer_index);
			}
			m_pop_index_clusters.clear();
			m_pop_index_clusters = pop_inx_clusters;
			//������ѡ����������ӵ������̭�����
			size_t select_rank = 0;
			while (select_index.size() < getPop()[i].size()) {
				//�����������ѡ��
				//�ȵõ���ǰ�������
				std::vector<size_t> candidate;
				for (size_t j = 0; j < cluster_layer_index.size(); ++j) {
					for (size_t k = 0; k < cluster_layer_index[j][select_rank].size(); ++k) {
						candidate.push_back(pop_inx_clusters[j][cluster_layer_index[j][select_rank][k]]);
					}
				}
				if (select_index.size() + candidate.size() <= getPop()[i].size()) {
					//��ǰ��ȫѡ
					for (size_t j = 0; j < candidate.size(); ++j) {
						select_index.push_back(candidate[j]);
					}
					select_rank++;
				}
				else {
					//�ӵ�ǰ��ѡ��һ������
					size_t num = getPop()[i].size() - select_index.size();
					std::vector<std::vector<Real>> temp_data;
					std::vector<size_t> temp_ranks(candidate.size(), 0);
					for (size_t j = 0; j < candidate.size(); ++j) {
						temp_data.emplace_back(getPop()[i].getOffspring()[candidate[j]].objective());
					}
					std::vector<size_t> sele_inx;
					sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, num);
					for (size_t j = 0; j < sele_inx.size(); ++j) {
						select_index.push_back(candidate[sele_inx[j]]);
					}
				}
			}
			for (size_t j = 0; j < select_index.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_index[j]];
			}
		}
	}

	void SPMOEA4_9::clusterSelection(size_t select_num, Problem* pro, Random* rnd, size_t m_num_neigh) {
		auto front_spaces = getFrontSpace();
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			std::map<size_t, std::vector<size_t>> pop2space;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto var = getPop()[i].getOffspring()[j].variable().vect();
				auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
				if (pop2space[space_idx].empty()) {
					std::vector<size_t> temp_inx;
					pop2space.insert(std::make_pair(space_idx, temp_inx));
				}
				pop2space[space_idx].emplace_back(j);
			}

			std::vector<std::vector<Real>> vars;
			for (size_t k = 0; k < getPop()[i].getOffspring().size(); ++k) {
				vars.emplace_back(getPop()[i].getOffspring()[k].variable().vect());
			}
			//��һ��
			dataNormalizeInBound(vars, bound);
			//���������ӿռ��ڽ����̾���ֵ
			std::vector<std::vector<Real>> all_dist;
			for (size_t j = 0; j < vars.size(); ++j) {
				std::vector<Real> ind_dist;//�ӿռ��ڽ�֮�����̾���
				Real min_dist = std::numeric_limits<Real>::max();//�ӿռ�֮�����̾����е���С����
				for (size_t k = 0; k < vars.size(); ++k) {
					Real dist = 0.;
					if (k > j) {
						dist = euclideanDistance(vars[j].begin(), vars[j].end(), vars[k].begin());
					}
					ind_dist.push_back(dist);
				}
				all_dist.emplace_back(ind_dist);
			}
			for (size_t j = 0; j < all_dist.size(); ++j) {
				for (size_t k = 0; k <= j; ++k) {
					if (k == j) {
						all_dist[j][k] = std::numeric_limits<Real>::max();
					}
					else {
						all_dist[j][k] = all_dist[k][j];
					}
				}
			}
			std::vector<std::vector<size_t>> ind_clusters;

			//���ݸ���֮�����С�����γɵıջ�����
			SpaceDirectedGraph ind_graph;
			for (size_t j = 0; j < all_dist.size(); ++j) {
				//��ӽڵ�ͱ�
				ind_graph.addNode(j);
				for (size_t k = 0; k < m_num_neigh; ++k) {
					size_t next_ind = std::distance(all_dist[j].begin(), std::min_element(all_dist[j].begin(), all_dist[j].end()));
					//��ӱ�
					Real dist;
					if (j == next_ind) {
						dist = 0.;
					}
					else {
						dist = all_dist[j][next_ind];
					}
					ind_graph.addEdge(j, next_ind, dist);
					all_dist[j][next_ind] = std::numeric_limits<Real>::max();
				}
			}

			std::vector<std::vector<size_t>> path;
			for (auto ss : ind_graph.getGraphNode()) {
				//�õ��Դ˽ڵ�����ıջ�·��
				std::vector<size_t> temp_path;
				temp_path.push_back(ss.first);
				std::vector<size_t> behind_node_inx;
				auto behind_node = ind_graph.getGraphNode()[ss.first]->getBehindNeighbors();
				for (size_t k = 0; k < behind_node.size(); ++k) {
					behind_node_inx.push_back(behind_node[k]->getSpaceInx());
				}
				std::vector<size_t> new_node_inx;
				for (size_t k = 0; k < behind_node_inx.size(); ++k) {
					if (std::find(temp_path.begin(), temp_path.end(), behind_node_inx[k]) == temp_path.end()) {
						temp_path.push_back(behind_node_inx[k]);
						new_node_inx.push_back(behind_node_inx[k]);
					}
				}
				while (new_node_inx.size() > 0) {
					std::vector<size_t> new_behind_inx;
					for (size_t k = 0; k < new_node_inx.size(); ++k) {
						auto behind_node = ind_graph.getGraphNode()[new_node_inx[k]]->getBehindNeighbors();
						for (size_t p = 0; p < behind_node.size(); ++p) {
							size_t index = behind_node[p]->getSpaceInx();
							if (std::find(temp_path.begin(), temp_path.end(), index) == temp_path.end()) {
								temp_path.push_back(index);
								new_behind_inx.push_back(index);
							}

						}
					}
					new_node_inx = new_behind_inx;
				}
				path.emplace_back(temp_path);
			}

			std::vector<std::vector<size_t>> common_flag(path.size());
			for (size_t j = 0; j < path.size(); ++j) {
				common_flag[j].resize(path.size());
				std::set<size_t> set1;
				for (size_t k = 0; k < path[j].size(); ++k) {
					set1.insert(path[j][k]);
				}
				for (size_t k = j + 1; k < path.size(); ++k) {
					std::set<size_t> set2 = set1;
					for (size_t p = 0; p < path[k].size(); ++p) {
						set2.insert(path[k][p]);
					}
					if (set2.size() < set1.size() + path[k].size()) {
						common_flag[j][k] = 1;
					}
				}
			}
			std::vector<size_t> merge_flag(path.size(), 0);
			while (std::find(merge_flag.begin(), merge_flag.end(), 0) != merge_flag.end()) {
				size_t begin_index;
				for (size_t j = 0; j < merge_flag.size(); ++j) {
					if (merge_flag[j] == 0) {
						begin_index = j;
						break;
					}
				}
				std::vector<size_t> temp_cluster = path[begin_index];
				merge_flag[begin_index] = 1;
				std::vector<size_t> merge_inx;
				for (size_t j = 0; j < common_flag[begin_index].size(); ++j) {
					if (common_flag[begin_index][j] == 1) {
						merge_inx.push_back(j);
					}
				}
				while (merge_inx.size() > 0) {
					for (size_t j = 0; j < merge_inx.size(); ++j) {
						if (merge_flag[merge_inx[j]] == 0) {
							for (size_t k = 0; k < path[merge_inx[j]].size(); ++k) {
								if (std::find(temp_cluster.begin(), temp_cluster.end(), path[merge_inx[j]][k]) == temp_cluster.end()) {
									temp_cluster.push_back(path[merge_inx[j]][k]);
								}
							}
							merge_flag[merge_inx[j]] = 1;
						}
					}
					std::vector<size_t> temp_merge_inx;
					for (size_t j = 0; j < merge_inx.size(); ++j) {
						for (size_t k = 0; k < common_flag[merge_inx[j]].size(); ++k) {
							if (merge_flag[k] == 0) {
								if (common_flag[merge_inx[j]][k] == 1) {
									temp_merge_inx.push_back(k);
								}
							}

						}
					}
					merge_inx = temp_merge_inx;
				}
				ind_clusters.emplace_back(temp_cluster);
			}

			//�Ӿ�����ӿռ���ѡ�����
			//ÿһ��ѡȡǰ�أ�Ȼ���������ѡ����ʷ��֧��⸨��ѡ��
			auto his_front_sols = getHisFrontSols();
			auto front_obj_range = getFrontObjRange();
			for (size_t i = 0; i < front_obj_range.size(); ++i) {
				auto span = front_obj_range[i].second - front_obj_range[i].first;
				front_obj_range[i].second += span * 0.2;
			}
			std::vector<size_t> cluster_ind_num;
			std::vector<size_t> survive_num;
			for (size_t j = 0; j < ind_clusters.size(); ++j) {
				cluster_ind_num.push_back(ind_clusters[j].size());
			}
			std::cout << "the number clusters in var space: " << ind_clusters.size() << std::endl;
			std::vector<size_t> select_index;
			std::vector<std::vector<size_t>> pop_inx_clusters;//�����ĸ�������
			std::vector<std::vector<size_t>> sele_in_clusters;
			std::vector<std::map<size_t, std::vector<size_t>>> cluster_layer_indexs;//����ֲ�����
			for (size_t j = 0; j < ind_clusters.size(); ++j) {
				Population<Solution<>> temp_pop;
				for (size_t k = 0; k < ind_clusters[j].size(); ++k) {
					temp_pop.append(getPop()[i].getOffspring()[ind_clusters[j][k]]);
				}
				pop_inx_clusters.emplace_back(ind_clusters[j]);
				std::map<size_t, std::vector<size_t>> layer_index;//�ֲ�λ������
				//��ʷǰ�ظ�������
				for (size_t k = 0; k < his_front_sols.size(); ++k) {
					temp_pop.append(*his_front_sols[k]);
				}
				auto layer_num = layerNDSort(temp_pop);
				std::vector<size_t> pop_rank;
				for (size_t k = 0; k < ind_clusters[j].size(); ++k) {
					pop_rank.push_back(temp_pop[k].fitness());
				}
				for (size_t k = 0; k < ind_clusters[j].size(); ++k) {
					auto& obj = temp_pop[k].objective();
					if (!indInRange(front_obj_range, obj)) {
						pop_rank[k] += layer_num;
					}
				}
				size_t max_rank = *std::max_element(pop_rank.begin(), pop_rank.end());
				for (size_t k = 0; k < 1 + max_rank; ++k) {
					std::vector<size_t> t_inx;
					for (size_t p = 0; p < pop_rank.size(); ++p) {
						if (pop_rank[p] == k) {
							t_inx.push_back(p);
						}
					}
					if (!t_inx.empty()) {
						layer_index.insert(std::make_pair<>(k, t_inx));
					}
				}
				cluster_layer_indexs.emplace_back(layer_index);
			}
			m_pop_index_clusters.clear();
			m_pop_index_clusters = pop_inx_clusters;
			//������ѡ����������ӵ������̭�����
			size_t select_rank = 0;
			while (select_index.size() < getPop()[i].size()) {
				//�����������ѡ��
				//�ȵõ���ǰ�������
				std::vector<size_t> candidate;
				for (size_t j = 0; j < cluster_layer_indexs.size(); ++j) {
					for (size_t k = 0; k < cluster_layer_indexs[j][select_rank].size(); ++k) {
						candidate.push_back(pop_inx_clusters[j][cluster_layer_indexs[j][select_rank][k]]);
					}
				}
				//���ظ�����
				std::set<size_t> inx;
				for (size_t j = 0; j < candidate.size(); ++j) {
					if (std::find(select_index.begin(), select_index.end(), candidate[j]) == select_index.end()) {
						inx.insert(candidate[j]);
					}
				}
				if (select_index.size() + inx.size() <= getPop()[i].size()) {
					//��ǰ��ȫѡ
					for (auto j : inx) {
						select_index.push_back(j);
					}
					select_rank++;
				}
				else {
					//�ӵ�ǰ��ѡ��һ������
					size_t num = getPop()[i].size() - select_index.size();
					std::vector<std::vector<Real>> temp_data;
					std::vector<size_t> temp_inx;
					for (auto j : inx) {
						temp_inx.push_back(j);
					}
					std::vector<size_t> temp_ranks(inx.size(), 0);
					for (auto j : inx) {
						temp_data.emplace_back(getPop()[i].getOffspring()[j].objective());
					}
					std::vector<size_t> sele_inx;
					sele_inx = selectMaxMinFromFront(temp_data, temp_ranks, num);
					for (size_t j = 0; j < sele_inx.size(); ++j) {
						select_index.push_back(temp_inx[sele_inx[j]]);
					}
				}
			}
			for (size_t j = 0; j < select_index.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[select_index[j]];
			}
		}
	}

	void SPMOEA4_9::clusterInVarspace(Problem* pro) {
		m_clusters_in_var.clear();
		auto front_spaces = getFrontSpace();
		auto bound = CAST_CONOP(pro)->boundary();
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				getPop()[i].getOffspring()[getPop()[i].size() + j] = getPop()[i][j];
			}
			//�����ӿռ����
			std::map<size_t, std::vector<size_t>> pop2space;
			for (size_t j = 0; j < getPop()[i].getOffspring().size(); ++j) {
				auto var = getPop()[i].getOffspring()[j].variable().vect();
				auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
				if (pop2space[space_idx].empty()) {
					std::vector<size_t> temp_inx;
					pop2space.insert(std::make_pair(space_idx, temp_inx));
				}
				pop2space[space_idx].emplace_back(j);
			}
			//ǰ���ӿռ����
			//������Щ�ӿռ�ǰ�ؽ���������
			std::vector<size_t> spaces;//ǰ���ӿռ�
			std::vector<size_t> behind_spaces;
			for (auto jj : pop2space) {
				if (std::find(front_spaces.begin(), front_spaces.end(), jj.first) != front_spaces.end()) {
					spaces.push_back(jj.first);
				}
				else {
					behind_spaces.push_back(jj.first);
				}
			}
			std::vector<std::vector<Real>> space_min_dist;
			for (size_t j = 0; j < spaces.size(); ++j) {
				std::vector<std::vector<Real>> vars;
				auto& front_sols1 = getMO_HLC().getSubspaceInfo(spaces[j]).m_subspace_front_sol;
				for (size_t i = 0; i < front_sols1.size(); ++i) {
					vars.emplace_back(front_sols1[i]->variable().vect());
				}
				std::vector<Real> temp_min_dist;//�ӿռ�֮�����̾���
				//temp_min_dist.resize(spaces.size());
				Real min_dist = std::numeric_limits<Real>::max();//�ӿռ�֮�����̾����е���С����
				for (size_t k = 0; k < spaces.size(); ++k) {
					if (k <= j) {
						temp_min_dist.push_back(0.);
					}
					else {
						std::vector<std::vector<Real>> vars2;
						auto& front_sols2 = getMO_HLC().getSubspaceInfo(spaces[k]).m_subspace_front_sol;
						for (size_t j = 0; j < front_sols2.size(); ++j) {
							vars2.emplace_back(front_sols2[j]->variable().vect());
						}
						//���������ӿռ�ǰ���е���̾���ֵ
						//��һ��
						dataNormalizeInBound(vars, bound);
						dataNormalizeInBound(vars2, bound);
						Real min_v = std::numeric_limits<Real>::max();
						size_t idx1 = 0, idx2 = 0;
						for (size_t p = 0; p < vars.size(); ++p) {
							for (size_t q = 0; q < vars2.size(); ++q) {
								auto dist = euclideanDistance(vars[p].begin(), vars[p].end(), vars2[q].begin());
								min_v = min_v < dist ? min_v : dist;
								idx1 = p;
								idx2 = q;
							}
						}
						if (min_v < min_dist) {
							min_dist = min_v;
						}
						temp_min_dist.push_back(min_v);
					}
				}
				space_min_dist.emplace_back(temp_min_dist);
			}
			for (size_t j = 0; j < space_min_dist.size(); ++j) {
				for (size_t k = 0; k <= j; ++k) {
					if (k == j) {
						space_min_dist[j][k] = std::numeric_limits<Real>::max();
					}
					else {
						space_min_dist[j][k] = space_min_dist[k][j];
					}
				}
			}
			//�����ӿռ�֮�����С�����γɵıջ�����
			std::vector<size_t> sele_flag(spaces.size(), 0);
			std::vector<std::vector<size_t>> space_clusters;
			while (std::find(sele_flag.begin(), sele_flag.end(), 0) != sele_flag.end()) {
				size_t begin_index = 0;
				for (size_t j = 0; j < sele_flag.size(); ++j) {
					if (sele_flag[j] == 0) {
						begin_index = j;
						break;
					}
				}
				std::vector<size_t> temp_cluster;
				temp_cluster.push_back(spaces[begin_index]);
				sele_flag[begin_index] = 1;
				size_t next_index = std::distance(space_min_dist[begin_index].begin(), std::min_element(space_min_dist[begin_index].begin(), space_min_dist[begin_index].end()));
				while (next_index != begin_index && sele_flag[next_index] == 0) {
					temp_cluster.push_back(spaces[next_index]);
					sele_flag[next_index] = 1;
					next_index = std::distance(space_min_dist[next_index].begin(), std::min_element(space_min_dist[next_index].begin(), space_min_dist[next_index].end()));
				}
				space_clusters.emplace_back(temp_cluster);
			}
			//�����ӿռ���ǰ�ž�
			std::vector<size_t> belong_cluster;
			for (size_t j = 0; j < behind_spaces.size(); ++j) {
				std::vector<std::vector<Real>> vars;
				auto& front_sols1 = getMO_HLC().getSubspaceInfo(behind_spaces[j]).m_subspace_front_sol;
				for (size_t i = 0; i < front_sols1.size(); ++i) {
					vars.emplace_back(front_sols1[i]->variable().vect());
				}
				std::vector<Real> temp_min_dist;//�����ӿռ䵽ǰ����֮�����̾���
				//temp_min_dist.resize(spaces.size());
				Real min_dist = std::numeric_limits<Real>::max();//�ӿռ�֮�����̾����е���С����
				for (size_t k = 0; k < space_clusters.size(); ++k) {
					std::vector<std::vector<Real>> vars2;
					for (size_t p = 0; p < space_clusters[k].size(); ++p) {
						auto& front_sols2 = getMO_HLC().getSubspaceInfo(space_clusters[k][p]).m_subspace_front_sol;
						for (size_t q = 0; q < front_sols2.size(); ++q) {
							vars2.emplace_back(front_sols2[q]->variable().vect());
						}
					}
					//���������ӿռ�ǰ���е���̾���ֵ
					//��һ��
					dataNormalizeInBound(vars, bound);
					dataNormalizeInBound(vars2, bound);
					Real min_v = std::numeric_limits<Real>::max();
					size_t idx1 = 0, idx2 = 0;
					for (size_t p = 0; p < vars.size(); ++p) {
						for (size_t q = 0; q < vars2.size(); ++q) {
							auto dist = euclideanDistance(vars[p].begin(), vars[p].end(), vars2[q].begin());
							min_v = min_v < dist ? min_v : dist;
							idx1 = p;
							idx2 = q;
						}
					}
					if (min_v < min_dist) {
						min_dist = min_v;
					}
					temp_min_dist.push_back(min_v);
				}
				auto index = std::distance(temp_min_dist.begin(), std::min_element(temp_min_dist.begin(), temp_min_dist.end()));
				belong_cluster.push_back(index);
			}
			for (size_t j = 0; j < belong_cluster.size(); ++j) {
				space_clusters[belong_cluster[j]].push_back(behind_spaces[j]);
			}
			m_clusters_in_var = space_clusters;
		}
	}

	void SPMOEA4_9::doubleClusterSelection(size_t select_num, Problem* pro, Random* rnd) {
		//�����ռ����
		clusterInVarspace(pro);
		//Ŀ��ռ����
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t num_obj = pro->numberObjectives();
			std::vector<std::vector<Real>> all_obj_data;
			for (size_t j = 0; j < getPop()[i].size(); ++i) {
				all_obj_data.emplace_back(getPop()[i][j].objective());
			}
			auto normal_obj_data = all_obj_data;
			dataNormalize(normal_obj_data);
			std::vector<Real> min_dist;
			std::vector<std::vector<Real>> all_obj_dist;
			Real max_dist = 0;
			for (size_t k = 0; k < normal_obj_data.size(); ++k) {
				std::vector<Real> temp_dist;
				for (size_t j = 0; j < normal_obj_data.size(); ++j) {
					if (j != k) {
						auto dist = euclideanDistance(normal_obj_data[k].begin(), normal_obj_data[k].end(), normal_obj_data[j].begin());
						temp_dist.push_back(dist);
					}
					else {
						temp_dist.push_back(0.);
					}
				}
				all_obj_dist.emplace_back(temp_dist);
				min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
				auto temp = *std::max_element(temp_dist.begin(), temp_dist.end());
				if (temp > max_dist) {
					max_dist = temp;
				}
			}
			Real sum_dist = 0.;
			for (size_t j = 0; j < min_dist.size(); ++j) {
				sum_dist += min_dist[j];
			}
			std::vector<std::vector<size_t>> ind2cluster;//�����������
			//dbscan����
			Real mean_dist = max_dist / normal_obj_data.size();
			//��ƽ��������Ϊ�ܶȾ���ľ���ֵ
			size_t minPts = 5;//ÿ��������С�ĸ�����
			Real epsilon = 5 * mean_dist;

			DBSCAN dscluster(minPts, epsilon, normal_obj_data);
			dscluster.run();
			std::vector<int> cluster_id;//ÿ������������,��ʵ���δ�����
			for (size_t i = 0; i < dscluster.m_points.size(); ++i) {
				cluster_id.push_back(dscluster.m_points[i]->clusterID);
			}
			std::vector<int> cluster_num;//���е���
			cluster_num.push_back(cluster_id[0]);
			for (size_t j = 1; j < cluster_id.size(); ++j) {
				if (std::find(cluster_num.begin(), cluster_num.end(), cluster_id[j]) == cluster_num.end()) {
					cluster_num.push_back(cluster_id[j]);
				}
			}
			//����ͬ����ĵ����һ��
			for (size_t k = 0; k < cluster_num.size(); ++k) {
				if (cluster_num[k] > 0) {
					std::vector<size_t> temp_cluster;
					for (size_t j = 0; j < cluster_id.size(); ++j) {
						if (cluster_id[j] == cluster_num[k]) {
							temp_cluster.push_back(j);
						}
					}
					ind2cluster.emplace_back(temp_cluster);
				}
			}

		}
	}

	void SPMOEA4_9::assignPop(Population<Solution<>>& pop, Problem* pro, Algorithm* alg, Random* rnd) {
		//�ȵõ���ǰ���ӿռ�
		std::vector<size_t> behind_spaces;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			if (getMO_HLC().getSubspaceInfo(j).m_best_rank > 0) {
				behind_spaces.push_back(j);
			}
		}
		std::vector<Real> box_span;
		for (size_t i = 0; i < behind_spaces.size(); ++i) {
			auto box = getMO_HLC().subspaceTree().getBox(behind_spaces[i]);
			Real v = getMO_HLC().subspaceTree().getBoxVolume(behind_spaces[i]);
			box_span.push_back(v);
		}
		std::vector<Real> space_density;
		for (size_t j = 0; j < behind_spaces.size(); ++j) {
			space_density.push_back(getMO_HLC().getSubspaceInfo(behind_spaces[j]).m_sub_freq / box_span[j]);
		}
		std::vector<Real> space_rank;
		for (size_t j = 0; j < behind_spaces.size(); ++j) {
			space_rank.push_back(1 + getMO_HLC().getSubspaceInfo(behind_spaces[j]).m_best_rank);
		}
		//�ӿռ�����
		std::vector<std::vector<Real>*> objs;
		std::vector<std::vector<Real>> all_metrics;
		for (size_t i = 0; i < space_rank.size(); ++i) {
			std::vector<Real> temp = { space_density[i],space_rank[i] };
			all_metrics.emplace_back(temp);
		}
		for (size_t i = 0; i < all_metrics.size(); ++i) {
			objs.emplace_back(&all_metrics[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		std::vector<size_t> first_layer;//��һ���ӿռ�
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				first_layer.push_back(behind_spaces[i]);
			}
		}
		//���ѡ��һ���ӿռ�
		if (first_layer.size() > 0 && pop.size() > 0) {
			size_t space_idx = (size_t)std::floor(first_layer.size() * rnd->uniform.next());
			auto& his_sols = getMO_HLC().getSubspaceInfo(first_layer[space_idx]).m_history_inds;
			size_t count = 0;
			if (his_sols.size() > pop.size()) {
				std::vector<size_t> sele_inx;
				while (sele_inx.size() < pop.size()) {
					size_t idx = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
					if (sele_inx.empty()) {
						sele_inx.push_back(idx);
					}
					else if (std::find(sele_inx.begin(), sele_inx.end(), idx) == sele_inx.end()) {
						sele_inx.push_back(idx);
					}
				}
				for (size_t i = 0; i < sele_inx.size(); ++i) {
					pop[i].variable() = his_sols[sele_inx[i]]->variable();
					pop[i].objective() = his_sols[sele_inx[i]]->objective();
				}
			}
			else {
				while (count < pop.size()) {
					size_t inx = (size_t)std::floor(behind_spaces.size() * rnd->uniform.next());
					if (getMO_HLC().getSubspaceInfo(behind_spaces[inx]).m_history_inds.size() > 0) {
						auto& his_sol = getMO_HLC().getSubspaceInfo(behind_spaces[inx]).m_history_inds;
						size_t idx = (size_t)std::floor(his_sol.size() * rnd->uniform.next());
						pop[count].variable() = his_sol[idx]->variable();
						pop[count].objective() = his_sol[idx]->objective();
						count++;
					}
				}


			}
		}

		m_immigrant_times += pop.size();
	}

	void SPMOEA4_9::testCoverage(Problem* pro) {
		//����ӿռ��Ƿ������ͬ��PSƬ��
		m_real_front_subspaces.clear();
		m_multi_segment_subspaces.clear();
		m_match_subspace.clear();
		m_error_subspace.clear();
		m_loss_subspace.clear();
		m_match_multi_seg_subspace.clear();
		auto front_sp = getFrontSpace();
		//��̬�ָ���ʹ������Ĵ��룬��ʵPS���ڵ��ӿռ�
		//std::vector<size_t> real_front_subspaces;
		////���ж�ε��ӿռ���Ŀ
		//std::vector<size_t> multi_segment_subspaces;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			auto& box = getMO_HLC().subspaceTree().getBox(i);
			//���ŵ�һά�����
			size_t num = 200;
			std::vector<size_t> flag(num - 1, 1);
			Real delta = (box[0].second - box[0].first) / num;
			for (size_t j = 0; j < num - 1; ++j) {
				std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
				auto sol = CAST_CONOP(pro)->createVar(s);
				//�ж�sol�Ƿ���box��
				for (size_t k = 0; k < box.size(); ++k) {
					if (sol[k] > box[k].second || sol[k] < box[k].first) {
						flag[j] = 0;
						break;
					}
				}
			}
			//ͳ�ƶ���
			size_t segment_time = 0;
			for (size_t j = 1; j < flag.size(); ++j) {
				if (flag[j - 1] == 0 && flag[j] == 1) {
					segment_time++;
				}
			}
			if (flag[0] == 1) {
				segment_time++;
			}

			if (segment_time > 0) {
				m_real_front_subspaces.push_back(i);
				if (segment_time > 1) {
					m_multi_segment_subspaces.push_back(i);
				}
			}
			//std::cout << "the number of segment of the " << i << " th subspace is: " << segment_time << std::endl;
			//1�����ӿռ�ķ����Գ̶�

			//2�����ӿռ��ڵĽ��Ƿ�ɷ�Ϊ����࣬����dbscan����

		}

		////���Լ�������ʶ���׼ȷ��
		//if (m_divide_iteration % m_switch_period >= m_num_push) {
		//	std::cout << "************************************************************** " << std::endl;
		//	std::cout << "neighs match result: " << std::endl;
		//	int total_match_num = 0;
		//	for (size_t i = 0; i < m_recog_neighs.size(); ++i) {
		//		int match_num = 0;
		//		auto center_space = m_recog_neighs[i][0];
		//		if (std::find(m_real_front_subspaces.begin(), m_real_front_subspaces.end(), center_space) != m_real_front_subspaces.end()) {
		//			//����Щλ��ƥ��
		//			std::vector<size_t> match_inx;
		//			for (size_t j = 0; j < m_real_all_front_subspaces.size(); ++j) {
		//				if (m_real_all_front_subspaces[j] == center_space) {
		//					match_inx.push_back(j);
		//				}
		//			}
		//			for (size_t j = 0; j < match_inx.size(); ++j) {
		//				if (match_inx[j] == 0) {
		//					size_t pos1 = m_real_all_front_subspaces[match_inx[j] + 1];
		//					for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
		//						if (m_recog_neighs[i][k] == pos1) {
		//							match_num++;
		//						}
		//					}
		//				}
		//				else if (match_inx[j] == m_real_all_front_subspaces.size() - 1) {
		//					size_t pos1 = m_real_all_front_subspaces[match_inx[j] - 1];
		//					for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
		//						if (m_recog_neighs[i][k] == pos1) {
		//							match_num++;
		//						}
		//					}
		//				}
		//				else {
		//					size_t pos1 = m_real_all_front_subspaces[match_inx[j] - 1];
		//					size_t pos2 = m_real_all_front_subspaces[match_inx[j] + 1];
		//					for (size_t k = 1; k < m_recog_neighs[i].size(); ++k) {
		//						if (m_recog_neighs[i][k] == pos1 || m_recog_neighs[i][k] == pos2) {
		//							match_num++;
		//						}
		//					}
		//				}
		//			}
		//			total_match_num += match_num;
		//		}
		//		else {
		//			match_num = -1;
		//		}
		//		std::cout << match_num << " ";
		//	}
		//	std::cout << std::endl;
		//	//����ܵ�ƥ����
		//	std::cout << "neighs match result: " << std::setprecision(2) << (Real)total_match_num / (2. * m_real_front_subspaces.size() - 2) << std::endl;
		//}

		//��ʵ�ӿռ���ƽ��ӿռ�Ĳ���
		for (size_t i = 0; i < front_sp.size(); ++i) {
			if (std::find(m_real_front_subspaces.begin(), m_real_front_subspaces.end(), front_sp[i]) != m_real_front_subspaces.end()) {
				m_match_subspace.push_back(front_sp[i]);
			}
			else {
				m_error_subspace.push_back(front_sp[i]);
			}
		}
		for (size_t i = 0; i < m_real_front_subspaces.size(); ++i) {
			if (std::find(front_sp.begin(), front_sp.end(), m_real_front_subspaces[i]) == front_sp.end()) {
				m_loss_subspace.push_back(m_real_front_subspaces[i]);
			}
		}
		for (size_t i = 0; i < m_match_subspace.size(); ++i) {
			if (std::find(m_multi_segment_subspaces.begin(), m_multi_segment_subspaces.end(), m_match_subspace[i]) != m_multi_segment_subspaces.end()) {
				m_match_multi_seg_subspace.push_back(m_match_subspace[i]);
			}
		}

		std::cout << "************************************************************** " << std::endl;
		std::cout << "the number of total   subspace is: " << getMO_HLC().numSubspace() << std::endl;
		std::cout << "the number of real-PS subspace is: " << m_real_front_subspaces.size() << std::endl;
		std::cout << "the number of mul-seg subspace is: " << m_multi_segment_subspaces.size() << std::endl;
		std::cout << "the number of matched subspace is: " << m_match_subspace.size() << std::endl;
		std::cout << "match multiseg subspace number is: " << m_match_multi_seg_subspace.size() << std::endl;
		std::cout << "the number of errored subspace is: " << m_error_subspace.size() << std::endl;
		std::cout << "the number of lossed  subspace is: " << m_loss_subspace.size() << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "space sequence: " << std::endl;
		for (size_t i = 0; i < m_real_all_front_subspaces.size(); ++i) {
			std::cout << m_real_all_front_subspaces[i] << " ";
		}
		std::cout << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "error index: " << std::endl;
		for (size_t i = 0; i < m_error_subspace.size(); ++i) {
			std::cout << m_error_subspace[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "frequency  : ";
		for (size_t i = 0; i < m_error_subspace.size(); ++i) {
			size_t fre = getMO_HLC().getSubspaceInfo(m_error_subspace[i]).m_sub_freq;
			std::cout << fre << " ";
		}
		std::cout << std::endl;
		std::cout << "rank       : ";
		for (size_t i = 0; i < m_error_subspace.size(); ++i) {
			size_t rank = getMO_HLC().getSubspaceInfo(m_error_subspace[i]).m_best_rank;
			std::cout << rank << " ";
		}
		std::cout << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "loss index: " << std::endl;
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			std::cout << m_loss_subspace[i] << " ";
		}
		std::cout << std::endl;
		std::cout << "frequency : ";
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			size_t fre = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_sub_freq;
			std::cout << fre << " ";
		}
		std::cout << std::endl;
		std::cout << "rank      : ";
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			size_t rank = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_best_rank;
			std::cout << rank << " ";
		}
		std::cout << std::endl;
		std::cout << "pop span  : ";
		for (size_t i = 0; i < m_loss_subspace.size(); ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(m_loss_subspace[i]).m_history_inds;
			auto& box = getMO_HLC().subspaceTree().getBox(m_loss_subspace[i]);
			std::vector<Real> span;
			if (his_sols.size() > 0) {
				for (size_t k = 0; k < pro->numberVariables(); ++k) {
					Real min_v = (Real)INT16_MAX;
					Real max_v = -1. * INT16_MAX;
					for (size_t j = 0; j < his_sols.size(); ++j) {
						if (his_sols[j]->variable()[k] > max_v) {
							max_v = his_sols[j]->variable()[k];
						}
						if (his_sols[j]->variable()[k] < min_v) {
							min_v = his_sols[j]->variable()[k];
						}
					}
					Real temp = max_v - min_v;
					span.push_back(temp / (box[k].second - box[k].first));
				}
			}
			else {
				for (size_t k = 0; k < pro->numberVariables(); ++k) {
					span.push_back(0.);
				}
			}

			auto max_inx = std::distance(span.begin(), max_element(span.begin(), span.end()));
			auto min_inx = std::distance(span.begin(), min_element(span.begin(), span.end()));
			std::cout << std::setprecision(2) << "(" << span[min_inx] << "," << span[max_inx] << ")";
		}
		std::cout << std::endl;

		std::cout << "************************************************************** " << std::endl;
		std::cout << "the   number   of   immigrant  is: " << m_immigrant_times << std::endl;
		std::cout << std::endl;

		std::cout << std::endl;
	}

	//�ӿռ��ǰ���ӿռ����򽻻��γ�һ���⵱���Ӵ�
	std::vector<std::vector<Real>> SPMOEA4_9::sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
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
			sele_inx = front_link_spaces[0];
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
			size_t first_space_num = 0;
			for (size_t j = 0; j < nei1.size(); ++j) {
				auto& front_sol = getMO_HLC().getSubspaceInfo(nei1[j]).m_front_sol_in_subspace;
				pop_size += front_sol.size();
				if (j == 1) {
					first_space_num = front_sol.size();
				}
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
			std::vector<size_t> candidate_inx;
			if (pop_size < 5) {
				size_t inx = (size_t)std::floor(first_space_num * rnd->uniform.next());
				candidate_inx.push_back(inx);
				auto off1 = sampleInRange(temp_pop[inx], sample_num, pro, alg, rnd);
				temp_pop[candidate_inx[0]].donor().variable() = off1[0];
				temp_pop[candidate_inx[0]].donor().objective() = getInteractiveSols().back().back()->objective();
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//�ӿռ����������
				candidate_inx.push_back((size_t)std::floor(first_space_num * rnd->uniform.next()));
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				size_t inx = 0;
				do {
					inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				} while (inx == candidate_inx[1]);
				candidate_inx.push_back(inx);
				temp_pop[candidate_inx[0]].mutate(rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(candidate_inx[0], rnd, pro);

				temp_pop[candidate_inx[0]].donor().evaluate(pro, this);
				Solution<> ind4(temp_pop[candidate_inx[0]].donor());
				/*ind4.variable() = temp_pop[candidate_inx[0]].trial().variable();
				ind4.objective() = temp_pop[candidate_inx[0]].trial().objective();*/
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].donor().variable().vect());
			}
			m_total_extend_times++;
			//�������Ӵ���֧���ϵ
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].donor().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_extend_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_extend_results[0]++;
			}
			else {
				m_extend_results[2]++;
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleLinearInFrontSpace(size_t space, size_t ind_inx, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		//����ͨ�ռ��ڲ���
		auto& sol = getPop()[0][ind_inx].variable().vect();
		for (size_t i = 0; i < sample_num; ++i) {
			//������ǰ�ظ��幹������Ⱥ
			size_t pop_size = 0;
			auto& front_sol = getMO_HLC().getSubspaceInfo(space).m_front_sol_in_subspace;
			pop_size += front_sol.size();
			PopDE<> temp_pop(pop_size, pro);
			for (size_t k = 0; k < front_sol.size(); ++k) {
				temp_pop[k].variable() = front_sol[k]->variable();
				temp_pop[k].objective() = front_sol[k]->objective();
			}
			std::vector<size_t> candidate_inx;
			if (pop_size < 5) {
				size_t inx = (size_t)std::floor(pop_size * rnd->uniform.next());
				candidate_inx.push_back(inx);
				auto off1 = sampleInRange(temp_pop[inx], sample_num, pro, alg, rnd);
				temp_pop[candidate_inx[0]].donor().variable() = off1[0];
				temp_pop[candidate_inx[0]].donor().objective() = getInteractiveSols().back().back()->objective();
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				//�ӿռ�����������
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				size_t inx = 0;
				do {
					inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				} while (inx == candidate_inx[0]);
				candidate_inx.push_back(inx);

				auto sol1 = temp_pop[candidate_inx[0]].variable().vect();
				auto sol2 = temp_pop[candidate_inx[1]].variable().vect();
				Real rand = rnd->uniform.next();

				std::vector<Real> new_sol;
				for (size_t k = 0; k < sol1.size(); ++k) {
					new_sol.push_back(sol1[k] + rand * (sol2[k] - sol1[k]));
				}
				temp_pop[candidate_inx[0]].donor().variable().vect() = new_sol;
				temp_pop[candidate_inx[0]].donor().evaluate(pro, this);
				Solution<> ind4(temp_pop[candidate_inx[0]].donor());

				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[0]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);

				out_off.emplace_back(temp_pop[candidate_inx[0]].donor().variable().vect());
			}
			m_total_extend_times++;
			//�������Ӵ���֧���ϵ
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].donor().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_extend_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_extend_results[0]++;
			}
			else {
				m_extend_results[2]++;
			}
		}
		return out_off;
	}

	//�����ӿռ�������������������ӿռ�ǰ�ؽ��ƽ�����������
	std::vector<std::vector<Real>> SPMOEA4_9::samplePush(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd) {
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
			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				//size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				temp_pop.append(IndDE(sol));
				candidate_inx.push_back(0);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				temp_pop[candidate_inx[0]].trial().variable().vect() = off[0];
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
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
				temp_pop.append(IndDE(*his_sols[sol_inx]));
				candidate_inx = { 0,1,2 };
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					if (rand < push_prob) {
						if (better_space.empty()) {
							//�ֲ������ӿռ䣬�ڲ�����
							if (subspace_behind_inx.empty()) {//���ɽ���
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
								m_push_behind_empty++;
							}
							else {
								//�ӿռ���ǰ���Ž���
								size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));

								////�ӿռ�����֧���
								//size_t a = m_evaluations;
								//size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								//std::vector<size_t> dominate_inx;
								//auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
								//for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
								//	auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
								//	Dominance ship = objectiveCompare(p2_obj,p1_obj,pro->optimizeMode());
								//	if (ship == Dominance::kDominant) {
								//		dominate_inx.push_back(j);
								//	}
								//}
								//size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
								//temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);
								//temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);

							}
						}
						else {//ѡ������һ�������ӿռ��ǰ�ؽ�
							size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
							//�����ӿռ���֧���

						}
					}
					else {
						if (equal_space.empty()) {
							if (subspace_front_inx.size() < 3) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {
								size_t inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[inx1]]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[inx2]]));
							}
							/*size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
							temp_pop.append(*his_sols[inx1]);
							size_t inx2 = 0;
							do {
								inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
							} while (inx2 == inx1);
							temp_pop.append(*his_sols[inx2]);*/
						}
						else {
							size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
							size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

							size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
						}
					}

				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					Real external_prob = 0.5;
					if (rand < 1 - push_prob) {//ǰ�ƣ���ǰ�Ž⽻��
						size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx1]]));
						size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));

						////�ӿռ�����֧���
						//size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						//std::vector<size_t> dominate_inx;
						//auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
						//for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
						//	auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
						//	Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
						//	if (ship == Dominance::kDominant) {
						//		dominate_inx.push_back(j);
						//	}
						//}
						//size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
						//temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

						//temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
					}
					else {//��������õ��ӿռ佻��
						size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[inx1]));
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(IndDE(*his_sols[inx2]));
						//if (better_space.empty()) {//�ֲ������ӿռ䣬�ڲ�����
						//	candidate_inx = { 0,1,2 };
						//	/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

						//	size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);*/

						//	//�ӿռ�����֧���
						//	size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
						//	std::vector<size_t> dominate_inx;
						//	auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
						//	for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
						//		auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
						//		Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
						//		if (ship == Dominance::kDominant) {
						//			dominate_inx.push_back(j);
						//		}
						//	}
						//	size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]);

						//	temp_pop.append(*his_sols[subspace_behind_inx[sele_inx3]]);
						//}
						//else {//ѡ������һ�������ӿռ��ǰ�ؽ�
						//	candidate_inx = { 0,1,2 };
						//	size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
						//	size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
						//	temp_pop.append(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]);

						//	size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
						//	temp_pop.append(*his_sols[subspace_front_inx[sele_inx2]]);
						//}
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
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
			m_total_push_times++;
			//�������Ӵ���֧���ϵ
			auto& pp = temp_pop[candidate_inx[0]].objective();
			auto& off = temp_pop[candidate_inx[0]].trial().objective();
			Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			if (dominance_ship == Dominance::kDominant) {
				m_push_results[1]++;
			}
			else if (dominance_ship == Dominance::kDominated) {
				m_push_results[0]++;
			}
			else {
				m_push_results[2]++;
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleByIndPos(Solution<>& sol, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto space = getMO_HLC().subspaceTree().getRegionIdx(sol.variable().vect());
		auto& box = getMO_HLC().subspaceTree().getBox(space);
		auto bound = CAST_CONOP(pro)->boundary();
		std::vector<size_t> equal_space;
		auto neighs = getMO_HLC().getSubspaceInfo(space).m_sub_neighbors;
		size_t space_rank = getMO_HLC().getSubspaceInfo(space).m_best_rank;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank == space_rank && space_rank != INT16_MAX) {
				equal_space.push_back(ii);
			}
		}
		std::vector<size_t> better_space;
		for (auto ii : neighs) {
			if (getMO_HLC().getSubspaceInfo(ii).m_best_rank < space_rank) {
				better_space.push_back(ii);
			}
		}
		for (size_t i = 0; i < sample_num; ++i) {
			auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			size_t flag = 0;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				//size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				temp_pop.append(IndDE(sol));
				candidate_inx.push_back(0);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				temp_pop[candidate_inx[0]].trial().variable() = sol.variable();
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
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
				temp_pop.append(IndDE(*his_sols[sol_inx]));
				candidate_inx = { 0,1,2 };
				if (std::find(subspace_front_inx.begin(), subspace_front_inx.end(), sol_inx) != subspace_front_inx.end()) {
					Real rand2 = rnd->uniform.next();
					if (rand < push_prob) {//ǰ��
						if (rand2 < 0.5) {//�ڲ�ǰ��
							flag = 1;
							if (subspace_behind_inx.empty()) {//���ɽ���
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
								m_push_behind_empty++;
							}
							else {
								//�ӿռ���ǰ���Ž���
								/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);
								*/
								//�ӿռ�����֧���
								size_t a = m_evaluations;
								size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								std::vector<size_t> dominate_inx;
								auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
								for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
									auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
									Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
									if (ship == Dominance::kDominant) {
										dominate_inx.push_back(j);
									}
								}
								size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]));

								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx3]]));
							}
						}
						else {//�ⲿǰ��
							flag = 2;
							if (better_space.empty()) {
								//�ֲ������ӿռ䣬�ڲ�����
								if (subspace_behind_inx.empty()) {//���ɽ���
									size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
									temp_pop.append(IndDE(*his_sols[inx1]));
									size_t inx2 = 0;
									do {
										inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
									} while (inx2 == inx1);
									temp_pop.append(IndDE(*his_sols[inx2]));
									m_push_behind_empty++;
								}
								else {
									//�ӿռ���ǰ���Ž���
									/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
									temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

									size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
									temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);
									*/
									//�ӿռ�����֧���
									size_t a = m_evaluations;
									size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
									std::vector<size_t> dominate_inx;
									auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
									for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
										auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
										Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
										if (ship == Dominance::kDominant) {
											dominate_inx.push_back(j);
										}
									}
									size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
									temp_pop.append(IndDE(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]));

									temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx3]]));
								}
							}
							else {//ѡ������һ�������ӿռ��ǰ�ؽ�
								size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_history_inds.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_history_inds[sele_ind_inx]));

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
								//�����ӿռ���֧���

							}
						}
					}
					else {//��չ
						if (rand2 < 0.5) {//�ڲ���չ
							flag = 3;
							if (subspace_front_inx.size() < 3) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {
								size_t inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[inx1]]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[inx2]]));
							}
						}
						else {//�ⲿ��չ
							flag = 4;
							if (equal_space.empty()) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {
								size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
							}
						}
					}
				}
				else {//ѡ���ӿռ���֧�����ǰ�ؽ������
					Real rand2 = rnd->uniform.next();
					if (rand < 1 - push_prob) {//ǰ�ƣ���ǰ�Ž⽻��
						if (rand2 < 0.5) {//�ڲ�ǰ��
							flag = 5;
							size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
							std::vector<size_t> dominate_inx;
							auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
							for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
								auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
								Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
								if (ship == Dominance::kDominant) {
									dominate_inx.push_back(j);
								}
							}
							size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
							temp_pop.append(IndDE(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]));

							temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx3]]));
						}
						else {//�ⲿǰ��
							flag = 6;
							if (better_space.empty()) {//�ֲ������ӿռ䣬�ڲ�����
								/*size_t sele_inx1 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_front_inx[sele_inx1]]);

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(*his_sols[subspace_behind_inx[sele_inx2]]);*/

								//�ӿռ�����֧���
								size_t sele_inx3 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								std::vector<size_t> dominate_inx;
								auto& p1_obj = his_sols[subspace_behind_inx[sele_inx3]]->objective();
								for (size_t j = 0; j < subspace_front_inx.size(); ++j) {
									auto& p2_obj = his_sols[subspace_front_inx[j]]->objective();
									Dominance ship = objectiveCompare(p2_obj, p1_obj, pro->optimizeMode());
									if (ship == Dominance::kDominant) {
										dominate_inx.push_back(j);
									}
								}
								size_t sele_inx2 = (size_t)std::floor(dominate_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[dominate_inx[sele_inx2]]]));

								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx3]]));
							}
							else {//ѡ������һ�������ӿռ��ǰ�ؽ�
								size_t sele_space = better_space[(size_t)std::floor(better_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_subspace_front_sol[sele_ind_inx]));

								size_t sele_inx2 = (size_t)std::floor(subspace_front_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_front_inx[sele_inx2]]));
							}
						}
					}
					else {
						if (rand2 < 0.5) {//�ڲ���չ
							flag = 7;
							if (subspace_behind_inx.size() < 3) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {
								size_t inx1 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[inx1]]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[inx2]]));
							}
						}
						else {//�ⲿ��չ
							flag = 8;
							if (equal_space.empty()) {
								size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[inx1]));
								size_t inx2 = 0;
								do {
									inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
								} while (inx2 == inx1);
								temp_pop.append(IndDE(*his_sols[inx2]));
							}
							else {
								size_t sele_space = equal_space[(size_t)std::floor(equal_space.size() * rnd->uniform.next())];
								size_t sele_ind_inx = (size_t)std::floor(getMO_HLC().getSubspaceInfo(sele_space).m_history_inds.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*getMO_HLC().getSubspaceInfo(sele_space).m_history_inds[sele_ind_inx]));

								size_t sele_inx2 = (size_t)std::floor(subspace_behind_inx.size() * rnd->uniform.next());
								temp_pop.append(IndDE(*his_sols[subspace_behind_inx[sele_inx2]]));
							}
						}
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
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
			if (flag != 0) {
				m_ind_operator_times[flag - 1]++;
				//�������Ӵ���֧���ϵ
				auto& pp = temp_pop[candidate_inx[0]].objective();
				auto& off = temp_pop[candidate_inx[0]].trial().objective();
				Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
				if (dominance_ship == Dominance::kDominant) {
					m_operator_results[flag - 1][1]++;
				}
				else if (dominance_ship == Dominance::kDominated) {
					m_operator_results[flag - 1][0]++;
				}
				else {
					m_operator_results[flag - 1][2]++;
				}
			}

		}
		return out_off;
	}

	//��������ѡ�񣬷�֧��ʱ���ѡ���Ӵ�����ʱ�����ݸ����Ƿ�Ϊǰ�ؾ����Ӵ������ĸ���
	void SPMOEA4_9::SolutionSelection(Problem* pro, Random* rnd) {
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
					m_total_results[0]++;
				}
				else if (dominance_ship == Dominance::kDominated) {//�Ӵ�֧�丸��
					/*if (flag2) {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					else {
						sele_inx.push_back(j);
					}*/
					sele_inx.push_back(j);
					m_total_results[1]++;
				}
				else if (dominance_ship == Dominance::kNonDominated) {
					Real rand = rnd->uniform.next();
					if (rand > 0.5) {
						sele_inx.push_back(j);
					}
					else {
						sele_inx.push_back(getPop()[i].size() + j);
					}
					//auto pp_sol = getPop()[i][j].variable().vect();
					//size_t p_space = getMO_HLC().subspaceTree().getRegionIdx(pp_sol);
					//auto off_sol = getPop()[i].getOffspring()[j].variable().vect();
					//size_t o_space = getMO_HLC().subspaceTree().getRegionIdx(off_sol);
					//Real p_span = 0.;
					//Real o_span = 0.;
					//auto p_box = getMO_HLC().subspaceTree().getBox(p_space);
					//for (size_t j = 0; j < p_box.size(); ++j) {
					//	p_span += (p_box[j].second - p_box[j].first);
					//}
					//p_span /= p_box.size();
					//auto o_box = getMO_HLC().subspaceTree().getBox(o_space);
					//for (size_t j = 0; j < o_box.size(); ++j) {
					//	o_span += (o_box[j].second - o_box[j].first);
					//}
					//o_span /= o_box.size();

					//if (flag1 && flag2) {
					//	//�Ƚ�ǰ�ؽ�ĸ���
					//	Real density1 = getMO_HLC().getSubspaceInfo(p_space).m_front_sol_in_subspace.size() / p_span;
					//	Real density2 = getMO_HLC().getSubspaceInfo(o_space).m_front_sol_in_subspace.size() / o_span;
					//	if (density1 < density2) {
					//		sele_inx.push_back(getPop()[i].size() + j);
					//	}
					//	else {
					//		sele_inx.push_back(j);
					//	}
					//}
					//else if (flag1 && (!flag2)) {
					//	sele_inx.push_back(j);
					//}
					//else if (flag2 && (!flag1)) {
					//	sele_inx.push_back(getPop()[i].size() + j);
					//}
					//else {
					//	/*Real density1 = getMO_HLC().getSubspaceInfo(p_space).m_history_inds.size() / p_span;
					//	Real density2 = getMO_HLC().getSubspaceInfo(o_space).m_history_inds.size() / o_span;
					//	if (density1 < density2) {
					//		sele_inx.push_back(getPop()[i].size() + j);
					//	}
					//	else {
					//		sele_inx.push_back(j);
					//	}*/
					//	size_t rank1 = getMO_HLC().getSubspaceInfo(p_space).m_best_rank;
					//	size_t rank2 = getMO_HLC().getSubspaceInfo(o_space).m_best_rank;
					//	if (rank1 < rank2) {
					//		Real rand = rnd->uniform.next();
					//		if (rand < 0.2) {
					//			sele_inx.push_back(j);
					//		}
					//		else {
					//			sele_inx.push_back(getPop()[i].size() + j);
					//		}
					//		//sele_inx.push_back(getPop()[i].size() + j);
					//	}
					//	else if (rank1 > rank2) {
					//		Real rand = rnd->uniform.next();
					//		if (rand > 0.2) {
					//			sele_inx.push_back(j);
					//		}
					//		else {
					//			sele_inx.push_back(getPop()[i].size() + j);
					//		}
					//		//sele_inx.push_back(j);
					//	}
					//	else {
					//		Real rand = rnd->uniform.next();
					//		if (rand > 0.5) {
					//			sele_inx.push_back(j);
					//		}
					//		else {
					//			sele_inx.push_back(getPop()[i].size() + j);
					//		}
					//	}
					//}
					m_total_results[2]++;
				}
			}
			//ʹ��ѡ�����½���¸����ݻ��켣
			for (size_t j = 0; j < sele_inx.size(); ++j) {
				getPop()[i][j] = getPop()[i].getOffspring()[sele_inx[j]];
			}
		}
	}


	void SPMOEA4_9::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			auto pop_state = getPop()[i].getPopState();
			auto search_box = getPop()[i].getPopSearchBox();
			//��Ⱥ��ʷǰ�ؽ���������ĸ���
			size_t in_count = 0;
			auto& pop_his_front_sols = getPop()[i].getPopHisFrontSols();
			for (size_t j = 0; j < pop_his_front_sols.size(); ++j) {
				if (indInRange(search_box, pop_his_front_sols[j]->variable().vect())) {
					in_count++;
				}
			}
			Real in_ratio = (Real)in_count / pop_his_front_sols.size();
			front_bound_ratio.push_back(in_ratio);
			if (in_ratio < 0.1 && pop_state == "active") {//����
				getPop()[i].setPopState("sleep");
			}
		}
	}

	bool SPMOEA4_9::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA4_9::updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd) {
		//��������Ⱥ����Ϣ����������Ⱥ������������ӿռ�����ֵ
		//�����ӿռ�������Ⱥ���������ڵ����rankֵ
		if (getPop()[pop_inx].getPopState() == "active") {
			Population<Solution<>> temp_off_pop;
			for (size_t j = 0; j < getPop()[pop_inx].getOffspring().size() - getPop()[pop_inx].size(); ++j) {
				temp_off_pop.append(getPop()[pop_inx].getOffspring()[j]);
			}
			getMO_HLC().updateSubspaceInfo(temp_off_pop, pro, rnd);

			std::vector<std::shared_ptr<Solution<>>> temp_pop1;
			for (size_t j = 0; j < getMO_HLC().getCluster(pop_inx)[0].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_represent_sol.size(); ++k) {
					temp_pop1.emplace_back(getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_represent_sol[k]);
				}
			}
			std::vector<std::vector<Real>*> objs;
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				objs.emplace_back(&temp_pop1[j]->objective());
			}
			std::vector<int> rank;
			int num_layer = ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				temp_pop1[j]->setFitness(rank[j]);
			}
			//������Ⱥ�������ӿռ��������ֵ
			for (size_t j = 0; j < getMO_HLC().getCluster(pop_inx)[0].size(); ++j) {
				getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = INT16_MAX;
			}
			for (size_t j = 0; j < temp_pop1.size(); ++j) {
				auto temp_sol = temp_pop1[j]->variable().vect();
				auto inx = getMO_HLC().subspaceTree().getRegionIdx(temp_sol);
				if (temp_pop1[j]->fitness() < getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank) {
					getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank = temp_pop1[j]->fitness();
				}
			}
			//����û�и�����ӿռ䣬��������õ�����ֵ
			std::vector<size_t> cluster_spaces = getMO_HLC().getCluster(pop_inx)[0];
			for (size_t j = 0; j < cluster_spaces.size(); ++j) {
				if (getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_best_rank == INT16_MAX) {
					auto neigh_space = getMO_HLC().getSubspaceInfo(cluster_spaces[j]).m_sub_neighbors;
					std::vector<int> neigh_ranks;
					for (auto iter = neigh_space.begin(); iter != neigh_space.end(); ++iter) {
						size_t inx = *iter;
						if (std::find(cluster_spaces.begin(), cluster_spaces.end(), inx) != cluster_spaces.end()) {
							if (getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank < INT16_MAX) {
								neigh_ranks.push_back(getMO_HLC().getSubspaceInfo(inx).m_sub_best_rank);
							}
						}
					}
					if (neigh_ranks.empty()) {
						getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = num_layer;
					}
					else if (neigh_ranks.size() == 1) {
						getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end()) + 1;
					}
					else {
						int min_rank = *std::min_element(neigh_ranks.begin(), neigh_ranks.end());
						int max_rank = *std::max_element(neigh_ranks.begin(), neigh_ranks.end());
						getMO_HLC().getSubspaceInfo(getMO_HLC().getCluster(pop_inx)[0][j]).m_sub_best_rank = std::ceil((min_rank + max_rank) / 2) + 1;
					}
				}
			}
		}
	}

	void SPMOEA4_9::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//������Ⱥ����ʷǰ�ؽ��ڱ߽��ڵı������������Դ
		size_t max_pop_size = getPopsize();
		//�����ӿռ������ǰ�ؽ�Ŀ��
		auto front_spaces = getFrontSpace();

		assign_pop_resource.push_back(max_pop_size);
		//updatePopResource(assign_pop_resource);
	}

	void SPMOEA4_9::generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, const std::vector<size_t>& pop_resource) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			size_t first_num = 0;
			for (size_t j = 0; j < getPop()[i].size(); ++j) {
				if (getPop()[i][j].fitness() == 0) {
					first_num++;
				}
			}
			if (first_num >= getPop()[i].size()) {
				m_exploit_flag = true;
			}
			Real rand = rnd->uniform.next();
			m_exploit_flag = true;
			if (m_exploit_flag && rand < m_greedy_rate) {
				Real prob_select_var = 1. - std::pow((Real)m_evaluations / m_maximum_evalutions, 0.3);
				prob_select_var = 0.2;
				biasSample(i, pro, alg, rnd, prob_select_var);
				m_local_selection = true;
			}
			else {
				uniformSample(i, pro, alg, rnd);
				m_local_selection = false;
			}

			auto& interactive_sols = getInteractiveSols();
			for (size_t k = 0; k < interactive_sols.size(); ++k) {
				auto& ind = interactive_sols[k].back();
				getPop()[i].getOffspring()[k].variable() = ind->variable();
				getPop()[i].getOffspring()[k].objective() = ind->objective();
				getPop()[i].getOffspring()[k].setTimeEvaluate(0);
			}
		}
	}

	void SPMOEA4_9::uniformSample(size_t pop_inx, Problem* pro, Algorithm* alg, Random* rnd) {
		//���ȷ�����Դ
		auto search_bound = CAST_CONOP(pro)->boundary();
		std::map<size_t, std::vector<size_t>> pop2space;
		for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
			auto var = getPop()[pop_inx][j].variable().vect();
			auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
			if (pop2space[space_idx].empty()) {
				std::vector<size_t> temp_inx;
				pop2space.insert(std::make_pair(space_idx, temp_inx));
			}
			pop2space[space_idx].emplace_back(j);
		}
		Real rand2 = rnd->uniform.next();
		if (rand2 < m_disturbance_rate) {
			//�����Ŷ��ķ�ʽ
			size_t num_var = pro->numberVariables();
			for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
				std::vector<std::pair<Real, Real>> box;
				auto& sol = getPop()[pop_inx][j].variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				auto space_box = getMO_HLC().subspaceTree().getBox(space_inx);
				for (size_t k = 0; k < num_var; ++k) {
					Real span = space_box[k].second - space_box[k].first;
					std::pair<Real, Real> temp1;
					Real lower1 = sol[k] - span / 40;
					temp1.first = lower1 > search_bound[k].first ? lower1 : search_bound[k].first;
					Real upper1 = sol[k] + span / 40;
					temp1.second = upper1 < search_bound[k].second ? upper1 : search_bound[k].second;
					box.emplace_back(temp1);
				}
				sampleRandom(getPop()[pop_inx][j], box, 1, pro, alg, rnd);
			}
		}
		else if (rand2 < m_disturbance_rate + m_subspace_interactive_rate) {
			//�ӿռ�֧���ϵ����
			for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
				auto space = getMO_HLC().subspaceTree().getRegionIdx(getPop()[pop_inx][j].variable().vect());
				auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
				sampleByPush(getPop()[pop_inx][j], his_sols, 1, 0.5, pro, alg, rnd);
			}
		}
		else {
			//����֮�����򽻻�,������ʷ��
			sampleAmongNeighs(1, pop2space, search_bound, getPop()[pop_inx].size(), pro, alg, rnd);
		}
		//if (rand2 < m_subspace_interactive_rate) {
		//	//�ӿռ�֧���ϵ����
		//	for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
		//		auto space = getMO_HLC().subspaceTree().getRegionIdx(getPop()[pop_inx][j].variable().vect());
		//		auto& his_sols = getMO_HLC().getSubspaceInfo(space).m_history_inds;
		//		sampleByPush(getPop()[pop_inx][j], his_sols, 1, 0.5, pro, alg, rnd);
		//	}
		//}
		//else if (rand2 >= m_subspace_interactive_rate && rand2 < m_subspace_interactive_rate + m_neigh_interactive_rate) {
		//	//����֮�����򽻻�,������ʷ��
		//	sampleAmongNeighs(1, pop2space, search_bound, getPop()[pop_inx].size(), pro, alg, rnd);
		//}
		//else {
		//	//����֮��ȫ�ֽ���
		//	for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
		//		Population<Solution<>> temp_pop;
		//		temp_pop.append(getPop()[pop_inx][j]);
		//		size_t sele_inx2 = (size_t)std::floor(getPop()[pop_inx].size() * rnd->uniform.next());
		//		auto& sol2 = getPop()[pop_inx][sele_inx2].variable().vect();
		//		temp_pop.append(getPop()[pop_inx][sele_inx2]);
		//		size_t sele_inx3 = 0;
		//		std::vector<Real> sol3;
		//		do {
		//			sele_inx3 = (size_t)std::floor(getPop()[pop_inx].size() * rnd->uniform.next());
		//			sol3 = getPop()[pop_inx][sele_inx3].variable().vect();
		//		} while (ifSame(sol2, sol3));
		//		temp_pop.append(getPop()[pop_inx][sele_inx3]);
		//		sampleAmongGlobal(temp_pop, search_bound, 1, pro, alg, rnd);
		//	}
		//}
	}

	void SPMOEA4_9::biasSample(size_t pop_inx, Problem* pro, Algorithm* alg, Random* rnd, Real prob) {
		m_exploit_ind_index.clear();
		m_exploit_space_bounds.clear();
		m_exploit_pos.clear();
		auto search_bound = CAST_CONOP(pro)->boundary();
		auto front_spaces = getFrontSpace();
		size_t total_spaces = getMO_HLC().numSubspace();
		auto his_front_sols = getHisFrontSols();
		std::map<size_t, std::vector<size_t>> pop2space;
		//Ǳ��������������Դ��1�������ռ����̽��������2��Ŀ��ռ��ϡ������
		//�ҵ�Ǳ�����򣬷�����Դ��Ŀ��ռ�������ռ��Ǳ�����򣬾ֲ������ڻ���ѡ��

		Real rand = rnd->uniform.next();
		if (rand < prob) {
			//�����ռ��ϡ��������������������ϡ�����
			//�ȵõ��ӿռ串����
			auto select_space = selectSubspaceFromVar(pro, rnd);
			size_t mean_num = getPop()[pop_inx].size() / select_space.size();
			if (mean_num == 0) {
				for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
					auto temp_box = getMO_HLC().subspaceTree().getBox(select_space[j]);
					sampleRandom(getPop()[pop_inx][j], temp_box, 1, pro, alg, rnd);
					m_exploit_space_bounds.emplace_back(temp_box);
				}
			}
			else {
				size_t total_num = mean_num * select_space.size();
				std::vector<size_t> assign_num(select_space.size(), mean_num);
				for (size_t j = 0; j < getPop()[pop_inx].size() - total_num; ++j) {
					assign_num[j]++;
				}
				for (size_t j = 0; j < select_space.size(); ++j) {
					auto temp_box = getMO_HLC().subspaceTree().getBox(select_space[j]);
					sampleRandom(getPop()[pop_inx][j], temp_box, assign_num[j], pro, alg, rnd);
					m_exploit_space_bounds.emplace_back(temp_box);
				}
			}
			/*auto temp_box = getMO_HLC().subspaceTree().getBox(select_space);
			sampleRandom(getPop()[pop_inx][0], temp_box, getPop()[pop_inx].size(), pro, alg, rnd);
			m_exploit_space_bounds.emplace_back(temp_box);*/

			//std::vector<size_t> expliot_spaces;
			//std::vector<Real> box_span;
			//for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//	Real v = getMO_HLC().subspaceTree().getBoxVolume(j);
			//	box_span.push_back(v);
			//}
			//std::vector<Real> space_density;
			//for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//	space_density.push_back(getMO_HLC().getSubspaceInfo(j).m_sub_freq / box_span[j]);
			//}
			//std::vector<Real> space_rank;
			//for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//	space_rank.push_back(getMO_HLC().getSubspaceInfo(j).m_best_rank);
			//}

			////�ӿռ�����
			//std::vector<std::vector<Real>*> objs;
			//std::vector<std::vector<Real>> all_metrics;
			//for (size_t i = 0; i < space_rank.size(); ++i) {
			//	std::vector<Real> temp = { space_density[i],space_rank[i] };
			//	all_metrics.emplace_back(temp);
			//}
			//for (size_t i = 0; i < all_metrics.size(); ++i) {
			//	objs.emplace_back(&all_metrics[i]);
			//}
			//std::vector<int> rank;
			//ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
			////std::vector<size_t> first_layer;//��һ���ӿռ�
			//for (size_t j = 0; j < rank.size(); ++j) {
			//	if (rank[j] == 0) {
			//		expliot_spaces.push_back(j);
			//	}
			//}
			//size_t mean_num = getPop()[pop_inx].size() / expliot_spaces.size();
			//if (mean_num == 0) {
			//	for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
			//		auto temp_box = getMO_HLC().subspaceTree().getBox(expliot_spaces[j]);
			//		sampleRandom(getPop()[pop_inx][j], temp_box, 1,pro,alg,rnd);
			//		m_exploit_space_bounds.emplace_back(temp_box);
			//	}
			//}
			//else {
			//	size_t total_num=mean_num* expliot_spaces.size();
			//	std::vector<size_t> assign_num(expliot_spaces.size(),mean_num);
			//	for (size_t j = 0; j < getPop()[pop_inx].size() - total_num; ++j) {
			//		assign_num[j]++;
			//	}
			//	for (size_t j = 0; j < expliot_spaces.size(); ++j) {
			//		auto temp_box = getMO_HLC().subspaceTree().getBox(expliot_spaces[j]);
			//		sampleRandom(getPop()[pop_inx][j], temp_box, assign_num[j], pro, alg, rnd);
			//		m_exploit_space_bounds.emplace_back(temp_box);
			//	}
			//}
		}
		else {
			//Ŀ��ռ��ϡ���������
			int source = 2;
			size_t his_front_num = his_front_sols.size();
			if (his_front_num > getPop()[pop_inx].size()) {
				source = 2;
			}
			else {
				source = 1;
			}
			auto sol_info = selectSubspaceFromObj(pro, rnd, source);//��Ⱥ��Ŀ��ռ�ֲ�ϡ��ĵ�
			auto max_pair = *std::max_element(sol_info.begin(), sol_info.end(), [](std::pair<Real, std::pair<size_t, size_t>> left, std::pair<Real, std::pair<size_t, size_t>> right) {return left.first < right.first; });
			//���ݵ��ȷ��bound
			auto point_pair = max_pair.second;
			//���ڲ����Ӵ��������ӿռ���ʷ���彻�������Ӵ�
			Real rand2 = rnd->uniform.next();
			size_t num_var = pro->numberVariables();
			std::vector<Real> sol1, sol2;
			size_t index1, index2;
			index1 = point_pair.first;
			index2 = point_pair.second;
			if (source == 1) {
				//��ǰ��Ⱥ
				sol1 = getPop()[pop_inx][index1].variable().vect();
				sol2 = getPop()[pop_inx][index2].variable().vect();
				for (size_t j = 0; j < getPop()[pop_inx].size(); ++j) {
					auto var = getPop()[pop_inx][j].variable().vect();
					auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
					if (pop2space[space_idx].empty()) {
						std::vector<size_t> temp_inx;
						pop2space.insert(std::make_pair(space_idx, temp_inx));
					}
					pop2space[space_idx].emplace_back(j);
				}
			}
			else if (source == 2) {
				//��ʷǰ��
				sol1 = his_front_sols[index1]->variable().vect();
				sol2 = his_front_sols[index2]->variable().vect();
				for (size_t j = 0; j < his_front_sols.size(); ++j) {
					auto var = his_front_sols[j]->variable().vect();
					auto space_idx = getMO_HLC().subspaceTree().getRegionIdx(var);
					if (pop2space[space_idx].empty()) {
						std::vector<size_t> temp_inx;
						pop2space.insert(std::make_pair(space_idx, temp_inx));
					}
					pop2space[space_idx].emplace_back(j);
				}
			}

			//ѡ����Զ�ĵ�Ը�����Ϊ�������ɾֲ���������
			std::vector<std::pair<Real, Real>> box1, box2;
			for (size_t k = 0; k < num_var; ++k) {
				Real span = search_bound[k].second - search_bound[k].first;
				std::pair<Real, Real> temp1, temp2;
				Real lower1 = sol1[k] - span / 100;
				temp1.first = lower1 > search_bound[k].first ? lower1 : search_bound[k].first;
				Real upper1 = sol1[k] + span / 100;
				temp1.second = upper1 < search_bound[k].second ? upper1 : search_bound[k].second;
				box1.emplace_back(temp1);
				Real lower2 = sol2[k] - span / 100;
				temp2.first = lower2 > search_bound[k].first ? lower2 : search_bound[k].first;
				Real upper2 = sol2[k] + span / 100;
				temp2.second = upper2 < search_bound[k].second ? upper2 : search_bound[k].second;
				box2.emplace_back(temp2);
			}
			size_t space1 = getMO_HLC().subspaceTree().getRegionIdx(sol1);
			size_t space2 = getMO_HLC().subspaceTree().getRegionIdx(sol2);
			auto& his_sols1 = getMO_HLC().getSubspaceInfo(space1).m_history_inds;
			auto& his_sols2 = getMO_HLC().getSubspaceInfo(space2).m_history_inds;

			std::vector<size_t> exploit_ind_index = { index1,index2 };
			std::vector<size_t> exploit_space_index = { space1,space2 };
			std::vector<size_t> m_exploit_space_index;
			std::vector<std::vector<Real>> exploit_pos = { sol1,sol2 };
			for (size_t i = 0; i < exploit_ind_index.size(); ++i) {
				if (std::find(m_exploit_ind_index.begin(), m_exploit_ind_index.end(), exploit_ind_index[i]) == m_exploit_ind_index.end()) {
					m_exploit_ind_index.push_back(exploit_ind_index[i]);
					m_exploit_pos.emplace_back(exploit_pos[i]);
				}
				if (std::find(m_exploit_space_index.begin(), m_exploit_space_index.end(), exploit_space_index[i]) == m_exploit_space_index.end()) {
					m_exploit_space_index.push_back(exploit_space_index[i]);
					m_exploit_space_bounds.emplace_back(getMO_HLC().subspaceTree().getBox(exploit_space_index[i]));
				}
			}
			/*box1 = getMO_HLC().subspaceTree().getBox(space1);
			box2 = getMO_HLC().subspaceTree().getBox(space2);*/
			if (rand2 < 0.2) {
				//�ֲ��������
				sampleRandom(getPop()[pop_inx][0], box1, getPop()[pop_inx].size() / 2, pro, alg, rnd);
				sampleRandom(getPop()[pop_inx][1], box2, getPop()[pop_inx].size() / 2, pro, alg, rnd);

				//Population<Solution<>> temp_pop1, temp_pop2;
				//size_t ind_inx1, ind_inx2;
				//for (size_t j = 0; j < his_sols1.size(); ++j) {
				//	temp_pop1.append(*his_sols1[j]);
				//	auto& sol = his_sols1[j]->variable().vect();
				//	if (ifSame(sol, sol1)) {
				//		ind_inx1 = j;
				//	}
				//}
				//if (temp_pop1.size() < m_neighs + 1) {
				//	//������Ⱥ�и���
				//	auto nearest_dist = getPop()[pop_inx].nearestNeighbour(index1, pro, getPop()[pop_inx].size() - 1);
				//	size_t add_num = m_neighs + 1 - temp_pop1.size();
				//	auto space_ind_inx = pop2space[space1];
				//	for (auto jj : nearest_dist) {
				//		if (std::find(space_ind_inx.begin(), space_ind_inx.end(), jj.second) == space_ind_inx.end()) {
				//			temp_pop1.append(getPop()[pop_inx][jj.second]);
				//		}
				//		if (temp_pop1.size() >= m_neighs + 1) {
				//			break;
				//		}
				//	}
				//}
				//for (size_t j = 0; j < his_sols2.size(); ++j) {
				//	temp_pop2.append(*his_sols2[j]);
				//	auto& sol = his_sols2[j]->variable().vect();
				//	if (ifSame(sol, sol2)) {
				//		ind_inx2 = j;
				//	}
				//}
				//if (temp_pop2.size() < m_neighs + 1) {
				//	//������Ⱥ�и���
				//	auto nearest_dist = getPop()[pop_inx].nearestNeighbour(index2, pro, getPop()[pop_inx].size() - 1);
				//	size_t add_num = m_neighs + 1 - temp_pop2.size();
				//	auto space_ind_inx = pop2space[space2];
				//	for (auto jj : nearest_dist) {
				//		if (std::find(space_ind_inx.begin(), space_ind_inx.end(), jj.second) == space_ind_inx.end()) {
				//			temp_pop2.append(getPop()[pop_inx][jj.second]);
				//		}
				//		if (temp_pop2.size() >= m_neighs + 1) {
				//			break;
				//		}
				//	}
				//}
				//sampleAmongNeighs(temp_pop1, ind_inx1, box1, getPop()[pop_inx].size() / 2, pro, alg, rnd);
				//sampleAmongNeighs(temp_pop2, ind_inx2, box2, getPop()[pop_inx].size() / 2, pro, alg, rnd);

			}
			else {
				//�ӿռ佻�������Ӵ�
				if (source == 1) {
					sampleByPush(getPop()[pop_inx][index1], his_sols1, getPop()[pop_inx].size() / 2, 0.5, pro, alg, rnd);
					sampleByPush(getPop()[pop_inx][index2], his_sols2, getPop()[pop_inx].size() / 2, 0.5, pro, alg, rnd);
				}
				else if (source == 2) {
					sampleByPush(*his_front_sols[index1], his_sols1, getPop()[pop_inx].size() / 2, 0.5, pro, alg, rnd);
					sampleByPush(*his_front_sols[index2], his_sols2, getPop()[pop_inx].size() / 2, 0.5, pro, alg, rnd);
				}
			}
		}
	}

	//�����ӿռ�������������������ӿռ�ǰ�ؽ��ƽ�����������
	std::vector<std::vector<Real>> SPMOEA4_9::sampleByPush(Solution<>& sol, const std::vector<std::shared_ptr<Solution<>>>& his_sols, size_t sample_num, Real push_prob, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		for (size_t i = 0; i < sample_num; ++i) {

			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			if (his_sols.size() < 5) {
				//auto off = sampleInRange(sol, 1, pro, alg, rnd);
				//size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
				temp_pop.append(IndDE(sol));
				candidate_inx.push_back(0);
				auto off = sampleBySubspace(sol, 1, pro, alg, rnd);
				temp_pop[candidate_inx[0]].trial().variable().vect() = off[0];
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
				for (auto ii : off) {
					out_off.emplace_back(ii);
				}
			}
			else {//���������ӿռ�λ�ø���ѡ���ƽ�����
				auto& sol1 = sol.variable().vect();
				//����ʷ������ࣺ֧����壬������֧�䣬��֧��
				std::vector<size_t> dominant_inx;
				std::vector<size_t> dominated_inx;
				std::vector<size_t> nondominated_inx;
				auto& obj1 = sol.objective();
				for (size_t j = 0; j < his_sols.size(); ++j) {
					auto& obj2 = his_sols[j]->objective();
					auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
					if (ship == Dominance::kDominant) {
						dominant_inx.push_back(j);
					}
					else if (ship == Dominance::kDominated) {
						dominated_inx.push_back(j);
					}
					else if (ship == Dominance::kNonDominated) {
						nondominated_inx.push_back(j);
					}
				}
				//���ݸ���ѡ��ǰ�ƻ���չ
				temp_pop.append(IndDE(sol));
				candidate_inx = { 0,1,2 };
				Real rand = rnd->uniform.next();
				if (rand < push_prob) {
					//ǰ��
					if (dominated_inx.size() > 0 && dominant_inx.size() > 0) {
						//ǰ�������һ��
						size_t inx1 = (size_t)std::floor(dominant_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[dominant_inx[inx1]]));
						size_t inx2 = (size_t)std::floor(dominated_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[dominated_inx[inx2]]));
					}
					else if (dominant_inx.size() > 0) {
						size_t inx1 = (size_t)std::floor(dominant_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[dominant_inx[inx1]]));
						temp_pop.append(IndDE(sol));
					}
					else if (dominated_inx.size() > 0) {
						temp_pop.append(IndDE(sol));
						size_t inx2 = (size_t)std::floor(dominated_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[dominated_inx[inx2]]));
					}
					else {
						//������ʷ��
						size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[inx1]));
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(IndDE(*his_sols[inx2]));
					}
				}
				else {
					//��չ
					if (nondominated_inx.size() < 3) {
						//������ʷ��
						size_t inx1 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[inx1]));
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(his_sols.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(IndDE(*his_sols[inx2]));
					}
					else {
						//�����֧���
						size_t inx1 = (size_t)std::floor(nondominated_inx.size() * rnd->uniform.next());
						temp_pop.append(IndDE(*his_sols[nondominated_inx[inx1]]));
						size_t inx2 = 0;
						do {
							inx2 = (size_t)std::floor(nondominated_inx.size() * rnd->uniform.next());
						} while (inx2 == inx1);
						temp_pop.append(IndDE(*his_sols[nondominated_inx[inx2]]));
					}
				}

				temp_pop[candidate_inx[0]].mutate(temp_pop.scalingFactor() * (2 * rnd->uniform.next() - 1), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
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
				getInteractiveSols().emplace_back(temp_pair);
				if (ifSame(temp_pop[candidate_inx[0]].variable().vect(), ind4.variable().vect())) {
					size_t a = 1;
				}

				out_off.emplace_back(temp_pop[candidate_inx[0]].trial().variable().vect());
			}
		}
		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleAmongGlobal(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		PopDE<> temp_pop;
		for (size_t i = 0; i < pop.size(); ++i) {
			temp_pop.append(IndDE(pop[i]));
		}
		for (size_t i = 0; i < sample_num; ++i) {
			temp_pop[0].mutate(temp_pop.scalingFactor() * rnd->uniform.next(), &temp_pop[0], &temp_pop[1], &temp_pop[2], pro);
			temp_pop.recombine(0, rnd, pro);

			//�޸���
			SPMOEA::repairSol(temp_pop[0].trial().variable().vect(), bound, rnd);

			temp_pop[0].trial().evaluate(pro, this);

			Solution<> ind4(temp_pop[0]);
			ind4.variable() = temp_pop[0].trial().variable();
			ind4.objective() = temp_pop[0].trial().objective();
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[0]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[1]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[2]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);

			out_off.emplace_back(ind4.variable().vect());
		}

		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleAmongNeighs(size_t method, std::map<size_t, std::vector<size_t>>& pop2space, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto front_spaces = getFrontSpace();
		for (auto& sub : pop2space) {
			size_t idx = sub.first;
			PopDE<> temp_pop;
			std::vector<size_t> candidate_inx;
			std::vector<size_t> act_space = { idx };
			std::vector<size_t> neis;
			if (method == 1) {//��ǰ��Ⱥ����
				for (auto& sub_nei : pop2space) {
					if (sub_nei.first != idx) {
						neis.push_back(sub_nei.first);
					}
				}
			}
			else if (method == 2) {//��ʷ�⽻��
				auto neighs = getMO_HLC().getSubspaceInfo(idx).m_sub_neighbors;
				for (auto& sub_nei : neighs) {
					neis.push_back(sub_nei);
				}
			}

			auto nei_space_prob = spaceLinkProbability(idx, neis);
			std::vector<size_t> pop_index;
			if (nei_space_prob.size() > 0) {
				//ֱ��ȡ���ֵ
				auto max_pair = *std::max_element(nei_space_prob.begin(), nei_space_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
				//�������ӿռ��ڵĸ��彻��
				act_space.push_back(max_pair.first);
			}
			//����ѡ���ӿռ��ڵ�ǰ�ؽ�
			if (method == 1) {
				for (size_t j = 0; j < act_space.size(); ++j) {
					for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
						temp_pop.append(IndDE(getPop()[0][pop2space[act_space[j]][k]]));
					}
				}
			}
			else if (method == 2) {
				for (size_t j = 0; j < act_space.size(); ++j) {
					for (size_t k = 0; k < pop2space[act_space[j]].size(); ++k) {
						temp_pop.append(IndDE(getPop()[0][pop2space[act_space[j]][k]]));
					}
				}
				for (size_t j = 0; j < act_space.size(); ++j) {
					if (std::find(front_spaces.begin(), front_spaces.end(), act_space[j]) == front_spaces.end()) {
						auto& front_sols = getMO_HLC().getSubspaceInfo(act_space[j]).m_subspace_front_sol;
						for (size_t j = 0; j < front_sols.size(); ++j) {
							temp_pop.append(IndDE(*front_sols[j]));
						}
					}
					else {
						auto& front_sols = getMO_HLC().getSubspaceInfo(act_space[j]).m_front_sol_in_subspace;
						for (size_t j = 0; j < front_sols.size(); ++j) {
							temp_pop.append(IndDE(*front_sols[j]));
						}
					}
				}
			}

			for (size_t j = 0; j < sub.second.size(); ++j) {
				candidate_inx.clear();
				candidate_inx.push_back(j);

				size_t sele_inx2 = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				auto& sol2 = temp_pop[sele_inx2].variable().vect();
				candidate_inx.push_back(sele_inx2);
				size_t sele_inx3 = 0;
				std::vector<Real> sol3;
				do {
					sele_inx3 = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
					sol3 = temp_pop[sele_inx3].variable().vect();
				} while (ifSame(sol2, sol3));
				candidate_inx.push_back(sele_inx3);

				temp_pop[j].mutate(temp_pop.scalingFactor() * rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
				temp_pop.recombine(j, rnd, pro);
				temp_pop[j].trial().evaluate(pro, this);

				Solution<> ind4(temp_pop[j]);
				ind4.variable() = temp_pop[j].trial().variable();
				ind4.objective() = temp_pop[j].trial().objective();
				std::vector<std::shared_ptr<Solution<>>> temp_pair;
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[j]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[1]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidate_inx[2]]));
				temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
				getInteractiveSols().emplace_back(temp_pair);

				out_off.emplace_back(ind4.variable().vect());
			}
		}

		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleAmongNeighs(Population<Solution<>>& pop, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		auto neighs = pop.nearestNeighbour(ind_inx, pro, m_neighs);
		std::vector<size_t> inx;
		for (auto& nei : neighs) {
			inx.push_back(nei.second);
		}
		PopDE<> temp_pop;
		temp_pop.append(IndDE(pop[ind_inx]));
		std::vector<size_t> candidate, ridx;
		for (size_t k = 0; k < inx.size(); ++k) {
			temp_pop.append(IndDE(pop[inx[k]]));
		}
		candidate.push_back(0);
		for (size_t k = 0; k < inx.size(); ++k) {
			candidate.push_back(k + 1);
		}
		for (size_t j = 0; j < sample_num; ++j) {
			ridx.clear();
			temp_pop.selectInCandidates(3, candidate, ridx, rnd);
			temp_pop[0].mutate(temp_pop.scalingFactor() * rnd->uniform.next(), &temp_pop[ridx[0]], &temp_pop[ridx[1]], &temp_pop[ridx[2]], pro);
			temp_pop.recombine(0, rnd, pro);

			//�޸���
			SPMOEA::repairSol(temp_pop[0].trial().variable().vect(), bound, rnd);

			temp_pop[0].trial().evaluate(pro, this);

			Solution<> ind4(pop[ind_inx]);
			ind4.variable() = temp_pop[0].trial().variable();
			ind4.objective() = temp_pop[0].trial().objective();
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[0]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);

			out_off.emplace_back(ind4.variable().vect());
		}

		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleAmongNeighs(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		for (size_t j = 0; j < sample_num; ++j) {
			auto neighs = pop.nearestNeighbour(j, pro, m_neighs);
			std::vector<size_t> ind_inx;
			for (auto& nei : neighs) {
				ind_inx.push_back(nei.second);
			}
			PopDE<> temp_pop;
			temp_pop.append(IndDE(pop[j]));
			std::vector<size_t> candidate, ridx;
			for (size_t k = 0; k < ind_inx.size(); ++k) {
				temp_pop.append(IndDE(pop[ind_inx[k]]));
			}
			candidate.push_back(0);
			for (size_t k = 0; k < ind_inx.size(); ++k) {
				candidate.push_back(k + 1);
			}
			temp_pop.selectInCandidates(3, candidate, ridx, rnd);
			temp_pop[0].mutate(temp_pop.scalingFactor(), &temp_pop[ridx[0]], &temp_pop[ridx[1]], &temp_pop[ridx[2]], pro);
			temp_pop.recombine(0, rnd, pro);

			//�޸���
			SPMOEA::repairSol(temp_pop[0].trial().variable().vect(), bound, rnd);

			temp_pop[0].trial().evaluate(pro, this);

			Solution<> ind4(pop[j]);
			ind4.variable() = temp_pop[0].trial().variable();
			ind4.objective() = temp_pop[0].trial().objective();
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[0]));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);

			out_off.emplace_back(ind4.variable().vect());
		}

		return out_off;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleAmongSubspace(Population<Solution<>>& pop, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;

		return out_off;
	}


	size_t SPMOEA4_9::selectSubspace(Problem* pro, Random* rnd) {
		//�ۺ��ӿռ串���ʺ�rankֵ����ѡ��ǰ���ӿռ�
		/*std::vector<Real> volumn;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			Real v = getMO_HLC().subspaceTree().getBoxVolume(j);
			volumn.push_back(v);
		}
		Real min_v = *std::min_element(volumn.begin(), volumn.end());*/
		std::vector<size_t> split_num;
		size_t var_num = pro->numberVariables();
		size_t min_num = 1000;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//size_t num = (size_t)std::floor(min_num*volumn[j]/min_v);

			split_num.push_back(min_num);
		}
		std::vector<Real> space_coverage;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//kd-treeϸ��
			auto box = getMO_HLC().subspaceTree().getBox(j);
			std::vector<Real> ratios(split_num[j], 1. / split_num[j]);
			KDTree temp_tree(ratios, box);
			temp_tree.buildIndex();
			auto& his_sols = getMO_HLC().getSubspaceInfo(j).m_history_inds;
			//������и�����ӿռ���
			std::vector<size_t> have_inds;
			for (size_t k = 0; k < his_sols.size(); ++k) {
				auto& sol = his_sols[k]->variable().vect();
				auto space_idx = temp_tree.getRegionIdx(sol);
				if (have_inds.empty()) {
					have_inds.push_back(space_idx);
				}
				else {
					if (std::find(have_inds.begin(), have_inds.end(), space_idx) == have_inds.end()) {
						have_inds.push_back(space_idx);
					}
				}
			}
			space_coverage.push_back((Real)have_inds.size() / split_num[j]);
		}
		//�ӿռ�����ֵ
		std::vector<Real> space_rank;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			space_rank.push_back(1 + getMO_HLC().getSubspaceInfo(j).m_best_rank);
		}
		//�ӿռ�����
		std::vector<std::vector<Real>*> objs;
		std::vector<std::vector<Real>> all_metrics;
		for (size_t i = 0; i < space_rank.size(); ++i) {
			std::vector<Real> temp = { space_coverage[i],space_rank[i] };
			all_metrics.emplace_back(temp);
		}
		for (size_t i = 0; i < all_metrics.size(); ++i) {
			objs.emplace_back(&all_metrics[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		std::vector<size_t> first_layer;//��һ���ӿռ�
		std::vector<Real> coverage;
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				first_layer.push_back(i);
				coverage.push_back(all_metrics[i][0]);
			}
		}
		//�ӵ�һ����ѡ�񸲸��ʵ͵��ӿռ�
		auto min_cover = std::distance(coverage.begin(), std::min_element(coverage.begin(), coverage.end()));
		//size_t space = (size_t)std::floor(first_layer.size() * rnd->uniform.next());
		return first_layer[min_cover];
	}

	//����Ŀ��ռ�ֲ����ҳ���ʷǰ���е�ϡ����
	std::map<Real, std::pair<size_t, size_t>> SPMOEA4_9::selectSubspaceFromObj(Problem* pro, Random* rnd, int source) {
		//��Ŀ��ռ����з�֧��������
		size_t num_obj = pro->numberObjectives();
		std::vector<std::vector<Real>> all_data;
		if (source == 1) {
			//��ǰ��Ⱥ����
			for (size_t i = 0; i < getPop()[0].size(); ++i) {
				all_data.emplace_back(getPop()[0][i].objective());
			}
		}
		else if (source == 2) {
			//��ʷǰ������
			auto his_front_sols = getHisFrontSols();
			for (size_t i = 0; i < his_front_sols.size(); ++i) {
				all_data.emplace_back(his_front_sols[i]->objective());
			}
		}
		auto normal_data = all_data;
		dataNormalize(normal_data);
		std::vector<Real> min_dist;
		std::vector<std::vector<Real>> all_dist;
		Real max_dist = 0;
		for (size_t i = 0; i < normal_data.size(); ++i) {
			std::vector<Real> temp_dist;
			for (size_t j = 0; j < normal_data.size(); ++j) {
				if (j != i) {
					auto dist = euclideanDistance(normal_data[i].begin(), normal_data[i].end(), normal_data[j].begin());
					temp_dist.push_back(dist);
				}
				else {
					temp_dist.push_back(0.);
				}
			}
			all_dist.emplace_back(temp_dist);
			min_dist.push_back(*std::min_element(temp_dist.begin(), temp_dist.end()));
			auto temp = *std::max_element(temp_dist.begin(), temp_dist.end());
			if (temp > max_dist) {
				max_dist = temp;
			}
		}
		Real sum_dist = 0.;
		for (size_t i = 0; i < min_dist.size(); ++i) {
			sum_dist += min_dist[i];
		}
		int method = 1;
		std::vector<std::vector<size_t>> ind2cluster;//�����������
		if (method == 1) {
			//dbscan����
			Real mean_dist = max_dist / normal_data.size();
			//��ƽ��������Ϊ�ܶȾ���ľ���ֵ
			size_t minPts = 3;//ÿ��������С�ĸ�����
			Real epsilon = 20 * mean_dist;

			DBSCAN dscluster(minPts, epsilon, normal_data);
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
			for (size_t i = 0; i < cluster_num.size(); ++i) {
				if (cluster_num[i] > 0) {
					std::vector<size_t> temp_cluster;
					for (size_t j = 0; j < cluster_id.size(); ++j) {
						if (cluster_id[j] == cluster_num[i]) {
							temp_cluster.push_back(j);
						}
					}
					ind2cluster.emplace_back(temp_cluster);
				}
			}
		}
		else if (method == 2) {
			//k-means����
			size_t K = 2;
			size_t max_iterations = 100;
			KMeans kmeans(K, normal_data.size(), num_obj, max_iterations);
			kmeans.run(normal_data, rnd);
			for (size_t i = 0; i < K; i++) {
				size_t total_points_cluster = kmeans.getClusters()[i].getTotalPoints();
				std::vector<size_t> ind_inx;
				for (size_t j = 0; j < total_points_cluster; j++) {
					ind_inx.push_back(kmeans.getClusters()[i].getPoint(j).getID());
				}
				ind2cluster.emplace_back(ind_inx);
			}
		}

		//�����������
		std::map<Real, std::pair<size_t, size_t>> output_pairs;
		if (ind2cluster.size() < 2) {
			//ѡϡ��㣬�Ƚ�ÿ���������n�����ƽ������
			size_t neigh_num = 5;
			std::vector<Real> mean_dist;
			for (size_t j = 0; j < all_dist.size(); ++j) {
				auto temp = all_dist[j];
				Real sum_d = 0;
				for (size_t k = 0; k < neigh_num; ++k) {
					Real dd = *std::min_element(temp.begin(), temp.end());
					auto inx = std::distance(temp.begin(), std::min_element(temp.begin(), temp.end()));
					sum_d += dd;
					temp[inx] = std::numeric_limits<Real>::max();
				}
				mean_dist.push_back(sum_d / neigh_num);
			}
			auto inx1 = std::distance(mean_dist.begin(), std::max_element(mean_dist.begin(), mean_dist.end()));
			std::pair<size_t, size_t> temp_pair;
			temp_pair.first = inx1;
			mean_dist[inx1] = 0.;
			auto inx2 = std::distance(mean_dist.begin(), std::max_element(mean_dist.begin(), mean_dist.end()));
			temp_pair.second = inx2;
			output_pairs.insert(std::make_pair<>(mean_dist[inx1], temp_pair));

			////���ѡ������
			//size_t inx1 = (size_t)std::floor(ind2cluster[0].size()*rnd->uniform.next());
			//size_t inx2 = 0;
			//do {
			//	inx2 = (size_t)std::floor(ind2cluster[0].size() * rnd->uniform.next());
			//} while (inx1 == inx2);
			//Real dd = all_dist[ind2cluster[0][inx1]][ind2cluster[0][inx2]];
			//std::pair<size_t, size_t> temp_pair;
			//temp_pair.first = ind2cluster[0][inx1];
			//temp_pair.second = ind2cluster[0][inx2];
			//output_pairs.insert(std::make_pair<>(dd,temp_pair));
		}
		else {
			for (size_t i = 0; i < ind2cluster.size(); ++i) {
				auto ind_inx1 = ind2cluster[i];
				std::map <Real, std::pair<size_t, size_t>> cluster_dists;
				for (size_t j = 0; j < ind2cluster.size(); ++j) {
					if (j != i) {
						auto ind_inx2 = ind2cluster[j];
						std::vector<std::vector<Real>> temp_d;
						std::vector<Real> min_d;
						for (size_t k = 0; k < ind_inx1.size(); ++k) {
							std::vector<Real> temp;
							for (size_t p = 0; p < ind_inx2.size(); ++p) {
								temp.push_back(all_dist[ind_inx1[k]][ind_inx2[p]]);
							}
							temp_d.emplace_back(temp);
							auto min_v = *std::min_element(temp.begin(), temp.end());
							min_d.push_back(min_v);
						}
						auto mm = *std::min_element(min_d.begin(), min_d.end());
						size_t inx1 = std::distance(min_d.begin(), std::min_element(min_d.begin(), min_d.end()));
						size_t inx2 = std::distance(temp_d[inx1].begin(), std::min_element(temp_d[inx1].begin(), temp_d[inx1].end()));
						std::pair<size_t, size_t> temp_pair;
						temp_pair.first = ind_inx1[inx1];
						temp_pair.second = ind_inx2[inx2];
						cluster_dists.insert(std::make_pair<>(mm, temp_pair));
					}
				}
				//����С����
				auto min_pair = *std::min_element(cluster_dists.begin(), cluster_dists.end(), [](std::pair<Real, std::pair<size_t, size_t>> left, std::pair<Real, std::pair<size_t, size_t>> right) {return left.first < right.first; });
				output_pairs.insert(min_pair);
			}
		}

		return output_pairs;
	}

	//���������ռ�ֲ����ҳ���ʷǰ���е�ϡ����
	std::vector<size_t> SPMOEA4_9::selectSubspaceFromVar(Problem* pro, Random* rnd) {
		//�ۺ��ӿռ串���ʺ�rankֵ����ѡ��ǰ���ӿռ�
		/*std::vector<Real> volumn;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			Real v = getMO_HLC().subspaceTree().getBoxVolume(j);
			volumn.push_back(v);
		}
		Real min_v = *std::min_element(volumn.begin(), volumn.end());*/
		std::vector<size_t> split_num;
		size_t var_num = pro->numberVariables();
		size_t min_num = 200;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//size_t num = (size_t)std::floor(min_num*volumn[j]/min_v);
			split_num.push_back(min_num);
		}
		std::vector<Real> space_coverage;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			//kd-treeϸ��
			auto box = getMO_HLC().subspaceTree().getBox(j);
			std::vector<Real> ratios(split_num[j], 1. / split_num[j]);
			KDTree temp_tree(ratios, box);
			temp_tree.buildIndex();
			auto& his_sols = getMO_HLC().getSubspaceInfo(j).m_history_inds;
			//������и�����ӿռ���
			std::vector<size_t> have_inds;
			for (size_t k = 0; k < his_sols.size(); ++k) {
				auto& sol = his_sols[k]->variable().vect();
				auto space_idx = temp_tree.getRegionIdx(sol);
				if (have_inds.empty()) {
					have_inds.push_back(space_idx);
				}
				else {
					if (std::find(have_inds.begin(), have_inds.end(), space_idx) == have_inds.end()) {
						have_inds.push_back(space_idx);
					}
				}
			}
			space_coverage.push_back((Real)have_inds.size() / split_num[j]);
		}
		//�ӿռ�����ֵ
		std::vector<Real> space_rank;
		for (size_t j = 0; j < getMO_HLC().numSubspace(); ++j) {
			space_rank.push_back(1 + getMO_HLC().getSubspaceInfo(j).m_best_rank);
		}
		//�ӿռ�����
		std::vector<std::vector<Real>*> objs;
		std::vector<std::vector<Real>> all_metrics;
		for (size_t i = 0; i < space_rank.size(); ++i) {
			std::vector<Real> temp = { space_coverage[i],space_rank[i] };
			all_metrics.emplace_back(temp);
		}
		for (size_t i = 0; i < all_metrics.size(); ++i) {
			objs.emplace_back(&all_metrics[i]);
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		std::vector<size_t> first_layer;//��һ���ӿռ�
		std::vector<Real> coverage;
		for (size_t i = 0; i < rank.size(); ++i) {
			if (rank[i] == 0) {
				first_layer.push_back(i);
				coverage.push_back(all_metrics[i][0]);
			}
		}
		size_t min_rank_space = 0;
		size_t temp_rank = INT16_MAX;
		for (size_t i = 0; i < first_layer.size(); ++i) {
			if (space_rank[first_layer[i]] < temp_rank) {
				temp_rank = space_rank[first_layer[i]];
				min_rank_space = first_layer[i];
			}
		}
		//�ӵ�һ����ѡ�񸲸��ʵ͵��ӿռ�
		auto min_cover_inx = std::distance(coverage.begin(), std::min_element(coverage.begin(), coverage.end()));
		size_t space = (size_t)std::floor(first_layer.size() * rnd->uniform.next());
		//return first_layer[min_cover_inx];
		//return first_layer[space];
		//return min_rank_space;
		return first_layer;
	}

	std::vector<std::vector<Real>> SPMOEA4_9::sampleRandom(Solution<>& ind, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		for (size_t i = 0; i < sample_num; ++i) {
			std::vector<Real> sol;
			for (size_t j = 0; j < bound.size(); ++j) {
				sol.push_back(bound[j].first + rnd->uniform.next() * (bound[j].second - bound[j].first));
			}
			Solution<> ind4(ind);
			ind4.variable().vect() = sol;
			ind4.evaluate(pro, alg);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			getInteractiveSols().emplace_back(temp_pair);

			out_off.emplace_back(sol);
		}
		return out_off;
	}

	bool SPMOEA4_9::spaceSubdivision(Problem* pro, Random* rnd) {
		//��������Ⱥ�����������ڵ�����Ⱥǰ���ӿռ�
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();

		auto front_spaces = getFrontSpace();
		//ϸ���ӿռ�
		for (auto& sp : front_spaces) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
			//size_t min_num = m_divide_granularity;
			size_t min_num = getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity;
			// min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
			//min_num = 20;
			if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
				//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
				int dim = findSplitDim(sp, pro);
				auto& space_bound = getMO_HLC().subspaceTree().getBox(sp);
				Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
				splitSpace(sp, 2, dim, pos, false, pro, rnd);
			}
			else if (getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity < m_divide_granularity) {
				getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity++;
			}
		}
		size_t cur_num_spaces = getMO_HLC().numSubspace();
		bool flag = false;
		if (cur_num_spaces > pre_num_spaces) {
			flag = true;
		}

		//�ӿռ����ռ��
		std::vector<Real> temp_ratio;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
		return flag;
	}

	void SPMOEA4_9::updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd) {
		auto& neigh = getMO_HLC().getSubspaceInfo(inx).m_sub_neighbors;
		std::vector<size_t> link_neighs;
		for (auto& i : neigh) {
			if (std::find(front_spaces.begin(), front_spaces.end(), i) != front_spaces.end()) {
				Real link_probability = ifContinuousSpace(inx, i);
				Real rand = rnd->uniform.next();
				if (rand <= link_probability) {
					link_neighs.push_back(i);
				}
			}
		}
		getMO_HLC().getSubspaceInfo(inx).m_link_subspaces = link_neighs;
	}

	void SPMOEA4_9::updateNeighSpace() {

	}

	void SPMOEA4_9::clusterSubspace() {
		//ÿ�ξ���ֻ����һ�㣬����������ص�
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			getMO_HLC().getSubspaceInfo(i).idx_cluster.clear();
		}
		getMO_HLC().getClusters().clear();
		/*for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
					auto inx = getMO_HLC().getClusters()[i][j][k];
					getMO_HLC().getSubspaceInfo(inx).idx_cluster.clear();
				}
			}
		}*/
		auto& front_clusters = getFrontRegionLinkSpace();

		for (size_t i = 0; i < front_clusters.size(); ++i) {
			std::vector<std::vector<size_t>> clusters;
			clusters.emplace_back(front_clusters[i]);
			getMO_HLC().getClusters().emplace_back(clusters);
		}
		for (size_t i = 0; i < getMO_HLC().getClusters().size(); ++i) {
			for (size_t j = 0; j < getMO_HLC().getClusters()[i].size(); ++j) {
				for (size_t k = 0; k < getMO_HLC().getClusters()[i][j].size(); ++k) {
					auto inx = getMO_HLC().getClusters()[i][j][k];
					getMO_HLC().getSubspaceInfo(inx).idx_cluster.push_back(i);
				}
			}
		}
	}

	std::vector<std::vector<size_t>> SPMOEA4_9::clusterFrontSpace(const std::vector<size_t>& frontspace) {
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
								////�����ӿռ���������ж�
								//if (subspaceLink(inx,frontspace[k])) {
								//	head_cluster.push_back(frontspace[k]);
								//	temp.push_back(frontspace[k]);
								//	select_flag[k] = 1;
								//	count++;
								//}
								head_cluster.push_back(frontspace[k]);
								temp.push_back(frontspace[k]);
								select_flag[k] = 1;
								count++;
							}
						}
					}
				}
				if (temp.empty()) {//û���µ��������
					clustered.emplace_back(head_cluster);
					break;
				}
				else {
					temp_cluster = temp;
				}
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
		}
		getMO_HLC().setFrontClusters(clustered);
		return clustered;
	}

	bool SPMOEA4_9::subspaceLink(size_t inx1, size_t inx2) {
		bool flag = false;
		//�����ӿռ��е�ǰ�ؽ�ķֲ�
		auto ind1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
		auto ind2 = getMO_HLC().getSubspaceInfo(inx2).m_subspace_front_sol;
		////���ӿռ�߽��ƽ������
		auto bound1 = getMO_HLC().subspaceTree().getBox(inx1);
		auto bound2 = getMO_HLC().subspaceTree().getBox(inx2);
		//����������е���̾���ֵ
		Real min_v = INT16_MAX;
		for (size_t i = 0; i < ind1.size(); ++i) {
			for (size_t j = 0; j < ind2.size(); ++j) {
				auto dist = euclideanDistance(ind1[i]->variable().vect().begin(), ind1[i]->variable().vect().end(), ind2[j]->variable().vect().begin());
				min_v = min_v < dist ? min_v : dist;
			}
		}
		//������и��Եĵ�����ֵ
		Real max_dist1 = 0.;
		Real max_dist2 = 0;
		std::vector<Real> dim_span1;
		std::vector<Real> dim_span2;
		//Ѱ��ÿһά������
		for (size_t i = 0; i < bound1.size(); ++i) {
			Real min_v = (Real)INT16_MAX;
			Real max_v = -1. * INT16_MAX;
			for (size_t j = 0; j < ind1.size(); ++j) {
				if (min_v > ind1[j]->variable().vect()[i]) {
					min_v = ind1[j]->variable().vect()[i];
				}
				if (max_v < ind1[j]->variable().vect()[i]) {
					max_v = ind1[j]->variable().vect()[i];
				}
			}
			dim_span1.push_back(max_v - min_v);
		}
		for (size_t i = 0; i < dim_span1.size(); ++i) {
			max_dist1 += std::pow(dim_span1[i], 2);
		}
		for (size_t i = 0; i < bound2.size(); ++i) {
			Real min_v = (Real)INT16_MAX;
			Real max_v = -1. * INT16_MAX;
			for (size_t j = 0; j < ind2.size(); ++j) {
				if (min_v > ind2[j]->variable().vect()[i]) {
					min_v = ind2[j]->variable().vect()[i];
				}
				if (max_v < ind2[j]->variable().vect()[i]) {
					max_v = ind2[j]->variable().vect()[i];
				}
			}
			dim_span2.push_back(max_v - min_v);
		}
		for (size_t i = 0; i < dim_span2.size(); ++i) {
			max_dist2 += std::pow(dim_span2[i], 2);
		}
		max_dist1 = std::sqrt(max_dist1);
		max_dist2 = std::sqrt(max_dist2);
		if (min_v < max_dist1 / 3 || min_v < max_dist2 / 3) {
			flag = true;
		}

		return flag;
	}

	void SPMOEA4_9::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	//ͳ���ӿռ�ǰ�ؽ�Ľ����ɹ���
	void SPMOEA4_9::adaptiveSplit(Problem* pro) {
		for (size_t i = 0; i < getPop().size(); ++i) {
			for (auto space : m_subspace_pop_index) {
				SPMOEA_pop off_pop(0, pro);
				for (size_t k = 0; k < space.second.size(); ++k) {
					off_pop.append(getPop()[i].getOffspring()[space.second[k]]);
				}
				//ͳ���ӿռ佻���ĳɹ���
				auto space_idx = space.first;
				//�Ӵ����ӿռ������
				size_t in_space = 0;
				size_t sol_in_space_front = 0;
				std::queue<int> off_flag;//0:���棬1��ǰ�أ�-1������
				auto& front_sols = getMO_HLC().getSubspaceInfo(space_idx).m_subspace_front_sol;
				for (size_t k = 0; k < off_pop.size(); ++k) {
					auto& sol = off_pop[k].variable().vect();
					auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
					in_space++;
					if (space_inx == space_idx) {
						//�жϽ��ǲ����ӿռ�ǰ��
						bool flag = false;
						for (size_t p = 0; p < front_sols.size(); ++p) {
							auto& sol2 = front_sols[p]->variable().vect();
							if (ifSame(sol, sol2)) {
								sol_in_space_front++;
								flag = true;
								break;
							}
						}
						if (flag) {
							off_flag.push(1);
						}
						else {
							off_flag.push(-1);
						}
					}
					else {
						off_flag.push(0);
					}
				}
				//
				if (auto search = m_subspace_succ_rate.find(space_idx); search != m_subspace_succ_rate.end()) {
					//�����ӿռ�
					m_subspace_succ_rate[space_idx].first += in_space;
					while (!off_flag.empty()) {
						m_subspace_succ_rate[space_idx].second.push(off_flag.front());
						off_flag.pop();
					}
					//���������100���Ӵ����ɽ��
					while (m_subspace_succ_rate[space_idx].second.size() > m_queue_size) {
						m_subspace_succ_rate[space_idx].second.pop();
					}
				}
				else {
					//û���ӿռ�
					std::pair<size_t, std::queue<int>> temp;
					temp.first = in_space;
					temp.second = off_flag;
					m_subspace_succ_rate.insert(std::make_pair<>(space_idx, temp));
				}
			}
		}
	}

	//ͳ���ӿռ�ǰ�ؽ�Ľ����ɹ���
	void SPMOEA4_9::adaptiveSplitFrontSpace(Problem* pro, Random* rnd) {
		std::map<size_t, Real> fail_rate;
		size_t num_var = pro->numberVariables();
		for (auto& space_info : m_subspace_interactive_result) {
			if (getMO_HLC().getSubspaceInfo(space_info.first).m_best_rank == 0) {
				auto que = space_info.second.second;
				size_t fail_times = 0;
				size_t out_times = 0;
				while (!que.empty()) {
					int flag = que.front();
					if (flag == -1) {
						fail_times++;
					}
					if (flag == 0) {
						out_times++;
					}
					que.pop();
				}
				Real fail = (Real)fail_times / ((Real)space_info.second.second.size() + 0.001);
				Real out = (Real)out_times / (Real)space_info.second.second.size();
				fail_rate.insert(std::make_pair<>(space_info.first, fail));
				auto volumn = getMO_HLC().subspaceTree().getBoxVolume(space_info.first);
				if (space_info.second.second.size() >= m_queue_size && fail > 0.4 && volumn / getVarSpaceVolume() > std::pow(0.1, num_var)) {
					//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
					int dim = findSplitDim(space_info.first, pro);
					auto& space_bound = getMO_HLC().subspaceTree().getBox(space_info.first);
					Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
					auto space_front_sols = getMO_HLC().getSubspaceInfo(space_info.first).m_subspace_front_sol;
					splitSpace(space_info.first, 2, dim, pos, false, pro, rnd);
					//���ݷָ���ӿռ�ǰ�ؽ��ռ�ȷ���
					if (space_front_sols.size() == 0) {
						size_t a = 1;
					}
					size_t keep_num = 0;
					for (size_t j = 0; j < space_front_sols.size(); ++j) {
						if (space_front_sols[j]->variable()[dim] < pos) {
							keep_num++;
						}
					}
					keep_num = (size_t)std::floor(m_queue_size * ((Real)keep_num / (space_front_sols.size() + 0.001)));
					size_t sample_num = m_subspace_interactive_result[space_info.first].first;
					m_subspace_interactive_result[space_info.first].first = (size_t)std::floor(sample_num * keep_num / m_queue_size);
					std::queue<int> que2;
					while (m_subspace_interactive_result[space_info.first].second.size() > keep_num) {
						que2.push(m_subspace_interactive_result[space_info.first].second.front());
						m_subspace_interactive_result[space_info.first].second.pop();
					}
					auto num_spaces = getMO_HLC().numSubspace();
					std::pair<size_t, std::queue<int>> temp;
					temp.first = sample_num - m_subspace_interactive_result[space_info.first].first;
					temp.second = que2;
					m_subspace_interactive_result.insert(std::make_pair<>(num_spaces - 1, temp));
				}
			}
		}
		size_t a = 1;
	}

	void SPMOEA4_9::periodSplitSpace(Problem* pro, Random* rnd) {

	}

	//ͳ���ӿռ�ǰ�ؽ�Ľ����ɹ���
	void SPMOEA4_9::adaptiveSplitSpace(Problem* pro, Random* rnd) {
		////test��������
		////std::vector<size_t> residual_space_inx;
		//for (size_t j = 0; j < front_spaces.size(); ++j) {
		//	Real v = getMO_HLC().subspaceTree().getBoxVolume(front_spaces[j]);
		//	if (v / getVarSpaceVolume() > std::pow(0.1, num_var)) {
		//		residual_space_inx.push_back(front_spaces[j]);
		//	}
		//}
		//if (residual_space_inx.size() > 0) {
		//	for (size_t j = 0; j < test_pop.size(); ++j) {
		//		size_t space_idx = (size_t)std::floor(residual_space_inx.size() * rnd->uniform.next());

		//		//ǰ�������ӿռ佻��
		//		auto off = sampleLinearInFrontSpace(residual_space_inx[space_idx], j, 1, pro, alg, rnd);
		//		for (auto ii : off) {
		//			all_off.emplace_back(ii);
		//		}

		//		getPop()[i].getOffspring()[main_pop.size() + j].variable() = getInteractiveSols().back().back()->variable();
		//		getPop()[i].getOffspring()[main_pop.size() + j].objective() = getInteractiveSols().back().back()->objective();

		//		//ͳ���ӿռ佻���ĳɹ���
		//		//�Ӵ����ӿռ������
		//		size_t sol_in_space_front = 0;
		//		int off_flag;//0:���棬1��ǰ�أ�-1������
		//		auto& front_sols = getMO_HLC().getSubspaceInfo(residual_space_inx[space_idx]).m_front_sol_in_subspace;
		//		auto& sol = off[0];
		//		auto obj1 = getInteractiveSols().back().back()->objective();
		//		auto space_index = getMO_HLC().subspaceTree().getRegionIdx(sol);
		//		if (residual_space_inx[space_idx] == space_index) {
		//			//�жϽ��ǲ����ӿռ�ǰ��
		//			bool flag = false;
		//			for (size_t p = 0; p < front_sols.size(); ++p) {
		//				auto& obj2 = front_sols[p]->objective();
		//				auto ship = objectiveCompare(obj1, obj2, pro->optimizeMode());
		//				if (ship == Dominance::kDominated) {
		//					break;
		//				}
		//				else if (p == front_sols.size() - 1) {
		//					sol_in_space_front++;
		//					flag = true;
		//				}
		//			}
		//			if (flag) {
		//				off_flag = 1;
		//			}
		//			else {
		//				off_flag = -1;
		//			}
		//		}
		//		else {
		//			off_flag = 0;
		//		}
		//		//
		//		auto search = m_subspace_interactive_result.find(residual_space_inx[space_idx]);
		//		if (search != m_subspace_interactive_result.end()) {
		//			//�����ӿռ�
		//			m_subspace_interactive_result[residual_space_inx[space_idx]].first += 1;
		//			m_subspace_interactive_result[residual_space_inx[space_idx]].second.push(off_flag);
		//			//�������һ���������Ӵ����ɽ��
		//			while (m_subspace_interactive_result[residual_space_inx[space_idx]].second.size() > m_queue_size) {
		//				m_subspace_interactive_result[residual_space_inx[space_idx]].second.pop();
		//			}
		//		}
		//		else {
		//			//û���ӿռ�
		//			std::pair<size_t, std::queue<int>> temp;
		//			temp.first = 1;
		//			temp.second.push(off_flag);
		//			m_subspace_interactive_result.insert(std::make_pair<>(residual_space_inx[space_idx], temp));
		//		}
		//	}
		//}
		std::map<size_t, Real> fail_rate;
		size_t num_var = pro->numberVariables();
		for (auto& space_info : m_subspace_interactive_result) {
			if (getMO_HLC().getSubspaceInfo(space_info.first).m_best_rank == 0) {
				auto que = space_info.second.second;
				size_t fail_times = 0;
				size_t out_times = 0;
				while (!que.empty()) {
					int flag = que.front();
					if (flag == -1) {
						fail_times++;
					}
					if (flag == 0) {
						out_times++;
					}
					que.pop();
				}
				Real fail = (Real)fail_times / ((Real)space_info.second.second.size() + 0.001);
				Real out = (Real)out_times / (Real)space_info.second.second.size();
				fail_rate.insert(std::make_pair<>(space_info.first, fail));
				auto volumn = getMO_HLC().subspaceTree().getBoxVolume(space_info.first);
				if (space_info.second.second.size() >= m_queue_size && fail > 0.4 && volumn / getVarSpaceVolume() > std::pow(0.1, num_var)) {
					//�ȵõ��ӿռ���ʷ�����ڸ���ϸ�ֵ��ӿռ����Ϣ
					int dim = findSplitDim(space_info.first, pro);
					auto& space_bound = getMO_HLC().subspaceTree().getBox(space_info.first);
					Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
					auto space_front_sols = getMO_HLC().getSubspaceInfo(space_info.first).m_subspace_front_sol;
					splitSpace(space_info.first, 2, dim, pos, false, pro, rnd);
					//���ݷָ���ӿռ�ǰ�ؽ��ռ�ȷ���
					if (space_front_sols.size() == 0) {
						size_t a = 1;
					}
					size_t keep_num = 0;
					for (size_t j = 0; j < space_front_sols.size(); ++j) {
						if (space_front_sols[j]->variable()[dim] < pos) {
							keep_num++;
						}
					}
					keep_num = (size_t)std::floor(m_queue_size * ((Real)keep_num / (space_front_sols.size() + 0.001)));
					size_t sample_num = m_subspace_interactive_result[space_info.first].first;
					m_subspace_interactive_result[space_info.first].first = (size_t)std::floor(sample_num * keep_num / m_queue_size);
					std::queue<int> que2;
					while (m_subspace_interactive_result[space_info.first].second.size() > keep_num) {
						que2.push(m_subspace_interactive_result[space_info.first].second.front());
						m_subspace_interactive_result[space_info.first].second.pop();
					}
					auto num_spaces = getMO_HLC().numSubspace();
					std::pair<size_t, std::queue<int>> temp;
					temp.first = sample_num - m_subspace_interactive_result[space_info.first].first;
					temp.second = que2;
					m_subspace_interactive_result.insert(std::make_pair<>(num_spaces - 1, temp));
				}
			}
		}
		size_t a = 1;
	}

	void SPMOEA4_9::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA4_9::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA4_9::findSplitDim(int inx, Problem* pro) {
		////��������һά��ƽ�����ߵ��ӿռ�ǰ�ؽ�Ĳ�������ȷ�����ֵ�ά��
		//auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
		//auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		//std::vector<Real> front_ind_ratio;
		//for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
		//	Real middle_v = (space_bound[j].first + space_bound[j].second) / 2;
		//	size_t upper_count = 0, lower_count = 0;
		//	for (size_t k = 0; k < front_sols.size(); ++k) {
		//		if (front_sols[k]->variable()[j] >= middle_v) {
		//			upper_count++;
		//		}
		//		else {
		//			lower_count++;
		//		}
		//	}
		//	size_t max_v = upper_count > lower_count ? upper_count : lower_count;
		//	size_t min_v = upper_count < lower_count ? upper_count : lower_count;
		//	if (max_v == 0) {
		//		max_v += 1;
		//	}
		//	front_ind_ratio.push_back((Real)min_v / max_v);
		//}
		//int dim = std::distance(front_ind_ratio.begin(), std::min_element(front_ind_ratio.begin(), front_ind_ratio.end()));

		//ʹ�ÿ��ȷ���ָ��ά��
		auto bound = CAST_CONOP(pro)->boundary();
		auto& space_bound = getMO_HLC().subspaceTree().getBox(inx);
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			//dim_span.push_back(space_bound[j].second - space_bound[j].first);
			//��һ��ά���з�
			dim_span.push_back((space_bound[j].second - space_bound[j].first) / (bound[j].second - bound[j].first));
		}
		int dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));

		////��ά�ȿ������ʱ�����ѡ��һά
		//Real max_span = *std::max_element(dim_span.begin(), dim_span.end());
		//std::vector<size_t> range_span_inx;
		//for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
		//	if (dim_span[j] / max_span >= 0.9) {
		//		range_span_inx.push_back(j);
		//	}
		//}
		//size_t select_inx = (size_t)std::floor(range_span_inx.size()*rnd->uniform.next());
		//int dim = range_span_inx[select_inx];
		return dim;
	}

	void SPMOEA4_9::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i) {
			objs.emplace_back(&pop[i]->objective());
		}
		std::vector<int> rank;
		ofec::nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(m_problem.get())->optimizeMode());
		for (size_t i = 0; i < pop.size(); ++i) {
			pop[i]->setFitness(rank[i]);
		}
	}

	void SPMOEA4_9::recordMetrics(Problem* pro, Algorithm* alg) {
		/************************************/
		/*            ����ָ�����          */
		/************************************/
		//ʹ��archive��������ָ��
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < m_archive.size(); ++i) {
			temp_pop.append(*m_archive[i]);
		}
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < m_archive.size(); ++i) {
			pop_objs.push_back(m_archive[i]->objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);

		//Real temp_IGD = CAST_CONOP(pro)->optima()->invertGenDist(temp_pop);
		getIGD().push_back(temp_IGD);
		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;
		//record();//store metrics data
	}

	void SPMOEA4_9::initiObjSpace(Problem* pro) {

	}
}