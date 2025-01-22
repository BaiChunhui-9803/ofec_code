#include "spmoea2.h"
#include "../../../../../utility/linear_algebra/matrix.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	void SPMOEA2::initialize_() {
		SPMOEA::initialize_();
		m_rl_pop.reset(new SPMOEA_pop(getPopsize(), m_problem.get()));
		size_t num_space = getMO_HLC().numSubspace();
		std::vector<Real> temp_ratio;
		Real total_volume = getVarSpaceVolume();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
		auto& v = *m_param;
		m_divide_granularity = v.get<int>("split granularity");
		m_num_rl = v.get<int>("number of push");
		m_num_extend = v.get<int>("number of extend");
		m_switch_period = m_num_rl + m_num_extend;
	}

	void SPMOEA2::run_() {
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

	void SPMOEA2::record() {
		//std::vector<Real> entry;
		//entry.push_back(m_evaluations);

		//Population<Solution<>> temp_pop;
		//for (size_t i = 0; i < m_archive.size(); ++i) {
		//	temp_pop.append(*m_archive[i]);
		//}
		////IGD
		//Real IGD = m_problem->optimaBase()->invertGenDist(temp_pop);
		//entry.push_back(IGD);
		////HV
		////先得到PF采样点中的最大值作为ref_point
		//auto& refPoint = getMaxRefPoint();
		//Real HV = hypervolumePop(temp_pop, refPoint, m_problem.get(), m_random.get());
		//entry.push_back(HV);
		//dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SPMOEA2::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(3);//第二个为子代
			for (size_t j = 0; j <m_rl_pop->size(); ++j) {
				m_solution[0].push_back(&m_rl_pop->at(j).phenotype());
			}
			for (size_t j = 0; j < m_rl_pop->size(); ++j) {
				m_solution[1].push_back(&m_rl_pop->getOffspring()[j].phenotype());
			}
			for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
				auto& subspace_front_sols = getMO_HLC().getSubspaceInfo(i).m_subspace_front_sol;
				for (size_t j = 0; j < subspace_front_sols.size(); ++j) {
					m_solution.back().push_back(&subspace_front_sols[j]->phenotype());
				}
			}
			/*auto& his_sols = getHisFrontSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_solution.back().push_back(&his_sols[i]->phenotype());
			}*/
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	void SPMOEA2::initiVarSpace(Problem* pro) {
		/*************************************/
		/*      MO_HLC搜索空间划分初始化     */
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

	void SPMOEA2::initPop(Problem* pro, Algorithm* alg, Random* rnd) {
		//在各个子空间生成一个个体
		size_t number_objectives = CAST_CONOP(pro)->numberObjectives();
		size_t num_vars = CAST_CONOP(pro)->numberVariables();
		if (number_objectives == 2) {
			if (num_vars == 2) {
				setArchiveNum(20);
			}
			else if (num_vars == 5) {
				setArchiveNum(50);
			}
			else if (num_vars == 10) {
				setArchiveNum(100);
			}
			else {
				setArchiveNum(100);
			}
		}
		else if (number_objectives == 3) {
			if (num_vars == 3) {
				setArchiveNum(100);
			}
			else if (num_vars == 5) {
				setArchiveNum(100);
			}
			else if (num_vars == 10) {
				setArchiveNum(100);
			}
			else {
				setArchiveNum(200);
			}
		}
		else {
			setArchiveNum(500);
		}
		size_t pop_num = getPopsize();
		size_t space_num = getMO_HLC().numSubspace();//初始化子空间数量
		if (pop_num < space_num) {
			throw("the number of Solutions should more than the number of subspaces");
		}
		//在每个子空间初始化一个个体
		m_rl_pop->initialize(pro,rnd);
		for (size_t i = 0; i < space_num; ++i) {
			auto& box = getMO_HLC().subspaceTree().getBox(i);
			std::vector<Real> temp_sol;
			for (size_t j = 0; j < box.size(); ++j) {
				temp_sol.push_back(box[j].first+rnd->uniform.next()*(box[j].second-box[j].first));
			}
			(*m_rl_pop)[i].variable().vect() = temp_sol;
		}
		m_rl_pop->evaluate(pro,alg);
		Population<Solution<>> new_pop;
		
		auto bound = CAST_CONOP(pro)->boundary();
		std::vector<std::vector<Real>> sols;
		for (size_t j = 0; j < pop_num; ++j) {
			sols.emplace_back(m_rl_pop->at(j).variable().vect());
		}
		SPMOEA::NDSort(*m_rl_pop);
		//初始化子代
		for (size_t j = 0; j < m_rl_pop->size(); ++j) {
			m_rl_pop->getOffspring()[j] = (*m_rl_pop)[j];
			m_rl_pop->getOffspring()[pop_num + j] = (*m_rl_pop)[j];
			new_pop.append((*m_rl_pop)[j]);
		}
		m_rl_pop->setPopState("active");
		//初始化子种群的历史和每代前沿解
		Population<Solution<>> front_pop;
		for (size_t j = 0; j < pop_num; ++j) {
			if (m_rl_pop->at(j).fitness() == 0) {
				front_pop.append((*m_rl_pop)[j]);
				m_rl_pop->getPopHisFrontSols().emplace_back(std::make_shared<Solution<>>((*m_rl_pop)[j]));
			}
		}
		m_rl_pop->getPopGenFrontSols().emplace_back(std::make_shared<Population<Solution<>>>(front_pop));

		//初始化Qtable
		for (size_t i = 0; i < pop_num; ++i) {
			std::map<size_t, std::map<size_t, Real>> temp_Qtable;
			for (size_t i = 0; i < space_num; ++i) {
				std::map<size_t, Real> temp_action;
				for (size_t j = 0; j < space_num; ++j) {
					temp_action.insert(std::make_pair<>(j, 0.));
				}
				temp_Qtable.insert(std::make_pair<>(i, temp_action));
			}
			m_ind_Qtable.emplace_back(std::make_unique<std::map<size_t, std::map<size_t, Real>>>(temp_Qtable));
		}
		//初始化回报值
		for (size_t i = 0; i < pop_num; ++i) {
			m_reward.emplace_back(std::make_unique<std::tuple<size_t,size_t,size_t,Real>>(std::make_tuple<>(0, 0, 0, 0.)));
		}

		//更新历史信息
		SPMOEA::NDSort(new_pop);
		updateHistoryInfo(new_pop, pro);
		SPMOEA::updateObjRange(new_pop, pro);
		SPMOEA::updateNewPop(new_pop);
		SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		//更新子空间信息
		Population<Solution<>> temp_off_pop;
		for (size_t j = 0; j < pop_num; ++j) {
			temp_off_pop.append(m_rl_pop->getOffspring()[j]);
		}
		SPMOEA::updateSubspaceFrontSol(temp_off_pop, pro, rnd);
		updateFrontSpace();
		updateEE(pop_num, 0);

		updateFrontRegionLinkSpace(pro, rnd);
		//根据子空间划分聚类
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
		m_divide_granularity = getSplitGranularity();
		//加入历史演化轨迹
		for (size_t j = 0; j < pop_num; ++j) {
			std::vector<std::shared_ptr<Solution<>>> temp_ind;
			temp_ind.emplace_back(std::make_shared<Solution<>>((*m_rl_pop)[j]));
			getHisEvolveLocus().emplace_back(temp_ind);
		}
		//算子比例
		std::vector<Real> operator_ratio(2, 1.);
		operator_ratio[1] = 0.;
		getOperatorRatio().emplace_back(operator_ratio);
		//个体所在的子空间
		std::vector<size_t> ind_space;
		for (size_t j = 0; j < pop_num; ++j) {
			auto sol = m_rl_pop->at(j).variable().vect();
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

		//初始化选择状态和动作转移概率


		//SPMOEA::recordMetrics(pro, alg);
		SPMOEA::record();
	}

	int SPMOEA2::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		/********************************************************************************
							   根据子种群前沿子空间的体积分配搜索资源
		********************************************************************************/
		//根据算子类型和子连通区域长度分配计算资源
		getInteractiveSols().clear();
		size_t max_eval = maximumEvalutions();
		m_exploration_prob = 1-std::sqrt((Real)m_evaluations / max_eval);
		/********************************************************************************
										强化学习控制交互子空间
		********************************************************************************/
		//1.子代中探索与开发的比例；2.子代中探索与开发的方式
		//auto front_clusters = getFrontRegionLinkSpace();
		//auto front_spaces = getFrontSpace();
		int tag = EvaluationTag::kNormalEval;
		size_t interactive_type;
		if (m_divide_iteration % m_switch_period < m_num_rl) {
			interactive_type=1;
		}
		else {
			interactive_type=2;
		}
		//generateOffspring(pro, alg, rnd, assign_pop_resource, interactive_type);
		tag=generateOffspring(pro, alg, rnd, interactive_type);

		/*Population<Solution<>> offspring_pop;
		int tag = EvaluationTag::kNormalEval;
		for (size_t i = 0; i < m_rl_pop->size(); ++i) {
			for (size_t j = 0; j < getPop()[i].getOffspring().size() - getPop()[i].size(); j++) {
				offspring_pop.append(getPop()[i].getOffspring()[j]);
			}
		}*/
		
		SPMOEA::NDSort(m_rl_pop->getOffspring());
		SPMOEA::updateHistoryInfo(m_rl_pop->getOffspring(), pro);
		updateSubspaceFrontSol(m_rl_pop->getOffspring(), pro, rnd);//子代更新子空间前沿解和代表解
		auto pre_front_space = getFrontSpace();
		updateFrontSpace();//使用子种群子代更新子空间信息
		/**********************************************************************************
					检测是否细分子空间:激活的子种群在其搜索空间的前沿子空间细分
		***********************************************************************************/
		auto pre_num_spaces = getMO_HLC().numSubspace();
		//分析子空间是否为单段来决定是否划分
		auto front_spaces = getFrontSpace();
		size_t switch_period2 = 10;
		//if (m_divide_iteration % switch_period2 == 0) {
		//	for (size_t i = 0; i < front_spaces.size(); ++i) {
		//		Real probability = ifSingleSegSpace(front_spaces[i]);
		//		Real rand = rnd->uniform.next();
		//		if (rand > probability) {//判定为多段
		//			spaceSubdivision(front_spaces[i], pro, rnd);
		//		}
		//	}
		//	updateFrontSpace();
		//	//更新子空间邻域
		//	for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
		//		getMO_HLC().subspaceTree().findNeighbor(i, getMO_HLC().getSubspaceInfo(i).m_sub_neighbors);
		//	}
		//}
		//更新子空间连续性
		auto cur_front_spaces = getFrontSpace();
		for (size_t i = 0; i < cur_front_spaces.size(); ++i) {
			updateLinkSubspace(cur_front_spaces[i], cur_front_spaces, rnd);
		}
		//更新状态空间和动作选择概率


		auto cur_num_spaces = getMO_HLC().numSubspace();
		setAddNumSpace(cur_num_spaces - pre_num_spaces);
		/**********************************************************************************
							   综合各个子种群信息，更新子空间信息
		**********************************************************************************/
		updateFrontRegionLinkSpace(pro, rnd);
		setFrontLastGens(1);
		clusterSubspace();
		SPMOEA::updateNewPop(m_rl_pop->getOffspring());
		SPMOEA::updateObjRange(m_rl_pop->getOffspring(), pro);
		//使用历史所有非支配解更新archive
		//SPMOEA::updateArchive(archiveNum(), pro);
		updateObjSpace();

		//testCoverage(pro);

		Real total_volume = getVarSpaceVolume();
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
										 种群更新，子种群内部进行淘汰选择
		*********************************************************************************************/
		updateInteractiveSols();
		//两种选择模式：1、基于个体迭代的选择；2、
		//采用局部比较的方式，选择个体
		//在进行目标空间排序选择时，考虑解在决策空间的距离
		//multiObjSelection(pro);
		//sparseSelection(pro);
		//localSelection(pro);
		//localCrowdSelection(pro,rnd);
		//localCrowdSelection2(pro, rnd);
		//ensembleSelection(getPop()[0].size(), pro, rnd);
		//在每个linkSpace内结合历史非支配解选择保留的个体
		//个体直接更新，或者根据状态转移概率决定是否更新
		if (m_divide_iteration % m_switch_period < m_num_rl) {
			for (size_t i = 0; i < m_rl_pop->size(); ++i) {
				m_rl_pop->at(i).variable() = m_rl_pop->getOffspring()[i].variable();
				m_rl_pop->at(i).objective() = m_rl_pop->getOffspring()[i].objective();
			}
		}
		
		/*********************************************************************************************
												 记录迭代信息
		**********************************************************************************************/
		//SPMOEA::recordMetrics(pro, alg);
		//SPMOEA::record();
		m_divide_iteration++;
		return tag;
	}

	void SPMOEA2::calculateReward(Problem *pro) {
		auto& his_front = getHisFrontSols();
		//目标归一化
		size_t M = pro->numberObjectives();
		std::vector<std::pair<Real, Real>> bound;
		for (size_t i = 0; i < M; ++i) {
			std::pair<Real, Real> temp(std::numeric_limits<Real>::max(), -1 * std::numeric_limits<Real>::max());
			for (size_t j = 0; j < his_front.size(); ++j) {
				if (temp.first > his_front[j]->objective()[i]) {
					temp.first = his_front[j]->objective()[i];
				}
				if (temp.second < his_front[j]->objective()[i]) {
					temp.second = his_front[j]->objective()[i];
				}
			}
			for (size_t j = 0; j < m_rl_pop->size(); ++j) {
				if (temp.first > m_rl_pop->getOffspring()[j].objective()[i]) {
					temp.first = m_rl_pop->getOffspring()[j].objective()[i];
				}
				if (temp.second < m_rl_pop->getOffspring()[j].objective()[i]) {
					temp.second = m_rl_pop->getOffspring()[j].objective()[i];
				}
			}
			bound.emplace_back(temp);
		}
		std::vector<std::vector<Real>> normalFront;
		std::vector<std::vector<Real>> normalOff;
		for (size_t i = 0; i < his_front.size(); ++i) {
			std::vector<Real> temp;
			for (size_t j = 0; j < M; ++j) {
				temp.push_back((his_front[i]->objective()[j]-bound[j].first) / (bound[j].second-bound[j].first));
			}
			normalFront.emplace_back(temp);
		}
		for (size_t i = 0; i < m_rl_pop->size(); ++i) {
			std::vector<Real> temp;
			for (size_t j = 0; j < M; ++j) {
				temp.push_back((m_rl_pop->getOffspring()[i].objective()[j] - bound[j].first) / (bound[j].second - bound[j].first));
			}
			normalOff.emplace_back(temp);
		}
		//计算个体的被支配数和离前沿最小距离
		auto opt_mode = CAST_CONOP(pro)->optimizeMode();
		std::vector<size_t> dominateNum(m_rl_pop->size(),0);
		for (size_t i = 0; i < normalOff.size(); ++i) {
			for (size_t j = 0; j < normalFront.size(); ++j) {
				Dominance ship = objectiveCompare(normalFront[j],normalOff[i],opt_mode);
				if (ship == Dominance::kDominant) {
					dominateNum[i]++;
				}
			}
		}
		std::vector<Real> minDistFromFront(m_rl_pop->size(), 0.);
		for (size_t i = 0; i < normalOff.size(); ++i) {
			std::vector<Real> all_dist;
			for (size_t j = 0; j < normalFront.size(); ++j) {
				Real dist = euclideanDistance(normalFront[j].begin(),normalFront[j].end(),normalOff[i].begin());
				all_dist.push_back(dist);
			}
			minDistFromFront[i] = *std::min_element(all_dist.begin(),all_dist.end());
		}
		Real max_d = *std::max_element(minDistFromFront.begin(),minDistFromFront.end());
		//根据设定公式计算回报
		for (size_t i = 0; i < m_reward.size(); ++i) {
			std::get<3>(*m_reward[i]) = std::exp(2*(1.-dominateNum[i])+minDistFromFront[i]/max_d);
		}
	}

	void SPMOEA2::updateQtable() {
		//若相同状态采取相同动作得到了不同状态，状态动作值按照常规连续更新
		for (size_t i = 0; i < m_reward.size(); ++i) {
			auto state = std::get<0>(*(m_reward[i]));
			auto action = std::get<1>(*(m_reward[i]));
			auto next_state= std::get<2>(*(m_reward[i]));
			auto reward = std::get<3>(*(m_reward[i]));
			auto maxQ = *std::max_element((*(m_ind_Qtable[i]))[next_state].begin(), (*(m_ind_Qtable[i]))[next_state].end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
			(*(m_ind_Qtable[i]))[state][action] = (*(m_ind_Qtable[i]))[state][action] + m_learning_rate * (reward+m_discount_factor*maxQ.second- (*(m_ind_Qtable[i]))[state][action]);
		}
	}

	void SPMOEA2::testCoverage(Problem* pro) {
		//检测子空间是否包含不同的PS片段
		m_real_front_subspaces.clear();
		m_multi_segment_subspaces.clear();
		m_match_subspace.clear();
		m_error_subspace.clear();
		m_loss_subspace.clear();
		m_match_multi_seg_subspace.clear();
		auto front_sp = getFrontSpace();
		//真实PS所在的子空间
		//std::vector<size_t> real_front_subspaces;
		////具有多段的子空间数目
		//std::vector<size_t> multi_segment_subspaces;
		for (size_t i = 0; i < getMO_HLC().numSubspace(); ++i) {
			auto& box = getMO_HLC().subspaceTree().getBox(i);
			//沿着第一维构造解
			size_t num = 100;
			std::vector<size_t> flag(num - 1, 1);
			Real delta = (box[0].second - box[0].first) / num;
			for (size_t j = 0; j < num - 1; ++j) {
				std::vector<Real> s(CAST_CONOP(pro)->numberVariables(), box[0].first + (j + 1) * delta);
				auto sol = CAST_CONOP(pro)->createVar(s);
				//判断sol是否在box内
				for (size_t k = 0; k < box.size(); ++k) {
					if (sol[k] > box[k].second || sol[k] < box[k].first) {
						flag[j] = 0;
						break;
					}
				}
			}
			//统计段数
			size_t segment_time = 0;
			for (size_t j = 1; j < flag.size(); ++j) {
				if (flag[0] == 0) {
					if (flag[j - 1] == 0 && flag[j] == 1) {
						segment_time++;
					}
				}
				else {
					if (flag[j - 1] == 1 && flag[j] == 0) {
						segment_time++;
					}
				}
			}
			if (segment_time > 0) {
				m_real_front_subspaces.push_back(i);
				if (segment_time > 1) {
					m_multi_segment_subspaces.push_back(i);
				}
			}
			//std::cout << "the number of segment of the " << i << " th subspace is: " << segment_time << std::endl;
			//1、看子空间的非线性程度

			//2、看子空间内的解是否可分为多个类，采用dbscan聚类

		}
		//真实子空间与逼近子空间的差异
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
		std::cout << "the number of total   subspace is: " << getMO_HLC().numSubspace() << std::endl;
		std::cout << "the number of real-PS subspace is: " << m_real_front_subspaces.size() << std::endl;
		std::cout << "the number of mul-seg subspace is: " << m_multi_segment_subspaces.size() << std::endl;
		std::cout << "the number of matched subspace is: " << m_match_subspace.size() << std::endl;
		std::cout << "match multiseg subspace number is: " << m_match_multi_seg_subspace.size() << std::endl;
		std::cout << "the number of errored subspace is: " << m_error_subspace.size() << std::endl;
		std::cout << "the number of lossed  subspace is: " << m_loss_subspace.size() << std::endl;

	}

	void SPMOEA2::updateInteractiveSols() {
		//加入代际演化轨迹
		std::vector<std::pair<std::shared_ptr<Solution<>>, std::shared_ptr<Solution<>>>> temp_ind;
		auto& front_clusters = getFrontRegionLinkSpace();
		auto front_spaces = getFrontSpace();
		auto& interactive_sols = getInteractiveSols();
		for (size_t k = 0; k < interactive_sols.size(); ++k) {
			temp_ind.emplace_back(std::make_pair<>(interactive_sols[k].front(), interactive_sols[k].back()));
		}
		//轨迹解排序
		for (size_t k = 0; k < m_rl_pop->size(); ++k) {
			m_rl_pop->getOffspring()[m_rl_pop->size() + k] = *(interactive_sols[k].front());
		}
		SPMOEA::NDSort(m_rl_pop->getOffspring());
		for (size_t k = 0; k < temp_ind.size(); ++k) {
			temp_ind[k].first->setFitness(m_rl_pop->getOffspring()[k + temp_ind.size()].fitness());
			temp_ind[k].second->setFitness(m_rl_pop->getOffspring()[k].fitness());
		}
		getGenEvolveLocus().emplace_back(temp_ind);
		//搜索空间解的分布
		std::vector<size_t> parents_in_space;
		std::vector<size_t> offs_in_space;
		for (size_t j = 0; j < m_rl_pop->size(); ++j) {
			auto sol = m_rl_pop->at(j).variable().vect();
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			parents_in_space.push_back(space_inx);
		}
		getParentSpaces() = parents_in_space;
		for (size_t j = 0; j < m_rl_pop->size(); ++j) {
			auto sol = m_rl_pop->getOffspring()[j].variable().vect();
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			offs_in_space.push_back(space_inx);
		}
		getOffspringSpaces() = offs_in_space;
		//个体所在的子空间
		std::vector<size_t> ind_space;
		for (size_t j = 0; j < m_rl_pop->size(); ++j) {
			auto sol = m_rl_pop->at(j).variable().vect();
			auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
			if (std::find(ind_space.begin(), ind_space.end(), space_inx) == ind_space.end()) {
				ind_space.push_back(space_inx);
			}
		}
		for (size_t j = 0; j < m_rl_pop->size(); ++j) {
			auto sol = m_rl_pop->getOffspring()[j].variable().vect();
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
		/*for (size_t i = 0; i < m_rl_pop->size(); ++i) {
			
		}*/
	}

	void SPMOEA2::updateLinkSubspace(size_t inx, const std::vector<size_t>& front_spaces, Random* rnd) {
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

	void SPMOEA2::calFrontBoundRatio(std::vector<Real>& front_bound_ratio) {
		for (size_t i = 0; i < getPop().size() - 1; ++i) {
			auto pop_state = getPop()[i].getPopState();
			auto search_box = getPop()[i].getPopSearchBox();
			//种群历史前沿解在搜索域的个数
			size_t in_count = 0;
			auto& pop_his_front_sols = getPop()[i].getPopHisFrontSols();
			for (size_t j = 0; j < pop_his_front_sols.size(); ++j) {
				if (indInRange(search_box, pop_his_front_sols[j]->variable().vect())) {
					in_count++;
				}
			}
			Real in_ratio = (Real)in_count / pop_his_front_sols.size();
			front_bound_ratio.push_back(in_ratio);
			if (in_ratio < 0.1 && pop_state == "active") {//参数
				getPop()[i].setPopState("sleep");
			}
		}
	}

	bool SPMOEA2::indInRange(std::vector<std::pair<Real, Real>>& bound, std::vector<Real>& sol) {
		bool flag = true;
		for (size_t i = 0; i < bound.size(); ++i) {
			if (sol[i]<bound[i].first || sol[i]>bound[i].second) {
				flag = false;
				break;
			}
		}
		return flag;
	}

	void SPMOEA2::updateSubPopSpace(size_t pop_inx, Problem* pro, Random* rnd) {
		//根据新种群的信息，更新子种群所在搜索域的子空间排序值
		//更新子空间在子种群中吸引域内的最好rank值
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
			//更新种群搜索域子空间最好排序值
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
			//对于没有个体的子空间，找邻域最好的排序值
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

	void SPMOEA2::PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro) {
		//根据种群的历史前沿解在边界内的比例分配计算资源
		size_t max_pop_size = getPopsize();
		//updateFrontRegionLinkSpace();
		for (size_t i = 0; i < getPop().size(); ++i) {
			assign_pop_resource.push_back(max_pop_size);
		}
		/*if (m_stage_last_time > 0) {
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
		}*/
		updatePopResource(assign_pop_resource);
	}

	int SPMOEA2::generateOffspring(Problem* pro, Algorithm* alg, Random* rnd, size_t active_type) {
		auto bound = CAST_CONOP(pro)->boundary();
		size_t space_num = getMO_HLC().numSubspace();
		auto front_space = getFrontSpace();
		size_t M = CAST_CONOP(pro)->numberObjectives();
		int tag = EvaluationTag::kNormalEval;
		m_recog_neighs.clear();
		m_all_front_nei_prob.clear();
		std::vector<std::map<size_t, Real>> all_front_nei_prob;
		for (size_t j = 0; j < front_space.size(); ++j) {
			auto space_prob = frontSpaceLinkProbability(front_space[j]);
			all_front_nei_prob.emplace_back(space_prob);
			/***************测试识别的邻域关系是否准确**************/
			std::vector<size_t> recog_neigs = { front_space[j] };
			//找最大的两个
			size_t num_neis = 2;
			auto temp_prob = space_prob;
			if (temp_prob.size() > 0) {
				size_t sele_num = temp_prob.size() > 2 ? num_neis : temp_prob.size();
				for (size_t k = 0; k < sele_num; ++k) {
					auto max_pair = *std::max_element(temp_prob.begin(), temp_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
					recog_neigs.push_back(max_pair.first);
					temp_prob[max_pair.first] = 0.;
				}
				m_recog_neighs.emplace_back(recog_neigs);
			}
		}
		m_all_front_nei_prob = all_front_nei_prob;
		if (active_type == 1) {
			for (size_t i = 0; i < m_rl_pop->size(); ++i) {
				//个体所在子空间
				auto& sol = m_rl_pop->at(i).variable().vect();
				auto space_inx = getMO_HLC().subspaceTree().getRegionIdx(sol);
				//根据状态和Qtable选择动作
				Real rand = rnd->uniform.next();
				size_t interactive_space;
				if (rand < m_exploration_prob || m_divide_iteration == 0) {
					//随机选择动作
					interactive_space = (size_t)std::floor(space_num * rnd->uniform.next());
				}
				else {
					//根据Qtable选择Q值最大的动作
					auto& action = (*(m_ind_Qtable[i]))[space_inx];
					auto action_pair = *std::max_element(action.begin(), action.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
					interactive_space = action_pair.first;
				}
				//根据动作产生新解，子空间内的个体交互
				PopDE<> temp_pop;
				std::vector<size_t> candidates;
				std::vector<std::vector<Real>> off;
				auto& space_front_sols = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
				if (space_inx == interactive_space) {
					if (space_front_sols.size() < 5) {
						//在选出的子空间采样
						off = sampleBySubspace(m_rl_pop->at(i), 1, pro, alg, rnd);
					}
					else {
						for (size_t j = 0; j < 3; ++j) {
							size_t ind_inx = (size_t)std::floor(space_front_sols.size() * rnd->uniform.next());
							temp_pop.append(*space_front_sols[ind_inx]);
							candidates.push_back(j);
						}
					}
				}
				else {
					if (std::find(front_space.begin(), front_space.end(), space_inx) == front_space.end()) {//非前沿子空间
						auto& front_sol1 = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
						size_t ind_inx = (size_t)std::floor(front_sol1.size() * rnd->uniform.next());
						temp_pop.append(*front_sol1[ind_inx]);
					}
					else {
						auto& front_sol1 = getMO_HLC().getSubspaceInfo(space_inx).m_front_sol_in_subspace;
						size_t ind_inx = (size_t)std::floor(front_sol1.size() * rnd->uniform.next());
						temp_pop.append(*front_sol1[ind_inx]);
					}
					candidates.push_back(0);
					if (std::find(front_space.begin(), front_space.end(), interactive_space) == front_space.end()) {//非前沿子空间
						auto& front_sol1 = getMO_HLC().getSubspaceInfo(interactive_space).m_subspace_front_sol;
						size_t ind_inx = (size_t)std::floor(front_sol1.size() * rnd->uniform.next());
						temp_pop.append(*front_sol1[ind_inx]);
					}
					else {
						auto& front_sol1 = getMO_HLC().getSubspaceInfo(interactive_space).m_front_sol_in_subspace;
						size_t ind_inx = (size_t)std::floor(front_sol1.size() * rnd->uniform.next());
						temp_pop.append(*front_sol1[ind_inx]);
					}
					candidates.push_back(1);
					if (std::find(front_space.begin(), front_space.end(), space_inx) == front_space.end()) {//非前沿子空间
						auto& front_sol1 = getMO_HLC().getSubspaceInfo(space_inx).m_subspace_front_sol;
						size_t ind_inx = (size_t)std::floor(front_sol1.size() * rnd->uniform.next());
						temp_pop.append(*front_sol1[ind_inx]);
					}
					else {
						auto& front_sol1 = getMO_HLC().getSubspaceInfo(space_inx).m_front_sol_in_subspace;
						size_t ind_inx = (size_t)std::floor(front_sol1.size() * rnd->uniform.next());
						temp_pop.append(*front_sol1[ind_inx]);
					}
					candidates.push_back(2);
				}
				if (temp_pop.size() >= 3) {
					temp_pop[candidates[0]].mutate(temp_pop.scalingFactor(), &temp_pop[candidates[0]], &temp_pop[candidates[1]], &temp_pop[candidates[2]], pro);
					temp_pop.recombine(candidates[0], rnd, pro);
					temp_pop[candidates[0]].donor().evaluate(pro, this);
					off.emplace_back(temp_pop[candidates[0]].donor().variable().vect());

					Solution<> ind4(m_rl_pop->at(i));
					ind4.variable() = temp_pop[candidates[0]].donor().variable();
					ind4.objective() = temp_pop[candidates[0]].donor().objective();
					std::vector<std::shared_ptr<Solution<>>> temp_pair;
					temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidates[0]]));
					temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidates[1]]));
					temp_pair.emplace_back(std::make_shared<Solution<>>(temp_pop[candidates[2]]));
					temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
					getInteractiveSols().emplace_back(temp_pair);
				}
				//off = sampleByDE(temp_pop, bound, 1, 2, pro, alg, rnd);


				if (off.size() < 1) {
					tag = EvaluationTag::kTerminate;
				}
				//评估新解，计算reward
				auto& ind = *getInteractiveSols().back().back();
				m_rl_pop->getOffspring()[i].variable() = ind.variable();
				m_rl_pop->getOffspring()[i].objective() = ind.objective();
				//父代与子代竞争，胜出的作为下一状态，并更新Qtable
				size_t new_state = getMO_HLC().subspaceTree().getRegionIdx(off[0]);
				//更新即时奖励相关的状态动作
				std::get<0>(*m_reward[i]) = space_inx;
				std::get<1>(*m_reward[i]) = interactive_space;
				std::get<2>(*m_reward[i]) = new_state;
			}
			calculateReward(pro);
			updateQtable();
		}
		else if (active_type == 2) {
			//通过前沿解的最大距离估计密度
			std::vector<Real> space_max_dist;
			for (size_t j = 0; j < front_space.size(); ++j) {
				auto& front_sols = getMO_HLC().getSubspaceInfo(front_space[j]).m_front_sol_in_subspace;
				Real max_dist = 0.;
				for (size_t k = 0; k < front_sols.size(); ++k) {
					auto& sol1 = front_sols[k]->variable().vect();
					for (size_t p = k + 1; p < front_sols.size(); ++p) {
						auto& sol2 = front_sols[p]->variable().vect();
						Real temp_dist = euclideanDistance(sol1.begin(), sol1.end(), sol2.begin());
						if (max_dist < temp_dist) {
							max_dist = temp_dist;
						}
					}
				}
				space_max_dist.push_back(max_dist);
			}
			std::vector<size_t> front_fre;
			for (size_t j = 0; j < front_space.size(); ++j) {
				front_fre.push_back(getMO_HLC().getSubspaceInfo(front_space[j]).m_sub_freq);
			}
			std::vector<Real> space_density;
			for (size_t j = 0; j < front_space.size(); ++j) {
				space_density.push_back(space_max_dist[j] / getMO_HLC().getSubspaceInfo(front_space[j]).m_front_sol_in_subspace.size());
			}
			std::vector<Real> front_ratio;
			for (size_t j = 0; j < front_space.size(); ++j) {
				auto& his_front_sol = getMO_HLC().getSubspaceInfo(front_space[j]).m_front_sol_in_subspace;
				auto& his_sol = getMO_HLC().getSubspaceInfo(front_space[j]).m_front_sol_in_subspace;
				//计算前沿解的比例
				front_ratio.push_back((Real)his_front_sol.size() / his_sol.size());
			}

			for (size_t j = 0; j < m_rl_pop->size(); ++j) {
				size_t space_inx;
				if (j < m_rl_pop->size() / 2) {
					space_inx = std::distance(front_fre.begin(), std::min_element(front_fre.begin(), front_fre.end()));
					front_fre[space_inx]++;
					//space_inx = std::distance(front_ratio.begin(), std::min_element(front_ratio.begin(), front_ratio.end()));
					//space_density[space_inx] = (space_density[space_inx] * box_span[space_inx] + 1) / box_span[space_inx];
					space_inx = std::distance(space_density.begin(), std::max_element(space_density.begin(), space_density.end()));
				}
				else {
					space_inx = (size_t)std::floor(front_space.size() * rnd->uniform.next());
				}
				space_inx = (size_t)std::floor(front_space.size() * rnd->uniform.next());
				// method 1
				/*auto link_neis = getMO_HLC().getSubspaceInfo(front_spaces[space_inx]).m_link_subspaces;
				link_neis.push_back(front_spaces[space_inx]);*/
				//method 2
				//auto space_prob=spaceLinkProbability(front_spaces[space_inx]);
				auto space_prob = all_front_nei_prob[space_inx];
				size_t num_neis = 1;
				std::vector<size_t> neis;
				neis.push_back(front_space[space_inx]);
				if (space_prob.size() > 0) {
					////直接取最大值
					//auto max_pair = *std::max_element(space_prob.begin(), space_prob.end(), [](std::pair<size_t, Real> left, std::pair<size_t, Real> right) {return left.second < right.second; });
					//neis.push_back(max_pair.first);
					//概率选择邻域子空间
					for (size_t k = 0; k < num_neis; ++k) {
						Real rand_num = rnd->uniform.next();
						Real temp_prob = 0.;
						for (auto kk : space_prob) {
							if (rand_num > temp_prob && rand_num < temp_prob + kk.second) {
								if (std::find(neis.begin(), neis.end(), kk.first) == neis.end()) {
									neis.push_back(kk.first);
								}
								break;
							}
							else {
								temp_prob += kk.second;
							}
						}
					}
				}

				auto off = sampleInFrontNeighSpace(neis, j, bound, 1, pro, alg, rnd);
				if (off.size() < 1) {
					tag = EvaluationTag::kTerminate;
				}
				auto& ind = *getInteractiveSols().back().back();
				m_rl_pop->getOffspring()[j].variable() = ind.variable();
				m_rl_pop->getOffspring()[j].objective() = ind.objective();
			}
		}
		
		//updateEE(0, num_exploit);
		return tag;
	}

	//子空间的前沿子空间邻域交互形成一个解当做子代
	std::vector<std::vector<Real>> SPMOEA2::sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd) {
		std::vector<std::vector<Real>> out_off;
		//子连通空间内采样
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
		auto& sol = m_rl_pop->at(ind_inx).variable().vect();
		for (size_t i = 0; i < sample_num; ++i) {
			//选出一个子空间
			size_t sele_inx = 0;
			if (method == 1) {//随机选择
				sele_inx = (size_t)std::floor(front_link_spaces.size() * rnd->uniform.next());
			}
			else if (method == 2) {//根据子空间内的前沿数选择
				//auto inx = std::distance(front_sol_density.begin(), std::min_element(front_sol_density.begin(), front_sol_density.end()));
				auto inx = std::distance(front_sol_num.begin(), std::min_element(front_sol_num.begin(), front_sol_num.end()));
				sele_inx = inx;
				//front_sol_num[inx]++;
				front_sol_density[inx] = front_sol_num[inx] / getMO_HLC().subspaceTree().getBoxVolume(front_link_spaces[inx]);
			}
			sele_inx = front_link_spaces[0];
			//找出其前沿邻域子空间
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
			//将邻域前沿个体构成新种群
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
				temp_pop[candidate_inx[0]].trial().variable() = off1[0];
				temp_pop[candidate_inx[0]].trial().objective() = getInteractiveSols().back().back()->objective();
				for (auto jj : off1) {
					out_off.emplace_back(jj);
					if (ifSame(sol, jj)) {
						size_t a = 1;
					}
				}
			}
			else {
				candidate_inx.push_back((size_t)std::floor(first_space_num * rnd->uniform.next()));
				candidate_inx.push_back((size_t)std::floor(temp_pop.size() * rnd->uniform.next()));
				size_t inx = 0;
				do {
					inx = (size_t)std::floor(temp_pop.size() * rnd->uniform.next());
				} while (inx == candidate_inx[1]);
				candidate_inx.push_back(inx);
				temp_pop[candidate_inx[0]].mutate(rnd->uniform.next(), &temp_pop[candidate_inx[0]], &temp_pop[candidate_inx[1]], &temp_pop[candidate_inx[2]], pro);
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
				getInteractiveSols().emplace_back(temp_pair);
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
			//m_total_extend_times++;
			////父代与子代的支配关系
			//auto& pp = temp_pop[candidate_inx[0]].objective();
			//auto& off = temp_pop[candidate_inx[0]].trial().objective();
			//Dominance dominance_ship = objectiveCompare(pp, off, CAST_CONOP(pro)->optimizeMode());
			//if (dominance_ship == Dominance::kDominant) {
			//	m_extend_results[1]++;
			//}
			//else if (dominance_ship == Dominance::kDominated) {
			//	m_extend_results[0]++;
			//}
			//else {
			//	m_extend_results[2]++;
			//}
		}
		return out_off;
	}

	bool SPMOEA2::spaceSubdivision(Problem* pro, Random* rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		auto front_spaces = getFrontSpace();
		//细分子空间
		for (auto& sp : front_spaces) {
			auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(sp);
			size_t min_num = m_divide_granularity;
			//size_t min_num = getMO_HLC().getSubspaceInfo(sp).m_subspace_granularity;
			// min_num = 50 - 40. / 8 * ((num_var >= 10 ? 10 : num_var) - 2);
			//min_num = 20;
			if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
				//先得到子空间历史解用于更新细分的子空间的信息
				int dim = findSplitDim(sp, pro);
				auto& space_bound = getMO_HLC().subspaceTree().getBox(sp);
				Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
				splitSpace(sp, 2, dim, pos, false, pro, rnd);
			}
		}

		size_t cur_num_spaces = getMO_HLC().numSubspace();
		bool flag = false;
		if (cur_num_spaces > pre_num_spaces) {
			flag = true;
		}

		//子空间体积占比
		std::vector<Real> temp_ratio;
		size_t num_space = getMO_HLC().numSubspace();
		for (size_t i = 0; i < num_space; ++i) {
			temp_ratio.push_back(getMO_HLC().subspaceTree().getBoxVolume(i) / total_volume);
		}
		updateSpaceRatio(temp_ratio);
		return flag;
	}

	bool SPMOEA2::spaceSubdivision(size_t space_inx, Problem* pro, Random* rnd) {
		//划分子种群所在搜索域内的子种群前沿子空间
		Real total_volume = getVarSpaceVolume();
		size_t pre_num_spaces = getMO_HLC().numSubspace();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		//细分子空间
		auto space_volume = getMO_HLC().subspaceTree().getBoxVolume(space_inx);
		size_t min_num = m_divide_granularity;
		if (space_volume / total_volume > std::pow(1. / (Real)min_num, num_var)) {
			//先得到子空间历史解用于更新细分的子空间的信息
			/*int dim = findSplitDim(space_inx, pro);
			auto& space_bound = getMO_HLC().subspaceTree().getBox(space_inx);
			Real pos = (space_bound[dim].first + space_bound[dim].second) / 2;
			splitSpace(space_inx, 2, dim, pos, false, pro, rnd);*/
			auto divide_info = findSplitDimAndPos(space_inx, pro);
			splitSpace(space_inx, 2, divide_info.first, divide_info.second, false, pro, rnd);
		}
		size_t cur_num_spaces = getMO_HLC().numSubspace();
		bool flag = false;
		if (cur_num_spaces > pre_num_spaces) {
			flag = true;
		}
		return flag;
	}

	void SPMOEA2::updateFrontRegionLinkSpace(Problem* pro, Random* rnd) {
		//先找到前沿子空间
		std::vector<size_t> front_space = getFrontSpace();
		//auto link_space = clusterRegionFrontSpace(front_space);
		auto link_space = linearClusterFrontSpace(front_space, pro, rnd);//线性子空间聚类
		getFrontRegionLinkSpace() = link_space;
	}

	std::vector<std::vector<size_t>> SPMOEA2::linearClusterFrontSpace(std::vector<size_t>& frontspace, Problem* pro, Random* rnd) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//标记
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			//首先选择子空间内个体数最多的子空间开始聚
			std::vector<size_t> num_front_sols;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 1) {
					num_front_sols.push_back(0);
				}
				else {
					num_front_sols.push_back(getMO_HLC().getSubspaceInfo(frontspace[i]).m_front_sol_in_subspace.size());
				}
			}
			auto inx = std::distance(num_front_sols.begin(), std::max_element(num_front_sols.begin(), num_front_sols.end()));
			head_cluster.push_back(frontspace[inx]);//先加入一个子空间
			select_flag[inx] = 1;

			for (size_t i = 0; i < select_flag.size(); ++i) {
				count += select_flag[i];
			}
			if (count == select_flag.size()) {
				clustered.emplace_back(head_cluster);
				break;
			}
			if (getMO_HLC().getSubspaceInfo(head_cluster[0]).m_linear_flag == 0) {
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
									//判断是否满足邻域子空间线性关系
									if (ifContinuousSpace(inx, frontspace[k])) {
										head_cluster.push_back(frontspace[k]);
										temp.push_back(frontspace[k]);
										select_flag[k] = 1;
										count++;
									}
								}
							}
						}
					}
					if (temp.empty()) {//没有新的邻域加入
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

	void SPMOEA2::clusterSubspace() {
		//每次聚类只加入一层，允许类间有重叠
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

	std::vector<std::vector<size_t>> SPMOEA2::clusterFrontSpace(const std::vector<size_t>& frontspace) {
		std::vector<std::vector<size_t>> clustered;
		std::vector<size_t> select_flag(frontspace.size(), 0);//标记
		while (std::find(select_flag.begin(), select_flag.end(), 0) != select_flag.end()) {
			size_t begin_space;
			std::vector<size_t> head_cluster;
			size_t count = 0;
			for (size_t i = 0; i < select_flag.size(); ++i) {
				if (select_flag[i] == 0) {
					begin_space = i;
					head_cluster.push_back(frontspace[begin_space]);//先加入一个子空间
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
								////邻域子空间的连续性判断
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
				if (temp.empty()) {//没有新的邻域加入
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

	bool SPMOEA2::subspaceLink(size_t inx1, size_t inx2) {
		bool flag = false;
		//两个子空间中的前沿解的分布
		auto ind1 = getMO_HLC().getSubspaceInfo(inx1).m_subspace_front_sol;
		auto ind2 = getMO_HLC().getSubspaceInfo(inx2).m_subspace_front_sol;
		////到子空间边界的平均距离
		auto bound1 = getMO_HLC().subspaceTree().getBox(inx1);
		auto bound2 = getMO_HLC().subspaceTree().getBox(inx2);
		//计算两类点中的最短距离值
		Real min_v = INT16_MAX;
		for (size_t i = 0; i < ind1.size(); ++i) {
			for (size_t j = 0; j < ind2.size(); ++j) {
				auto dist = euclideanDistance(ind1[i]->variable().vect().begin(), ind1[i]->variable().vect().end(), ind2[j]->variable().vect().begin());
				min_v = min_v < dist ? min_v : dist;
			}
		}
		//两类点中各自的点间最大值
		Real max_dist1 = 0.;
		Real max_dist2 = 0;
		std::vector<Real> dim_span1;
		std::vector<Real> dim_span2;
		//寻找每一维的切面
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

	void SPMOEA2::findClusterCenterSsp() {
		getMO_HLC().findClusterCenterSsp();
	}

	void SPMOEA2::splitSpace(size_t inx, size_t num, int dim, Real pos, bool flag, Problem* pro, Random* rnd) {
		auto his_ind = getMO_HLC().getSubspaceInfo(inx).m_history_inds;
		splitSubspace(inx, num, dim, pos, flag);
		Population<Solution<>> temp_pop;
		for (size_t j = 0; j < his_ind.size(); ++j) {
			temp_pop.append(*his_ind[j]);
		}
		//NDSort(temp_pop);
		SPMOEA::updateSubspaceFrontSol(temp_pop, pro, rnd);
	}

	void SPMOEA2::splitSubspace(size_t inx, size_t num, int dim, Real pos, bool flag) {
		if (flag) {
			SPMOEA::divideSubspace(inx, num);
		}
		else {
			SPMOEA::splitSubspace(inx, dim, pos);
		}
	}

	int SPMOEA2::findSplitDim(int inx, Problem* pro) {
		//根据在哪一维上平分两边的子空间前沿解的差额最大来确定划分的维度
		auto& front_sols = getMO_HLC().getSubspaceInfo(inx).m_subspace_front_sol;
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

		//使用跨度确定分割的维度
		std::vector<Real> dim_span;
		for (size_t j = 0; j < CAST_CONOP(pro)->numberVariables(); ++j) {
			dim_span.push_back(space_bound[j].second - space_bound[j].first);
		}
		dim = std::distance(dim_span.begin(), std::max_element(dim_span.begin(), dim_span.end()));

		return dim;
	}

	void SPMOEA2::NDSort(std::vector<std::shared_ptr<Solution<>>>& pop) {
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

	void SPMOEA2::recordMetrics(Problem* pro, Algorithm* alg) {
		/************************************/
		/*            性能指标计算          */
		/************************************/
		//使用archive计算性能指标
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < m_archive.size(); ++i) {
			temp_pop.append(*m_archive[i]);
		}
		Real temp_IGD = CAST_CONOP(pro)->optima()->invertGenDist(temp_pop);
		getIGD().push_back(temp_IGD);
		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;
		//record();//store metrics data
	}

	void SPMOEA2::initiObjSpace(Problem* pro) {

	}

	std::vector<std::vector<Real>> SPMOEA2::qLearningNonStationaryMDP(const NonStationaryMDP& mdp, Random* rnd,
		Real learning_rate = 0.1,
		Real discount_factor = 0.9,
		Real exploration_prob = 0.2,
		int num_episodes = 1000) {
		int num_states = mdp.getTransitionProbabilities()[0].size();
		int num_actions = mdp.getTransitionProbabilities()[0][0].size();
		int horizon = mdp.getTransitionProbabilities().size();

		std::vector<std::vector<Real>> Q(num_states, std::vector<Real>(num_actions, 0.0));

		for (int episode = 0; episode < num_episodes; ++episode) {
			int state = (int)std::floor(num_states*rnd->uniform.next());

			for (int t = 0; t < horizon; ++t) {
				// 选择动作
				int action;
				if (rnd->uniform.next() < exploration_prob) {
					action = (int)std::floor(num_actions * rnd->uniform.next());
				}
				else {
					auto it = std::max_element(Q[state].begin(), Q[state].end());
					action = std::distance(Q[state].begin(), it);
				}

				// 获取下一个状态和奖励
				const auto& next_state_probabilities = mdp.getTransitionProbabilities()[t][state][action];
				int next_state = generateRandomChoice(next_state_probabilities, rnd);
				Real reward = mdp.getRewardFunction()[t][state][action];

				// 更新Q值
				Q[state][action] = (1 - learning_rate) * Q[state][action] +
					learning_rate * (reward + discount_factor * *std::max_element(Q[next_state].begin(), Q[next_state].end()));

				// 更新当前状态
				state = next_state;
			}
		}

		return Q;
	}

	int SPMOEA2::generateRandomChoice(const std::vector<Real>& probabilities, Random* rnd) {
		Real rand_num = rnd->uniform.next();
		Real cumulative_prob = 0.0;
		for (int i = 0; i < probabilities.size(); ++i) {
			cumulative_prob += probabilities[i];
			if (rand_num <= cumulative_prob) {
				return i;
			}
		}
		return probabilities.size() - 1;
	}

	void NonStationaryMDP::generateNonStationaryTransitions() {
		// 生成随时间变化的状态转移概率
		//transition_probabilities.resize(horizon, std::vector<std::vector<std::vector<Real>>>(num_states, std::vector<std::vector<Real>>(num_actions, std::vector<Real>(num_states, 0.0))));
		action_selection_probabilities.resize(horizon, std::vector<std::vector<Real>>(num_states, std::vector<Real>(num_actions, 0.)));
		for (int t = 0; t < horizon; ++t) {
			for (int s = 0; s < num_states; ++s) {
				for (int a = 0; a < num_actions; ++a) {
					// 在这里可以添加对状态转移概率的具体变化规律
					//定义动作选择概率
					//transition_probabilities[t][s][a] = generateDirichletDistribution(num_states);
					action_selection_probabilities[t][s][a] = 1./num_states;
				}
			}
		}
	}

	void NonStationaryMDP::generateNonStationaryRewards() {
		// 生成随时间变化的即时奖励
		reward_function.resize(horizon, std::vector<std::vector<Real>>(num_states, std::vector<Real>(num_actions, 0.0)));
		for (int t = 0; t < horizon; ++t) {
			for (int s = 0; s < num_states; ++s) {
				for (int a = 0; a < num_actions; ++a) {
					// 在这里可以添加对奖励函数的具体变化规律
					//reward_function[t][s][a] = generateNormalDistribution(0, 1);
					reward_function[t][s][a] = 1./num_actions;
				}
			}
		}
	}

	std::vector<Real> NonStationaryMDP::generateDirichletDistribution(int size) {
		std::vector<Real> distribution(size, 0.0);
		Real sum = 0.0;
		for (int i = 0; i < size - 1; ++i) {
			distribution[i] = generateUniformDistribution(0, 1 - sum);
			sum += distribution[i];
		}
		distribution[size - 1] = 1 - sum;
		return distribution;
	}

	Real NonStationaryMDP::generateUniformDistribution(Real min, Real max) {
		std::uniform_real_distribution<Real> distribution(min, max);
		static std::default_random_engine generator;
		return distribution(generator);
	}

	Real NonStationaryMDP::generateNormalDistribution(Real mean, Real stddev) {
		std::normal_distribution<Real> distribution(mean, stddev);
		static std::default_random_engine generator;
		return distribution(generator);
	}
}