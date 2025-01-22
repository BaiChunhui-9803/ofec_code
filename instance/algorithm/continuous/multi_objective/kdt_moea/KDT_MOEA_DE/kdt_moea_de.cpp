#include "KDT_MOEA_DE.h"
#include "../../../../../record/multi_objective/rcr_vec_real_moea.h"
#include "../../../../../../utility/functional.h"
#include "../../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../../utility/metricsMOP/IGD.h"
#include <fstream>
#include <iomanip>
#include <chrono>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif


namespace ofec {
	//KDT_MOEA_DE::KDT_MOEA_DE(param_map& v) :algorithm(v.at("algName")), m_pop(v.at("popSize")), m_pop_MOEAD(v.at("popSize")), m_offspring(2 * v.at("popSize")), m_Q(v.at("popSize")), m_kdt(v) {}
	void KDT_MOEA_DE::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_pop.reset();
	}

	void KDT_MOEA_DE::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			pop_objs.push_back(m_pop->at(i).objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);
		//Real IGD = CAST_CONOP(m_problem.get())->optimaBase()->invertGenDist(*m_pop);
		entry.push_back(temp_IGD);
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

	void KDT_MOEA_DE::initPop() {
		auto& v = *m_param;
		auto size_var = v.get<int>("number of variables");
		auto size_obj = v.get<int>("number of objectives");
		m_pop.reset(new KDT_MOEA_DE_pop(m_pop_size, m_problem.get(),v));
		m_pop->initialize_(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
	}

	void KDT_MOEA_DE::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
			/* store the population every generation for prediction*/
			//record_obj(this->iteration());
			/* store the population every generation for prediction*/
/*#ifdef OFEC_DEMO
			Real IGD = CONTINUOUS_CAST->get_optima().IGD_to_PF(m_pop);
			vector<vector<solution<>*>> pops(2);
			for (size_t i = 0; i < m_offspring.size(); ++i)
				pops[0].emplace_back(&m_offspring[i]);
			for (size_t i = 0; i < this->m_individuals.size(); ++i)
				pops[1].emplace_back(&m_pop[i]);
			dynamic_cast<ofec_demo::buffer_KDT<DE::MOEA_DE_pop<>, DE::Solution>*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, m_kdt, IGD);
#endif*/
#ifdef OFEC_DEMO
			updateBuffer();
            //vector<vector<solution<>*>> pops(3);
			//Real IGD;
			//vector<vector<solution<>*>> pops;
			//if (this->iteration() == 27) {
			//	IGD = CONTINUOUS_CAST->get_optima().IGD_to_PF(m_pop_MOEAD);
			//	pops.resize(2);
			//	//for (size_t i = 0; i < m_pop_MOEAD.size(); ++i)
			//	//	pops[0].emplace_back(&m_pop_MOEAD[i].solut());
			//	for (size_t i = 0; i < m_Q.size(); ++i)
			//		pops[0].emplace_back(&m_Q[i].solut());
			//	for (size_t i = 0; i < m_pop_MOEAD.size(); ++i)
			//		pops[1].emplace_back(&m_pop_MOEAD[i].solut());
			//}
			//else {
			//	IGD = CONTINUOUS_CAST->get_optima().IGD_to_PF(m_pop);
			//	pops.resize(4);
			//	int predict = global::ms_arg.at("predict");
			//	if (predict == 1) {
			//		m_prediction.resize(m_prediction_kalman_filter.size());
			//		for (size_t i = 0; i < m_prediction_kalman_filter.size(); ++i) {
			//			m_prediction[i] = m_Q[0];
			//			m_prediction[i].objective() = m_prediction_kalman_filter[i];
			//			pops[3].emplace_back(&m_prediction[i].solut());
			//		}
			//		//size_predict = m_prediction_kalman_filter.size();
			//	}
			//	else if (predict == 2) {
			//		m_prediction.resize(m_prediction_least_square.size());
			//		for (size_t i = 0; i < m_prediction_least_square.size(); ++i) {
			//			m_prediction[i] = m_Q[0];
			//			m_prediction[i].objective() = m_prediction_least_square[i];
			//			pops[3].emplace_back(&m_prediction[i].solut());
			//		}
			//		//size_predict = m_prediction_least_square.size();
			//	}
			//	if (m_prediction.size() == 0) {
			//		m_prediction.resize(m_ref_point.size());
			//		for (size_t i = 0; i < m_ref_point.size(); ++i) {
			//			m_prediction[i] = m_Q[0];
			//			m_prediction[i].objective() = m_ref_point[i];
			//			pops[3].emplace_back(&m_prediction[i].solut());
			//		}
			//	}

			//	if (m_count_F1 <= this->m_individuals.size()) {
			//		for (size_t i = 0; i < m_Q.size(); ++i)	//m_A
			//			pops[0].emplace_back(&m_Q[i].solut());
			//		for (size_t i = 0; i < m_count_F1; ++i)
			//			pops[1].emplace_back(&m_pop[i].solut());
			//		/*for (size_t i = 0; i < m_pop_second_select.size(); ++i)
			//			pops[2].emplace_back(&m_pop_second_select[i].solut());*/
			//		for (size_t i = m_count_F1; i < this->m_individuals.size(); ++i)
			//			pops[2].emplace_back(&m_pop[i].solut());

			//	}
			//	else {
			//		for (size_t i = 0; i < m_F1.size(); ++i)
			//			pops[0].emplace_back(&m_F1[i].solut());
			//		/*for (size_t i = 0; i < this->m_individuals.size(); ++i)
			//			pops[1].emplace_back(&m_pop[i].solut());*/
			//		for (size_t i = 0; i < m_pop_first_select.size(); ++i)
			//			pops[1].emplace_back(&m_pop_first_select[i].solut());
			//		for (size_t i = m_pop_first_select.size(); i < this->m_individuals.size(); ++i)
			//			pops[2].emplace_back(&m_pop[i].solut());
			//		/*pops.resize(1);
			//		for (size_t i = 0; i < this->m_individuals.size(); ++i)
			//			pops[0].emplace_back(&m_pop[i].solut());*/
			//	}
			//}

			//dynamic_cast<ofec_demo::buffer_KDT<DE::MOEA_DE_pop<>, DE::Solution>*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, m_kdt, IGD, m_prediction_size, m_prediction.size(), this->iteration());


#endif

			// update m_A
			m_pop->setAllNondomiSolu();	// m_A store nondominated Solutions
			// uodate t
			m_pop->m_iteration++;
		}

		//std::cout << "1: " << m_pop.get_diff1() << std::endl;
		//std::cout << "2: " << m_pop.get_diff2() << std::endl;
	}

#ifdef OFEC_DEMO
	void KDT_MOEA_DE::updateBuffer() {
		/*if(!m_initialized) return;
		Real IGD = CONTINUOUS_CAST->get_optima().IGD_to_PF(m_pop);;
		vector<vector<solution<>*>> pops(3);
		for (size_t i = 0; i < this->m_individuals.size(); ++i)
			pops[0].emplace_back(&m_offspring[i]);
		for (size_t i = this->m_individuals.size(); i < 2*this->m_individuals.size(); ++i)
			pops[1].emplace_back(&m_offspring[i]);
		for (size_t i = 0; i < this->m_individuals.size(); ++i)
			pops[2].emplace_back(&m_pop[i]);
		dynamic_cast<ofec_demo::buffer_KDT<DE::MOEA_DE_pop<>, DE::Solution>*>(ofec_demo::msp_buffer.get())->updateBuffer_(&pops, m_kdt, IGD);*/
	}
#endif

	//void KDT_MOEA_DE::initialize_() {
	//	m_pop->initialize_();
	//	m_pop.evaluate();
	//	m_pop.set_parameter(1.0, 0.5); //set parmeter: CR=1.0, F=0.5
	//	m_kdt->update_min_max(m_pop);
	//	m_kdt->set_kdtree(m_kdt->min_max());
	//	m_Kalman_Filter.resize(this->m_individuals.size());
	//	m_prediction_kalman_filter.reserve(m_offspring.size());
	//	m_prediction_least_square.reserve(m_offspring.size());
	//	m_initialized = true;
	//}

	KDT_MOEA_DE_pop::KDT_MOEA_DE_pop(size_t size_pop, Problem *pro,const ParameterMap& v) : \
		PopMODE<IndDE>(size_pop, pro),\
		m_kdt(std::make_unique<KDT_MOEA<PopMODE<>, IndDE>>(v)), \
		m_pop_MOEAD(PopMODE<IndDE>(size_pop, pro)){ }

	void KDT_MOEA_DE_pop::initialize_(Problem *pro, Random *rnd) {
		Population<IndDE>::initialize(pro, rnd);
		//MOEAD<IndDE>::initialize(this->m_individuals, pro);
	}

	/***********************MOEA/D-start*************************/
	int KDT_MOEA_DE_pop::evolve_MOEAD(PopMODE<>& pop_MOEAD,Problem *pro, Random *rnd) {
		initialize_MOEAD(pop_MOEAD.m_individuals,pro);
		int tag = EvaluationTag::kNormalEval;
		tag=evolve_Mo_MOEAD(pop_MOEAD,pro,rnd);
		return tag;
	}

	int KDT_MOEA_DE_pop::evolve_Mo_MOEAD(PopMODE<>& pop_MOEAD,Problem *pro, Random *rnd) {
		int tag = EvaluationTag::kNormalEval;
		/*std::vector<int> perm(m_pop.m_individuals.size());
		for (int i(0); i < perm.size(); ++i) {
			perm[i] = i;
		}
		global::ms_global->m_uniform[caller::Algorithm]->shuffle(perm.begin(), perm.end());*/
		for (int i = 0; i < pop_MOEAD.m_individuals.size(); i++) {
			//int n = perm[i];
			//// or int n = i;
			//int type;
			//double rnd = global::ms_global->m_uniform[caller::Algorithm]->next();
			//// mating selection based on probability
			//if (rnd < m_Realb)
			//	type = 1;   // neighborhood
			//else
			//	type = 2;   // whole population
			//// select the indexes of mating parents
			//std::vector<int> p;
			//matingselection(p, n, 2, type, m_pop.m_individuals.size());  // neighborhood selection
			//// produce a child Solution
			//DE::Solution child;
			//std::vector<size_t> index(3);
			//index[0] = n; index[1] = p[0]; index[2] = p[1];
			//this->cross_mutate(index, child);
			//tag = child.evaluate();
			//if (tag != kNormalEval) break;
			//// update the reference points and other TypeIndivs in the neighborhood or the whole population
			//update_reference(child);
			update_problem_MOEAD(pop_MOEAD.m_individuals, m_Q[i], i, m_type[i],pro,rnd);
		}
		return tag;
	}

	void KDT_MOEA_DE_pop::initialize_MOEAD(std::vector<std::unique_ptr<ofec::IndDE>>& parent,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		mv_sol_arr.resize(M);
		mv_ideal_point.resize(M, 1.0e+30);
		init_uniformweight_MOEAD(parent.size(),pro);
		init_neighbourhood_MOEAD();
		for (int i = 0; i < parent.size(); i++)
			update_reference_MOEAD(*(parent[i]),pro);
	}

	void KDT_MOEA_DE_pop::init_uniformweight_MOEAD(int parent_size,Problem *pro) {
		if (CAST_CONOP(pro)->numberObjectives() == 2) {
			mv_namda.resize(parent_size);
			for (int n = 0; n < parent_size; n++) {
				Real a = 1.0 * n / (parent_size - 1);
				mv_namda[n].pushBack(a);
				mv_namda[n].pushBack(1 - a);
			}
		}
		else {
			int n = 0;
			for (int i = 0; i <= m_unit; i++) {
				for (int j = 0; j <= m_unit; j++) {
					if (i + j <= m_unit) {
						std::vector<int> arr;
						arr.push_back(i);
						arr.push_back(j);
						arr.push_back(m_unit - i - j);
						mv_namda.push_back(std::vector<Real>(0));
						for (int k = 0; k < arr.size(); k++)
							mv_namda[n].pushBack(1.0 * arr[k] / m_unit);
						n++;
						arr.clear();
					}
				}
			}
		}
	}

	void KDT_MOEA_DE_pop::init_neighbourhood_MOEAD() {
		int pops = mv_namda.size();
		mvv_neigh.resize(pops);
		std::vector<Real> dis(pops);
		std::vector<int> index(pops);
		for (int i = 0; i < pops; i++) {
			// calculate the distances based on weight vectors
			for (int j = 0; j < pops; j++) {
				dis[j] = mv_namda[i].distance(mv_namda[j]);
				index[j] = j;
			}
			//find 'niche' nearest neighboring subproblems			
			min_max_sort(dis, index);
			for (int k = 0; k < m_niche; k++)
				mvv_neigh[i].push_back(index[k]);
		}
		dis.clear();
		index.clear();
	}

	void KDT_MOEA_DE_pop::min_max_sort(std::vector<Real>& v, std::vector<int>& idx) {
		for (int i = 0; i < v.size(); ++i) {
			int min = i;
			for (int j = i + 1; j < v.size(); ++j) {
				if (v[j] < v[min])
					min = j;
			}
			auto temp = v[i];
			v[i] = v[min];
			v[min] = temp;
			auto t = idx[i];
			idx[i] = idx[min];
			idx[min] = t;
		}
	}

	void KDT_MOEA_DE_pop::update_problem_MOEAD(std::vector<std::unique_ptr<IndDE>>& parent, IndDE& sol, int id, int type,Problem *pro, Random *rnd) {
		// sol: child Solution
		// id:  the id of current subproblem
		// type:update Solutions in - neighborhood (1) or whole population (otherwise)
		int size, time = 0;
		if (type == 1)
			size = m_niche;
		else
			size = parent.size();
		std::vector<int> perm(size);
		for (int i(0); i < perm.size(); ++i) {
			perm[i] = i;
		}
		
		rnd->uniform.shuffle(perm.begin(), perm.end());

		for (int i = 0; i < size; i++) {
			int k;
			if (type == 1)
				k = mvv_neigh[id][perm[i]];
			else
				k = perm[i];
			// calculate the values of objective function regarding the current subproblem
			Real f1, f2;
			f1 = fitnessfunction_MOEAD(parent[k]->objective(), k,pro);
			f2 = fitnessfunction_MOEAD(sol.objective(), k,pro);
			if (f2 < f1) {
				*parent[k] = sol;
				time++;
			}
			// the maximal number of Solution updated is not allowed to exceed 'limit'
			if (time >= m_limit)
				return;
		}
		perm.clear();
	}

	Real KDT_MOEA_DE_pop::fitnessfunction_MOEAD(std::vector<Real>& obj, int k,Problem *pro) {
		// Chebycheff Scalarizing Function
		Real fitness = 0;
		int numObj = CAST_CONOP(pro)->numberObjectives();
		//int numObj = obj.size();
		if (m_decom_function == _TCHE1) {
			Real max_fun = -1.0e+30;
			for (int n = 0; n < numObj; n++) {
				//Real diff = fabs(y_obj[n] - idealpoint[n] + scale[n]);
				//Real diff = fabs(y_obj[n] - idealpoint[n] + 0.05);
				Real diff = fabs(obj[n] - mv_ideal_point[n]);
				//Real diff = fabs(y_obj[n] - 0);
				Real feval;
				if (mv_namda[k][n] == 0)
					feval = 0.0001 * diff;
				else
					feval = diff * mv_namda[k][n];
				if (feval > max_fun) max_fun = feval;

			}
			fitness = max_fun;
		}

		if (m_decom_function == _TCHE2) {
			// reference point in the CHIM
			std::vector<int> scale(numObj);
			//throw myException("Please initialize the scale @MOEAD<Poppulation,Solution>::fitnessfuction");
			Real max_fun = -1.0e+30;
			for (int n = 0; n < numObj; n++) {
				Real diff = (obj[n] - mv_ideal_point[n]) / scale[n];  //note: the scale is not initialized, there has no knowledge
				Real feval;
				if (mv_namda[k][n] == 0)
					feval = 0.0001 * diff;
				else
					feval = diff * mv_namda[k][n];
				if (feval > max_fun) max_fun = feval;

			}
			fitness = max_fun;
		}

		// CHIM + Tchebycheff
		// CHIM is not available in 3 objectives
		if (m_decom_function == _NBI1) {
			// quasi normal direction
			Vector norm;
			for (int i = 0; i < numObj; i++) {
				norm.pushBack(0.0);
				for (int j = 0; j < numObj; j++) {
					norm[i] += -mv_sol_arr[j].objective()[i];
				}
			}

			// normalization
			norm.normalize();

			// reference point in the CHIM
			std::vector <Real> base;
			for (int i = 0; i < numObj; i++) {
				Real tp2 = 0;
				for (int j = 0; j < numObj; j++)
					tp2 += mv_sol_arr[j].objective()[i] * mv_namda[k][j];
				base.push_back(tp2);
			}

			// Tchebycheff function
			Real max_fun = -1.0e+30;
			for (int n = 0; n < numObj; n++) {
				Real diff = obj[n] - base[n];
				Real feval = -diff * norm[n];
				if (feval > max_fun) max_fun = feval;

			}
			fitness = max_fun;
		}

		//* Boundary intersection approach
		//* reference point is chosen as the ideal point
		//* the direction is independent of CHIM
		if (m_decom_function == _NBI2) {

			mv_namda[k].normalize();

			// penalty method 
			// temporary vectors NBI method
			Vector RealA(numObj);
			Vector RealB(numObj);

			// difference beween current point and reference point
			for (int n = 0; n < numObj; n++)
				RealA[n] = (obj[n] - mv_ideal_point[n]);

			// distance along the search direction norm
			Real d1 = fabs(RealA * mv_namda[k]);

			// distance to the search direction norm
			for (int n = 0; n < numObj; n++)
				RealB[n] = (obj[n] - (mv_ideal_point[n] + d1 * mv_namda[k][n]));
			Real d2 = RealB.norm();

			fitness = (d1 + 5 * d2);

			//t2 = clock();
			//total_sec+=(t2 - t1);
		}

		// NBI method
		if (m_decom_function == _NBI3) {

			// quasi normal direction
			Vector norm;
			for (int i = 0; i < numObj; i++) {
				norm.pushBack(0.0);
				for (int j = 0; j < numObj; j++) {
					norm[i] += -mv_sol_arr[j].objective()[i];
				}
			}

			// normalization
			norm.normalize();

			// reference point in the CHIM
			std::vector<Real> base;
			for (int i = 0; i < numObj; i++) {
				Real tp2 = 0;
				for (int j = 0; j < numObj; j++)
					tp2 += mv_sol_arr[j].objective()[i] * mv_namda[k][j];
				base.push_back(tp2);
			}

			// penalty method 
			// temporary vectors NBI method
			Vector RealA;
			Vector RealB;

			// difference beween current point and reference point
			for (int n = 0; n < numObj; n++)
				RealA.pushBack(obj[n] - base[n]);

			// distance along the search direction norm
			Real d1 = RealA * norm;

			// distance to the search direction norm
			for (int n = 0; n < numObj; n++)
				RealB.pushBack(obj[n] - (base[n] + d1 * norm[n]));
			Real d2 = RealB.norm();

			fitness = -d1 + 2 * d2;
		}
		return fitness;
	}

	void KDT_MOEA_DE_pop::update_reference_MOEAD(IndDE& sol,Problem *pro) {
		//sol: child Solution
		int numObj = CAST_CONOP(pro)->numberObjectives();
		for (int n = 0; n < numObj; n++) {
			if (sol.objective()[n] < mv_ideal_point[n]) {
				mv_ideal_point[n] = sol.objective()[n];
				mv_sol_arr[n] = sol;
			}
		}
	}
	/***********************MOEA/D-end***************************/


	void KDT_MOEA_DE_pop::recordObj(int iter,Algorithm *alg,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		std::ofstream out;
		out.open("./result_rnn/KDT_obj.txt", std::ios::app);
		if (out && M == 2) {
			out << "the " << m_iteration << " generation" << std::endl;
			out << "f1\t" << "f2" << std::endl;
			for (int i = 0; i < m_archive[iter - 1].size(); ++i) {
				for (int j = 0; j < M; ++j) {
					if (j == 0) {
						out << m_archive[iter - 1][i].objective()[j];
					}
					else {
						out << "\t" << m_archive[iter - 1][i].objective()[j];
					}

				}
				out << std::endl;
			}
			if (m_archive.size() >= 3) {
				out << "the " << this->iteration() << " generation" << std::endl;
				out << "predict result" << std::endl;
				for (int i = 0; i < m_prediction_least_square.size(); ++i) {
					for (int j = 0; j < M; ++j) {
						if (j == 0)
							out << m_prediction_least_square[i][j];
						else
							out << "\t" << m_prediction_least_square[i][j];
					}
					out << std::endl;
				}
			}
		}
		else if (out && M > 2) {
			out << "the " << this->iteration() << " generation" << std::endl;
			out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::setw(10) << "f3" << std::endl;
			for (int i = 0; i < m_archive[iter - 1].size(); ++i) {
				for (int j = 0; j < M; ++j) {
					if (j == 0) {
						out << m_archive[iter - 1][i].objective()[j];
					}
					else {
						out << "\t" << m_archive[iter - 1][i].objective()[j];
					}
					//out << std::setw(10) << m_pop[i].objective()[j];
				}
				out << std::endl;
			}
			if (m_archive.size() >= 3) {
				out << "the " << this->iteration() << " generation" << std::endl;
				out << "predict result" << std::endl;
				for (int i = 0; i < m_prediction_least_square.size(); ++i) {
					for (int j = 0; j < M; ++j) {
						if (j == 0)
							out << m_prediction_least_square[i][j];
						else
							out << "\t" << m_prediction_least_square[i][j];
					}
					out << std::endl;
				}
			}
		}
		out.close();

		//if (iter == 0) {
		//	std::ofstream out;
		//	out.open("./result_rnn/KDT_obj.txt", std::ios::app);
		//	if (out && m == 2) {
		//		out << "the " << this->iteration() << " generation" << std::endl;
		//		//out << "f1\t" << "f2" << std::endl;
		//		for (int i = 0; i < this->m_individuals.size(); ++i) {
		//			for (int j = 0; j < m; ++j) {
		//				out << "\t" << m_pop[i].objective()[j];
		//			}
		//			out << std::endl;
		//		}
		//	}
		//	else if (out && m == 3) {
		//		out << "the " << this->iteration() << " generation" << std::endl;
		//		//out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::setw(10) << "f3" << std::endl;
		//		for (int i = 0; i < this->m_individuals.size(); ++i) {
		//			for (int j = 0; j < m; ++j) {
		//				out << std::setw(10) << m_pop[i].objective()[j];
		//			}
		//			out << std::endl;
		//		}
		//	}
		//	out.close();
		//}
		//else {
		//	std::ofstream out;
		//	out.open("./result_rnn/KDT_obj.txt", std::ios::app);
		//	if (out && m == 2) {
		//		out << std::endl << "the " << this->iteration() << " generation" << std::endl;
		//		//out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::endl;
		//		for (int i = 0; i < this->m_individuals.size(); ++i) {
		//			for (int j = 0; j < m; ++j) {
		//				out << std::setw(10) << m_pop[i].objective()[j];
		//			}
		//			out << std::endl;
		//		}
		//	}
		//	else if (out && m == 3) {
		//		out << std::endl << "the " << this->iteration() << " generation" << std::endl;
		//		//out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::setw(10) << "f3" << std::endl;
		//		for (int i = 0; i < this->m_individuals.size(); ++i) {
		//			for (int j = 0; j < m; ++j) {
		//				out << std::setw(10) << m_pop[i].objective()[j];
		//			}
		//			out << std::endl;
		//		}
		//	}
		//	out.close();
		//}
	}

	int KDT_MOEA_DE_pop::evolve_mo(Problem *pro, Algorithm *alg, Random *rnd) {
		int pop_size = this->m_individuals.size();
		/*consider neighborhood*/
		std::vector<std::vector<size_t>> neighborhood(pop_size);
		for (int i = 0; i < pop_size; ++i) {
			//neighborhood[i].push_back(i);
			int k1 = m_kdt->get_kdt().getRegionIdx(this->m_individuals[i]->objective());
			for (int j = i + 1; j < pop_size; ++j) {
				int k2 = m_kdt->get_kdt().getRegionIdx(this->m_individuals[j]->objective());
				if (m_kdt->get_kdt().checkAdjacency(k1, k2)) {
					neighborhood[i].push_back(j);
					neighborhood[j].push_back(i);
				}
			}
		}
		//test
		std::cout << "neighborhood:" << std::endl;
		for (int i = 0; i < pop_size; ++i) {
			std::cout << i << ": ";
			for (int j = 0; j < neighborhood[i].size(); ++j)
				std::cout << neighborhood[i][j] << " ";
			std::cout << std::endl;
		}

		int tag = EvaluationTag::kNormalEval;
		double delta = 0.9;
		int m = 0;
		for (int i = 0; i < pop_size; ++i) {
			double rand = rnd->uniform.next();	//get a ramdom value in (0,1)				
			std::vector<size_t> result(2);
			if (rand < delta && neighborhood[i].size() >= 2)
				selectInCandidates(2, neighborhood[i], result,rnd);
			else
				select(i, 2, result,rnd);
			result.push_back(i);
			/*adopt DE operator*/
			m_offspring[m] = *(this->m_individuals[i]);
			m_offspring[m].setType(-1);
			m++;
			crossMutate(result, m_offspring[m],pro,rnd);
			tag = m_offspring[m].evaluate(pro,alg);
			m_offspring[m].setType(-1);

			if (tag != EvaluationTag::kNormalEval) break;
			m++;
		}
		return tag;
	}

	int KDT_MOEA_DE_pop::evolve_mo_new(Problem *pro, Algorithm *alg, Random *rnd) {
		m_type.resize(this->m_individuals.size(), 1); //for comparing MOEA/D

		int pop_size = this->m_individuals.size();
		/*consider neighborhood*/
		std::vector<std::vector<size_t>> neighborhood(pop_size);
		for (int i = 0; i < pop_size; ++i) {
			//neighborhood[i].push_back(i);
			int k1 = m_kdt->get_kdt().getRegionIdx(this->m_individuals[i]->objective());
			for (int j = i + 1; j < pop_size; ++j) {
				int k2 = m_kdt->get_kdt().getRegionIdx(this->m_individuals[j]->objective());
				if (m_kdt->get_kdt().checkAdjacency(k1, k2)) {
					neighborhood[i].push_back(j);
					neighborhood[j].push_back(i);
				}
			}
		}
		//test
		/*std::cout << "neighborhood:" << std::endl;
		for (int i = 0; i < pop_size; ++i) {
			std::cout << i << ": ";
			for (int j = 0; j < neighborhood[i].size(); ++j)
				std::cout << neighborhood[i][j] << " ";
			std::cout << std::endl;
		}*/

		int tag = kNormalEval;
		double delta = 0.9;
		int m = 0;
		for (int i = 0; i < pop_size; ++i) {
			double rand = rnd->uniform.next();	//get a ramdom value in (0,1)				
			std::vector<size_t> result(2);
			if (rand < delta && neighborhood[i].size() >= 2)
				selectInCandidates(2, neighborhood[i], result,rnd);
			else {
				select(i, 2, result,rnd);
				m_type[i] = 2;
			}
			result.push_back(i);
			/*adopt DE operator*/
			this->crossMutate(result, m_Q[m],pro,rnd);	//m_Q only store Q
			tag = m_Q[m].evaluate(pro,alg);

			if (tag != kNormalEval) break;
			m++;
		}
		return tag;
	}

	int KDT_MOEA_DE_pop::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		if (this->m_individuals.size() < 5) {
			throw MyExcept("the population size cannot be smaller than 5@DE::population<Solution>::evolve()");
		}
		int tag = kNormalEval;

		/****test LeastSquare****/
		/*LeastSquare<DE::Solution> test;
		std::vector<Real> result;
		std::vector<Real> vt(3);
		vt[0] = 2;
		vt[1] = 1;
		vt[2] = 0.5;
		Real r = test.least_square_estimation_e(vt);
		std::vector<Real> result_kalman;
		std::vector<Real> result_least;
		std::vector<DE::Solution> v(4);
		v[0].objective() = { 3,3,3 };
		v[1].objective() = { 2.5,1,1.5 };
		v[2].objective() = { 1,0,0.5 };
		v[3].objective() = { 0.5,-0.2,0 };
		LeastSquare<DE::Solution> least;
		result_least = least.predict(v);
		Kalman_Filter<DE::Solution> kalman;
		result_kalman = kalman.predict(v);*/
		/*****test judge_point*****/
		/*std::vector<std::vector<Real>> test(3, std::vector<Real>(3));
		test[0] = { 1,2,3};
		test[1] = { 0.1,0.2,0.3 };
		test[2] = { 0.15,0.2, 0.5};
		std::vector<Real> point = { 0.1,0.2,0.3};
		bool result = judge_point(point, test);*/

		std::cout << std::endl;
		std::cout << "the " << this->iteration() << " generation" << std::endl;

		//tag = evolve_mo();
		tag = evolve_mo_new(pro,alg,rnd);//for compare MOEA/D

		if (this->iteration() == 0) {
			for (int i = 0; i < this->m_individuals.size(); ++i) {
				m_A.push_back(*(this->m_individuals[i]));
			}
			for (int i = 0; i < m_Q.size(); ++i) {
				m_A.push_back(m_Q[i]);
			}
		}
		else {
			for (int i = 0; i < this->m_individuals.size(); ++i) {
				if (this->m_individuals[i]->fitness() != 0) {
					m_A.push_back(*(this->m_individuals[i]));
				}
			}
			for (int i = 0; i < m_Q.size(); ++i) {
				m_A.push_back(m_Q[i]);
			}
		}

		for (int i = 0; i < m_A.size(); ++i) {
			m_A[i].setType(-1);
		}

		/*****sort, get m_F1 and m_F2_Fn*****/
		m_count_F1 = find_non_dominance(m_A,pro);
		m_F1.resize(m_count_F1);
		m_F2_Fn.resize(m_A.size() - m_count_F1);
		int index1, index2;
		index1 = index2 = 0;
		for (int i = 0; i < m_A.size(); ++i) {
			if (m_A[i].fitness() == 0) {
				m_F1[index1++] = m_A[i];
			}
			else {
				m_F2_Fn[index2++] = m_A[i];
			}
		}
		if ((index1 != m_count_F1) || ((index1 + index2) != m_A.size()))
			std::cout << "number of m_F1 and m_F2_Fn is error!" << std::endl;


		/*****sort, get m_F1 and m_F2_Fn*****/
		/*m_count_F1 = find_non_dominance(m_offspring);
		m_F1.resize(m_count_F1);
		m_F2_Fn.resize(m_offspring.size() - m_count_F1);
		int index1, index2;
		index1 = index2 = 0;
		for (int i = 0; i < m_offspring.size(); ++i) {
			if (m_offspring[i].fitness() == 0) {
				m_F1[index1++] = m_offspring[i];
			}
			else {
				m_F2_Fn[index2++] = m_offspring[i];
			}
		}
		if ((index1 != m_count_F1) || ((index1 + index2) != m_offspring.size()))
			std::cout << "number of m_F1 and m_F2_Fn is error!" << std::endl;*/

			//test, display objective of m_F1
			/*int M = global::ms_global->m_problem->objective_size();
			std::cout << "m_F1: " << std::endl;
			for (int i = 0; i < m_F1.size(); ++i) {
				std::cout << "[" << i << "]:";
				for (int j = 0; j < M; ++j) {
					if (j == 0)
						std::cout << m_F1[i].objective(j);
					else
						std::cout << "\t" << m_F1[i].objective(j);
				}
				std::cout << std::endl;
			}*/

			/*****update m_all_nondominance*****/
			/*if (m_all_nondominance.size() == 0) {
				m_all_nondominance.reserve(m_F1.size());
				for (int i = 0; i < m_F1.size(); ++i) {
					m_all_nondominance.push_back(m_F1[i].objective());
				}
			}
			else {
				std::vector<int> add;
				std::vector<int> del;
				for (int i = 0; i < m_F1.size(); ++i) {
					bool flag = true;
					for (int j = 0; j < m_all_nondominance.size(); ++j) {
						dominationship rship = objective_compare(m_all_nondominance[j], m_F1[i].objective(), global::ms_global->m_problem->opt_mode());
						if (dominationship::Equal == rship) {
							flag = false;
							break;
						}
						else if (dominationship::Dominating == rship) {
							flag = false;
						}
						else if (dominationship::Dominated == rship) {
							if (std::find(del.begin(), del.end(), j) == del.end()) {
								del.push_back(j);
							}
						}
					}
					if (flag) {
						add.push_back(i);
					}
				}
				if (del.size() > 0) {
					std::sort(del.begin(), del.end());
					for (int i = del.size() - 1; i >= 0; --i) {
						m_all_nondominance.erase(m_all_nondominance.begin() + del[i]);
					}
				}
				if (add.size() > 0) {
					m_all_nondominance.reserve(m_all_nondominance.size() + add.size());
					for (int i = 0; i < add.size(); ++i) {
						m_all_nondominance.push_back(m_F1[add[i]].objective());
					}
				}
			}*/

			//eval_selection();
			//eval_selection_();  //select an Solution every cycle
			//eval_selection__(); //select an Solution every cycle, and consider F1 of nondominated-sorting


		m_prediction_kalman_filter.resize(0);
		m_prediction_least_square.resize(0);
		m_ref_point.resize(0);

		if (this->iteration() == 27) {
			for (int i = 0; i < this->m_individuals.size(); ++i) {
				m_pop_MOEAD[i] = *(this->m_individuals[i]);
				m_Q.push_back(*(this->m_individuals[i]));//for test
			}
			evolve_MOEAD(m_pop_MOEAD,pro,rnd);
		}

		eval_selection_new(pro); //select an Solution every cycle, and consider F1 of nondominated-sorting, and sample reference point, and not use cdi indicator

		//eval_selection_new(); //select an Solution every cycle, and consider F1 of nondominated-sorting, and sample reference point, and not use cdi indicator


		/*****using NSGAII to select*****/
		/*nondominated_sorting(m_offspring);
		eval_selection_NSGAII(m_pop, m_offspring);*/

		/*****handle m_archive*****/
		find_non_dominance(pro); //find non-dominated Solutions from m_pop, and store them into m_arhive
		if (this->iteration() != 0)
			update_index_archive(); //update index of inviduals according to their distance

		/*****update m_all_nondominance*****/
		//m_A = m_F1;	// m_A store nondominated Solutions

		//this->iteration()++;
		return tag;
	}

	//select an Solution every cycle, and consider F1 of nondominated-sorting, and sample reference point, and not use cdi indicator
	void KDT_MOEA_DE_pop::eval_selection_new(Problem *pro) {
		m_kdt->set_F1(m_count_F1);
		//clear m_pop_first_select and m_pop_second_select
		m_pop_first_select.resize(0);
		m_pop_second_select.resize(0);
		// step3 : select Solutions according to number of m_count_F1
		int pops = 0;
		if (m_count_F1 <= this->m_individuals.size()) {
			m_kdt->update_min_max(m_A,pro); // partition by m_A
			m_kdt->set_kdtree(m_kdt->min_max(),pro);
			for (int i = 0; i < m_F1.size(); ++i) {
				*(this->m_individuals[pops]) = m_F1[i];
				m_F1[i].setType(1);
				++pops;
			}
			if (pops < this->m_individuals.size()) {
				select_F2_Fn_new(m_F2_Fn, pops,pro);
			}
		}
		else {
			m_kdt->update_min_max(m_F1,pro);	//partition by F1
			m_kdt->set_kdtree(m_kdt->min_max(),pro);
			select_F1(m_F1, pops,pro);
		}
	}

	void KDT_MOEA_DE_pop::select_F2_Fn(std::vector<IndDE>& F2_Fn, int& pops,Problem *pro) {
		m_kdt->updateRegionLabel(F2_Fn,pro);
		int size_region = m_kdt->size_region();
		/*calculate Solutions for each subregion*/
		std::vector<std::vector<int>> region_to_sort_ind(size_region);
		int number = update_region_to_sort_ind(F2_Fn, region_to_sort_ind);
		if (number != F2_Fn.size())
			std::cout << "The matching subregion of Solution is error!" << std::endl;

		/*prediction*/
		if (m_archive.size() >= 2) {
			predict();
			add_pop_to_pred(pro);		//add the non-domianted Solutions of pre-generation population into prediction (its left point in box)
			m_kdt->cal_dds(F2_Fn, region_to_sort_ind,pro);
			add_subregion_to_pred(m_kdt->region(), 1,pro);//add the non-dominated(type=0) or all(type=1, sub_dds>0) subregions into prediction (its left lower point in box)
			//add_shrink_F1_to_pred();			
			//observe m_prediction_kalman_filter
			/*std::cout << "prediction_kalman_filter:" << std::endl;
			for (int i = 0; i < m_prediction_kalman_filter.size(); ++i) {
				for (int j = 0; j < m_prediction_kalman_filter[0].size(); ++j) {
					if (j == 0)
						std::cout << m_prediction_kalman_filter[i][j];
					else
						std::cout << "\t" << m_prediction_kalman_filter[i][j];
				}
				std::cout << std::endl;
			}*/
			//record number of prediction in pre and post of handling
			auto para_map = pro->getParam();
			int pred_method = para_map->get<int>("predict");
			if (pred_method == 1)
				m_prediction_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_size = m_prediction_least_square.size();
			handle_prediction(m_offspring,pro);//erase these dominated-prediction-point by m_A
			if (pred_method == 1)
				m_prediction_handled_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_handled_size = m_prediction_least_square.size();
			add_shrink_all_nondominance_to_pred(m_F1);
		}
		else {
			add_pop_to_pred(pro);
			m_kdt->cal_dds(m_offspring, region_to_sort_ind,pro);
			add_subregion_to_pred(m_kdt->region(), 1,pro);
			//add_shrink_F1_to_pred();
			handle_prediction(m_offspring,pro);//erase these dominated-prediction-point
			add_shrink_all_nondominance_to_pred(m_F1);
		}

		/*select offspring*/
		int pops_F1 = pops;	//record the number of the first selection 
		int test = pops;
		std::vector<std::pair<int, std::pair<int, int>>> change_of_pops;	//for test change of pops
		change_of_pops.reserve(F2_Fn.size() / 2);
		change_of_pops.push_back(std::pair<int, std::pair<int, int>>(pops_F1, std::pair<int, int>(-1, -1)));
		m_kdt->m_observe = true;
		while (pops < this->m_individuals.size()) {
			std::vector<int> rank_dds_sa;
			//m_kdt->sort_region(m_pop, rank_dds_sa, m_offspring, region_to_sort_ind, pops, pops_F1);	//non-dominated sort of subspaces 
			auto para_map = CAST_CONOP(pro)->getParam();
			int predict = para_map->get<int>("predict");
			if (m_archive.size() >= 2 && predict == 1)
				m_kdt->sort_region_predict(m_prediction_kalman_filter, *this, rank_dds_sa, F2_Fn, region_to_sort_ind, pops, pops_F1,pro);	//non-dominated sort of subspaces 
			else if (m_archive.size() >= 3 && predict == 2)
				m_kdt->sort_region_predict(m_prediction_least_square, *this, rank_dds_sa, F2_Fn, region_to_sort_ind, pops, pops_F1,pro);
			else
				m_kdt->sort_region_predict(m_ref_point, *this, rank_dds_sa, F2_Fn, region_to_sort_ind, pops, pops_F1,pro);
			//m_kdt->sort_region(m_pop, rank_dds_sa, F2_Fn, region_to_sort_ind, pops, pops_F1);	//non-dominated sort of subspaces

		// step4 : select Solution from subspace of the first rank	
			int num = 0;
			std::vector<std::pair<int, std::pair<int, int>>> sub_inds; //the first store index of Solution, the second store index of subspace
			sub_inds.reserve(F2_Fn.size());
			for (int k = 0; k < m_kdt->size_region(); ++k) {
				if (rank_dds_sa[k] == 0 && pops < this->m_individuals.size() && region_to_sort_ind[k].size() > 0) {
					for (int i = 0; i < region_to_sort_ind[k].size(); ++i) {
						sub_inds.push_back(std::pair<int, std::pair<int, int>>(region_to_sort_ind[k][i], std::pair<int, int>(k, i)));
					}
					num++;
				}
			}
			if (sub_inds.size() == 0)
				std::cout << "sub_inds is error!" << std::endl;
			select_an_ind_max_min(pops, F2_Fn, sub_inds, region_to_sort_ind);	//select an Solution from sub_inds		
			change_of_pops.push_back(std::pair<int, std::pair<int, int>>(pops, std::pair<int, int>(num, m_kdt->size_region())));//for test change of pops
			if (test == pops_F1 && pops_F1 < pops) {	//for test the second selection
				m_pop_second_select.resize(pops - pops_F1);
				for (; test < pops; ++test) {
					m_pop_second_select[test - pops_F1] = *(this->m_individuals[test]);
				}
			}

			/*update m_region_sort of kd-tree*/
			for (int k = 0; k < m_kdt->size_region();) {
				if (region_to_sort_ind[k].size() == 0) {
					region_to_sort_ind.erase(region_to_sort_ind.begin() + k);
					m_kdt->region().erase(m_kdt->region().begin() + k);
				}
				else
					++k;
			}
		}

		/*for test change of pops*/
		std::cout << "selection process: " << change_of_pops[0].first << "  ";
		for (int i = 1; i < change_of_pops.size(); ++i) {
			//output Solutions number every cycle
			std::cout << change_of_pops[i].first - change_of_pops[i - 1].first << "(" << change_of_pops[i].second.first << "/" << change_of_pops[i].second.second << ")  ";
		}
		std::cout << std::endl;
	}

	void KDT_MOEA_DE_pop::select_F2_Fn_new(std::vector<IndDE>& F2_Fn, int& pops,Problem *pro) {
		m_kdt->updateRegionLabel(F2_Fn,pro);	//update m_kdt->m_region_sort
		int size_region = m_kdt->size_region();
		/*calculate selected Solutions for each subregion*/
		std::vector<std::vector<int>> region_to_sort_selected_ind(size_region);
		int number1 = 0;
		/*for (int k = 0; k < size_region; ++k) {
			for (int i = 0; i < pops; ++i) {
				int regionIdx = m_kdt->get_kdt().getRegionIdx(m_pop[i].objective());
				if (m_kdt->region()[k].first == regionIdx) {
					region_to_sort_selected_ind[k].push_back(i);
				}
			}
			number1 += region_to_sort_selected_ind[k].size();
		}*/
		int number1_ = 0;
		for (int i = 0; i < pops; ++i) {
			int regionIdx = m_kdt->get_kdt().getRegionIdx(this->m_individuals[i]->objective());
			bool flag = true;
			for (int k = 0; k < size_region; ++k) {
				if (m_kdt->region()[k].first == regionIdx) {
					region_to_sort_selected_ind[k].push_back(i);
					flag = false;
					break;
				}
			}
			if (flag)
				number1_++;
		}
		for (int k = 0; k < size_region; ++k) {
			number1 += region_to_sort_selected_ind[k].size();
		}

		if ((number1 + number1_) != pops)
			std::cout << "The matching subregion of selected-Solution is error!" << std::endl;
		/*calculate Solutions for each subregion*/
		std::vector<std::vector<int>> region_to_sort_ind(size_region);
		int number2 = update_region_to_sort_ind(F2_Fn, region_to_sort_ind);
		if (number2 != F2_Fn.size())
			std::cout << "The matching subregion of Solution is error!" << std::endl;

		//test
		/*int M = global::ms_global->m_problem->objective_size();
		std::cout << "subspace" << std::endl;
		for (int k = 0; k < size_region; ++k) {
			std::cout << "[" << k << "] \t selected:" << region_to_sort_selected_ind[k].size() << "\t candicate: " << region_to_sort_ind[k].size() << std::endl;
			auto box = m_kdt->get_kdt().get_box(m_kdt->region()[k].first);
			for (int j = 0; j < M; ++j) {
				std::cout << std::setw(10) << box[j].first << std::setw(10) << box[j].second << std::endl;
			}
		}
		std::cout << "selected-Solutions: " << std::endl;
		for (int k = 0; k < size_region; ++k) {
			std::cout << "[" << k << "]: ";
			for (int i = 0; i < region_to_sort_selected_ind[k].size(); ++i)
				std::cout << region_to_sort_selected_ind[k][i] << " ";
			std::cout << std::endl;
		}*/

		/*prediction*/
		if (m_archive.size() >= 2) {
			predict();
			add_pop_to_pred(pro);		//add the non-domianted Solutions of pre-generation population into prediction (its left point in box)
			m_kdt->cal_dds(F2_Fn, region_to_sort_ind,pro);
			add_subregion_to_pred(m_kdt->region(), 1,pro);//add the non-dominated(type=0) or all(type=1, sub_dds>0) subregions into prediction (its left lower point in box)		
			//record number of prediction in pre and post of handling
			auto para_map = pro->getParam();
			int pred_method = para_map->get<int>("predict");
			if (pred_method == 1)
				m_prediction_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_size = m_prediction_least_square.size();
			handle_prediction(m_A,pro);//erase these dominated-prediction-point by m_A
			if (pred_method == 1)
				m_prediction_handled_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_handled_size = m_prediction_least_square.size();
			add_shrink_all_nondominance_to_pred(m_F1);

			if (m_archive.size() < 2)
				m_prediction_handled_size = m_ref_point.size();
			else if (pred_method == 1)
				m_prediction_handled_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_handled_size = m_prediction_least_square.size();
		}
		else {
			add_pop_to_pred(pro);
			m_kdt->cal_dds(F2_Fn, region_to_sort_ind,pro);
			add_subregion_to_pred(m_kdt->region(), 1,pro);
			//add_shrink_F1_to_pred();
			handle_prediction(m_A,pro);//erase these dominated-prediction-point by m_A
			add_shrink_all_nondominance_to_pred(m_F1);
		}

		/*select offspring*/
		int pops_F1 = pops;	//record the number of the first selection 
		int test = pops;
		std::vector<std::pair<int, std::pair<int, int>>> change_of_pops;	//for test change of pops
		change_of_pops.reserve(F2_Fn.size() / 2);
		change_of_pops.push_back(std::pair<int, std::pair<int, int>>(pops_F1, std::pair<int, int>(-1, -1)));
		m_kdt->m_observe = true;
		m_observe = true;
		m_observe_angle_subspace_inner.resize(0);
		while (pops < this->m_individuals.size()) {
			std::vector<int> rank_norm_angle;
			auto para_map = pro->getParam();
			int predict = para_map->get<int>("predict");
			if (m_archive.size() >= 2 && predict == 1)
				m_kdt->sort_region_norm_angle(m_prediction_kalman_filter, *this, rank_norm_angle, F2_Fn, region_to_sort_ind, pops,pro);
			else if (m_archive.size() >= 3 && predict == 2)
				m_kdt->sort_region_norm_angle(m_prediction_least_square, *this, rank_norm_angle, F2_Fn, region_to_sort_ind, pops,pro);
			else
				m_kdt->sort_region_norm_angle(m_ref_point, *this, rank_norm_angle, F2_Fn, region_to_sort_ind, pops,pro);

			// step4 : select Solution from subspace of the first rank	
			int num = 0;
			//std::vector<std::pair<int, std::pair<int, int>>> sub_inds; //the first store index of Solution, the second store index of subspace
			//sub_inds.reserve(F2_Fn.size());
			for (int k = 0; k < m_kdt->size_region(); ++k) {
				if (rank_norm_angle[k] == 0 && pops < this->m_individuals.size() && region_to_sort_ind[k].size() > 0) {
					select_an_ind_subspace_diversity(pops, F2_Fn, region_to_sort_ind, k,pro);
					/*for (int i = 0; i < region_to_sort_ind[k].size() && pops < this->m_individuals.size(); ++i) {
						sub_inds.push_back(std::pair<int, std::pair<int, int>>(region_to_sort_ind[k][i], std::pair<int, int>(k, i)));
					}*/
					num++;
				}
			}

			if (m_observe) {
				m_kdt->m_observe_angle_sub_inner = m_observe_angle_subspace_inner;
				m_observe = false;
			}

			//if (sub_inds.size() == 0)
				//std::cout << "sub_inds is error!" << std::endl;
			//select_an_ind_max_min(pops, F2_Fn, sub_inds, region_to_sort_ind);	//select an Solution from sub_inds		
			change_of_pops.push_back(std::pair<int, std::pair<int, int>>(pops, std::pair<int, int>(num, m_kdt->size_region())));//for test change of pops
			if (test == pops_F1 && pops_F1 < pops) {	//for test the second selection
				m_pop_second_select.resize(pops - pops_F1);
				for (; test < pops; ++test) {
					m_pop_second_select[test - pops_F1] = *(this->m_individuals[test]);
				}
			}

			/*update m_region_sort of kd-tree*/
			for (int k = 0; k < m_kdt->size_region();) {
				if (region_to_sort_ind[k].size() == 0) {
					region_to_sort_ind.erase(region_to_sort_ind.begin() + k);
					m_kdt->region().erase(m_kdt->region().begin() + k);
				}
				else
					++k;
			}
		}

		/*for test change of pops*/
		std::cout << "selection process: " << change_of_pops[0].first << "  ";
		for (int i = 1; i < change_of_pops.size(); ++i) {
			//output Solutions number every cycle
			std::cout << change_of_pops[i].first - change_of_pops[i - 1].first << "(" << change_of_pops[i].second.first << "/" << change_of_pops[i].second.second << ")  ";
		}
		std::cout << std::endl;
	}

	void KDT_MOEA_DE_pop::select_F1(std::vector<IndDE>& F1, int& pops,Problem *pro) {
		m_kdt->updateRegionLabel(F1,pro);	//update m_kdt->m_region_sort
		int size_region = m_kdt->size_region();
		/*calculate Solutions for each subregion*/
		std::vector<std::vector<int>> region_to_sort_ind(size_region);
		int number = update_region_to_sort_ind(F1, region_to_sort_ind);
		if (number != F1.size())
			std::cout << "The matching subregion of Solution is error!" << std::endl;
		/*prediction*/
		if (m_archive.size() >= 2) {
			predict();
			add_pop_to_pred(pro);		//add the non-domianted Solutions of pre-generation population into prediction (its left point in box)
			m_kdt->cal_dds(F1, region_to_sort_ind,pro);
			add_subregion_to_pred(m_kdt->region(), 1,pro);//add the non-dominated(type=0) or all(type=1, sub_dds>0) subregions into prediction (its left lower point in box)
			add_shrink_all_nondominance_to_pred(m_F1);
			//observe m_prediction_kalman_filter
			/*std::cout << "prediction_kalman_filter:" << std::endl;
			for (int i = 0; i < m_prediction_kalman_filter.size(); ++i) {
				for (int j = 0; j < m_prediction_kalman_filter[0].size(); ++j) {
					if (j == 0)
						std::cout << m_prediction_kalman_filter[i][j];
					else
						std::cout << "\t" << m_prediction_kalman_filter[i][j];
				}
				std::cout << std::endl;
			}*/
			//record number of prediction in pre and post of handling
			auto para_map = pro->getParam();
			int pred_method = para_map->get<int>("predict");
			if (pred_method == 1)
				m_prediction_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_size = m_prediction_least_square.size();
			handle_prediction(m_F1,pro);//erase these dominated-prediction-point
			if (pred_method == 1)
				m_prediction_handled_size = m_prediction_kalman_filter.size();
			else if (pred_method == 2)
				m_prediction_handled_size = m_prediction_least_square.size();
		}

		/*select offspring : the first stage*/
		if (size_region < this->m_individuals.size()) {
			/*select an Solution from each subapces, let them near the center of subspace*/
			for (int k = 0; k < size_region; ++k) {
				if (pops < this->m_individuals.size()) {
					select_an_ind_subspace_center(pops, F1, region_to_sort_ind, k,pro);
				}
			}
		}
		else {	//select M Solution from F1, let them near the objective axis
			select_M_ind_F1(pops, F1,pro);
		}

		/*record the first selection*/
		m_pop_first_select.resize(pops);
		for (int i = 0; i < pops; ++i) {
			m_pop_first_select[i] = *(this->m_individuals[i]);
		}

		/*select offspring : the second stage, Max-Min method*/
		while (pops < this->m_individuals.size()) {
			std::vector<Real> min_distance(F1.size(), -1);
			std::pair<Real, int> max_min_id = std::make_pair<Real, int>(-1, -1);
			for (int i = 0; i < F1.size(); ++i) {
				if (F1[i].getType() == -1) {
					min_distance[i] = 1e10;
					for (int p = 0; p < pops; ++p) {
						Real distance = euclideanDistance(F1[i].objective().begin(), F1[i].objective().end(), this->m_individuals[p]->objective().begin());
						if (distance < min_distance[i]) {
							min_distance[i] = distance;
						}
					}
					if (max_min_id.first < min_distance[i]) {
						max_min_id.first = min_distance[i];
						max_min_id.second = i;
					}
				}
			}
			if (max_min_id.second == -1)
				std::cout << "find max_min_id is error!" << std::endl;
			else {
				*(this->m_individuals[pops]) = F1[max_min_id.second];
				F1[max_min_id.second].setType(1);
				pops++;
			}
		}
	}

	//select an Solution from each subapces, let them near the center of subspace
	void KDT_MOEA_DE_pop::select_an_ind_subspace_center(int& pops, std::vector<IndDE>& F1, std::vector<std::vector<int>>& region_to_sort_ind, int k,Problem *pro) {
		if (region_to_sort_ind[k].size() == 1) {
			*(this->m_individuals[pops]) = F1[region_to_sort_ind[k][0]];
			F1[region_to_sort_ind[k][0]].setType(1);
			++pops;
			region_to_sort_ind[k].erase(region_to_sort_ind[k].begin());	//update region_to_sort_ind
			return;
		}
		int M = CAST_CONOP(pro)->numberObjectives();
		std::vector<Real> center(M);
		std::vector<Real> ref(M);
		/*if (!m_kdt->find_reference_point(ref_point, x_k, ref)) {
			for (int j = 0; j < M; ++j) {
				z_ref[j] = m_min_max[j].first;
			}
		}*/
		auto box = m_kdt->get_kdt().getBox(m_kdt->region()[k].first);
		for (int j = 0; j < M; ++j) {
			center[j] = (box[j].first + box[j].second) / 2;
			ref[j] = box[j].first;
		}
		//std::vector<std::vector<Real, Real>> data(region_to_sort_ind[k].size());		
		//for (int i = 0; i < data.size(); ++i) {
		//	data[i][0] = euclidean_distance(ref.begin(), ref.end(), F1[region_to_sort_ind[k][i]].objective().begin());
		//	Real dot_product = 0;
		//	Real ref_center = 0;
		//	Real ref_obj = 0;
		//	for (int j = 0; j < M; ++j) {
		//		dot_product += (center[j] - ref[j]) * (F1[region_to_sort_ind[k][i]].objective(j) - ref[j]);
		//		ref_center += pow(center[j] - ref[j], 2);
		//		ref_obj += pow(F1[region_to_sort_ind[k][i]].objective(j) - ref[j], 2);
		//	}
		//	data[i][1] = dot_product / (ref_center * ref_obj);
		//	data[i][1] = std::acos(data[i][1]);			
		//}
		//std::vector<optimization_mode> opt_mode({ optimization_mode::Minimization, optimization_mode::Minimization }); //min norm and min angle

		std::vector<Real> data(region_to_sort_ind[k].size());
		std::pair<Real, int> small_angle(1e10, -1);
		for (int i = 0; i < data.size(); ++i) {
			Real dot_product = 0;
			Real ref_center = 0;
			Real ref_obj = 0;
			for (int j = 0; j < M; ++j) {
				dot_product += (center[j] - ref[j]) * (F1[region_to_sort_ind[k][i]].objective(j) - ref[j]);
				ref_center += pow(center[j] - ref[j], 2);
				ref_obj += pow(F1[region_to_sort_ind[k][i]].objective(j) - ref[j], 2);
			}
			ref_center = sqrt(ref_center);
			ref_obj = sqrt(ref_obj);
			data[i] = dot_product / (ref_center * ref_obj);
			if ((fabs(data[i]) - 1.0) > 1e-6)
				std::cout << "calculate cos(theta) in select_an_ind_subspace_center is error!" << std::endl;
			if (fabs(data[i] - 1.0) < 1e-6)
				data[i] = acos(1.0);
			else if (fabs(data[i] + 1.0) < 1e-6)
				data[i] = acos(-1.0);
			else
				data[i] = acos(data[i]);

			if (data[i] < small_angle.first) {
				small_angle.first = data[i];
				small_angle.second = i;
			}
		}
		*(this->m_individuals[pops]) = F1[region_to_sort_ind[k][small_angle.second]];
		F1[region_to_sort_ind[k][small_angle.second]].setType(1);
		++pops;
		region_to_sort_ind[k].erase(region_to_sort_ind[k].begin() + small_angle.second);	//update region_to_sort_ind
	}

	//select M Solution from F1, let them near the objective axis
	void KDT_MOEA_DE_pop::select_M_ind_F1(int& pops, std::vector<IndDE>& F1,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		std::vector<Real> z_star(M);
		for (int j = 0; j < M; ++j) {
			z_star[j] = m_kdt->min_max()[j].first;
		}
		std::vector<bool> flag(F1.size(), true);
		for (int j = 0; j < M; ++j) {
			std::pair<Real, int> smallest_angle = std::make_pair<Real, int>(10, -1);
			for (int i = 0; i < F1.size(); ++i) {
				if (flag[i]) {
					Real norm = 0;
					norm = euclideanDistance(z_star.begin(), z_star.end(), F1[i].objective().begin());
					Real angle = std::acos((F1[i].objective()[j] - z_star[j]) / norm); //cos<a,b>=(a*b)/(||a||*||b||)

					if (angle < smallest_angle.first) {	//for minimize, F-z*
						smallest_angle.first = angle;
						smallest_angle.second = i;
					}
				}
			}
			flag[smallest_angle.second] = false;
			*(this->m_individuals[pops]) = F1[smallest_angle.second];
			F1[smallest_angle.second].setType(1);
			++pops;
		}
	}

	//select an Solution from non-dominated subapces, let them keep diversity with selected population
	void KDT_MOEA_DE_pop::select_an_ind_subspace_diversity(int& pops, std::vector<IndDE>& F2_Fn, std::vector<std::vector<int>>& region_to_sort_ind, int k,Problem *pro) {
		if (region_to_sort_ind[k].size() == 1) {
			*(this->m_individuals[pops]) = F2_Fn[region_to_sort_ind[k][0]];
			F2_Fn[region_to_sort_ind[k][0]].setType(1);
			++pops;
			region_to_sort_ind[k].erase(region_to_sort_ind[k].begin());	//update region_to_sort_ind
			return;
		}
		//int selected_size = region_to_sort_selected_ind[k].size();
		int candidate_size = region_to_sort_ind[k].size();

		std::vector<Real> angle(candidate_size);
		std::pair<Real, int> max_angle(-1e10, -1);
		int M = CAST_CONOP(pro)->numberObjectives();
		for (int i = 0; i < candidate_size; ++i) {
			std::pair<int, int> angle_idx(-1, -1);
			std::vector<std::pair<Real, int>> distance(pops, std::pair<Real, int>(-1, -1));
			for (int p = 0; p < pops; ++p) {
				if (Dominance::kDominant == objectiveCompare(F2_Fn[region_to_sort_ind[k][i]].objective(), this->m_individuals[p]->objective(), CAST_CONOP(pro)->optimizeMode()))
					distance[p].first = 1e10;
				else {
					std::vector<Real> ind_obj = this->m_individuals[p]->objective();
					std::vector<Real> candidate_obj = F2_Fn[region_to_sort_ind[k][i]].objective();
					for (int j = 0; j < M; ++j) {	//normalizate
						ind_obj[j] = (ind_obj[j] - m_kdt->min_max()[j].first) / (m_kdt->min_max()[j].second - m_kdt->min_max()[j].first);
						candidate_obj[j] = (candidate_obj[j] - m_kdt->min_max()[j].first) / (m_kdt->min_max()[j].second - m_kdt->min_max()[j].first);
					}
					distance[p].first = euclideanDistance(ind_obj.begin(), ind_obj.end(), candidate_obj.begin());
				}
				distance[p].second = p;
			}
			std::sort(distance.begin(), distance.end(), [](std::pair<Real, int>a, std::pair<Real, int>b) {return a.first < b.first; });
			angle_idx.first = distance[0].second;
			angle_idx.second = distance[1].second;
			std::vector<Real> x_first = this->m_individuals[angle_idx.first]->objective();
			std::vector<Real> x_second = this->m_individuals[angle_idx.second]->objective();
			//record angle in subspace inner
			if (m_observe) {
				std::vector<std::vector<Real>> temp(3, std::vector<Real>(M));
				temp[0] = F2_Fn[region_to_sort_ind[k][i]].objective();
				temp[1] = x_first;	//the first nearest Solution
				temp[2] = x_second;	//the second nearest Solution
				m_observe_angle_subspace_inner.push_back(temp);
			}

			////normalization x_first and x_second
			//for (int j = 0; j < M; ++j) {
			//	x_first[j] = (x_first[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			//	x_second[j] = (x_second[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
			//}
			Real point_product, first_norm, second_norm;//point product (x_k-x_first) and (x_p-x_second), modules of (x_k-x_first), modules of (x_p-x_second)
			point_product = 0;
			for (int j = 0; j < M; ++j) {
				point_product += (x_first[j] - F2_Fn[region_to_sort_ind[k][i]].objective()[j]) * (x_second[j] - F2_Fn[region_to_sort_ind[k][i]].objective()[j]);
			}
			first_norm = euclideanDistance(x_first.begin(), x_first.end(), F2_Fn[region_to_sort_ind[k][i]].objective().begin());
			second_norm = euclideanDistance(x_second.begin(), x_second.end(), F2_Fn[region_to_sort_ind[k][i]].objective().begin());
			Real cos_theta = point_product / (first_norm * second_norm);
			if ((fabs(cos_theta) - 1.0) > 1e-6)
				std::cout << "calculate angle is error!" << std::endl;
			if (fabs(cos_theta - 1.0) < 1e-6)
				angle[i] = acos(1.0);
			else if (fabs(cos_theta + 1.0) < 1e-6)
				angle[i] = acos(-1.0);
			else
				angle[i] = acos(cos_theta);
			if (max_angle.first < angle[i]) {
				max_angle.first = angle[i];
				max_angle.second = i;
			}
		}
		if (max_angle.second == -1) {
			std::cout << "the select_an_ind_subspace is error, or it has special condition!" << std::endl;
		}
		else {
			*(this->m_individuals[pops]) = F2_Fn[region_to_sort_ind[k][max_angle.second]];
			F2_Fn[region_to_sort_ind[k][max_angle.second]].setType(1);
			++pops;
			region_to_sort_ind[k].erase(region_to_sort_ind[k].begin() + max_angle.second);	//update region_to_sort_ind

		}

		/*if (selected_size == 0) {
			select_an_ind_subspace_center(pops, F2_Fn, region_to_sort_ind, k);
		}
		else if (selected_size == 1) {
			std::vector<Real> obj = m_pop[region_to_sort_selected_ind[k][0]].objective();
			std::vector<Real> angle(candidate_size);
			std::pair<Real, int> max_angle(-1e10, -1);
			auto box = m_kdt->get_kdt().get_box(m_kdt->region()[k].first);
			std::vector<Real> left_ref(box.size());
			for (int j = 0; j < box.size(); ++j) {
				left_ref[j] = box[j].first;
			}
			for (int i = 0; i < candidate_size; ++i) {
				Real dot_product = 0;
				Real ref_center = 0;
				Real ref_obj = 0;
				for (int j = 0; j < obj.size(); ++j) {
					dot_product += (obj[j] - left_ref[j]) * (F2_Fn[region_to_sort_ind[k][i]].objective(j) - left_ref[j]);
					ref_center += pow(obj[j] - left_ref[j], 2);
					ref_obj += pow(F2_Fn[region_to_sort_ind[k][i]].objective(j) - left_ref[j], 2);
				}
				ref_center = sqrt(ref_center);
				ref_obj = sqrt(ref_obj);
				angle[i] = dot_product / (ref_center * ref_obj);
				if ((fabs(angle[i]) - 1.0) > 1e-10)
					std::cout << "calculate cos(theta) is error!" << std::endl;
				angle[i] = std::acos(angle[i]);
				if (max_angle.first < angle[i]) {
					max_angle.first < angle[i];
					max_angle.second = i;
				}
			}
			m_pop[pops] = F2_Fn[region_to_sort_ind[k][max_angle.second]];
			F2_Fn[region_to_sort_ind[k][max_angle.second]].set_type(1);
			++pops;
			region_to_sort_ind[k].erase(region_to_sort_ind[k].begin() + max_angle.second);	//update region_to_sort_ind
		}
		else {
			std::vector<Real> angle(candidate_size);
			std::pair<Real, int> max_angle(-1e10, -1);
			int M = global::ms_global->m_problem->objective_size();
			for (int i = 0; i < candidate_size; ++i) {
				std::pair<int, int> angle_idx(-1, -1);
				std::vector<std::pair<Real, int>> distance(selected_size, std::pair<Real, int>(-1, -1));
				for (int p = 0; p < selected_size; ++p) {
					if (dominationship::Dominating == objective_compare(F2_Fn[region_to_sort_ind[k][i]].objective(), m_pop[region_to_sort_selected_ind[k][p]].objective(), global::ms_global->m_problem->opt_mode()))
						distance[p].first = 1e10;
					else {
						std::vector<Real> ind_obj = m_pop[region_to_sort_selected_ind[k][p]].objective();
						//for (int j = 0; j < ind_obj.size(); ++j)	//normalizate ind_obj, i.e., pop[i].objective()
						//	ind_obj[j] = (ind_obj[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
						distance[p].first = euclidean_distance(ind_obj.begin(),ind_obj.end(), F2_Fn[region_to_sort_ind[k][i]].objective().begin());
					}
					distance[p].second = p;
				}
				std::sort(distance.begin(), distance.end(), [](std::pair<Real, int>a, std::pair<Real, int>b) {return a.first < b.first; });
				angle_idx.first = distance[0].second;
				angle_idx.second = distance[1].second;
				std::vector<Real> x_first = m_pop[region_to_sort_selected_ind[k][angle_idx.first]].objective();
				std::vector<Real> x_second = m_pop[region_to_sort_selected_ind[k][angle_idx.second]].objective();
				////normalization x_first and x_second
				//for (int j = 0; j < M; ++j) {
				//	x_first[j] = (x_first[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
				//	x_second[j] = (x_second[j] - m_min_max[j].first) / (m_min_max[j].second - m_min_max[j].first);
				//}
				Real point_product, first_norm, second_norm;//point product (x_k-x_first) and (x_p-x_second), modules of (x_k-x_first), modules of (x_p-x_second)
				point_product = 0;
				for (int j = 0; j < M; ++j) {
					point_product += (x_first[j] - F2_Fn[region_to_sort_ind[k][i]].objective()[j]) * (x_second[j] - F2_Fn[region_to_sort_ind[k][i]].objective()[j]);
				}
				first_norm = euclidean_distance(x_first.begin(),x_first.end(),F2_Fn[region_to_sort_ind[k][i]].objective().begin());
				second_norm = euclidean_distance(x_second.begin(), x_second.end(), F2_Fn[region_to_sort_ind[k][i]].objective().begin());
				Real cos_theta = point_product / (first_norm * second_norm);
				if ((fabs(cos_theta) - 1.0) > 1e-6)
					std::cout << "calculate angle is error!" << std::endl;
				angle[i] = acos(cos_theta);
				if (max_angle.first < angle[i]) {
					max_angle.first = angle[i];
					max_angle.second = i;
				}
			}
			m_pop[pops] = F2_Fn[region_to_sort_ind[k][max_angle.second]];
			F2_Fn[region_to_sort_ind[k][max_angle.second]].set_type(1);
			++pops;
			region_to_sort_ind[k].erase(region_to_sort_ind[k].begin() + max_angle.second);	//update region_to_sort_ind*/

	}

	int KDT_MOEA_DE_pop::update_region_to_sort_ind(std::vector<IndDE>& offspring, std::vector<std::vector<int>>& region_to_sort_ind) {
		int size_region = region_to_sort_ind.size();
		int number = 0;
		for (int k = 0; k < size_region; ++k) {
			for (int i = 0; i < offspring.size(); ++i) {
				int regionIdx = m_kdt->get_kdt().getRegionIdx(offspring[i].objective());
				if (m_kdt->region()[k].first == regionIdx) {
					if (offspring[i].getType() == -1)
						region_to_sort_ind[k].push_back(i); //record only using not selected Solutions  					
				}
			}
			number += region_to_sort_ind[k].size();
		}
		return number;
	}

	int KDT_MOEA_DE_pop::find_non_dominance(std::vector<IndDE>& offspring,Problem *pro) {
		int count_F1 = 0;
		for (size_t i = 0; i < offspring.size(); ++i)
			offspring[i].setFitness(-1);

		std::vector<bool> flag(offspring.size(), true);

		for (size_t j = 0; j < offspring.size(); j++) {
			for (size_t i = 0; i < offspring.size(); i++) {
				if (i == j || !flag[j] || !flag[i]) continue;
				if (offspring[j].dominate(offspring[i], pro)) {
					flag[i] = false;
				}
			}
		}
		for (size_t i = 0; i < offspring.size(); i++) {
			if (flag[i]) {
				offspring[i].setFitness(0);
				count_F1++;
			}
		}
		return count_F1;
	}

	void KDT_MOEA_DE_pop::eval_selection_dd_dominance(Problem *pro) {
		m_pop_rank_is_0.resize(m_number_test);
		//m_kdt->sort_dd(m_offspring);
		m_kdt->sort_dd_unify(m_offspring,pro);
		m_kdt->update_min_max(m_offspring,pro);
		m_kdt->set_kdtree(m_kdt->min_max(),pro);
		int count = 0;
		int rank = -1;
		while (count < m_number_test) {
			rank++;
			for (int i = 0; i < m_offspring.size() && count < m_number_test; ++i) {
				if (m_offspring[i].fitness() == rank)
					m_pop_rank_is_0[count++] = m_offspring[i];
			}
		}

	}
	void KDT_MOEA_DE_pop::eval_selection_1_k_dominance(Problem *pro) {
		m_pop_1_k_dominance.resize(m_number_test);
		m_kdt->update_min_max(m_offspring,pro);
		m_kdt->set_kdtree(m_kdt->min_max(),pro);
		/*(1-k)-dominance------start*/
		//for (int i = 0; i < m_offspring.size(); ++i) {	//test
		//	std::cout << i << "\t";
		//	for (int j = 0; j < m_offspring[i].objective_size(); ++j) {
		//		std::cout << m_offspring[i].objective()[j] << "\t";
		//	}
		//	std::cout << std::endl;
		//}

		//test (1-k)-dominamce
		/*std::vector<std::vector<Real>> test(5, std::vector<Real>(3));
		std::vector<std::vector<Real>*> test_ptr(5);
		std::vector<int> test_rank(5);
		test = { {2,1,3},{1,0,4},{4,2,0},{3,3,2},{0,4,1} };
		for (int i = 0; i < 5; ++i) {
			test_ptr[i] = &test[i];
		}
		std::vector<optimization_mode> test_mode = { optimization_mode::Minimization, optimization_mode::Minimization, optimization_mode::Minimization };
		dominance_1_k(test_ptr, test_rank, test_mode);*/

		m_kdt->sort_1_k_dominance(m_offspring, m_k,pro);
		/*(1-k)-dominance------end*/
		int count = 0;
		int rank1 = -1;
		while (count < m_number_test) {
			rank1++;
			for (int i = 0; i < m_offspring.size() && count < m_number_test; ++i) {
				if (m_offspring[i].fitness() == rank1)
					m_pop_1_k_dominance[count++] = m_offspring[i];
			}
		}
	}

	void KDT_MOEA_DE_pop::eval_selection_ranking_dominance(Problem *pro) {
		m_pop_ranking_dominance.resize(m_number_test);
		m_kdt->update_min_max(m_offspring,pro);
		m_kdt->set_kdtree(m_kdt->min_max(),pro);
		/*ranking-dominance------start*/
		//test ranking-dominamce
		/*std::vector<std::vector<Real>> test(5, std::vector<Real>(3));
		std::vector<std::vector<Real>*> test_ptr(5);
		std::vector<int> test_rank(5);
		test = { {5,200,15},{9,350,14},{4,100,13},{6,270,12},{3,400,11} };
		for (int i = 0; i < 5; ++i) {
			test_ptr[i] = &test[i];
		}
		std::vector<optimization_mode> test_mode = { optimization_mode::Minimization, optimization_mode::Minimization, optimization_mode::Minimization };
		ranking_dominance_R_sum(test_ptr, test_rank, test_mode);*/

		std::vector<std::vector<Real>*> objs;
		for (auto& i : m_offspring)
			objs.emplace_back(&i.objective());
		std::vector<int> rank;
		ranking_dominance_R_sum(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < m_offspring.size(); ++i)
			m_offspring[i].setFitness(rank[i]);
		/*ranking-dominance------end*/
		int count = 0;
		int rank1 = -1;
		while (count < m_number_test) {
			rank1++;
			for (int i = 0; i < m_offspring.size() && count < m_number_test; ++i) {
				if (m_offspring[i].fitness() == rank1)
					m_pop_ranking_dominance[count++] = m_offspring[i];
			}
		}
	}

	int KDT_MOEA_DE_pop::ranking_dominance_R_sum(const std::vector<std::vector<Real>*>& data, std::vector<int>& rank, const std::vector<OptimizeMode>& opt_mode) {
		int popsize = data.size();
		if (popsize == 0) return 0;
		if (rank.size() != popsize)
			rank.resize(popsize);
		std::vector<std::vector<Real>> data_(popsize, std::vector<Real>(data[0]->size()));
		for (int i = 0; i < popsize; ++i) {
			data_[i] = *data[i];
		}
		std::vector<std::vector<int>> Ranks(popsize, std::vector<int>(data[0]->size()));
		std::vector<int> rank1(popsize);

		for (int j = 0; j < data[0]->size(); ++j) {
			for (int i = 0; i < popsize; ++i) {
				rank1[i] = i;
			}
			for (int i = 0; i < data_.size(); ++i) {
				int min = i;
				for (int k = i + 1; k < data_.size(); ++k) {
					if (data_[k][j] < data_[min][j])
						min = k;
				}
				if (min != i) {
					Real temp1 = data_[i][j];
					data_[i][j] = data_[min][j];
					data_[min][j] = temp1;
					int temp2 = rank1[i];
					rank1[i] = rank1[min];
					rank1[min] = temp2;
				}
				//				Ranks[min][j] = i;
			}
			for (int i = 0; i < popsize; ++i) {
				Ranks[rank1[i]][j] = i + 1;
			}
		}
		std::vector<int> rank_(data.size());
		for (int i = 0; i < data.size(); ++i) {
			for (int j = 0; j < data[0]->size(); ++j) {
				rank_[i] += Ranks[i][j];
			}
		}

		for (int i = 0; i < popsize; i++)
			rank[i] = -1;

		int m_curRank = 0;
		//std::vector<int> number;
		int number_;
		int number_sum = 0;
		while (1) {
			int stop_count = 0;
			int min_rank = 1e6;
			number_ = 0;
			for (int k = 0; k < popsize; k++) {
				if (rank[k] == -1 && rank_[k] < min_rank)
					min_rank = rank_[k];
			}
			for (int k = 0; k < popsize; k++) {
				if (rank[k] == -1 && rank_[k] == min_rank) {
					rank[k] = m_curRank;
					number_++;
				}
			}
			number_sum += number_;
			m_curRank++;
			if (number_sum == popsize)
				return m_curRank;
		}
	}

	void KDT_MOEA_DE_pop::find_non_dominance(Problem *pro) {	//find non-dominated Solutions from m_pop, and store them into m_arhive
		m_archive.resize(this->iteration() + 1);
		m_archive[this->iteration()].reserve(this->m_individuals.size());
		std::vector<bool> flag(this->m_individuals.size(), true);
		for (size_t j = 0; j < this->m_individuals.size(); j++) {
			for (size_t i = 0; i < this->m_individuals.size(); i++) {
				if (i == j || !flag[j] || !flag[i]) continue;
				if (this->m_individuals[j]->dominate(*(this->m_individuals[i]),pro)) {
					flag[i] = false;
				}
			}
		}
		for (size_t i = 0; i < this->m_individuals.size(); i++) {
			if (flag[i]) {
				m_archive[this->iteration()].push_back(*(this->m_individuals[i]));
			}
		}
	}

	void KDT_MOEA_DE_pop::update_index_archive() {
		int size = m_archive.size();
		for (int i = 0; i < m_archive[size - 2].size() && i < m_archive[size - 1].size(); ++i) {
			std::pair<int, Real> distance_min;
			distance_min.first = i;
			std::vector<Real> v_obj1 = m_archive[size - 2][0].objective();
			std::vector<Real> v_obj2 = m_archive[size - 1][i].objective();
			distance_min.second = euclideanDistance(v_obj1.begin(), v_obj1.end(), v_obj2.begin());
			for (int j = i + 1; j < m_archive[size - 1].size(); ++j) {
				Real distance_i_j = euclideanDistance(m_archive[size - 2][i].objective().begin(), m_archive[size - 2][i].objective().end(), m_archive[size - 1][j].objective().begin());
				if (distance_i_j < distance_min.second) {
					distance_min.second = distance_i_j;
					distance_min.first = j;
				}
			}
			if (distance_min.first != i) {
				IndDE temp = m_archive[size - 1][i];
				m_archive[size - 1][i] = m_archive[size - 1][distance_min.first];
				m_archive[size - 1][distance_min.first] = temp;
			}
		}
	}

	void KDT_MOEA_DE_pop::predict() { //predict some points according to m_archive by kalman filter and least square
		/*******predict by least square method********/
		if (m_archive.size() >= 3) {
			/*int size = m_archive[m_archive.size() - 3].size() < m_archive[m_archive.size() - 2].size() ? m_archive[m_archive.size() - 3].size() : m_archive[m_archive.size() - 2].size();
			size = size < m_archive[m_archive.size() - 1].size() ? size : m_archive[m_archive.size() - 1].size();*/
			int size = std::min(std::min(m_archive[m_archive.size() - 3].size(), m_archive[m_archive.size() - 2].size()), m_archive[m_archive.size() - 1].size());
			m_prediction_least_square.resize(size);
			std::vector<std::vector<IndDE>> data(size, std::vector<IndDE>(3));
			for (int i = 0; i < size; ++i) {
				data[i][0] = m_archive[m_archive.size() - 3][i];
				data[i][1] = m_archive[m_archive.size() - 2][i];
				data[i][2] = m_archive[m_archive.size() - 1][i];
				m_prediction_least_square[i] = m_least_square.predict_liner_e(data[i]);
				/*std::cout << "i = " << i << std::endl;
				std::cout << data[i][0].objective(0) << "\t" << data[i][0].objective(1) << std::endl;
				std::cout << data[i][1].objective(0) << "\t" << data[i][1].objective(1) << std::endl;
				std::cout << data[i][2].objective(0) << "\t" << data[i][2].objective(1) << std::endl;
				std::cout << "least square prediction result:" << std::endl;
				std::cout << m_prediction_least_square[i][0] << "\t" << m_prediction_least_square[i][1] << std::endl;*/
			}
		}
		/*******predict by Kalman Filter method********/
		if (m_archive.size() >= 2) {
			int size = std::min(m_archive[m_archive.size() - 2].size(), m_archive[m_archive.size() - 1].size());
			m_prediction_kalman_filter.resize(size);
			std::vector<std::vector<IndDE>> data(size, std::vector<IndDE>(2));
			for (int i = 0; i < size; ++i) {
				data[i][0] = m_archive[m_archive.size() - 2][i];
				data[i][1] = m_archive[m_archive.size() - 1][i];
				m_prediction_kalman_filter[i] = m_Kalman_Filter[i].predict(data[i]);
				/*std::cout << "i = " << i << std::endl;
				std::cout << data[i][0].objective(0) << "\t" << data[i][0].objective(1) << std::endl;
				std::cout << data[i][1].objective(0) << "\t" << data[i][1].objective(1) << std::endl;
				std::cout << "kalman filter prediction result:" << std::endl;
				std::cout << m_prediction_kalman_filter[i][0] << "\t" << m_prediction_kalman_filter[i][1] << std::endl;*/
			}
		}
	}

	void KDT_MOEA_DE_pop::add_pop_to_pred(Problem *pro) {	//add the non-dominated Solutions of pre-generation population into prediction (its left point in box)
		if (m_archive.size() >= 3) {
			std::vector<Real> left_box(m_prediction_least_square.back().size());
			for (int i = 0; i < m_archive.back().size(); ++i) {
				int k = m_kdt->get_kdt().getRegionIdx(m_archive.back()[i].objective());
				auto box = m_kdt->get_kdt().getBox(k);
				for (int j = 0; j < box.size(); ++j) {
					left_box[j] = box[j].first;
				}
				if (!(judge_point(left_box, m_prediction_least_square)))
					m_prediction_least_square.push_back(left_box);
			}
		}
		if (m_archive.size() >= 2) {
			std::vector<Real> left_box(m_prediction_kalman_filter.back().size());
			for (int i = 0; i < m_archive.back().size(); ++i) {
				int k = m_kdt->get_kdt().getRegionIdx(m_archive.back()[i].objective());
				auto box = m_kdt->get_kdt().getBox(k);
				for (int j = 0; j < box.size(); ++j) {
					left_box[j] = box[j].first;
				}
				if (!(judge_point(left_box, m_prediction_kalman_filter)))
					m_prediction_kalman_filter.push_back(left_box);
			}
		}

		if (m_archive.size() == 1) {
			std::vector<Real> left_box(CAST_CONOP(pro)->numberObjectives());
			for (int i = 0; i < m_archive.back().size(); ++i) {
				int k = m_kdt->get_kdt().getRegionIdx(m_archive.back()[i].objective());
				auto box = m_kdt->get_kdt().getBox(k);
				for (int j = 0; j < box.size(); ++j) {
					left_box[j] = box[j].first;
				}
				if (!(judge_point(left_box, m_ref_point)))
					m_ref_point.push_back(left_box);
			}
		}
	}

	//add the non-dominated(type=0) or all(type=1, sub_dds>0) subregions into prediction (its left lower point in box)		
	void KDT_MOEA_DE_pop::add_subregion_to_pred(const std::vector<std::pair<int, Real>>& region_dds, int type,Problem *pro) {
		//type = 0, select nondominated subregions
		//type = 1, select subregions of its dds > 0
		Real max_dds = -1;
		if (m_archive.size() >= 2) {
			for (int i = 0; i < region_dds.size(); ++i) {
				if (max_dds < region_dds[i].second)
					max_dds = region_dds[i].second;
			}
		}
		if (m_archive.size() >= 3) {
			std::vector<Real> left_box(m_prediction_least_square.back().size());
			for (int i = 0; i < region_dds.size(); ++i) {
				bool flag = false;
				if (type == 0)
					flag = (fabs(max_dds - region_dds[i].second) < 1e-10);
				else if (type == 1)
					flag = (region_dds[i].second > 0);
				if (flag) {
					auto box = m_kdt->get_kdt().getBox(region_dds[i].first);
					for (int j = 0; j < box.size(); ++j) {
						left_box[j] = box[j].first;
					}
					if (!(judge_point(left_box, m_prediction_least_square)))
						m_prediction_least_square.push_back(left_box);
				}
			}
		}
		if (m_archive.size() >= 2) {
			std::vector<Real> left_box(m_prediction_kalman_filter.back().size());
			for (int i = 0; i < region_dds.size(); ++i) {
				bool flag = false;
				if (type == 0)
					flag = (fabs(max_dds - region_dds[i].second) < 1e-10);
				else if (type == 1)
					flag = (region_dds[i].second > 0);
				if (flag) {
					auto box = m_kdt->get_kdt().getBox(region_dds[i].first);
					for (int j = 0; j < box.size(); ++j) {
						left_box[j] = box[j].first;
					}
					if (!(judge_point(left_box, m_prediction_kalman_filter)))
						m_prediction_kalman_filter.push_back(left_box);
				}
			}
		}
		if (m_archive.size() < 2) {
			std::vector<Real> left_box(CAST_CONOP(pro)->numberObjectives());
			for (int i = 0; i < region_dds.size(); ++i) {
				bool flag = false;
				if (type == 0)
					flag = (fabs(max_dds - region_dds[i].second) < 1e-10);
				else if (type == 1)
					flag = (region_dds[i].second > 0);
				if (flag) {
					auto box = m_kdt->get_kdt().getBox(region_dds[i].first);
					for (int j = 0; j < box.size(); ++j) {
						left_box[j] = box[j].first;
					}
					if (!(judge_point(left_box, m_ref_point)))
						m_ref_point.push_back(left_box);
				}
			}
		}
	}

	//add shrinked non-dominated Solutions into prediction
	void KDT_MOEA_DE_pop::add_shrink_F1_to_pred() {
		std::vector<Real> temp;
		if (m_archive.size() >= 3) {
			for (int i = 0; i < m_F1.size(); ++i) {
				temp = m_F1[i].objective();
				for (int j = 0; j < temp.size(); ++j) {
					//temp[j] = temp[j] - 0.1 * (temp[j] - m_kdt->min_max()[j].first);
					temp[j] = 0.9 * temp[j];
				}
				if (!(judge_point(temp, m_prediction_least_square)))
					m_prediction_least_square.push_back(temp);
			}
		}
		if (m_archive.size() >= 2) {
			for (int i = 0; i < m_F1.size(); ++i) {
				temp = m_F1[i].objective();
				for (int j = 0; j < temp.size(); ++j) {
					//temp[j] = temp[j] - 0.1 * (temp[j] - m_kdt->min_max()[j].first);
					temp[j] = 0.9 * temp[j];
				}
				if (!(judge_point(temp, m_prediction_kalman_filter)))
					m_prediction_kalman_filter.push_back(temp);
			}
		}
		if (m_archive.size() < 2) {
			for (int i = 0; i < m_F1.size(); ++i) {
				temp = m_F1[i].objective();
				for (int j = 0; j < temp.size(); ++j) {
					//temp[j] = temp[j] - 0.1 * (temp[j] - m_kdt->min_max()[j].first);
					temp[j] = 0.9 * temp[j];
				}
				if (!(judge_point(temp, m_ref_point)))
					m_ref_point.push_back(temp);
			}
		}
	}

	//add shrinked m_all_nondominance into prediction
	void KDT_MOEA_DE_pop::add_shrink_all_nondominance_to_pred(const std::vector<IndDE>& all_nondominance) {
		std::vector<Real> temp;
		if (m_archive.size() >= 3) {
			for (int i = 0; i < all_nondominance.size(); ++i) {
				temp = all_nondominance[i].objective();
				for (int j = 0; j < temp.size(); ++j) {
					//temp[j] = temp[j] - 0.1 * (temp[j] - m_kdt->min_max()[j].first);
					temp[j] = 0.9 * temp[j];
				}
				if (!(judge_point(temp, m_prediction_least_square)))
					m_prediction_least_square.push_back(temp);
			}
		}
		if (m_archive.size() >= 2) {
			for (int i = 0; i < all_nondominance.size(); ++i) {
				temp = all_nondominance[i].objective();
				for (int j = 0; j < temp.size(); ++j) {
					//temp[j] = temp[j] - 0.1 * (temp[j] - m_kdt->min_max()[j].first);
					temp[j] = 0.9 * temp[j];
				}
				if (!(judge_point(temp, m_prediction_kalman_filter)))
					m_prediction_kalman_filter.push_back(temp);
			}
		}
		if (m_archive.size() < 2) {
			for (int i = 0; i < all_nondominance.size(); ++i) {
				temp = all_nondominance[i].objective();
				for (int j = 0; j < temp.size(); ++j) {
					//temp[j] = temp[j] - 0.1 * (temp[j] - m_kdt->min_max()[j].first);
					temp[j] = 0.9 * temp[j];
				}
				if (!(judge_point(temp, m_ref_point)))
					m_ref_point.push_back(temp);
			}
		}
	}


	void KDT_MOEA_DE_pop::handle_prediction(const std::vector<IndDE>& offspring,Problem *pro) {	//erase these dominated-prediction-point

		auto para_map = pro->getParam();
		int predict = para_map->get<int>("predict");
		if (m_prediction_kalman_filter.size() > 0) {
			if (predict != 0) {
				std::cout << "kalman_filter: " << std::endl;
				std::cout << "prediction number: " << m_prediction_kalman_filter.size() << "\t";
			}
			//erase these dominated-point in m_prediction_kalman_filter
			/*std::vector<std::vector<Real>*> prediction_kalman_filter(m_prediction_kalman_filter.size());
			for (int i = 0; i < m_prediction_kalman_filter.size(); ++i) {
				prediction_kalman_filter[i] = &m_prediction_kalman_filter[i];
			}
			std::vector<int> rank_prediction_kal;
			nd_sort::fast_sort<Real>(prediction_kalman_filter, rank_prediction_kal, global::ms_global->m_problem->opt_mode());
			for (int i = m_prediction_kalman_filter.size() - 1; i >= 0; --i) {
				if (rank_prediction_kal[i] != 0)
					m_prediction_kalman_filter.erase(m_prediction_kalman_filter.begin() + i);
			}
			if (predict != 0) {
				std::cout << "non-dominance Solution: " << m_prediction_kalman_filter.size() << "\t";
			}*/
			//test whether Solution of offspring dominate Solution of m_prediction, and erase them
			int count = 0;
			for (int i = 0; i < m_prediction_kalman_filter.size(); ) {
				bool flag1 = false;
				for (int j = 0; j < offspring.size(); ++j) {
					Dominance x = objectiveCompare(offspring[j].objective(), m_prediction_kalman_filter[i], CAST_CONOP(pro)->optimizeMode());
					if (x == Dominance::kDominant) {
						flag1 = true;
						break;
					}
				}
				if (flag1) {
					count++;
					m_prediction_kalman_filter.erase(m_prediction_kalman_filter.begin() + i);//erase dominated-prediction point
				}
				else
					i++;
			}
			if (predict != 0) {
				std::cout << "dominated-Solution by offspring: " << count << std::endl;
			}
		}
		if (m_prediction_least_square.size() > 0) {
			if (predict != 0) {
				std::cout << "least_square: " << std::endl;
				std::cout << "prediction number: " << m_prediction_least_square.size() << "\t";
			}
			//erase these dominated-point in m_prediction_least_square
			/*std::vector<std::vector<Real>*> prediction_least_square(m_prediction_least_square.size());
			for (int i = 0; i < m_prediction_least_square.size(); ++i) {
				prediction_least_square[i] = &m_prediction_least_square[i];
			}
			std::vector<int> rank_prediction_ls;
			nd_sort::fast_sort<Real>(prediction_least_square, rank_prediction_ls, global::ms_global->m_problem->opt_mode());
			for (int i = m_prediction_least_square.size() - 1; i >= 0; --i) {
				if (rank_prediction_ls[i] != 0)
					m_prediction_least_square.erase(m_prediction_least_square.begin() + i);
			}
			if (predict != 0) {
				std::cout << "non-dominance Solution: " << m_prediction_least_square.size() << "\t";
			}*/
			//test whether Solution of offspring dominate Solution of m_prediction, and erase them
			int count = 0;
			for (int i = 0; i < m_prediction_least_square.size(); ) {
				bool flag2 = false;
				for (int j = 0; j < offspring.size(); ++j) {
					Dominance x = objectiveCompare(offspring[j].objective(), m_prediction_least_square[i], CAST_CONOP(pro)->optimizeMode());
					if (x == Dominance::kDominant) {
						flag2 = true;
						break;
					}
				}
				if (flag2) {
					count++;
					m_prediction_least_square.erase(m_prediction_least_square.begin() + i); //erase dominated-prediction point
				}
				else
					i++;
			}
			if (predict != 0) {
				std::cout << "dominated-Solution by offspring: " << count << std::endl;
			}
		}

		if (m_ref_point.size() > 0) {
			if (predict != 0) {
				std::cout << "ref_point: " << std::endl;
				std::cout << "refernce point number: " << m_ref_point.size() << "\t";
			}
			//test whether Solution of offspring dominate Solution of m_prediction, and erase them
			int count = 0;
			for (int i = 0; i < m_ref_point.size(); ) {
				bool flag2 = false;
				for (int j = 0; j < offspring.size(); ++j) {
					Dominance x = objectiveCompare(offspring[j].objective(), m_ref_point[i], CAST_CONOP(pro)->optimizeMode());
					if (x == Dominance::kDominant) {
						flag2 = true;
						break;
					}
				}
				if (flag2) {
					count++;
					m_ref_point.erase(m_ref_point.begin() + i); //erase dominated-prediction point
				}
				else
					i++;
			}
			if (predict != 0) {
				std::cout << "dominated-Solution by offspring: " << count << std::endl;
			}
		}

	}

	bool KDT_MOEA_DE_pop::judge_point(const std::vector<Real>& point, const std::vector<std::vector<Real>>& prediction) {//judge whether prediction include point
		bool result = false;
		for (int i = 0; i < prediction.size(); ++i) {
			bool flag = true;
			for (int j = 0; j < point.size(); ++j) {
				if (!(fabs(point[j] - prediction[i][j]) < 1e-10)) {
					flag = false;
					break;
				}
			}
			if (flag) {
				result = true;
				break;
			}
		}
		return result;
	}


	void KDT_MOEA_DE_pop::nondominatedSorting(std::vector<IndDE>& offspring,Problem *pro) {
		std::vector<std::vector<Real>*> objs;
		for (auto& i : offspring)
			objs.emplace_back(&i.objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < offspring.size(); ++i)
			offspring[i].setFitness(rank[i]);
	}

	void KDT_MOEA_DE_pop::eval_selection_NSGAII(PopMODE<>& parent, std::vector<IndDE>& offspring,Problem *pro) {
		std::vector<IndDE*> pop;
		for (auto& i : offspring)
			pop.emplace_back(&i);

		int numobj = CAST_CONOP(pro)->numberObjectives();
		int pops = 0;  //indicate parent population size be 0
		int size = pop.size();
		int rank = 0;
		while (true) {
			int count = 0;
			for (size_t i = 0; i < size; i++)
				if (pop[i]->fitness() == rank)
					count++;
			int size2 = pops + count;
			if (size2 > parent.size()) {
				break;
			}
			for (size_t i = 0; i < size; i++)
				if (pop[i]->fitness() == rank)
				{
					parent[pops] = offspring[i];
					++pops;
				}
			rank++;
			if (pops >= parent.size()) break;
		}
		if (pops < parent.size()) {
			std::vector<int> list;
			// save the Solutions in the overflowed front
			for (size_t i = 0; i < size; i++)
				if (pop[i]->fitness() == rank)
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
			for (size_t j = 0; j < numobj; j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop[list[i]]->objective()[j];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = parent.size() - pops;
			mergeSort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				parent[pops] = offspring[list[idx[i]]];
				++pops;
			}
			density.clear();
			idx.clear();
			list.clear();
		}
	}

	void KDT_MOEA_DE_pop::select_ind_F1(int& pops, std::vector<int>& F1,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		//std::vector<std::vector<Real>> v_f(M);	//cos<a,b>=(a*b)/(||a||*||b||)
		std::vector<Real> z_star(M);
		for (int j = 0; j < M; ++j) {
			z_star[j] = m_kdt->min_max()[j].first;
		}
		std::vector<bool> flag(F1.size(), true);
		for (int j = 0; j < M; ++j) {
			std::pair<Real, int> smallest_angle = std::make_pair<Real, int>(10, -1);
			for (int i = 0; i < F1.size(); ++i) {
				if (flag[i]) {
					Real norm = 0;
					/*for (int k = 0; k < M; ++k) {
						norm += std::pow((m_offspring[F1[i]].objective()[k] - m_kdt->min_max()[k].first), 2);	//for minimize, F-z*
					}
					norm = std::sqrt(norm);*/
					norm = euclideanDistance(z_star.begin(), z_star.end(), m_offspring[F1[i]].objective().begin());
					//Real angle = std::acos((m_offspring[F1[i]].objective()[j] - m_kdt->min_max()[j].first) / norm); //for minimize, F-z*
					Real angle = std::acos((m_offspring[F1[i]].objective()[j] - z_star[j]) / norm);

					if (angle < smallest_angle.first) {
						smallest_angle.first = angle;
						smallest_angle.second = i;
					}
				}
			}
			flag[smallest_angle.second] = false;
			*(this->m_individuals[pops]) = m_offspring[F1[smallest_angle.second]];
			m_offspring[F1[smallest_angle.second]].setType(1);
			++pops;
		}
	}

	void KDT_MOEA_DE_pop::select_ind_F1(int& pops, std::vector<IndDE>& F1,Problem *pro) {
		int M = CAST_CONOP(pro)->numberObjectives();
		//std::vector<std::vector<Real>> v_f(M);	//cos<a,b>=(a*b)/(||a||*||b||)
		std::vector<Real> z_star(M);
		for (int j = 0; j < M; ++j) {
			z_star[j] = m_kdt->min_max()[j].first;
		}
		std::vector<bool> flag(F1.size(), true);
		for (int j = 0; j < M; ++j) {
			std::pair<Real, int> smallest_angle = std::make_pair<Real, int>(10, -1);
			for (int i = 0; i < F1.size(); ++i) {
				if (flag[i]) {
					Real norm = 0;
					/*for (int k = 0; k < M; ++k) {
						norm += std::pow((m_offspring[F1[i]].objective()[k] - m_kdt->min_max()[k].first), 2);	//for minimize, F-z*
					}
					norm = std::sqrt(norm);*/
					norm = euclideanDistance(z_star.begin(), z_star.end(), F1[i].objective().begin());
					//Real angle = std::acos((m_offspring[F1[i]].objective()[j] - m_kdt->min_max()[j].first) / norm); //for minimize, F-z*
					Real angle = std::acos((F1[i].objective()[j] - z_star[j]) / norm);

					if (angle < smallest_angle.first) {
						smallest_angle.first = angle;
						smallest_angle.second = i;
					}
				}
			}
			flag[smallest_angle.second] = false;
			*(this->m_individuals[pops]) = F1[smallest_angle.second];
			F1[smallest_angle.second].setType(1);
			++pops;
		}
	}


	void KDT_MOEA_DE_pop::select_ind_subregion(int& pops, std::vector<int>& vec_inds, int num_max,Problem *pro) {
		if (vec_inds.size() == 1) {
			*(this->m_individuals[pops]) = m_offspring[vec_inds[0]];
			m_offspring[vec_inds[0]].setType(1);
			vec_inds.erase(vec_inds.begin());
			pops++;
			return;
		}

		/*sort for Solutions in vec_inds*/
		std::vector<std::vector<Real>> D(vec_inds.size(), std::vector<Real>(vec_inds.size()));
		std::vector<Real> rank(vec_inds.size(), -1);
		//caculate dominance matrix D
		for (int i = 0; i < vec_inds.size(); ++i) {
			for (int j = i + 1; j < vec_inds.size(); ++j) {
				D[i][j] = dd(m_offspring[vec_inds[i]], m_offspring[vec_inds[j]],pro);
				D[j][i] = -D[i][j];
			}
		}
		//calculate rank of Solution
		std::vector<double> ratio_dd(vec_inds.size());
		for (int i = 0; i < vec_inds.size(); ++i) {
			double val_dominate = 0, val_dominated = 0;
			for (int j = 0; j < vec_inds.size(); ++j) {
				if (D[i][j] > 0)
					val_dominate += D[i][j];
				else
					val_dominated += std::fabs(D[i][j]);
			}
			if (std::fabs(val_dominated) < 1e-6) {
				ratio_dd[i] = 1e6;
				//max_dd = 1e6;
			}
			else if (std::fabs(val_dominate) < 1e-6)
				ratio_dd[i] = -1e6;
			else {
				ratio_dd[i] = val_dominate / val_dominated;
				/*if (max_dd < ratio_dd[i])
				max_dd = ratio_dd[i];*/
			}
		}
		int rank_ = -1;
		int count = 0;
		while (count < vec_inds.size()) {
			rank_++;
			double max_dd = -1e6;
			for (int i = 0; i < vec_inds.size(); ++i) {
				if (rank[i] == -1 && max_dd < ratio_dd[i])
					max_dd = ratio_dd[i];
			}
			for (int i = 0; i < vec_inds.size(); ++i) {
				if (std::fabs(ratio_dd[i] - max_dd) < 1e-6) {
					rank[i] = rank_;
					count++;
					//ratio_dd[i] = 0;
				}
			}
		}
		std::vector<int> inds_rank0;
		for (int i = 0; i < vec_inds.size(); ++i) {
			if (rank[i] == 0) {
				inds_rank0.push_back(i);
				//inds_rank0.push_back(vec_inds[i]);
			}
		}
		if (inds_rank0.size() == 0) {
			std::cout << "error" << std::endl;
			throw "error";
		}
		/*else if (inds_rank0.size() == 1) {
			*(this->m_pop[pops]) = m_offspring[inds_rank0[0]];
			m_offspring[inds_rank0[0]].set_type(1);
			pops++;
		}
		else {
			select an Solution from inds_rank0;
		}*/
		if ((pops + inds_rank0.size()) <= this->m_individuals.size()) {
			int i;
			for (i = 0; i < inds_rank0.size() && num_max > 0; ++i) {
				*(this->m_individuals[pops]) = m_offspring[vec_inds[inds_rank0[i]]];
				m_offspring[vec_inds[inds_rank0[i]]].setType(1);
				pops++;
				num_max--;
			}
			for (int j = i - 1; j >= 0; --j) {
				vec_inds.erase(vec_inds.begin() + inds_rank0[j]);
			}
		}
		else {
			int i;
			for (i = 0; i < inds_rank0.size() && pops < this->m_individuals.size() && num_max > 0; ++i) {
				*(this->m_individuals[pops]) = m_offspring[vec_inds[inds_rank0[i]]];
				m_offspring[vec_inds[inds_rank0[i]]].setType(1);
				pops++;
				num_max--;
			}
			for (int j = i - 1; j >= 0; --j) {
				vec_inds.erase(vec_inds.begin() + inds_rank0[j]);
			}
		}
	}

	//select an Solution by Max-Min method, for sub_inds, the first store index of Solution, the second store index of subspace
	void KDT_MOEA_DE_pop::select_an_ind_max_min(int& pops, std::vector<IndDE>& offspring, std::vector<std::pair<int, std::pair<int, int>>>& sub_inds, std::vector<std::vector<int>>& region_to_sort_ind) {
		if (pops < this->m_individuals.size()) {
			if (sub_inds.size() == 1) {
				*(this->m_individuals[pops]) = offspring[sub_inds[0].first];
				offspring[sub_inds[0].first].setType(1);
				pops++;
				int k1 = sub_inds[0].second.first;
				int k2 = sub_inds[0].second.second;
				region_to_sort_ind[k1].erase(region_to_sort_ind[k1].begin() + k2);
				return;
			}

			std::vector<Real> min_distance(sub_inds.size(), 0);
			std::pair<Real, int> max_id = std::make_pair<Real, int>(0, 0);
			for (int t = 0; t < sub_inds.size(); ++t) {
				min_distance[t] = 10e10;
				for (int i = 0; i < pops; ++i) {
					Real distance = minkowski_distance(offspring[sub_inds[t].first].objective().begin(), offspring[sub_inds[t].first].objective().end(), this->m_individuals[i]->objective().begin(), 2);//p=2, the distance is euclidean distance
					if (distance < min_distance[t])
						min_distance[t] = distance;
				}
				if (max_id.first < min_distance[t]) {
					max_id.first = min_distance[t];
					max_id.second = t;
				}
			}
			*(this->m_individuals[pops++]) = offspring[sub_inds[max_id.second].first];
			offspring[sub_inds[max_id.second].first].setType(1);
			int k1 = sub_inds[max_id.second].second.first;
			int k2 = sub_inds[max_id.second].second.second;
			region_to_sort_ind[k1].erase(region_to_sort_ind[k1].begin() + k2);
		}
	}

	Real KDT_MOEA_DE_pop::dd(IndDE& ind1, IndDE& ind2,Problem *pro) {
		int better, worse;
		better = worse = 0;
		int M = pro->numberObjectives();
		for (int j = 0; j < M; ++j) {
			if (ind1.objective(j) < ind2.objective(j))
				better++;
			else if (ind2.objective(j) < ind1.objective(j))
				worse++;
		}
		return (Real)(better - worse) / M;
	}

	/*void KDT_MOEA_DE::normalization(const std::vector<DE::Solution> &offspring, std::vector<DE::Solution> &offspring_map, const std::vector<std::pair<Real, Real>> &min_max) {
		int M = global::ms_global->m_problem->objective_size();
		for (int i = 0; i < offspring.size(); ++i) {
			offspring_map[i] = offspring[i];
			for (int j = 0; j < M; ++j) {
				offspring_map[i].objective(j) = (offspring_map[i].objective(j) - min_max[j].first) / (min_max[j].second - min_max[j].first);
			}
		}
	}*/


}
