#include "GA_DVRP.h"

namespace ofec {
	GA_DVRP::GA_DVRP() :m_cr(0.8), m_mr(0.05) {}

	void GA_DVRP::crossover(DCMOEA_ind<Solution<DVRP::routes>> & ind1, DCMOEA_ind<Solution<DVRP::routes>> & ind2, Random *rnd) {
		if (rnd->uniform.nextNonStd<Real>(0, 1) > m_cr) return;
		auto &c1 = ind1.variable().m_cus_order;
		//auto &c2 = ind2.variable().m_cus_order;
		int car_id_1 = rnd->uniform.nextNonStd<int>(0, c1.size());
		int car_id_2 = car_id_1;
		while(car_id_2==car_id_1)
			car_id_2 = rnd->uniform.nextNonStd<int>(0, c1.size());
		int pos1_1 = rnd->uniform.nextNonStd<int>(1, c1[car_id_1].size() - 1);
		int pos1_2 = pos1_1;
		if (c1[car_id_1].size() > 3)
			while (pos1_1 == pos1_2)
				pos1_2 = rnd->uniform.nextNonStd<int>(1, c1[car_id_1].size() - 1);

		int pos2_1 = rnd->uniform.nextNonStd<int>(1, c1[car_id_2].size() - 1);
		int pos2_2 = pos2_1;
		if (c1[car_id_2].size() > 3)
			while (pos2_1 == pos2_2)
				pos2_2 = rnd->uniform.nextNonStd<int>(1, c1[car_id_2].size() - 1);
		std::vector<int> customers1_, customers2_;
		if (pos1_1 > pos1_2) {
			auto t = pos1_2; pos1_2 = pos1_1;
			pos1_1 = t;
		}
		if (pos1_1 == pos1_2) {
			customers1_.push_back(c1[car_id_1][pos1_1]);
			c1[car_id_1].erase(c1[car_id_1].begin() + pos1_1);
		}
		else {
			for (int i = pos1_1; i < pos1_2; ++i) {
				customers1_.push_back(c1[car_id_1][i]);
			}
			c1[car_id_1].erase(c1[car_id_1].begin() + pos1_1, c1[car_id_1].begin() + pos1_2);
		}

		if (pos2_1 > pos2_2) {
			auto t = pos2_2; pos2_2 = pos2_1;
			pos2_1 = t;
		}
		if (pos2_1 == pos2_2) {
			customers2_.push_back(c1[car_id_2][pos2_1]);
			c1[car_id_2].erase(c1[car_id_2].begin() + pos2_1);
		}
		else {
			for (int i = pos2_1; i < pos2_2; ++i) {
				customers2_.push_back(c1[car_id_2][i]);
			}
			c1[car_id_2].erase(c1[car_id_2].begin() + pos2_1, c1[car_id_2].begin() + pos2_2);
		}

		c1[car_id_1].insert(c1[car_id_1].begin() + pos1_1, customers2_.begin(), customers2_.end());
		c1[car_id_2].insert(c1[car_id_2].begin() + pos2_1, customers1_.begin(), customers1_.end());
	}

	void GA_DVRP::mutate(DCMOEA_ind<Solution<DVRP::routes>> &ind, Random *rnd) {
		if (rnd->uniform.nextNonStd<Real>(0, 1) > m_mr) return;
		auto &c1 = ind.variable().m_cus_order;
		int car_1 = rnd->uniform.nextNonStd<int>(0, c1.size());
		int car_2 = car_1;
		while(car_2==car_1)
			car_2 = rnd->uniform.nextNonStd<int>(0, c1.size());
		int pos1 = rnd->uniform.nextNonStd<int>(1, c1[car_1].size() - 1);
		int pos2 = rnd->uniform.nextNonStd<int>(1, c1[car_2].size() - 1);
		
		int customer = c1[car_1][pos1];
		c1[car_1][pos1] = c1[car_2][pos2];
		c1[car_2][pos2] = customer;
		
	}
	

}