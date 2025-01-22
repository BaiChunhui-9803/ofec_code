#include "reference_point.h"
#include <limits>


namespace ofec {
	RefPoint::~RefPoint() {
		mv_potential_members_.clear();
		mv_position_.clear();
	}

	void RefPoint::clear() {
		m_member_size_ = 0;
		mv_potential_members_.clear();
	}

	void RefPoint::addMember() {
		m_member_size_ += 1;
	}

	void RefPoint::addPotentialMember(size_t member_ind, Real distance) {
		mv_potential_members_.push_back(std::move(std::make_pair(member_ind, distance)));
	}

	int RefPoint::findClosestMember() const {
		Real min_dist = std::numeric_limits<Real>::max();
		int min_indv = -1;
		for (size_t i = 0; i < mv_potential_members_.size(); i += 1)
		{
			if (mv_potential_members_[i].second < min_dist)
			{
				min_dist = mv_potential_members_[i].second;
				min_indv = mv_potential_members_[i].first;
			}
		}
		return min_indv;
	}

	int RefPoint::randomMember(Random *rnd) const {
		if (mv_potential_members_.size() > 0)
			return mv_potential_members_[rnd->uniform.nextNonStd(0, (int)mv_potential_members_.size())].first;
		else
			return -1;
	}

	void RefPoint::removePotentialMember(size_t member_ind) {
		for (size_t i = 0; i < mv_potential_members_.size(); i += 1) {
			if (mv_potential_members_[i].first == member_ind) {
				mv_potential_members_.erase(mv_potential_members_.begin() + i);
				return;
			}
		}
	}

	namespace reference_point {
		void generateRecursive(std::vector<RefPoint>* rps, RefPoint* pt, size_t number_objectives, size_t left, size_t total, size_t element) {
			if (element == number_objectives - 1) {
				pt->pos()[element] = static_cast<Real>(left) / total;
				rps->push_back(*pt);
			}
			else {
				for (size_t i = 0; i <= left; i += 1) {
					pt->pos()[element] = static_cast<Real>(i) / total;
					generateRecursive(rps, pt, number_objectives, left - i, total, element + 1);
				}
			}
		}

		void generateRefPoints(std::vector<RefPoint>* rps, size_t M, const std::vector<size_t>& p) {
			RefPoint pt(M);
			generateRecursive(rps, &pt, M, p[0], p[0], (size_t)0);
			if (p.size() > 1) { // two layers of reference points (Check Fig. 4 in NSGA-III paper)
				std::vector<RefPoint> inside_rps;
				generateRecursive(&inside_rps, &pt, M, p[1], p[1], (size_t)0);
				Real center = 1.0 / M;
				for (size_t i = 0; i < inside_rps.size(); i += 1) {
					for (size_t j = 0; j < inside_rps[i].pos().size(); j += 1)
						inside_rps[i].pos()[j] = (center + inside_rps[i].pos()[j]) / 2; // (k=num_divisions/M, k, k, ..., k) is the center point
					rps->push_back(inside_rps[i]);
				}
			}
		}

		void associate(std::vector<RefPoint>* prps, const std::vector<std::vector<Real>>& conv_obj, const std::vector<std::vector<int>>& fronts) {
			std::vector<RefPoint>& rps = *prps;
			for (size_t t = 0; t < fronts.size(); t += 1) {
				for (size_t i = 0; i < fronts[t].size(); i += 1) {
					size_t min_rp = rps.size();
					Real min_dist = std::numeric_limits<Real>::max();
					for (size_t r = 0; r < rps.size(); r += 1) {
						Real d = math_aux::perpendicularDistance(rps[r].pos(), conv_obj[fronts[t][i]]);
						if (d < min_dist) {
							min_dist = d;
							min_rp = r;
						}
					}
					if (t + 1 != fronts.size()) // associating members in St/Fl (only counting)
						rps[min_rp].addMember();
					else
						rps[min_rp].addPotentialMember(fronts[t][i], min_dist);
				} // for - members in front
			} // for - fronts
		}
	}
}