/* -------------------------------ReferencePoint-------------------------------------   
    Reference points play very important roles in NSGA-III. Solutions in the population
    are associated with reference points, and the survivors in the environmental selection
    are determined based on the niche count of the reference points.
   
    Check Algorithms 1-4 in the orignal paper for the usage of reference points.
------------------------------------------------------------------------------------*/

#ifndef OFEC_REF_POINT_H
#define OFEC_REF_POINT_H

#include <vector>
#include <utility>
#include "math_aux.h"
#include "../../../../../../utility/random/newran.h"

namespace ofec {
	class RefPoint {
	public:
		explicit RefPoint(size_t s) : mv_position_(s), m_member_size_(0) {}
		RefPoint() {}
		~RefPoint();
		const std::vector<Real>& pos() const { return mv_position_; }
		std::vector<Real>& pos() { return mv_position_; }
		size_t memberSize() const { return m_member_size_; }
		bool hasPotentialMember() const { return !mv_potential_members_.empty(); }
		void clear();
		void addMember();
		void addPotentialMember(size_t member_ind, Real distance);
		int findClosestMember() const;
		int randomMember(Random *rnd) const;
		void removePotentialMember(size_t member_ind);
	private:
		std::vector<Real> mv_position_;
		// pair<indices of Solutions in the population, distance>
		// note. only the data of Solutions in the last considered front will be stored.
		std::vector<std::pair<size_t, Real> > mv_potential_members_;
		size_t m_member_size_;
	};

	namespace reference_point {

		void generateRecursive(std::vector<RefPoint>* rps, RefPoint* pt, size_t number_objectives, size_t left, size_t total, size_t element);

		// GenerateReferencePoints():
		// Given the number of objectives (M) and the number of divisions (p), generate the set of 
		// reference points. Check Section IV-B and equation (3) in the original paper.
		void generateRefPoints(std::vector<RefPoint>* rps, size_t M, const std::vector<size_t>& p);

		// Associate():
		// Associate Solutions in the population with reference points.
		// Check Algorithm 3 in the original paper.
		void associate(std::vector<RefPoint>* prps, const std::vector<std::vector<Real> >& conv_obj, const std::vector<std::vector<int> >& fronts);
	}
}

#endif // !OFEC_REF_POINT_H

