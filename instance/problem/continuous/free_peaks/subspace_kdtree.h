#ifndef OFEC_FREE_PEAKS_SUBSPACE_KDTREE_H
#define OFEC_FREE_PEAKS_SUBSPACE_KDTREE_H

#include "subproblem/subproblem.h"
#include "../../../../utility/kd-tree/kdtree_space.h"

namespace ofec::free_peaks {
	struct SubspaceKDTree {
		struct Box {
			std::list<std::shared_ptr<Box>> children;
			std::string name;
			Real ratio = 1.0;
			Box(const std::string &name, Real ratio = 1.0) : name(name), ratio(ratio) {}
		};
		std::shared_ptr<Box> root_box;
		std::map<std::string, std::pair<std::shared_ptr<Box>, std::unique_ptr<Subproblem>>> name_box_subproblem;
		std::unique_ptr<nanoflann::KDTreeSpace<Real>> tree;

		void createTree(const std::vector<std::pair<Real, Real>> &ranges, const std::vector<std::pair<std::string, std::vector<std::pair<std::string, Real>>>> &tree_name);
		std::vector<std::pair<std::string, std::vector<std::pair<std::string, Real>>>> treeData() const;
		void clear();
	};
}

#endif // !OFEC_FREE_PEAKS_SUBSPACE_KDTREE_H