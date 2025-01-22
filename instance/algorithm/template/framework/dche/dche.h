#ifndef OFEC_DCHE_H
#define OFEC_DCHE_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../../../utility/kd-tree/kdtree_space.h"
#include "../../../../../core/problem/solution.h"

namespace ofec {
	class DCHE : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(DCHE)
	public:
		using SPTree = nanoflann::KDTreeSpace<Real>;

		class Subspace;
		class Hill {
			friend class Subspace;
		protected:
			std::list<Subspace*> m_subspaces;
			Real m_volume = 0;
		public:
			void addSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace = false);
			void removeSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace = false);
			void merge(Hill *hill, const SPTree *sp_tree, bool update_subspace = false);
			void clear(const SPTree *sp_tree, bool update_subspace = false);
			Subspace* rouletteWheelSelection(Random *rnd, const SPTree *sp_tree) const;
			const std::list<Subspace*>& subspaces() const { return m_subspaces; }
			Real volume() const { return m_volume; }
		};

		class Subspace {
		protected:
			const size_t m_id;
			std::list<Hill*> m_assigned_hills;
			std::list<Subspace*> m_adjacent_subspaces;
		public:
			Subspace(size_t id) : m_id(id) {}
			size_t ID() const { return m_id; }
			void joinHill(Hill *hill, const SPTree *sp_tree, bool update_hill = false);
			void quitHill(Hill *hill, const SPTree *sp_tree, bool update_hill = false);
			Subspace* bisect(SPTree *sp_tree, int *dim = nullptr, Real *pivot = nullptr);
			void clearAssignedHills(const SPTree *sp_tree, bool update_hill = false);
			const std::list<Hill*>& assignedHills() const { return m_assigned_hills; }
			const std::list<Subspace*>& adjacentSubspaces() const { return m_adjacent_subspaces; }
		};

	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;

		virtual void updateSeeds(Environment *env);
		virtual void updateHills(Environment *env);

		void randomVariablesInHill(VariableVector<> &var, const Hill *hill) const;
		void randomSolutionInHill(Solution<> &sol, const Hill *hill) const;
		bool isVariablesInHill(const VariableVector<> &var, const Hill *hill) const;
		bool isSolutionInHill(const Solution<> &sol, const Hill *hill) const;

		virtual void identifyCandidates(std::list<const Solution<> *> &candidates, Environment *env) = 0;

		bool m_omniscient_hills;
		bool m_omniscient_seeds;
		bool m_heuristic_split;
		
		std::vector<Solution<>> m_seeds;
		std::unique_ptr<SPTree> m_sp_tree;
		std::list<std::unique_ptr<Hill>> m_hills;
		std::vector<std::unique_ptr<Subspace>> m_subspaces;

	public:
		const std::unique_ptr<SPTree>& spacePartitionTree() const { return m_sp_tree; }
		const std::list<std::unique_ptr<Hill>>& hills()const { return m_hills; }
		const std::vector<std::unique_ptr<Subspace>>& subspaces() const { return m_subspaces; }

	private:
		std::tuple<int, Real> dimensionWithSmallestIntersection(const std::vector<Solution<>*> &sols1, const std::vector<Solution<>*> &sols2, size_t num_dims);
	};
}

#endif // !OFEC_DCHE_H
