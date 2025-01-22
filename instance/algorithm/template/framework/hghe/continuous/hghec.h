#ifndef OFEC_HGHEC_H
#define OFEC_HGHEC_H

#include "../hghe.h"
#include "../../../../../../utility/kd-tree/kdtree_space.h"

namespace ofec {
	class HGHEC : public HGHE<VariableVector<Real>> {
		OFEC_ABSTRACT_INSTANCE(HGHEC)
	public:
		using SPTree = nanoflann::KDTreeSpace<Real>;

		class Subspace;
		class Hill {
			friend class Subspace;
		protected:
			std::list<Subspace*> m_subspaces;
			Real m_volume = 0;

		public:
			Real potential_explore;

			void addSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace = false);
			void removeSubspace(Subspace *subspace, const SPTree *sp_tree, bool update_subspace = false);
			void merge(Hill *hill, const SPTree *sp_tree, bool update_subspace = false);
			void clear(const SPTree *sp_tree, bool update_subspace = false);
			const std::list<Subspace*>& subspaces() const { return m_subspaces; }
			Real volume() const { return m_volume; }
			virtual Subspace* rouletteWheelSelection(Random *rnd, const SPTree *sp_tree) const;
		};

		class Subspace {
		protected:
			std::list<Hill*> m_assigned_hills;
			std::set<Hill*> m_set_assigned_hills;
			std::list<Subspace*> m_adjacent_subspaces;

		public:
			const size_t id;
			std::vector<const SolutionType*> his_sols;
			std::list<const SolutionType*> his_exploit;
			std::vector<const SolutionType*> his_explore;
			const SolutionType *best_sol = nullptr;
			bool best_in_exploit = false;
			Real obj = 0; // for grouping subspace
			int num_learn = 1;

			Subspace(size_t id) : id(id) {}
			void joinHill(Hill *hill, const SPTree *sp_tree, bool update_hill = false);
			void quitHill(Hill *hill, const SPTree *sp_tree, bool update_hill = false);
			void clearAssignedHills(const SPTree *sp_tree, bool update_hill = false);
			const std::list<Hill*>& assignedHills() const { return m_assigned_hills; }
			const std::set<Hill*>& setassignedHills() const { return m_set_assigned_hills; }
			const std::list<Subspace*>& adjacentSubspaces() const { return m_adjacent_subspaces; }
			virtual Subspace* halve(SPTree *sp_tree);
		};

		void randomVarInHill(VariableVector<Real> &var, Hill *hill);
		bool isVarInHill(const VariableVector<Real> &var, Hill *hill) const;
		const std::list<HGHEC::Hill*>& theHillsLocated(const VariableVector<Real> &var) const;

		void randomSolInHill(SolutionType &sol, Hill *hill);
		bool isSolInHill(const SolutionType &sol, Hill *hill) const;
		const std::list<HGHEC::Hill*>& theHillsLocated(const SolutionType &sol) const;

	protected:
		std::unique_ptr<SPTree> m_sp_tree;
		std::list<std::unique_ptr<Hill>> m_hills;
		std::map<Hill*, size_t> m_ptr_to_id;
		std::map<size_t, Hill*> m_id_to_ptr;
		std::vector<std::unique_ptr<Subspace>> m_subspaces;

		Real m_phi;
		bool m_use_acceleration_mode;
		size_t m_num_sols_thold, m_num_sols_intvl;
		bool m_use_his_explore_only;
		bool m_keep_subspace_valley;

		void addInputParameters();
		void initialize_(Environment *env) override;
		const SolutionType* archiveSolution(const SolutionType &sol, TaskEval task, Environment *env) override;
		void updateInfoSSP(const SolutionType *sol, TaskEval task, Environment *env);

		virtual void initHills(Environment *env);
		virtual void updateHills(Environment *env);
		virtual void getHisSolsInHill(std::vector<const SolutionType *> &his_sols, Hill *hill);
		virtual void identifySeedSolutions(const std::vector<const SolutionType *> &his_sols,
			std::vector<const SolutionType *> &seed_sols, Environment *env);
		virtual void subdivideSpace(const std::vector<const SolutionType *> &seed_sols, Environment *env);
		virtual void assignSubspaceFitness();
		virtual void groupSubspace(Environment *env);

		Subspace* halveSubspace(Subspace *subspace, Environment *env);
		void bilinearInterpolation();
		void findInferiorNeighbors1(size_t center, std::set<size_t> &cluster,
			std::vector<std::map<size_t, size_t>> &min_num_step_to_cluster, size_t id_cluster, size_t num_step);
		void findInferiorNeighbors2(size_t center, std::set<size_t> &cluster,
			std::vector<std::map<size_t, size_t>> &min_num_step_to_cluster, size_t id_cluster, size_t num_step);

		bool isTimeToLearn(Hill *hill);
		void updateNumLearn(size_t num_seed_sols, Hill *hill);

	public:
		const std::unique_ptr<SPTree>& spacePartitionTree() const { return m_sp_tree; }
		const std::list<std::unique_ptr<Hill>>& hills() const { return m_hills; }
		const std::vector<std::unique_ptr<Subspace>>& subspaces() const { return m_subspaces; }
	};
}

#endif // !OFEC_HGHEC_H
