// Include directory "../"

#ifndef _LKH_H
#define _LKH_H

/*
 * This header is used by almost all functions of the program. It defines 
 * macros and specifies data structures and function prototypes.
 */

 /*************************************************************************
 * Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
 *************************************************************************
 * Author: Changhe Li
 * Email: changhe.lw@google.com
 * Language: C++
 *************************************************************************
 *  This file is part of OFEC. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.

 K. Helsgaun,
 "An Effective Implementation of the Lin-Kernighan Traveling Salesman Heuristic."
 DATALOGISKE SKRIFTER (Writings on Computer Science), No. 81, 1998,
 Roskilde University.

 The source code comes from K. Helsgaun's homepage.
 *************************************************************************/
 // Created: 29 Oct 2019
 // Last modified:XiaoLong

#undef NDEBUG
#define _CRT_SECURE_NO_WARNINGS




//#include <fcntl.h>
#include <assert.h>
//#include <ctype.h>
//#include <float.h>
//#include <limits.h>
//#include <math.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <time.h>
#include "GainType.h"
#include "Hashing.h"

#include "../../../../../core/exception.h"
#include "../../../../../core/random/newran.h"
#include "../../../../../core/global.h"

/* Macro definitions */


#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

//

namespace LKH {

	enum Types {
		TSP, ATSP, SOP, TSPTW, HCP, CCVRP, CVRP, ACVRP, CVRPTW, RCTVRP,
		RCTVRPTW, VRPB, VRPBTW, PDPTW, PDTSP, PDTSPF, PDTSPL, VRPPD, OVRP,
		ONE_PDTSP, MLP, M_PDTSP, M1_PDTSP, TSPDL, TSPPD, TOUR, HPP, BWTSP, TRP,
		CTSP
	};
	enum CoordTypes { TWOD_COORDS, THREED_COORDS, NO_COORDS };
	enum EdgeWeightTypes {
		EXPLICIT, EUC_2D, EUC_3D, MAX_2D, MAX_3D,
		MAN_2D, MAN_3D, CEIL_2D, CEIL_3D, FLOOR_2D, FLOOR_3D,
		GEO, GEOM, GEO_MEEUS, GEOM_MEEUS, ATT, TOR_2D, TOR_3D,
		XRAY1, XRAY2, SPECIAL
	};
	enum EdgeWeightFormats {
		FUNCTION, FULL_MATRIX, UPPER_ROW, LOWER_ROW,
		UPPER_DIAG_ROW, LOWER_DIAG_ROW, UPPER_COL, LOWER_COL,
		UPPER_DIAG_COL, LOWER_DIAG_COL
	};
	enum CandidateSetTypes { ALPHA, DELAUNAY, NN, POPMUSIC, QUADRANT };
	enum InitialTourAlgorithms {
		BORUVKA, CTSP_ALG, CVRP_ALG, GREEDY, MOORE,
		MTSP_ALG, NEAREST_NEIGHBOR, QUICK_BORUVKA, SIERPINSKI, SOP_ALG,
		TSPDL_ALG, WALK
	};
	enum Objectives { MINMAX, MINMAX_SIZE, MINSUM };
	enum RecombinationTypes { IPT, GPX2 };

	class LKHAlg {
	public:

		/*	typedef struct Node Node;
	typedef struct Candidate Candidate;
	typedef struct Segment Segment;
	typedef struct SSegment SSegment;
	typedef struct SwapRecord SwapRecord;
	typedef struct Constraint Constraint;*/
		struct Node;//typedef struct Node Node;
		struct Candidate;
		struct Segment;
		struct SSegment;
		struct SwapRecord;
		struct Constraint;//typedef struct Constraint Constraint;
		//
		typedef Node* (LKHAlg::* MoveFunction) (Node* t1, Node* t2, GainType* G0,
			GainType* Gain);
		typedef int(LKHAlg::* CostFunction) (Node* Na, Node* Nb);
		typedef GainType(LKHAlg::* PenaltyFunction) (void);
		typedef GainType(LKHAlg::* MergeTourFunction) (void);


		/* The Node structure is used to represent nodes (cities) of the problem */

		struct Node {
			int Id;     /* Number of the node (1...Dimension) */
			int Loc;    /* Location of the node in the heap
						   (zero, if the node is not in the heap) */
			int Rank;   /* During the ascent, the priority of the node.
						   Otherwise, the ordinal number of the node in
						   the tour */
			int V;      /* During the ascent the degree of the node minus 2.
						   Otherwise, the variable is used to mark nodes */
			int LastV;  /* Last value of V during the ascent */
			int Cost;   /* "Best" cost of an edge emanating from the node */
			int NextCost;  /* During the ascent, the next best cost of an edge
							  emanating from the node */
			int PredCost,  /* The costs of the neighbor edges on the current tour */
				SucCost;
			int SavedCost;
			int Pi;     /* Pi-value of the node */
			int BestPi; /* Currently best pi-value found during the ascent */
			int Beta;   /* Beta-value (used for computing alpha-values) */
			int Subproblem; /* Number of the subproblem the node is part of */
			int Sons;   /* Number of sons in the minimum spanning tree */
			int* C = 0;     /* A row in the cost matrix */
			int Special;       /* Is the node a special node in Jonker and
								  Volgenant's mTSP to TSP transformation? */
			int Demand;        /* Customer demand in a CVRP or 1-PDTSP instance */
			int* M_Demand = 0;     /* Table of demands in an M-PDTSP instance */
			int Seq;           /* Sequence number in the current tour */
			int DraftLimit;    /* Draft limit in a TSPDL instance */
			int Load;          /* Accumulated load in the current route */
			Node* Pred = 0, * Suc = 0;  /* Predecessor and successor node in
								  the two-way list of nodes */
			Node* OldPred = 0, * OldSuc = 0; /* Previous values of Pred and Suc */
			Node* BestSuc = 0,     /* Best and next best successor node in the */
				* NextBestSuc = 0; /* currently best tour */
			Node* Dad = 0;         /* Father of the node in the minimum 1-tree */
			Node* Nearest = 0;     /* Nearest node (used in the greedy heuristics) */
			Node* Next = 0; /* Auxiliary pointer, usually to the next node in a list
						   of nodes (e.g., the list of "active" nodes) */
			Node* Prev = 0; /* Auxiliary pointer, usually to the previous node
						   in a list of nodes */
			Node* Mark = 0; /* Visited mark */
			Node* FixedTo1 = 0,    /* Pointers to the opposite end nodes of fixed edges. */
				* FixedTo2 = 0;    /* A maximum of two fixed edges can be incident
								 to a node */
			Node* FixedTo1Saved = 0, /* Saved values of FixedTo1 and FixedTo2 */
				* FixedTo2Saved = 0;
			Node* Head = 0; /* Head of a segment of fixed or common edges */
			Node* Tail = 0; /* Tail of a segment of fixed or common edges */
			Node* InputSuc = 0;    /* Successor in the INPUT_TOUR file */
			Node* InitialSuc = 0;  /* Successor in the INITIAL_TOUR file */
			Node* SubproblemPred = 0; /* Predecessor in the SUBPROBLEM_TOUR file */
			Node* SubproblemSuc = 0;  /* Successor in the SUBPROBLEM_TOUR file */
			Node* SubBestPred = 0; /* The best predecessor node in a subproblem */
			Node* SubBestSuc = 0;  /* The best successor node in a subproblem */
			Node* MergePred = 0;   /* Predecessor in the first MERGE_TOUR file */
			Node** MergeSuc = 0;   /* Successors in the MERGE_TOUR files */
			Node* Added1 = 0, * Added2 = 0; /* Pointers to the opposite end nodes
									  of added edges in a submove */
			Node* Deleted1 = 0, * Deleted2 = 0;  /* Pointers to the opposite end nodes
										   of deleted edges in a submove */
			Node* SucSaved = 0;             /* Saved pointer to successor node */
			Candidate* CandidateSet = 0;    /* Candidate array */
			Candidate* BackboneCandidateSet = 0; /* Backbone candidate array */
			Segment* Parent = 0;   /* Parent segment of a node when the two-level
								  tree representation is used */
			Constraint* FirstConstraint = 0;
			double ServiceTime;
			int Pickup, Delivery;
			int DepotId;     /* Equal to Id if the node is a depot; otherwize 0 */
			double Earliest, Latest;
			int Backhaul;
			int Serial;
			int Color;
			double X, Y, Z;     /* Coordinates of the node */
			double Xc, Yc, Zc;  /* Converted coordinates */
			char Axis;  /* The axis partitioned when the node is part of a KDTree */
			char OldPredExcluded, OldSucExcluded;  /* Booleans used for indicating
													  whether one (or both) of the
													  adjoining nodes on the old tour
													  has been excluded */
		};

		// function for node:
		bool Fixed(Node* a, Node* b) { return ((a)->FixedTo1 == (b) || (a)->FixedTo2 == (b)); }
		bool FixedOrCommon(Node* a, Node* b) { return (Fixed(a, b) || IsCommonEdge(a, b)); }
		bool InBestTour(Node* a, Node* b) { return ((a)->BestSuc == (b) || (b)->BestSuc == (a)); }
		bool InNextBestTour(Node* a, Node* b) {
			return ((a)->NextBestSuc == (b) || (b)->NextBestSuc == (a));
		}
		bool InInputTour(Node* a, Node* b) { return ((a)->InputSuc == (b) || (b)->InputSuc == (a)); }
		bool InInitialTour(Node* a, Node* b) {
			return ((a)->InitialSuc == (b) || (b)->InitialSuc == (a));
		}
		bool Near(Node* a, Node* b) {
			return ((a)->BestSuc ? InBestTour(a, b) : (a)->Dad == (b) || (b)->Dad == (a));
		}

		void Link(Node* a, Node* b) { ((a)->Suc = (b))->Pred = (a); }
		void Follow(Node* b, Node* a)
		{
			Link((b)->Pred, (b)->Suc); Link(b, b); Link(b, (a)->Suc); Link(a, b);
		}
		void Precede(Node* a, Node* b)
		{
			Link((a)->Pred, (a)->Suc); Link(a, a); Link((b)->Pred, a); Link(a, b);
		}
		template<typename T>
		void SLink(T a, T b) { (a)->Suc = (b); (b)->Pred = (a); }


		/* The Candidate structure is used to represent candidate edges */

		struct Candidate {
			Node* To = nullptr;   /* The end node of the edge */
			int Cost=0;   /* Cost (distance) of the edge */
			int Alpha=0;  /* Its alpha-value */
		};

		/* The Segment structure is used to represent the segments in the two-level
		   representation of tours */

		struct Segment {
			char Reversed=0;       /* Reversal bit */
			Node* First = nullptr, * Last = nullptr;  /* First and last node in the segment */
			Segment* Pred = nullptr, * Suc = nullptr; /* Predecessor and successor in the two-way
									list of segments */
			int Rank=0;   /* Ordinal number of the segment in the list */
			int Size=0;   /* Number of nodes in the segment */
			SSegment* Parent = nullptr;    /* The parent super segment */
		};

		struct SSegment {
			char Reversed=0;         /* Reversal bit */
			Segment* First = nullptr, * Last = nullptr; /* The first and last node in the segment */
			SSegment* Pred = nullptr, * Suc = nullptr;  /* The predecessor and successor in the
									  two-way list of super segments */
			int Rank=0;   /* The ordinal number of the segment in the list */
			int Size=0;   /* The number of nodes in the segment */
		};

		/* The SwapRecord structure is used to record 2-opt moves (swaps) */

		struct SwapRecord {
			/*SwapRecord()
			{
				t1 = 0; t2 = 0; t3 = 0; t4 = 0;
			}*/
			Node* t1 = nullptr, * t2 = nullptr, * t3 = nullptr, * t4 = nullptr;    /* The 4 nodes involved in a 2-opt move */
		};

		/* The Constraint structure ts used to represent a precedence constraint
		 * for a SOP instance */

		struct Constraint {
			Node* t1 = nullptr, * t2 = nullptr;    /* Node t1 has to precede node t2 */
			Constraint* Suc = nullptr;  /* The next constraint in the list of all constraints */
			Constraint* Next = nullptr; /* The next constraint in the list of constraints
								 for a node */
		};



		struct Opt3Struct {
			int n__=0;
			int** dist= nullptr;
			int* tour_= nullptr, * pos= nullptr;
			int** neighbor= nullptr;
			int* neighbors= nullptr;
			int reversed=0;
			char* dontLook= nullptr;
			GainType tourLength=0;
		};

		//Sequence related info


		Node** t= nullptr;       /* The sequence of nodes to be used in a move */
		Node** T= nullptr;       /* The currently best t's */
		Node** tSaved= nullptr;  /* For saving t when using the BacktrackKOptMove function */
		std::vector<int> p;         /* The permutation corresponding to the sequence in which
						   the t's occur on the tour */
		std::vector<int> q;         /* The inverse permutation of p */
		int* incl= nullptr;      /* Array: incl[i] == j, if (t[i], t[j]) is an inclusion edge */
		int* cycle= nullptr;     /* Array: cycle[i] is cycle number of t[i] */
		GainType* G= nullptr;    /* For storing the G-values in the BestKOptMove function */
		int K=0;

		// BIT related

		struct BITstruct {
			int n=0, p=0;
			std::vector<int> MinTree ;
			std::vector<int> MaxTree ;
			std::vector<int> L ;
		};

		
		BITstruct m_BITinfo;



		typedef int(*ContainsFunction) (LKHAlg::Node* T, int Q, LKHAlg::Node* N, LKHAlg* Alg);
		typedef int(*BoxOverlapsFunction) (LKHAlg::Node* T, int Q, LKHAlg::Node* N, LKHAlg* Alg);


		struct QuadrantInfo {
			Node** KDTree = nullptr;
			Candidate* CandidateSet = nullptr;
			double* XMin = nullptr, * XMax = nullptr, * YMin = nullptr, * YMax = nullptr, * ZMin = 0, * ZMax = 0;
			int Candidates=0, Radius=0;
			int Level = 0;

			ContainsFunction Contains;
			BoxOverlapsFunction BoxOverlaps;
		};

		QuadrantInfo m_quaInfo;
		

		struct PopMusicInfo {
			Node** node = nullptr;
			Node** node_path = nullptr;
		};

		PopMusicInfo m_popmusic_info;



		typedef struct CVPRSaving {
			int Value=0, i=0, j=0;
		} Saving;


		typedef struct CVPRSaveInfo {
			Saving* S = nullptr;
			int Sets=0, SSize=0;
		};
		

		CVPRSaveInfo m_cvpr_info;




		struct ERXTinfo {
			Node* FirstFree = nullptr;
			int Tabu=0;
		};
		ERXTinfo m_erxtInfo;


		GainType OrdinalTourCost;

		struct Gain32Info {
			Node* s1 = nullptr;
			short OldReversed = 0;
		};
		Gain32Info m_gain32Info;


		typedef void(*CrossoverFunction) (LKH::LKHAlg* Alg);

		struct GenecticInfo {
			int MaxPopulationSize=0; /* The maximum size of the population */
			int PopulationSize=0;    /* The current size of the population */

			CrossoverFunction Crossover;
			int** Population= nullptr;      /* Array of individuals (solution tours) */
			GainType* PenaltyFitness= nullptr;  /* The f itnesslty  (tour penalty) of each i                             individual */
			GainType* Fitness=0;     /* The fitness (tour cost) of each individual */
		};

		GenecticInfo m_geneticInfo;



		struct gpxInfo {
			int n_cities=0, n_cand = 0;
			int n_partitions_size2 = 0, n_partitions_before_fusion = 0,
				n_partitions_after_fusion1 = 0, n_partitions_after_fusion2 = 0,
				n_partitions_after_fusion3 = 0;
			int n_partitions_after_fusion4 = 0, n_partitions_after_fusion5 = 0,
				n_partitions_after_fusionB = 0;
			LKHAlg::Node** Map2Node = nullptr;
		};


		gpxInfo m_gpxInfo;


		struct HeapInfo {
			int HeapCount = 0;           /* Its current number of elements */
			int HeapCapacity = 0;        /* Its capacity */

		};

		HeapInfo m_heapInfo;


		struct PatchCyclesInfo {

			std::unique_ptr<int> CurrentCycle, Patchwork, RecLevel;
		};
		PatchCyclesInfo m_patchInfo;


		struct RandomInfo {
			int a = 0, b = 24, arr[55], initialized = 0;
		};

		RandomInfo m_randomInfo;

		struct BufferInfo {
			char* Buffer = nullptr;
			int MaxBuffer = 0;
		};

		BufferInfo m_bufferInfo;
		
		int PenaltiesRead = 0;


		struct StatisticInfo {
			int TrialsMin = 0, TrialsMax = 0, TrialSum = 0, Successes = 0;
			GainType CostMin = 0, CostMax = 0, CostSum = 0;
			GainType PenaltyMin = 0, PenaltyMax = 0, PenaltySum = 0;
			double TimeMin = 0, TimeMax = 0, TimeSum = 0;

		};

		StatisticInfo m_statistic_info;


		int m_DelaunayLevel = 0;
		

		std::string m_filename;
		char* m_filestream_p = 0;
		//unique_ptr<LKHAlg> mp_algorithm;

		MergeTourFunction MergeWithTour;
		

		GainType m_bestCostInside = 0;
	//	std::vector<int> m_bestSolInside;
	//	std::vector<std::vector<char>> mvv_best_diff_matrix;
		std::vector<std::array<int, 2>> m_bestSolInsideRecord;
		int m_bestSolSameEdge = 0;
		int m_curSameEdge = 0;
		GainType m_curCost = 0;
		std::vector<std::array<int, 2>> m_curSolRecord;
		std::vector<std::vector<char>> mvv_diff_matrix;

		
		int m_max_radius = 0;
		bool m_flag_divideBySol = false;

		unsigned long long m_swap_step = 0;

		void initBestInsideRecorder(const std::vector<int>& initSol, GainType cost);

		void initDiffMatrix(const std::vector<int>& initSol);

		/*
				* Edges (t1,t2) and (t3,t4)
 * are exchanged with edges (t2,t3) and (t4,t1).*/
		void updateCurSolRecord(int t1,int t2,int t3,int t4);

		/*
				* Edges (t1,t2) and (t3,t4)
 * are exchanged with edges (t2,t3) and (t4,t1).*/

		void updateDiffDis(int t1, int t2, int t3, int t4);
		/*
		* 
		* Edges (t1,t2) and (t3,t4) 
 * are exchanged with edges (t2,t3) and (t4,t1).
		* 		SwapStack[Swaps].t1 = t1;
		SwapStack[Swaps].t2 = t2;
		SwapStack[Swaps].t3 = t3;
		SwapStack[Swaps].t4 = t4;
		Swaps++;
		Hash ^= (Rand[t1->Id] * Rand[t2->Id]) ^
			(Rand[t3->Id] * Rand[t4->Id]) ^
			(Rand[t2->Id] * Rand[t3->Id]) ^ (Rand[t4->Id] * Rand[t1->Id]);
		*/
		//void updateBestSolInside(int t1, int t2, int t3, int t4);

		// created by DYY 
		int runGA_tsp(const std::string& filename,
			const std::string& parFile,std::vector<int> & tour);
		int run(const std::string &fileName, std::vector<size_t> &tour);


		int LKlocalSearch(std::vector<int>& cursol,int moveType ) {
			MoveType = moveType;
			auto cost = setCurrentSolution(cursol);
			auto cost2 = LinKernighanLocalSearch();
			RecordBetterTour();
			transferTourType(BetterTour, cursol);
			return cost2;
		}

		void assignedMemory() {
			AllocateStructures();
			CreateCandidateSet();
			InitializeStatistics();
		}

		void setInitialTourAlgorithm(InitialTourAlgorithms type) {
			InitialTourAlgorithm = type;
		}

		void set2optLocalSearchParamenters();
		void readProblem(const std::string& tspDir,
			const std::string& tspname);

	

		char* lkh_strtok(char* str, const char* delim, char** saveptr);
		
	//	int lon_sampling2(const std::string& tspFile);
		bool m_flag = true;

		GainType generateSolRandom();
		GainType setCurrentSolution(const std::vector<int>& sol);
		void setCurrentSolutionFromBest();
		


		void transferSol(std::vector<int>& sol);
		
		//GainType generateSolGreedy();
		GainType LinKernighanLocalSearch(void);
		unsigned calHash(const std::vector<int> &sol);


		void freeERXT();
		void freeGain23();
	//	void freeSequence();
		void freeCreateDelaunay();
		void freeCreateQuadrant();
		void freeGreedyTour();
		void freePatchCycle();
		void freeReadPenalties();
		void freeBuffer();

		void freeAll();



		int AscentCandidates;   /* Number of candidate edges to be associated
						   with each node during the ascent */
		int Asymmetric;         /* Is the instance asymmetric? */
		int BackboneTrials;     /* Number of backbone trials in each run */
		int Backtracking;       /* Specifies whether backtracking is used for
								   the first move in a sequence of moves */
		GainType BestCost=0;      /* Cost of the tour in BestTour */
		GainType BestPenalty=0;   /* Penalty of the tour in BestTour */
		int *BestTour = 0;          /* Table containing best tour found */
		GainType BetterCost=0;    /* Cost of the tour stored in BetterTour */
		GainType BetterPenalty = 0; /* Penalty of the tour stored in BetterTour */
		int *BetterTour = 0;        /* Table containing the currently best tour
								   in a run */
		int BWTSP_B;            /* Number of black nodes in a BWTSP instance */
		int BWTSP_Q;            /* Maximum number of subsequent white nodes in a
								   BWTSP instance */
		int BWTSP_L;            /* Maximum length of any path between two black nodes
								   in a BTWSP instance */
		int CacheMask;  /* Mask for indexing the cache */
		int *CacheVal = 0;  /* Table of cached distances */
		int *CacheSig = 0;  /* Table of the signatures of cached
						   distances */
		int CandidateFiles;     /* Number of CANDIDATE_FILEs */
		int EdgeFiles;          /* Number of EDGE_FILEs */
		int *CostMatrix = 0;        /* Cost matrix */
		GainType CurrentGain = 0;
		GainType CurrentPenalty = 0;
		int DemandDimension;    /* Number of commodities in a M-PDTSP instance */
		Node *Depot = 0;
		int Dimension;          /* Number of nodes in the problem */
		int DimensionSaved;     /* Saved value of Dimension */
		int Dim;                /* DimensionSaved - Salesmen + 1 */
		double DistanceLimit;   /* Maxixim route distance for a CVRP instance */
		double Excess;          /* Maximum alpha-value allowed for any
								   candidate edge is set to Excess times the
								   absolute value of the lower bound of a
								   solution tour */
		int ExtraCandidates;    /* Number of extra neighbors to be added to
								   the candidate set of each node */
		Node *FirstActive = 0, *LastActive = 0; /* First and last node in the list
										   of "active" nodes */
		Constraint *FirstConstraint = 0; /* First constraint in the list of SOP
										precedence constraints */
		Node *FirstNode = 0;        /* First node in the list of nodes */
		Segment *FirstSegment = 0;  /* A pointer to the first segment in the cyclic
								   list of segments */
		SSegment *FirstSSegment = 0;        /* A pointer to the first super segment in
										   the cyclic list of segments */
		int Gain23Used; /* Specifies whether Gain23 is used */
		int GainCriterionUsed;  /* Specifies whether L&K's gain criterion is
								   used */
		double GridSize;        /* Grid size used in toroidal instances */
		int GroupSize;  /* Desired initial size of each segment */
		int SGroupSize; /* Desired initial size of each super segment */
		int Groups;     /* Current number of segments */
		int SGroups;    /* Current number of super segments */
		unsigned Hash;  /* Hash value corresponding to the current tour */
		Node **Heap = 0;    /* Heap used for computing minimum spanning trees */
		HashTable *HTable = 0;      /* Hash table used for storing tours */
		int InitialPeriod;      /* Length of the first period in the ascent */
		int InitialStepSize;    /* Initial step size used in the ascent */
		double InitialTourFraction;     /* Fraction of the initial tour to be
										   constructed by INITIAL_TOUR_FILE edges */
		int Kicks;      /* Specifies the number of K-swap-kicks */
		int KickType;   /* Specifies K for a K-swap-kick */
		char *LastLine=0; /* Last input line */

		double LowerBound;      /* Lower bound found by the ascent */
		int M;          /* The M-value is used when solving an ATSP-
						   instance by transforming it to a STSP-instance */
		int MaxBreadth; /* The maximum number of candidate edges
						   considered at each level of the search for
						   a move */
		int MaxCandidates;      /* Maximum number of candidate edges to be
								   associated with each node */
								   //double DistanceLimit; /* Maxixim route distance for a CVRP instance */
		int MaxMatrixDimension; /* Maximum dimension for an explicit cost matrix */
		int MaxSwaps;   /* Maximum number of swaps made during the
						   search for a move */
		int MaxTrials;  /* Maximum number of trials in each run */
		int MergeTourFiles;     /* Number of MERGE_TOUR_FILEs */
		int MoveType;   /* Specifies the sequantial move type to be used
						   in local search. A value K >= 2 signifies
						   that a k-opt moves are tried for k <= K */
		int MoveTypeSpecial; /* A special (3- or 5-opt) move is used */
		Node *NodeSet = 0;  /* Array of all nodes */
		int Norm;       /* Measure of a 1-tree's discrepancy from a tour */
		int NonsequentialMoveType;      /* Specifies the nonsequential move type to
										   be used in local search. A value
										   L >= 4 signifies that nonsequential
										   l-opt moves are tried for l <= L */
		GainType Optimum = 0;       /* Known optimal tour length.
								   If StopAtOptimum is 1, a run will be
								   terminated as soon as a tour length
								   becomes equal this value */
		int PatchingA;  /* Specifies the maximum number of alternating
						   cycles to be used for patching disjunct cycles */
		int PatchingC;  /* Specifies the maximum number of disjoint cycles to be
						   patched (by one or more alternating cycles) */
		GainType PenaltyGain = 0;
		int Precision;  /* Internal precision in the representation of
						   transformed distances */
		int PredSucCostAvailable; /* PredCost and SucCost are available */
		int POPMUSIC_InitialTour;  /* Specifies whether the first POPMUSIC tour
									  is used as initial tour for LK */
		int POPMUSIC_MaxNeighbors; /* Maximum number of nearest neighbors used
									  as candidates in iterated 3-opt */
		int POPMUSIC_SampleSize;   /* The sample size */
		int POPMUSIC_Solutions;    /* Number of solutions to generate */
		int POPMUSIC_Trials;       /* Maximum trials used for iterated 3-opt */
		unsigned *Rand = 0;           /* Table of random values */
		int Recombination;        /* IPT or GPX2 */
		int RestrictedSearch;     /* Specifies whether the choice of the first
									 edge to be broken is restricted */
		short Reversed; /* Boolean used to indicate whether a tour has
						   been reversed */
		int Run;        /* Current run number */
		int Runs;       /* Total number of runs */
		int Scale;      /* Scale factor for Euclidean and ATT instances */
		double ServiceTime; /* Service time for a CVRP instance */
		int Serial;
		unsigned Seed;  /* Initial seed for random number generation */
		int StopAtOptimum;      /* Specifies whether a run will be terminated if
								   the tour length becomes equal to Optimum */
		int Subgradient;        /* Specifies whether the Pi-values should be
								   determined by subgradient optimization */
		int SubproblemSize;     /* Number of nodes in a subproblem */
		int SubsequentMoveType; /* Specifies the move type to be used for all
								   moves following the first move in a sequence
								   of moves. The value K >= 2 signifies that a
								   K-opt move is to be used */
		int SubsequentMoveTypeSpecial; /* A special (3- or 5-opt) sunbsequent move
										  is used */
		int SubsequentPatching; /* Species whether patching is used for
								   subsequent moves */
		SwapRecord *SwapStack = 0;  /* Stack of SwapRecords */
		int Swaps;      /* Number of swaps made during a tentative move */
		int OldSwaps;   /* Saved number of swaps */
		double TimeLimit;       /* The time limit in seconds for each run */
		int TotalDemand;        /* Sum of demands for a CVRP instance */
		int TraceLevel; /* Specifies the level of detail of the output
						   given during the solution process.
						   The value 0 signifies a minimum amount of
						   output. The higher the value is the more
						   information is given */
		int Trial;      /* Ordinal number of the current trial */
		GainType TSPTW_CurrentMakespanCost = 0;
		int TSPTW_Makespan;

		/* The following variables are read by the functions ReadParameters and
		   ReadProblem: */

		char *ParameterFileName = 0, *ProblemFileName = 0, *PiFileName = 0,
			*TourFileName = 0, *OutputTourFileName = 0, *InputTourFileName = 0,
			**CandidateFileName = 0, **EdgeFileName = 0, *InitialTourFileName = 0,
			*SubproblemTourFileName = 0, **MergeTourFileName = 0,
			*MTSPSolutionFileName = 0, *SINTEFSolutionFileName = 0;
		char *Name = 0, *Type = 0, *EdgeWeightType = 0, *EdgeWeightFormat = 0,
			*EdgeDataFormat = 0, *NodeCoordType = 0, *DisplayDataType = 0;
		int CandidateSetSymmetric, CandidateSetType, Capacity,
			CoordType, DelaunayPartitioning, DelaunayPure,
			ExternalSalesmen,
			ExtraCandidateSetSymmetric, ExtraCandidateSetType,
			InitialTourAlgorithm,
			KarpPartitioning, KCenterPartitioning, KMeansPartitioning,
			MTSPDepot, MTSPMinSize, MTSPMaxSize, MTSPObjective,
			MoorePartitioning,
			PatchingAExtended, PatchingARestricted,
			PatchingCExtended, PatchingCRestricted,
			ProblemType, RiskThreshold,
			RohePartitioning, Salesmen, SierpinskiPartitioning,
			SubproblemBorders, SubproblemsCompressed, WeightType, WeightFormat;

		FILE *ParameterFile = 0, *ProblemFile = 0, *PiFile = 0, *InputTourFile = 0,
			*TourFile = 0, *InitialTourFile = 0, *SubproblemTourFile = 0, **MergeTourFile = 0;
		CostFunction Distance = 0, D = 0,
			C = 0 , c = 0, OldDistance = 0;
		MoveFunction BestMove = 0, BacktrackMove = 0, BestSubsequentMove = 0;
		PenaltyFunction Penalty=0;

		/* Function prototypes: */

		int Distance_1(Node * Na, Node * Nb);
		int Distance_LARGE(Node * Na, Node * Nb);
		int Distance_Asymmetric(Node * Na, Node * Nb);
		int Distance_ATSP(Node * Na, Node * Nb);
		int Distance_ATT(Node * Na, Node * Nb);
		int Distance_CEIL_2D(Node * Na, Node * Nb);
		int Distance_CEIL_3D(Node * Na, Node * Nb);
		int Distance_EXPLICIT(Node * Na, Node * Nb);
		int Distance_EUC_2D(Node * Na, Node * Nb);
		int Distance_EUC_3D(Node * Na, Node * Nb);
		int Distance_FLOOR_2D(Node * Na, Node * Nb);
		int Distance_FLOOR_3D(Node * Na, Node * Nb);
		int Distance_GEO(Node * Na, Node * Nb);
		int Distance_GEOM(Node * Na, Node * Nb);
		int Distance_GEO_MEEUS(Node * Na, Node * Nb);
		int Distance_GEOM_MEEUS(Node * Na, Node * Nb);
		int Distance_MAN_2D(Node * Na, Node * Nb);
		int Distance_MAN_3D(Node * Na, Node * Nb);
		int Distance_MAX_2D(Node * Na, Node * Nb);
		int Distance_MAX_3D(Node * Na, Node * Nb);
		int Distance_MTSP(Node * Na, Node * Nb);
		int Distance_SOP(Node * Na, Node * Nb);
		int Distance_SPECIAL(Node * Na, Node * Nb);
		int Distance_TOR_2D(Node * Na, Node * Nb);
		int Distance_TOR_3D(Node * Na, Node * Nb);
		int Distance_XRAY1(Node * Na, Node * Nb);
		int Distance_XRAY2(Node * Na, Node * Nb);

		int D_EXPLICIT(Node * Na, Node * Nb);
		int D_FUNCTION(Node * Na, Node * Nb);

		int C_EXPLICIT(Node * Na, Node * Nb);
		int C_FUNCTION(Node * Na, Node * Nb);

		int c_ATT(Node * Na, Node *Nb);
		int c_CEIL_2D(Node * Na, Node * Nb);
		int c_CEIL_3D(Node * Na, Node * Nb);
		int c_EUC_2D(Node * Na, Node * Nb);
		int c_EUC_3D(Node * Na, Node * Nb);
		int c_FLOOR_2D(Node * Na, Node * Nb);
		int c_FLOOR_3D(Node * Na, Node * Nb);
		int c_GEO(Node * Na, Node * Nb);
		int c_GEOM(Node * Na, Node * Nb);
		int c_GEO_MEEUS(Node * Na, Node * Nb);
		int c_GEOM_MEEUS(Node * Na, Node * Nb);

		void Activate(Node * t);
		int AddCandidate(Node * From, Node * To, int Cost, int Alpha);
		void AddExtraCandidates(int K, int CandidateSetType, int Symmetric);
		void AddExtraDepotCandidates(void);
		void AddTourCandidates(void);
		void AdjustCandidateSet(void);
		void AdjustClusters(int K, Node ** Center);
		void AllocateSegments(void);
		void AllocateStructures(void);
		GainType Ascent(void);
		Node *Best2OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);
		Node *Best3OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);
		Node *Best4OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);
		Node *Best5OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);
		Node *BestKOptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);
		Node *BestSpecialOptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);

		GainType BestKOptMoveRec(int k, GainType G0, GainType& BestG2);

		int Between(const Node * ta, const Node * tb, const Node * tc);
		int Between_SL(const Node * ta, const Node * tb, const Node * tc);
		int Between_SSL(const Node * ta, const Node * tb, const Node * tc);
		GainType BridgeGain(Node * s1, Node * s2, Node * s3, Node * s4,
			Node * s5, Node * s6, Node * s7, Node * s8,
			int Case6, GainType G);
		Node **BuildKDTree(int Cutoff);
		void ChooseInitialTour(void);
		void Connect(Node * N1, int Max, int Sparse);
		void CandidateReport(void);
		void CreateCandidateSet(void);
		void CreateDelaunayCandidateSet(void);
		void CreateNearestNeighborCandidateSet(int K);
		void CreateNNCandidateSet(int K);
		void Create_POPMUSIC_CandidateSet(int K);
		void CreateQuadrantCandidateSet(int K);
		GainType CTSP_InitialTour(void);
		GainType CVRP_InitialTour(void);
		void eprintf(const char *fmt, ...);
		int Excludable(Node * ta, Node * tb);
		void Exclude(Node * ta, Node * tb);
		int FixedOrCommonCandidates(Node * N);
		GainType FindTour(void);
		void Flip(Node * t1, Node * t2, Node * t3);
		void Flip_SL(Node * t1, Node * t2, Node * t3);
		void Flip_SSL(Node * t1, Node * t2, Node * t3);
		int Forbidden(Node * Na, Node * Nb);
		void FreeCandidateSets(void);
		void FreeSegments(void);
		void FreeStructures(void);
		char *FullName(char * Name, GainType Cost);
		int fscanint(FILE *f, int *v);
		GainType Gain23(void);
		void GenerateCandidates(int MaxCandidates, GainType MaxAlpha, int Symmetric);
		double GetTime(void);
		GainType GreedyTour(void);
		int Improvement(GainType  * Gain, Node * t1, Node * SUCt1);
		void InitializeStatistics(void);

		int IsBackboneCandidate(const Node * ta, const Node * tb);
		int IsCandidate(const Node * ta, const Node * tb);
		int IsCommonEdge(const Node * ta, const Node * tb);
		int IsPossibleCandidate(Node * From, Node * To);
		void KSwapKick(int K);
		GainType LinKernighan(void);


		void Make2OptMove(Node * t1, Node * t2, Node * t3, Node * t4);
		void Make3OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
			Node * t5, Node * t6, int Case);
		void Make4OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
			Node * t5, Node * t6, Node * t7, Node * t8,
			int Case);
		void Make5OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
			Node * t5, Node * t6, Node * t7, Node * t8,
			Node * t9, Node * t10, int Case);
		void MakeKOptMove(int K);
		GainType MergeTourWithBestTour(void);
		GainType MergeWithTourIPT(void);
		GainType MergeWithTourGPX2(void);
		GainType Minimum1TreeCost(int Sparse);
		void MinimumSpanningTree(int Sparse);
		void MTSP2TSP(void);
		GainType MTSP_InitialTour(void);
		int MTSP_DiffSize(void);
		void MTSP_WriteSolution(char * FileName, GainType Penalty, GainType Cost);
		void MTSP_Report(GainType Penalty, GainType Cost);
		void NormalizeNodeList(void);
		void NormalizeSegmentList(void);
		void OrderCandidateSet(int MaxCandidates,
			GainType MaxAlpha, int Symmetric);
		GainType PatchCycles(int k, GainType Gain);
		GainType Penalty_ACVRP(void);
		GainType Penalty_BWTSP(void);
		GainType Penalty_CCVRP(void);
		GainType Penalty_CTSP(void);
		GainType Penalty_CVRP(void);
		GainType Penalty_CVRPTW(void);
		GainType Penalty_MTSP_MINMAX(void);
		GainType Penalty_MTSP_MINMAX_SIZE(void);
		GainType Penalty_MTSP_MINSUM(void);
		GainType Penalty_VRPPD(void);
		GainType Penalty_1_PDTSP(void);
		GainType Penalty_MLP(void);
		GainType Penalty_M_PDTSP(void);
		GainType Penalty_M1_PDTSP(void);
		GainType Penalty_OVRP(void);
		GainType Penalty_PDPTW(void);
		GainType Penalty_PDTSP(void);
		GainType Penalty_PDTSPF(void);
		GainType Penalty_PDTSPL(void);
		GainType Penalty_SOP(void);
		GainType Penalty_SCVRPTW(void);
		GainType Penalty_RCTVRP(void);
		GainType Penalty_TRP(void);
		GainType Penalty_TSPDL(void);
		GainType Penalty_TSPPD(void);
		GainType Penalty_TSPTW(void);
		GainType Penalty_VRPB(void);
		GainType Penalty_VRPBTW(void);
		//GainType Penalty_VRPPD(void);
		void PDPTW_Reduce(void);
		void printff(const char * fmt, ...);
		void PrintParameters(void);
		void PrintStatistics(void);
		unsigned Random(void);
		int ReadCandidates(int MaxCandidates);
		//int ReadCandidates(int MaxCandidates);
		int ReadEdges(int MaxCandidates);
		char *ReadLine(FILE * InputFile);
		void ReadParameters(void);
		
		void setDefaultParaments();
		int ReadPenalties(void);
		void ReadProblem(const std::string& filepath);
		void ReadTour(char * FileName, FILE ** File);
		void RecordBestTour(void);
		void RecordCurTour(std::vector<int>& sol)const;
	
		void RecordBetterTour(void);
		Node *RemoveFirstActive(void);
		void ResetCandidateSet(void);
		void RestoreTour(void);
		int SegmentSize(Node * ta, Node * tb);
		void SINTEF_WriteSolution(char * FileName, GainType Cost);
		GainType SFCTour(int CurveType);
		void SolveCompressedSubproblem(int CurrentSubproblem, int Subproblems,
			GainType * GlobalBestCost);
		void SolveDelaunaySubproblems(void);
		void SolveKarpSubproblems(void);
		void SolveKCenterSubproblems(void);
		void SolveKMeansSubproblems(void);
		void SolveRoheSubproblems(void);
		void SolveSFCSubproblems(void);
		int SolveSubproblem(int CurrentSubproblem, int Subproblems,
			GainType * GlobalBestCost);
		void SolveSubproblemBorderProblems(int Subproblems, GainType * GlobalCost);
		void SolveTourSegmentSubproblems(void);
		GainType SOP_InitialTour(void);
		GainType SOP_RepairTour(void);
		void SOP_Report(GainType Cost);
		int SpecialMove(Node * t1, Node * t2, GainType * G0, GainType * Gain);
		void StatusReport(GainType Cost, double EntryTime, char * Suffix);
		void StoreTour(void);
		void SRandom(unsigned seed);
		void SymmetrizeCandidateSet(void);
		GainType TSPDL_InitialTour(void);
		GainType TSPTW_MakespanCost(void);
		void TSPTW_Reduce(void);
		void TrimCandidateSet(int MaxCandidates);
		void UpdateStatistics(GainType Cost, double Time);
		void VRPB_Reduce(void);
		void WriteCandidates(void);
		void WritePenalties(void);
		void WriteTour(char * FileName, int * Tour, GainType Cost);

		void transferTourType(int* Tour, std::vector<int>& vTour)const;
		void getBetterSolution(std::vector<int>& vTour) const {
			transferTourType(BetterTour, vTour);
		}
	};
//	extern thread_local LKHAlg *lkh_ptr;


	// for functions in lkh



}

#endif
