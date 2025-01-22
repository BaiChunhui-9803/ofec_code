#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Segment.h"
#include "./INCLUDE/BIT.h"

/* An implementation of the binary indexed tree data structure (BIT)
 * proposed in 
 *
 *      N. Mladenovic et al.,
 *      A general variable neighborhood search for the one-commodity
 *      pickup-and-delivery travelling salesman problem.
 *      European Journal of Operational Research, 220(1), 270–285, 2012.
 *
 * The BIT structure can be used to calculate the maximum (or minimum) of
 * an array of info.n elements in O(logn) time. 
 *
 * The following operations are provided
 *     
 *     BIT_make:         Create empty min and max trees.
 *     BIT_update:       Update the trees according to the current tour.
 *     BIT_LoadDiff3Opt,
 *     BIT_LoadDiff4Opt,
 *     BIT_LoadDiff5Opt,
 *     BIT_LoadDiff6Opt: Return the load difference for a proposed
 *                       3-, 4-, 5- or 6-opt move.
 */
//static thread_local unique_ptr<int> info.n, info.p;
//static thread_local unique_ptr<int*> info.MinTree, info.MaxTree, info.L;
namespace LKH {


#undef min
#undef max
	static int min(int a, int b)
	{
		return a < b ? a : b;
	}

	static int max(int a, int b)
	{
		return a > b ? a : b;
	}

	static int minSeq(int a, int b, LKH::LKHAlg *Alg)
	{
		return abs(b - a) == Alg->m_BITinfo.n - 1 ? Alg->Dim : a < b ? a : b;
	}

	static int compare(const void *a, const void *b)
	{
		return *(int *)a - *(int *)b;
	}

	void BIT_Make(int Size, LKH::LKHAlg *Alg)
	{
		if (Alg->ProblemType != LKH::ONE_PDTSP)
			return;
		auto& info = Alg->m_BITinfo;
		info.n = Size;
		for (info.p = 1; (1 << info.p) < info.n; info.p++);
		info.MinTree.clear();
		info.MinTree.resize(1 << (info.p + 1),0);
		info.MaxTree.clear();
		info.MaxTree.resize(1 << (info.p + 1), 0);
		info.L.clear();
		info.L.resize(info.n + 1, 0);

		//assert(info.info.MinTree = (int *)calloc(1 << (info.info.p + 1), sizeof(int)));
		//assert(info.info.MaxTree = (int *)calloc(1 << (info.info.p + 1), sizeof(int)));
		//assert(info.info.L = (int *)malloc((info.info.n + 1) * sizeof(int)));
	}

	static void Build(LKH::LKHAlg::BITstruct& info)
	{
		int i, i1;
		for (i = 1; i <= info.n; i++) {
			i1 = i + (1 << info.p) - 1;
			(info.MinTree)[i1] = (info.MaxTree)[i1] = (info.L)[i];
		}
		for (i = (1 << info.p) - 1; i >= 1; i--) {
			info.MinTree[i] = min(info.MinTree[2 * i], info.MinTree[2 * i + 1]);
			info.MaxTree[i] = max(info.MaxTree[2 * i], info.MaxTree[2 * i + 1]);
		}
	}

	static int BIT_Min(int i, int j, LKH::LKHAlg::BITstruct& info)
	{
		int vmin = std::numeric_limits<int>::max();
		if (i > j)
			return vmin;
		i += (1 << info.p) - 2;
		j += (1 << info.p);
		for (; i / 2 != j / 2; i /= 2, j /= 2) {
			if ((i & 1) == 0)
				vmin = min(vmin, info.MinTree[i + 1]);
			if ((j & 1) != 0)
				vmin = min(vmin, info.MinTree[j - 1]);
		}
		return vmin;
	}

	static int BIT_Max(int i, int j, LKH::LKHAlg::BITstruct& info)
	{
		int vmax = std::numeric_limits<int>::min();
		if (i > j)
			return vmax;
		i += (1 << info.p) - 2;
		j += (1 << info.p);
		for (; i / 2 != j / 2; i /= 2, j /= 2) {
			if ((i & 1) == 0)
				vmax = max(vmax, info.MaxTree[i + 1]);
			if ((j & 1) != 0)
				vmax = max(vmax, info.MaxTree[j - 1]);
		}
		return vmax;
	}

	void BIT_Update(LKH::LKHAlg *Alg)
	{

		/*if (!info.n.get())
			info.n.reset(new int(0));
		if (!info.p.get())
			info.p.reset(new int(0));
		if (!info.MinTree.get())
			info.MinTree.reset(new int*(0));
		if (!info.MaxTree.get())
			info.MaxTree.reset(new int*(0));
		if (!info.L.get())
			info.L.reset(new int*(0));*/

		if (Alg->ProblemType != LKH::ONE_PDTSP)
			return;

		auto& info = Alg->m_BITinfo;
		int Forward = (Alg->Reversed == (Alg->Depot)->Parent->Reversed ? (Alg->Depot)->Suc : (Alg->Depot)->Pred)->Id != Alg->Depot->Id + Alg->DimensionSaved;
		int Load = 0, Seq = 0;
		LKH::LKHAlg::Node *N = Alg->Depot;
		do {
			if (N->Id <= Alg->Dim) {
				N->Seq = ++Seq;
				info.L[Seq] = N->Load = Load += N->Demand;
				Alg->NodeSet[N->Id + Alg->DimensionSaved].Seq = Seq;
				Alg->NodeSet[N->Id + Alg->DimensionSaved].Load = Load;
			}
			N = Forward ? (Alg->Reversed == (N)->Parent->Reversed ? (N)->Suc : (N)->Pred) : (Alg->Reversed == (N)->Parent->Reversed ? (N)->Pred : (N)->Suc);
		} while (N != Alg->Depot);
		Build(info);
	}

	static int LoadDiffKOpt(int *t, int K, LKH::LKHAlg *Alg)
	{
		auto& info = Alg->m_BITinfo;
		int MinLoad = min(BIT_Min(1, t[0],info), BIT_Min(t[2 * K - 1] + 1, Alg->Dim,info));
		int MaxLoad = max(BIT_Max(1, t[0],info), BIT_Max(t[2 * K - 1] + 1, Alg->Dim,info));
		int Diff = 0, i, j;
		for (i = 0; i <= 2 * K - 4; i += 2) {
			Diff += info.L[t[i]] - info.L[t[i + 1]];
			j = t[i + 1] % Alg->Dim + 1;
			MinLoad = min(MinLoad, Diff + min(info.L[j], BIT_Min(j, t[i + 2],info)));
			MaxLoad = max(MaxLoad, Diff + max(info.L[j], BIT_Max(j, t[i + 2],info)));
		}
		return MaxLoad - MinLoad;
	}

	int BIT_LoadDiff3Opt(LKH::LKHAlg::Node * t1, LKH::LKHAlg::Node * t2, LKH::LKHAlg::Node * t3, LKH::LKHAlg::Node * t4,
		LKH::LKHAlg::Node * t5, LKH::LKHAlg::Node * t6, LKH::LKHAlg *Alg)
	{
		if (Alg->ProblemType != LKH::ONE_PDTSP || Alg->Swaps > 0)
			return Alg->Capacity;
		int s[3] = { minSeq(t1->Seq, t2->Seq,Alg),
			minSeq(t3->Seq, t4->Seq, Alg),
			minSeq(t5->Seq, t6->Seq, Alg)
		};
		qsort(s, 3, sizeof(int), compare);
		int t[6] = { s[0], s[1], s[2], s[0], s[1], s[2] };
		return LoadDiffKOpt(t, 3, Alg);
	}

	int BIT_LoadDiff4Opt(LKH::LKHAlg::Node * t1, LKH::LKHAlg::Node * t2, LKH::LKHAlg::Node * t3, LKH::LKHAlg::Node * t4,
		LKH::LKHAlg::Node * t5, LKH::LKHAlg::Node * t6, LKH::LKHAlg::Node * t7, LKH::LKHAlg::Node * t8, LKH::LKHAlg *Alg)
	{
		if (Alg->ProblemType != LKH::ONE_PDTSP || Alg->Swaps > 0)
			return Alg->Capacity;
		int s[4] = { minSeq(t1->Seq, t2->Seq, Alg),
			minSeq(t3->Seq, t4->Seq, Alg),
			minSeq(t5->Seq, t6->Seq, Alg),
			minSeq(t7->Seq, t8->Seq, Alg)
		};
		qsort(s, 4, sizeof(int), compare);
		int t[8] = { s[0], s[2], s[3], s[1], s[2], s[0], s[1], s[3] };
		return LoadDiffKOpt(t, 4, Alg);
	}

	int BIT_LoadDiff5Opt(LKH::LKHAlg::Node * t1, LKH::LKHAlg::Node * t2, LKH::LKHAlg::Node * t3, LKH::LKHAlg::Node * t4,
		LKH::LKHAlg::Node * t5, LKH::LKHAlg::Node * t6, LKH::LKHAlg::Node * t7, LKH::LKHAlg::Node * t8,
		LKH::LKHAlg::Node * t9, LKH::LKHAlg::Node * t10, int Case10, LKH::LKHAlg *Alg)
	{
		if (Alg->ProblemType != LKH::ONE_PDTSP || Alg->Swaps > 0)
			return Alg->Capacity;
		int Forward = (Alg->Reversed == (Alg->Depot)->Parent->Reversed ? (Alg->Depot)->Suc : (Alg->Depot)->Pred)->Id != Alg->Depot->Id + Alg->DimensionSaved;
		int s[5] = { minSeq(t1->Seq, t2->Seq, Alg),
			minSeq(t3->Seq, t4->Seq, Alg),
			minSeq(t5->Seq, t6->Seq, Alg),
			minSeq(t7->Seq, t8->Seq, Alg),
			minSeq(t9->Seq, t10->Seq, Alg)
		};
		qsort(s, 5, sizeof(int), compare);
		if (Case10 == 4) {
			int t[10] = { s[0], s[3], s[4], s[2], s[3],
				s[1], s[2], s[0], s[1], s[4]
			};
			return LoadDiffKOpt(t, 5, Alg);
		}
		if (Case10 == 5) {
			if (Alg->Between_SL(t6, Alg->Depot, t1)) {
				int t[10] = { s[0], s[1], s[2], s[0], s[1],
					s[3], s[4], s[2], s[3], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward ? Alg->Between_SL(t8, Alg->Depot, t5) : Alg->Between_SL(t6, Alg->Depot, t1)) {
				int t[10] = { s[0], s[3], s[4], s[0], s[1],
					s[2], s[3], s[1], s[2], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward ? Alg->Between_SL(t4, Alg->Depot, t7) : Alg->Between_SL(t10, Alg->Depot, t3)) {
				int t[10] = { s[0], s[1], s[2], s[3], s[4],
					s[2], s[3], s[0], s[1], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward ? Alg->Between_SL(t10, Alg->Depot, t3) : Alg->Between_SL(t4, Alg->Depot, t7)) {
				int t[10] = { s[0], s[3], s[4], s[1], s[2],
					s[0], s[1], s[2], s[3], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward || Alg->Between_SL(t8, Alg->Depot, t5)) {
				int t[10] = { s[0], s[2], s[3], s[1], s[2],
					s[3], s[4], s[0], s[1], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			else {
				int t[10] = { s[0], s[3], s[4], s[0], s[1],
					s[2], s[3], s[1], s[2], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
		}
		if (Case10 == 13) {
			if (Forward ? Alg->Between_SL(t8, Alg->Depot, t1) : Alg->Between_SL(t4, Alg->Depot, t7)) {
				int t[10] = { s[0], s[1], s[2], s[3], s[4],
					s[2], s[3], s[0], s[1], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward ? Alg->Between_SL(t4, Alg->Depot, t7) : Alg->Between_SL(t8, Alg->Depot, t1)) {
				int t[10] = { s[0], s[3], s[4], s[1], s[2],
					s[0], s[1], s[2], s[3], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward ? Alg->Between_SL(t6, Alg->Depot, t3) : Alg->Between_SL(t2, Alg->Depot, t9)) {
				int t[10] = { s[0], s[2], s[3], s[1], s[2],
					s[3], s[4], s[0], s[1], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Alg->Between_SL(t10, Alg->Depot, t5)) {
				int t[10] = { s[0], s[1], s[2], s[0], s[1],
					s[3], s[4], s[2], s[3], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			else {
				int t[10] = { s[0], s[3], s[4], s[0], s[1],
					s[2], s[3], s[1], s[2], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
		}
		if (Case10 == 14) {
			int t[10] = { s[0], s[2], s[3], s[0], s[1],
				s[3], s[4], s[1], s[2], s[4]
			};
			return LoadDiffKOpt(t, 5, Alg);
		}
		if (Case10 == 15) {
			if (Forward ? Alg->Between_SL(t8, Alg->Depot, t1) : Alg->Between_SL(t2, Alg->Depot, t5)) {
				int t[10] = { s[0], s[3], s[4], s[1], s[2],
					s[0], s[1], s[2], s[3], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Forward ? Alg->Between_SL(t10, Alg->Depot, t7) : Alg->Between_SL(t6, Alg->Depot, t3)) {
				int t[10] = { s[0], s[2], s[3], s[1], s[2],
					s[3], s[4], s[0], s[1], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			if (Alg->Between_SL(t4, Alg->Depot, t9)) {
				int t[10] = { s[0], s[1], s[2], s[0], s[1],
					s[3], s[4], s[2], s[3], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			else if (Forward ? Alg->Between_SL(t6, Alg->Depot, t3) :
				Alg->Between_SL(t10, Alg->Depot, t7)) {
				int t[10] = { s[0], s[3], s[4], s[0], s[1],
					s[2], s[3], s[1], s[2], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
			else {
				int t[10] = { s[0], s[1], s[2], s[3], s[4],
					s[2], s[3], s[0], s[1], s[4]
				};
				return LoadDiffKOpt(t, 5, Alg);
			}
		}
		return 1;
	}

	int BIT_LoadDiff6Opt(LKH::LKHAlg::Node * t1, LKH::LKHAlg::Node * t2, LKH::LKHAlg::Node * t3, LKH::LKHAlg::Node * t4,
		LKH::LKHAlg::Node * t5, LKH::LKHAlg::Node * t6, LKH::LKHAlg::Node * t7, LKH::LKHAlg::Node * t8,
		LKH::LKHAlg::Node * t9, LKH::LKHAlg::Node * t10, LKH::LKHAlg::Node * t11, LKH::LKHAlg::Node * t12, LKH::LKHAlg *Alg)
	{
		if (Alg->ProblemType != LKH::ONE_PDTSP || Alg->Swaps > 0)
			return Alg->Capacity;
		int s[6] = { minSeq(t1->Seq, t2->Seq,Alg),
			minSeq(t3->Seq, t4->Seq,Alg),
			minSeq(t5->Seq, t6->Seq,Alg),
			minSeq(t7->Seq, t8->Seq,Alg),
			minSeq(t9->Seq, t10->Seq,Alg),
			minSeq(t11->Seq, t12->Seq,Alg)
		};
		qsort(s, 6, sizeof(int), compare);
		int t[12] = { s[0], s[4], s[5], s[3], s[4], s[2], s[3], s[1], s[2],
			s[0], s[1], s[5]
		};
		int r = LoadDiffKOpt(t, 6, Alg);
		return r;
	}
}