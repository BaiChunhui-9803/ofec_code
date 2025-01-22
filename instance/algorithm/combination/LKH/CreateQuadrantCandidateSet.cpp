#include "./INCLUDE/LKH.h"

/*
 * This file contains the CreateQuadrantCandidateSet function
 * and the CreateNearestNeighborCandidateSet function.
 */
namespace LKH {
	static void NearestQuadrantNeighbors(LKHAlg::Node * N, int Q, int K, LKHAlg *Alg);
	static int Contains2D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg);
	static int Contains3D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg);
	static int BoxOverlaps2D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg);
	static int BoxOverlaps3D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg);
	static void ComputeBounds(int start, int end, LKHAlg *Alg);



	//static thread_local unique_ptr<LKHAlg::Node *> info.KDTree;
	//static thread_local unique_ptr<LKHAlg::Candidate> info.CandidateSet;
	//static thread_local unique_ptr<double> info.XMin, info.XMax, info.YMin, info.YMax, info.ZMin, info.ZMax;
	//static thread_local unique_ptr<int> info.Candidates, info.Radius, Level_;



	/*
	 * The CreateQuadrantCandidateSet function creates for each node
	 * a candidate set consisting of the K/L least costly neighbor edges
	 * in each of the L geometric quadrants around the node, where L
	 * is 4 for 2-D instances, and 8 for 3-D instances.
	 * If these totals less than K nodes, the candidate set is
	 * augmented by the nearest remaining nodes overall to bring
	 * the total up to K.
	 *
	 * The function is called from the CreateCandidateSet and
	 * AddExtraCandidates functions.
	 */

	void LKHAlg::freeCreateQuadrant()
	{

		auto& info = m_quaInfo;
		/*info.Candidates.reset();
		info.Radius.reset();
		Level_.reset();*/
		free(info.KDTree);
		free(info.CandidateSet);
		free(info.XMin);
		free(info.XMax);
		free(info.YMin);
		free(info.YMax);
		free(info.ZMin);
		free(info.ZMax);
		info.KDTree = 0;
		info.CandidateSet = 0;
		info.Candidates = 0;
		info.XMax = 0;
		info.YMin = 0;
		info.YMax = 0;
		info.ZMin = 0;
		info.ZMax = 0;
		info.Radius = 0;
		info.Level = 0;
	}

	void LKHAlg::CreateQuadrantCandidateSet(int K)
	{

		auto& info = m_quaInfo;



		Node *From, *To;
		Candidate *NFrom;
		int L, Q, CandPerQ, Added, Count, i;

		if (K <= 0)
			return;
		if (TraceLevel >= 2)
			printff("Creating quadrant candidate set ... ");
		/*if (!Level_.get())
		{
			info.Candidates.reset(new int(0));
			info.Radius.reset(new int(0));
			Level_.reset(new int(0));
		}*/
		info.KDTree = BuildKDTree(1);
		/*info.XMin.reset(new double[1 + DimensionSaved]);
		info.XMax.reset(new double[1 + DimensionSaved]);
		info.YMin.reset(new double[1 + DimensionSaved]);
		info.YMax.reset(new double[1 + DimensionSaved]);*/
		assert(info.XMin =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		assert(info.XMax =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		assert(info.YMin =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		assert(info.YMax =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		if (CoordType == THREED_COORDS) {
			/*info.ZMin.reset(new double[1 + DimensionSaved]);
			info.ZMax.reset(new double[1 + DimensionSaved]);*/
			assert(info.ZMin =
				(double *)malloc((1 + DimensionSaved) * sizeof(double)));
			assert(info.ZMax =
				(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		}
		ComputeBounds(0, Dimension - 1, this);
		/*info.Contains.reset(new ContainsFunction(CoordType == THREED_COORDS ? Contains3D : Contains2D));
		info.BoxOverlaps.reset(new BoxOverlapsFunction(CoordType == THREED_COORDS ? BoxOverlaps3D : BoxOverlaps2D));*/
		info.Contains = CoordType == THREED_COORDS ? Contains3D : Contains2D;
		info.BoxOverlaps =
			CoordType == THREED_COORDS ? BoxOverlaps3D : BoxOverlaps2D;
		L = CoordType == THREED_COORDS ? 8 : 4;
		CandPerQ = K / L;
		//info.CandidateSet.reset(new Candidate[K + 1]);
		assert(info.CandidateSet =
			(Candidate *)malloc((K + 1) * sizeof(Candidate)));

		From = FirstNode;
		do {
			Count = 0;
			for (NFrom = From->CandidateSet; NFrom && NFrom->To; NFrom++)
				if (FixedOrCommon(From, NFrom->To) && ++Count == 2)
					break;
			if (Count == 2)
				continue;
			Added = 0;
			for (Q = 1; Q <= L; Q++) {
				NearestQuadrantNeighbors(From, Q, CandPerQ, this);
				for (i = 0; i < info.Candidates; i++) {
					To = info.CandidateSet[i].To;
					if (AddCandidate(From, To, (this->*D)(From, To), 1))
						Added++;
				}
			}
			if (K > Added) {
				NearestQuadrantNeighbors(From, 0, K - Added, this);
				for (i = 0; i < info.Candidates; i++) {
					To = info.CandidateSet[i].To;
					AddCandidate(From, To, (this->*D)(From, To), 2);
				}
			}
		} while ((From = From->Suc) != FirstNode);

		/*info.CandidateSet.reset();
		info.KDTree.reset();
		info.XMin.reset();
		info.XMax.reset();
		info.YMin.reset();
		info.YMax.reset();*/
		free(info.CandidateSet);
		free(info.KDTree);
		free(info.XMin);
		free(info.XMax);
		free(info.YMin);
		free(info.YMax);
		if (CoordType == THREED_COORDS) {
			/*info.ZMin.reset();
			info.ZMax.reset();*/
			free(info.ZMin);
			free(info.ZMax);
		}
		if (info.Level == 0 &&
			(WeightType == GEO || WeightType == GEOM ||
				WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS)) {
			Candidate **SavedCandidateSet;
			assert(SavedCandidateSet =
				(Candidate **)malloc((1 + DimensionSaved) *
					sizeof(Candidate *)));
			if (TraceLevel >= 2)
				printff("done\n");
			From = FirstNode;
			while ((From = From->Suc) != FirstNode)
				if ((From->Y > 0) != (FirstNode->Y > 0))
					break;
			if (From != FirstNode) {
				/* Transform longitude (180 and -180 map to 0) */
				From = FirstNode;
				do {
					SavedCandidateSet[From->Id] = From->CandidateSet;
					From->CandidateSet = 0;
					From->Zc = From->Y;
					if (WeightType == GEO || WeightType == GEO_MEEUS)
						From->Y =
						(int)From->Y + 5.0 * (From->Y -
						(int)From->Y) / 3.0;
					From->Y += From->Y > 0 ? -180 : 180;
					if (WeightType == GEO || WeightType == GEO_MEEUS)
						From->Y =
						(int)From->Y + 3.0 * (From->Y -
						(int)From->Y) / 5.0;
				} while ((From = From->Suc) != FirstNode);
				info.Level++;
				CreateQuadrantCandidateSet(K);
				info.Level--;
				From = FirstNode;
				do
					From->Y = From->Zc;
				while ((From = From->Suc) != FirstNode);
				do {
					Candidate *QCandidateSet = From->CandidateSet;
					From->CandidateSet = SavedCandidateSet[From->Id];
					for (NFrom = QCandidateSet; (To = NFrom->To); NFrom++)
						AddCandidate(From, To, NFrom->Cost, NFrom->Alpha);
					free(QCandidateSet);
				} while ((From = From->Suc) != FirstNode);
				free(SavedCandidateSet);
			}
		}
		if (info.Level == 0) {
			ResetCandidateSet();
			AddTourCandidates();
			if (CandidateSetSymmetric)
				SymmetrizeCandidateSet();
			if (TraceLevel >= 2)
				printff("done\n");
		}
	}

	/*
	 * The CreateNearestNeighborCandidateSet function creates for each node
	 * a candidate set consisting of the K least costly neighbor edges.
	 *
	 * The function is called from the CreateCandidateSet and
	 * AddExtraCandidates functions.
	 */

	void LKHAlg::CreateNearestNeighborCandidateSet(int K)
	{
		auto& info = m_quaInfo;
		Node *From, *To;
		int i;

		if (TraceLevel >= 2)
			printff("Creating nearest neighbor candidate set ... ");
		/*if (!Level_.get())
		{
			info.Candidates.reset(new int(0));
			info.Radius.reset(new int(0));
			Level_.reset(new int(0));
		}*/
		/*info.KDTree = BuildKDTree(1);*/
		/*info.XMin.reset(new double[1 + DimensionSaved]);
		info.XMax.reset(new double[1 + DimensionSaved]);
		info.YMin.reset(new double[1 + DimensionSaved]);
		info.YMax.reset(new double[1 + DimensionSaved]);*/
		info.KDTree = BuildKDTree(1);
		assert(info.XMin =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		assert(info.XMax =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		assert(info.YMin =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		assert(info.YMax =
			(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		if (CoordType == THREED_COORDS) {
			/*info.ZMin.reset(new double[1 + DimensionSaved]);
			info.ZMax.reset(new double[1 + DimensionSaved]);*/
			assert(info.ZMin =
				(double *)malloc((1 + DimensionSaved) * sizeof(double)));
			assert(info.ZMax =
				(double *)malloc((1 + DimensionSaved) * sizeof(double)));
		}
		ComputeBounds(0, Dimension - 1, this);
		/*info.Contains.reset(new ContainsFunction(CoordType == THREED_COORDS ? Contains3D : Contains2D));
		info.BoxOverlaps.reset(new BoxOverlapsFunction(CoordType == THREED_COORDS ? BoxOverlaps3D : BoxOverlaps2D));
		info.CandidateSet.reset(new Candidate[K + 1]);*/
		info.Contains = CoordType == THREED_COORDS ? Contains3D : Contains2D;
		info.BoxOverlaps =
			CoordType == THREED_COORDS ? BoxOverlaps3D : BoxOverlaps2D;
		assert(info.CandidateSet =
			(Candidate *)malloc((K + 1) * sizeof(Candidate)));

		From = FirstNode;
		do {
			NearestQuadrantNeighbors(From, 0, K, this);
			for (i = 0; i < info.Candidates; i++) {
				To = info.CandidateSet[i].To;
				AddCandidate(From, To, (this->*D)(From, To), 1);
			}
		} while ((From = From->Suc) != FirstNode);

		/*info.CandidateSet.reset();
		info.KDTree.reset();
		info.XMin.reset();
		info.XMax.reset();
		info.YMin.reset();
		info.YMax.reset();*/
		free(info.CandidateSet);
		free(info.KDTree);
		free(info.XMin);
		free(info.XMax);
		free(info.YMin);
		free(info.YMax);
		if (CoordType == THREED_COORDS) {
			/*info.ZMin.reset();
			info.ZMax.reset();*/
			free(info.ZMin);
			free(info.ZMax);
		}
		if (info.Level == 0 && (WeightType == GEOM || WeightType == GEOM_MEEUS)) {
			Candidate **SavedCandidateSet;
			assert(SavedCandidateSet =
				(Candidate **)malloc((1 + DimensionSaved) *
					sizeof(Candidate *)));
			if (TraceLevel >= 2)
				printff("done\n");
			/* Transform longitude (180 and -180 map to 0) */
			From = FirstNode;
			do {
				SavedCandidateSet[From->Id] = From->CandidateSet;
				From->CandidateSet = 0;
				From->Yc = From->Y;
				From->Y += From->Y > 0 ? -180 : 180;
			} while ((From = From->Suc) != FirstNode);
			info.Level++;
			CreateNearestNeighborCandidateSet(K);
			info.Level--;
			From = FirstNode;
			do
				From->Y = From->Yc;
			while ((From = From->Suc) != FirstNode);
			do {
				Candidate *QCandidateSet = From->CandidateSet;
				Candidate *NFrom;
				From->CandidateSet = SavedCandidateSet[From->Id];
				for (NFrom = QCandidateSet; (To = NFrom->To); NFrom++)
					AddCandidate(From, To, NFrom->Cost, NFrom->Alpha);
				free(QCandidateSet);
			} while ((From = From->Suc) != FirstNode);
			free(SavedCandidateSet);
		}
		if (info.Level == 0) {
			ResetCandidateSet();
			AddTourCandidates();
			if (CandidateSetSymmetric)
				SymmetrizeCandidateSet();
			if (TraceLevel >= 2)
				printff("done\n");
		}
	}

	/*
	 * The ComputeBounds function computes the bounding boxes
	 * for the K-d tree nodes.
	 */

	static void ComputeBounds(int start, int end, LKHAlg*Alg)
	{
		auto& info = Alg->m_quaInfo;
		if (start <= end) {
			int mid = (start + end) / 2, i;
			LKHAlg::Node *T = info.KDTree[mid];
			info.XMin[T->Id] = info.YMin[T->Id] = std::numeric_limits<double>::max();
			info.XMax[T->Id] = info.YMax[T->Id] = -std::numeric_limits<double>::max();
			if (Alg->CoordType == LKH::THREED_COORDS) {
				info.ZMin[T->Id] = std::numeric_limits<double>::max();
				info.ZMax[T->Id] = -std::numeric_limits<double>::max();
			}
			for (i = start; i <= end; i++) {
				LKHAlg::Node *N = info.KDTree[i];
				if (N == T)
					continue;
				if (N->X < info.XMin[T->Id])
					info.XMin[T->Id] = N->X;
				if (N->X > info.XMax[T->Id])
					info.XMax[T->Id] = N->X;
				if (N->Y < info.YMin[T->Id])
					info.YMin[T->Id] = N->Y;
				if (N->Y > info.YMax[T->Id])
					info.YMax[T->Id] = N->Y;
				if (Alg->CoordType == LKH::THREED_COORDS) {
					if (N->Z < info.ZMin[T->Id])
						info.ZMin[T->Id] = N->Z;
					if (N->Z > info.ZMax[T->Id])
						info.ZMax[T->Id] = N->Z;
				}
			}
			ComputeBounds(start, mid - 1, Alg);
			ComputeBounds(mid + 1, end, Alg);
		}
	}


#define Coord(N, axis) (axis == 0 ? (N)->X : axis == 1 ? (N)->Y : (N)->Z)

	/*
	 * The Contains2D function returns 1 if T belongs to 2-D quadrant Q
	 * relative to N; otherwise 0.
	 *
	 *          Q = 2 | Q = 1
	 *          ===== N =====
	 *          Q = 3 | Q = 4
	 *
	 */

	static int Contains2D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg)
	{
		switch (Q) {
		case 1:
			return T->X >= N->X && T->Y >= N->Y;
		case 2:
			return T->X <= N->X && T->Y >= N->Y;
		case 3:
			return T->X <= N->X && T->Y <= N->Y;
		case 4:
			return T->X >= N->X && T->Y <= N->Y;
		default:
			return 1;
		}
	}

	/*
	 * The Contains3D function returns 1 if T belongs to 3-D
	 * quadrant Q relative to N; otherwise 0.
	 *
	 *          Q = 2 | Q = 1
	 *          ===== N =====   for T.Z >= N.Z
	 *          Q = 3 | Q = 4
	 *
	 *          Q = 6 | Q = 5
	 *          ===== N =====   for T.Z <= N.Z
	 *          Q = 7 | Q = 8
	 */

	static int Contains3D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg)
	{
		switch (Q) {
		case 1:
			return T->X >= N->X && T->Y >= N->Y && T->Z >= N->Z;
		case 2:
			return T->X <= N->X && T->Y >= N->Y && T->Z >= N->Z;
		case 3:
			return T->X <= N->X && T->Y <= N->Y && T->Z >= N->Z;
		case 4:
			return T->X >= N->X && T->Y <= N->Y && T->Z >= N->Z;
		case 5:
			return T->X >= N->X && T->Y >= N->Y && T->Z <= N->Z;
		case 6:
			return T->X <= N->X && T->Y >= N->Y && T->Z <= N->Z;
		case 7:
			return T->X <= N->X && T->Y <= N->Y && T->Z <= N->Z;
		case 8:
			return T->X >= N->X && T->Y <= N->Y && T->Z <= N->Z;
		default:
			return 1;
		}
	}

	/*
	 * The BoxOverlaps2D function returns 1 if T's bounding box
	 * overlaps the 2-D quadrant Q relative to N; otherwise 0.
	 */

	static int BoxOverlaps2D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg)
	{
		auto& info = Alg->m_quaInfo;
		switch (Q) {
		case 1:
			return info.XMax[T->Id] >= N->X && info.YMax[T->Id] >= N->Y;
		case 2:
			return info.XMin[T->Id] <= N->X && info.YMax[T->Id] >= N->Y;
		case 3:
			return info.XMin[T->Id] <= N->X && info.YMin[T->Id] <= N->Y;
		case 4:
			return info.XMax[T->Id] >= N->X && info.YMin[T->Id] <= N->Y;
		default:
			return 1;
		}
	}

	/*
	 * The BoxOverlaps3D function returns 1 if T's bounding box
	 * overlaps the 3-D quadrant Q relative to N; otherwise 0.
	 */

	static int BoxOverlaps3D(LKHAlg::Node * T, int Q, LKHAlg::Node * N, LKHAlg* Alg)
	{
		auto& info = Alg->m_quaInfo;
		switch (Q) {
		case 1:
			return info.XMax[T->Id] >= N->X && info.YMax[T->Id] >= N->Y &&
				info.ZMax[T->Id] >= N->Z;
		case 2:
			return info.XMin[T->Id] <= N->X && info.YMax[T->Id] >= N->Y &&
				info.ZMax[T->Id] >= N->Z;
		case 3:
			return info.XMin[T->Id] <= N->X && info.YMin[T->Id] <= N->Y &&
				info.ZMax[T->Id] >= N->Z;
		case 4:
			return info.XMax[T->Id] >= N->X && info.YMin[T->Id] <= N->Y &&
				info.ZMax[T->Id] >= N->Z;
		case 5:
			return info.XMax[T->Id] >= N->X && info.YMax[T->Id] >= N->Y &&
				info.ZMin[T->Id] <= N->Z;
		case 6:
			return info.XMin[T->Id] <= N->X && info.YMax[T->Id] >= N->Y &&
				info.ZMin[T->Id] <= N->Z;
		case 7:
			return info.XMin[T->Id] <= N->X && info.YMin[T->Id] <= N->Y &&
				info.ZMin[T->Id] <= N->Z;
		case 8:
			return info.XMax[T->Id] >= N->X && info.YMin[T->Id] <= N->Y &&
				info.ZMin[T->Id] <= N->Z;
		default:
			return 1;
		}
	}

	/*
	 * The Overlaps function returns 1 if High is zero and the half
	 * plane to the left of a point T overlaps quadrant Q relative
	 * to a point N.
	 * If High is not zero the Overlaps function returns 1 if the
	 * half plane to the right of T overlaps quadrant Q relative to N.
	 * Otherwise the function returns 0.
	 *
	 * The directions are relative to the given coordinate axis.
	 * The parameter diff is <= 0 if T is to the left of N,
	 * and >= 0 if T is to the right of N.
	 */

	static int Overlaps(int Q, double diff, int High, int axis)
	{
		switch (Q) {
		case 1:
			return High || diff >= 0;
		case 2:
			return axis == 0 ? !High || diff <= 0 : High || diff >= 0;
		case 3:
			return axis <= 1 ? !High || diff <= 0 : High || diff >= 0;
		case 4:
			return axis == 1 ? !High || diff <= 0 : High || diff >= 0;
		case 5:
			return axis <= 1 ? High || diff >= 0 : !High || diff <= 0;
		case 6:
			return axis == 1 ? High || diff >= 0 : !High || diff <= 0;
		case 7:
			return !High || diff <= 0;
		case 8:
			return axis == 0 ? High || diff >= 0 : !High || diff <= 0;
		default:
			return 1;
		}
	}

	static int InCandidateSet(LKHAlg::Node * N, LKHAlg::Node * T, LKHAlg *Alg)
	{

		auto& info = Alg->m_quaInfo;
		int i;
		for (i = 0; i < info.Candidates; i++)
			if (info.CandidateSet[i].To == T)
				return 1;
		return Alg->IsCandidate(N, T);
	}

	/*
	 * The NQN function searches info.KDTree[start:end] in an attempt to
	 * find the K quad-nearest neighbors in quadrant Q relative to N.
	 *
	 * The function is called from the NearestQuadrantNeighbors
	 * and the CreateNearestNeighborCandidateSet function.
	 */

	static void NQN(LKHAlg::Node * N, int Q, int start, int end, int K, LKHAlg *Alg)
	{

		auto& info = Alg->m_quaInfo;
		int mid = (start + end) / 2, d;
		LKHAlg::Node *T = info.KDTree[mid], P = { 0 };
		int axis = T->Axis;

		if (K == 0 || N->FixedTo2)
			return;
		if (start <= end && T != N && !T->FixedTo2 &&
			Alg->IsPossibleCandidate(N, T) &&
			info.Contains(T, Q, N, Alg) &&
			!InCandidateSet(N, T, Alg) &&
			(!Alg->c || (Alg->*(Alg->c))(N, T) - N->Pi - T->Pi <= info.Radius) &&
			(d = (Alg->*(Alg->Distance))(N, T) * Alg->Precision) <= info.Radius) {
			int i = info.Candidates;
			while (--i >= 0 && d < info.CandidateSet[i].Cost)
				info.CandidateSet[i + 1] = info.CandidateSet[i];
			info.CandidateSet[i + 1].To = T;
			info.CandidateSet[i + 1].Cost = d;
			if (info.Candidates < K)
				(info.Candidates)++;
			if (info.Candidates == K)
				info.Radius = info.CandidateSet[info.Candidates - 1].Cost;
		}
		if (start < end && info.BoxOverlaps(T, Q, N, Alg)) {
			double diff = Coord(T, axis) - Coord(N, axis);
			P.X = axis == 0 ? T->X : N->X;
			P.Y = axis == 1 ? T->Y : N->Y;
			P.Z = axis == 2 ? T->Z : N->Z;
			P.Pi = 0;
			if (diff >= 0) {
				if (Overlaps(Q, diff, 0, axis))
					NQN(N, Q, start, mid - 1, K, Alg);
				if (Overlaps(Q, diff, 1, axis) &&
					(!Alg->c || (Alg->*(Alg->c))(N, &P) - N->Pi <= info.Radius) &&
					(Alg->*(Alg->Distance))(N, &P) * Alg->Precision <= info.Radius)
					NQN(N, Q, mid + 1, end, K, Alg);
			}
			else {
				if (Overlaps(Q, diff, 1, axis))
					NQN(N, Q, mid + 1, end, K, Alg);
				if (Overlaps(Q, diff, 0, axis) &&
					(!Alg->c || (Alg->*(Alg->c))(N, &P) - N->Pi <= info.Radius) &&
					(Alg->*(Alg->Distance))(N, &P) * Alg->Precision <= info.Radius)
					NQN(N, Q, start, mid - 1, K, Alg);
			}
		}
	}

	/*
	 * The NearestQuadrantNeighbors function searches the K-d tree
	 * in an attempt to find the K quad-nearest neighbors in
	 * quadrant Q relative to N.
	 * If Q = 0, the funtion computes the K nearest neighbors to N.
	 */

	static void NearestQuadrantNeighbors(LKHAlg::Node * N, int Q, int K, LKHAlg *Alg)
	{
		auto& info = Alg->m_quaInfo;
		info.Candidates = 0;
		info.Radius = std::numeric_limits<int>::max();
		NQN(N, Q, 0, Alg->Dimension - 1, K, Alg);
	}
}