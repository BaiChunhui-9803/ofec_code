#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Delaunay.h"

/*
 * The CreateDelaunayCandidateSet function determines for each node its set 
 * of incident candidate edges. The edges are found by Delaunay triangulation. 
 *
 * The function is called from CreateCandidateSet.
 */
namespace LKH {

	//static thread_local unique_ptr<int> m_DelaunayLevel;

	void LKHAlg::freeCreateDelaunay()
	{
		//m_DelaunayLevel.reset();
		m_DelaunayLevel = 0;
	}
	void LKHAlg::CreateDelaunayCandidateSet()
	{
		Node *From, *To;
		point *u, *v;
		edge *e_start, *e;
		int d, i, Count;

		/*if (!m_DelaunayLevel.get())
			m_DelaunayLevel.reset(new int(0));*/
		if (TraceLevel >= 2)
			printff("Creating Delaunay candidate set ... ");
		if (m_DelaunayLevel == 0 && MaxCandidates == 0) {
			AddTourCandidates();
			From = FirstNode;
			do {
				if (!From->CandidateSet)
					eprintf("MAX_CANDIDATES = 0: No candidates");
			} while ((From = From->Suc) != FirstNode);
			if (TraceLevel >= 2)
				printff("done\n");
			return;
		}

		/* Find the Delaunay edges */
		delaunay(Dimension, this);

		/* Add the Delaunay edges to the candidate set */
		for (i = 0; i < Dimension; i++) {
			u = &p_array[i];
			From = &NodeSet[u->id];
			e_start = e = u->entry_pt;
			Count = 0;
			do {
				v = Other_point(e, u);
				if (u < v) {
					To = &NodeSet[v->id];
					d = (this->*D)(From, To);
					AddCandidate(From, To, d, 1);
					AddCandidate(To, From, d, 1);
				}
			} while ((e = Next(e, u)) != e_start && ++Count < Dimension);
		}
		free_memory();
		if (m_DelaunayLevel == 0 &&
			(WeightType == GEO || WeightType == GEOM ||
				WeightType == GEO_MEEUS || WeightType == GEOM_MEEUS)) {
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
				m_DelaunayLevel++;
				CreateDelaunayCandidateSet();
				m_DelaunayLevel--;
				From = FirstNode;
				do
					From->Y = From->Zc;
				while ((From = From->Suc) != FirstNode);
			}
		}
		if (m_DelaunayLevel == 0) {
			AddTourCandidates();
			/* Add quadrant neighbors if any node has less than two candidates.
			   That is, if it should happen that delaunay_edges fails. */
			if (ExtraCandidates < 2) {
				/* Add quadrant neighbors if any node has less than two candidates.
				   That is, if it should happen that delaunay_edges fails. */
				From = FirstNode;
				do {
					if (From->CandidateSet == 0 ||
						From->CandidateSet[0].To == 0 ||
						From->CandidateSet[1].To == 0) {
						if (TraceLevel >= 2)
							printff("*** Not complete ***\n");
						AddExtraCandidates(CoordType == THREED_COORDS ? 8 : 4,
							QUADRANT, 1);
						break;
					}
				} while ((From = From->Suc) != FirstNode);
			}
			if (TraceLevel >= 2)
				printff("done\n");
		}
	}
}