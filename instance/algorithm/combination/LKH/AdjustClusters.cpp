#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Heap.h"

/* 
 * The AdjustClusters function adjusts K clusters with centers Center
 * so that they all have a size less than SubproblemSize.
 *
 * The function is called from KMeansClustering and KCenterClustering.
 */

void LKH::LKHAlg::AdjustClusters(int K, Node ** Center)
{
	int d, i, j;
	Node *N = 0;
	int *Size = 0;

    assert(Size = (int *) calloc((K + 1), sizeof(int)));
    N = FirstNode;
    do
        Size[N->Subproblem]++;
    while ((N = N->Suc) != FirstNode);

    /* Remove nodes from clusters with size > SubproblemSize */
    for (i = 1; i <= K; i++) {
        if (Size[i] > SubproblemSize) {
            N = FirstNode;
            do {
                if (N->Subproblem == i) {
                    N->Rank = -N->Cost;
                    HeapLazyInsert(N,this);
                }
            } while ((N = N->Suc) != FirstNode);
            Heapify(this);
            for (j = 1; j <= SubproblemSize; j++)
                HeapDeleteMin(this);
            while ((N = HeapDeleteMin(this)))
                N->Subproblem = 0;
            Size[i] = SubproblemSize;
        }
    }
    /* Insert removed nodes into cluster with size < SubproblemSize */
    N = FirstNode;
    do {
        if (N->Subproblem == 0) {
            N->Cost = std::numeric_limits<int>::max();
            j = 0;
            for (i = 1; i <= K; i++) {
				if (Size[i] < SubproblemSize &&
					(d = (this->*Distance)(N, Center[i])) < N->Cost) {
                    N->Cost = d;
                    j = i;
                }
            }
            Size[N->Subproblem = j]++;
        }
    } while ((N = N->Suc) != FirstNode);
    free(Size);
}
