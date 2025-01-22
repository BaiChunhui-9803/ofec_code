#include "./INCLUDE/LKH.h"

void LKH::LKHAlg::VRPB_Reduce()
{
    int i, j, n = Dim;
    const int M = std::numeric_limits<int>::max() / 2 / Precision;

    for (i = 1; i <= n; i++) {
        if (NodeSet[i].Backhaul) {
            for (j = 1; j <= n; j++)
                if (j != i && j != MTSPDepot && !NodeSet[j].Backhaul)
                    NodeSet[i].C[j] = M;
        }
    }
}
