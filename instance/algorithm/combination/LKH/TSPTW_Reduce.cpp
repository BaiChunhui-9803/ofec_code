#include "./INCLUDE/LKH.h"

void LKH::LKHAlg::TSPTW_Reduce()
{
    int i, j, k, n = DimensionSaved;
    char **R;
    const int M = std::numeric_limits<int>::max() / 2 / Precision;

    if (Salesmen > 1)
        return;
    assert(R = (char **) malloc(((n + 1) * sizeof(char *))));
    for (i = 1; i <= n; i++)
        assert(R[i] = (char *) calloc(n + 1, sizeof(char)));
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            if (j != i &&
                NodeSet[j].Earliest + NodeSet[j].C[i] > NodeSet[i].Latest)
                R[i][j] = 1;
    /* Compute the transitive closure */
    for (k = 1; k <= n; k++)
        for (i = 1; i <= n; i++)
            if (R[i][k])
                for (j = 1; j <= n; j++)
                    R[i][j] |= R[k][j];
    /* Eliminate edges */
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            if (j != i && R[i][j])
                NodeSet[j].C[i] = M;
    /* Generate constraints */
    for (i = 1; i <= n; i++) {
        Node *Ni = &NodeSet[i];
        for (j = 1; j <= n; j++) {
            if (i != j && R[i][j]) {
                Node *Nj = &NodeSet[j];
                Constraint *Con;
                assert(Con = (Constraint *) malloc(sizeof(Constraint)));
                Con->t1 = Ni;
                Con->t2 = Nj;
                Con->Suc = FirstConstraint;
                FirstConstraint = Con;
                Con->Next = Ni->FirstConstraint;
                Ni->FirstConstraint = Con;
            }
        }
    }
    for (i = 1; i <= n; i++)
        free(R[i]);
    free(R);
}
