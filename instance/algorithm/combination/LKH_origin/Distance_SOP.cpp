#include "./INCLUDE/LKH.h"

/*
 * The Distance_SOP function computes the distance for a SOP instance.
 */

int LKH::LKHAlg::Distance_SOP(Node * Na, Node * Nb)
{
    int d = (this->*(this->OldDistance))(Na, Nb);
    return d >= 0 ? d : M;
}
