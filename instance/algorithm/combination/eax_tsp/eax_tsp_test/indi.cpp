#include "indi.h"
using namespace ofec::eax_tsp2;

TIndi::TIndi()
{
    fN = 0;
    fLink.clear();
   // fLink = NULL;
    fEvaluationValue = 0;
}
//
//TIndi::~TIndi()
//{
//    //for ( int i = 0; i < fN; ++i ) delete[] fLink[ i ];
//    //delete[] fLink;
//}

void TIndi::define( int N )
{
    fN = N;
    fLink.resize(N);
    for (auto& it : fLink) it.resize(2);
    //fLink = new int* [ fN ];
    //for( int i = 0; i < fN; ++i ) fLink[ i ] = new int [ 2 ];
    variable().resize(N);
    
}

void TIndi::transferSol(std::vector<int>& sol) const{
    sol.resize(fN);
    //std::vector<int> Array(Ncity);
    auto& Array(sol);
    int curr = 0, st = 0, count = 0, pre = -1, next;
    while (1)
    {
        Array[count++] = curr;
        if (count > fN)
        {
            //  printf("Invalid\n");
            return;
        }
        if (fLink[curr][0] == pre)
            next = fLink[curr][1];
        else
            next = fLink[curr][0];

        pre = curr;
        curr = next;
        if (curr == st)
            break;
    }
}

void ofec::eax_tsp2::TIndi::toCurSol(const std::vector<int>& sol)
{
    auto& fGene = sol;
    //for (auto& it : fGene) {
    //    --it;
    //}
    for (int j2 = 1; j2 < fN - 1; ++j2) {
        fLink[fGene[j2]][0] = fGene[j2 - 1];
        fLink[fGene[j2]][1] = fGene[j2 + 1];
    }
    fLink[fGene[0]][0] = fGene[fN - 1];
    fLink[fGene[0]][1] = fGene[1];
    fLink[fGene[fN - 1]][0] = fGene[fN - 2];
    fLink[fGene[fN - 1]][1] = fGene[0];
}

TIndi& TIndi::operator = ( const TIndi& src )
{
    //if (*this != src) {
    //       
    //}
    //else {
    //    *this = src;
    //    return *this;
    //}
    fN = src.fN;
    fLink = src.fLink;
    //for ( int i = 0; i < fN; ++i )
    //    for ( int j = 0; j < 2; ++j ) fLink[i][j] = src.fLink[i][j];
    fEvaluationValue = src.fEvaluationValue;
    return *this;
}

double TIndi::distanceTo(const TIndi& indi2)const {
    double dis = 0;
    for (int idx(0); idx < fN; ++idx) {
        auto& v1= fLink[idx];
        auto& v2 = indi2.fLink[idx];
        for (auto& it1 : v1) {
            for (auto& it2 : v2) {
                if (it1 == it2) {
                    ++dis;
                    break;
                }
            }
        }
    }
    //auto tmpdis = (fN * 2 - dis) / 2.0;
    //if (tmpdis < 0) {
    //    int stop = -1;
    //}
    return (fN * 2 - dis) / 2.0;
}

bool TIndi::operator == ( const TIndi& src ){



    int curr, next, pre, flag_identify;

    if( fN != src.fN ) return false;
    if( fEvaluationValue != src.fEvaluationValue ) return false;

    curr = 0;
    pre = -1;
    for( int i = 0; i < fN; ++i ){
        if( fLink[curr][0] == pre ) next = fLink[curr][1];
        else next = fLink[curr][0];

        if( src.fLink[curr][0] != next && src.fLink[curr][1] != next ) return false;
        pre = curr;
        curr = next;
    }
    return true;
}

