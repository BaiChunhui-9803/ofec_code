#include "evaluator.h"
#include<vector>
#include<fstream>
#include<string>
#include "../../../../../instance/problem/combination/travelling_salesman/travelling_salesman.h"

using namespace ofec::eax_tsp_mt;

TEvaluator::TEvaluator()
{
    Ncity = 0;
    fNearNumMax = 50;
}

TEvaluator::~TEvaluator(){
}


void TEvaluator::setInstance(ofec::Environment *env) {
    using namespace ofec;
    Ncity = env->problem()->numberVariables();


    //x = new double[Ncity];
    //y = new double[Ncity];
    //int* checkedN = new int[Ncity];
    int n(0);
    std::vector<int> checkedN(Ncity, 0);
    //  fclose(fp);


    //fEdgeDis = new int* [Ncity];
    //for (int i = 0; i < Ncity; ++i) fEdgeDis[i] = new int[Ncity];
    fEdgeDis.resize(Ncity);
    for (auto& it : fEdgeDis) it.resize(Ncity);

    //fNearCity = new int* [Ncity];
    //for (int i = 0; i < Ncity; ++i) fNearCity[i] = new int[fNearNumMax + 1];
    auto& cost(CAST_TSP(env->problem())->cost());
    for (int idx(0); idx < Ncity;++idx) {
        for (int idy(0); idy < Ncity; ++idy) {
            fEdgeDis[idx][idy] = cost[idx][idy];
        }
    }

    fNearCity.resize(Ncity);
    for (auto& it : fNearCity) it.resize(fNearNumMax + 1);

    int ci, j1, j2, j3;
    int cityNum = 0;
    double minDis;
    for (ci = 0; ci < Ncity; ++ci) {
        for (j3 = 0; j3 < Ncity; ++j3) checkedN[j3] = 0;
        checkedN[ci] = 1;
        fNearCity[ci][0] = ci;
        for (j1 = 1; j1 <= fNearNumMax; ++j1) {
            minDis = std::numeric_limits<double>::max();
            for (j2 = 0; j2 < Ncity; ++j2) {
                if (fEdgeDis[ci][j2] <= minDis && checkedN[j2] == 0) {
                    cityNum = j2;
                    minDis = fEdgeDis[ci][j2];
                }
            }
            fNearCity[ci][j1] = cityNum;
            checkedN[cityNum] = 1;
        }
    }


}
void TEvaluator::doIt( TIndi& indi )const{
    double d = 0;
    for( int i = 0; i < Ncity; ++i ) d += fEdgeDis[ i ][ indi.fLink[i][0] ] + fEdgeDis[ i ][ indi.fLink[i][1] ];
    indi.fEvaluationValue = d/2;
}


void TEvaluator::transferSol(const TIndi& indi, std::vector<int>& sol)const {
    sol.resize(Ncity);
    //std::vector<int> Array(Ncity);
    auto& Array(sol);
    int curr = 0, st = 0, count = 0, pre = -1, next;
    while (1)
    {
        Array[count++] = curr + 1;
        if (count > Ncity)
        {
            //  printf("Invalid\n");
            return;
        }
        if (indi.fLink[curr][0] == pre)
            next = indi.fLink[curr][1];
        else
            next = indi.fLink[curr][0];

        pre = curr;
        curr = next;
        if (curr == st)
            break;
    }
}


bool TEvaluator::checkValid( int* array, double value )const {
    //int *check=new int[Ncity];
    std::vector<int> check(Ncity);
    for( int i = 0; i < Ncity; ++i ) check[ i ] = 0;
    for( int i = 0; i < Ncity; ++i ) ++check[ array[ i ]-1 ];
    for( int i = 0; i < Ncity; ++i )
        if( check[ i ] != 1 ) return false;
    int distance = 0;
    for( int i = 0; i < Ncity-1; ++i )
        distance += fEdgeDis[ array[ i ]-1 ][ array[ i+1 ]-1 ];

    distance += fEdgeDis[ array[ Ncity-1 ]-1 ][ array[ 0 ]-1 ];

 //   delete [] check;
    if( distance != value ) return false;
    return true;
}

