#include "evaluator.h"

#include<vector>
#include<fstream>
#include<string>

#include "../../../../../instance/problem/combination/travelling_salesman/travelling_salesman.h"

using namespace ofec::eax_tsp;

TEvaluator::TEvaluator()
{
  //  fEdgeDis = NULL;
  //  fNearCity = NULL;
    Ncity = 0;
    fNearNumMax = 50;
}

TEvaluator::~TEvaluator(){
    //for ( int i = 0; i < Ncity; ++i ) delete[] fEdgeDis[ i ];
    //delete[] fEdgeDis;

    //for ( int i = 0; i < Ncity; ++i ) delete[] fNearCity[ i ];
    //delete[] fNearCity;

  //  delete [] x;
  //  delete [] y;
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

    //delete[] checkedN;
    //checkedN = nullptr;

}
//
//void TEvaluator::setInstance(const std::string& filename)
//{
//
//    std::fstream in(filename);
//    if(!in.is_open())     exit(0);
//
//    {
//        std::string word, type;
//        while (in >> word) {
//            if (word == "DIMENSION") {
//                in >> word >> Ncity;
//            }
//            if (word == "EDGE_WEIGHT_TYPE") {
//                in >> word >> type;
//            }
//            if (word == "NODE_COORD_SECTION") break;
//        }
//
//        if (word != "NODE_COORD_SECTION") {
//            printf("Error in reading the instance\n");
//            exit(0);
//        }
//
//        x = new double[Ncity];
//        y = new double[Ncity];
//        int* checkedN = new int[Ncity];
//        int n(0);
//        for (int i = 0; i < Ncity; ++i) {
//            in >> n >> word;
//            x[i] = std::stod(word);
//            in >> word;
//            y[i] = std::stod(word);
//            //fscanf(fp, "%d", &n);
//            //fscanf(fp, "%s", word);
//            //x[i] = atof(word);
//            //fscanf(fp, "%s", word);
//            //y[i] = atof(word);
//        }
//      //  fclose(fp);
//        fEdgeDis = new int* [Ncity];
//        for (int i = 0; i < Ncity; ++i) fEdgeDis[i] = new int[Ncity];
//        fNearCity = new int* [Ncity];
//        for (int i = 0; i < Ncity; ++i) fNearCity[i] = new int[fNearNumMax + 1];
//
//        if (type =="EUC_2D") {
//            for (int i = 0; i < Ncity; ++i)
//                for (int j = 0; j < Ncity; ++j)
//                    fEdgeDis[i][j] = (int)(sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j])) + 0.5);
//        }
//        else if (type == "ATT") {
//            for (int i = 0; i < Ncity; ++i) {
//                for (int j = 0; j < Ncity; ++j) {
//                    double r = (sqrt(((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j])) / 10.0));
//                    int t = (int)r;
//                    if ((double)t < r) fEdgeDis[i][j] = t + 1;
//                    else fEdgeDis[i][j] = t;
//                }
//            }
//        }
//        else if (type == "CEIL_2D" ) {
//            for (int i = 0; i < Ncity; ++i)
//                for (int j = 0; j < Ncity; ++j)
//                    fEdgeDis[i][j] = (int)ceil(sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j])));
//        }
//        else {
//            printf("EDGE_WEIGHT_TYPE is not supported\n");
//            exit(1);
//        }
//        int ci, j1, j2, j3;
//        int cityNum = 0;
//        int minDis;
//        for (ci = 0; ci < Ncity; ++ci) {
//            for (j3 = 0; j3 < Ncity; ++j3) checkedN[j3] = 0;
//            checkedN[ci] = 1;
//            fNearCity[ci][0] = ci;
//            for (j1 = 1; j1 <= fNearNumMax; ++j1) {
//                minDis = 100000000;
//                for (j2 = 0; j2 < Ncity; ++j2) {
//                    if (fEdgeDis[ci][j2] <= minDis && checkedN[j2] == 0) {
//                        cityNum = j2;
//                        minDis = fEdgeDis[ci][j2];
//                    }
//                }
//                fNearCity[ci][j1] = cityNum;
//                checkedN[cityNum] = 1;
//            }
//        }
//
//        delete[] checkedN;
//        checkedN = nullptr;
//    }
//
//
//    in.close();
//
//    //FILE* fp;
//    //int n;
//    //char word[ 80 ], type[ 80 ];
//    //
//
//
//    //fp = fopen( filename, "r" );
//
//    ///* read instance */
//    //while( 1 ){
//    //    if( fscanf( fp, "%s", word ) == EOF ) break;
//    //    if( strcmp( word, "DIMENSION" ) == 0 ){
//    //        fscanf( fp, "%s", word );
//    //        fscanf( fp, "%d", &Ncity );
//    //    }
//    //    if( strcmp( word, "EDGE_WEIGHT_TYPE" ) == 0 ){
//    //        fscanf( fp, "%s", word );
//    //        fscanf( fp, "%s", type );
//    //    }
//    //    if( strcmp( word, "NODE_COORD_SECTION" ) == 0 ) break;
//    //}
//    //if( strcmp( word, "NODE_COORD_SECTION" ) != 0 ){
//    //    printf( "Error in reading the instance\n" );
//    //    exit(0);
//    //}
//
//
//    //x = new double [ Ncity ];
//    //y = new double [ Ncity ];
//    //int *checkedN = new int[Ncity];
//
//    //for( int i = 0; i < Ncity; ++i ){
//    //    fscanf( fp, "%d", &n );
//    //    fscanf( fp, "%s", word );
//    //    x[ i ] = atof( word );
//    //    fscanf( fp, "%s", word );
//    //    y[ i ] = atof( word );
//    //}
//    //fclose(fp);
//    //fEdgeDis = new int* [ Ncity ];
//    //for( int i = 0; i < Ncity; ++i ) fEdgeDis[ i ] = new int [ Ncity ];
//    //fNearCity = new int* [ Ncity ];
//    //for( int i = 0; i < Ncity; ++i ) fNearCity[ i ] = new int [ fNearNumMax+1 ];
//
//    //if( strcmp( type, "EUC_2D" ) == 0  ) {
//    //    for( int i = 0; i < Ncity ; ++i )
//    //        for( int j = 0; j < Ncity ; ++j )
//    //            fEdgeDis[ i ][ j ]=(int)(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))+0.5);
//    //}
//    //else if( strcmp( type, "ATT" ) == 0  ) {
//    //    for( int i = 0; i < Ncity; ++i ){
//    //        for( int j = 0; j < Ncity; ++j ) {
//    //            double r = (sqrt(((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))/10.0));
//    //            int t = (int)r;
//    //            if( (double)t < r ) fEdgeDis[ i ][ j ] = t+1;
//    //            else fEdgeDis[ i ][ j ] = t;
//    //        }
//    //    }
//    //}
//    //else if( strcmp( type, "CEIL_2D" ) == 0  ){
//    //for( int i = 0; i < Ncity ; ++i )
//    //    for( int j = 0; j < Ncity ; ++j )
//    //        fEdgeDis[ i ][ j ]=(int)ceil(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
//    //}
//    //else{
//    //    printf( "EDGE_WEIGHT_TYPE is not supported\n" );
//    //    exit( 1 );
//    //}
//    //int ci, j1, j2, j3;
//    //int cityNum = 0;
//    //int minDis;
//    //for( ci = 0; ci < Ncity; ++ci ){
//    //    for( j3 = 0; j3 < Ncity; ++j3 ) checkedN[ j3 ] = 0;
//    //    checkedN[ ci ] = 1;
//    //    fNearCity[ ci ][ 0 ] = ci;
//    //    for( j1 = 1; j1 <= fNearNumMax; ++j1 ) {
//    //        minDis = 100000000;
//    //        for( j2 = 0; j2 < Ncity; ++j2 ){
//    //            if( fEdgeDis[ ci ][ j2 ] <= minDis && checkedN[ j2 ] == 0 ){
//    //                cityNum = j2;
//    //                minDis = fEdgeDis[ ci ][ j2 ];
//    //            }
//    //        }
//    //        fNearCity[ ci ][ j1 ] = cityNum;
//    //        checkedN[ cityNum ] = 1;
//    //    }
//    //}
//}

void TEvaluator::doIt( TIndi& indi ){
    double d = 0;
    for( int i = 0; i < Ncity; ++i ) d += fEdgeDis[ i ][ indi.fLink[i][0] ] + fEdgeDis[ i ][ indi.fLink[i][1] ];
    indi.fEvaluationValue = d/2;
}

//void TEvaluator::writeTo( FILE* fp, TIndi& indi )
//{
//    Array=new int[Ncity];
//    int curr=0, st=0, count=0, pre=-1, next;
//    while( 1 ){
//        Array[ count++ ] = curr + 1;
//        if( count > Ncity ){
//            printf( "Invalid\n" );
//            return;
//        }
//        if( indi.fLink[ curr ][ 0 ] == pre ) next = indi.fLink[ curr ][ 1 ];
//        else next = indi.fLink[ curr ][ 0 ];
//
//        pre = curr;
//        curr = next;
//        if( curr == st ) break;
//    }
//    if( this->checkValid( Array, indi.fEvaluationValue ) == false )
//        printf( "Individual is invalid \n" );
//
//    fprintf( fp, "%d %d\n", indi.fN, indi.fEvaluationValue );
//    for( int i = 0; i < indi.fN; ++i )
//        fprintf( fp, "%d ", Array[ i ] );
//    fprintf( fp, "\n" );
//}
//
//void TEvaluator::writeToStdout(TIndi &indi)
//{
//    Array = new int[Ncity];
//    int curr = 0, st = 0, count = 0, pre = -1, next;
//    while (1)
//    {
//        Array[count++] = curr + 1;
//        if (count > Ncity)
//        {
//            printf("Invalid\n");
//            return;
//        }
//        if (indi.fLink[curr][0] == pre)
//            next = indi.fLink[curr][1];
//        else
//            next = indi.fLink[curr][0];
//
//        pre = curr;
//        curr = next;
//        if (curr == st)
//            break;
//    }
//    if (this->checkValid(Array, indi.fEvaluationValue) == false)
//        printf("Individual is invalid \n");
//
//    printf("%d %d\n", indi.fN, indi.fEvaluationValue);
//    for (int i = 0; i < indi.fN; ++i)
//        printf("%d ", Array[i]);
//    printf("\n");
//}
//
//void TEvaluator::writeTo(std::ofstream& out, const TIndi& indi)
//{
//
//    std::vector<int> Array (Ncity);
//    int curr = 0, st = 0, count = 0, pre = -1, next;
//    while (1)
//    {
//        Array[count++] = curr + 1;
//        if (count > Ncity)
//        {
//          //  printf("Invalid\n");
//            return;
//        }
//        if (indi.fLink[curr][0] == pre)
//            next = indi.fLink[curr][1];
//        else
//            next = indi.fLink[curr][0];
//
//        pre = curr;
//        curr = next;
//        if (curr == st)
//            break;
//    }
//
//
//
//    if (this->checkValid(Array.data(), indi.fEvaluationValue) == false) {
//        out << "invalid" << std::endl;
//    }
//    else {
//        out << "valid" << std::endl;
//    }
//   //     printf("Individual is invalid \n");
//    out << indi.fN << "\t" << indi.fEvaluationValue << std::endl;
//   // printf("%d %d\n", indi.fN, indi.fEvaluationValue);
//    for (int i = 0; i < indi.fN; ++i) {
//        out << Array[i] << "\t";
//    }
//    out << std::endl;
//    //    printf("%d ", Array[i]);
//    //printf("\n");
//}
void TEvaluator::transferSol(const TIndi& indi, std::vector<int>& sol) {
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


bool TEvaluator::checkValid( int* array, double value ){
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

