#include "cross.h"
#include "environment.h"
using namespace ofec::eax_tsp2;

TCross::TCross( int N){

//    Random *rnd = ofec::ADD_RND(0.5);

    tRand.reset(new TRandom());

    tSort.reset(new TSort());

    fMaxNumOfABcycle = 2000; /* Set an appropriate value (2000 is usually enough) */
    fN = N;
    tBestTmp.define( fN );

    nearData.resize(fN);
    for (auto& it : nearData) it.resize(5);

    //nearData = new int* [ fN ];
    //for ( int j = 0; j < fN; ++j )
    //    nearData[j] = new int [ 5 ];

    
    fABcycle.resize(fMaxNumOfABcycle);
    for (auto& it : fABcycle) it.resize(2 * fN + 4);
    //fABcycle = new int* [ fMaxNumOfABcycle ];
    //for ( int j = 0; j < fMaxNumOfABcycle; ++j ) fABcycle[j] = new int [ 2*fN + 4 ];

    koritsu.resize(fN);
//    koritsu = new int [ fN ];
    bunki.resize(fN);
 //   bunki = new int [ fN ];
    koriInv.resize(fN);
  //  koriInv = new int [ fN ];
    bunInv.resize(fN);
  //  bunInv = new int [ fN ];
    checkKoritsu.resize(fN);
   // checkKoritsu = new int [ fN ];
    fRoute.resize(2 * fN + 1);
  //  fRoute = new int [ 2*fN + 1 ];
    fPermu.resize(fMaxNumOfABcycle);
//    fPermu = new int [ fMaxNumOfABcycle ];

    fC.resize(2 * fN + 4);
  //  fC = new int [ 2*fN+4 ];
    fJun.resize(fN + 1);
 //   fJun = new int[ fN+ 1 ];
    fOrd1.resize(fN );
 //   fOrd1 = new int [ fN ];
    fOrd2.resize(fN);
 //   fOrd2 = new int [ fN ];

    // Speed Up Start
    fOrder.resize(fN);
   // fOrder = new int [ fN ];
    fInv.resize(fN);
  //  fInv = new int [ fN ];
    fSegment.resize(fN);
    for (auto& it : fSegment) it.resize(2);
 //   fSegment = new int* [ fN ];
 //   for ( int j = 0; j < fN; ++j ) fSegment[ j ] = new int [ 2 ];

    fSegUnit.resize(fN);
   // fSegUnit = new int [ fN ];
    fSegPosiList.resize(fN);
  //  fSegPosiList = new int[ fN ];

    LinkAPosi.resize(fN);
    //LinkAPosi = new int [ fN ];

    LinkBPosi.resize(fN);
    for (auto& it : LinkBPosi) it.resize(2);
  //  LinkBPosi = new int* [ fN ];
  //  for ( int j = 0; j < fN; ++j ) LinkBPosi[ j ] = new int [ 2 ];

    fPosiSeg.resize(fN);
   // fPosiSeg = new int [ fN ];
    fNumOfElementInUnit.resize(fN);
   // fNumOfElementInUnit = new int [ fN ];

    fCenterUnit.resize(fN);
    std::fill(fCenterUnit.begin(), fCenterUnit.end(), 0);
  //  fCenterUnit = new int [ fN ];
//    for ( int j = 0; j < fN; ++j ) fCenterUnit[ j ] = 0;

    fListOfCenterUnit.resize(fN + 2);
//    fListOfCenterUnit = new int [ fN+2 ];
    fSegForCenter.resize(fN );
//    fSegForCenter = new int [ fN ];
    fGainAB.resize(fN);
  //  fGainAB = new int [ fN ];
    fModiEdge.resize(fN);
    for (auto& it : fModiEdge) it.resize(4);
  //  fModiEdge = new int* [ fN ];
 //   for ( int j = 0; j < fN; ++j ) fModiEdge[ j ] = new int [ 4 ];

    fBestModiEdge.resize(fN);
    for (auto& it : fBestModiEdge) it.resize(4);
 //   fBestModiEdge = new int* [ fN ];
 //   for ( int j = 0; j < fN; ++j ) fBestModiEdge[ j ] = new int [ 4 ];
    fAppliedCylce.resize(fN);
  //  fAppliedCylce = new int [ fN ];
    fBestAppliedCylce.resize(fN);
 //   fBestAppliedCylce = new int [ fN ];
    // Speed Up End

    // Block2
    fNumOfElementINAB.resize(fMaxNumOfABcycle);
   // fNumOfElementINAB = new int [ fMaxNumOfABcycle ];
    fInEffectNode.resize(fN);
    for (auto& it : fInEffectNode) it.resize(2);
  //  fInEffectNode = new int* [ fN ];
 //   for( int i = 0; i < fN; ++i ) fInEffectNode[ i ] = new int [ 2 ];

    fWeightRR.resize(fMaxNumOfABcycle);
    for (auto& it : fWeightRR) it.resize(fMaxNumOfABcycle);

 //   fWeightRR = new int* [ fMaxNumOfABcycle ];
 //   for( int i = 0; i < fMaxNumOfABcycle; ++i ) fWeightRR[ i ] = new int [ fMaxNumOfABcycle ];

    fWeightSR.resize(fMaxNumOfABcycle);
  //  fWeightSR = new int [ fMaxNumOfABcycle ];
    fWeightC.resize(fMaxNumOfABcycle);
   // fWeightC = new int [ fMaxNumOfABcycle ];
    fUsedAB.resize(fN);
  //  fUsedAB = new int [ fN ];
    fMovedAB.resize(fN);
  //  fMovedAB = new int [ fN ];
    fABcycleInEset.resize(fMaxNumOfABcycle);
 //   fABcycleInEset = new int [ fMaxNumOfABcycle ];


    
}

void TCross::bindEnvironment(TEnvironment* env) {
    m_tmp_environemt = env;
    tRand->setEnv(m_tmp_environemt);
    tSort->setEnv(m_tmp_environemt);
}
//
//TCross::~TCross()
//{
//
//
//    //delete [] koritsu;
//    //delete [] bunki;
//    //delete [] koriInv;
//    //delete [] bunInv;
//    //delete [] checkKoritsu;
//    //delete [] fRoute;
//    //delete [] fPermu;
//    //for ( int j = 0; j < fN; ++j ) delete[] nearData[ j ];
//
//    //delete[] nearData;
//    //for ( int j = 0; j < fMaxNumOfABcycle; ++j ) delete[] fABcycle[ j ];
//
//    //delete[] fABcycle;
//    //delete [] fC;
//    //delete [] fJun;
//    //delete [] fOrd1;
//    //delete [] fOrd2;
//
//    //// Speed Up Start
//    //delete [] fOrder;
//    //delete [] fInv;
//
//    //for ( int j = 0; j < fN; ++j ) delete[] fSegment[ j ];
//
//    //delete[] fSegment;
//    //delete[] fSegUnit;
//    //delete [] fSegPosiList;
//    //delete [] LinkAPosi;
//    //for ( int j = 0; j < fN; ++j ) delete[] LinkBPosi[ j ];
//
//    //delete [] LinkBPosi;
//    //delete [] fPosiSeg;
//    //delete [] fNumOfElementInUnit;
//    //delete [] fCenterUnit;
//    //delete [] fListOfCenterUnit;
//    //delete [] fSegForCenter;
//    //delete [] fGainAB;
//
//    //for ( int j = 0; j < fN; ++j ) delete[] fModiEdge[ j ];
//    //delete [] fModiEdge;
//    //for ( int j = 0; j < fN; ++j ) delete[] fBestModiEdge[ j ];
//    //delete [] fBestModiEdge;
//
//    //delete [] fAppliedCylce;
//    //delete [] fBestAppliedCylce;
//    //// Speed Up End
//
//    //// Block2
//    //delete [] fNumOfElementINAB;
//    //for ( int j = 0; j < fN; ++j ) delete [] fInEffectNode[ j ];
//
//    //delete [] fInEffectNode;
//    //for( int i = 0; i < fMaxNumOfABcycle; ++i ) delete [] fWeightRR[ i ];
//
//    //delete [] fWeightSR;
//    //delete [] fWeightC;
//    //delete [] fUsedAB;
//    //delete [] fMovedAB;
//    //delete [] fABcycleInEset;
//
//    //delete tRand;
//    //tRand = nullptr;
//    //delete tSort;
//    //tSort = nullptr;
//
//  //  ofec::DEL_RND(0.5);
//}

void TCross::setParents( const TIndi& tPa1, const TIndi& tPa2, const std::vector<int>& flagC, int numOfKids, Random* rnd)
{
    this->setABcycle( tPa1, tPa2, flagC, numOfKids , rnd);

    m_tmp_environemt->testRandom(rnd, "TCross::setParents end");
    fDisAB = 0;
    int curr=0, next=0, st=0, pre=0;
    st = 0;
    curr=-1;
    next = st;
    for( int i = 0; i < fN; ++i ){
        pre=curr;
        curr=next;
        if( tPa1.fLink[curr][0] != pre ) next = tPa1.fLink[ curr ][ 0 ];
        else next=tPa1.fLink[curr][1];

        if( tPa2.fLink[ curr ][ 0 ] != next && tPa2.fLink[ curr ][ 1 ] != next ) ++fDisAB;
        fOrder[ i ] = curr;
        fInv[ curr ] = i;
    }
    if (flagC[1] == 2) /* Block2 */
    {
        fTmax = 10;
        fMaxStag = 20; //1:Greedy LS, 20:Tabu Search
        this->setWeight( tPa1, tPa2 );
    }

}

void TCross::doIt( TIndi& tKid, TIndi& tPa2, int numOfKids, int flagP, const std::vector<int>& flagC, 
    std::vector<std::vector<int>>& fEdgeFreq, Random* rnd)
{

    m_tmp_environemt->testRandom(rnd, "before TCross doIt");
    int Num;
    int jnum, centerAB;
    double gain;
    double BestGain;
    double pointMax, point;
    double DLoss(0);

    fEvalType = flagC[ 0 ]; //1:Greedy, 2:---, 3:Distance, 4:Entropy
    fEsetType = flagC[ 1 ]; // 1:Single-AB, 2:Block2

    if ( numOfKids <= fNumOfABcycle ) Num = numOfKids;
    else Num = fNumOfABcycle;

    if (fEsetType == 1) /* Single-AB */
    {
        m_tmp_environemt->testRandom(rnd, "before tRand->permutation at doIt");
        tRand->permutation( fPermu, fNumOfABcycle, fNumOfABcycle , rnd);
        m_tmp_environemt->testRandom(rnd, "after tRand->permutation at doIt");
    }

    else if( fEsetType == 2 )
    {  // Block2
        for( int k =0; k< fNumOfABcycle; ++k ) fNumOfElementINAB[ k ] = fABcycle[ k ][ 0 ];
        tSort->indexB( fNumOfElementINAB, fNumOfABcycle, fPermu, fNumOfABcycle );
    }
    fNumOfGeneratedCh = 0;
    pointMax = 0.0;
    BestGain = 0;
    fFlagImp = 0;

///    m_localsearch_indis.clear();
    for( int j =0; j < Num; ++j )
    {
        fNumOfABcycleInEset = 0;
        if (fEsetType == 1) /* Single-AB */
        {
            jnum = fPermu[ j ];
            fABcycleInEset[ fNumOfABcycleInEset++ ] = jnum;
        }
        else if (fEsetType == 2) /* Block2 */
        {
            jnum = fPermu[ j ];
            centerAB = jnum;
            for( int s = 0; s < fNumOfABcycle; ++s ){
                if( s == centerAB ) fABcycleInEset[ fNumOfABcycleInEset++ ] = s;
                else{

                    m_tmp_environemt->testRandom(rnd, " fABcycleInEset 289 at doIt before");
                    if( fWeightRR[ centerAB ][ s ] > 0 && fABcycle[ s ][ 0 ] < fABcycle[ centerAB ][ 0 ] )
                    if(tRand->Double(0,1, rnd) <=0.5 ) fABcycleInEset[ fNumOfABcycleInEset++ ] = s;

                    m_tmp_environemt->testRandom(rnd, " fABcycleInEset 289 at doIt");
                }
            }

            m_tmp_environemt->testRandom(rnd, " searchEset 292 at doIt before");
        this->searchEset( centerAB,rnd );

        m_tmp_environemt->testRandom(rnd, " searchEset 292 at doIt");
        }
        //fNumOfSPL = 0;
        fSegPosiList.clear();
        gain = 0;
        fNumOfAppliedCycle = 0;
        fNumOfModiEdge = 0;

        fNumOfAppliedCycle = fNumOfABcycleInEset;
        for( int k = 0; k < fNumOfAppliedCycle; ++k ){
            fAppliedCylce[ k ] = fABcycleInEset[ k ];
            jnum = fAppliedCylce[ k ];
            this->changeSol( tKid, jnum, flagP );
            gain += fGainAB[ jnum ];
        }
        m_tmp_environemt->testRandom(rnd, " makeUnit at doIt before");
        this->makeUnit(rnd);
        m_tmp_environemt->testRandom(rnd, " makeUnit 310 at doIt");
        this->makeCompleteSol( tKid ,rnd);

        m_tmp_environemt->testRandom(rnd, " makeCompleteSol 312 at doIt");
        gain += fGainModi;

        ++fNumOfGeneratedCh;

        if (fEvalType == 1) /* Greedy */
            DLoss = 1.0;
        else if (fEvalType == 3) /* Distance preservation */
            DLoss = this->calAdpLoss( fEdgeFreq );
        else if( fEvalType == 4 )  /* Entropy preservation */
            DLoss = this->calEntLoss( fEdgeFreq );

        if( DLoss <= 0.0 ) DLoss = 0.00000001;

        point = (double)gain / DLoss;
        tKid.fEvaluationValue = tKid.fEvaluationValue - gain;
//        m_localsearch_indis.push_back(tKid);
        



        if( pointMax < point && (2 * fBestNumE < fDisAB || tKid.fEvaluationValue != tPa2.fEvaluationValue ) )
        {
            pointMax = point;
            BestGain = gain;
            fFlagImp = 1;

            fNumOfBestAppliedCycle = fNumOfAppliedCycle;
            for( int s = 0; s < fNumOfBestAppliedCycle; ++s ) fBestAppliedCylce[ s ] = fAppliedCylce[ s ];

            fNumOfBestModiEdge = fNumOfModiEdge;
            for( int s = 0; s < fNumOfBestModiEdge; ++s ){
                fBestModiEdge[ s ][ 0 ] = fModiEdge[ s ][ 0 ];
                fBestModiEdge[ s ][ 1 ] = fModiEdge[ s ][ 1 ];
                fBestModiEdge[ s ][ 2 ] = fModiEdge[ s ][ 2 ];
                fBestModiEdge[ s ][ 3 ] = fModiEdge[ s ][ 3 ];
            }

        }
        this->backToPa1( tKid );
        tKid.fEvaluationValue = tKid.fEvaluationValue + gain;
    }
    if( fFlagImp == 1 ){
        this->goToBest( tKid );
        tKid.fEvaluationValue = tKid.fEvaluationValue - BestGain;

      //  m_localsearch_indis.push_back(tKid);
        this->incrementEdgeFreq( fEdgeFreq );
    }
}

void TCross::doItWithoutParent(TIndi& tKid, int numOfKids, int flagP, const std::vector<int>& flagC,
    std::vector<std::vector<int>>& fEdgeFreq, Random* rnd) {
    int Num;
    int jnum, centerAB;
    double gain;
    double BestGain;
    double pointMax, point;
    double DLoss(0);

    fEvalType = flagC[0]; //1:Greedy, 2:---, 3:Distance, 4:Entropy
    fEsetType = flagC[1]; // 1:Single-AB, 2:Block2

    if (numOfKids <= fNumOfABcycle) Num = numOfKids;
    else Num = fNumOfABcycle;

    if (fEsetType == 1) /* Single-AB */
    {
        tRand->permutation(fPermu, fNumOfABcycle, fNumOfABcycle, rnd);
    }

    else if (fEsetType == 2)
    {  // Block2
        for (int k = 0; k < fNumOfABcycle; ++k) fNumOfElementINAB[k] = fABcycle[k][0];
        tSort->indexB(fNumOfElementINAB, fNumOfABcycle, fPermu, fNumOfABcycle);
    }
    fNumOfGeneratedCh = 0;
    pointMax = 0.0;
    BestGain = 0;
    fFlagImp = 0;

    ///    m_localsearch_indis.clear();
    for (int j = 0; j < Num; ++j)
    {
        fNumOfABcycleInEset = 0;
        if (fEsetType == 1) /* Single-AB */
        {
            jnum = fPermu[j];
            fABcycleInEset[fNumOfABcycleInEset++] = jnum;
        }
        else if (fEsetType == 2) /* Block2 */
        {
            jnum = fPermu[j];
            centerAB = jnum;
            for (int s = 0; s < fNumOfABcycle; ++s) {
                if (s == centerAB) fABcycleInEset[fNumOfABcycleInEset++] = s;
                else {
                    if (fWeightRR[centerAB][s] > 0 && fABcycle[s][0] < fABcycle[centerAB][0])
                        if (tRand->Double(0, 1, rnd) <= 0.5) fABcycleInEset[fNumOfABcycleInEset++] = s;
                }
            }
            this->searchEset(centerAB,rnd);
        }
        //fNumOfSPL = 0;
        fSegPosiList.clear();
        gain = 0;
        fNumOfAppliedCycle = 0;
        fNumOfModiEdge = 0;

        fNumOfAppliedCycle = fNumOfABcycleInEset;
        for (int k = 0; k < fNumOfAppliedCycle; ++k) {
            fAppliedCylce[k] = fABcycleInEset[k];
            jnum = fAppliedCylce[k];
            this->changeSol(tKid, jnum, flagP);
            gain += fGainAB[jnum];
        }

        this->makeUnit(rnd);
        this->makeCompleteSol(tKid,rnd);
        gain += fGainModi;

        ++fNumOfGeneratedCh;

        if (fEvalType == 1) /* Greedy */
            DLoss = 1.0;
        else if (fEvalType == 3) /* Distance preservation */
            DLoss = this->calAdpLoss(fEdgeFreq);
        else if (fEvalType == 4)  /* Entropy preservation */
            DLoss = this->calEntLoss(fEdgeFreq);

        if (DLoss <= 0.0) DLoss = 0.00000001;

        point = (double)gain / DLoss;
        tKid.fEvaluationValue = tKid.fEvaluationValue - gain;
        //        m_localsearch_indis.push_back(tKid);




        if (pointMax < point )
        {
            pointMax = point;
            BestGain = gain;
            fFlagImp = 1;

            fNumOfBestAppliedCycle = fNumOfAppliedCycle;
            for (int s = 0; s < fNumOfBestAppliedCycle; ++s) fBestAppliedCylce[s] = fAppliedCylce[s];

            fNumOfBestModiEdge = fNumOfModiEdge;
            for (int s = 0; s < fNumOfBestModiEdge; ++s) {
                fBestModiEdge[s][0] = fModiEdge[s][0];
                fBestModiEdge[s][1] = fModiEdge[s][1];
                fBestModiEdge[s][2] = fModiEdge[s][2];
                fBestModiEdge[s][3] = fModiEdge[s][3];
            }

        }
        this->backToPa1(tKid);
        tKid.fEvaluationValue = tKid.fEvaluationValue + gain;
    }
    if (fFlagImp == 1) {
        this->goToBest(tKid);
        tKid.fEvaluationValue = tKid.fEvaluationValue - BestGain;

        //  m_localsearch_indis.push_back(tKid);
        this->incrementEdgeFreq(fEdgeFreq);
    }
}

void TCross::setABcycle( const TIndi& tPa1, const TIndi& tPa2, const std::vector<int>& flagC,
    int numOfKids , Random* rnd)
{

    for (auto& it : nearData) {
        std::fill(it.begin(), it.end(), 0);
    }

    bunkiMany=0; koritsuMany=0;
    for( int j = 0; j < fN ; ++j ){
        nearData[j][1]=tPa1.fLink[j][0];
        nearData[j][3]=tPa1.fLink[j][1];
        nearData[j][0] = 2;

        koritsu[koritsuMany]=j;
        koritsuMany++;

        nearData[j][2]=tPa2.fLink[j][0];
        nearData[j][4]=tPa2.fLink[j][1];
    }
    for(int j = 0; j < fN; ++j ){
        checkKoritsu[j]=-1;
        koriInv[koritsu[j]]=j;
    }
    fNumOfABcycle=0;
    flagSt=1;
    while(koritsuMany!=0){
        if(flagSt==1){
            fPosiCurr=0;
            //r=rand()%koritsuMany;

            m_tmp_environemt->testRandom(rnd, "setABcycle 506");
            r = rnd->uniform.nextNonStd<int>(0, koritsuMany);
            m_tmp_environemt->testRandom(rnd, "setABcycle 508");
           // r = tRand->Integer(0, koritsuMany);
            st=koritsu[r];
            checkKoritsu[st]=fPosiCurr;
            fRoute[fPosiCurr]=st;
            ci=st;
            prType=2;
        }
        else if(flagSt==0) ci=fRoute[fPosiCurr];

        flagCycle=0;
        while(flagCycle==0){
            fPosiCurr++;
            pr=ci;
            switch(prType){
            case 1:
                ci=nearData[pr][fPosiCurr%2+1];
            break;
            case 2:
                //r=rand()%2;
             //   r = tRand->Integer(0, 3);
                m_tmp_environemt->testRandom(rnd, "setABcycle 529");
                //r = rnd->uniform.nextNonStd<int>(0, 2);
                r = 1;
                m_tmp_environemt->testRandom(rnd, "setABcycle 531");
                ci=nearData[pr][fPosiCurr%2+1+2*r];
                if(r==0) this->swap(nearData[pr][fPosiCurr%2+1],nearData[pr][fPosiCurr%2+3]);
            break;
            case 3:
                ci=nearData[pr][fPosiCurr%2+3];
            }
            fRoute[fPosiCurr]=ci;
            if(nearData[ci][0]==2){
                if(ci==st){
                    if(checkKoritsu[st]==0){
                    if((fPosiCurr-checkKoritsu[st])%2==0){
                        if(nearData[st][fPosiCurr%2+1]==pr) this->swap(nearData[ci][fPosiCurr%2+1],nearData[ci][fPosiCurr%2+3]);

                        stAppear = 1;
                        this->formABcycle();
                        if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto RETURN;
                        if( fNumOfABcycle == fMaxNumOfABcycle ) goto RETURN;

                        flagSt=0;
                        flagCycle=1;
                        prType=1;
                    }
                    else{
                        this->swap(nearData[ci][fPosiCurr%2+1],nearData[ci][fPosiCurr%2+3]);
                        prType=2;
                    }
                    checkKoritsu[st]=fPosiCurr;
                    }
                    else{
                    stAppear = 2;
                    this->formABcycle();
                    if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto RETURN;
                    if( fNumOfABcycle == fMaxNumOfABcycle ) goto RETURN;

                    flagSt=1;
                    flagCycle=1;
                    }
                }
                else if(checkKoritsu[ci]==-1) {
                    checkKoritsu[ci]=fPosiCurr;
                    if(nearData[ci][fPosiCurr%2+1]==pr) this->swap(nearData[ci][fPosiCurr%2+1],nearData[ci][fPosiCurr%2+3]);
                    prType=2;
                }
                else if(checkKoritsu[ci]>0){
                    this->swap(nearData[ci][fPosiCurr%2+1],nearData[ci][fPosiCurr%2+3]);
                    if((fPosiCurr-checkKoritsu[ci])%2==0){
                        stAppear = 1;
                        this->formABcycle();
                        if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto RETURN;
                        if( fNumOfABcycle == fMaxNumOfABcycle ) goto RETURN;

                        flagSt=0;
                        flagCycle=1;
                        prType=1;
                    }
                    else{
                        this->swap(nearData[ci][(fPosiCurr+1)%2+1],nearData[ci][(fPosiCurr+1)%2+3]);
                        prType=3;
                    }
                }
            }
            else if(nearData[ci][0]==1){
                if(ci==st){
                    stAppear = 1;
                    this->formABcycle();
                    if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto RETURN;
                    if( fNumOfABcycle == fMaxNumOfABcycle ) goto RETURN;
                    flagSt=1;
                    flagCycle=1;
                }
                else prType=1;
            }
        }
    }
    while(bunkiMany!=0){
        fPosiCurr=0;
        //r=rand()%bunkiMany;
        m_tmp_environemt->testRandom(rnd, "setABcycle 609");
    //    r = rnd->uniform.nextNonStd<int>(0, bunkiMany);
        r = bunkiMany / 2;
        m_tmp_environemt->testRandom(rnd, "setABcycle 611");
        //r = tRand->Integer(0, bunkiMany+1);
        st=bunki[r];
        fRoute[fPosiCurr]=st;
        ci=st;

        flagCycle=0;
        while(flagCycle==0){
            pr=ci;
            fPosiCurr++;
            ci=nearData[pr][fPosiCurr%2+1];
            fRoute[fPosiCurr]=ci;
            if(ci==st){
                stAppear = 1;
                this->formABcycle();
                if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto RETURN;
                if( fNumOfABcycle == fMaxNumOfABcycle ) goto RETURN;

                flagCycle=1;
            }
        }
    }
RETURN:
    if( fNumOfABcycle == fMaxNumOfABcycle ){
        m_tmp_environemt->testRandom(rnd, "exit???");
        printf( "fMaxNumOfABcycle(%d) must be increased\n", fMaxNumOfABcycle );
        exit( 1 );
    }

    m_tmp_environemt->testRandom(rnd, "tcross after set parent");
}

void TCross::formABcycle(){
    int j;
    int st_count;
    int edge_type;
    int st,ci, stock;
    int cem;
    double diff;

    if(fPosiCurr%2==0) edge_type=1;
    else edge_type=2;
    st=fRoute[fPosiCurr];
    cem=0;
    fC[cem]=st;

    st_count=0;
    while(1){
        cem++;
        fPosiCurr--;
        ci=fRoute[fPosiCurr];
        if(nearData[ci][0]==2){
            koritsu[koriInv[ci]]=koritsu[koritsuMany-1];
            koriInv[koritsu[koritsuMany-1]]=koriInv[ci];
            koritsuMany--;
            bunki[bunkiMany]=ci;
            bunInv[ci]=bunkiMany;
            bunkiMany++;
        }
        else if(nearData[ci][0]==1){
            bunki[bunInv[ci]]=bunki[bunkiMany-1];
            bunInv[bunki[bunkiMany-1]]=bunInv[ci];
            bunkiMany--;
        }

        nearData[ci][0]--;
        if(ci==st) st_count++;
        if(st_count==stAppear) break;
        fC[cem]=ci;
    }

    if(cem==2) return;

    fABcycle[fNumOfABcycle][0]=cem;

    if(edge_type==2){
        stock=fC[0];
        for( int j=0;j<cem-1;j++) fC[j]=fC[j+1];
        fC[cem-1]=stock;
    }

    for( int j=0;j<cem;j++) fABcycle[fNumOfABcycle][j+2]=fC[j];

    fABcycle[fNumOfABcycle][1]=fC[cem-1];
    fABcycle[fNumOfABcycle][cem+2]=fC[0];
    fABcycle[fNumOfABcycle][cem+3]=fC[1];

    fC[ cem ] = fC[ 0 ];
    fC[ cem+1 ] = fC[ 1 ];
    diff = 0;
    for( j = 0; j < cem/2; ++j ) diff = diff + eval->fEdgeDis[fC[2*j]][fC[1+2*j]] - eval->fEdgeDis[fC[1+2*j]][fC[2+2*j]];

    fGainAB[fNumOfABcycle] = diff;
    ++fNumOfABcycle;
}

void TCross::swap(int &x, int &y){
    int s=x;
    x=y;
    y=s;
}

void TCross::changeSol( TIndi& tKid, int ABnum, int type ){
    int j;
    int cem, r1, r2, b1, b2;
    int po_r1, po_r2, po_b1, po_b2;

    cem=fABcycle[ABnum][0];
    fC[0]=fABcycle[ABnum][0];

    if(type==2) for(j=0;j<cem+3;j++) fC[cem+3-j]=fABcycle[ABnum][j+1];
    else for(j=1;j<=cem+3;j++) fC[j]=fABcycle[ABnum][j];

    for(j=0;j<cem/2;j++){
        r1=fC[2+2*j];r2=fC[3+2*j];
        b1=fC[1+2*j];b2=fC[4+2*j];

        if(tKid.fLink[r1][0]==r2) tKid.fLink[r1][0]=b1;
        else tKid.fLink[r1][1]=b1;
        if(tKid.fLink[r2][0]==r1) tKid.fLink[r2][0]=b2;
        else tKid.fLink[r2][1]=b2;

        po_r1 = fInv[ r1 ];
        po_r2 = fInv[ r2 ];
        po_b1 = fInv[ b1 ];
        po_b2 = fInv[ b2 ];

        if (po_r1 == 0 && po_r2 == fN - 1) fSegPosiList.push_back(po_r1);
        else if (po_r1 == fN - 1 && po_r2 == 0) fSegPosiList.push_back(po_r2);
        else if (po_r1 < po_r2) fSegPosiList.push_back(po_r2);
        else if (po_r2 < po_r1) fSegPosiList.push_back(po_r1);


    //    if( po_r1 == 0 && po_r2 == fN-1 ) fSegPosiList[ fNumOfSPL++ ] = po_r1;
    //    else if( po_r1 == fN-1 && po_r2 == 0 ) fSegPosiList[ fNumOfSPL++ ] = po_r2;
    //    else if( po_r1 < po_r2 ) fSegPosiList[ fNumOfSPL++ ] = po_r2;
    //    else if( po_r2 < po_r1 ) fSegPosiList[ fNumOfSPL++ ] = po_r1;

        LinkBPosi[ po_r1 ][ 1 ] = LinkBPosi[ po_r1 ][ 0 ];
        LinkBPosi[ po_r2 ][ 1 ] = LinkBPosi[ po_r2 ][ 0 ];
        LinkBPosi[ po_r1 ][ 0 ] = po_b1;
        LinkBPosi[ po_r2 ][ 0 ] = po_b2;
    }
}

void TCross::makeCompleteSol( TIndi& tKid, Random* rnd){
    int j, j1, j2;
    int st ,pre, curr, next, a, b, c, d, aa, bb, a1, b1;
    int min_unit_city;
    int center_un, select_un;
    double diff, max_diff;
    int near_num, nearMax;

   // m_opt2_circles.clear();

    std::array<int, 4> cur_circle;

    fGainModi = 0;
    while( fNumOfUnit != 1 ){
        min_unit_city = fN + 12345;
        for( int u = 0; u < fNumOfUnit; ++u )
            if( fNumOfElementInUnit[ u ] < min_unit_city ){
                center_un = u;
                min_unit_city = fNumOfElementInUnit[ u ];
            }


        st = -1;
        fNumOfSegForCenter = 0;
        for( int s = 0; s < fNumOfSeg; ++s )
            if( fSegUnit[ s ] == center_un ){
                int posi = fSegment[ s ][ 0 ];
                st = fOrder[ posi ];
                fSegForCenter[  fNumOfSegForCenter++ ] = s;
            }
        curr = -1;
        next = st;
        fNumOfElementInCU = 0;
        while(1){
            pre = curr;
            curr = next;
            fCenterUnit[ curr ] = 1;
            fListOfCenterUnit[ fNumOfElementInCU ] = curr;
            ++fNumOfElementInCU;
            if( tKid.fLink[ curr ][ 0 ] != pre ) next = tKid.fLink[ curr ][ 0 ];
            else next = tKid.fLink[ curr ][ 1 ];
            if( next == st ) break;
        }
        fListOfCenterUnit[ fNumOfElementInCU ] = fListOfCenterUnit[ 0 ];
        fListOfCenterUnit[ fNumOfElementInCU+1 ] = fListOfCenterUnit[ 1 ];

        max_diff = -999999999;
        a1 = -1; b1 = -1;
        nearMax = 10; /* N_near (see Step 5.3 in Section 2.2 of the Online Supplement) */
        /* nearMax must be smaller than or equal to eva->fNearNumMax (kopt.cpp ) */

    RESTART:
        for( int s = 1; s <= fNumOfElementInCU; ++s ){
            a = fListOfCenterUnit[ s ];
            cur_circle[0] = a;
            for( near_num = 1; near_num <= nearMax; ++near_num ){
                c = eval->fNearCity[ a ][ near_num ];
                cur_circle[2] = c;
                if( fCenterUnit[ c ] == 0 ){
                    for( j1 = 0; j1 < 2; ++j1 ){
                        b = fListOfCenterUnit[ s-1+2*j1 ];
                        cur_circle[1] = b;
                        for( j2 = 0; j2 < 2; ++j2 ){
                            d = tKid.fLink[ c ][ j2 ];
                            cur_circle[3] = d;
                           // m_opt2_circles.push_back(cur_circle);
                            diff = eval->fEdgeDis[a][b] + eval->fEdgeDis[c][d] - eval->fEdgeDis[a][c] - eval->fEdgeDis[b][d];
                            if( diff > max_diff ){
                                aa = a; bb = b; a1 = c; b1 = d;
                                max_diff = diff;
                            }
                            diff = eval->fEdgeDis[a][b] + eval->fEdgeDis[d][c] -
                                eval->fEdgeDis[a][d] - eval->fEdgeDis[b][c];
                            if( diff > max_diff ){
                                aa = a; bb = b; a1 = d; b1 = c;
                                max_diff = diff;
                            }
                        }
                    }
                }
            }
        }

        if (a1 == -1 && nearMax == 10) /* This value must also be changed if nearMax is changed above */
        {
            nearMax = 50;
            goto RESTART;
        }
        else if( a1 == -1 && nearMax == 50  ){
            //int r = rand() % ( fNumOfElementInCU - 1 );
            int r = rnd->uniform.nextNonStd<int>(0, fNumOfElementInCU - 1);
            //int r = tRand->Integer(0, fNumOfElementInCU);
            a = fListOfCenterUnit[ r ];
            b = fListOfCenterUnit[ r+1 ];
            for( j = 0; j < fN; ++j ){
                if( fCenterUnit[ j ] == 0 ){
                    aa = a; bb = b;
                    a1 = j;
                    b1 = tKid.fLink[ j ][ 0 ];
                    break;
                }
            }
            max_diff = eval->fEdgeDis[aa][bb] + eval->fEdgeDis[a1][b1] - eval->fEdgeDis[a][a1] - eval->fEdgeDis[b][b1];
        }

        if( tKid.fLink[aa][0] == bb ) tKid.fLink[aa][0]=a1;
        else tKid.fLink[aa][1] = a1;
        if( tKid.fLink[bb][0] == aa ) tKid.fLink[bb][0] = b1;
        else tKid.fLink[bb][1] = b1;
        if( tKid.fLink[a1][0] == b1 ) tKid.fLink[a1][0] = aa;
        else tKid.fLink[a1][1] = aa;
        if( tKid.fLink[b1][0] == a1 ) tKid.fLink[b1][0] = bb;
        else tKid.fLink[b1][1] = bb;

        fModiEdge[ fNumOfModiEdge ][ 0 ] = aa;
        fModiEdge[ fNumOfModiEdge ][ 1 ] = bb;
        fModiEdge[ fNumOfModiEdge ][ 2 ] = a1;
        fModiEdge[ fNumOfModiEdge ][ 3 ] = b1;
        ++fNumOfModiEdge;

        fGainModi += max_diff;

        int posi_a1 = fInv[ a1 ];
        select_un = -1;
        for( int s = 0; s < fNumOfSeg; ++s )
            if( fSegment[ s ][ 0 ] <= posi_a1 && posi_a1 <=  fSegment[ s ][ 1 ] ){
                select_un = fSegUnit[ s ];
                break;
            }


        for( int s = 0; s < fNumOfSeg; ++s )
            if( fSegUnit[ s ] == select_un ) fSegUnit[ s ] = center_un;

        fNumOfElementInUnit[ center_un ] += fNumOfElementInUnit[ select_un ];

        for( int s = 0; s < fNumOfSeg; ++s )
            if( fSegUnit[ s ] == fNumOfUnit - 1 ) fSegUnit[ s ] = select_un;

        fNumOfElementInUnit[ select_un ] = fNumOfElementInUnit[ fNumOfUnit - 1 ];
        --fNumOfUnit;

        for( int s = 0; s < fNumOfElementInCU; ++s ){
            c = fListOfCenterUnit[ s ];
            fCenterUnit[ c ] = 0;
        }
    }
}

void TCross::makeUnit(Random* rnd){
    int flag = 1;
    for( int s = 0; s < fSegPosiList.size(); ++s ){
        if( fSegPosiList[ s ] == 0 ){
            flag = 0;
            break;
        }
    }
    if( flag == 1 ){
        fSegPosiList.push_back(0);
        //fSegPosiList[ fNumOfSPL++ ] = 0;
        LinkBPosi[ fN-1 ][ 1 ]  = LinkBPosi[ fN-1 ][ 0 ];
        LinkBPosi[ 0 ][ 1 ] = LinkBPosi[ 0 ][ 0 ];
        LinkBPosi[ fN-1 ][ 0 ] = 0;
        LinkBPosi[ 0 ][ 0 ] = fN-1;
    }

    tSort->sort( fSegPosiList, fSegPosiList.size(), rnd);
   // tSort->sort(fSegPosiList, fNumOfSPL);
  //  fNumOfSeg = fNumOfSPL;
    fNumOfSeg = fSegPosiList.size();
    for( int s = 0; s < fNumOfSeg-1; ++s ){
        fSegment[ s ][ 0 ] = fSegPosiList[ s ];
        fSegment[ s ][ 1 ] = fSegPosiList[ s+1 ]-1;
    }

    fSegment[ fNumOfSeg-1 ][ 0 ] = fSegPosiList[ fNumOfSeg-1 ];
    fSegment[ fNumOfSeg-1 ][ 1 ] = fN - 1;

    for( int s = 0; s < fNumOfSeg; ++s ){
        LinkAPosi[ fSegment[ s ][ 0 ] ] = fSegment[ s ][ 1 ];
        LinkAPosi[ fSegment[ s ][ 1 ] ] = fSegment[ s ][ 0 ];
        fPosiSeg[ fSegment[ s ][ 0 ] ] = s;
        fPosiSeg[ fSegment[ s ][ 1 ] ] = s;
    }

    for( int s = 0; s < fNumOfSeg; ++s ) fSegUnit[ s ] = -1;
    fNumOfUnit = 0;

    int p_st, p1, p2, p_next, p_pre;
    int segNum;
    while(1){
        flag = 0;
        for( int s = 0; s < fNumOfSeg; ++s ){
            if( fSegUnit[ s ] == -1 ){
                p_st = fSegment[ s ][ 0 ];
                p_pre = -1;
                p1 = p_st;
                flag = 1;
                break;
            }
        }
        if( flag == 0 ) break;

        while(1){
            segNum = fPosiSeg[ p1 ];
            fSegUnit[ segNum ] = fNumOfUnit;

            p2 = LinkAPosi[ p1 ];
            p_next = LinkBPosi[ p2 ][ 0 ];
            if( p1 == p2 )
                if( p_next == p_pre ) p_next = LinkBPosi[ p2 ][ 1 ];

            if( p_next == p_st ){
                ++fNumOfUnit;
                break;
            }

            p_pre = p2;
            p1 = p_next;
        }
    }

    for( int s = 0; s < fNumOfUnit; ++s ) fNumOfElementInUnit[ s ] = 0;

    int unitNum = -1;
    int tmpNumOfSeg = -1;
    for( int s = 0; s < fNumOfSeg; ++s ){
        if( fSegUnit[ s ] != unitNum ){
            ++tmpNumOfSeg;
            fSegment[ tmpNumOfSeg ][ 0 ] = fSegment[ s ][ 0 ];
            fSegment[ tmpNumOfSeg ][ 1 ] = fSegment[ s ][ 1 ];
            unitNum = fSegUnit[ s ];
            fSegUnit[ tmpNumOfSeg ] = unitNum;
            fNumOfElementInUnit[ unitNum ] +=
            fSegment[ s ][ 1 ] - fSegment[ s ][ 0 ] + 1;
        }
        else{
            fSegment[ tmpNumOfSeg ][ 1 ] = fSegment[ s ][ 1 ];
            fNumOfElementInUnit[ unitNum ] +=
            fSegment[ s ][ 1 ] - fSegment[ s ][ 0 ] + 1;
        }
    }
    fNumOfSeg = tmpNumOfSeg + 1;
}

void TCross::backToPa1( TIndi& tKid ){
    int aa, bb, a1, b1;
    int jnum;

    for( int s = fNumOfModiEdge -1; s >= 0; --s ){
        aa = fModiEdge[ s ][ 0 ];
        a1 = fModiEdge[ s ][ 1 ];
        bb = fModiEdge[ s ][ 2 ];
        b1 = fModiEdge[ s ][ 3 ];

        if( tKid.fLink[aa][0] == bb ) tKid.fLink[aa][0] = a1;
        else tKid.fLink[aa][1] = a1;
        if( tKid.fLink[b1][0] == a1 ) tKid.fLink[b1][0] = bb;
        else tKid.fLink[b1][1] = bb;
        if( tKid.fLink[bb][0] == aa ) tKid.fLink[bb][0] = b1;
        else tKid.fLink[bb][1] = b1;
        if( tKid.fLink[a1][0] == b1 ) tKid.fLink[a1][0] = aa;
        else tKid.fLink[a1][1] = aa;
    }

    for( int s = 0; s < fNumOfAppliedCycle; ++s ){
        jnum = fAppliedCylce[ s ];
        this->changeSol( tKid, jnum, 2 );
    }
}

void TCross::goToBest( TIndi& tKid ){
    int aa, bb, a1, b1;
    int jnum;

    for( int s = 0; s < fNumOfBestAppliedCycle; ++s ){
        jnum = fBestAppliedCylce[ s ];
        this->changeSol( tKid, jnum, 1 );
    }

     for( int s = 0; s < fNumOfBestModiEdge; ++s ){
        aa = fBestModiEdge[ s ][ 0 ];
        bb = fBestModiEdge[ s ][ 1 ];
        a1 = fBestModiEdge[ s ][ 2 ];
        b1 = fBestModiEdge[ s ][ 3 ];

        if( tKid.fLink[aa][0] == bb ) tKid.fLink[aa][0]=a1;
        else tKid.fLink[aa][1] = a1;
        if( tKid.fLink[bb][0] == aa ) tKid.fLink[bb][0] = b1;
        else tKid.fLink[bb][1] = b1;
        if( tKid.fLink[a1][0] == b1 ) tKid.fLink[a1][0] = aa;
        else tKid.fLink[a1][1] = aa;
        if( tKid.fLink[b1][0] == a1 ) tKid.fLink[b1][0] = bb;
        else tKid.fLink[b1][1] = bb;
    }
}

void TCross::incrementEdgeFreq(std::vector<std::vector<int>>& fEdgeFreq ){
    int j, jnum, cem;
    int r1, r2, b1, b2;
    int aa, bb, a1;

    for( int s = 0; s < fNumOfBestAppliedCycle; ++s ){
        jnum = fBestAppliedCylce[ s ];

        cem = fABcycle[ jnum ][ 0 ];
        fC[ 0 ] = fABcycle[ jnum ][ 0 ];

        for( j = 1; j <= cem+3; ++j )
            fC[ j ] = fABcycle[ jnum ][ j ];

        for( j = 0; j <cem/2; ++j ){
            r1 = fC[2+2*j]; r2 = fC[3+2*j];
            b1 = fC[1+2*j]; b2 = fC[4+2*j];

            // r1 - b1 add
            // r1 - r2 remove
            // r2 - r1 remove
            // r2 - b2 add

            ++fEdgeFreq[ r1 ][ b1 ];
            --fEdgeFreq[ r1 ][ r2 ];
            --fEdgeFreq[ r2 ][ r1 ];
            ++fEdgeFreq[ r2 ][ b2 ];
        }
    }
    for( int s = 0; s < fNumOfBestModiEdge; ++s ){
        aa = fBestModiEdge[ s ][ 0 ];
        bb = fBestModiEdge[ s ][ 1 ];
        a1 = fBestModiEdge[ s ][ 2 ];
        b1 = fBestModiEdge[ s ][ 3 ];

        --fEdgeFreq[ aa ][ bb ];
        --fEdgeFreq[ a1 ][ b1 ];
        ++fEdgeFreq[ aa ][ a1 ];
        ++fEdgeFreq[ bb ][ b1 ];
        --fEdgeFreq[ bb ][ aa ];
        --fEdgeFreq[ b1 ][ a1 ];
        ++fEdgeFreq[ a1 ][ aa ];
        ++fEdgeFreq[ b1 ][ bb ];
    }
}

int TCross::calAdpLoss(std::vector<std::vector<int>>& fEdgeFreq )
{
    int j, jnum, cem;
    int r1, r2, b1, b2;
    int aa, bb, a1;
    double DLoss;

    DLoss = 0;
    for( int s = 0; s < fNumOfAppliedCycle; ++s )
    {
        jnum = fAppliedCylce[ s ];

        cem = fABcycle[ jnum ][ 0 ];
        fC[ 0 ] = fABcycle[ jnum ][ 0 ];

        for( j = 1; j <= cem+3; ++j ) fC[ j ] = fABcycle[ jnum ][ j ];

        for( j = 0; j <cem/2; ++j ){
            r1 = fC[2+2*j]; r2 = fC[3+2*j];
            b1 = fC[1+2*j]; b2 = fC[4+2*j];

            // r1 - b1 add
            // r1 - r2 remove
            // r2 - r1 remove
            // r2 - b2 add

            DLoss -= (fEdgeFreq[ r1 ][ r2 ]-1);
            DLoss -= (fEdgeFreq[ r2 ][ r1 ]-1);
            DLoss += fEdgeFreq[ r2 ][ b2 ];
            DLoss += fEdgeFreq[ b2 ][ r2 ];

            // Remove
            --fEdgeFreq[ r1 ][ r2 ];
            --fEdgeFreq[ r2 ][ r1 ];

            // Add
            ++fEdgeFreq[ r2 ][ b2 ];
            ++fEdgeFreq[ b2 ][ r2 ];
        }
    }

    for( int s = 0; s < fNumOfModiEdge; ++s )
    {
        aa = fModiEdge[ s ][ 0 ];
        bb = fModiEdge[ s ][ 1 ];
        a1 = fModiEdge[ s ][ 2 ];
        b1 = fModiEdge[ s ][ 3 ];

        DLoss -= (fEdgeFreq[ aa ][ bb ]-1);
        DLoss -= (fEdgeFreq[ bb ][ aa ]-1);
        DLoss -= (fEdgeFreq[ a1 ][ b1 ]-1);
        DLoss -= (fEdgeFreq[ b1 ][ a1 ]-1);

        DLoss += fEdgeFreq[ aa ][ a1 ];
        DLoss += fEdgeFreq[ a1 ][ aa ];
        DLoss += fEdgeFreq[ bb ][ b1 ];
        DLoss += fEdgeFreq[ b1 ][ bb ];

        // Remove
        --fEdgeFreq[ aa ][ bb ];
        --fEdgeFreq[ bb ][ aa ];
        --fEdgeFreq[ a1 ][ b1 ];
        --fEdgeFreq[ b1 ][ a1 ];

        // Add
        ++fEdgeFreq[ aa ][ a1 ];
        ++fEdgeFreq[ a1 ][ aa ];
        ++fEdgeFreq[ bb ][ b1 ];
        ++fEdgeFreq[ b1 ][ bb ];
    }

    for( int s = 0; s < fNumOfAppliedCycle; ++s ){
        jnum = fAppliedCylce[ s ];
        cem = fABcycle[ jnum ][ 0 ];
        fC[ 0 ] = fABcycle[ jnum ][ 0 ];
        for( j = 1; j <= cem+3; ++j ) fC[ j ] = fABcycle[ jnum ][ j ];

        for( j = 0; j <cem/2; ++j ){
            r1 = fC[2+2*j]; r2 = fC[3+2*j];
            b1 = fC[1+2*j]; b2 = fC[4+2*j];

            ++fEdgeFreq[ r1 ][ r2 ];
            ++fEdgeFreq[ r2 ][ r1 ];
            --fEdgeFreq[ r2 ][ b2 ];
            --fEdgeFreq[ b2 ][ r2 ];
        }
    }

    // Modification
    for( int s = 0; s < fNumOfModiEdge; ++s )
    {
        aa = fModiEdge[ s ][ 0 ];
        bb = fModiEdge[ s ][ 1 ];
        a1 = fModiEdge[ s ][ 2 ];
        b1 = fModiEdge[ s ][ 3 ];

        // Remove
        ++fEdgeFreq[ aa ][ bb ];
        ++fEdgeFreq[ bb ][ aa ];
        ++fEdgeFreq[ a1 ][ b1 ];
        ++fEdgeFreq[ b1 ][ a1 ];

        --fEdgeFreq[ aa ][ a1 ];
        --fEdgeFreq[ a1 ][ aa ];
        --fEdgeFreq[ bb ][ b1 ];
        --fEdgeFreq[ b1 ][ bb ];
    }
    return int(DLoss / 2);
}

double TCross::calEntLoss(std::vector<std::vector<int>>& fEdgeFreq ){
    int j, jnum, cem;
    int r1, r2, b1, b2;
    int aa, bb, a1;
    double DLoss;
    double h1, h2;

    DLoss = 0;
    // AB-cycle
    for( int s = 0; s < fNumOfAppliedCycle; ++s ){
        jnum = fAppliedCylce[ s ];
        cem = fABcycle[ jnum ][ 0 ];
        fC[ 0 ] = fABcycle[ jnum ][ 0 ];

        for( j = 1; j <= cem+3; ++j ) fC[ j ] = fABcycle[ jnum ][ j ];

        for( j = 0; j <cem/2; ++j )
        {
            r1 = fC[2+2*j]; r2 = fC[3+2*j];
            b1 = fC[1+2*j]; b2 = fC[4+2*j];

            // r1 - b1 add
            // r1 - r2 remove
            // r2 - r1 remove
            // r2 - b2 add

            // Remove
            h1 = (double)(fEdgeFreq[r1][r2] - 1) / (double)Npop;
            h2 = (double)(fEdgeFreq[r1][r2]) / (double)Npop;
            if (fEdgeFreq[r1][r2] - 1 != 0)
                DLoss -= h1 * log(h1);
            DLoss += h2 * log(h2);
            --fEdgeFreq[r1][r2];
            --fEdgeFreq[r2][r1];

            // Add
            h1 = (double)(fEdgeFreq[r2][b2] + 1) / (double)Npop;
            h2 = (double)(fEdgeFreq[r2][b2]) / (double)Npop;
            DLoss -= h1 * log(h1);
            if (fEdgeFreq[r2][b2] != 0)
                DLoss += h2 * log(h2);
            ++fEdgeFreq[r2][b2];
            ++fEdgeFreq[b2][r2];
        }
    }

    // Modification
    for( int s = 0; s < fNumOfModiEdge; ++s ){
        aa = fModiEdge[ s ][ 0 ];
        bb = fModiEdge[ s ][ 1 ];
        a1 = fModiEdge[ s ][ 2 ];
        b1 = fModiEdge[ s ][ 3 ];

        h1 = (double)( fEdgeFreq[ aa ][ bb ] - 1 )/(double)Npop;
        h2 = (double)( fEdgeFreq[ aa ][ bb ] )/(double)Npop;
        if( fEdgeFreq[ aa ][ bb ] - 1 != 0 )
          DLoss -= h1 * log( h1 );
        DLoss += h2 * log( h2 );
        --fEdgeFreq[ aa ][ bb ];
        --fEdgeFreq[ bb ][ aa ];

        h1 = (double)( fEdgeFreq[ a1 ][ b1 ] - 1 )/(double)Npop;
        h2 = (double)( fEdgeFreq[ a1 ][ b1 ] )/(double)Npop;
        if( fEdgeFreq[ a1 ][ b1 ] - 1 != 0 )
          DLoss -= h1 * log( h1 );
        DLoss += h2 * log( h2 );
        --fEdgeFreq[ a1 ][ b1 ];
        --fEdgeFreq[ b1 ][ a1 ];

        h1 = (double)( fEdgeFreq[ aa ][ a1 ] + 1 )/(double)Npop;
        h2 = (double)( fEdgeFreq[ aa ][ a1 ])/(double)Npop;
        DLoss -= h1 * log( h1 );
        if( fEdgeFreq[ aa ][ a1 ] != 0 )
          DLoss += h2 * log( h2 );
        ++fEdgeFreq[ aa ][ a1 ];
        ++fEdgeFreq[ a1 ][ aa ];

        h1 = (double)( fEdgeFreq[ bb ][ b1 ] + 1 )/(double)Npop;
        h2 = (double)( fEdgeFreq[ bb ][ b1 ])/(double)Npop;
        DLoss -= h1 * log( h1 );
        if( fEdgeFreq[ bb ][ b1 ] != 0 )
          DLoss += h2 * log( h2 );
        ++fEdgeFreq[ bb ][ b1 ];
        ++fEdgeFreq[ b1 ][ bb ];
    }
    DLoss = -DLoss;

    for( int s = 0; s < fNumOfAppliedCycle; ++s ){
        jnum = fAppliedCylce[ s ];

        cem = fABcycle[ jnum ][ 0 ];
        fC[ 0 ] = fABcycle[ jnum ][ 0 ];

        for( j = 1; j <= cem+3; ++j ) fC[ j ] = fABcycle[ jnum ][ j ];

        for( j = 0; j <cem/2; ++j ){
          r1 = fC[2+2*j]; r2 = fC[3+2*j];
          b1 = fC[1+2*j]; b2 = fC[4+2*j];

          ++fEdgeFreq[ r1 ][ r2 ];
          ++fEdgeFreq[ r2 ][ r1 ];
          --fEdgeFreq[ r2 ][ b2 ];
          --fEdgeFreq[ b2 ][ r2 ];
        }
    }
    for( int s = 0; s < fNumOfModiEdge; ++s ){
        aa = fModiEdge[ s ][ 0 ];
        bb = fModiEdge[ s ][ 1 ];
        a1 = fModiEdge[ s ][ 2 ];
        b1 = fModiEdge[ s ][ 3 ];

        ++fEdgeFreq[ aa ][ bb ];
        ++fEdgeFreq[ bb ][ aa ];
        ++fEdgeFreq[ a1 ][ b1 ];
        ++fEdgeFreq[ b1 ][ a1 ];

        --fEdgeFreq[ aa ][ a1 ];
        --fEdgeFreq[ a1 ][ aa ];
        --fEdgeFreq[ bb ][ b1 ];
        --fEdgeFreq[ b1 ][ bb ];
    }
  return DLoss;
}

void TCross::setWeight( const TIndi& tPa1, const TIndi& tPa2 ){
    int cem;
    int r1, r2, v1, v2, v_p;
    int AB_num;

    for( int i = 0; i < fN; ++i ){
        fInEffectNode[ i ][ 0 ] = -1;
        fInEffectNode[ i ][ 1 ] = -1;
    }

    // Step 1:
    for( int s = 0; s < fNumOfABcycle; ++s ){
        cem = fABcycle[ s ][ 0 ];
        for( int j = 0; j < cem/2; ++j ){
            r1 = fABcycle[ s ][ 2*j+2 ];  // red edge
            r2 = fABcycle[ s ][ 2*j+3 ];

            if( fInEffectNode[ r1 ][ 0 ] == -1 ) fInEffectNode[ r1 ][ 0 ] = s;
            else if ( fInEffectNode[ r1 ][ 1 ] == -1 ) fInEffectNode[ r1 ][ 1 ] = s;


            if( fInEffectNode[ r2 ][ 0 ] == -1 ) fInEffectNode[ r2 ][ 0 ] = s;
            else if ( fInEffectNode[ r2 ][ 1 ] == -1 ) fInEffectNode[ r2 ][ 1 ] = s;

        }
    }

    // Step 2:
    for( int i = 0; i < fN; ++i ){
        if( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] == -1 ){
            AB_num = fInEffectNode[ i ][ 0 ];
            v1 = i;

            if( tPa1.fLink[ v1 ][ 0 ] != tPa2.fLink[ v1 ][ 0 ] && tPa1.fLink[ v1 ][ 0 ] != tPa2.fLink[ v1 ][ 1 ] )
                v_p = tPa1.fLink[ v1 ][ 0 ];
            else if( tPa1.fLink[ v1 ][ 1 ] != tPa2.fLink[ v1 ][ 0 ] && tPa1.fLink[ v1 ][ 1 ] != tPa2.fLink[ v1 ][ 1 ] )
                v_p = tPa1.fLink[ v1 ][ 1 ];


            while( 1 ){
                fInEffectNode[ v1 ][ 1 ] = AB_num;

                if( tPa1.fLink[ v1 ][ 0 ] != v_p ) v2 = tPa1.fLink[ v1 ][ 0 ];
                else if( tPa1.fLink[ v1 ][ 1 ] != v_p ) v2 = tPa1.fLink[ v1 ][ 1 ];


                if( fInEffectNode[ v2 ][ 0 ] == -1 ) fInEffectNode[ v2 ][ 0 ] = AB_num;
                else if( fInEffectNode[ v2 ][ 1 ] == -1 ) fInEffectNode[ v2 ][ 1 ] = AB_num;


                if( fInEffectNode[ v2 ][ 1 ] != -1 ) break;
                v_p = v1;
                v1 = v2;
            }
        }
    }

    // Step 3:

    for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
        fWeightC[ s1 ] = 0;
        for( int s2 = 0; s2 < fNumOfABcycle; ++s2 ) fWeightRR[ s1 ][ s2 ] = 0;
    }

    for( int i = 0; i < fN; ++i ){

        if( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] != -1 ){
            ++fWeightRR[ fInEffectNode[ i ][ 0 ] ][ fInEffectNode[ i ][ 1 ] ];
            ++fWeightRR[ fInEffectNode[ i ][ 1 ] ][ fInEffectNode[ i ][ 0 ] ];
        }
        if( fInEffectNode[ i ][ 0 ] != fInEffectNode[ i ][ 1 ] ){
            ++fWeightC[ fInEffectNode[ i ][ 0 ] ];
            ++fWeightC[ fInEffectNode[ i ][ 1 ] ];
        }
    }
    for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ) fWeightRR[ s1 ][ s1 ] = 0;
}

int TCross::calCNaive(){
    int count_C;
    int tt;

    count_C = 0;

    for( int i = 0; i < fN; ++i ){
        if( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] != -1 ){
            tt = 0;
            if( fUsedAB[ fInEffectNode[ i ][ 0 ] ] == 1 ) ++tt;
            if( fUsedAB[ fInEffectNode[ i ][ 1 ] ] == 1 ) ++tt;
            if( tt == 1 ) ++count_C;
        }
    }
    return count_C;
}

void TCross::searchEset( int centerAB, Random* rnd){
    int nIter, stagImp;
    int delta_weight, min_delta_weight_nt;
    int flag_AddDelete, flag_AddDelete_nt;
    int selected_AB, selected_AB_nt;
    int jnum;

    fNumC = 0;  // Number of C nodes in E-set
    fNumE = 0;  // Number of Edges in E-set

    fNumOfUsedAB = 0;
    for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
        fUsedAB[ s1 ] = 0;
        fWeightSR[ s1 ] = 0;
        fMovedAB[ s1 ] = 0;
    }

    for( int s = 0; s < fNumOfABcycleInEset; ++s ){
        jnum = fABcycleInEset[ s ];
        this->addAB( jnum );
    }
    fBestNumC = fNumC;
    fBestNumE = fNumE;

    stagImp = 0;
    nIter = 0;
    while( 1 ){
        ++nIter;
        min_delta_weight_nt = 99999999;
        flag_AddDelete = 0;
        flag_AddDelete_nt = 0;
        for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
            if( fUsedAB[ s1 ] == 0 && fWeightSR[ s1 ] > 0 ){
                delta_weight = fWeightC[ s1 ] - 2 * fWeightSR[ s1 ];
                if( fNumC + delta_weight < fBestNumC ){
                    selected_AB = s1;
                    flag_AddDelete = 1;
                    fBestNumC = fNumC + delta_weight;
                }
                if( delta_weight < min_delta_weight_nt && nIter > fMovedAB[ s1 ] ){
                    selected_AB_nt = s1;
                    flag_AddDelete_nt = 1;
                    min_delta_weight_nt = delta_weight;
                }
            }
            else if( fUsedAB[ s1 ] == 1 && s1 != centerAB ){
                delta_weight = - fWeightC[ s1 ] + 2 * fWeightSR[ s1 ];
                if( fNumC + delta_weight < fBestNumC ){
                    selected_AB = s1;
                    flag_AddDelete = -1;
                    fBestNumC = fNumC + delta_weight;
                }
                if( delta_weight < min_delta_weight_nt && nIter > fMovedAB[ s1 ] ){
                    selected_AB_nt = s1;
                    flag_AddDelete_nt = -1;
                    min_delta_weight_nt = delta_weight;
                }
            }
        }

        if( flag_AddDelete != 0 ){
            if( flag_AddDelete == 1 ) this->addAB( selected_AB );
            else if( flag_AddDelete == -1 ) this->deleteAB( selected_AB );

            fMovedAB[ selected_AB ] = nIter + tRand->Integer( 1, fTmax , rnd);

            fBestNumE = fNumE;

            fNumOfABcycleInEset = 0;
            for( int s1 = 0; s1 < fNumOfABcycle; ++s1 )
                if( fUsedAB[ s1 ] == 1 ) fABcycleInEset[ fNumOfABcycleInEset++ ] = s1;

            stagImp = 0;
        }
        else if( flag_AddDelete_nt != 0 ) {
            if( flag_AddDelete_nt == 1 ) this->addAB( selected_AB_nt );
            else if( flag_AddDelete_nt == -1 )

            this->deleteAB( selected_AB_nt );
            fMovedAB[ selected_AB_nt ] = nIter + tRand->Integer( 1, fTmax, rnd);
        }
        if( flag_AddDelete == 0 ) ++stagImp;
        if( stagImp == fMaxStag ) break;
    }
}

void TCross::addAB( int num ){
    fNumC += fWeightC[ num ] - 2 * fWeightSR[ num ];
    fNumE += fABcycle[ num ][ 0 ] / 2;

    fUsedAB[ num ] = 1;
    ++fNumOfUsedAB;
    for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ) fWeightSR[ s1 ] += fWeightRR[ s1 ][ num ];
}

void TCross::deleteAB( int num ){
    fNumC -= fWeightC[ num ] - 2 * fWeightSR[ num ];
    fNumE -= fABcycle[ num ][ 0 ] / 2;

    fUsedAB[ num ] = 0;
    --fNumOfUsedAB;
    for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ) fWeightSR[ s1 ] -= fWeightRR[ s1 ][ num ];
}
