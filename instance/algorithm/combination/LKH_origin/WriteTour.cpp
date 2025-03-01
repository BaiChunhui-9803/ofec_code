#include "./INCLUDE/LKH.h"
#include<cstring>
/*
 * The WriteTour function writes a tour to file. The tour 
 * is written in TSPLIB format to file FileName. 
 * 
 * The tour is written in "normal form": starting at node 1,
 * and continuing in direction of its lowest numbered
 * neighbor.
 * 
 * Nothing happens if FileName is 0. 
 */

void LKH::LKHAlg::WriteTour(char *FileName, int *Tour, GainType Cost)
{
    FILE *TourFile;
    int i, j, n, Forward;
    char *FullFileName;
    time_t Now;

    if (CurrentPenalty != 0 && MTSPObjective == -1 &&
        ProblemType != CCVRP && ProblemType != TRP &&
        ProblemType != MLP)
        return;
    if (FileName == 0)
        return;
    FullFileName = FullName(FileName, Cost);
    Now = time(&Now);
    if (TraceLevel >= 1)
        printff("Writing%s: \"%s\" ... ",
                FileName == TourFileName ? " TOUR_FILE" :
                FileName == OutputTourFileName ? " OUTPUT_TOUR_FILE" : "",
                FullFileName);
    assert(TourFile = fopen(FullFileName, "w"));
    if (CurrentPenalty == 0) {
        fprintf(TourFile, "NAME : %s." GainFormat ".tour\n", Name, Cost);
        fprintf(TourFile, "COMMENT : Length = " GainFormat "\n", Cost);
    } else {
        fprintf(TourFile, "NAME : %s." GainFormat "_" GainFormat ".tour\n",
                Name, CurrentPenalty, Cost);
        fprintf(TourFile,
                "COMMENT : Cost = " GainFormat "_" GainFormat "\n",
                CurrentPenalty, Cost);
    }
    fprintf(TourFile, "COMMENT : Found by LKH [Keld Helsgaun] %s",
            ctime(&Now));
    fprintf(TourFile, "TYPE : TOUR\n");
    fprintf(TourFile, "DIMENSION : %d\n", DimensionSaved);
    fprintf(TourFile, "TOUR_SECTION\n");

    n = DimensionSaved;
    for (i = 1; i < n && Tour[i] != 1; i++);
    Forward = Asymmetric ||
        Tour[i < n ? i + 1 : 1] < Tour[i > 1 ? i - 1 : Dimension];
    for (j = 1; j <= n; j++) {
        if (Tour[i] <= n)
            fprintf(TourFile, "%d\n", Tour[i]);
        if (Forward) {
            if (++i > n)
                i = 1;
        } else if (--i < 1)
            i = n;
    }
    fprintf(TourFile, "-1\nEOF\n");
    fclose(TourFile);
    free(FullFileName);
    if (TraceLevel >= 1)
        printff("done\n");
}

void LKH::LKHAlg::transferTourType(int* Tour, std::vector<int>& vTour)const {
    
    vTour.clear();
    int i, j, n, Forward;
    

    n = DimensionSaved;
    for (i = 1; i < n && Tour[i] != 1; i++);
    Forward = Asymmetric ||
        Tour[i < n ? i + 1 : 1] < Tour[i > 1 ? i - 1 : Dimension];
    for (j = 1; j <= n; j++) {
        if (Tour[i] <= n)
            vTour.push_back(Tour[i]);
        //    fprintf(TourFile, "%d\n", Tour[i]);
        if (Forward) {
            if (++i > n)
                i = 1;
        }
        else if (--i < 1)
            i = n;
    }

   // for (auto& it : vTour) --it;
}

/*
 * The FullName function returns a copy of the string Name where all 
 * occurrences of the character '$' have been replaced by Cost.        
 */

char *LKH::LKHAlg::FullName(char *Name, GainType Cost)
{
    char *NewName = 0, *CostBuffer, *Pos;

    if (!(Pos = strstr(Name, "$"))) {
        assert(NewName = (char *) calloc(strlen(Name) + 1, 1));
        strcpy(NewName, Name);
        return NewName;
    }
    assert(CostBuffer = (char *) malloc(400));
    if (CurrentPenalty != 0)
        sprintf(CostBuffer, GainFormat "_" GainFormat,
                CurrentPenalty, Cost);
    else
        sprintf(CostBuffer, GainFormat, Cost);
    do {
        free(NewName);
        assert(NewName =
               (char *) calloc(strlen(Name) + strlen(CostBuffer) + 1, 1));
        strncpy(NewName, Name, Pos - Name);
        strcat(NewName, CostBuffer);
        strcat(NewName, Pos + 1);
        Name = NewName;
    }
    while ((Pos = strstr(Name, "$")));
    free(CostBuffer);
    return NewName;
}
