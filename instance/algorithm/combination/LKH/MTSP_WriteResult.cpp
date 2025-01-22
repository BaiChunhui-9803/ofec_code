#include "./INCLUDE/LKH.h"
namespace LKH {
	void MTSP_WriteResult(char * FileName, GainType Penalty, GainType Cost, LKHAlg *Alg)
	{
		FILE *ResultFile;
		LKHAlg::Node *N, *NextN;
		int Size, Forward;
		char *FullFileName;
		GainType Sum;

		if (FileName == 0)
			return;
		FullFileName = Alg->FullName(FileName, Cost);
		if (Alg->TraceLevel >= 1)
			Alg->printff("Writing MTSP_SOLUTION_FILE: \"%s\" ... ", FullFileName);
		assert(ResultFile = fopen(FullFileName, "w"));
		fprintf(ResultFile, "%s, Cost: " GainFormat "_" GainFormat "\n",
			Alg->Name, Penalty, Cost);
		fprintf(ResultFile, "The tours traveled by the %d salesmen are:\n",
			Alg->Salesmen);
		N = Alg->Depot;
		Forward = N->Suc->Id != N->Id + Alg->DimensionSaved;
		do {
			Sum = 0;
			Size = -1;
			do {
				fprintf(ResultFile, "%d ", N->Id <= Alg->Dim ? N->Id : Alg->Depot->Id);
				NextN = Forward ? N->Suc : N->Pred;
				Sum += (Alg->*(Alg->C))(N, NextN) - N->Pi - NextN->Pi;
				Size++;
				if (NextN->Id > Alg->DimensionSaved)
					NextN = Forward ? NextN->Suc : NextN->Pred;
				N = NextN;
			} while (N->DepotId == 0);
			fprintf(ResultFile, "%d (#%d)  Cost: " GainFormat "\n",
				Alg->Depot->Id, Size, Sum / Alg->Precision);
		} while (N != Alg->Depot);
		fclose(ResultFile);
		if (Alg->TraceLevel >= 1)
			Alg->printff("done\n");
	}
}