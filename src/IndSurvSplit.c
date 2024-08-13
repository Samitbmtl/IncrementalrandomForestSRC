#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  "IndSurvSplit.h"
#ifdef _OPENMP
#include <omp.h>
#endif
 
double getDelta(int thread_id,uint treeID,uint nodeid,uint covariate,double split,uint log, int writeLogHeader, double * userTime,double * userEvent,double * userGroup,char * userSplitIndicator,int NodeNbrObrservation,int incTailType, int incMinObs, int incMinObsControl,int incMinObsTreatment,double incTreatmentLowBound,double incTreatmentUpperBound)
{
    int nbrLeft = 0;
    int nbrRight = 0;
    double errorTotalLeft = 0;
    double errorTotalRight = 0;
    int SNbrLeft_Control = 0;
    int SNbrLeft_Traitement = 0;
    int SNbrRight_Control = 0;
    int SNbrRight_Traitement = 0;
    double CumulTLeft = 0;
    double CumulTRight = 0;
    double CumulDiffBrutLeft = 0;
    double CumulDiffBrutRight = 0;

    double MaxT_Left = 0;
    double UsedMaxT_Left = 0;
    double MaxT_Right = 0;
    double UsedMaxT_Right = 0;
        
    int nbrLeftTreatment = 0;
    int nbrLeftControl = 0;
    int nbrRightTreatment = 0;
    int nbrRightControl = 0;

    int nbrEventLeftTreatment = 0;
    int nbrEventLeftControl = 0;
    int nbrEventRightTreatment = 0;
    int nbrEventRightControl = 0;

    int res =  SortedCalculateSplit(userTime,
                                    userEvent,
                                    userGroup,
                                    userSplitIndicator,
                                    NodeNbrObrservation,
                                    incTailType, 
									incMinObs,
									incMinObsControl,
									incMinObsTreatment,
									incTreatmentLowBound,
									incTreatmentUpperBound,
                                    &nbrLeft,
                                    &nbrLeftTreatment,
                                    &nbrLeftControl,
                                    &nbrRight,
                                    &nbrRightTreatment,
                                    &nbrRightControl,
                                    &errorTotalLeft,
                                    &errorTotalRight,
                                    &SNbrLeft_Control,
                                    &SNbrLeft_Traitement,
                                    &SNbrRight_Control,
                                    &SNbrRight_Traitement,
                                    &CumulTLeft,
                                    &CumulTRight,
                                    &MaxT_Left,
                                    &UsedMaxT_Left,
                                    &MaxT_Right,
                                    &UsedMaxT_Right,
                                    &CumulDiffBrutLeft,
                                    &CumulDiffBrutRight,
                                    &nbrEventLeftTreatment,
                                    &nbrEventLeftControl,
                                    &nbrEventRightTreatment,
                                    &nbrEventRightControl);
    if (log && writeLogHeader)
    {
        LogSplitHeader("SplitVal.csv");
    }
    if (res == 1)
    {
        double ert = (errorTotalLeft - errorTotalRight);
        ert = (ert < 0 ? -ert : ert);
        double errorTotal = sqrt(nbrLeft * nbrRight) * ert;
		if (log)  
		{
			#ifdef _OPENMP
				#pragma omp critical (writelog)
			#endif
			{                          
        		LogSplitInFile("SplitVal.csv",thread_id, treeID,nodeid,covariate,split, nbrLeft, nbrLeftTreatment, nbrLeftControl, nbrRight, nbrRightTreatment, nbrRightControl, errorTotalLeft, errorTotalRight, errorTotal, errorTotal, errorTotal, SNbrLeft_Control, SNbrLeft_Traitement, SNbrRight_Control, SNbrRight_Traitement, CumulTLeft, CumulTRight, MaxT_Left, UsedMaxT_Left, MaxT_Right, UsedMaxT_Right, CumulDiffBrutLeft, CumulDiffBrutRight, nbrEventLeftTreatment, nbrEventLeftControl, nbrEventRightTreatment, nbrEventRightControl);
			}
		}
        return errorTotal;
    }
    {
        return 0;
    }
}

int SortedCalculateSplit(double * dataT, double * dataC, double * dataG, char * indCut, int NodeNbrObrservation,int incTailType, int incMinObs, int incMinObsControl,int incMinObsTreatment,double incTreatmentLowBound,double incTreatmentUpperBound, int * p_nbrLeft, int * p_nbrLeftTreatment, int * p_nbrLeftControl, int * p_nbrRight, int * p_nbrRightTreatment, int * p_nbrRightControl, double * p_errorTotalLeft, double * p_errorTotalRight, int * SNbrLeft_Control, int * SNbrLeft_Traitement, int * SNbrRight_Control, int * SNbrRight_Traitement, double * CumulTLeft, double * CumulTRight, double * MaxT_Left, double * UsedMaxT_Left, double * MaxT_Right, double * UsedMaxT_Right, double * CumulDiffBrutLeft, double * CumulDiffBrutRight, int * p_nbrEventLeftTreatment, int * p_nbrEventLeftControl, int * p_nbrEventRightTreatment, int * p_nbrEventRightControl)
{
	int nbrLeft = 0;
	int nbrLeftTreatment = 0;
	int nbrLeftControl = 0;
	int nbrRight = 0;
	int nbrRightTreatment = 0;
	int nbrRightControl = 0;
	int nbrEventLeftTreatment = 0;
	int nbrEventLeftControl = 0;
	int nbrEventRightTreatment = 0;
	int nbrEventRightControl = 0;
    
	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			nbrLeft++;
			if (dataG[i] == TREATMENT)
			{
				nbrLeftTreatment++;
				if (dataC[i] == EVENT)
					nbrEventLeftTreatment++;
			}
			else
			{
				nbrLeftControl++;
				if (dataC[i] == EVENT)
					nbrEventLeftControl++;
			}
		}
		else
		{
			nbrRight++;
			if (dataG[i] == TREATMENT)
			{
				nbrRightTreatment++;
				if (dataC[i] == EVENT)
					nbrEventRightTreatment++;
			}
			else
			{
				nbrRightControl++;
				if (dataC[i] == EVENT)
					nbrEventRightControl++;
			}
		}
	}
    
	*p_nbrLeft = nbrLeft;
	*p_nbrRight = nbrRight;
	*p_nbrLeftTreatment = nbrLeftTreatment;
	*p_nbrLeftControl = nbrLeftControl;
	*p_nbrRightTreatment = nbrRightTreatment;
	*p_nbrRightControl = nbrRightControl;

	*p_nbrEventLeftTreatment = nbrEventLeftTreatment;
	*p_nbrEventLeftControl = nbrEventLeftControl;
	*p_nbrEventRightTreatment = nbrEventRightTreatment;
	*p_nbrEventRightControl = nbrEventRightControl;

    if (nbrLeft < incMinObs || 
		nbrRight < incMinObs || 
		nbrLeftTreatment < incMinObsTreatment || 
		nbrRightTreatment < incMinObsTreatment || 
		nbrLeftControl < incMinObsControl || 
		nbrRightControl < incMinObsControl ||
		(double)nbrLeftTreatment / (double)nbrLeft < incTreatmentLowBound ||
		(double)nbrLeftTreatment / (double)nbrLeft > incTreatmentUpperBound ||
		(double)nbrRightTreatment / (double)nbrRight < incTreatmentLowBound ||
		(double)nbrRightTreatment / (double)nbrRight > incTreatmentUpperBound ||
		nbrEventLeftTreatment < 1 ||
		nbrEventLeftControl < 1 ||
		nbrEventLeftTreatment < 1 ||
		nbrEventRightControl < 1)
	{
		*p_errorTotalRight = 0.0;
		*p_errorTotalLeft = 0.0;
		return 0;
	}
    

	double * TLeft_Control = NULL;
	double * SLeft_Control = NULL;

	double * TLeft_Traitement = NULL;
	double * SLeft_Traitement = NULL;

	double * TRight_Control = NULL;
	double * SRight_Control = NULL;

	double * TRight_Traitement = NULL;
	double * SRight_Traitement = NULL;

	double IntegralLeft_Treatment = 0;
	double IntegralLeft_Control = 0;
	double IntegralRight_Treatment = 0;
	double IntegralRight_Control = 0;
	
	SortedKaplanMaier(dataT, dataC, dataG, CONTROL, NodeNbrObrservation, indCut, LEFT, &TLeft_Control, &SLeft_Control, SNbrLeft_Control);
	SortedKaplanMaier(dataT, dataC, dataG, TREATMENT, NodeNbrObrservation, indCut, LEFT, &TLeft_Traitement, &SLeft_Traitement, SNbrLeft_Traitement);
	SortedKaplanMaier(dataT, dataC, dataG, CONTROL, NodeNbrObrservation, indCut, RIGHT, &TRight_Control, &SRight_Control, SNbrRight_Control);
	SortedKaplanMaier(dataT, dataC, dataG, TREATMENT, NodeNbrObrservation, indCut, RIGHT, &TRight_Traitement, &SRight_Traitement, SNbrRight_Traitement);

	NewCompare2KaplanMaier(TLeft_Control, SLeft_Control, *SNbrLeft_Control, nbrLeftControl, TLeft_Traitement, SLeft_Traitement, *SNbrLeft_Traitement, nbrLeftTreatment,incTailType, p_errorTotalLeft, CumulDiffBrutLeft, CumulTLeft, MaxT_Left, UsedMaxT_Left, &IntegralLeft_Control, &IntegralLeft_Treatment);
	NewCompare2KaplanMaier(TRight_Control, SRight_Control, *SNbrRight_Control, nbrRightControl, TRight_Traitement, SRight_Traitement, *SNbrRight_Traitement, nbrRightTreatment,incTailType, p_errorTotalRight, CumulDiffBrutRight, CumulTRight, MaxT_Right, UsedMaxT_Right, &IntegralRight_Control, &IntegralRight_Treatment);
	
	free(TLeft_Control);
	TLeft_Control = NULL;
	free(SLeft_Control);
	SLeft_Control = NULL;
	free(TLeft_Traitement);
	TLeft_Traitement = NULL;
	free(SLeft_Traitement);
	SLeft_Traitement = NULL;

	free(TRight_Control);
	TRight_Control = NULL;
	free(SRight_Control);
	SRight_Control = NULL;
	free(TRight_Traitement);
	TRight_Traitement = NULL;
	free(SRight_Traitement);
	SRight_Traitement = NULL;

	return 1;
}


//V18b
//Modification pour ajouter le Tail.
//Doit être lu comme (int * TLeft, double * SLeft, int SNbrLeft, int NbrLeft) => Partie control 
//et (int * TRight, double * SRight, int SNbrRight, int NbrRight) => Partie traitement
int NewCompare2KaplanMaier(double * TLeft, double * SLeft, int SNbrLeft, int NbrLeft, double * TRight, double * SRight, int SNbrRight, int NbrRight,int Tail, double * CumulDiff, double * CumulDiffBrut, double * CumulT, double * MaxT, double * UsedMaxT, double * IntegralLeft, double * IntegralRight)
{
	double IntLeft = TLeft[1];
	for (int i = 1; i < SNbrLeft - 1; i++)
	{
		IntLeft += SLeft[i] * (TLeft[i + 1] - TLeft[i]);
	}
	if (SLeft[SNbrLeft - 1] != 0.0 && SLeft[SNbrLeft - 1] != 1.0)
	{
        if (Tail == 1)
        {
		    IntLeft -= TLeft[SNbrLeft - 1] * SLeft[SNbrLeft - 1] / log(SLeft[SNbrLeft - 1]);
        }
	}

	double IntRight = TRight[1];
	for (int i = 1; i < SNbrRight - 1; i++)
	{
		IntRight += SRight[i] * (TRight[i + 1] - TRight[i]);
	}
	if (SRight[SNbrRight - 1] != 0.0 && SRight[SNbrRight - 1] != 1.0)
	{
        if (Tail == 1)
        {
		    IntRight -= TRight[SNbrRight - 1] * SRight[SNbrRight - 1] / log(SRight[SNbrRight - 1]);
        }
	}

	(*IntegralLeft) = IntLeft;
	(*IntegralRight) = IntRight;

	*CumulDiffBrut = IntRight - IntLeft;
	*CumulDiff = (*CumulDiffBrut);

	*CumulT = -1;
	*MaxT = -1;
	*UsedMaxT = -1;

	return 1;
}


int SortedKaplanMaier(double * tTra, double * cTra, double * gTra, int IsTreatmentGroup, int nbr, char * indObs, char BelowCut, double ** T, double ** S, int * SNbr)
{
	//tTra should be sorted (of lenth nbr)

	double * tTraSorted = (double *) malloc(nbr* sizeof(double));
    //(double *) malloc((size_t) ((nh+1) * (sizeof(double))));
	double * cTraSorted = (double *) malloc(nbr * sizeof(double));
	//u contient le count exact des valeur qui sont dans chaque tableau (incluant les doublons)
	int u = 0;

	for (int i = 0; i < nbr; i++)
	{
		if (indObs[i] == BelowCut && gTra[i] == IsTreatmentGroup)
		{
			tTraSorted[u] = tTra[i];
			cTraSorted[u] = cTra[i];
			u++;
		}
	}

	double * D  = (double *) malloc((u + 1) * sizeof(double));
	double * TT = (double *) malloc((u + 1) * sizeof(double));
	double * SS = (double *) malloc((u + 1) * sizeof(double));

	TT[0] = 0;
	SS[0] = 1;

	double LastS = SS[0];
	int n = u;
	int i = 0;
	int k = 1;
	int j;
	int m;

	while (i < u) {
		D[i] = 0.0;
		if (cTraSorted[i] == 1)
		{
			D[i] = 1;
			j = i + 1;
			m = 0;
			while ((tTraSorted[j] == tTraSorted[i]) && (j < u)) {
				if (cTraSorted[j] == 1)
				{
					D[i] = D[i] + 1;
				}
				j++;
				m++;
			}
			LastS = LastS * (n - D[i]) / n;
			n = n - m - 1;
			TT[k] = tTraSorted[i];
			SS[k] = LastS;
		}
		else
		{
			j = i + 1;
			m = 0;
			while ((tTraSorted[j] == tTraSorted[i]) && (j < u)) {
				if (cTraSorted[j] == 1)
				{
					D[i] = D[i] + 1;
				}
				j++;
				m++;
			}
			LastS = LastS * (n - D[i]) / n;
			n = n - m - 1;
			TT[k] = tTraSorted[i];
			SS[k] = LastS;
		}
		i = i + m + 1;
		k++;
	}

	*T = (double *) malloc(k * sizeof(double));
	*S = (double *) malloc(k * sizeof(double));

	for (int i = 0; i < k; i++)
	{
		(*T)[i] = TT[i];
		(*S)[i] = SS[i];
	}

	*SNbr = k;

	free((double *) TT);
	TT = NULL;
	free((double *) SS);
	SS = NULL;
	free((double *) D);
	D = NULL;
	free((double *) tTraSorted);
	tTraSorted = NULL;
	free((double *) cTraSorted);
	cTraSorted = NULL;

	return 1;
}

int LogSplitInFile(char * FileName,int thread_id, uint treeID,uint nodeid,uint covariate,double split, int nbrLeft, int nbrLeftTreatment, int nbrLeftControl, int nbrRight, int nbrRightTreatment, int nbrRightControl, double errorLeft, double errorRight, double error, double error1, double error2, int SNbrLeft_Control, int SNbrLeft_Traitement, int SNbrRight_Control, int SNbrRight_Traitement, int CumulTLeft, int CumulTRight, int MaxT_Left, int UsedMaxT_Left, int MaxT_Right, int UsedMaxT_Right, double CumulDiffBrutLeft, double CumulDiffBrutRight, int nbrEventLeftTreatment, int nbrEventLeftControl, int nbrEventRightTreatment, int nbrEventRightControl)
{
	FILE *fp;

	fp = fopen(FileName, "a");
	if (fp == NULL)
	{
		return -1;
	}

	fprintf(fp, "%i;%u;%u;%u;%f;%i;%i;%i;%i;%i;%i;%f;%f;%f;%f;%f;%f;%f;%i;%i;%i;%i;%i;%i;%i;%i;%i;%i;%i;%i;%i;%i\r",thread_id, treeID,nodeid,covariate,split,nbrLeft, nbrLeftTreatment, nbrLeftControl, nbrRight, nbrRightTreatment, nbrRightControl, errorLeft, errorRight, CumulDiffBrutLeft, CumulDiffBrutRight, error, error1, error2, SNbrLeft_Control, SNbrLeft_Traitement, SNbrRight_Control, SNbrRight_Traitement, CumulTLeft, CumulTRight, MaxT_Left, UsedMaxT_Left, MaxT_Right, UsedMaxT_Right, nbrEventLeftTreatment, nbrEventLeftControl, nbrEventRightTreatment, nbrEventRightControl);

	fclose(fp);
	return 1;
}

int LogSplitHeader(char * FileName)
{
	FILE *fp;

	fp = fopen(FileName, "w");
	if (fp == NULL)
	{
		return -1;
	}

	fprintf(fp, "threadID;treeID;nodeid;covariate;split;NbrLeftNode;NbrTreatmentLeftNode;NbrControlLeftNode;NbrRightNode;NbrTreatmentRightNode;NbrControlRightNode;ErrorLeft;ErrorRight;CumulDiffBrutLeft;CumulDiffBrutRight;TotalError;TotalError1;TotalError2;NbrTLeft_Control;NbrTLeft_Traitement;NbrTRight_Control;NbrTRight_Traitement;CumulTLeft;CumulTRight;MaxT_Left; UsedMaxT_Left; MaxT_Right; UsedMaxT_Right; nbrEventLeftTreatment; nbrEventLeftControl; nbrEventRightTreatment; nbrEventRightControl\r");

	fclose(fp);
	return 1;
}
