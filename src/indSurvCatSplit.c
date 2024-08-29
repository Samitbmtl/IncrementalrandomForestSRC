#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  "indSurvCatSplit.h"
#include  "IndSurvSplit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

double getDeltaSurvCatSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime,double * userEvent,double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incTailType, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double TauMaxLeft = 0;
	double TauMaxRight = 0;
	int maxClassLeft = 0;
	int maxClassRigth = 0;
	
	//Il faut au moins 2 classe, une contol et une traitement. 
	if (nbrClass < 1)
	{
	  return 0;
	}
	
	int res = CalculateSurvCatSplit(userTime,userEvent, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, nbrClass, nbrObsMinInControlClass, nbrObsMinInTreatmentClass,incTailType, &nbrLeft, &nbrRight, &TauMaxLeft, &TauMaxRight,&maxClassLeft,&maxClassRigth);

	if (log && writeLogHeader)
	{
	  LogSplitHeaderSurvCatSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		double ert = (TauMaxLeft - TauMaxRight);
		ert = (ert < 0 ? -ert : ert);
		double errorTotal = sqrt(nbrLeft * nbrRight) * ert;
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
        LogSplitInFileSurvCatSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, TauMaxLeft, TauMaxRight, errorTotal,maxClassLeft,maxClassRigth);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}

int CalculateSurvCatSplit(double * dataT,double * dataC, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass, int incTailType, int * nbrLeft, int * nbrRight, double * TauMaxLeft, double * TauMaxRight,int * maxClassLeft, int * maxClassRigth)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*TauMaxLeft = 0;
	*TauMaxRight = 0;
	*maxClassLeft = 0;
	*maxClassRigth = 0;
	
	double * TLeft_Control = NULL;
	double * SLeft_Control = NULL;
	int SNbrLeft_Control = 0;
	
	double * TRight_Control = NULL;
	double * SRight_Control = NULL;
	int SNbrRight_Control = 0;
	
	double * TLeft_Tr = NULL;
	double * SLeft_Tr = NULL;
	int SNbrLeft_Tr = 0;
	
	double * TRight_Tr = NULL;
	double * SRight_Tr = NULL;
	int SNbrRight_Tr = 0;
	
	double CumulDiff = 0;
	double CumulDiffBrut = 0;
	double CumulT = 0;
	double MaxT = 0;
	double UsedMaxT = 0;
	double IntegralLeft = 0;
	double IntegralRight = 0;
	
	//I=0 pour control
	int ClassI = 0;
	int res = 0;
	res = CalculateSurvCatOneClassSplit(dataT, dataC, dataG, indCut, NodeNbrObrservation, LEFT, 0, nbrObsMinInControlClass, &TLeft_Control, &SLeft_Control, &SNbrLeft_Control);
	if (res != 1)
	{
	  //pas assez de traitement control gauche pour faire le calcul.
	  return 0;
	}
	res = CalculateSurvCatOneClassSplit(dataT, dataC, dataG, indCut, NodeNbrObrservation, RIGHT, 0, nbrObsMinInControlClass, &TRight_Control, &SRight_Control, &SNbrRight_Control);
	if (res != 1)
	{
	  //pas assez de traitement control droite pour faire le calcul.
	  return 0;
	}
	
	for (ClassI = 1; ClassI <= nbrClass ; ClassI++)
	{
	  res = CalculateSurvCatOneClassSplit(dataT, dataC, dataG, indCut, NodeNbrObrservation, LEFT, ClassI, nbrObsMinInTreatmentClass, &TLeft_Tr, &SLeft_Tr, &SNbrLeft_Tr);
	  if (res == 1)
	  {
		  double val = NewCompare2KaplanMaier(TLeft_Control, SLeft_Control, SNbrLeft_Control, SNbrLeft_Control, TLeft_Tr, SLeft_Tr, SNbrLeft_Tr, SNbrLeft_Tr,incTailType, &CumulDiff, &CumulDiffBrut, &CumulT, &MaxT, &UsedMaxT, &IntegralLeft, &IntegralRight);
	    if ((*maxClassLeft) == 0)
	    {	
	      (*TauMaxLeft) = CumulDiff;
		  (*maxClassLeft) = ClassI;
	    }
	    else
	    {
	      if ((*TauMaxLeft) < (CumulDiff))
	      {
	        (*TauMaxLeft) = CumulDiff;
	        (*maxClassLeft) = ClassI;
	      }
	    }
	    if (val == 1) {}
	  }
	  
	  res = CalculateSurvCatOneClassSplit(dataT, dataC, dataG, indCut, NodeNbrObrservation, RIGHT, ClassI, nbrObsMinInTreatmentClass, &TRight_Tr, &SRight_Tr, &SNbrRight_Tr);
	  if (res == 1)
	  {
		  double val = NewCompare2KaplanMaier(TRight_Control, SRight_Control, SNbrRight_Control, SNbrRight_Control, TRight_Tr, SRight_Tr, SNbrRight_Tr, SNbrRight_Tr,incTailType, &CumulDiff, &CumulDiffBrut, &CumulT, &MaxT, &UsedMaxT, &IntegralLeft, &IntegralRight);
	    if ((*maxClassRigth) == 0)
	    {
	      (*TauMaxRight) = CumulDiff;
	      (*maxClassRigth) = ClassI;
	    }
	    else
	    {
	      if ((*TauMaxRight) < (CumulDiff))
	      {
	        (*TauMaxRight) = CumulDiff;
	        (*maxClassRigth) = ClassI;
	      }
	    }
	    if (val == 1) {}
	  }
	  
		free(TLeft_Tr);
		TLeft_Tr = NULL;
		free(SLeft_Tr);
		SLeft_Tr = NULL;
		free(TRight_Tr);
		TRight_Tr = NULL;
		free(SRight_Tr);
		SRight_Tr = NULL;
	
	}

	free(TLeft_Control);
	TLeft_Control = NULL;
	free(SLeft_Control);
	SLeft_Control = NULL;
	free(TRight_Control);
	TRight_Control = NULL;
	free(SRight_Control);
	SRight_Control = NULL;
		
	if ((*maxClassLeft) == 0 || (*maxClassRigth) == 0)
	{
		return 0;
	}
  
  //Pour retourner le nombre total d'observtion a gauche et a droite.
	for (int i = 0; i < NodeNbrObrservation; i++)
	{
	  if (indCut[i] == LEFT)
	  {
	    (*nbrLeft)++;
	  }
	  if (indCut[i] == RIGHT)
	  {
	    (*nbrRight)++;
	  }
	}
	
	if ((*nbrLeft) < incMinObs ||
     (*nbrRight) < incMinObs)
	{
	  return 0;
	}

	return 1;
}

int CalculateSurvCatOneClassSplit(double * dataT,double * dataC, double * dataG, char * indCut, int NodeNbrObrservation, char LeftOrRight, int TheClassI, int nbrObsMinInEachClass, double ** T, double ** S, int * ClassNbr)
{	
  SortedKaplanMaier(dataT, dataC, dataG, TheClassI, NodeNbrObrservation, indCut, LeftOrRight, T, S, ClassNbr);
  
  //ClassNbr toujour supérieur au nombre reel d'observation de 1 (qui provient de SortedKaplanMaier)
  //d'ou (*ClassNbr) <= nbrObsMinInEachClass et non pas (*ClassNbr) < nbrObsMinInEachClass
  
  if ((*ClassNbr) <= 0 ||
      (*ClassNbr) <= nbrObsMinInEachClass)
  {
    return 0;
  }
  
  return 1;
}

int LogSplitInFileSurvCatSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double TauMaxLeft, double TauMaxRight, double errorTotal,int maxClassLeft, int maxClassRigth)
{
	FILE *fp;

	fp = fopen(FileName, "a");
	if (fp == NULL)
	{
		return -1;
	}
	//thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight,maxClassLeft,maxClassRigth, TauMaxLeft, TauMaxRight, errorTotal);
	fprintf(fp, "%i;%u;%u;%u;%f;%i;%i;%i;%i;%f;%f;%f\r", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight,maxClassLeft,maxClassRigth, TauMaxLeft, TauMaxRight, errorTotal);

	fclose(fp);
	return 1;
}

int LogSplitHeaderSurvCatSplit(char * FileName)
{
	FILE *fp;

	fp = fopen(FileName, "w");
	if (fp == NULL)
	{
		return -1;
	}

	fprintf(fp, "threadID;treeID;nodeid;covariate;split;NbrLeftNode;NbrRightNode;MaxClassLeft;axClassRight;ErrorLeft;ErrorRight;TotalError\r");

	fclose(fp);
	return 1;
}
