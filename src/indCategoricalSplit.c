#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  "indCategoricalSplit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

double getDeltailCategoricalSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userY, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass)
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
	
	int res = CalculateCategoricalSplit(userY, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, nbrClass, nbrObsMinInControlClass, nbrObsMinInTreatmentClass, &nbrLeft, &nbrRight, &TauMaxLeft, &TauMaxRight,&maxClassLeft,&maxClassRigth);

	if (log && writeLogHeader)
	{
	  LogSplitHeaderCategoricalCSplit("SplitVal.csv");
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
        LogSplitInFileCategoricalSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, TauMaxLeft, TauMaxRight, errorTotal,maxClassLeft,maxClassRigth);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}

int CalculateCategoricalSplit(double * dataY, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int nbrClass, int nbrObsMinInControlClass, int nbrObsMinInTreatmentClass, int * nbrLeft, int * nbrRight, double * TauMaxLeft, double * TauMaxRight,int * maxClassLeft, int * maxClassRigth)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*TauMaxLeft = 0;
	*TauMaxRight = 0;
	*maxClassLeft = 0;
	*maxClassRigth = 0;
	
	double MeanLeftControl = 0;
	double MeanRightControl = 0;
	double MeanLeft = 0;
	double MeanRight = 0;
	int nbrLeftClass = 0;
	int nbrRightClass = 0;
	
	//I=0 pour control
	int ClassI = 0;
	int res = 0;
	res = CalculateCategoricalOneClassSplit(dataY, dataG, indCut, NodeNbrObrservation, LEFT, 0, nbrObsMinInControlClass, &nbrLeftClass, &MeanLeftControl);
	if (res != 1)
	{
	  //pas assez de traitement control gauche pour faire le calcul.
	  return 0;
	}
	res = CalculateCategoricalOneClassSplit(dataY, dataG, indCut, NodeNbrObrservation, RIGHT, 0, nbrObsMinInControlClass, &nbrRightClass, &MeanRightControl);
	if (res != 1)
	{
	  //pas assez de traitement control droite pour faire le calcul.
	  return 0;
	}
	
	for (ClassI = 1; ClassI <= nbrClass ; ClassI++)
	{
	  res = CalculateCategoricalOneClassSplit(dataY, dataG, indCut, NodeNbrObrservation, LEFT, ClassI, nbrObsMinInTreatmentClass, &nbrLeftClass, &MeanLeft);
	  if (res == 1)
	  {
	    if ((*maxClassLeft) == 0)
	    {
	      (*TauMaxLeft) = MeanLeft - MeanLeftControl;
        (*maxClassLeft) = ClassI;
	    }
	    else
	    {
	      if ((*TauMaxLeft) < (MeanLeft - MeanLeftControl))
	      {
	        (*TauMaxLeft) = MeanLeft - MeanLeftControl;
	        (*maxClassLeft) = ClassI;
	      }
	    }
	  }
	  
	  res = CalculateCategoricalOneClassSplit(dataY, dataG, indCut, NodeNbrObrservation, RIGHT, ClassI, nbrObsMinInTreatmentClass, &nbrRightClass, &MeanRight);
	  if (res == 1)
	  {
	    if ((*maxClassRigth) == 0)
	    {
	      (*TauMaxRight) = MeanRight - MeanRightControl;
	      (*maxClassRigth) = ClassI;
	    }
	    else
	    {
	      if ((*TauMaxRight) < (MeanRight - MeanRightControl))
	      {
	        (*TauMaxRight) = MeanRight - MeanRightControl;
	        (*maxClassRigth) = ClassI;
	      }
	    }
	  }
	}

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

int CalculateCategoricalOneClassSplit(double * dataY, double * dataG, char * indCut, int NodeNbrObrservation, char LeftOrRight, int TheClassI, int nbrObsMinInEachClass, int * nbrClass, double * Mean)
{
  *nbrClass = 0;
  *Mean = 0;
  
  double sumY = 0;
  
  for (int i = 0; i < NodeNbrObrservation; i++)
  {
    if (indCut[i] == LeftOrRight && dataG[i] == TheClassI)
    {
      (*nbrClass)++;
      sumY += dataY[i];
    }
  }
  
  if ((*nbrClass) <= 0 ||
      (*nbrClass) < nbrObsMinInEachClass)
  {
    return 0;
  }
  
  (*Mean) = sumY / (*nbrClass);
  
  return 1;
}

// int CalculateCategoricalOneClassSplit(double * dataY, double * dataG, char * indCut, int NodeNbrObrservation, int TheClassI, int nbrObsMinInEachClass, int * nbrLeftClass, int * nbrRightClass, double * MeanLeft, double * MeanRight)
// {
//   *nbrLeftClass = 0;
//   *nbrRightClass = 0;
//   *MeanLeft = 0;
//   *MeanRight = 0;
// 
//   double sumLeftY = 0;
//   double sumRightY = 0;
// 
//   for (int i = 0; i < NodeNbrObrservation; i++)
//   {
//   	if (indCut[i] == LEFT && dataG[i] == TheClassI)
//   	{
//   		(*nbrLeftClass)++;
//   		sumLeftY += dataY[i];
//   	}
//   	if (indCut[i] == RIGHT && dataG[i] == TheClassI)
//   	{
//   		(*nbrRightClass)++;
//   		sumRightY += dataY[i];
//   	}
//   }
// 
//   if ((*nbrLeftClass) <= 0 ||
//       (*nbrRightClass) <= 0 ||
//       (*nbrLeftClass) < nbrObsMinInEachClass ||
//       (*nbrRightClass) < nbrObsMinInEachClass)
//   {
//     return 0;
//   }
//   
//   (*MeanLeft) = sumLeftY / (*nbrLeftClass);
//   (*MeanRight) = sumRightY / (*nbrRightClass);
//   
//   return 1;
// }

int LogSplitInFileCategoricalSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double TauMaxLeft, double TauMaxRight, double errorTotal,int maxClassLeft, int maxClassRigth)
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

int LogSplitHeaderCategoricalCSplit(char * FileName)
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
