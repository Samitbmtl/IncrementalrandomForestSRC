#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  "indLinearCSplit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

double getDeltaindLinearCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;

	int res = CalculateindLinearCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta1Left, &Beta1Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindLinearCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		double ert = (Beta1Left - Beta1Right);
		ert = (ert < 0 ? -ert : ert);
		double errorTotal = sqrt(nbrLeft * nbrRight) * ert;
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindLinearCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindTSTCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;

	int res = CalculateindTSTCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta1Left, &Beta1Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindLinearCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		//maximizing a function is equivalent to minimizing its negative
		double errorTotal = -(Beta1Left + Beta1Right);
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindLinearCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindNSRCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;

	int res = CalculateindLinearCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta1Left, &Beta1Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindLinearCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		double errorTotal = (nbrLeft * Beta1Left * Beta1Left) + (nbrRight * Beta1Right * Beta1Right);
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindLinearCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindTST2CSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;

	int res = CalculateindTST2CSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta1Left, &Beta1Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindLinearCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		double errorTotal = ((Beta1Left * Beta1Left / nbrLeft) + (Beta1Right * Beta1Right / nbrRight));
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindLinearCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindMaxLCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;

	int res = CalculateindLinearCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta1Left, &Beta1Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindLinearCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		double ert = (Beta1Left > Beta1Right ? Beta1Left : Beta1Right);
		double errorTotal = ert;
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindLinearCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}

//Version optimisé
int CalculateindLinearCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*Beta1Left = 0;
	*Beta1Right = 0;

	double x_bar_left = 0;
	double x_bar_right = 0;
	double y_bar_left = 0;
	double y_bar_right = 0;

	double sumLeftX = 0;
	double sumRightX = 0;
	double sumLeftY = 0;
	double sumRightY = 0;

	double sumLeftX2 = 0;
	double sumLeftXY = 0;
	double sumRightX2 = 0;
	double sumRightXY = 0;


	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*nbrLeft)++;
			sumLeftX += dataG[i];
			sumLeftY += dataT[i];
			sumLeftX2 += dataG[i] * dataG[i];
			sumLeftXY += dataT[i] * dataG[i];
		}
		else
		{
			(*nbrRight)++;
			sumRightX += dataG[i];
			sumRightY += dataT[i];
			sumRightX2 += dataG[i] * dataG[i];
			sumRightXY += dataT[i] * dataG[i];
		}
	}

	x_bar_left = sumLeftX / (*nbrLeft);
	y_bar_left = sumLeftY / (*nbrLeft);

	x_bar_right = sumRightX / (*nbrRight);
	y_bar_right = sumRightY / (*nbrRight);

	if ((*nbrLeft) < incMinObs ||
		(*nbrRight) < incMinObs ||
		(sumLeftX2 - (*nbrLeft) * x_bar_left*x_bar_left) == 0 ||
		(sumRightX2 - (*nbrRight) * x_bar_right*x_bar_right) == 0)
	{
		return 0;
	}

	(*Beta1Left) = (sumLeftXY - (*nbrLeft) * x_bar_left * y_bar_left) / (sumLeftX2 - (*nbrLeft) * x_bar_left*x_bar_left);
	(*Beta1Right) = (sumRightXY - (*nbrRight) * x_bar_right * y_bar_right) / (sumRightX2 - (*nbrRight) * x_bar_right*x_bar_right);

	return 1;
}

//int CalculateindLinearCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right)
//{
// //Remplace par la version optimise.
//	*nbrLeft = 0;
//	*nbrRight = 0;
//	*Beta1Left = 0;
//	*Beta1Right = 0;
//
//	double x_bar_left = 0;
//	double x_bar_right = 0;
//	double y_bar_left = 0;
//	double y_bar_right = 0;
//
//	double sumLeftX = 0;
//	double sumRightX = 0;
//	double sumLeftY = 0;
//	double sumRightY = 0;
//
//	double sumNominatorLeft = 0;
//	double sumDenominatorLeft = 0;
//	double sumNominatorRight = 0;
//	double sumDenominatorRight = 0;
//
//	for (int i = 0; i < NodeNbrObrservation; i++)
//	{
//		if (indCut[i] == LEFT)
//		{
//			(*nbrLeft)++;
//			sumLeftX += dataG[i];
//			sumLeftY += dataT[i];
//		}
//		else
//		{
//			(*nbrRight)++;
//			sumRightX += dataG[i];
//			sumRightY += dataT[i];
//		}
//	}
//
//	x_bar_left = sumLeftX / (*nbrLeft);
//	y_bar_left = sumLeftY / (*nbrLeft);
//
//	x_bar_right = sumRightX / (*nbrRight);
//	y_bar_right = sumRightY / (*nbrRight);
//
//	for (int i = 0; i < NodeNbrObrservation; i++)
//	{
//		if (indCut[i] == LEFT)
//		{
//			sumNominatorLeft += (dataG[i] - x_bar_left) * (dataT[i] - y_bar_left);
//			sumDenominatorLeft += (dataG[i] - x_bar_left) * (dataG[i] - x_bar_left);
//		}
//		else
//		{
//			sumNominatorRight += (dataG[i] - x_bar_right) * (dataT[i] - y_bar_right);
//			sumDenominatorRight += (dataG[i] - x_bar_right) * (dataG[i] - x_bar_right);
//		}
//	}
//
//
//	if ((*nbrLeft) < incMinObs ||
//		(*nbrRight) < incMinObs ||
//		sumDenominatorLeft == 0 ||
//		sumDenominatorRight == 0)
//	{
//		return 0;
//	}
//
//
//	(*Beta1Left) = sumNominatorLeft / sumDenominatorLeft;
//	(*Beta1Right) = sumNominatorRight / sumDenominatorRight;
//
//	return 1;
//}
int CalculateindTSTCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*Beta1Left = 0;
	*Beta1Right = 0;

	double y_bar_left = 0;
	double y_bar_right = 0;

	double sumLeftY = 0;
	double sumRightY = 0;

	double sumErrorLeft = 0;
	double sumErrorRight = 0;
	
	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*nbrLeft)++;
			sumLeftY += dataT[i];
		}
		else
		{
			(*nbrRight)++;
			sumRightY += dataT[i];
		}
	}

	y_bar_left = sumLeftY / (*nbrLeft);
	y_bar_right = sumRightY / (*nbrRight);

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			sumErrorLeft += (dataT[i] - y_bar_left) * (dataT[i] - y_bar_left);
		}
		else
		{
			sumErrorRight += (dataT[i] - y_bar_right) * (dataT[i] - y_bar_right);
		}
	}


	if ((*nbrLeft) < incMinObs ||
		(*nbrRight) < incMinObs)
	{
		return 0;
	}

	(*Beta1Left) = sumErrorLeft;
	(*Beta1Right) = sumErrorRight;

	return 1;
}
//int CalculateindNSRCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right)
//{
//	*nbrLeft = 0;
//	*nbrRight = 0;
//	*Beta1Left = 0;
//	*Beta1Right = 0;
//
//	double x_bar_left = 0;
//	double x_bar_right = 0;
//	double y_bar_left = 0;
//	double y_bar_right = 0;
//
//	double sumLeftX = 0;
//	double sumRightX = 0;
//	double sumLeftY = 0;
//	double sumRightY = 0;
//
//	double sumNominatorLeft = 0;
//	double sumDenominatorLeft = 0;
//	double sumNominatorRight = 0;
//	double sumDenominatorRight = 0;
//
//	for (int i = 0; i < NodeNbrObrservation; i++)
//	{
//		if (indCut[i] == LEFT)
//		{
//			(*nbrLeft)++;
//			sumLeftX += dataG[i];
//			sumLeftY += dataT[i];
//		}
//		else
//		{
//			(*nbrRight)++;
//			sumRightX += dataG[i];
//			sumRightY += dataT[i];
//		}
//	}
//
//	x_bar_left = sumLeftX / (*nbrLeft);
//	y_bar_left = sumLeftY / (*nbrLeft);
//
//	x_bar_right = sumRightX / (*nbrRight);
//	y_bar_right = sumRightY / (*nbrRight);
//
//	for (int i = 0; i < NodeNbrObrservation; i++)
//	{
//		if (indCut[i] == LEFT)
//		{
//			sumNominatorLeft += (dataG[i] - x_bar_left) * (dataT[i] - y_bar_left);
//			sumDenominatorLeft += (dataG[i] - x_bar_left) * (dataG[i] - x_bar_left);
//		}
//		else
//		{
//			sumNominatorRight += (dataG[i] - x_bar_right) * (dataT[i] - y_bar_right);
//			sumDenominatorRight += (dataG[i] - x_bar_right) * (dataG[i] - x_bar_right);
//		}
//	}
//
//
//	if ((*nbrLeft) < incMinObs ||
//		(*nbrRight) < incMinObs ||
//		sumDenominatorLeft == 0 || 
//		sumDenominatorRight == 0)
//	{
//		return 0;
//	}
//
//
//	(*Beta1Left) = sumNominatorLeft / sumDenominatorLeft;
//	(*Beta1Right) = sumNominatorRight / sumDenominatorRight;
//
//	return 1;
//}
int CalculateindTST2CSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta1Left, double * Beta1Right)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*Beta1Left = 0;
	*Beta1Right = 0;

	double sumLeftY = 0;
	double sumRightY = 0;

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*nbrLeft)++;
			sumLeftY += dataT[i];
		}
		else
		{
			(*nbrRight)++;
			sumRightY += dataT[i];
		}
	}

	if ((*nbrLeft) < incMinObs ||
		(*nbrRight) < incMinObs)
	{
		return 0;
	}

	(*Beta1Left) = sumLeftY;
	(*Beta1Right) = sumRightY;

	return 1;
}

int LogSplitInFileindLinearCSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double Beta1Left, double Beta1Right, double errorTotal)
{
	FILE *fp;

	fp = fopen(FileName, "a");
	if (fp == NULL)
	{
		return -1;
	}
	//thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);
	fprintf(fp, "%i;%u;%u;%u;%f;%i;%i;%f;%f;%f\r", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta1Left, Beta1Right, errorTotal);

	fclose(fp);
	return 1;
}

int LogSplitHeaderindLinearCSplit(char * FileName)
{
	FILE *fp;

	fp = fopen(FileName, "w");
	if (fp == NULL)
	{
		return -1;
	}

	fprintf(fp, "threadID;treeID;nodeid;covariate;split;NbrLeftNode;NbrRightNode;ErrorLeft;ErrorRight;TotalError\r");

	fclose(fp);
	return 1;
}
