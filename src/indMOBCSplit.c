#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  "indMOBCSplit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

double getDeltaindMOB0CSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double SumLeft = 0;
	double SumRight = 0;
	double Error2Left = 0;
	double Error2Right = 0;

	int res = CalculateindMOB0CSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &SumLeft, &SumRight, &Error2Left, &Error2Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindMOBLCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		//On veux minimiser
		double errorTotal = -(Error2Left + Error2Right);
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindMOBLCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, SumLeft, SumRight, Error2Left, Error2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindMOBLCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double SumLeft = 0;
	double SumRight = 0;
	double Error2Left = 0;
	double Error2Right = 0;

	int res = CalculateindMOBLCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &SumLeft, &SumRight,&Error2Left,&Error2Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindMOBLCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		//On veux minimiser
		double errorTotal = -(Error2Left + Error2Right);
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindMOBLCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, SumLeft, SumRight, Error2Left, Error2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindMOBQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double SumLeft = 0;
	double SumRight = 0;
	double Error2Left = 0;
	double Error2Right = 0;

	int res = CalculateindMOBQCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &SumLeft, &SumRight, &Error2Left, &Error2Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindMOBLCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		//On veux minimiser
		double errorTotal = -(Error2Left + Error2Right);
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindMOBLCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, SumLeft, SumRight, Error2Left, Error2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}

int CalculateindMOB0CSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * SumLeft, double * SumRight, double * ErrorLeft, double * ErrorRight)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*SumLeft = 0;
	*SumRight = 0;
	*ErrorLeft = 0;
	*ErrorRight = 0;

	double Beta0Left = 0;
	double Beta0Right = 0;

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

	Beta0Left = sumLeftY / (*nbrLeft);
	Beta0Right = sumRightY / (*nbrRight);

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*ErrorLeft) += (dataT[i] - Beta0Left) * (dataT[i] - Beta0Left);
		}
		else
		{
			(*ErrorRight) += (dataT[i] - Beta0Right) * (dataT[i] - Beta0Right);
		}
	}

	(*SumLeft) = sumLeftY;
	(*SumRight) = sumRightY;

	return 1;
}
int CalculateindMOBLCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * SumLeft, double * SumRight, double * ErrorLeft, double * ErrorRight)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*SumLeft = 0;
	*SumRight = 0;
	*ErrorLeft = 0;
	*ErrorRight = 0;

	double Beta0Left = 0;
	double Beta0Right = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;

	double x_bar_left = 0;
	double x_bar_right = 0;
	double y_bar_left = 0;
	double y_bar_right = 0;

	double sumLeftX = 0;
	double sumRightX = 0;
	double sumLeftY = 0;
	double sumRightY = 0;

	double sumNominatorLeft = 0;
	double sumDenominatorLeft = 0;
	double sumNominatorRight = 0;
	double sumDenominatorRight = 0;

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*nbrLeft)++;
			sumLeftX += dataG[i];
			sumLeftY += dataT[i];
		}
		else
		{
			(*nbrRight)++;
			sumRightX += dataG[i];
			sumRightY += dataT[i];
		}
	}

	x_bar_left = sumLeftX / (*nbrLeft);
	y_bar_left = sumLeftY / (*nbrLeft);

	x_bar_right = sumRightX / (*nbrRight);
	y_bar_right = sumRightY / (*nbrRight);

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			sumNominatorLeft += (dataG[i] - x_bar_left) * (dataT[i] - y_bar_left);
			sumDenominatorLeft += (dataG[i] - x_bar_left) * (dataG[i] - x_bar_left);
		}
		else
		{
			sumNominatorRight += (dataG[i] - x_bar_right) * (dataT[i] - y_bar_right);
			sumDenominatorRight += (dataG[i] - x_bar_right) * (dataG[i] - x_bar_right);
		}
	}

	if ((*nbrLeft) < incMinObs ||
		(*nbrRight) < incMinObs ||
		sumDenominatorLeft == 0 ||
		sumDenominatorRight == 0)
	{
		return 0;
	}

	Beta1Left = sumNominatorLeft / sumDenominatorLeft;
	Beta1Right = sumNominatorRight / sumDenominatorRight;

	Beta0Left = y_bar_left - (Beta1Left * x_bar_left);
	Beta0Right = y_bar_right - (Beta1Right * x_bar_right);

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*ErrorLeft) += (dataT[i] - Beta0Left - Beta1Left * dataG[i]) * (dataT[i] - Beta0Left - Beta1Left * dataG[i]);
		}
		else
		{
			(*ErrorRight) += (dataT[i] - Beta0Right - Beta1Right * dataG[i]) * (dataT[i] - Beta0Right - Beta1Right * dataG[i]);
		}
	}

	(*SumLeft) = sumLeftY;
	(*SumRight) = sumRightY;

	return 1;
}
int CalculateindMOBQCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * SumLeft, double * SumRight, double * ErrorLeft, double * ErrorRight)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*SumLeft = 0;
	*SumRight = 0;
	*ErrorLeft = 0;
	*ErrorRight = 0;

	double Beta0Left = 0;
	double Beta0Right = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;
	double Beta2Left = 0;
	double Beta2Right = 0;

	double x_bar_Left = 0;
	double x_bar_Right = 0;
	double y_bar_Left = 0;
	double y_bar_Right = 0;

	double sumLeftX = 0;
	double sumRightX = 0;
	double sumLeftY = 0;
	double sumRightY = 0;

	double sumLeftX2 = 0;
	double sumRightX2 = 0;
	double sumLeftXX2 = 0;
	double sumRightXX2 = 0;
	double sumLeftX22 = 0;
	double sumRightX22 = 0;
	double sumLeftXY = 0;
	double sumRightXY = 0;
	double sumLeftYX2 = 0;
	double sumRightYX2 = 0;
	double x2_bar_Left = 0;
	double x2_bar_Right = 0;

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*nbrLeft)++;
			sumLeftX += dataG[i];
			sumLeftY += dataT[i];

			sumLeftX2 += dataG[i] * dataG[i];
			sumLeftXX2 += dataG[i] * dataG[i] * dataG[i];
			sumLeftX22 += dataG[i] * dataG[i] * dataG[i] * dataG[i];

			sumLeftXY += dataG[i] * dataT[i];
			sumLeftYX2 += dataG[i] * dataG[i] * dataT[i];

		}
		else
		{
			(*nbrRight)++;
			sumRightX += dataG[i];
			sumRightY += dataT[i];

			sumRightX2 += dataG[i] * dataG[i];
			sumRightXX2 += dataG[i] * dataG[i] * dataG[i];
			sumRightX22 += dataG[i] * dataG[i] * dataG[i] * dataG[i];

			sumRightXY += dataG[i] * dataT[i];
			sumRightYX2 += dataG[i] * dataG[i] * dataT[i];
		}
	}

	x_bar_Left = sumLeftX / *nbrLeft;
	x2_bar_Left = sumLeftX2 / *nbrLeft;
	y_bar_Left = sumLeftY / *nbrLeft;

	x_bar_Right = sumRightX / *nbrRight;
	x2_bar_Right = sumRightX2 / *nbrRight;
	y_bar_Right = sumRightY / *nbrRight;

	double s11Left = sumLeftX2 - (sumLeftX * sumLeftX / *nbrLeft);
	double s12Left = sumLeftXX2 - (sumLeftX * sumLeftX2 / *nbrLeft);
	double s22Left = sumLeftX22 - (sumLeftX2 * sumLeftX2 / *nbrLeft);
	double sy1Left = sumLeftXY - (sumLeftY * sumLeftX / *nbrLeft);
	double sy2Left = sumLeftYX2 - (sumLeftY * sumLeftX2 / *nbrLeft);

	double s11Right = sumRightX2 - (sumRightX * sumRightX / *nbrRight);
	double s12Right = sumRightXX2 - (sumRightX * sumRightX2 / *nbrRight);
	double s22Right = sumRightX22 - (sumRightX2 * sumRightX2 / *nbrRight);
	double sy1Right = sumRightXY - (sumRightY * sumRightX / *nbrRight);
	double sy2Right = sumRightYX2 - (sumRightY * sumRightX2 / *nbrRight);


	if ((*nbrLeft) < incMinObs ||
		(*nbrRight) < incMinObs ||
		(s22Left * s11Left - s12Left * s12Left) == 0 ||
		(s22Right * s11Right - s12Right * s12Right) == 0)
	{
		return 0;
	}

	Beta1Left = (sy1Left * s22Left - sy2Left * s12Left) / (s22Left * s11Left - s12Left * s12Left);
	Beta2Left = (sy2Left * s11Left - sy1Left * s12Left) / (s22Left * s11Left - s12Left * s12Left);
	Beta0Left = y_bar_Left - Beta1Left * x_bar_Left - Beta2Left * x2_bar_Left;

	Beta1Right = (sy1Right * s22Right - sy2Right * s12Right) / (s22Right * s11Right - s12Right * s12Right);
	Beta2Right = (sy2Right * s11Right - sy1Right * s12Right) / (s22Right * s11Right - s12Right * s12Right);
	Beta0Right = y_bar_Right - Beta1Right * x_bar_Right - Beta2Right * x2_bar_Right;

	for (int i = 0; i < NodeNbrObrservation; i++)
	{
		if (indCut[i] == LEFT)
		{
			(*ErrorLeft) += (dataT[i] - Beta0Left - Beta1Left * dataG[i] - Beta2Left * dataG[i] * dataG[i]) * (dataT[i] - Beta0Left - Beta1Left * dataG[i] - Beta2Left * dataG[i] * dataG[i]);
		}
		else
		{
			(*ErrorRight) += (dataT[i] - Beta0Right - Beta1Right * dataG[i] - Beta2Right * dataG[i] * dataG[i]) * (dataT[i] - Beta0Right - Beta1Right * dataG[i] - Beta2Right * dataG[i] * dataG[i]);
		}
	}

	(*SumLeft) = sumLeftY;
	(*SumRight) = sumRightY;

	return 1;
}

int LogSplitInFileindMOBLCSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double SumLeft, double SumRight, double ErrorLeft, double ErrorRight, double errorTotal)
{
	FILE *fp;

	fp = fopen(FileName, "a");
	if (fp == NULL)
	{
		return -1;
	}
	fprintf(fp, "%i;%u;%u;%u;%f;%i;%i;%f;%f;%f;%f;%f\r", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, SumLeft, SumRight, ErrorLeft, ErrorRight, errorTotal);

	fclose(fp);
	return 1;
}

int LogSplitHeaderindMOBLCSplit(char * FileName)
{
	FILE *fp;

	fp = fopen(FileName, "w");
	if (fp == NULL)
	{
		return -1;
	}

	fprintf(fp, "threadID;treeID;nodeid;covariate;split;NbrLeftNode;NbrRightNode;SumYLeft;SumYRight;ErrorLeft;ErrorRight;TotalError\r");

	fclose(fp);
	return 1;
}
