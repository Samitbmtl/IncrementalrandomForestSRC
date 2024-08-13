#include <stdio.h>
#include <stdlib.h>
#include  <math.h>
#include  "indQuadCSplit.h"
#ifdef _OPENMP
#include <omp.h>
#endif

double getDeltaindQuadCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta0Left = 0;
	double Beta0Right = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;
	double Beta2Left = 0;
	double Beta2Right = 0;

	int res = CalculateindQuadCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta0Left, &Beta0Right, &Beta1Left, &Beta1Right, &Beta2Left, &Beta2Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindQuadCSplit("SplitVal.csv");
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
				LogSplitInFileindQuadCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta0Left, Beta0Right, Beta1Left, Beta1Right, Beta2Left, Beta2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindQuad2CSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta0Left = 0;
	double Beta0Right = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;
	double Beta2Left = 0;
	double Beta2Right = 0;


	double errorLeft = 0;
	double errorRight = 0;

	int res = CalculateindQuadCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta0Left, &Beta0Right, &Beta1Left, &Beta1Right, &Beta2Left, &Beta2Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindQuadCSplit("SplitVal.csv");
	}
	if (res == 1 && Beta2Left != 0 && Beta2Right != 0)
	{
//#let m = -b1 / (2 * b2)
//#if(b2<0 and m<=1 and m>=0) then the maximum is b0+b1*m+b2*m^2
//#else
//#the maximum is max(b0, b0 + b1 + b2)

		errorLeft = -Beta1Left / (2 * Beta2Left);
		if (Beta2Left < 0 && errorLeft <= 1 && errorLeft >= 0)
		{
			errorLeft = Beta0Left + (Beta1Left * errorLeft)  + (Beta2Left * errorLeft*errorLeft);
		}
		else
		{
			errorLeft = (Beta0Left >= Beta0Left + Beta1Left + Beta2Left ? Beta0Left : Beta0Left + Beta1Left + Beta2Left);
		}
		errorRight = -Beta1Right / (2 * Beta2Right);
		if (Beta2Right < 0 && errorRight <= 1 && errorRight >= 0)
		{
			errorRight = Beta0Right + (Beta1Right * errorRight) + (Beta2Right * errorRight * errorRight);
		}
		else
		{
			errorRight = (Beta0Right >= Beta0Right + Beta1Right + Beta2Right ? Beta0Right : Beta0Right + Beta1Right + Beta2Right);
		}
	
		double ert = (errorLeft - errorRight);
		ert = (ert < 0 ? -ert : ert);
		double errorTotal = sqrt(nbrLeft * nbrRight) * ert;
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindQuadCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta0Left, Beta0Right, Beta1Left, Beta1Right, Beta2Left, Beta2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}

// double getDeltaindMaxQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
// {
	// int nbrLeft = 0;
	// int nbrRight = 0;
	// double Beta0Left = 0;
	// double Beta0Right = 0;
	// double Beta1Left = 0;
	// double Beta1Right = 0;
	// double Beta2Left = 0;
	// double Beta2Right = 0;

	// int res = CalculateindQuadCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta0Left, &Beta0Right, &Beta1Left, &Beta1Right, &Beta2Left, &Beta2Right);

	// if (log && writeLogHeader)
	// {
		// LogSplitHeaderindQuadCSplit("SplitVal.csv");
	// }
	// if (res == 1)
	// {
		// double ert = (Beta1Left > Beta1Right ? Beta1Left : Beta1Right);
		// double errorTotal = ert;
		// if (log)
		// {
// #ifdef _OPENMP
// #pragma omp critical (writelog)
// #endif
			// {
				// LogSplitInFileindQuadCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta0Left, Beta0Right, Beta1Left, Beta1Right, Beta2Left, Beta2Right, errorTotal);
			// }
		// }
		// return errorTotal;
	// }
	// {
		// return 0;
	// }
// }

//Prise 2 : apres changement formule max effect pour article.
double getDeltaindMaxQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta0Left = 0;
	double Beta0Right = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;
	double Beta2Left = 0;
	double Beta2Right = 0;
	double tauLeft = 0;
	double tauRight = 0;

	int res = CalculateindQuadCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta0Left, &Beta0Right, &Beta1Left, &Beta1Right, &Beta2Left, &Beta2Right);
	if (Beta2Left < 0 && (-Beta1Left / (2 * Beta2Left)) >= 0 && (-Beta1Left / (2 * Beta2Left)) <= 1)
	{
		tauLeft = Beta0Left - ((Beta1Left*Beta1Left) / (4 * Beta2Left));
	}
	else
	{
		if (Beta0Left > (Beta0Left + Beta1Left + Beta2Left))
		{
			tauLeft = Beta0Left;
		}
		else
		{
			tauLeft = (Beta0Left + Beta1Left + Beta2Left);
		}
	}
	
	if (Beta2Right < 0 && (-Beta1Right / (2 * Beta2Right)) >= 0 && (-Beta1Right / (2 * Beta2Right)) <= 1)
	{
		tauRight = Beta0Right - ((Beta1Right*Beta1Right) / (4 * Beta2Right));
	}
	else
	{
		if (Beta0Right > (Beta0Right + Beta1Right + Beta2Right))
		{
			tauRight = Beta0Right;
		}
		else
		{
			tauRight = (Beta0Right + Beta1Right + Beta2Right);
		}
	}
	
	if (log && writeLogHeader)
	{
		LogSplitHeaderindQuadCSplit("SplitVal.csv");
	}
	if (res == 1)
	{
		double ert = (tauLeft - tauRight);
		ert = (ert < 0 ? -ert : ert);
		double errorTotal = sqrt(nbrLeft * nbrRight) * ert;
		
		//double ert = (Beta1Left > Beta1Right ? Beta1Left : Beta1Right);
		//double errorTotal = ert;
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindQuadCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta0Left, Beta0Right, Beta1Left, Beta1Right, Beta2Left, Beta2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}
double getDeltaindMaxMaxQCSplit(int thread_id, uint treeID, uint nodeid, uint covariate, double split, uint log, int writeLogHeader, double * userTime, double * userGroup, char * userSplitIndicator, int NodeNbrObrservation, int incMinObs)
{
	int nbrLeft = 0;
	int nbrRight = 0;
	double Beta0Left = 0;
	double Beta0Right = 0;
	double Beta1Left = 0;
	double Beta1Right = 0;
	double Beta2Left = 0;
	double Beta2Right = 0;

	double errorLeft = 0;
	double errorRight = 0;

	int res = CalculateindQuadCSplit(userTime, userGroup, userSplitIndicator, NodeNbrObrservation, incMinObs, &nbrLeft, &nbrRight, &Beta0Left, &Beta0Right, &Beta1Left, &Beta1Right, &Beta2Left, &Beta2Right);

	if (log && writeLogHeader)
	{
		LogSplitHeaderindQuadCSplit("SplitVal.csv");
	}
	if (res == 1 && Beta2Left != 0 && Beta2Right != 0)
	{
		//#let m = -b1 / (2 * b2)
		//#if(b2<0 and m<=1 and m>=0) then the maximum is b0+b1*m+b2*m^2
		//#else
		//#the maximum is max(b0, b0 + b1 + b2)

		errorLeft = -Beta1Left / (2 * Beta2Left);
		if (Beta2Left < 0 && errorLeft <= 1 && errorLeft >= 0)
		{
			errorLeft = Beta0Left + (Beta1Left * errorLeft) + (Beta2Left * errorLeft*errorLeft);
		}
		else
		{
			errorLeft = (Beta0Left >= Beta0Left + Beta1Left + Beta2Left ? Beta0Left : Beta0Left + Beta1Left + Beta2Left);
		}
		errorRight = -Beta1Right / (2 * Beta2Right);
		if (Beta2Right < 0 && errorRight <= 1 && errorRight >= 0)
		{
			errorRight = Beta0Right + (Beta1Right * errorRight) + (Beta2Right * errorRight * errorRight);
		}
		else
		{
			errorRight = (Beta0Right >= Beta0Right + Beta1Right + Beta2Right ? Beta0Right : Beta0Right + Beta1Right + Beta2Right);
		}

		double ert = (errorLeft > errorRight ? errorLeft : errorRight);
		double errorTotal = ert;
		if (log)
		{
#ifdef _OPENMP
#pragma omp critical (writelog)
#endif
			{
				LogSplitInFileindQuadCSplit("SplitVal.csv", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta0Left, Beta0Right, Beta1Left, Beta1Right, Beta2Left, Beta2Right, errorTotal);
			}
		}
		return errorTotal;
	}
	{
		return 0;
	}
}

int CalculateindQuadCSplit(double * dataT, double * dataG, char * indCut, int NodeNbrObrservation, int incMinObs, int * nbrLeft, int * nbrRight, double * Beta0Left, double * Beta0Right, double * Beta1Left, double * Beta1Right, double * Beta2Left, double * Beta2Right)
{
	*nbrLeft = 0;
	*nbrRight = 0;
	*Beta0Left = 0;
	*Beta0Right = 0;
	*Beta1Left = 0;
	*Beta1Right = 0;
	*Beta2Left = 0;
	*Beta2Right = 0;

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


	*Beta1Left = (sy1Left * s22Left - sy2Left * s12Left) / (s22Left * s11Left - s12Left * s12Left);
	*Beta2Left = (sy2Left * s11Left - sy1Left * s12Left) / (s22Left * s11Left - s12Left * s12Left);
	*Beta0Left = y_bar_Left - *Beta1Left * x_bar_Left - *Beta2Left * x2_bar_Left;

	*Beta1Right = (sy1Right * s22Right - sy2Right * s12Right) / (s22Right * s11Right - s12Right * s12Right);
	*Beta2Right = (sy2Right * s11Right - sy1Right * s12Right) / (s22Right * s11Right - s12Right * s12Right);
	*Beta0Right = y_bar_Right - *Beta1Right * x_bar_Right - *Beta2Right * x2_bar_Right;

	return 1;
}

int LogSplitInFileindQuadCSplit(char * FileName, int thread_id, uint treeID, uint nodeid, uint covariate, double split, int nbrLeft, int nbrRight, double Beta0Left, double Beta0Right, double Beta1Left, double Beta1Right, double Beta2Left, double Beta2Right, double errorTotal)
{
	FILE *fp;

	fp = fopen(FileName, "a");
	if (fp == NULL)
	{
		return -1;
	}
	fprintf(fp, "%i;%u;%u;%u;%f;%i;%i;%f;%f;%f;%f;%f;%f;%f\r", thread_id, treeID, nodeid, covariate, split, nbrLeft, nbrRight, Beta0Left, Beta0Right, Beta1Left, Beta1Right, Beta2Left, Beta2Right, errorTotal);

	fclose(fp);
	return 1;
}
int LogSplitHeaderindQuadCSplit(char * FileName)
{
	FILE *fp;

	fp = fopen(FileName, "w");
	if (fp == NULL)
	{
		return -1;
	}

	fprintf(fp, "threadID;treeID;nodeid;covariate;split;NbrLeftNode;NbrRightNode;Beta0Left; Beta0Right; Beta1Left; Beta1Right; Beta2Left; Beta2Right;TotalError\r");

	fclose(fp);
	return 1;
}
