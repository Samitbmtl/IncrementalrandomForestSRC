//#define LEFT      0x01
//#define RIGHT     0x00
//#define TREATMENT    1
//#define CONTROL      0
//#define EVENT        1
//#define NOEVENT      0

#ifndef LEFT
#define LEFT      0x01
#endif

#ifndef RIGHT
#define RIGHT     0x00
#endif

#ifndef TREATMENT
#define TREATMENT    1
#endif

#ifndef CONTROL
#define CONTROL      0
#endif

#ifndef EVENT
#define EVENT        1
#endif

#ifndef NOEVENT
#define NOEVENT      0
#endif

#ifndef uint
typedef unsigned int  uint;
#endif

double getDelta(int thread_id,uint treeID,uint nodeid,uint covariate,double split,uint log, int writeLogHeader, double * userTime,double * userEvent,double * userGroup,char * userSplitIndicator,int NodeNbrObrservation,int incTailType, int incMinObs, int incMinObsControl,int incMinObsTreatment,double incTreatmentLowBound,double incTreatmentUpperBound);
int SortedCalculateSplit(double * dataT, double * dataC, double * dataG, char * indCut, int NodeNbrObrservation,int incTailType, int incMinObs, int incMinObsControl,int incMinObsTreatment,double incTreatmentLowBound,double incTreatmentUpperBound, int * p_nbrLeft, int * p_nbrLeftTreatment, int * p_nbrLeftControl, int * p_nbrRight, int * p_nbrRightTreatment, int * p_nbrRightControl, double * p_errorTotalLeft, double * p_errorTotalRight, int * SNbrLeft_Control, int * SNbrLeft_Traitement, int * SNbrRight_Control, int * SNbrRight_Traitement, double * CumulTLeft, double * CumulTRight, double * MaxT_Left, double * UsedMaxT_Left, double * MaxT_Right, double * UsedMaxT_Right, double * CumulDiffBrutLeft, double * CumulDiffBrutRight, int * p_nbrEventLeftTreatment, int * p_nbrEventLeftControl, int * p_nbrEventRightTreatment, int * p_nbrEventRightControl);
int NewCompare2KaplanMaier(double * TLeft, double * SLeft, int SNbrLeft, int NbrLeft, double * TRight, double * SRight, int SNbrRight, int NbrRight,int Tail, double * CumulDiff, double * CumulDiffBrut, double * CumulT, double * MaxT, double * UsedMaxT, double * IntegralLeft, double * IntegralRight);
int SortedKaplanMaier(double * tTra, double * cTra, double * gTra, int IsTreatmentGroup, int nbr, char * indObs, char BelowCut, double ** T, double ** S, int * SNbr);
int LogSplitInFile(char * FileName,int thread_id, uint treeID,uint nodeid,uint covariate,double split, int nbrLeft, int nbrLeftTreatment, int nbrLeftControl, int nbrRight, int nbrRightTreatment, int nbrRightControl, double errorLeft, double errorRight, double error, double error1, double error2, int SNbrLeft_Control, int SNbrLeft_Traitement, int SNbrRight_Control, int SNbrRight_Traitement, int CumulTLeft, int CumulTRight, int MaxT_Left, int UsedMaxT_Left, int MaxT_Right, int UsedMaxT_Right, double CumulDiffBrutLeft, double CumulDiffBrutRight, int nbrEventLeftTreatment, int nbrEventLeftControl, int nbrEventRightTreatment, int nbrEventRightControl);
int LogSplitHeader(char * FileName);
//void SortRX(int n, double * dataX, int * index);