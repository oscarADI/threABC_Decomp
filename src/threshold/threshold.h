/**HFile****************************************************************

  FileName    [threshold.h]

  SystemName  [ABC]

  PackageName [threshold]

  Synopsis    [Header for threshold gate]

  Author      [Nian-Ze Lee]
  
  Affiliation [NTU]

  Date        [May 11, 2016]

***********************************************************************/

#ifndef THRESHOLD_H_
#define THRESHOLD_H_

////////////////////////////////////////////////////////////////////////
///                          INCLUDES                                ///
////////////////////////////////////////////////////////////////////////

#include "base/abc/abc.h"
#include "stdio.h"
#include <math.h>

////////////////////////////////////////////////////////////////////////
///                         PARAMETERS                               ///
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
///                         BASIC TYPES                              ///
////////////////////////////////////////////////////////////////////////

typedef enum {
	Th_Unknown = 0,
	Th_Pi,
	Th_Po,
	Th_Node,
	Th_CONST1 // const1 gate always has Id = 0
}Th_Gate_Type;

typedef struct Thre_S_       Thre_S;   
typedef struct Pair_S_       Pair_S;
typedef struct Th_Stat_      Th_Stat;
typedef struct Th_Set_       Th_Set;

struct Thre_S_
{
	int         thre;
  int         Type;
  int         Ischoose;
  int         Id;
	int         oId;
  int         nId;
	int         cost;
	int         level;
  char      * pName;
  Vec_Int_t * weights;
  Vec_Int_t * Fanins;
  Vec_Int_t * Fanouts;
  Abc_Obj_t * pCopy;
  Vec_Int_t * dtypes;    // 1 for 1->1 , 0 for 0->0 
  Vec_Int_t * thres;
  Vec_Int_t * FaninCs;
};

struct Pair_S_
{
	int IntK;
  int IntL;
};

struct Th_Stat_
{
  int numTotal;
	int numSingleFoutIf;
	int numSingleFoutIff;
	int numSingleFoutFail;
	int numMultiFout;
	int numMultiFoutOk;
	int numTwoFout;
	int numTwoFoutOk;
	int numNotThNode;
  // threshold --> mux redundancy check
  int numRedundancy;
  int redund[50];
};

// for extraction
struct Th_Set_
{
  Vec_Int_t * set;
  int Dtype; // 1 for 1->1 , 0 for 0->0 
  int gate;
};

////////////////////////////////////////////////////////////////////////
///                         GLOABAL VAR                              ///
////////////////////////////////////////////////////////////////////////

Vec_Ptr_t * current_TList;
//Vec_Ptr_t * another_TList;
Vec_Ptr_t * cut_TList;
int         globalRef; // for iterative collapse
Th_Stat     thProfiler;

// for approx
Vec_Ptr_t * Golden_Vec;
Vec_Int_t * Critical;
int nPi;
int nPo;
int nThres;

////////////////////////////////////////////////////////////////////////
///                    FUNCTION DECLARATIONS                         ///
////////////////////////////////////////////////////////////////////////

//===threFileIo.c================================//

extern void         DeleteTList(Vec_Ptr_t*);
extern Vec_Ptr_t *  func_readFileOAO(char*);
extern Thre_S    *  newThre(int , char*);
extern int          getIdByName(Vec_Ptr_t* , char*);
extern Thre_S    *  getPtrByName(Vec_Ptr_t* , char*);
extern void         printNodeStats(Vec_Ptr_t*);
extern void         dumpTh2FileNZ(Vec_Ptr_t * , char*);

//===threEC.c===================================//

extern void        func_EC_writePB(Abc_Ntk_t *, Vec_Ptr_t *, char*);
extern Vec_Ptr_t*  thre_PB(FILE*, Vec_Ptr_t *);
extern Vec_Ptr_t*  alan_PB(FILE*, Abc_Ntk_t *);
extern void        miter_PB(FILE*, Vec_Ptr_t*, Vec_Ptr_t *);
extern int         Max_Thre(Thre_S *);
extern int         min_Thre(Thre_S *);

//===threAig2Th.c===================================//

extern void aigThCreateConst1( Abc_Ntk_t * , Vec_Ptr_t * , Vec_Int_t * );
extern void aigThCreateCi( Abc_Ntk_t * , Vec_Ptr_t * , Vec_Int_t *  );
extern void aigThCreateCo( Abc_Ntk_t * , Vec_Ptr_t * , Vec_Int_t * );
extern void aigThCreateNode( Abc_Ntk_t * , Vec_Ptr_t * , Vec_Int_t * );
extern void aigThConnectFanin( Abc_Ntk_t * , Vec_Ptr_t * , Vec_Int_t * );
extern void aigThConnectFanout( Abc_Ntk_t * , Vec_Ptr_t *  , Vec_Int_t * );
extern Thre_S*    Th_CreateObj( Vec_Ptr_t * , Th_Gate_Type );
extern void       Th_DeleteObj( Thre_S * );
extern Vec_Ptr_t* aig2Th( Abc_Ntk_t * ); 

//===threKLMerge.c===================================//

extern Thre_S*    Th_GetObjById          ( Vec_Ptr_t * , int );
extern Thre_S*    Th_MergeNodes          ( Thre_S *    , int );
extern Thre_S*    Th_CalKLMerge          ( const Thre_S * , const Thre_S * , int );
extern Pair_S*    Th_CalKL               ( const Thre_S * , const Thre_S * , int , int , int );
extern Thre_S*    Th_KLMerge             ( const Thre_S * , const Thre_S * , const Pair_S * , int , int );
extern Thre_S*    Th_InvertObj           ( const Thre_S* );
extern Thre_S*    Th_CopyObj             ( const Thre_S* );
extern int        MaxF                   ( Vec_Int_t * , int );
extern int        MinF                   ( Vec_Int_t * , int );
extern Thre_S*    Th_CreateObjNoInsert   ( Th_Gate_Type );
extern void       Th_DeleteObjNoInsert   ( Thre_S * );
extern void       Th_DumpObj             ( const Thre_S * );
extern void       Th_DumpMergeObj        ( const Thre_S * , const Thre_S * , const Thre_S *);
extern void       mergeThreNtk_Iter      ( Vec_Ptr_t * , int );
extern void       deleteNode             ( Vec_Ptr_t * , Thre_S * );

//===threStat.c===================================//

extern void       Th_PrintStat          ( Vec_Ptr_t * );
extern int        Th_CountGate          ( Vec_Ptr_t * , Th_Gate_Type );
extern int        Th_CountLevel         ( Vec_Ptr_t * );
extern void       Th_SetLevel           ( Vec_Ptr_t * );
extern void       Th_ComputeLevel       ( Vec_Ptr_t * );
extern void       Th_ComputeLevel_rec   ( Vec_Ptr_t * , Thre_S * );
extern int        Th_FindMaxLevel       ( Vec_Ptr_t * );
extern Vec_Int_t* Th_CountFanio         ( Vec_Ptr_t * , int );
extern void       Th_PrintData          ( int , int , int , int , Vec_Int_t * , Vec_Int_t * );
extern void       Th_PrintNode          ( Thre_S *);

//===threProfile.c================================//

extern void       Th_ProfileInit        ();
extern void       Th_ProfilePrint       ();

//===threCNF.c================================//

extern void       func_EC_writeCNF      ( Abc_Ntk_t * , Vec_Ptr_t * , char * );

//===thresholdNZ.c================================//

extern void       func_EC_compareTH     ( Vec_Ptr_t * , Vec_Ptr_t * );
extern void       func_CNF_compareTH    ( Vec_Ptr_t * , Vec_Ptr_t * );

//===threKLCollapse.c================================//

extern void       Th_CollapseNtk         ( Vec_Ptr_t * , int , int );
extern void       Th_NtkDfs              ();

//===threTh2Blif.c===================================//

extern void       Th_WriteBlif           ( Vec_Ptr_t * thre_list , const char * );

//===threTh2Mux.c====================================//

extern Abc_Ntk_t* Th_Ntk2Mux             ( Vec_Ptr_t * , int , int );

//=================== New Added ====================//


//===threTransWeight.c=============================//

extern void   Th_NtkTransWeight(Vec_Ptr_t *);
extern void   Th_NtkTransBack(Vec_Ptr_t* );
extern void   Th_ObjSortWeight(Thre_S* );

//===threTh2Dlist.c================================//
extern void   Th_NtkTh2DList(Vec_Ptr_t *);
extern int    Th2Dlist(Thre_S *);
void          Th_DLBuildSets(Vec_Ptr_t * vThres);
extern int    compare(Th_Set** i1, Th_Set** i2);
extern void   Th_Decompose(Vec_Ptr_t *, Thre_S *);
//===threCost.c====================================//
extern double Th_NtkCost(Vec_Ptr_t *);
extern double Th_NtkCost2(Vec_Ptr_t *);
extern double Th_NtkCost3(Vec_Ptr_t *);
extern double Th_NtkCostThre(Vec_Ptr_t *);

//===threReWeight.c================================//
extern void   Th_NtkReWeight(Vec_Ptr_t *);
extern void   Th_ObjReWeight(Thre_S *);

//===threDLExtract.c================================//
extern  void      Th_IterativeUpdate(Vec_Ptr_t *, Vec_Ptr_t *);
extern  int       FindMaxSet(Vec_Ptr_t * , Vec_Ptr_t * , int , int);
extern  void      Th_ObjUpdate(Vec_Ptr_t * , Thre_S *);
extern  void      Th_RemoveDummyNode(Vec_Ptr_t *, Thre_S *);
extern  Thre_S *  Th_CreateTempObj(Vec_Ptr_t *, Th_Gate_Type ); 
extern  void      Th_PrintSet();
extern  int       Th_FindCritical(int);
extern  void      Th_FindCritical_rec(Thre_S* tObj, int level);
extern  int       Th_CountLevelWithDummy( Vec_Ptr_t * );
extern  void      Th_ComputeLevelWithDummy( Vec_Ptr_t * );
extern  void      Th_ComputeLevelWithDummy_rec( Vec_Ptr_t *  , Thre_S *  );
extern  void      Th_RemoveDummy(); 
extern  void      Th_IterativeUpdatePQ(Vec_Ptr_t* vThres, Vec_Ptr_t* vSet, int levelbound);



#endif

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////
