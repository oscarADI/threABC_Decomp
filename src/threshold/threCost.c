/**CFile****************************************************************
 
  FileName    [threApprox.c] 

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [threshold.]
  
  Synopsis    [Threshold network approximation.]

  Author      [Hao Chen]
   
  Affiliation [NTU]

  Date        [Jan 11, 2018.]

  Revision    [$Id: abc.c,v 1.00 2005/06/20 00:00:00 alanmi Exp $]

***********************************************************************/

////////////////////////////////////////////////////////////////////////
///                          INCLUDES                                ///
////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include "base/abc/abc.h"
#include "map/if/if.h"
#include "map/if/ifCount.h"
#include "bdd/extrab/extraBdd.h"
#include "misc/extra/extra.h"
#include "threshold.h"

////////////////////////////////////////////////////////////////////////
///                        DECLARATIONS                              ///
////////////////////////////////////////////////////////////////////////

double Th_NtkCost(Vec_Ptr_t *);
double Th_NtkCost2(Vec_Ptr_t *); // Weight
double Th_NtkCost3(Vec_Ptr_t *);
double Th_NtkCostThre(Vec_Ptr_t *);
int Th_ObjCost(Thre_S *);
int Th_ObjCost2(Thre_S *);
int Th_ObjCost3(Thre_S *);

/**Function*************************************************************

  Synopsis    []

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/

double Th_NtkCost(Vec_Ptr_t * vThres) {
    int i;
    Thre_S* tObj;
    int cost = 0;
    int nTh = 0;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
      if (tObj->Type != Th_Node) continue;
      cost += Th_ObjCost(tObj);
      nTh += 1;
    }
    return cost;
    // return 0.5 * cost + 0.5 * nTh;
}
double Th_NtkCostThre(Vec_Ptr_t * vThres) {
  int cost = 0;
  Thre_S* tObj;
  int i;
  Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
    if (tObj->Type != Th_Node) continue;
    cost += tObj->thre;
  }
  return cost;
}

int Th_ObjCost(Thre_S * tObj) {
    return Vec_IntSize(tObj->Fanins);
}

double Th_NtkCost2(Vec_Ptr_t * vThres) {
    int i;
    Thre_S* tObj;
    int cost = 0;
    int nTh = 0;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
      if (tObj->Type != Th_Node) continue;
      cost += Th_ObjCost2(tObj);
    }
    return cost;
}
double Th_NtkCost3(Vec_Ptr_t * vThres) {
    int i;
    int cost = 0;
    Thre_S* tObj;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
      if (tObj->Type != Th_Node) continue;
      cost += Th_ObjCost3(tObj);
    }
    return cost;
}

int Th_ObjCost2(Thre_S * tObj) {
    int cost = 0;
    int i, w;
    Vec_IntForEachEntry(tObj->weights, w, i) {
      cost += w;
    }
    // cost += tObj->thre;
    return cost;
}

int Th_ObjCost3(Thre_S * tObj) {
  int s = Vec_IntSize(tObj->Fanins);
  if (s >= 2) return 12+8*(s-2);
  else return 0;
 /* switch (s) {
    case 2: return 3;
    case 3: return 5;
    case 4: return 5;
    case 5: return 7;
    case 6: return 8;
    case 7: return 9;
    case 8: return 14;
    case 9: return 14;
    case 10: return 15;
    case 11: return 19;
    case 12: return 20;
    case 13: return 21;
    case 14: return 21;
    case 15: return 27;
    



    ddouble Th_NtkCostThre(Vec_Ptr_t *)
     double Th_NtkCostThre(Vec_Ptr_t *)ins : %d\n", s);
     double Th_NtkCostThre(Vec_Ptr_t *)
    }double Th_NtkCostThre(Vec_Ptr_t *)
  }*/
}

