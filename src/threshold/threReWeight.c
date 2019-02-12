/**CFile****************************************************************
 
  FileName    [threReWeight.c] 

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [threshold.]
  
  Synopsis    [Threshold network approximation.]

  Author      [Hao Chen]
   
  Affiliation [NTU]

  Date        [Jan 13, 2018.]

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

void Th_NtkReWeight(Vec_Ptr_t* vThres);
void Th_ObjReWeight(Thre_S* tObj);

/**Function*************************************************************

  Synopsis    [ReWeight 1-DL threshold gates to optimal weights and threshold]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Th_NtkReWeight(Vec_Ptr_t* vThres) {
    int i;
    Thre_S* tObj;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
        if (tObj->Type == Th_CONST1 || tObj->Type == Th_Pi || tObj->Type == Th_Po)
            continue;
        if (Vec_IntFind(tObj->dtypes, 2) != -1)  {
            // printf("Not 1-DL\n");
            continue;
        }
        Th_ObjReWeight(tObj);
    }
}

/**Function*************************************************************

  Synopsis    [ReWeight a 1-DL threshold gate to optimal weights and threshold]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Th_ObjReWeight(Thre_S* tObj) {
    int i, id;
    int sumW = 0;
    int oldThre = 0;
    int newThre = 1;
    Vec_IntForEachEntryReverse(tObj->Fanins, id, i) {
        if (i == Vec_IntSize(tObj->Fanins) - 1) {
            Vec_IntWriteEntry(tObj->weights, i, 1);
            sumW = 1;
            oldThre = 1;
        }
        else if (Vec_IntEntry(tObj->dtypes, i) == 1) {
            newThre = oldThre;
            Vec_IntWriteEntry(tObj->weights, i, oldThre);
            sumW += oldThre;
        }
        else if (Vec_IntEntry(tObj->dtypes, i) == 0) {
            newThre = sumW + 1;
            Vec_IntWriteEntry(tObj->weights, i, newThre - oldThre);
            sumW += (newThre - oldThre);
            oldThre = newThre;
        }
        else assert(0);
    }
    tObj->thre = newThre;
}

////////////////////////////////////////////////////////////////////////
///                       END OF FILE                                ///
////////////////////////////////////////////////////////////////////////
