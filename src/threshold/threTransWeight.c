/**CFile****************************************************************
 
  FileName    [threTransWeight.c] 

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [threshold.]
  
  Synopsis    [Threshold network approximation.]

  Author      [Hao Chen]
   
  Affiliation [NTU]

  Date        [Jan 9, 2018.]

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

void Th_NtkTransWeight(Vec_Ptr_t* vThres);
void Th_NtkTransBack(Vec_Ptr_t* vThres);
void Th_ObjSortWeight(Thre_S* tObj);

/**Function*************************************************************

  Synopsis    [transform threshold gate with negative weight to postive]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/

void Th_NtkTransWeight(Vec_Ptr_t* vThres) {
    int i,j;
    Thre_S* tObj;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
        if (tObj->Type == Th_CONST1 || tObj->Type == Th_Pi) continue;
        tObj->FaninCs = Vec_IntStart(Vec_IntSize(tObj->Fanins));
        int weight;
        int sum = 0;
        Vec_IntForEachEntry(tObj->weights, weight, j) {
            if (weight < 0) {
                Vec_IntSetEntry(tObj->weights, j, -weight);
                Vec_IntSetEntry(tObj->FaninCs, j, 1);
                sum += -weight;
            }
            else {
                Vec_IntSetEntry(tObj->FaninCs, j, 0);
            }
        }
        tObj->thre += sum;
    }
    // // sort weights and fanins from big to small
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
            Th_ObjSortWeight(tObj);
    //     if (tObj->Type == Th_CONST1 || tObj->Type == Th_Pi) continue;
    //     // bubble sort
    //     for (j = 0; j < Vec_IntSize(tObj->weights); ++j) {
    //         for (k = j + 1; k < Vec_IntSize(tObj->weights); ++k) {
    //             if (Vec_IntEntry(tObj->weights, j) < Vec_IntEntry(tObj->weights, k)) {
    //                 // swap weight
    //                 int wj = Vec_IntEntry(tObj->weights, j);
    //                 int wk = Vec_IntEntry(tObj->weights, k);
    //                 Vec_IntSetEntry(tObj->weights, k, wj);
    //                 Vec_IntSetEntry(tObj->weights, j, wk);
    //                 // swap fanin
    //                 int idj = Vec_IntEntry(tObj->Fanins, j);
    //                 int idk = Vec_IntEntry(tObj->Fanins, k);
    //                 Vec_IntSetEntry(tObj->Fanins, k, idj);
    //                 Vec_IntSetEntry(tObj->Fanins, j, idk);
    //                 // swap faninC
    //                 int cj = Vec_IntEntry(tObj->FaninCs, j);
    //                 int ck = Vec_IntEntry(tObj->FaninCs, k);
    //                 Vec_IntSetEntry(tObj->FaninCs, k, cj);
    //                 Vec_IntSetEntry(tObj->FaninCs, j, ck);
    //             }
    //         }
    //     }
    }
}

void Th_NtkTransBack(Vec_Ptr_t* vThres) {
    int i, j;
    Thre_S* tObj;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
        if (tObj->Type == Th_CONST1 || tObj->Type == Th_Pi) continue;
        int comp;
        int weight;
        Vec_IntForEachEntry(tObj->FaninCs, comp, j) {
            if (comp == 1) {
                weight = Vec_IntEntry(tObj->weights,j);
                Vec_IntSetEntry(tObj->weights, j, -weight);
                Vec_IntSetEntry(tObj->FaninCs, j, 0);
                tObj->thre -= weight;
            }
        }
    }
}

void Th_ObjSortWeight(Thre_S* tObj) {
    int j,k;
    // sort weights and fanins from big to small
    for (j = 0; j < Vec_IntSize(tObj->weights); ++j) {
        for (k = j + 1; k < Vec_IntSize(tObj->weights); ++k) {
            if (Vec_IntEntry(tObj->weights, j) < Vec_IntEntry(tObj->weights, k)) {
                // swap weight
                int wj = Vec_IntEntry(tObj->weights, j);
                int wk = Vec_IntEntry(tObj->weights, k);
                Vec_IntWriteEntry(tObj->weights, k, wj);
                Vec_IntWriteEntry(tObj->weights, j, wk);
                // swap fanin
                int idj = Vec_IntEntry(tObj->Fanins, j);
                int idk = Vec_IntEntry(tObj->Fanins, k);
                Vec_IntWriteEntry(tObj->Fanins, k, idj);
                Vec_IntWriteEntry(tObj->Fanins, j, idk);
                // swap faninC
                int cj = Vec_IntEntry(tObj->FaninCs, j);
                int ck = Vec_IntEntry(tObj->FaninCs, k);
                Vec_IntWriteEntry(tObj->FaninCs, k, cj);
                Vec_IntWriteEntry(tObj->FaninCs, j, ck);
            }
        }
    }
}