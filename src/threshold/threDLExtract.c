/**CFile****************************************************************
 
  FileName    [threDLExtract.c] 

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [threshold.]
  
  Synopsis    [Threshold network approximation.]

  Author      [Oscar Hung, Hao Chen]
   
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
void     Th_IterativeUpdate(Vec_Ptr_t* vThres, Vec_Ptr_t* vSet);
int      FindMaxSet(Vec_Ptr_t * vThres, Vec_Ptr_t * vSet, int bound, int n);
void     Th_ObjUpdate(Vec_Ptr_t* vThres, Thre_S * tObj);
void     Th_RemoveNode(Vec_Ptr_t* vThres, Thre_S* tObj);
void     Th_Decompose(Vec_Ptr_t * vThres, Thre_S *);
Thre_S * Th_CreateTempObj(Vec_Ptr_t *TList, Th_Gate_Type Type);
int      Th_FindCritical(int l);
void     Th_FindCritical_rec(Thre_S* tObj, int level);
int      Th_CountLevelWithDummy( Vec_Ptr_t * tList );
void     Th_ComputeLevelWithDummy( Vec_Ptr_t * tList );
void     Th_ComputeLevelWithDummy_rec( Vec_Ptr_t * tList , Thre_S * tObj );
void     Th_RemoveDummy(); 
void     Th_IterativeUpdatePQ(Vec_Ptr_t* vThres, Vec_Ptr_t* vSet, int levelbound);
/**Function*************************************************************

  Synopsis    [Find the max set]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void Th_IterativeUpdate(Vec_Ptr_t* vThres, Vec_Ptr_t* vSet) {
    Th_Set* tObj;
    int level_ori = Th_CountLevel(current_TList) + 3;
    int l = 0;

    // find max fanio number
    // Th_Set* tObj;
    int i;
    int n = 0;
    // int numFan;
	Vec_PtrForEachEntry( Th_Set * , vSet , tObj , i )
	{
        if (Vec_IntSize(tObj->set) > n) n = Vec_IntSize(tObj->set);
	}

    // Thre_S* tObj;
    while (1) {
        if (n < 2) break;
        FindMaxSet(vThres, vSet, 0, n);
        if (n == 2) break;
        ++l;
        --n;
        // break;
    }

}

void Th_IterativeUpdatePQ(Vec_Ptr_t* vThres, Vec_Ptr_t* vSet, int levelbound) {
    Th_Set* tObj;
    Thre_S* original;
    int i;
    int j, id, k;
    int size = 100;
    Vec_Ptr_t* newSet = Vec_PtrStart(0);
    Vec_PtrSort(vSet, (int (*)(void))compare);
    /*printf("vSet size = %d\n", Vec_PtrSize(vSet));*/
    // Vec_PtrForEachEntry(Th_Set*, vSet, tObj, i) {
    //     printf("#set: %d, gate:%d\n", Vec_IntSize(tObj->set), tObj->gate);
    // }
    // return;
    /*int bound = Th_FindCritical(levelbound);*/
    int change = 1;
    int bound = 0;
   /*while(Vec_PtrSize(vSet) != 0) {     */
     /*tObj = Vec_PtrPop(vSet);          */
     /*printf("gate = %d\n", tObj->gate);*/
     /*Vec_IntForEachEntry(tObj->set, id, k) {*/
       /*printf("%d\t",id);*/
     /*}*/
     /*printf("\n");*/
   /*}                                   */
   /*return;                             */

    while(Vec_PtrSize(vSet) > 1) {
        if (change) {
          bound = Th_FindCritical(levelbound);
          change = 0;
          /*printf("Change!\n");*/
        }
        /*printf("bound: %d\n", bound);*/
        tObj = Vec_PtrPop(vSet);
        if (bound == 1) {
          if (abs(Vec_IntEntry(Critical, tObj->gate)) == 1) continue;
        }

        int unavailable = 0;
        original = Th_GetObjById(vThres, tObj->gate);
        Vec_IntForEachEntry(tObj->set, id, k) {
            if (Vec_IntFind(original->Fanins, abs(id)) == -1) {
                unavailable = 1;
                break;
            }
          //  if (bound == 1) {
          //    if (id >= 0) {
          //      if (Vec_IntEntry(Critical, id) == 1) {
          //        unavailable = 1;
          //        break;
          //      }
          //    }
          //    else {
          //      if (Vec_IntEntry(Critical, -id) == -1) {
          //        unavailable = 1;
          //        break;
          //      }
          //    }   
          //  }
        }
        if (unavailable) {
            /*printf("\tUnavailable!!!\n");*/
            continue;
        }

        // tObj set
        Vec_Int_t * iValue = Vec_IntStart(nPi+nPo+nThres+1);
        Vec_Int_t * Max_gate = Vec_IntStart(0);
        Vec_IntForEachEntry(tObj->set,id,j)
        {
            if (id < 0) 
            {
                Vec_IntWriteEntry(iValue,-id,-1);
                // if(bound==1)
                //     if(Vec_IntEntry(Critical,-id) == 1) 
                //         cp = 1;
            }
            else 
            {
                Vec_IntWriteEntry(iValue,id,1);
                // if(bound==1)
                //     if(Vec_IntEntry(Critical,id) == 1) 
                //         cp = 1;
            }
        }
        // set of every Th_Set
        Th_Set* jObj;
        Vec_Int_t * temp_fanin = Vec_IntStart(0);
        Vec_Ptr_t* overlapSet = Vec_PtrStart(0);  // as new obj fanout
        Vec_Ptr_t* criticalSet = Vec_PtrStart(0);  // as new obj fanout
        Vec_PtrPush(overlapSet, tObj); // Add tObj
        Vec_PtrForEachEntry(Th_Set * ,vSet, jObj, j)
        {
            if (jObj->gate == tObj->gate) continue;                
            if (jObj->Dtype != tObj->Dtype) continue;
            if (bound == 1) {
              if (abs(Vec_IntEntry(Critical, jObj->gate)) == 1) {
                /*Vec_PtrPush(criticalSet, jObj);*/
                continue;
              }
            }
            /*int critic = 0;*/
            Vec_IntClear(temp_fanin);
            Vec_IntForEachEntry(jObj->set,id,k)
            {
                /*if (bound == 1) {                           */
                /*  if (id >= 0) {                            */
                /*    if (Vec_IntEntry(Critical, id) == 1) {  */
                /*      critic = 1;                           */
                /*      break;                                */
                /*    }                                       */
                /*  }                                         */
                /*  else {                                    */
                /*    if (Vec_IntEntry(Critical, -id) == -1) {*/
                /*      critic = 1;                           */
                /*      break;                                */
                /*    }                                       */
                /*  }                                         */
                /*}                                           */
                if (id < 0 && Vec_IntEntry(iValue,-id) == -1)
                {
                    Vec_IntPush(temp_fanin,id);
                }
                else if (id >=0 && Vec_IntEntry(iValue,id) == 1)
                {
                    Vec_IntPush(temp_fanin,id);
                }
            }
            /*if (critic) continue;*/
            if (Vec_IntSize(temp_fanin) > Vec_IntSize(Max_gate))
            {
                Vec_IntClear(Max_gate);
                Vec_IntForEachEntry(temp_fanin,id,k) Vec_IntPush(Max_gate,id);
                Vec_PtrClear(overlapSet);
                Vec_PtrPush(overlapSet, jObj);
                Vec_PtrPush(overlapSet, tObj);
            }
            else if (Vec_IntSize(temp_fanin) == Vec_IntSize(Max_gate)) {
                int lu = 0;
                Vec_IntForEachEntry(temp_fanin, id, k) {
                    if (Vec_IntFind(Max_gate, id) == -1) {
                        lu = 1;
                        break;
                    }
                }
                if (lu == 0) Vec_PtrPush(overlapSet, jObj);
            }
        }
        /*Vec_PtrForEachEntry(Th_Set*, criticalSet, jObj, i) {*/
        /*  Vec_PtrRemove(vSet, jObj);                        */
        /*}                                                   */
        /*Vec_PtrFree(criticalSet);*/
        Vec_IntFree(temp_fanin);
        Vec_IntFree(iValue);
        if (Vec_IntSize(Max_gate) < 2) {
            /*printf("\t#Fanin = %d < 2!!!\n", Vec_IntSize(Max_gate));*/
            Vec_PtrFree(overlapSet);
            Vec_IntFree(Max_gate);
            continue;
        }

        // Find overlap set with max_gate
        // assert(Vec_IntSize(Max_gate) != Vec_IntSize(tObj->set));
        // Vec_Ptr_t* overlapSet = Vec_PtrStart(0);  // as new obj fanout
        // Vec_PtrPush(overlapSet, tObj); // Add tObj

        // Vec_PtrForEachEntry(Th_Set *, vSet, jObj, j) {
        //     int RetValue = 0;
        //     if (jObj->Dtype != tObj->Dtype) continue;
        //     Vec_Int_t* maxValue = Vec_IntStart(nPi+nPo+nThres+1);  // check overlap
        //     Vec_IntForEachEntry(jObj->set, id, k) {
        //         if (id >= 0) {
        //             Vec_IntWriteEntry(maxValue, id, 1);
        //         }
        //         else {
        //             Vec_IntWriteEntry(maxValue, -id, -1);
        //         }
        //     }
        //     Vec_IntForEachEntry(Max_gate, id, k) {
        //         if (id >= 0) {
        //             if (Vec_IntEntry(maxValue, id) != 1){
        //                 RetValue = 1;
        //                 break;
        //             }
        //         }
        //         else {
        //             if (Vec_IntEntry(maxValue, -id) != -1){
        //                 RetValue = 1;
        //                 break;
        //             }
        //         }
        //     }
        //     if (RetValue == 0) Vec_PtrPush(overlapSet, jObj);
        //     Vec_IntFree(maxValue);
        // }
        // Create new Obj
        if (Vec_IntSize(overlapSet) > 2 || Vec_IntSize(Max_gate) > 2){
            Thre_S* rObj = Th_CreateTempObj(vThres, Th_Node);
            rObj->dtypes = Vec_IntStart(0);
            Vec_IntForEachEntry(Max_gate,id,k)
            {
                if(id >= 0) 
                {
                    Vec_IntPush(rObj->Fanins,id);
                    Vec_IntPush(rObj->FaninCs,0);
                }
                else
                {
                    Vec_IntPush(rObj->Fanins,-id);
                    Vec_IntPush(rObj->FaninCs,1);
                }
                Vec_IntPush(rObj->weights, 1);
                Vec_IntPush(rObj->dtypes, tObj->Dtype);
            }
            if (tObj->Dtype == 1) rObj->thre = 1;
            else rObj->thre = Vec_IntSize(Max_gate);
            // Put in Golden_Vec if #fanin > 2
            if (Vec_IntSize(Max_gate) > 2)
            {
                Th_Set * t = ABC_ALLOC(Th_Set , 1);
                t->set = Vec_IntDup(Max_gate);
                t->Dtype = tObj->Dtype;
                t->gate = rObj->Id;
                Vec_PtrPush(newSet, t);
                // pqueue_enqueue(pq, &t);
            }
            // Renew set & rObj->fanout
            Vec_PtrForEachEntry(Th_Set*, overlapSet, jObj, j) {
                Vec_IntPush(rObj->Fanouts, jObj->gate);
                if(Vec_IntSize(jObj->set) - Vec_IntSize(Max_gate) < 2) 
                {
                    if (jObj->gate != tObj->gate) Vec_PtrRemove(vSet,jObj);
                }
                else
                {
                    Vec_IntForEachEntry(Max_gate,id,k) Vec_IntRemove(jObj->set,id);
                    Vec_PtrPush(newSet, jObj);
                    if (jObj->gate != tObj->gate) Vec_PtrRemove(vSet,jObj);
                }
            }
            Vec_PtrPush(vThres, rObj);
            Th_ObjUpdate(vThres, rObj);
            /*return;*/
            change = 1;
            /*printf("\tFanin_num : %d, Fanout_num : %d \n",Vec_IntSize(rObj->Fanins),Vec_IntSize(rObj->Fanouts));*/
        }
        /*else {                         */
        /*    printf("\tNot #Fanout!\n");*/
        /*}                              */
    }
    Golden_Vec = Vec_PtrDup(newSet);
    Vec_PtrFree(newSet);
    /*printf("NewSet size = %d\n", Vec_PtrSize(Golden_Vec));*/
}

int FindMaxSet(Vec_Ptr_t * vThres, Vec_Ptr_t * vSet, int levelbound, int n)
{  
    // printf("n = %d\n", n);
    int i,j;
    // int temp_type;
    Th_Set * iObj;
    Th_Set * jObj;
    Thre_S * rObj;

    // Critical
    int bound = Th_FindCritical(levelbound);

    // Vec_Int_t * Max_gate = Vec_IntStart(0);//gate
    // printf("#set = %d\n", Vec_PtrSize(vSet));
    Vec_Ptr_t* store = Vec_PtrStart(0);
    Vec_Int_t* type = Vec_IntStart(0);

    Vec_PtrForEachEntry(Th_Set *, vSet, iObj, i)
    {
        int id,k,cp;
        cp = 0;
        Vec_Int_t * iValue = Vec_IntStart(nPi+nPo+nThres+1);
        Vec_IntForEachEntry(iObj->set,id,k)
        {
            if (id < 0) 
            {
                Vec_IntWriteEntry(iValue,-id,-1);
                if(bound==1)
                    if(Vec_IntEntry(Critical,-id) == 1) 
                        cp = 1;
            }
            else 
            {
                Vec_IntWriteEntry(iValue,id,1);
                if(bound==1)
                    if(Vec_IntEntry(Critical,id) == 1) 
                        cp = 1;
            }
        }
        if(cp == 1) continue;
        Vec_Int_t * temp_fanin = Vec_IntStart(0);
        Vec_PtrForEachEntryStart(Th_Set * ,vSet, jObj, j , i+1 )
        {
            if(jObj->Dtype != iObj->Dtype) continue;
            Vec_IntClear(temp_fanin);
            Vec_IntForEachEntry(jObj->set,id,k)
            {
                if (id < 0 && Vec_IntEntry(iValue,-id) == -1)
                {
                    Vec_IntPush(temp_fanin,id);
                }
                else if (id >=0 && Vec_IntEntry(iValue,id) == 1)
                {
                    Vec_IntPush(temp_fanin,id);
                }
            }
            // if (Vec_IntSize(temp_fanin) > Vec_IntSize(Max_gate))
            // {
            //     Vec_IntClear(Max_gap_fanin);
            //     Vec_IntForEachEntryp_fanin);ntPush(Max_gate,id);
            //     temp_type = jObj->Dp_fanin);
            //     Vec_PtrClear(store)p_fanin);
            //     Vec_PtrPush(store, p_fanin);;
            //     Vec_IntClear(type);p_fanin);
            //     Vec_IntPush(type, jObj->Dtype);
            //     MaxFanin = Vec_IntSize(temp_fanin);
            // }
            // else if (Vec_IntSize(temp_fanin) == Vec_IntSize(Max_gate) && Vec_IntSize(temp_fanin) != 0) {
            //     Vec_PtrPush(store, Vec_IntDup(temp_fanin));
            //     Vec_IntPush(type, jObj->Dtypp_fanin);
            // }
            if (Vec_IntSize(temp_fanin) == n) {        
                Vec_PtrPush(store, Vec_IntDup(temp_fanin));
                Vec_IntPush(type, jObj->Dtype); 
            }
        }
        Vec_IntFree(iValue);
        Vec_IntFree(temp_fanin);
    }
    if (Vec_PtrSize(store) < 2) {
        // Vec_IntFree(Max_gate);
        return 0;
    }
    Vec_Int_t* tt;
    int available = 0;
    // return 0;
    int u,p;
    Vec_PtrForEachEntry(Vec_Int_t*, store, tt, u) {
        // Vec_IntPrint(tt);
        int id,k;
        Vec_Ptr_t* mod = Vec_PtrStart(0);
        bound = Th_FindCritical(levelbound);
        Vec_PtrForEachEntry(Th_Set *, vSet, iObj, p)
        {
            if(iObj->Dtype != Vec_IntEntry(type, u)) continue;
            Vec_Int_t * iValue = Vec_IntStart(nPi+nPo+nThres+1);
            int cp = 0;
            Vec_IntForEachEntry(iObj->set,id,k)
            {
                if (id < 0) 
                {
                    Vec_IntWriteEntry(iValue,-id,-1);
                    if(bound==1)
                        if(Vec_IntEntry(Critical,-id) == 1) 
                            cp = 1;
                }
                else 
                {
                    Vec_IntWriteEntry(iValue,id,1);
                    if(bound==1)
                        if(Vec_IntEntry(Critical,id) == 1) 
                            cp = 1;
                }
            }
            if(cp == 1) {
                Vec_IntFree(iValue);
                continue;
            }
            int RetValue = 0;
            Vec_IntForEachEntry(tt,id,k)
            {
                if (id >= 0 && Vec_IntEntry(iValue,id) != 1)
                {
                    RetValue = 0;
                    break;
                }
                else if (id < 0 && Vec_IntEntry(iValue,-id) != -1)
                {
                    RetValue = 0;
                    break;
                }
                else ++RetValue;
            }
            if(RetValue) 
            {
                Vec_PtrPush(mod, iObj);
            }
            Vec_IntFree(iValue);
        }
        // printf("Mod size = %d\n",Vec_PtrSize(mod)  );
        if (Vec_PtrSize(mod) > 1) {
            // if ((Vec_PtrSize(tt)+3) < Vec_PtrSize(mod)*(Vec_PtrSize(tt)-1)) { 
                // printf("\t%d + 3 < %d * (%d - 1)\n", Vec_PtrSize(tt), Vec_PtrSize(mod), Vec_PtrSize(tt));
                rObj = Th_CreateTempObj(vThres, Th_Node);
                // Add new node
                if (Vec_IntSize(tt) > 2 || Vec_PtrSize(mod) > 2) {
                    available = 1;
                    rObj->dtypes = Vec_IntStart(0);
                    Vec_IntForEachEntry(tt,id,k)
                    {
                        if(id >= 0) 
                        {
                            Vec_IntPush(rObj->Fanins,id);
                            Vec_IntPush(rObj->FaninCs,0);
                        }
                        else
                        {
                            Vec_IntPush(rObj->Fanins,-id);
                            Vec_IntPush(rObj->FaninCs,1);
                        }
                        Vec_IntPush(rObj->weights,1);
                        Vec_IntPush(rObj->dtypes,Vec_IntEntry(type, u));
                    }
                    if (Vec_IntEntry(type, u) == 1) rObj->thre = 1;
                    else rObj->thre = Vec_IntSize(tt);

                    if (Vec_IntSize(tt) > 2)
                    {
                        Th_Set * t = ABC_ALLOC(Th_Set , 1);
                        t->set = Vec_IntDup(tt);
                        t->Dtype = Vec_IntEntry(type, u);
                        t->gate = rObj->Id;
                        Vec_PtrPush(vSet, t);
                    }
                }
                // renew vSet
                Vec_PtrForEachEntry(Th_Set*, mod, iObj, j) {
                    if(Vec_IntSize(rObj->Fanins) > 0) Vec_IntPush(rObj->Fanouts,iObj->gate);
                    if(Vec_IntSize(iObj->set) - Vec_IntSize(tt) < 2) 
                    {
                        Vec_PtrRemove(vSet,iObj);
                        --i;
                    }
                    else
                    {
                        Vec_IntForEachEntry(tt,id,k) Vec_IntRemove(iObj->set,id);
                    }
                }
                if (Vec_IntSize(rObj->Fanins) > 0) {
                    assert(available);
                    Vec_PtrPush(vThres,rObj);
                    Th_ObjUpdate(vThres, rObj);
                    printf("\tFanin_num : %d, Fanout_num : %d \n",Vec_IntSize(rObj->Fanins),Vec_IntSize(rObj->Fanouts));
                    Vec_IntForEachEntry(rObj->Fanins, id, k) printf("%d ", id);
                    printf("\n");
                    Vec_IntForEachEntry(rObj->Fanouts, id, k) printf("%d ", id);
                    printf("\n");
                    // printf("%d < %d\n", (Vec_PtrSize(tt)+3), Vec_PtrSize(mod)*(Vec_PtrSize(tt)-1));
                }
            // }
        }

    }
    if(available) return 1;
    else return 0;
    
}

void Th_ObjUpdate(Vec_Ptr_t * vThres, Thre_S * tObj) {
    int faninId, fanoutId;
    int i, j;
    // update faninGs
    Vec_IntForEachEntry(tObj->Fanins, faninId, i) {
        Thre_S* faninG = Th_GetObjById(vThres, faninId);
        Vec_IntForEachEntry(tObj->Fanouts, fanoutId, j) {
            Vec_IntRemove(faninG->Fanouts, fanoutId);
        }
        Vec_IntPush(faninG->Fanouts, tObj->Id);
    }
    // update fanoutGs
    int check = 0;
    Vec_Int_t* de = Vec_IntStart(0);
    Vec_IntForEachEntry(tObj->Fanouts, fanoutId, i) {
        Thre_S* fanoutG = Vec_PtrEntry(vThres, fanoutId);
        int index;
        int weight_old;
        check = 0;
        int ccc = 0;
        Vec_IntClear(de);
        Vec_IntForEachEntry(tObj->Fanins, faninId, j) {
            index = Vec_IntFind(fanoutG->Fanins, faninId);
            if(j == 0)
            {
                Vec_IntPush(de,Vec_IntEntry(fanoutG->dtypes, index));
                weight_old = Vec_IntEntry(fanoutG->weights, index);
            }
            else
                if(weight_old != Vec_IntEntry(fanoutG->weights, index)) 
                {
                    check = 1;
                    if(Vec_IntEntry(fanoutG->weights, index) > weight_old)
                        weight_old = Vec_IntEntry(fanoutG->weights, index);
                }
            ccc += Vec_IntEntry(fanoutG->weights, index);
            Vec_IntPush(de,Vec_IntEntry(fanoutG->weights, index));
            // printf("%d ", weight_old);
            Vec_IntDrop(fanoutG->weights, index);
            Vec_IntDrop(fanoutG->Fanins, index);
            Vec_IntDrop(fanoutG->FaninCs, index);
            Vec_IntDrop(fanoutG->dtypes, index);
        }
         /*printf("%d\n",fanoutG->Id);*/
         Vec_IntInsert(fanoutG->Fanins, index, tObj->Id);
         Vec_IntInsert(fanoutG->FaninCs, index, 0); // must be 0
         Vec_IntInsert(fanoutG->dtypes, index, Vec_IntEntry(tObj->dtypes, 0));
         Vec_IntInsert(fanoutG->weights, index, weight_old); 
        /*Vec_IntPush(fanoutG->Fanins, tObj->Id);                     */
        /*Vec_IntPush(fanoutG->FaninCs, 0); // must be 0              */
        /*Vec_IntPush(fanoutG->dtypes, Vec_IntEntry(tObj->dtypes, 0));*/
        /*Vec_IntPush(fanoutG->weights, weight_old);                  */
        if(Vec_IntEntry(tObj->dtypes, 0) == 0)
            fanoutG->thre -= ccc - weight_old;
        Th_ObjSortWeight(fanoutG);
    }
    Vec_IntFree(de);
}

void Th_RemoveDummyNode(Vec_Ptr_t* vThres, Thre_S* tObj) {
    assert(tObj->Type == Th_Node);
    assert(Vec_IntSize(tObj->Fanins) == 1);
    /*printf("Remove Gate %d\n", tObj->Id);*/

    int i;
    int fanoutId;
    Thre_S* faninG = Vec_PtrEntry(vThres, Vec_IntEntry(tObj->Fanins, 0));
    Vec_IntRemove(faninG->Fanouts, tObj->Id);
    Vec_IntForEachEntry(tObj->Fanouts, fanoutId, i) {
        Vec_IntPush(faninG->Fanouts, fanoutId);
    }
    Vec_IntForEachEntry(tObj->Fanouts, fanoutId, i) {
        Thre_S* fanoutG = Vec_PtrEntry(vThres, fanoutId);
        int index = Vec_IntFind(fanoutG->Fanins, tObj->Id);
        Vec_IntSetEntry(fanoutG->Fanins, index, faninG->Id);
    }
    Th_DeleteObj(tObj);
    /*Th_NtkDfs(current_TList);*/
    /*printf("Done!\n");*/
}

Thre_S *
Th_CreateTempObj(Vec_Ptr_t *TList, Th_Gate_Type Type)
{
	Thre_S *tObj 	= ABC_ALLOC(Thre_S, 1);
	tObj->thre 		= 0;
	tObj->Type 		= Type;
	tObj->Ischoose 	= 0;
	tObj->Id 		= Vec_PtrSize(TList);
	tObj->oId 		= 0;
	tObj->nId 		= 0;
	tObj->cost 		= 0;
	tObj->level 	= 0;
	tObj->pName 	= NULL;
	tObj->weights 	= Vec_IntAlloc(16);
	tObj->Fanins 	= Vec_IntAlloc(16);
	tObj->Fanouts 	= Vec_IntAlloc(16);
	tObj->FaninCs   = Vec_IntAlloc(16);
	tObj->pCopy 	= NULL;

	return tObj;
}

/**Function*************************************************************

  Synopsis    [Find critical path]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int      
Th_FindCritical(int l)
{
    if(l == 0) return 0;
    int level;
    // level = Th_CountLevelWithDummy(current_TList)-1;
    level = Th_CountLevelWithDummy(current_TList);
    /*printf("level = %d, level bound = %d\n",level-1,l-1);*/
    if(level < l) return 0;
    
    //     int i;
    //     Thre_S* kObj;
    //     Vec_PtrForEachEntry(Thre_S*, current_TList, kObj, i)
    //     {
    //         printf("Gate %d level = %d\n", kObj->Id, kObj->level);
    //     }
    //     return 2;
    // }

    
    // if (Critical) Vec_IntFree(Critical);
    Critical = Vec_IntStart( Vec_PtrSize(current_TList) );
    Thre_S * tObj;
    int i;
    Vec_PtrForEachEntry(Thre_S *, current_TList, tObj, i)
    {
        if(tObj && tObj->Type == Th_Pi)
        {
            // printf("Fanin  %d, level = %d\n",tObj->Id, tObj->level);
            assert(tObj->level <= level);
            if(tObj->level == level)
            { 
                Vec_IntWriteEntry(Critical, tObj->Id, 1);
                Th_FindCritical_rec(tObj, level-1);
            }
        }
    }
    /*if (level > l) return 2;*/
    /*printf("LULU\n");*/
    return 1;
}

void
Th_FindCritical_rec(Thre_S* tObj, int level)
{
    /*printf("Current level : %d, GateId : %d\n", level+1, tObj->Id);*/
    // if(level == -1) assert(tObj->Type == Th_Po);
    // if(tObj) if(tObj->Type == Th_Po) printf("lulu level = %d\n",level);
    if(level <= 0) return;
    int id,i;
    Vec_IntForEachEntry(tObj->Fanouts, id, i)
    {
        Thre_S * fObj = Th_GetObjById(current_TList, id);
        if(Vec_IntSize(fObj->Fanins) == 1 && fObj->Type == Th_Node && fObj->level == level) 
        {
            // level++;
            // printf("Dummy! ID = %d\n", id);
            // printf("Parent level : %d, Dummy node level : %d\n", tObj->level, fObj->level);
            if (Vec_IntEntry(Critical, id) == 1) continue;
            Vec_IntWriteEntry(Critical, id, 1);
            Th_FindCritical_rec(fObj, level);
        }
        else if (fObj->level == level)
        {
            if (Vec_IntEntry(Critical, id) == 1) continue;
            Vec_IntWriteEntry(Critical, id, 1);
            Th_FindCritical_rec(fObj, level-1);
        }
    }
}

int Th_CountLevelWithDummy( Vec_Ptr_t * tList )
{
	Th_SetLevel  ( tList );
	Th_ComputeLevelWithDummy( tList );
	return Th_FindMaxLevel( tList );
}

void Th_ComputeLevelWithDummy( Vec_Ptr_t * tList )
{
	Thre_S * tObj;
	Vec_Ptr_t * vPo;
	int i;

	vPo = Vec_PtrAlloc(16);

	Vec_PtrForEachEntry( Thre_S * , tList , tObj , i )
	{
		if ( tObj && tObj->Type == Th_Po )
			Vec_PtrPush( vPo , tObj );
	}
   
	Vec_PtrForEachEntry( Thre_S * , vPo , tObj , i )
	{
		Th_ComputeLevelWithDummy_rec( tList , tObj );
	}

	Vec_PtrFree(vPo);
}

void Th_ComputeLevelWithDummy_rec( Vec_Ptr_t * tList , Thre_S * tObj )
{
   Thre_S * tObjFanin;
	int i , Entry;
  
	//printf( "Dump tObj: \n");
	//Th_DumpObj( tObj );
	if ( tObj->Type == Th_Po || tObj->Type == Th_Node ) {
	   Vec_IntForEachEntry( tObj->Fanins , Entry , i ) 
		{
			tObjFanin = Th_GetObjById( tList , Entry );
	      //printf( "Dump tObjFanin: \n");
	      //Th_DumpObj( tObjFanin );
			assert(tObjFanin);
            if(Vec_IntSize(tObjFanin->Fanins) == 1 && tObjFanin->Type == Th_Node) { 
                if(tObjFanin->level < tObj->level) {
                    tObjFanin->level = tObj->level;
                    // Abc_Print(ABC_STANDARD, "LuLu %d gate\n", tObjFanin->Id);
                    // printf("Dummy node level %d\n", tObjFanin->level);
                    Th_ComputeLevelWithDummy_rec( tList , tObjFanin );
                }
            }
			else if ( tObjFanin-> level < tObj->level + 1 ) {
                // if (Vec_IntSize(tObj->Fanins) == 1 && tObj->Type == Th_Node) 
                //     printf("Dummy level %d, level %d\n", tObj->level, tObjFanin->level);
			   tObjFanin->level = tObj->level + 1;
			   Th_ComputeLevelWithDummy_rec( tList , tObjFanin );
			}
		}
	}
}

void Th_RemoveDummy() {
    // remove dummy nodes
    int counter = 0;
    Thre_S* tObj;
    /*while(1) {*/
        /*printf("No.%d remove!\n", counter);*/
        /*int check = 0;*/
        int i;
        Vec_PtrForEachEntry(Thre_S*, current_TList, tObj, i) {
            if (tObj->Type != Th_Node) continue;
            if (Vec_IntSize(tObj->Fanins) == 1) {
                Th_RemoveDummyNode(current_TList, tObj);
                /*check = 1;*/
                /*break;*/
                /*printf("Remove Dummy: %d\n", ++counter);*/
            }
        }
        /*if (check == 0) break;*/
        /*++counter;*/
    /*}*/
        Th_NtkDfs(current_TList);
}


void   Th_Decompose(Vec_Ptr_t * vThres, Thre_S * tObj) {
  int num = 5;
  int id, i;
  int count = 0;
  Thre_S* rObj = Th_CreateTempObj(vThres, Th_Node);
  
  printf("ID = %d\n", tObj->Id);
  Th_PrintNode(tObj);
  rObj->dtypes = Vec_IntStart(0);
  rObj->thres = Vec_IntStart(0);

  Vec_IntForEachEntryReverse(tObj->Fanins, id, i) {
    Vec_IntPushFirst(rObj->Fanins, id);
    Vec_IntPushFirst(rObj->dtypes, Vec_IntEntry(tObj->dtypes, i));
    Vec_IntPushFirst(rObj->FaninCs, Vec_IntEntry(tObj->FaninCs, i));
    Vec_IntPushFirst(rObj->weights, Vec_IntEntry(tObj->weights, i));
    Vec_IntPushFirst(rObj->thres, Vec_IntEntry(tObj->thres, i));
    ++count;
    if(count == num) 
    {
      rObj->thre = Vec_IntEntry(tObj->thres, i);
      break;
    }
    Vec_IntPop(tObj->thres);
  }
  Vec_IntPush(rObj->Fanouts, tObj->Id);

  Vec_PtrPush(vThres, rObj);
  Th_ObjUpdate(vThres, rObj);
  Th_PrintNode(rObj);
  Th_PrintNode(tObj);
}
