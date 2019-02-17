/**CFile****************************************************************
 
  FileName    [threTh2Dlist.c] 

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [threshold.]
  
  Synopsis    [Threshold network approximation.]

  Author      [Oscar Hung]
   
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
int    Th2DList(Thre_S * tObj);
void   Th_NtkTh2DList(Vec_Ptr_t * vThres);
void   Th_DLBuildSets(Vec_Ptr_t * vThres);
int    compare(Th_Set** i1, Th_Set** i2);

/*Function*************************************************************

  Synopsis    [Initialize 1-DL]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void   Th_NtkTh2DList(Vec_Ptr_t * vThres)
{
    int i;
    Thre_S * tObj;
    Vec_PtrForEachEntry(Thre_S*, vThres, tObj, i) {
        if (tObj->Type == Th_CONST1 || tObj->Type == Th_Pi) continue;
        /*if (!Th2DList(tObj)){                                                          */
        /*    Abc_Print(ABC_STANDARD, "Threshold gate %s is not 1-DL...\n", tObj->pName);*/
        /*    retValue = 0;                                                              */
        /*}*/
        Th2DList(tObj);
        /*Th_PrintNode(tObj);*/
    }
    Th_DLBuildSets(vThres);
}
/**Function*************************************************************

  Synopsis    [Transform a threshold gate to G-DL]

  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
int Th2DList(Thre_S * tObj)
{
    /*if(tObj->thre < 0) return 0;*/
    int fanin; 
    int w;
    int thre = tObj->thre;
    int sumWeight = Vec_IntSum(tObj->weights);
    int general = 0;
    int upper = -1;
    int lower = -1;
    Vec_Int_t * type2Node = Vec_IntStart(0);
    /*int debug = 0;*/
    tObj->dtypes = Vec_IntStart(0);
    tObj->thres = Vec_IntStart(0);
    Vec_IntForEachEntry(tObj->weights, w, fanin) {
        // type-2 node, ILP
        if(general) 
        {
            Vec_IntPush(tObj->dtypes, 2);
            Vec_IntPush(tObj->thres, -1);
            sumWeight -= w;
            
            Vec_IntPush(type2Node, w);
            upper = constructLP(type2Node, thre, -1, 0);
            /*printf("upper = %d, done!\n\n", upper);*/
            lower = constructLP(type2Node, thre-sumWeight, 0, -1);
            /*printf("lower = %d, done!\n\n", lower);*/
            if(upper == lower)
            {
              general = 0;
              upper = -1;
              lower = -1;
            }
            else
            {
              int start, weight;
              Vec_Int_t * following = Vec_IntStart(0);
              Vec_IntForEachEntryStart(tObj->weights, weight, start, fanin+1)
              {
                Vec_IntPush(following, weight);
              }
              if(constructLP(following, 0, upper, lower) == 0)
              {
                general = 0;
                upper = -1;
                lower = -1;
              }
              /*printf("equivalency done!\n");*/
              Vec_IntFree(following);
            }
            continue;
        }
        
        //  dummy inputs
        if (sumWeight < thre) {                                                  
         Thre_S* fi = Th_GetObjById(current_TList, Vec_IntEntry(tObj->Fanins,fanin));
         Vec_IntPop(tObj->weights);                                            
         Vec_IntPop(tObj->Fanins);                                             
         Vec_IntPop(tObj->FaninCs);                                            
         Vec_IntRemove(fi->Fanouts, tObj->Id);                                 
         --fanin;                                                                  
         continue;                                                             
        }

        //  last input
        if (fanin == Vec_IntSize(tObj->weights)-1)
        {
            // dummy
            /*if (w < thre) {                                                                 */
            /*    Thre_S* fi = Th_GetObjById(current_TList, Vec_IntEntry(tObj->Fanins,fanin));*/
            /*    Vec_IntPop(tObj->weights);                                                  */
            /*    Vec_IntPop(tObj->Fanins);                                                   */
            /*    Vec_IntPop(tObj->FaninCs);                                                  */
            /*    Vec_IntRemove(fi->Fanouts, tObj->Id);                                       */
            /*    break;                                                                      */
            /*}                                                                               */

            // assign dtype as the previous input
            if (fanin > 0) {
                Vec_IntPush(tObj->dtypes,Vec_IntEntryLast(tObj->dtypes));
                Vec_IntPush(tObj->thres, w);      
            }
        }

        // case1: 1->1
        else if (w >= thre) { 
            Vec_IntPush(tObj->thres, thre);      
            sumWeight -= w;
            Vec_IntPush(tObj->dtypes,1);
        }

        // case0: 0->0
        else if (sumWeight - w < thre) { 
            Vec_IntPush(tObj->thres, thre);      
            thre -= w;
            sumWeight -= w;
            Vec_IntPush(tObj->dtypes,0);
        }

        // case2: type-2
        else 
        {
            Vec_IntPush(tObj->thres, thre);      
            Vec_IntPush(tObj->dtypes,2);
            sumWeight -= w;
            general = 1;

            Vec_IntZero(type2Node);
            Vec_IntPush(type2Node, w);
            // check if there exists a majority gate
            /*if(w == Vec_IntEntry(tObj->weights,fanin+1))                                                           */
            /*{                                                                                                      */
            /*    if (2*w >= thre && sumWeight-2*w < thre) {                                                         */
            /*        Vec_IntPush(tObj->thres, -1);                                                                  */
            /*        Vec_IntPush(tObj->dtypes,2);                                                                   */
            /*        thre -= w;                                                                                     */
            /*        sumWeight -= 2*w;                                                                              */
            /*        fanin++;                                                                                       */
            /*        if (w == 1) general = 1; // remaining weights are all equal to 1, so they are all type-2 inputs*/
            /*    }                                                                                                  */
            /*    else general = 1;                                                                                  */
            /*}                                                                                                      */
            /*else general = 1;                                                                                      */
        }
    }
    return 1;
}
/**Function*************************************************************

  Synopsis    [Find sets of 1-DL for extraction]
1000000
  Description []
               
  SideEffects []

  SeeAlso     []

***********************************************************************/
void   Th_DLBuildSets(Vec_Ptr_t * vThres)
{
    int i;
    Thre_S * tObj;
    Golden_Vec = Vec_PtrStart(0);
    // int bound = Th_CountLevelWithDummy(vThres);
    // for (i = 0; i < bound; ++i) {
    // for (i = 0; i < 2; ++i) {
    //     Vec_Ptr_t* temp = Vec_PtrStart(0);
    //     Vec_PtrPush(Golden_Vec, temp);
    // }
    // int bound = Th_CountLevelWithDummy(vThres);
    /*printf("bound = %d\n", bound);*/
    int cnt = 0;
    Vec_PtrForEachEntry(Thre_S *, vThres, tObj, i)    
    {
        if (tObj->thre < 0) {++cnt;continue;}
        if (tObj->Type == Th_CONST1 || tObj->Type == Th_Pi) continue;
        if (!Vec_IntSize(tObj->dtypes)) continue;
        int type,k;
        int temp_type = 0;
        Vec_Int_t * temp_set = Vec_IntStart(0);
        Vec_IntForEachEntry(tObj->dtypes , type , k)
        {
            /*if (tObj->Id == 1331) printf("temp_type = %d, current type = %d\n", temp_type, type);*/
            if(type >= 2 && temp_type == 2)
            {
                Vec_IntClear(temp_set);
                continue;
            }
            if(k == Vec_IntSize(tObj->dtypes)-1 && k > 1)
            {
                // if(temp_type == 2 || type == 2) break;
                if(Vec_IntEntry(tObj->FaninCs,k) == 1)
                    Vec_IntPush(temp_set,((-1)*Vec_IntEntry(tObj->Fanins,k)));
                else
                    Vec_IntPush(temp_set,Vec_IntEntry(tObj->Fanins,k));
                // if(Vec_IntSize(temp_set) > 2) {
                //     printf("BUS\n");
                //     if(temp_type == 0) break;
                // }

                // else printf("hi %d\n",i);
                // assert(Vec_IntSize(temp_set) >= 2 && Vec_IntSize(temp_set) <= Vec_IntSize(tObj->dtypes));

                Th_Set * t = ABC_ALLOC(Th_Set, 1);
                t->set = Vec_IntDup(temp_set);
                Vec_IntSort(t->set,0);
                t->Dtype = temp_type;
                t->gate = tObj->Id;
                // push to space  
                // printf("level = %d\n", tObj->level);
                // if (tObj->level <= bound/2)
                //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 0),t);
                // else 
                //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 1),t);
                // else if (tObj->level <= bound*3/4)  
                //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 2),t); 
                // else
                //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 3),t);
                Vec_PtrPush(Golden_Vec,t);
                // Vec_PtrPush(Vec_PtrEntry(Golden_Vec, tObj->level),t);
                break;  
            }
            if(type != temp_type)
            {
                /*if (tObj->Id == 1331) printf("LuLU\n");*/
                if(Vec_IntSize(temp_set) > 1 )
                {
                    Th_Set * t = ABC_ALLOC(Th_Set, 1);
                    t->set = Vec_IntDup(temp_set);
                    //t->set = Vec_IntDup(temp_set);
                    Vec_IntSort(t->set,0);
                    t->Dtype = temp_type;
                    t->gate = tObj->Id; 
                    // push to space
                    // printf("level-1 = %d\n", tObj->level);
                    // if (tObj->level <= bound/2)
                    //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 0),t);
                    // else 
                    //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 1),t);
                    // else if (tObj->level <= bound*3/4)  
                    //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 2),t); 
                    // else
                    //     Vec_PtrPush(Vec_PtrEntry(Golden_Vec, 3),t);
                    Vec_PtrPush(Golden_Vec,t);
                    // Vec_PtrPush(Vec_PtrEntry(Golden_Vec, tObj->level),t);
                }
                Vec_IntClear(temp_set);
                temp_type = type;
            }
            if(Vec_IntEntry(tObj->FaninCs,k) == 1)
                Vec_IntPush(temp_set,((-1)*Vec_IntEntry(tObj->Fanins,k)));
            else
                Vec_IntPush(temp_set,Vec_IntEntry(tObj->Fanins,k));
        }
        Vec_IntFree(temp_set);
    }
    printf("Redundant = %d\n", cnt);
    // printf("Build Done!\n");
}

int compare(Th_Set** i1, Th_Set** i2) {
  Th_Set* p1 = *i1;
  Th_Set* p2 = *i2;
  int value1 = Vec_IntSize(p1->set);
  int value2 = Vec_IntSize(p2->set);
  if (value1 > value2) return 1;
  else if (value1 == value2) return 0;
  else return -1;
}

