#include <iostream>
#include <vector>
#include "../../lib/lp_solve_5.5.2.5_dev/lp_lib.h"
#include "base/abc/abc.h"

using namespace std;
#define M 1000
#define wRestrict 255

extern "C" int constructLP(Vec_Int_t * var, int thre, int upper, int lower);
extern "C" int solveLp(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, int max, int min);

int constructLP(Vec_Int_t * var, int thre, int upper, int lower) 
{
  vector< vector<int> > A;
  vector<int> B;
  vector<int> ans;
  for(int i = 0; i < Vec_IntSize(var); ++i)
  {
    //  x >= 0
    vector<int> temp;
    for(int j = 0; j < Vec_IntSize(var); ++j)
    {
      if(i==j) temp.push_back(1);
      else temp.push_back(0); 
    }
    A.push_back(temp);
    B.push_back(0);

    //  x <= 1 (-x >= -1)
    temp.clear();
    for(int j = 0; j < Vec_IntSize(var); ++j)
    {
      if(i==j) temp.push_back(-1);
      else temp.push_back(0); 
    }
    A.push_back(temp);
    B.push_back(-1);
  }

  int k, weight;
  vector<int> objective;
 
  // calculate upper bound
  if(upper == -1) 
  {
    //cout << "upper\n";
    Vec_IntForEachEntry(var, weight, k) 
    {
      objective.push_back(-1*weight);
    }
    A.push_back(objective);
    B.push_back(-1*thre+1);
    solveLp(A, B, ans, 1, 0);
    int ret = 0;
    Vec_IntForEachEntry(var, weight, k)
    {
      ret += weight * ans[k];
    }
    return ret;
  }
  
  // calculate lower bound
  else if(lower == -1)
  {
    //cout << "lower\n";
    Vec_IntForEachEntry(var, weight, k) 
    {
      objective.push_back(weight);
    }
    A.push_back(objective);
    B.push_back(thre);
    solveLp(A, B, ans, 0, 1); 
    int ret = 0;
    Vec_IntForEachEntry(var, weight, k)
    {
      ret += weight * ans[k];
    }
    return ret;
  }

  // calculate equivalency 
  else 
  {
    //cout << "equivalency\n";
    // sum(w*x) < upper
    Vec_IntForEachEntry(var, weight, k) 
    {
      objective.push_back(-1*weight);
    }
    A.push_back(objective);
    B.push_back(-1*upper+1);

    // sum(w*x) >= lower
    objective.clear();
    Vec_IntForEachEntry(var, weight, k) 
    {
      objective.push_back(weight);
    }
    A.push_back(objective);
    B.push_back(lower);
    return solveLp(A, B, ans, 0, 0);
  }

  //  debug
  //for(int i = 0; i < A.size(); ++i)
  //{
  //  for(int j = 0;j < A[i].size(); ++j)
  //  {
  //    cout << A[i][j] << " ";
  //  }
  //  cout << "| " << B[i] << endl;
  //}
  //solveLp(A, B, ans);
  //int l;
  //cin >> l;
  return 0;
}
int solveLp(vector< vector<int> >& A, vector<int>& B, vector<int>& ans, int max, int min)
{
   lprec *lp;
   int Ncol, *colno = NULL, j, ret = 0;
   REAL *row = NULL;
   
   Ncol = A[0].size();
   lp = make_lp(0, Ncol);
   if(lp == NULL)
      ret = 1; /* couldn't construct a new model... */
   
   if(ret == 0) {
      for (int i=0; i<Ncol; ++i)
         set_int(lp,i+1,true);
         //set_binary(lp,i+1,true);
      
      /* create space large enough for one row */
      colno = (int *) malloc(Ncol * sizeof(*colno));
      row = (REAL *) malloc(Ncol * sizeof(*row));
      if((colno == NULL) || (row == NULL))
         ret = 2;
      set_add_rowmode(lp, TRUE);  /* makes building the model faster if it is done rows by row */
   }
   
   int k=0;
   while ((ret == 0)&&(k<A.size())) {
      j = 0;
      for (int i=0; i<A[k].size(); ++i)
      {
         if (A[k][i]==0) continue;
         colno[j]=i+1;
         row[j++]=A[k][i];
      }
      
      /* add the row to lpsolve */
      if(!add_constraintex(lp, j, row, colno, GE, B[k]))
         ret = 3;
      k++;
   }

   if(ret == 0) {
      set_add_rowmode(lp, FALSE); /* rowmode should be turned off again when done building the model */
      
      /* set the objective function */
      j = 0;
      
      for (int i=0; i<A[0].size(); ++i)
      {
         colno[j]=i+1;
         if(max)      row[j++] = -1*A[A.size()-1][i];
         else if(min) row[j++] = A[A.size()-1][i];
         else         row[j++] = 0;
      }
      
      /* set the objective in lpsolve */
      if(!set_obj_fnex(lp, j, row, colno))
         ret = 4;
   }
   
   if(ret == 0) {
      /* set the object direction to minimize */
      if(max) set_maxim(lp);
      else set_minim(lp);
      
      /* just out of curioucity, now show the model in lp format on screen */
      //write_LP(lp, stdout);
      
      /* I only want to see important messages on screen while solving */
      set_verbose(lp,NEUTRAL);// IMPORTANT);
      
      /* Now let lpsolve calculate a solution */
      ret = solve(lp);

      if(ret == OPTIMAL)
         ret = 0;
      else
         ret = 5;
   }
   if (ret==0) {
   
      /* objective value */
      //printf("Objective value: %f\n", get_objective(lp));
      
      /* variable values */
      get_variables(lp, row);
      for(j = 0; j < Ncol; j++)
      {
         //printf("%s: %f\n", get_col_name(lp, j + 1), row[j]);
         ans.push_back(row[j]);
         //cout<<row[j]<<" ";
         assert(ans.back()>=0);
      }
      //cout<<endl;
      /* we are done now */
   }
   
   /* free allocated memory */
   if(row != NULL)
      free(row);
   if(colno != NULL)
      free(colno);
   
   if(lp != NULL) {
      /* clean up such that all used memory by lpsolve is freed */
      delete_lp(lp);
   }
   if (ret != 0) return 0;
   //if (!min && !max) {
     //cout << "ans.size() = " << ans.size() << endl;
   //}
   //cout << "LU\n";
   return(ans.size());
}
