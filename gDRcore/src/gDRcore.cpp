#include <algorithm>
#include <string>
#include <vector>
#include <numeric>
#include <iostream>
#define R_NO_REMAP 
#include <R.h>
#include <Rdefines.h>

extern "C" SEXP matches(SEXP a, SEXP b)
{
  
  int alength = LENGTH(a);
  int blength = LENGTH(b);
  SEXP sortedA = PROTECT(Rf_allocVector(INTSXP,alength));
  int* apoint = INTEGER(sortedA);
  internalOrder(apoint,a);
  SEXP sortedB = PROTECT(Rf_allocVector(INTSXP,blength));
  int* bpoint = INTEGER(sortedB);
  internalOrder(bpoint,b);
  std::vector<int> indexsA;
  indexsA.reserve(alength);
  std::vector<int> indexsB;
  indexsB.reserve(blength);
  
  switch(TYPEOF(a))
  {
    case INTSXP:
      {
        int* astart = &INTEGER(a)[0];
        int* bstart = &INTEGER(b)[0];
        nmatch(astart,bstart,indexsA,indexsB,apoint,bpoint,alength,blength);
        break;
      }
    case REALSXP:
      {
        double* astart = &REAL(a)[0];
        double* bstart = &REAL(b)[0];
        nmatch(astart,bstart,indexsA,indexsB,apoint,bpoint,alength,blength);
        break;
      }
    case LGLSXP:
      {
        int* astart = &INTEGER(a)[0];
        int* bstart = &INTEGER(b)[0];
        nmatch(astart,bstart,indexsA,indexsB,apoint,bpoint,alength,blength);
        
        break;
      }
    case STRSXP:
      {
        SEXP* astart = &STRING_PTR(a)[0];
        SEXP* bstart = &STRING_PTR(b)[0];
        cmatch(astart,bstart,indexsA,indexsB,apoint,bpoint,alength,blength);
        break;
      }
    default:
      UNPROTECT(2);
    Rf_error("Unsupported type for matching.")
    ;
  }
  
  //Unfortunate overhead needed to convert vector to SEXP
  SEXP result1 = PROTECT(Rf_allocVector(INTSXP,indexsA.size()));
  SEXP result2 = PROTECT(Rf_allocVector(INTSXP,indexsB.size()));
  int* respoint1 = INTEGER(result1);
  int* respoint2 = INTEGER(result2);
  //Will this work?  If so, then we could probably dispense with allocating new vectors above (and instead allocate of size 0)
  //respoint1 = &indexsA.front();
  //respoint1 = &indexsB.front();
  std::copy(indexsA.begin(),indexsA.end(),respoint1);
  std::copy(indexsB.begin(),indexsB.end(),respoint2);
  SEXP combined = PROTECT(Rf_allocVector(VECSXP,2));
  SET_VECTOR_ELT(combined,0,result1);
  SET_VECTOR_ELT(combined,1,result2);
  UNPROTECT(5);
  return(combined);
  
}