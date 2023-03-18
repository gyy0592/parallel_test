/*! \file ptapsotestfunc.c
\brief Example of a fitness function.
 */
#include <stdio.h>
#include <string.h>
#include "ptapso_omp.h"
#include "ptapsotestfunc.h"
//#include "maxphase.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <fnmatch.h>

/*
fitFuncVal = ptapsotestfunc_gsl(xVec,P)
A benchmark test function for PTAPSO
that computes the Rastrigin fitness function for
each row of xVec.  The fitness value is returned in fitFuncVal.
xVec is standardized, that is 0<=xVec(i,j)<=1. 
The values used to convert xVec(i,j)
internally before computing fitness are given in P.rmin and
P.rangeVec: 
xVec(j) -> xVec(j)*rangevec(j)+rmin(j).
fitFuncVal = infty if the point xVec falls
outside the hypercube defined by 0<=xVec(j)<=1.
The real coordinates are returned in P.realCoord. 

Soumya D. Mohanty, Jan 2016
- Derived from ptapsotestfunc.c. Converts to gsl_vector inputs and 
  uses the interface needed by GSL multi-dimensional local minimization routines.
*/



double ptapsotestfunc(gsl_vector *xVec, void  *inParamsPointer){
	
	unsigned int validPt;
        unsigned int lpc;
	//! [Cast fit func params]
	struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;
	//! [Cast fit func params]
	unsigned int ncols = inParams->nDim;
	double rangeVec, rmin, x;
	double fitFuncVal;
	//! [Cast special params]
	/* Shows how to retrieve special parameters (dummy ones are supplied for this particular
	fitness function).
	*/
	struct ptapsotestfunc_params *splParams = (struct ptapsotestfunc_params *)inParams->splParams;
		
	/* This fitness function knows what fields are given in the special parameters struct */
        int dummy = splParams -> dummyParam;
        double *aaa  = splParams -> aa;
	//printf("PTAPSOTESTFUNC called with special parameters %d\n",dummy);
	//! [Cast special params]
	
	s2rvector(xVec,inParams->rmin,inParams->rangeVec,inParams->realCoord);
	
	validPt = chkstdsrchrng(xVec);
	
	if (validPt)
        {
		inParams->fitEvalFlag = 1;
		fitFuncVal = 0;
//
//give your fitness function here, do not change the others.
//
//
//#---------------------------------fitness  function example 1-------------------------------------------------------
//		for (lpc = 0; lpc < ncols; lpc++)
//              {
//			x = gsl_vector_get(inParams->realCoord,lpc);
//			fitFuncVal = fitFuncVal + gsl_pow_int(x,2) - 10.0*cos(2*M_PI*x) + 10;
//              }
//#---------------------------------fitness  function example 2----------------------------------------------------------------------

// set up your parameters
double x0 = gsl_vector_get(inParams->realCoord,0);
double x1 = gsl_vector_get(inParams->realCoord,1);
double x2 = gsl_vector_get(inParams->realCoord,2);

double sum=0.0;
// set up your fitness functions (pso search the minimum)
for(int i=0;i<100000;i++)
{
    sum = sum + pow(x0-aaa[0]-aaa[1],2.0) + pow(x1-aaa[2]-aaa[3],2.0) + pow(x2-aaa[4]-aaa[5],2.0);  
    //fitFuncVal = (x0-1.0)*(x0-1.0) + (x1-2.0)*(x1-2.0) + (x2-3.0)*(x2-3.0);  
}
//printf("x0=%f  x1=%f  x2=%f  fitFuncVal=%f\n",x0,x1,x2,fitFuncVal);
//fitFuncVal = pow(x0-1.0,2.0) + pow(x1-2.0,2.0); 
fitFuncVal = sum;
//#--------------------------------------------------------------------------------------------------------
        }
	else{
		fitFuncVal=GSL_POSINF;
		inParams->fitEvalFlag = 0;
	}
   return fitFuncVal;
}

/*! 
Allocate a fitness function parameter structure for a given dimensionality of the fitness function.
*/
//ffparam_alloc is a struct pointer.
// double *a = (double*)calloc(Nl,sizeof(double));   
struct fitFuncParams *ffparam_alloc(size_t nDim)
{
	struct fitFuncParams *ffp = (struct fitFuncParams *)malloc(sizeof(struct fitFuncParams));
	ffp->nDim = nDim;
	ffp->rmin = gsl_vector_calloc(nDim);
	ffp->rangeVec = gsl_vector_calloc(nDim);
	ffp->realCoord = gsl_vector_calloc(nDim);
	return ffp;
}

/*!
Deallocate a fitness function parameter structure.
Any special parameter structure contained in the 
fitness function parameter structure must be freed 
before calling this function.
*/
void ffparam_free(struct fitFuncParams *ffp)
{
	if (ffp == NULL)
        {
		printf("Invalid FitFuncParams supplied\n");
		return;
	}
	gsl_vector_free(ffp->rmin);
	gsl_vector_free(ffp->rangeVec);
	gsl_vector_free(ffp->realCoord);
	free(ffp); // free struct pointer
}

/*!					   			   
Takes standardized coordinates and
returns real coordinates using the supplied range and minimum limits. 
The range and limits can be different for the different coordinates. 

Notes:
 -  Derived from sr2vector.c
 -  Shifted to using gsl_vector
*/
void s2rvector(const gsl_vector *xVec, /*!< Standardized coordinates of the point.*/
 			   const gsl_vector *rmin,       /*!< Minimum value of each coordinate.*/
			   const gsl_vector *rangeVec,   /*!< Range of each coordinate. */
               gsl_vector *realCoord  /*!< Returns real coordinate values.*/)
{
			   
	size_t lpc;
	double x, rv, rmi, rc;
	size_t ncols = xVec->size;
	
    for(lpc=0; lpc < ncols; lpc++)
    {
		x = gsl_vector_get(xVec,lpc);
		rv = gsl_vector_get(rangeVec,lpc);
		if(rv<0)
                {
			printf("Invalid range\n");
			abort();
		}
		rmi = gsl_vector_get(rmin,lpc);
		rc = x*rv+rmi;
    	gsl_vector_set(realCoord,lpc,rc);
    }
}

/*! 
Returns 0 or 1 corresponding to invalid/valid
coordinates. The invalid flag is set if any of the coordinates
 fall outside the closed interval [0,1].

Notes:
 - Derived from chkstdsrchrng.c
*/
size_t chkstdsrchrng(const gsl_vector *xVec/*!< Standardized coordinates of the point.*/)
{
    double x;
    unsigned int validPt,lpc;
    size_t ncols = xVec->size;;
	
    /*printf("from chkstdsrchrng\n");*/
    validPt=1;
    for (lpc=0; lpc <ncols; lpc++)
    {
        x = gsl_vector_get(xVec,lpc);
        if (x<0 || x>1)
        {
            validPt = 0;
            break;
        }						
    }
    return validPt;
}

/*! Limit each compononent of a vector to a specified range */
void limitVecComponent(gsl_vector *xVec, double min, double max)
{
	size_t lpc;
	size_t nDim = xVec->size;
	double x;
	for (lpc = 0; lpc < nDim; lpc++)
        {
		x = gsl_vector_get(xVec,lpc);
		if (x < min)
		    gsl_vector_set(xVec,lpc,min);
		else if (x > max)
		    gsl_vector_set(xVec,lpc,max);
	}
}

/*! List all files in a directory with a given extension. 
   - First parameter: Extension (without a '.'). Example: "mat" and not ".mat".
   - Second parameter: Directory name
   - Third parameter:  number of files found.
   - Fourth parameter: length of longest file name string.
   - Output: List of filenames. Dynamically allocated memory: must be freed by 
           calling function.
 */
char ** listfileswext (const char *ext, const char *dirName, size_t *nFiles, size_t *maxFileNameLen)
{
	/*printf("------- %s -------\n",ext);
	printf("Extension string length %zu\n",strlen(ext));*/
  char *pattern = (char *)malloc((3+strlen(ext))*sizeof(char));
  DIR *dp;
  struct dirent *ep;
  char **fileList;
  size_t lpc;
 
  pattern[0] = '*';
  pattern[1] = '.';
  pattern[2] = '\0';
  strcat(pattern,ext);
  
  /*printf("Pattern %s\n",pattern);*/
  /* Step 1: Count the number of files with the required 
     extension.
  */
  size_t countValidFiles = 0;
  dp = opendir(dirName);
  if (dp != NULL)
  {
      while ((ep = readdir(dp)))
      {
	   /*printf("Checking file %s\n",ep->d_name);*/
	   //if(!fnmatch(pattern, ep->d_name, (FNM_FILE_NAME|FNM_PERIOD))){
          if(!fnmatch(pattern, ep->d_name, (FNM_PATHNAME|FNM_PERIOD)))
          {
              /* Match. Increment counter */
              /*printf("found match with pattern %s\n",pattern);*/
              countValidFiles++;
          }
      }
      (void) closedir (dp);
      /*Apparently, there is no need to free ep as it is declared to be 'static' in the readdir function */
  }
  else
  {
    printf ("Couldn't open the directory %s\n",dirName);
    free(pattern);
    return NULL;
  }
  *nFiles = countValidFiles;
  /* Step 2: Create storage for list of filenames */
  fileList = (char **)malloc(countValidFiles*sizeof(char *));
  /* Can't find a better way than to repeat the whole loop again */
  countValidFiles = 0;
  dp = opendir(dirName);
  if (dp != NULL)
  {
      while ((ep = readdir(dp)))
      {
          if(!fnmatch(pattern, ep->d_name, (FNM_PATHNAME|FNM_PERIOD)))
          {
              fileList[countValidFiles] = (char *)malloc((strlen(ep->d_name)+1)*sizeof(char));
              strcpy(fileList[countValidFiles],ep->d_name);
              /* Match. Increment counter */
              countValidFiles++;
          }
      }
      (void) closedir (dp);
  }
  else
  {
    printf ("Couldn't open the directory %s\n",dirName);
    return NULL;
  }
  
	/*Find longest filename */
	size_t fileNameLen;
	*maxFileNameLen = 0;
	size_t lpc1;
	for (lpc1 = 0; lpc1 < *nFiles; lpc1++)
        {
		fileNameLen = strlen(fileList[lpc1]);
		if ( fileNameLen > *maxFileNameLen)
			*maxFileNameLen = fileNameLen;
	}
  
  /* Wrap up */
  free(pattern);
  return fileList;
}


