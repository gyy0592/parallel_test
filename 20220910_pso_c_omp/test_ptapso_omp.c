/* test function for ptapso.c */

#include <stdio.h>
#include "ptapso_omp.h"
/* The header file tells us what fitness function we are calling
and what the parameter structure for this function is.
*/
#include "ptapsotestfunc.h"
#include <omp.h>
#include<time.h>

int main(){
	/* Test set up */
	unsigned int nDim = 3, lpc;
	double rmin[3] = {-10,-10,-10};
	double rmax[3] = {10,  10,  10};
	double rangeVec[3];
	
	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();
	
	/* Initialize random number generator */
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(rngGen,2571971);
	
    /* Allocate fitness function parameter struct. 
	 */
	struct fitFuncParams *inParams = ffparam_alloc(nDim);
	
	/* Load fitness function parameter struct */
	for (lpc = 0; lpc < nDim; lpc++)
        {
		rangeVec[lpc]=rmax[lpc]-rmin[lpc];
		gsl_vector_set(inParams->rmin,lpc,rmin[lpc]);
		gsl_vector_set(inParams->rangeVec,lpc,rangeVec[lpc]);		
	}
	/* Set up pointer to fitness function. Use the prototype
	declaration given in the header file for the fitness function. */
	double (*fitfunc)(gsl_vector *, void *) = ptapsotestfunc;
	/* Set up special parameters, if any, needed by the fitness function used.
	   These should be provided in a structure that should be defined in 
	   the fitness function's header file.
	 */
	struct ptapsotestfunc_params splParams;
	splParams.dummyParam = 0;
        splParams.aa = (double *)calloc(6, sizeof(double));
        for(int i=0; i<6; i++)
        {
            splParams.aa[i] = (double)i;
        }
	/* Pass on the special parameters through the generic fitness function parameter
	struct */
	inParams->splParams = &splParams;
	
	/* Set up storage for output from ptapso. */
	struct returnData *psoResults = returnData_alloc(nDim);
	
	/* Set up the pso parameter structure.*/
	struct psoParamStruct psoParams;
	psoParams.popsize=40;
	psoParams.maxSteps= 100; 
	psoParams.c1=2;
	psoParams.c2=2;
	psoParams.max_velocity = 0.2;
	psoParams.dcLaw_a = 0.9;
	psoParams.dcLaw_b = 0.4;
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = 0.2;
	psoParams.locMinIter = 0;
	psoParams.locMinStpSz = 0.01;
	psoParams.rngGen = rngGen;
	psoParams.debugDumpFile = fopen("test_ptapso_Dump.txt","w");


	/* Run PSO multiple times to check for memory leaks */
        /*
	double start = omp_get_wtime();
        for(lpc = 0; lpc < 1; lpc++){
	printf("Loop # %zu\n", lpc);
		ptapso_omp(nDim, fitfunc, inParams, &psoParams, psoResults);	
	}
        double end = omp_get_wtime();  
        */
        /*ZXB: we need to run pso lpc times at the same time for global bestfit.
        assuming probability that one pso run catching the global bestfit is p (e.g.  p ~ 0.5),
        then probability that all lp pso runs missing the global bestfit is (1-p)^lp,
        probability that at least one pso run catching the global bestfit is 1-(1-p)^lp,
        1-0.5^10 is close to 1
        */


        /* Call PSO*/
        double start = omp_get_wtime();
        ptapso_omp(nDim, fitfunc, inParams, &psoParams, psoResults);
        double end = omp_get_wtime(); 
        printf("\n");
        printf("using pso_openmp, elpase time is %f sec  by omp_get_wtime.\n", end-start);
        printf("\n");
	fclose(psoParams.debugDumpFile);

	/* Information returned by PSO */
	printf("Total number of iterations %zu\n", psoResults->totalIterations);
	printf("Total number of function evaluations %zu\n", psoResults->totalFuncEvals);
	printf("Best Location found: \n");
	for (lpc = 0; lpc < nDim; lpc++){
		printf("%f, ",rmin[lpc]+gsl_vector_get(psoResults->bestLocation,lpc)*(rmax[lpc]-rmin[lpc]));		
	}
	printf("\n");
	printf("Best Fitness Value: %f\n", psoResults->bestFitVal);
        printf("\n");
	
	/* Free allocated memory */
	ffparam_free(inParams);
	returnData_free(psoResults);
	gsl_rng_free(rngGen);
        free(splParams.aa);
}



