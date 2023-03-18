#include "ptapso_omp.h"
#include "ptapsotestfunc.h"
#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include"omp.h"

// global pointer, used for updating velocity/position update of particles in openmp
// size_t *ringNbrs is for id of local ring, e.g. [1,2,3]
gsl_vector *accVecPbest;
gsl_vector *accVecLbest;
gsl_vector *chi1Vec;
gsl_vector *chi2Vec;
size_t *ringNbrs;
#pragma omp threadprivate(accVecPbest,accVecLbest,chi1Vec,chi2Vec,ringNbrs)
/* 2022/04/24*/


/*! \file
\brief Particle Swarm Optimization (PSO) and support functions.

The function \ref ptapso is the main one. The rest are support functions 
that do jobs such as initialization, memory allocation/deallocation, 
and handling output.

\author Soumya D. Mohanty
*/
/*! 
This function accepts a pointer to a fitness function and searches 
for its global optimum using Particle Swarm Optimization. The fitness
function must accept coordinates in the range [0,1]. See how_to_code_fitnessFunc.txt for
an example of the interface required for a fitness function.

Notes on the PSO implementation used:
   - Follows the prescription of Bratton, Kennedy, 2007.
   - Local best (lbest) PSO with three nearest neighbors in a ring topology.
   - Linearly deacreasing inertia weight.
   - Velocity clamping
*/
void ptapso_omp(size_t nDim, double (*fitfunc)(gsl_vector *, void *), void *ffParams,struct psoParamStruct *psoParams, struct returnData *psoResults)
//size_t nDim, /*!< Number of search dimensions */
//double (*fitfunc)(gsl_vector *, void *), /*!< Pointer to Fitness function */
//void *ffParams, /*!< Fitness function parameter structure */
//struct psoParamStruct *psoParams, /*!< PSO parameter structure */
//struct returnData *psoResults /*!< Output structure */)
{				
	/* 
	  Random numbers can be generated on the fly or read from a file. If a
	  valid file is specified, generation of random numbers is overriden.
	*/
	gsl_rng *rngGen = psoParams->rngGen;

	/* Initialize local minimizer of gbest */
	gsl_multimin_function func2minimz;
	func2minimz.n = nDim; /*dimensionality of function to minimize */
	func2minimz.f = dummyfitfunc;//fitfunc; /* Name of function to minimize */

	struct dummyFitFuncParam dffp;
        // give fitfunc to dffp.trufuncPr
	dffp.trufuncPr = fitfunc;
        // give ffParams to dffp.trufuncParam
	dffp.trufuncParam = ffParams;


	func2minimz.params = &dffp;//ffParams; /* Parameters needed by this function */
	/* Local Minimization method: Nelder Mead */
	gsl_multimin_fminimizer *minimzrState = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, nDim);
	/* Initial step vector of local minimization method */
	gsl_vector *locMinStp = gsl_vector_alloc(nDim);
	gsl_vector_set_all(locMinStp,psoParams->locMinStpSz);
	size_t lpLocMin;/* Local minimization iteration counter */
	int status;
	
	/* PSO loop counters */
	size_t lpParticles, lpPsoIter;
	/* Number of particles */
	const size_t popsize = psoParams->popsize;
	/* Number of iterations */
	const size_t maxSteps = psoParams->maxSteps;
	/* Information about a particles is stored in a struct array.
	*/

        struct particleInfo pop[popsize];

	/* initialize particles */
	for (lpParticles = 0; lpParticles < popsize; lpParticles++)
        {
            initPsoParticles(&pop[lpParticles], nDim, rngGen);
        }
	
	/* Variables needed to find and track gbest */
	double gbestFitVal = GSL_POSINF;
	gsl_vector *gbestCoord = gsl_vector_alloc(nDim);
	gsl_vector *partSnrCurrCol = gsl_vector_alloc(popsize);
	size_t bestfitParticle;
	double currBestFitVal;
	/* Variables needed to keep track of number of fitness function evaluations */
	unsigned char computeOK;
	size_t funcCount;


	/* Variables needed in PSO dynamical equation update */
	size_t lpNbrs; /* Loop counter over nearest neighbors */
	size_t nNbrs = 3;
	//size_t ringNbrs[nNbrs];
	double nbrFitVal; /*Fitness of a neighbor */
	size_t lbestPart; /* local best particle */
	double lbestFit; /* Fitness of local best particle */
	/* Variables needed for PSO dynamical equations */
	size_t lpCoord;
        //gsl_vector *local_ring_fitness_3 = gsl_vector_alloc (3); // 2022/04/24 use this gsl vector to get localbest of each local ring.

        // constant vector, c1,c2 for gsl_vector maultiplication
	//gsl_vector *accVecPbest = gsl_vector_alloc(nDim); // private for each particle, should be allocate within omp region  2022/04/24
	//gsl_vector *accVecLbest = gsl_vector_alloc(nDim);
	//gsl_vector *chi1Vec = gsl_vector_alloc(nDim);
	//gsl_vector *chi2Vec = gsl_vector_alloc(nDim);
	
	/* 
	   Start PSO iterations from the second iteration since the first is used
	   above for initialization.
	*/
        int tid;

        /* ----------------------------------------------------------------------------------------- */
        /* open matlab parpool  take  ~10 sec, open omp by #pragma omp parallel is cheap */
        /*open them 1e7 times by  for i-1:1:1e7  takes*/
        /*    #pragma omp parallel
              {
              }
        ~10 sec.
        */
        /*    #pragma omp parallel
              {
                  #pragma omp parallel for
                  for(int j=0;j<1;j++)
                  {
                  }
              }
         ~30 sec 
        */
        /*    #pragma omp parallel for
              for(int j=0;j<1;j++)
              {
              }
        ~ 15 sec
        */
        /*for a few thounsand pso steps, cost from init omp multiple times (within pso step loop) is ignorable, 1e-3 sec level.*/
        /* 2022/04/23 */
        /* ----------------------------------------------------------------------------------------- */
	for (lpPsoIter = 1; lpPsoIter <= maxSteps-1; lpPsoIter++)
        {
            //printf("...%zu/%zu-step ... \n",lpPsoIter,maxSteps-1);

            if (psoParams->debugDumpFile != NULL)
            {
                fprintf(psoParams->debugDumpFile,"Loop %zu \n",lpPsoIter);
                particleInfoDump(psoParams->debugDumpFile,pop,popsize);
            }        	
            /*----------------------------------------------------------------------------------------*/	
            /*----------------------------------------------------------------------------------------*/	
            /*----------------------------------------------------------------------------------------*/	
            /* Calculate fitness values */
            /* ZXB 0815 test openmp for Mohanty's pso*/
            /* 2022/03/05 FIXME: when they(different cpu) call the same fitnessfunction, is there any fight? */
            /* 2022/03/05 FIXME: when they(different cpu) call the same gsl_vector_set, is there any fight? */
            /* 2022/03/05 FIXME: when they(different cpu) share the same computeOK/funcCount, is there any fight? */
            // 2022/04/22 try fitness in omp //
            // pop[popsize] is a struct array,which should be shared by particles.
            // ffParams is a pointer,which should be shared by particles.
            /*----------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/* openmp region 1: evaluating fitness function of particles in parallel  */
/* -----------------------2022/04/23------------------------------------------*/
/*----------------------------------------------------------------------------*/
#pragma omp parallel private(lpParticles,tid,computeOK,funcCount) shared(pop,partSnrCurrCol,ffParams)
{
                #pragma omp for
                for (lpParticles = 0; lpParticles < popsize; lpParticles++)
                {
                    printf("\n");                    
                    tid=omp_get_thread_num(); // Obtain thread id
                    //printf("tid=%d, lpParticles=%zu\n", tid,lpParticles);

                    /* Evaluate fitness */
                    pop[lpParticles].partSnrCurr = fitfunc(pop[lpParticles].partCoord,ffParams);
                    /* Separately store all fitness values -- needed to find best particle */
                    gsl_vector_set(partSnrCurrCol,lpParticles,pop[lpParticles].partSnrCurr);
                    /* Check if fitness function was actually evaluated or not */
	            computeOK = ((struct fitFuncParams *)ffParams)->fitEvalFlag;
	            funcCount = 0;
	        
                    if (computeOK)
                    {
                        /* Increment fitness function evaluation count */
	                funcCount = 1;
                    }
	            pop[lpParticles].partFitEvals += funcCount;
                    /* Update pbest fitness and coordinates if needed */
	            if (pop[lpParticles].partSnrPbest > pop[lpParticles].partSnrCurr)
                    {
	                pop[lpParticles].partSnrPbest = pop[lpParticles].partSnrCurr;
                        //copies element, use pop[lpParticles].partCoord to replace pop[lpParticles].partPbest
	                gsl_vector_memcpy(pop[lpParticles].partPbest,pop[lpParticles].partCoord);
	            }
                } //end omp for
}// end of omp parallel region 1                   
            /*----------------------------------------------------------------------------------------*/	
            /*----------------------------------------------------------------------------------------*/	
            /*----------------------------------------------------------------------------------------*/					
            /* Find the best particle in the current iteration */
            bestfitParticle = gsl_vector_min_index(partSnrCurrCol);
	    currBestFitVal = pop[bestfitParticle].partSnrCurr; 
		
            if (gbestFitVal > currBestFitVal)
            {
             //Do local minimization iterations since gbest has changed.
		
            /* Merely calling the gsl_multimin_fiminizer_set function leads to nDim+1  function evaluations! */
                
            //gsl_multimin_fminimizer_set(minimzrState,&func2minimz,pop[bestfitParticle].partCoord,locMinStp);			
            //funcCount = nDim+1;
			
            //for (lpLocMin = 0; lpLocMin < psoParams->locMinIter; lpLocMin++)
            //  {
            //	status = gsl_multimin_fminimizer_iterate(minimzrState);
            //	  //A non-zero value of status indicates some type of failure  
            //	if (status)
            //		break;
            //	 /*Note that the function evaluation count is only an approximate
            //	   one for the nmsimplex2 algorithm as GSL routines 
            //	   do not return this information.*/
            //	funcCount += nDim+1;
            //	pop[bestfitParticle].partSnrCurr = gsl_multimin_fminimizer_minimum(minimzrState);
            //      gsl_vector_memcpy(pop[bestfitParticle].partCoord, minimzrState->x);
            // }
            //pop[bestfitParticle].partFitEvals += funcCount;
                
            /* Update particle pbest */
            pop[bestfitParticle].partSnrPbest = pop[bestfitParticle].partSnrCurr;
            gsl_vector_memcpy(pop[bestfitParticle].partPbest,pop[bestfitParticle].partCoord);
            /* Update gbest */
            gbestFitVal = pop[bestfitParticle].partSnrCurr;
            gsl_vector_memcpy(gbestCoord,pop[bestfitParticle].partCoord);
            /* Update list of fitness Values */
            gsl_vector_set(partSnrCurrCol,bestfitParticle,pop[bestfitParticle].partSnrCurr);
            }

/*----------------------------------------------------------------------------------*/
/* openmp region 2: update localbest by 3-ring structure for particles in parallel  */
/* -----------------------2022/04/23------------------------------------------------*/
/*----------------------------------------------------------------------------------*/
// nNbrs=3, shared
// ringNbrs[3], fixed-size array, private
// lbestPart, id of local ring best
// lbestFit, fitness of local ring best
// partSnrCurrCol,vector of current snr, size=popsize
// nbrFitVal, fitness vaue of ringNbrs[1] and ringNbrs[2]
// lpNbrs  =1:3
//
//double t0,t1,t2;
#pragma omp parallel private(lpParticles,lbestPart,lbestFit,lpNbrs,nbrFitVal) shared(nNbrs,partSnrCurrCol,pop)
{	
            /* Get lbest updated, using ring topology */
            #pragma omp for
            for (lpParticles = 0; lpParticles < popsize; lpParticles++)
            {
                //size_t ringNbrs[nNbrs]; // ringNbrs is a private array, need allocation for wach particle, same as TDI temp_h, 2022/04/24
                ringNbrs = (size_t *)calloc(nNbrs,sizeof(size_t));  
 
                if (lpParticles == 0)
                {
                    ringNbrs[0]=popsize-1; ringNbrs[1]=lpParticles; ringNbrs[2]=lpParticles+1;
                }
                else if (lpParticles == popsize -1)
                {
                    ringNbrs[0]=lpParticles-1; ringNbrs[1]=lpParticles; ringNbrs[2]=0;
                }
                else
                {
                    ringNbrs[0]=lpParticles-1; ringNbrs[1]=lpParticles; ringNbrs[2]=lpParticles+1;
                }					  
                /* Get best particle in neighborhood */
                // use ringNbrs[0] first
                lbestPart = ringNbrs[0];
                lbestFit = gsl_vector_get(partSnrCurrCol,lbestPart);

               // t0 = gsl_vector_get(partSnrCurrCol,ringNbrs[0]);
               // t1 = gsl_vector_get(partSnrCurrCol,ringNbrs[1]);
               // t2 = gsl_vector_get(partSnrCurrCol,ringNbrs[2]);

               // gsl_vector_set(local_ring_fitness_3,0,t0);
               // gsl_vector_set(local_ring_fitness_3,1,t1);
               // gsl_vector_set(local_ring_fitness_3,2,t2);

               // lbestFit  = gsl_vector_min(local_ring_fitness_3);
               // lbestPart = gsl_vector_min_index(local_ring_fitness_3);

               // compare ringNbrs[1] and ringNbrs[2] with ringNbrs[0], get local ring best
                for (lpNbrs = 1; lpNbrs < nNbrs; lpNbrs++)
                {
                    nbrFitVal = gsl_vector_get(partSnrCurrCol,ringNbrs[lpNbrs]);
                    if (nbrFitVal < lbestFit)
                    {
                        lbestPart = ringNbrs[lpNbrs];
                        lbestFit = nbrFitVal;
                    }
                }
                //compare local ring best woth pbest
                if (lbestFit < pop[lpParticles].partSnrLbest)
                {
                    pop[lpParticles].partSnrLbest = lbestFit;
                    gsl_vector_memcpy(pop[lpParticles].partLocalBest,pop[lbestPart].partCoord);
                }

                // free ptivate array ringNbrs[3]
                free(ringNbrs);
            } //end omp for
}// end of omp parallel region 2     
        
/*----------------------------------------------------------------------------*/
/* openmp region 3: update velocity,position  for particles in parallel  */
/* -----------------------2022/04/23------------------------------------------*/
/*----------------------------------------------------------------------------*/
// accVecPbest, private vector 1, pbest with size=nDim
// accVecLbest, private vector 2, pbest with size=nDim
// chi1Vec,private vector 3, randon coefficients with size=nDim, it is r1 that we get used to vector multiplication as (x-x_lbest_local_ring)*r1.
// chi2Vec,private vector 4, randon coefficients with size=nDim, it is r2 that we get used to vector multiplication as (x-x_pbest)*r2.
// lpCoord=1:nDim
// rngGen is a pointer, gsl_rng *rngGen = psoParams->rngGen;
//
// note the private vector need allocation/de-alloaction within the parallel region like temp_h in TDI omp. 2022/04/24
//
#pragma omp parallel private(lpParticles,lpCoord) shared(nDim,pop,psoParams,lpPsoIter,rngGen)
{
            #pragma omp for
            for (lpParticles = 0; lpParticles < popsize; lpParticles++)
            {
                /* 2022/04/24 allocation for private vector*/
	        accVecPbest = gsl_vector_alloc(nDim);
	        accVecLbest = gsl_vector_alloc(nDim);
	        chi1Vec = gsl_vector_alloc(nDim);
	        chi2Vec = gsl_vector_alloc(nDim);

                /* Update inertia Weight */
                pop[lpParticles].partInertia = psoParams->dcLaw_a-(psoParams->dcLaw_b/psoParams->dcLaw_c)*lpPsoIter;
                if (pop[lpParticles].partInertia < psoParams->dcLaw_d)
                {
                    pop[lpParticles].partInertia = psoParams->dcLaw_d;
                }
                /* Acceleration vector to pbest */
                gsl_vector_memcpy(accVecPbest,pop[lpParticles].partPbest);
                gsl_vector_sub(accVecPbest,pop[lpParticles].partCoord);
                /* Acceleration vector to lbest */
                gsl_vector_memcpy(accVecLbest,pop[lpParticles].partLocalBest);
                gsl_vector_sub(accVecLbest,pop[lpParticles].partCoord); //x-x_p_local_ring_best
                /* Random weights for acceleration components */
                for (lpCoord = 0; lpCoord < nDim; lpCoord++)
                {
                    gsl_vector_set(chi1Vec,lpCoord,gsl_rng_uniform(rngGen));
                }
                for (lpCoord = 0; lpCoord < nDim; lpCoord++)
                {
                    gsl_vector_set(chi2Vec,lpCoord,gsl_rng_uniform(rngGen));
                }
                /* Multiply random weights */
                gsl_vector_mul(accVecPbest,chi1Vec);
                gsl_vector_mul(accVecLbest,chi2Vec);
                /* Scale with acceleration constants */
                gsl_vector_scale(accVecPbest,psoParams->c1);
                gsl_vector_scale(accVecLbest,psoParams->c2);
                /* Velocity update 			
                pop(k,partVelCols)=partInertia*pop(k,partVelCols)+...
                                   c1*(pop(k,partPbestCols)-pop(k,partCoordCols))*chi1+...
                                   c2*(pop(k,partLocalBestCols)-pop(k,partCoordCols))*chi2;
                */
                gsl_vector_scale(pop[lpParticles].partVel,pop[lpParticles].partInertia);
                gsl_vector_add(pop[lpParticles].partVel,accVecPbest);
                gsl_vector_add(pop[lpParticles].partVel,accVecLbest);
                /*
			  Apply max. velocity threshold
		        maxvBustCompPos = find(pop(k,partVelCols) > max_velocity);
		        maxvBustCompNeg = find(pop(k,partVelCols) < -max_velocity);
		        if ~isempty(maxvBustCompPos)
		            pop(k,partVelCols(maxvBustCompPos))= max_velocity;
		        end
		        if ~isempty(maxvBustCompNeg)
		            pop(k,partVelCols(maxvBustCompNeg))= -max_velocity(1);
		        end
                */
               limitVecComponent(pop[lpParticles].partVel, -psoParams->max_velocity, psoParams->max_velocity);
               /*Position Update*/
                gsl_vector_add(pop[lpParticles].partCoord,pop[lpParticles].partVel); 


                /* 2022/04/24 de-allocation for private vector*/      
	        gsl_vector_free(accVecPbest);
                gsl_vector_free(accVecLbest); 
	        gsl_vector_free(chi1Vec);
	        gsl_vector_free(chi2Vec); 
	    } //end omp for
}// end of omp parallel region 3
		
            if (psoParams->debugDumpFile != NULL)
            {
                fprintf(psoParams->debugDumpFile,"After dynamical update\n");   
                particleInfoDump(psoParams->debugDumpFile,pop,popsize);
                fprintf(psoParams->debugDumpFile,"--------\n");			      
            }
            printf("\n");
        } // end of pso steps


	/* Prepare output */
	psoResults->totalIterations = lpPsoIter-1;
	/* 	actualEvaluations = sum(pop(:,partFitEvalsCols)); */
	psoResults->totalFuncEvals = 0;
	for (lpParticles = 0; lpParticles < popsize; lpParticles ++)
        {
		psoResults->totalFuncEvals += pop[lpParticles].partFitEvals;
	}
	gsl_vector_memcpy(psoResults->bestLocation, gbestCoord);
	psoResults->bestFitVal = gbestFitVal;
	
	/* Free function minimizer state */
	gsl_multimin_fminimizer_free(minimzrState);
	/* Deallocate vectors */
	gsl_vector_free(locMinStp);
	gsl_vector_free(gbestCoord);
	gsl_vector_free(partSnrCurrCol);
	//gsl_vector_free(accVecPbest);
        //gsl_vector_free(accVecLbest); 
	//gsl_vector_free(chi1Vec);
	//gsl_vector_free(chi2Vec); 
	//gsl_vector_free(local_ring_fitness_3);//2022/04/24 use gsl vector to get localbest of the local ring 
	/* Deallocate members of pop */
	for(lpParticles = 0; lpParticles < popsize; lpParticles++)
        {
		particleinfo_free(&pop[lpParticles]);
	}
} //end of ptapso_omp function

/*Dummy function to wrap supplied function so that the 
  call to the GSL local optimizer is compatible */
double dummyfitfunc(const gsl_vector *xVec, void *dffParams)
{
	struct dummyFitFuncParam *dfp = (struct dummyFitFuncParam *)dffParams;
	double (*trueFitFunc)(gsl_vector *, void *) = dfp->trufuncPr;
	gsl_vector *xVec2 = (gsl_vector *)xVec;
	double funcVal = trueFitFunc(xVec2,dfp->trufuncParam);
}

/*! Initializer of particle position, velocity, and other properties. */
void initPsoParticles(struct particleInfo *p, size_t nDim, gsl_rng *rngGen)
{
	
	double rngNum;
	size_t lpCoord;
	
	particleinfo_alloc(p,nDim);
	
	for (lpCoord = 0; lpCoord < nDim; lpCoord++){
		rngNum = gsl_rng_uniform(rngGen);
		gsl_vector_set(p->partCoord,lpCoord,rngNum);
	}
	
	for (lpCoord = 0; lpCoord < nDim; lpCoord++){
		rngNum = gsl_rng_uniform(rngGen);
		rngNum = - gsl_vector_get(p->partCoord,lpCoord)
			     + rngNum;
		gsl_vector_set(p->partVel,lpCoord,rngNum);
	}
	
	if(gsl_vector_memcpy(p->partPbest, p->partCoord))
		printf("Error in copying vectors\n");	
	
	p->partSnrPbest = GSL_POSINF;
	p->partSnrCurr = 0;
	p->partSnrLbest = GSL_POSINF;
	p->partInertia = 0;
	p->partFitEvals = 0;
}


/*! Allocate storage for members of particleInfo structure */
void particleinfo_alloc(struct particleInfo *p, size_t nDim){
	p->partCoord = gsl_vector_alloc(nDim); /* Current coordinates */
	p->partVel = gsl_vector_alloc(nDim);  /* Current velocity */
	p->partPbest = gsl_vector_alloc(nDim); /* Coordinates of pbest */
	p->partLocalBest = gsl_vector_alloc(nDim); /* Coordinates of neighborhood best */
}

/*! Free the storage assigned to members of particleInfo structure */
void particleinfo_free(struct particleInfo *p){
	gsl_vector_free(p->partCoord);
	gsl_vector_free(p->partVel);
	gsl_vector_free(p->partPbest);
	gsl_vector_free(p->partLocalBest);	
}

/*! Allocate storage for returnData struct members */
struct returnData * returnData_alloc(size_t nDim){
	struct returnData *psoResults = (struct returnData *)malloc(sizeof(struct returnData));
	gsl_vector *bestLocation = gsl_vector_alloc(nDim);
	psoResults->bestLocation = bestLocation;
	return psoResults;
};

/*! Free storage assigned to returnData struct members */
void returnData_free(struct returnData *psoResults){
	gsl_vector_free(psoResults->bestLocation);
	free(psoResults);
};


/*! Dump information stored in particleInfo struct */
void particleinfo_fwrite(FILE *outF, struct particleInfo *p){
	
	size_t nDim = p->partCoord->size;
	size_t lpc;
	
	fprintf(outF,"Particle locations in standardized coordinates\n");
	for (lpc = 0; lpc < nDim; lpc++){
		fprintf(outF,"%f ",gsl_vector_get(p->partCoord,lpc));
	}
	fprintf(outF,"\n");
	fprintf(outF,"Particle velocities in standardized coordinates\n");
	for (lpc = 0; lpc < nDim; lpc++){
		fprintf(outF,"%f ",gsl_vector_get(p->partVel,lpc));
	}
	fprintf(outF,"\n -------- \n");
}

/*! Dump particleInfo struct array information as a matrix
with all information pertaining to one particle in a row.
*/
void particleInfoDump(FILE *outF, struct particleInfo *p, size_t popsize){	
	size_t nDim = p[0].partCoord->size;
	size_t lpParticles, lpCoord;
	
	for (lpParticles = 0; lpParticles < popsize; lpParticles++){
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partCoord,lpCoord));
		}
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partVel,lpCoord));
		}
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partPbest,lpCoord));
		}
		fprintf(outF,"%lf ",p[lpParticles].partSnrPbest);
		fprintf(outF,"%lf ",p[lpParticles].partSnrCurr);	
		fprintf(outF,"%lf ",p[lpParticles].partSnrLbest);
		fprintf(outF,"%lf ",p[lpParticles].partInertia);		
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partLocalBest,lpCoord));
		}
		fprintf(outF,"X ");	
		fprintf(outF,"%zu ",p[lpParticles].partFitEvals);		
		fprintf(outF,"\n");				
	}
}



