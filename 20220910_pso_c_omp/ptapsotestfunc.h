/*! 
\file ptapsotestfunc.h

\brief This is an example of the header file that is required for a fitness function.
*/
#if !defined(PTAPSOTESTFUNCHDR)
#define PTAPSOTESTFUNCHDR
#include <gsl/gsl_vector.h>

//! [Special Params struct] 
/*! 
\brief Struct to pass fitness function specific parameters.

   If a fitness function foo.c needs special parameters, define them through a struct in its header file foo.h. 
   The fitness function here does not need special parameters but we define a 
   dummy one as an example. Special parameters are provided as fields of a struct. The name of
this struct should be unique and indicate the fitness function with which it is associated.
 This struct should be initialized before a call to the fitness function and a pointer to it should be assigned to the splParams field
   of the initialized fitFuncParams struct.
*/
struct ptapsotestfunc_params
{
    int dummyParam; /*!< Just a dummy special parameter */
    double *aa;   // vector in the model　　
};
//! [Special Params struct] 

/*! \brief A test fitness function.

All fitness functions must have the same input and output arguments.
*/



/*! [Declaration of fitness function] */
double ptapsotestfunc(gsl_vector *, void *);


/*! 
A fitness function must use this
struct to ferry in parameter values needed to compute the fitness value, and
ferry out additional outputs, if any, from the fitness function besides the fitness value.
A fitness function can use the parameters supplied in this structure to translate
standardized input coordinates, each of which belongs to [0,1], into their real values.
*/
struct fitFuncParams
{
    size_t nDim;    /*!< Dimensionality of the fitness function. */
    gsl_vector *rmin;     /*!< Minimum value of each coordinate. */
    gsl_vector *rangeVec; /*!< Range of each coordinate. */
    gsl_vector *realCoord;/*!< The unstandardized value of each coordinate is returned in this vector.*/
    unsigned char fitEvalFlag; /*!< Set to 0 if fitness is infinity, else to 1.*/

/*! Pointer to a struct that carries additional parameters that are 
specific to a fitness function.
The header file of a fitness function can define this struct.
*/
void *splParams;
};


struct fitFuncParams * ffparam_alloc(size_t );

void ffparam_free(struct fitFuncParams *);

size_t chkstdsrchrng(const gsl_vector *);

void s2rvector(const gsl_vector *, const gsl_vector *, const gsl_vector *, gsl_vector *);

void limitVecComponent(gsl_vector *, double, double);

char **listfileswext(const char *, const char *, size_t *, size_t *);
#endif
