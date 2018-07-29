#include <R_ext/Rdynload.h>

#include <Rinternals.h>
#include <R_ext/Random.h>
#include <R_ext/Print.h>

#include "fastlib/fastlib.h"
#include "fscm.h"


extern "C"
void do_fscm(const char **datafile, 
             const char **neighborfile, 
             const char **clustfile,
             const char **phifile,
             const char **psifile,
             const int *c,
             const double *tol_em,
             const int *max_iter_gamma,
             const int *max_iter_z,
             const int *max_iter_em,
             const int *max_iter_hmrf,
             const int *verbose,
	         int *ok) {
             
	Dataset dataset,neighborset,clustset,phiset, psiset;

	*ok = 0;

	if (!PASSED(dataset.InitFromFile(*datafile))) {    
		
		//fprintf(stderr, "main: Couldn't open file '%s'.\n", *datafile);    
        error("fscm: couldnt open datafile");
		return;

	}

	if (!PASSED(neighborset.InitFromFile(*neighborfile))) {    
		
		//fprintf(stderr, "main: Couldn't open file '%s'.\n", *neighborfile);    
        error("fscm: couldnt open neighborfile");
		return;

	}

	if (!PASSED(clustset.InitFromFile(*clustfile))) {    
		
		//fprintf(stderr, "main: Couldn't open file '%s'.\n", *clustfile);    
        error("fscm: couldnt open clustfile");
		return;

	}

	if (!PASSED(phiset.InitFromFile(*phifile))) {    
	
		//fprintf(stderr, "main: Couldn't open file '%s'.\n", *phifile);    
        error("fscm: couldn't open phifile");
		return;
	
	}


	if (!PASSED(psiset.InitFromFile(*psifile))) {    
		
		//fprintf(stderr, "main: Couldn't open file '%s'.\n", *psifile);    
        error("fscm: couldn't open psifile");
		return;

	}

    GetRNGstate();


    FSCM fit;
   	Vector alpha,gamma,sigma;
	Matrix beta,pi,pi_y,gamma_all,Vgamma;

    fit.verbose = *verbose;




	fit.InitTrain(dataset,neighborset,clustset,phiset,psiset,alpha,beta,sigma,gamma,Vgamma,pi,pi_y,*c,*max_iter_hmrf);

    int iter=0, error=0;

    while(iter < *max_iter_em) {
        if(*verbose) Rprintf("3. Update f(Z|Y;theta) \n");
		fit.HMRF(alpha,beta,sigma,gamma,Vgamma,pi,pi_y,*max_iter_gamma); //given alpha, beta, sigma, pi, (gamma, Vgamma); update pi_y

		if(hasnan(fit.psi)) {
			if(*verbose) Rprintf("nan occured during the estimation (psi)\n");
			error=1;
			break;
		}
	
		if(*verbose) Rprintf("2. Update f(gamma|Y;sigma_s)\n");
		fit.RandomEffect(alpha,beta,pi,sigma,gamma,Vgamma,*max_iter_z); //given alpha, beta, sigma, pi; update gamma, Vgamma;
	
		
		if(*verbose) Rprintf("1. Update f(Y|gamma,z; alpha,beta,sigma_epsilon) \n");
		fit.MStep(gamma, Vgamma, pi_y, alpha, beta , sigma, pi, *max_iter_hmrf); //given gamma, Vgamma,pi_y; update alpha, beta, sigma, pi

		if(yet_another_isnan(sigma[0]) || yet_another_isnan(sigma[1])) {
			if(*verbose) Rprintf("nan occured during the estimation (sigma)\n");
			error=1;
			break;
		}

		//printf("4. Calculate the improvement of log-likelihood\n");
		//fit.LogLikelihood(alpha, beta, sigma, pi, log_likelihood_new); //given alpha, beta, sigma, pi,theta, update log_likelihood_new

		iter++;

    
    }

	if(!error) {
		if (*verbose) Rprintf("5. Save the results\n");

		//Matrix probsAll;
		//probsAll.Copy(pi);
		data::Save("probsAll.csv", pi);


		GenVector<fl__index_t> index;
		index.Init(2);
		for (fl__index_t i = 0; i < 2; i++)
			index[i] = i;
		data::Save("parameter.csv", index, sigma);

		fit.Results(alpha, beta, gamma, pi_y);
		*ok = 1;
	}

    PutRNGstate();

}

static R_NativePrimitiveArgType do_fscm_t[] = {
	STRSXP, STRSXP, STRSXP, STRSXP, STRSXP, 
	INTSXP, REALSXP, 
	INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_CMethodDef cMethods[] = {
	{"do_fscm", (DL_FUNC) &do_fscm, 13, do_fscm_t},
	{NULL, NULL, 0}
};

extern "C"
void R_init_funcy(DllInfo *info) {
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
}

