#ifndef U_FSCM_FSCM
#define U_FSCM_FSCM

#include "fastlib/fastlib.h"

class FSCM{
    
    public: //cannot be accessed by functions outside FSCM

        int verbose;
    
        Matrix y,neighbors,phi,psi,z;
    
        void InitTrain(const Dataset& dataset, const Dataset& neighborset, const Dataset& clustset, const Dataset& phiset, const Dataset& psiset, Vector& alpha, Matrix& beta, Vector& sigma, Vector& gamma, Matrix& Vgamma, Matrix& pi, Matrix& pi_y, int c,double max_iter_hmrf); 
        void HMRF(Vector& alpha, Matrix& beta, Vector& sigma, Vector& gamma, Matrix& Vgamma, Matrix& pi,Matrix& pi_y, double max_iter_mc);
        void RandomEffect(Vector& alpha, Matrix& beta, Matrix& pi_y, Vector& sigma, Vector& gamma, Matrix& Vgamma, double max_iter_mc);       
        void MStep(Vector& gamma, Matrix& Vgamma, Matrix& pi_y, Vector& alpha, Matrix& beta ,Vector& sigma, Matrix& pi,double max_iter_hmrf);    
        //void LogLikelihood(Vector& alpha, Matrix& beta, Vector& sigma, Matrix& pi, double log_likelihood_new);
        void Results(Vector& alpha, Matrix& beta, Vector& gamma, Matrix& pi_y);
        //void STCM::Predict(Matrix& beta, Vector& u, Vector& sigma_t, Vector& sigma_s, double& sigma, Matrix& y_pred);
};

int hasnan(Matrix x);


#endif
