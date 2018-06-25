#include "fscm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R_ext/Print.h>


void FSCM::InitTrain(const Dataset& dataset, const Dataset& neighborset, const Dataset& clustset, const Dataset& phiset, const Dataset& psiset, Vector& alpha, Matrix& beta, Vector& sigma, Vector& gamma, Matrix& Vgamma, Matrix& pi, Matrix& pi_y, int c,double max_iter_hmrf){

    //printf("1. read input datasets\n");
    y.Copy(dataset.matrix()); //m-by-n matrix
    neighbors.Copy(neighborset.matrix()); //n-by-m matrix
    z.Copy(clustset.matrix());
    phi.Copy(phiset.matrix()); //m-by-p//
    psi.Copy(psiset.matrix()); //n-by-q//

    fl__index_t n=psi.n_rows(), m=phi.n_rows(), p=phi.n_cols(), q=psi.n_cols(), k=neighbors.n_cols();
    
    //printf("2. Initialize alpha\n");
    Matrix y_tmp,phi_tran_phi,phi_tran,psi_tran_psi,psi_tran;
    y_tmp.Copy(y);
    la::MulTransAInit(phi,phi,&phi_tran_phi);
    la::Inverse(&phi_tran_phi);//(phi'phi)^{-1}
    la::TransposeInit(phi,&phi_tran); //phi'
    la::MulTransAInit(psi,psi,&psi_tran_psi); 
    la::Inverse(&psi_tran_psi); //(psi'psi)^{-1}
    la::TransposeInit(psi,&psi_tran); //psi'

    Vector phi_tran_y,psi_tran_y;
    phi_tran_y.Init(p);
    phi_tran_y.SetZero();

    psi_tran_y.Init(q);
    psi_tran_y.SetZero();
    
    for(fl__index_t s_j=0;s_j<n;s_j++){
        Vector y_j;
        y_tmp.MakeColumnVector(s_j,&y_j);
        la::MulExpert(1.0/n,phi_tran,y_j,1,&phi_tran_y);//sum_j 1/n*phi'y_j

    }

    la::MulInit(phi_tran_phi,phi_tran_y,&alpha); //(phi'phi)^(-1) sum_j 1/n*phi'y_j

    //printf("3. Initialize pi_y\n"); 
    pi_y.Init(n,c);
    pi_y.SetZero();
    for(fl__index_t s_j=0;s_j<n;s_j++){
    
        pi_y.set(s_j,z.get(s_j,0),1);
    
    }

    //printf("4. Initialize beta_k\n");
    beta.Init(p,c);
    beta.SetZero();

    Vector n_c; //n_k=sum_j z_jk
    n_c.Init(c);
    n_c.SetZero(); 

    double sum_nc_inv=0;

    for(fl__index_t c_i=0;c_i<c;c_i++){
        
        for(fl__index_t s_j=0;s_j<n;s_j++){
                
            n_c[c_i]+=pi_y.get(s_j,c_i);
            
        } //n_c[c_i] = sum_j z_jc 
        
        sum_nc_inv+=1.0/n_c[c_i]; //sum_k 1/n_k
    
    }
    
    //printf("2.1 Compute the Lagrange multiplier lambda\n");
	
    Vector lambda;
    lambda.Init(p);
    lambda.SetZero();

    for(fl__index_t s_j=0;s_j<n;s_j++){

        //remove global effects//
        Vector y_j;
        y_tmp.MakeColumnVector(s_j,&y_j); 
        la::MulExpert(-1,phi,alpha,1,&y_j);

        for(fl__index_t c_i=0;c_i<c;c_i++){

			//lambda=sum_j sum_k z_jk/n_k*phi'(y_j-phi*alpha)
            la::MulExpert(pi_y.get(s_j,c_i)/n_c[c_i],phi_tran,y_j,1,&lambda); //Caution! if n_c[c_i]=0, it will produce NA!

        }//end of c_i

    }//end of s_j

    la::Scale(1/sum_nc_inv,&lambda);

    //printf("2.2 Compute beta_k\n");
    for(fl__index_t c_i=0;c_i<c;c_i++){

        Vector beta_k;
        beta.MakeColumnVector(c_i,&beta_k);
        phi_tran_y.SetZero(); //clear up for each beta_k

        for(fl__index_t s_j=0;s_j<n;s_j++){

            Vector y_j;
            y_tmp.MakeColumnVector(s_j,&y_j);
            la::MulExpert(pi_y.get(s_j,c_i),phi_tran,y_j,1,&phi_tran_y); //sum_j z_jk*phi'*(y_j-phi*alpha)
            
        }//end of s_j

            la::SubFrom(lambda,&phi_tran_y); //sum_j z_jk*phi'*(y_j-phi*alpha)-lambda
            la::MulOverwrite(phi_tran_phi,phi_tran_y,&beta_k); //(phi'phi)^{-1}*[sum_j z_jk*phi'*(y_j-phi*alpha)-lambda]
            la::Scale(1.0/n_c[c_i],&beta_k); //1/n_k*(phi'phi)^{-1}*[sum_j z_jk*phi'*(y_j-phi*alpha)-lambda]

    }//end of c_i
	
    //printf("5. Initialize gamma\n");

    gamma.Init(q);

    //remove cluster effects//
    for(fl__index_t s_j=0;s_j<n;s_j++){

        Vector y_j;
        y_tmp.MakeColumnVector(s_j,&y_j);

        for(fl__index_t c_i=0;c_i<c;c_i++){
            
            Vector beta_k;
            beta.MakeColumnVector(c_i,&beta_k);   
            la::MulExpert(-pi_y.get(s_j,c_i),phi,beta_k,1,&y_j); //y_j-phi*alpha-sum_k z_jk*phi*beta_k
            
        }
            		
    }
   
    Matrix y_tran;
    la::TransposeInit(y_tmp,&y_tran);

    for(fl__index_t t_i=0;t_i<m;t_i++){
    
        Vector y_i;
        y_tran.MakeColumnVector(t_i,&y_i);
        la::MulExpert(1.0/m, psi_tran,y_i,1,&psi_tran_y); // 1/m*sum_i psi'(y_i-phi_i*alpha-E*phi_i*beta_k)

    }

    la::MulInit(psi_tran_psi,psi_tran_y,&gamma); // 1/m*sum_i (psi'psi)^{-1}*psi'(y_i-phi_i*alpha-E*phi_i*beta_k)
	
    //printf("6. Initial sigma_s and sigma_epsilon\n");

    sigma.Init(2); //sigma_s + sigma_epsilon
    sigma.SetZero();

    Vector tau,one;
    la::MulInit(psi,gamma,&tau);
    one.Init(m);
    one.SetAll(1);
    
    for(fl__index_t s_j=0;s_j<n;s_j++){
    
        Vector y_j;
        y_tmp.MakeColumnVector(s_j,&y_j);
        la::AddExpert(-tau[s_j],one,&y_j); //y_j-phi*alpha-sum_k z_j[k]*phi*beta_k - psi_j*gamma
        sigma[0]+=la::Dot(y_j,y_j)/m/n;
       
    }
    
    if(verbose) Rprintf("Initial sigma_epsilon is %f\n", sigma[0]);

    one.Destruct();
    one.Init(q);
    one.SetAll(1);

    double gamma_mean=0;

    for(fl__index_t q_i=0;q_i<q;q_i++){
    
        gamma_mean+=gamma[q_i]/q; 
    
    }

    la::AddExpert(-gamma_mean,one,&gamma);

    sigma[1]=la::Dot(gamma,gamma)/q;	
	if(verbose) Rprintf("Initial sigma_s is %f\n", sigma[1]);

	la::AddExpert(gamma_mean,one,&gamma);
	
	la::MulTransAInit(psi,psi,&Vgamma);
    //printf("%f\n",Vgamma.get(0,0));
    la::Scale(m,&Vgamma); //Vgamma=mPsi'Psi
    //printf("%f\n",Vgamma.get(0,0));

    for(fl__index_t q_i=0;q_i<q;q_i++){

        Vgamma.set(q_i,q_i,Vgamma.get(q_i,q_i)+sigma[0]/sigma[1]); //Vgamma=mPsi'Psi+sigma_e/sigma_s I

    }

    la::Inverse(&Vgamma); //Vgamma=(m psi'psi+sigma/tau)^{-1} =1/m(psi'psi+sigma_e/sigma_s/m I)^{-1}
	
	la::Scale(sigma[0],&Vgamma);	
	
	//data::Save("Vgamma.csv",Vgamma);
    
    //printf("7 Initialize pi \n");
    //printf("7.1 Initialize the Gibbs parameter theta: find theta which minimize the neg_log_f using Bisection Method\n");

    pi.Init(n,c);
    
    //log_f=sum_j log_f_j//
    double theta_ub=1,theta_lb=0, neg_log_f_ub=0,neg_log_f_lb=0; //initial bounds

    //printf("check if the inital bounds bracket the root\n");
    
    for (fl__index_t s_j=0;s_j<n; s_j++){
    
        double der_Vj_lb=0,der_Vj_ub=0,Vj_lb=0,Vj_ub=0; //vary with s_j
        Vector n_j;
        n_j.Init(c);
        n_j.SetZero();

        for(fl__index_t c_i=0;c_i<c; c_i++){
        
            for(fl__index_t k_i=0;k_i<k;k_i++){

                fl__index_t s_i=neighbors.get(s_j,k_i);
                n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
            } //end of k_i
        
            Vj_lb +=exp(theta_lb*n_j[c_i]); // sum_c U_jc(theta_lb)=sum_c exp(theta_lb n_jc)
            der_Vj_lb +=n_j[c_i]*exp(theta_lb*n_j[c_i]); //sum_c derivative of U_jc(theta_lb)= sum_c n_jc exp(theta_lb n_jc)
        
            Vj_ub +=exp(theta_ub*n_j[c_i]); // sum_c U_jc(theta_ub)=sum_c exp(theta_ub n_jc)
            der_Vj_ub +=n_j[c_i]*exp(theta_ub*n_j[c_i]); //sum_c derivative of U_jc(theta_ub)= sum_c n_jc exp(theta_ub n_jc)
            
        }//end of c_i

        //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
        for(fl__index_t c_i=0;c_i<c;c_i++){
            
            neg_log_f_lb += pi_y.get(s_j,c_i)*(der_Vj_lb/Vj_lb-n_j[c_i]); //derivative of -log p(z_jc=1)
            neg_log_f_ub += pi_y.get(s_j,c_i)*(der_Vj_ub/Vj_ub-n_j[c_i]);
    
        }//end of c_i

    }//end of s_j

    while (neg_log_f_lb > 0){
    
        neg_log_f_lb=0;  
        theta_lb--;
  
        for (fl__index_t s_j=0;s_j<n; s_j++){
    
            double der_Vj_lb=0,Vj_lb=0; 
            Vector n_j;
            n_j.Init(c);
            n_j.SetZero();

            for(fl__index_t c_i=0;c_i<c; c_i++){
        
                for(fl__index_t k_i=0;k_i<k;k_i++){

                    fl__index_t s_i=neighbors.get(s_j,k_i);
                    n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                    //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
                } //end of k_i
        
                Vj_lb +=exp(theta_lb*n_j[c_i]); // sum_c U_jc(theta_lb)=sum_c exp(theta_lb n_jc)
                der_Vj_lb +=n_j[c_i]*exp(theta_lb*n_j[c_i]); //sum_c derivative of U_jc(theta_lb)= sum_c n_jc exp(theta_lb n_jc)
                    
            }//end of c_i

            //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
            for(fl__index_t c_i=0;c_i<c;c_i++){
        
                neg_log_f_lb += pi_y.get(s_j,c_i)*(der_Vj_lb/Vj_lb-n_j[c_i]); //derivative of -log p(z_jc=1)
            
            }//end of c_i

        }//end of s_j
    
    }//end of while

    while (neg_log_f_ub < 0){
    
        neg_log_f_ub=0;  
        theta_ub++;

        for (fl__index_t s_j=0;s_j<n; s_j++){
    
            double der_Vj_ub=0,Vj_ub=0; 
            Vector n_j;
            n_j.Init(c);
            n_j.SetZero();

            for(fl__index_t c_i=0;c_i<c; c_i++){
        
                for(fl__index_t k_i=0;k_i<k;k_i++){

                    fl__index_t s_i=neighbors.get(s_j,k_i);
                    n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                    //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
                } //end of k_i
        
                Vj_ub +=exp(theta_ub*n_j[c_i]); // sum_c U_jc(theta_ub)=sum_c exp(theta_ub n_jc)
                der_Vj_ub +=n_j[c_i]*exp(theta_ub*n_j[c_i]); //sum_c derivative of U_jc(theta_ub)= sum_c n_jc exp(theta_ub n_jc)
            
            }//end of c_i

            //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
            for(fl__index_t c_i=0;c_i<c;c_i++){
        
                neg_log_f_ub += pi_y.get(s_j,c_i)*(der_Vj_ub/Vj_ub-n_j[c_i]);
    
            }//end of c_i

        }//end of s_j

    }//end of while

    //printf("find the theta within the bounds\n");

    int iter=0;

    double theta=0.5*(theta_ub+theta_lb);

    while(iter<max_iter_hmrf){
    
        double neg_log_f=0;

        for(fl__index_t s_j=0;s_j<n; s_j++){
    
            double der_Vj=0,Vj=0; //
            Vector n_j;
            n_j.Init(c);
            n_j.SetZero();

            for(fl__index_t c_i=0;c_i<c; c_i++){
        
                for(fl__index_t k_i=0;k_i<k; k_i++){

                    fl__index_t s_i=neighbors.get(s_j,k_i);
                    n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                    //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
                } //end of k_i
        
                Vj+=exp(theta*n_j[c_i]); // sum_c U_jc(theta_ub)=sum_c exp(theta_ub n_jc)
                der_Vj+=n_j[c_i]*exp(theta*n_j[c_i]); //sum_c derivative of U_jc(theta_ub)= sum_c n_jc exp(theta_ub n_jc)
            
            }//end of c_i

            //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
            for(fl__index_t c_i=0;c_i<c;c_i++){
        
                neg_log_f += pi_y.get(s_j,c_i)*(der_Vj/Vj-n_j[c_i]); //derivative of -log p(z_jc=1)
            
            }//end of c_i

        }//end of s_j

        if(neg_log_f>0) theta_ub=theta;
        if(neg_log_f<0) theta_lb=theta;

            theta=0.5*(theta_ub+theta_lb);
            iter++;

        //printf("%dth iteration: neg_log_f=%f, phi=%f, phi_ub=%f, phi_lb=%f\n",j,neg_log_f,phi,phi_ub,phi_lb);
    }

    if(verbose) Rprintf("Final theta is %f\n",theta);
     FILE *fp1 = fopen("theta.csv", "w");
    fprintf(fp1, "%f\n", theta);
    fclose(fp1);

    //printf("7.2 Compute pi \n");
    for(fl__index_t s_j=0;s_j<n;s_j++){

        Vector n_j;
        n_j.Init(c);
        n_j.SetZero(); //sum_{i in neighbor_j} z_ic
        double V_j=0; // the normalizing parameter V_j

        //printf("1.1 Compute the normalizing paramter of s_j: sum_k theta*sum_{i in neighbor j} z_ic \n");
        for(fl__index_t c_i=0;c_i<c;c_i++){
                
            for(fl__index_t k_i=0;k_i<k;k_i++){
                    
                fl__index_t s_i=neighbors.get(s_j,k_i); //neighbors of s_j
                n_j[c_i]+=pi_y.get(s_i,c_i); //sum over k
                
            }//end of k_i
        
            V_j+=exp(theta*n_j[c_i]); //sum over c
            
        }//end of c_i

        //printf("1.2. Compute Prob(z_jc)\n");      
        for(fl__index_t c_i=0;c_i<c;c_i++){
            
            pi.set(s_j,c_i,exp(theta*n_j[c_i])/V_j);
            
        }

    }
	
	//8. Save initial estimates
	Matrix spatial;
    Matrix temporal;

    spatial.Init(n,2);
    temporal.Init(m,c+1);

    Vector temp;
    spatial.MakeColumnVector(0,&temp);
    la::MulOverwrite(psi,gamma, &temp);

    for(fl__index_t s_j=0;s_j<n;s_j++){
        
        fl__index_t z_j=0;
    
        for(fl__index_t c_i=1;c_i<c;c_i++){
        
            if(pi_y.get(s_j,c_i)>pi_y.get(s_j,z_j)){
            
                z_j=c_i;
            
            }//end of if        
            
        }//end of c_i
    
        spatial.set(s_j,1,z_j);

    }//end of s_j

    temp.Destruct();
    temporal.MakeColumnVector(0,&temp);
    la::MulOverwrite(phi,alpha, &temp);

    for(fl__index_t c_i=0;c_i<c;c_i++){
    
        Vector beta_k,temporal_k;
        beta.MakeColumnVector(c_i,&beta_k);
        temporal.MakeColumnVector(c_i+1,&temporal_k);
        la::MulOverwrite(phi,beta_k,&temporal_k);

    }

    //data::Save("spatial_init.csv",spatial);
    //data::Save("temporal_init.csv",temporal);
	
}

void FSCM::HMRF(Vector& alpha, Matrix& beta, Vector& sigma, Vector& gamma, Matrix& Vgamma, Matrix& pi, Matrix& pi_y, double max_iter_mc){

    fl__index_t n=psi.n_rows(), m=phi.n_rows(), c=beta.n_cols(), q=psi.n_cols();

    //#pragma GCC diagnostic ignored "-Wunused-variable"
    //fl__index_t k = neighbors.n_cols();
    //fl__index_t p = phi.n_cols();
    
    //ArrayList<fl__index_t> random_index;
    //math::MakeRandomPermutation(n,&random_index);

    pi_y.SetZero(); //f(z_ik=1|y_1,...,y_n)

    Matrix y_tmp,psi_tran;
    la::TransposeInit(psi,&psi_tran);
    y_tmp.Copy(y);

    Vector one_t;
    one_t.Init(m);
    one_t.SetAll(1);
    
    //printf("Remove global effects\n");
    for(fl__index_t s_j=0;s_j<n;s_j++){
        
        Vector y_j;
        y_tmp.MakeColumnVector(s_j,&y_j);
        la::MulExpert(-1,phi,alpha,1,&y_j); //y_j-phi*alpha
    
    }
        
	Matrix H;
	la::CholeskyInit(Vgamma,&H);
	la::TransposeSquare(&H);
	
    for(fl__index_t r=0;r<max_iter_mc;r++){

        //2.1 Generate gamma_r
        Vector ind_norm,gamma_r;
        ind_norm.Init(q);
    
        //Generate independent standard normal distributed variables
        //printf("generate independent standard normal distributed variables\n");
        for(fl__index_t q_i=0;q_i<q;q_i++){

            double u=math::Random(); 
            double v=math::Random();

            ind_norm[q_i]=sqrt(-2*log(u))*cos(2*math::PI*v);
        
        }

        la::MulInit(H,ind_norm,&gamma_r);
		la::AddTo(gamma,&gamma_r);
		
		Matrix pi_y_r;
		pi_y_r.Init(n,c);
		pi_y_r.SetZero();
              
        for(fl__index_t s_j=0;s_j<n;s_j++){
        
            //printf("Remove spatial effect from rth gamma sample\n");
            Vector y_j,psi_j;
            y_tmp.MakeColumnVector(s_j,&y_j);
            psi_tran.MakeColumnVector(s_j,&psi_j);
            
            double tau_s=la::Dot(psi_j,gamma_r); 
            la::AddExpert(-tau_s,one_t,&y_j); //y_j-phi*alpha-psi_j*gamma

            for(fl__index_t c_i=0;c_i<c;c_i++){

                //remove cluster effects
                Vector beta_k;
                beta.MakeColumnVector(c_i,&beta_k);
                la::MulExpert(-1,phi,beta_k,1,&y_j); //y_j-phi*alpha-phi*beta_k-psi_j*gamma
              
                //printf("2.5 compute Prob(y_j|gamma_r)\n");
                double prob_y_c=exp(-0.5/sigma[0]*la::Dot(y_j,y_j)); //f(y_j|z_jk=1,gamma_r)       
                pi_y_r.set(s_j,c_i,pi.get(s_j,c_i)*prob_y_c); //f(y_j,z_jk|gamma_r)=pi_jk*f(y_j|gamma_r,z_jk)       

                //add back cluster effect
                la::MulExpert(1,phi,beta_k,1,&y_j); //y_j-phi*alpha-psi_j*gamma
            
            }//end of c_i

            //add back the spatial effects from rth gamma sample
            la::AddExpert(tau_s,one_t,&y_j);   //y_j-phi*alpha      

        }//end of s_j

        //2.3 Compute f(y_1,...,y_n, z_jk=1)

        for(fl__index_t s_j=0;s_j<n;s_j++){

			double prob_y_r=0;
			
            for(fl__index_t c_i=0;c_i<c;c_i++){
            
                prob_y_r+=pi_y_r.get(s_j,c_i);                       
                
            }//end of c_i
			
			for(fl__index_t c_i=0;c_i<c;c_i++){
			 
				pi_y.set(s_j,c_i,pi_y.get(s_j,c_i)+pi_y_r.get(s_j,c_i)/prob_y_r);
			
			}//end of c_i

        }//end of s_j

    }//end of r    

    //printf("2.4 Compute f(y_1,...,y_n)\n"); 
   la::Scale(1/max_iter_mc,&pi_y);
   
   //apply hard clustering
   for(fl__index_t s_j=0;s_j<n;s_j++){
	
		double z_j=0;
		
		for(fl__index_t c_i=(z_j+1);c_i<c;c_i++){
		
			if(pi_y.get(s_j,c_i)>pi_y.get(s_j,z_j)){
			
				z_j=c_i;
			
			}
		
		}
	
		z.set(s_j,0,z_j);
	
	}
	
	pi_y.SetZero();
    for(fl__index_t s_j=0;s_j<n;s_j++){
    
        pi_y.set(s_j,z.get(s_j,0),1);
    
    }
       
} //end of function

void FSCM::RandomEffect(Vector& alpha, Matrix& beta, Matrix& pi_y, Vector& sigma, Vector& gamma, Matrix& Vgamma, double max_iter_mc){

    fl__index_t n=psi.n_rows(), m=phi.n_rows(), /*c=beta.n_cols(), p=phi.n_cols(),*/ q=psi.n_cols();


    //printf("Compute Vgamma=Cov(gamma|Y)\n");
	la::MulTransAOverwrite(psi,psi,&Vgamma);
    la::Scale(m,&Vgamma); //Vgamma=mPsi'Psi

    for(fl__index_t q_i=0;q_i<q;q_i++){

        Vgamma.set(q_i,q_i,Vgamma.get(q_i,q_i)+sigma[0]/sigma[1]); //Vgamma=mPsi'Psi+sigma_e/sigma_s I

    }



    la::Inverse(&Vgamma); //Vgamma=(m psi'psi+sigma/tau)^{-1} =1/m(psi'psi+sigma_e/sigma_s/m I)^{-1}
    
    Matrix y_tmp;
    y_tmp.Copy(y);
	
		
    for(fl__index_t s_j=0;s_j<n;s_j++){
	
		//printf("Remove global effects\n");
		Vector y_j, beta_k;
        y_tmp.MakeColumnVector(s_j,&y_j);         
        la::MulExpert(-1,phi,alpha,1,&y_j); //y_j-phi*alpha
	
        beta.MakeColumnVector(z.get(s_j,0),&beta_k);        
        la::MulExpert(-1,phi,beta_k,1,&y_j); //y_j-phi*alpha-phi*beta_k

    }//end of s_j
        
    Matrix y_tran, psi_tran;
    la::TransposeInit(y_tmp,&y_tran); // n-by-m matrix
    la::TransposeInit(psi,&psi_tran);

    Vector psi_tran_y;
    psi_tran_y.Init(q);
    psi_tran_y.SetZero();

    for(fl__index_t t_i=0;t_i<m;t_i++){
    
        Vector y_i;
        y_tran.MakeColumnVector(t_i,&y_i);
        la::MulExpert(1, psi_tran,y_i,1,&psi_tran_y); //sum_i psi'(y_i - phi_i alpha - phi_i beta_k)
    
    }

    la::MulOverwrite(Vgamma,psi_tran_y,&gamma);  //E[gamma|z]=(Psi'Psi+sigma/tau/m I)^{-1}/m*sum_i psi'(y_i - phi_i alpha - phi_i beta_k)

    la::Scale(sigma[0],&Vgamma);

}//end of function

void FSCM::MStep(Vector& gamma, Matrix& Vgamma, Matrix& pi_y, Vector& alpha, Matrix& beta ,Vector& sigma, Matrix& pi,double max_iter_hmrf){

    fl__index_t n=psi.n_rows(), m=phi.n_rows(), c=beta.n_cols(), p=phi.n_cols(), q=psi.n_cols(), k=neighbors.n_cols();
    
    Matrix phi_tran_phi,phi_tran,psi_tran;
    la::MulTransAInit(phi,phi,&phi_tran_phi); //Phi'Phi
    la::Inverse(&phi_tran_phi); //(Phi'Phi)^{-1}
    la::TransposeInit(phi,&phi_tran);
    la::TransposeInit(psi,&psi_tran);
    
    Vector phi_tran_y;
    //printf("1. Update alpha\n");
    phi_tran_y.Init(p);
    phi_tran_y.SetZero();

    Matrix y_tmp;
    y_tmp.Copy(y);

    for(fl__index_t s_j=0;s_j<n;s_j++){
        
        Vector y_j,psi_j,effect;
        y_tmp.MakeColumnVector(s_j,&y_j);
		psi_tran.MakeColumnVector(s_j,&psi_j);
		
		//remove spatial effects 
		double tau=la::Dot(psi_j,gamma);            
        effect.Init(m); 
		effect.SetAll(tau);
        la::AddExpert(-1,effect,&y_j); //y_j-psi_j*gamma
        
		//printf("remove cluster effects\n");
        Vector beta_k;
        beta.MakeColumnVector(z.get(s_j,0),&beta_k);
        la::MulExpert(-1,phi,beta_k,1,&y_j); //y_j-phi*beta_k-psi_j*gamma
            
        //printf("compute sum_j sum_k z_jk*phi'(y_j-phi*beta_k-psi_j*gamma)\n");
        la::MulExpert(1,phi_tran,y_j,1,&phi_tran_y); 
		
		//printf("add back cluster effects\n");
		la::MulExpert(1,phi,beta_k,1,&y_j); 
       
    }//end of s_j

    la::MulOverwrite(phi_tran_phi,phi_tran_y,&alpha); //  (phi'phi)^{-1}[sum_j sum_k z_jk*phi'(y_j-phi*beta_k-psi_j*gamma)]

    la::Scale(1.0/n,&alpha);
    //   alpha.SetZero();
	
    //printf("2. Update beta\n");
    Vector n_c; //n_k=sum_j z_jk
    n_c.Init(c);
    n_c.SetZero(); 

    double sum_nc_inv=0;

    for(fl__index_t c_i=0;c_i<c;c_i++){
        
        for(fl__index_t s_j=0;s_j<n;s_j++){
                
            n_c[c_i]+=pi_y.get(s_j,c_i);
            
        } //n_c[c_i] = sum_j z_jc 
        
        sum_nc_inv+=1/n_c[c_i]; //sum_k 1/n_k
    
    }
    
    //printf("2.1 Compute the Lagrange multiplier lambda\n");
    Vector lambda;
    lambda.Init(p);
    lambda.SetZero();

    for(fl__index_t s_j=0;s_j<n;s_j++){

        //remove global effects//
        Vector y_j;
        y_tmp.MakeColumnVector(s_j,&y_j);
        la::MulExpert(-1,phi,alpha,1,&y_j); //y_j-phi*alpha-psi_j*gamma

        for(fl__index_t c_i=0;c_i<c;c_i++){
		
			//sum_j sum_k z_jk/n_k*phi'(y_j-phi*alpha-psi_j*gamma)
            la::MulExpert(pi_y.get(s_j,c_i)/n_c[c_i],phi_tran,y_j,1,&lambda); //Caution! if n_c[c_i]=0, it will produce NA!

        }//end of c_i

    }//end of s_j

    la::Scale(1/sum_nc_inv,&lambda);

    //printf("2.2 Compute beta_k\n");
    for(fl__index_t c_i=0;c_i<c;c_i++){

        Vector beta_k;
        beta.MakeColumnVector(c_i,&beta_k);
        phi_tran_y.SetZero(); //clear up for each beta_k

        for(fl__index_t s_j=0;s_j<n;s_j++){

            Vector y_j;
            y_tmp.MakeColumnVector(s_j,&y_j);
            la::MulExpert(pi_y.get(s_j,c_i),phi_tran,y_j,1,&phi_tran_y); //sum_j z_jk/n_k*phi'(y_j-phi*alpha-psi_j*gamma_jk)
           
        }//end of s_j

            la::SubFrom(lambda,&phi_tran_y); ////sum_j z_jk/n_k*phi'(y_j-phi*alpha-psi_j*gamma_jk)-lambda
            la::MulOverwrite(phi_tran_phi,phi_tran_y,&beta_k);
            la::Scale(1/n_c[c_i],&beta_k); 

    }//end of c_i

    //printf("3. Update sigma_epsilon\n");
    sigma[0]=0;

    for(fl__index_t s_j=0; s_j<n; s_j++){

        Vector y_j,beta_k;
        y_tmp.MakeColumnVector(s_j,&y_j);
        beta.MakeColumnVector(z.get(s_j,0),&beta_k);
        la::MulExpert(-1,phi,beta_k,1,&y_j); //y_j-phi*alpha-phi*beta_k-psi_j*gamma_jk
            
        sigma[0]+=la::Dot(y_j,y_j)/n/m;

    }

	if(verbose) Rprintf("sigma is %f \n",sigma[0]);
	
    Matrix temp,psi_Vgamma_psi;
    la::MulInit(psi,Vgamma,&temp);
    la::MulInit(temp,psi_tran,&psi_Vgamma_psi);

    for(fl__index_t s_j=0;s_j<n;s_j++){
    
        sigma[0]+=psi_Vgamma_psi.get(s_j,s_j)/n;

    }

    if(verbose) Rprintf("sigma is %f \n",sigma[0]);

    //printf("4. Update sigma_s");
    //compute E[gamma|y]
	
    sigma[1]=la::Dot(gamma,gamma)/q;
	
	if(verbose) Rprintf("sigma_s is %f \n",sigma[1]);

    for(fl__index_t q_i=0;q_i<q;q_i++){
    
        sigma[1]+=Vgamma.get(q_i,q_i)/q;

    }

	if(verbose) Rprintf("sigma_s is %f \n",sigma[1]);
	
    //printf("5. Update the Gibbs parameter theta: find theta which minimize the neg_log_f using Bisection Method\n");
    
    //log_f=sum_j log_f_j//
    double theta_ub=1,theta_lb=0, neg_log_f_ub=0,neg_log_f_lb=0; //initial bounds

    //printf("5.1 Check if the inital bounds bracket the root\n");
    
    for (fl__index_t s_j=0;s_j<n; s_j++){
    
        double der_Vj_lb=0,der_Vj_ub=0,Vj_lb=0,Vj_ub=0; //vary with s_j
        Vector n_j;
        n_j.Init(c);
        n_j.SetZero();

        for(fl__index_t c_i=0;c_i<c; c_i++){
        
            for(fl__index_t k_i=0;k_i<k;k_i++){

                fl__index_t s_i=neighbors.get(s_j,k_i);
                n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
            } //end of k_i
        
            Vj_lb +=exp(theta_lb*n_j[c_i]); // sum_c U_jc(theta_lb)=sum_c exp(theta_lb n_jc)
            der_Vj_lb +=n_j[c_i]*exp(theta_lb*n_j[c_i]); //sum_c derivative of U_jc(theta_lb)= sum_c n_jc exp(theta_lb n_jc)
        
            Vj_ub +=exp(theta_ub*n_j[c_i]); // sum_c U_jc(theta_ub)=sum_c exp(theta_ub n_jc)
            der_Vj_ub +=n_j[c_i]*exp(theta_ub*n_j[c_i]); //sum_c derivative of U_jc(theta_ub)= sum_c n_jc exp(theta_ub n_jc)
            
        }//end of c_i

        //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
        for(fl__index_t c_i=0;c_i<c;c_i++){
            
            neg_log_f_lb += pi_y.get(s_j,c_i)*(der_Vj_lb/Vj_lb-n_j[c_i]); //derivative of -log p(z_jc=1)
            neg_log_f_ub += pi_y.get(s_j,c_i)*(der_Vj_ub/Vj_ub-n_j[c_i]);
    
        }//end of c_i

    }//end of s_j

    while (neg_log_f_lb > 0){
    
        neg_log_f_lb=0;  
        theta_lb--;
  
        for (fl__index_t s_j=0;s_j<n; s_j++){
    
            double der_Vj_lb=0,Vj_lb=0; 
            Vector n_j;
            n_j.Init(c);
            n_j.SetZero();

            for(fl__index_t c_i=0;c_i<c; c_i++){
        
                for(fl__index_t k_i=0;k_i<k;k_i++){

                    fl__index_t s_i=neighbors.get(s_j,k_i);
                    n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                    //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
                } //end of k_i
        
                Vj_lb +=exp(theta_lb*n_j[c_i]); // sum_c U_jc(theta_lb)=sum_c exp(theta_lb n_jc)
                der_Vj_lb +=n_j[c_i]*exp(theta_lb*n_j[c_i]); //sum_c derivative of U_jc(theta_lb)= sum_c n_jc exp(theta_lb n_jc)
                    
            }//end of c_i

            //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
            for(fl__index_t c_i=0;c_i<c;c_i++){
        
                neg_log_f_lb += pi_y.get(s_j,c_i)*(der_Vj_lb/Vj_lb-n_j[c_i]); //derivative of -log p(z_jc=1)
            
            }//end of c_i

        }//end of s_j
    
    }//end of while

    while (neg_log_f_ub < 0){
    
        neg_log_f_ub=0;  
        theta_ub++;

        for (fl__index_t s_j=0;s_j<n; s_j++){
    
            double der_Vj_ub=0,Vj_ub=0; 
            Vector n_j;
            n_j.Init(c);
            n_j.SetZero();

            for(fl__index_t c_i=0;c_i<c; c_i++){
        
                for(fl__index_t k_i=0;k_i<k;k_i++){

                    fl__index_t s_i=neighbors.get(s_j,k_i);
                    n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                    //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
                } //end of k_i
        
                Vj_ub +=exp(theta_ub*n_j[c_i]); // sum_c U_jc(theta_ub)=sum_c exp(theta_ub n_jc)
                der_Vj_ub +=n_j[c_i]*exp(theta_ub*n_j[c_i]); //sum_c derivative of U_jc(theta_ub)= sum_c n_jc exp(theta_ub n_jc)
            
            }//end of c_i

            //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
            for(fl__index_t c_i=0;c_i<c;c_i++){
        
                neg_log_f_ub += pi_y.get(s_j,c_i)*(der_Vj_ub/Vj_ub-n_j[c_i]);
    
            }//end of c_i

        }//end of s_j

    }//end of while

    //printf("5.2 Find the theta within the bounds\n");

    int iter=0;

    double theta=0.5*(theta_ub+theta_lb);

    while(iter<max_iter_hmrf){
    
        double neg_log_f=0;

        for(fl__index_t s_j=0;s_j<n; s_j++){
    
            double der_Vj=0,Vj=0; //
            Vector n_j;
            n_j.Init(c);
            n_j.SetZero();

            for(fl__index_t c_i=0;c_i<c; c_i++){
        
                for(fl__index_t k_i=0;k_i<k; k_i++){

                    fl__index_t s_i=neighbors.get(s_j,k_i);
                    n_j[c_i]+=pi_y.get(s_i,c_i);    //sum_{i in neighbor_j} pi_ic
                    //if (z[neighbors.get(s_i,k_i)]==c_i) n_eq[c_i]++; 
                } //end of k_i
        
                Vj+=exp(theta*n_j[c_i]); // sum_c U_jc(theta_ub)=sum_c exp(theta_ub n_jc)
                der_Vj+=n_j[c_i]*exp(theta*n_j[c_i]); //sum_c derivative of U_jc(theta_ub)= sum_c n_jc exp(theta_ub n_jc)
            
            }//end of c_i

            //double -log Prob(z_1,...,z_n)=- sum_j log Prob(z_j1,...,z_jc) = - sum_j log (pi_j1^z_j1...pi_jc^z_jc) = - sum_j sum_c E[z_jc|Y_j] log pi_jc
            for(fl__index_t c_i=0;c_i<c;c_i++){
        
                neg_log_f += pi_y.get(s_j,c_i)*(der_Vj/Vj-n_j[c_i]); //derivative of -log p(z_jc=1)
            
            }//end of c_i

        }//end of s_j

        if(neg_log_f>0) theta_ub=theta;
        if(neg_log_f<0) theta_lb=theta;

            theta=0.5*(theta_ub+theta_lb);
            iter++;

        //printf("%dth iteration: neg_log_f=%f, phi=%f, phi_ub=%f, phi_lb=%f\n",j,neg_log_f,phi,phi_ub,phi_lb);
    }

    if(verbose) Rprintf("5.3 Final theta is %f\n",theta);

    //printf("6. Compute p(z_j)\n");
    for(fl__index_t s_j=0;s_j<n;s_j++){

        Vector n_j;
        n_j.Init(c);
        n_j.SetZero(); //sum_{i in neighbor_j} z_ic
        double V_j=0; // the normalizing parameter V_j

        //printf("1.1 Compute the normalizing paramter of s_j: sum_k theta*sum_{i in neighbor j} z_ic \n");
        for(fl__index_t c_i=0;c_i<c;c_i++){
                
            for(fl__index_t k_i=0;k_i<k;k_i++){
                    
                fl__index_t s_i=neighbors.get(s_j,k_i); //neighbors of s_j
                n_j[c_i]+=pi_y.get(s_i,c_i); //sum over k
                
            }//end of k_i
        
            V_j+=exp(theta*n_j[c_i]); //sum over c
            
        }//end of c_i

        //printf("1.2. Compute Prob(z_jc)\n");      
        for(fl__index_t c_i=0;c_i<c;c_i++){
            
            pi.set(s_j,c_i,exp(theta*n_j[c_i])/V_j);
            
        
}
    }

    
}

void FSCM::Results(Vector& alpha, Matrix& beta, Vector& gamma, Matrix& pi_y){

  fl__index_t n=psi.n_rows(),m=phi.n_rows(),c=beta.n_cols();

    Matrix spatial;
    Matrix temporal;
    Matrix probs;

    spatial.Init(n,2);
    temporal.Init(m,c+1);

    probs.Copy(pi_y);
    //probs.Init(n,c);
    data::Save("probsy.csv",probs);

    Vector temp;
    spatial.MakeColumnVector(0,&temp);
    la::MulOverwrite(psi,gamma, &temp);

    for(fl__index_t s_j=0;s_j<n;s_j++){
        
        fl__index_t z_j=0;
    
        for(fl__index_t c_i=1;c_i<c;c_i++){
        
            if(pi_y.get(s_j,c_i)>pi_y.get(s_j,z_j)){
            
                z_j=c_i;
            
            }//end of if        
            
        }//end of c_i
    
        spatial.set(s_j,1,z_j);

    }//end of s_j

    temp.Destruct();
    temporal.MakeColumnVector(0,&temp);
    la::MulOverwrite(phi,alpha, &temp);
    
    for(fl__index_t c_i=0;c_i<c;c_i++){
    
        Vector beta_k,temporal_k;
        beta.MakeColumnVector(c_i,&beta_k);
        temporal.MakeColumnVector(c_i+1,&temporal_k);
        la::MulOverwrite(phi,beta_k,&temporal_k);
    }
   
    data::Save("spatial.csv",spatial);
    data::Save("temporal.csv",temporal);
    data::Save("beta.csv",beta);
    
    GenVector<fl__index_t> index;
    index.Init(gamma.length());
    for(fl__index_t i=0; i<gamma.length();i++)
      index[i]=i;
    data::Save("gamma.csv", index, gamma);
    
    index.Init(alpha.length());
    for(fl__index_t i=0; i<alpha.length();i++)
      index[i]=i;
    data::Save("alpha.csv", index, alpha);
   
}


int hasnan(Matrix x) {
    int i;
    double *p = x.ptr();
    for(i=0; i<x.n_rows()*x.n_cols(); i++) {
        if(yet_another_isnan(p[i])) return 1;
    }
    return 0;
}


