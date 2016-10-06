//
//  hos_model.c
//  UCD_WAVE_HOS
//
//  Created by Joseph on 02/09/2016.
//  Copyright © 2016 Joseph Brennan. All rights reserved.
//

#include "hos_model.h"


#ifdef USE_DOUBLES
static fftw_complex*   hwn;
static fftw_complex*   hwM;
static fftw_complex*   hwM2;
static fftw_complex*   hw2M;
static fftw_complex*   hw2M2;
#else
static fftwf_complex*   hwn;
static fftwf_complex*   hwM;
static fftwf_complex*   hwM2;
static fftwf_complex*   hw2M;
static fftwf_complex*   hw2M2;
#endif
static float*         Coeff;

#ifdef USE_DOUBLES
void rhs_test(fftw_complex* hrhs, fftw_complex* hu){
#else
    void rhs_test(fftwf_complex* hrhs, fftwf_complex* hu){
#endif
        
        // Linear advection
        Dx(hu, hrhs);
        Dy(&hu[alloc_local], &hrhs[alloc_local]);
        
    }
    
    void rhs_hos_setup(){
        
        int n;
        
#ifdef USE_DOUBLES
        
        temp1 = fftw_alloc_real(2 * alloc_local);
        temp2 = fftw_alloc_real(2 * alloc_local);
        htemp1 = fftw_alloc_complex(2 * alloc_local);
        htemp2 = &htemp1[alloc_local];
        
        hetan = fftw_alloc_complex((NLevs+1)*alloc_local); //Nlevs (expansion coefficients) of heta
        hphin = fftw_alloc_complex((NLevs+1)*alloc_local); //Nlevs pf hphi
        hwn   = fftw_alloc_complex((NLevs+1)*alloc_local); //Nlevs of hwn
        
        
        hwM   = fftw_alloc_complex(alloc_local);    // W^(M)
        hwM2  = fftw_alloc_complex(alloc_local);   // W^(M-2)
        hw2M  = fftw_alloc_complex(alloc_local);   // (W^2)^(M)
        hw2M2 = fftw_alloc_complex(alloc_local);  // (W^2)^(M-2)
        
#else
        
        temp1  = fftwf_alloc_real(2 * alloc_local);
        temp2  = fftwf_alloc_real(2 * alloc_local);
        htemp1 = fftwf_alloc_complex(alloc_local);
        htemp2 = fftwf_alloc_complex(alloc_local);
        
        hetan = fftwf_alloc_complex((NLevs+1)*alloc_local);
        hphin = fftwf_alloc_complex((NLevs+1)*alloc_local);
        
        hwn   = fftwf_alloc_complex((NLevs+1)*alloc_local);
        hwM   = fftwf_alloc_complex(alloc_local);
        hwM2  = fftwf_alloc_complex(alloc_local);
        hw2M  = fftwf_alloc_complex(alloc_local);
        hw2M2 = fftwf_alloc_complex(alloc_local);
        
#endif
        
        Coeff = malloc(sizeof(double)*(NLevs+1)); //factorial componant of expansion for phi and W
        Coeff[0] = 1.0;
        Coeff[1] = 1.0;
        
        for (n=2; n<=NLevs; n++)
        {
            Coeff[n] = Coeff[n-1]*n;
        }
        
        
    }
    
    void rhs_hos_clean(){
        
        fftw_free(temp1);
        fftw_free(temp2);
        fftw_free(htemp1);
        fftw_free(htemp2);
        fftw_free(hetan);
        fftw_free(hphin);
        fftw_free(hwn);
        fftw_free(hwM);
        fftw_free(hwM2);
        fftw_free(hw2M);
        fftw_free(hw2M2);
        
        free(Coeff);
        
    }
    /* HOS scheme (West et al. 1987) */
#ifdef USE_DOUBLES
void rhs_hos(fftw_complex* hrhs, fftw_complex* hu, float t){
#else
void rhs_hos(fftwf_complex* hrhs, fftwf_complex* hu, float t){
#endif
            
            
            /* NOTE: hu = [heta; hphi] */
            
            ptrdiff_t index, index_shift;
            
            
            /*------------------*/
            /* RHS for eta part */
            
            /* eta_x*phi_x + eta_y*phi_y */
            Dx(hu, htemp1);
            Dx(&hu[alloc_local], htemp2);
            Mult(htemp1, htemp2, hrhs);
            
            Dy(hu, htemp1);
            Dy(&hu[alloc_local], htemp2);
            Mult(htemp1, htemp2, htemp1);
            
            for (i=0; i<local_N; i++)
            {
                hrhs[i] = - hrhs[i] - htemp1[i];
            }
            
            /* Vertical velocity */
            Zvel(hu, hwM, hwM2, hw2M, hw2M2, t);
            
            Dx(hu, htemp1);                /* eta_x --> htemp1 */
            Mult(htemp1, htemp1, htemp1);  /* eta_x*eta_x --> htemp1 */
            
            Dy(hu, htemp2);                /* eta_y --> htemp2 */
            Mult(htemp2, htemp2, htemp2);  /* eta_y*eta_y --> htemp2 */
            
            Sum(htemp1, htemp2, htemp1);   /* eta_x*eta_x + eta_y*eta_y --> htemp1  */
            Mult(htemp1, hwM2, htemp2);    /* W^(M-2)*(eta_x*eta_x + eta_y*eta_y) --> htemp2  */
            Sum(hwM, htemp2, htemp2);      /* W^(M-2)*(eta_x*eta_x + eta_y*eta_y) + W^(M) --> htemp2  */
            
            Sum(hrhs, htemp2, hrhs);
            
            
            /*------------------*/
            /* RHS for phi part */
            Mult(htemp1, hw2M2, htemp2);   /* W2^(M-2)*(eta_x*eta_x + eta_y*eta_y) --> htemp2  */
            Sum(hw2M, htemp2, htemp2);     /* W2^(M-2)*(eta_x*eta_x + eta_y*eta_y) + W2^(M) --> htemp2  */
            
            for (i=0; i<local_N; i++)
            {
                index_shift = i + alloc_local;
                hrhs[index_shift] = -g*hu[i];
                hrhs[index_shift] = hrhs[index_shift] + 0.5*htemp2[i];
            }
            
            /* -0.5*(phi_x*phi_x + phi_y*phi_y) */
            Dx(&hu[alloc_local], htemp1);
            Mult(htemp1, htemp1, htemp1);
            
            Dy(&hu[alloc_local], htemp2);
            Mult(htemp2, htemp2, htemp2);
            
            /* phi_x*phi_x + phi_y*phi_y --> htemp1  */
            Sum(htemp1, htemp2, htemp1);
            
            for (i=0; i<local_N; i++)
            {
                index = i;
                index_shift = index + alloc_local;
                hrhs[index_shift] = hrhs[index_shift] - 0.5*htemp1[index];
            }
            
            
            /* Dommermuth initialiation scheme -----------------------*/
            /* Multiply nonlinear part of rhs by ramping coefficient. */
            if (rampflg==1)
            {
                ZvelLinear(hu, htemp1);
                
                for (i=0; i<local_N; i++)
                {
                    hrhs[i] -= htemp1[i];
                    hrhs[i] *= RampFun(t);
                    hrhs[i] += htemp1[i];
                }
                
                for (i=0; i<local_N; i++)
                {
                    index_shift = i + alloc_local;
                    hrhs[index_shift] += g*hu[i];
                    hrhs[index_shift] *= RampFun(t);
                    hrhs[index_shift] -= g*hu[i];
                }
            }
            
            /* Dealiasing --------------------------------------------*/
            Dealias(hrhs);
            Dealias(&hrhs[alloc_local]);
            
        }
        
        /**
         * @brief computes the vertical velocity on the free surface using the HOS scheme (West et al., 1987)
         *
         * @param hu a fftw_complex array pointer
         * @param hZvelM a fftw_complex array ponter
         * @param hZvelM2 a fftw_complex array pointer
         * @param hZvel2M a fftw_complex array pointer
         * @param hZvel2M2 a fftw_complex array pointer
         * @param t a TYPE_REAL
         */
        
        /* Compute vertical velocity on the free surface. */
        /* HOS scheme (West et al. 1987) */
#ifdef USE_DOUBLES
void Zvel(fftw_complex* hu, fftw_complex* hZvelM, fftw_complex* hZvelM2, fftw_complex* hZvel2M, fftw_complex* hZvel2M2, float t){
#else
void Zvel(fftwf_complex* hu, fftwf_complex* hZvelM, fftwf_complex* hZvelM2, fftwf_complex* hZvel2M, fftwf_complex* hZvel2M2, float t){
#endif
                
                // hZvelM   -> W^(M)
                // hZvelM2  -> W^(M-2)
                // hZvel2M  -> (W^2)^(M)
                // hZvel2M2 -> (W^2)^(M-2)
                
                ptrdiff_t Nindex, Mindex;
                ptrdiff_t n, m;
                
                
                
                for (n=0; n<=NLevs; n++) {
                    
                    for (i=0; i<local_N; i++) {
                        
                        hwn[n*alloc_local + i] = 0.0;
                        hphin[n*alloc_local + i] = 0.0;
                        
                    }
                    
                }
                
                for (i=0; i<local_N; i++) {
                    
                    hZvelM[i]   = 0.0;
                    hZvelM2[i]  = 0.0;
                    hZvel2M[i]  = 0.0;
                    hZvel2M2[i] = 0.0;
                    hetan[i]    = 0.0;
                    
                }
                
                if (mpi_rank==0) {
                    
                    /* heta^0 = 1*/
                    hetan[0] = 1.0;
                }
                
                Dz( &hu[alloc_local], hphin); /* hphin[0] -> d/dz psi */
                
                
                for (i=0; i<local_N; i++) {
                    
                    hwn[i] = hphin[i];
                    
                    /* hwn[0] -> d/dz psi */
                    
                    
                    
                }
                
                for (n=1; n<=NLevs; n++) {
                    
                    /* eta^n */
                    /* Construct array containing successive powers of heta in sequence */
                    
                    Nindex = n*alloc_local;
                    
                    Mult(hu, &hetan[(n-1)*alloc_local], &hetan[Nindex]);
                    
                    if (mpi_rank == 0){
                        
                        
                    }
                    
                    /* Phi_n */
                    /* Taylor Series for Phi (n -> nth term in series) */
                    
                    for (m=0; m<=n-1; m++) {
                        
                        Mindex = m*alloc_local;
                        
                        
                        Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
                        
                        
                        for (i=0; i<local_N; i++) {
                            
                            
                            hphin[Nindex + i] = hphin[Nindex + i] - htemp1[i]/Coeff[n-m];
                            
                            
                        }
                        
                    }
                    
                    /* W_n */
                    /* Taylor Series for W (n -> nth term in series) */
                    
                    for (m=0; m<=n; m++) {
                        
                        Mindex = m*alloc_local;
                        
                        
                        
                        /* odd derivatives (1st/3rd/5th...) dz hu -> kn*tanh(kn*h)*hu */
                        /* even derivatives (2nd/4th/6th...) dz hu -> kn*hu */
                        
                        if ( n % 2){                   //(n) odd
                            Dz_trem(&hphin[Mindex], &hphin[Mindex]);
                        }
                        
                        else{
                            Dz( &hphin[Mindex], &hphin[Mindex]);
                        }
                        
                        
                        Mult( &hphin[Mindex], &hetan[Nindex - Mindex], htemp1);
                        
                        
                        /* Taylor series computation */
                        for (i=0; i<local_N; i++) {
                            
                            hwn[Nindex + i] = hwn[Nindex + i] + htemp1[i]/Coeff[n-m];
                            
                            
                        }
                        
                        
                    }
                    
                }
                
                /* hZvelM2 */
                for (n=0; n<=(NLevs-2); n++) {  //W*(Deta)^2 3rd order, included for Nlev >= 2
                    
                    for (i=0; i<local_N; i++) {
                        
                        hZvelM2[i] = hZvelM2[i] + hwn[n*alloc_local + i];
                        
                    }
                    
                }
                
                /* hZvelM */
                for (i=0; i<local_N; i++) {
                    
                    hZvelM[i] = hZvelM2[i];
                    
                }
                
                
                for (n=(NLevs-1); n<=NLevs; n++) {
                    
                    if (n >= 0){  //this makes sure W^(-n) isn't accessed
                        
                        for (i=0; i<local_N; i++) {
                            
                            
                            hZvelM[i] = hZvelM[i] + hwn[n*alloc_local + i];
                            
                        }
                        
                        
                        
                    }
                }
                
                /* hZvel2M2 */
                for (n=0; n<=(NLevs-3); n++) { //W^2*(Deta)^2 4th order, nonzero for Nlev >= 3
                    
                    
                    /* W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0) */
                    /* West et al (1987)'s formulation for keeping ordering consistent wrt epsilon */
                    
                    for (m=0; m<=n; m++) {
                        
                        
                        
                        Mult(&hwn[m*alloc_local], &hwn[(n-m)*alloc_local], htemp1);
                        
                        
                        Sum(hZvel2M2, htemp1, hZvel2M2);
                        
                    }
                    
                    //Sum(hZvel2M2, htemp2, hZvel2M2);
                    
                }
                
                
                /* hZvel2M */
                for (i=0; i<local_N; i++) { //W^2 2nd order, nonzero for Nlev >= 1
                    
                    hZvel2M[i] = hZvel2M2[i];
                    
                }
                
                for (n=(NLevs-2); n<=(NLevs-1); n++) {
                    
                    
                    
                    /* W2^(n) = W^(0)*W^(n) + W^(1)*W^(n-1) + ... W^(n)*W^(0) */
                    for (m=0; m<=n; m++) { //if n < 0, this loop won't run (ie, Nlevs < 1)
                        //This forces the exclussion of W^(-n)
                        
                        Mult(&hwn[m*alloc_local], &hwn[(n-m)*alloc_local], htemp1);
                        
                        Sum(hZvel2M, htemp1, hZvel2M);
                        
                    }
                    
                    //        Sum(hZvel2M, htemp2, hZvel2M);
                    
                }
                
                
            }
            /* Compute vertical velocity on the free surface at linear order. */
#ifdef USE_DOUBLES
    void ZvelLinear(const fftw_complex* hu, fftw_complex* hZvelLinear){
#else
    void ZvelLinear(const fftwf_complex* hu, fftwf_complex* hZvelLinear){
#endif
                
        Dz(&hu[alloc_local], hZvelLinear);
                
    }
            
            /* Compute the Hamiltonian for the HOS system, as given in West et al. 1987, eq. 8.*/
#ifdef USE_DOUBLES
float Hamiltonian(const fftw_complex* heta, const fftw_complex* heta_t, const fftw_complex* hphi){
#else
float Hamiltonian(const fftwf_complex* heta, const fftwf_complex* heta_t, const fftwf_complex* hphi){
#endif
                    
                    ptrdiff_t index;
                    double H;
                    
                    H = 0.0;
                    
                    /* j=0 */
                    for (i=0; i < Nx; i++)
                    {
                        H += g*creal(heta[i])*creal(heta[i]);
                        H += g*cimag(heta[i])*cimag(heta[i]);
                        
                        H += creal(hphi[i])*creal(heta_t[i]);
                        H += cimag(hphi[i])*cimag(heta_t[i]);
                    }
                    
                    for (j=1; j<local_Nyhpo; j++)
                    {
                        for (i=0; i < Nx; i++)
                        {
                            index = Nx*j + i;
                            
                            H += 2*g*creal(heta[index])*creal(heta[index]);
                            H += 2*g*cimag(heta[index])*cimag(heta[index]);
                            
                            H += 2*creal(hphi[index])*creal(heta_t[index]);
                            H += 2*cimag(hphi[index])*cimag(heta_t[index]);
                        }
                    }
                    
                    return H;
                    
                }
                
/* Ramping function used for Dommermuth adjusted initialization (D. Dommermuth, Wave Motion, Vol 32, 2000). */
#ifdef USE_DOUBLES
    double RampFun(const double t){
#else
    float RampFun(const float t){
#endif
                    
                    return ( 1 - exp(-t*t/Tramp/Tramp) );
                    
                }
                
                
