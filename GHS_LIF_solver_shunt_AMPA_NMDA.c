/*=================================================================
 *
 * GHS_LIF_solver_shunt.C
 *
 *
 * Mark Humphries & Rob Stewart (6/1/2005)
 *=================================================================*/

#include <math.h>
#include "mex.h"

/* Input Arguments */
#define IN_N    prhs[0] 
#define IN_T    prhs[1]
#define LN_IN   prhs[2]
#define AMPA_IN   prhs[3]
#define NMDA_IN   prhs[4]
#define GABAA_IN   prhs[5]
#define LT_IN   prhs[6]
#define MAX_P   prhs[7]
#define MAX_S   prhs[8]
#define IG      prhs[9]
#define DELAYS  prhs[10]
#define BOUNDS  prhs[11]
#define N_PATHS prhs[12]
#define SPON    prhs[13]
#define T_ONSET prhs[14]
#define T_OFF   prhs[15]
#define P_CELLS prhs[16]
#define B_T1    prhs[17]
#define B_T2    prhs[18]
#define A_CA    prhs[19]
#define THETA_CA prhs[20]
#define C_CELLS prhs[21]
#define THETA   prhs[22]
#define MLIMIT  prhs[23]
#define TAU_AMPA   prhs[24]
#define TAU_NMDA   prhs[25]
#define TAU_GABAA   prhs[26]
#define TAU_M   prhs[27]
#define R_IN    prhs[28]

#define SCALAR_WS  prhs[29] /* weights = [SD1_w SD2_w STN_GPew STN_GPiw GPe_STNw GPe_GPiw GPe_GPew GPi_GPiw EXT_w STN_ext_ratio]; */
                            /* recall elements numbered from 0... For future expansion */

#define DPS     prhs[30]
#define IPS     prhs[31]
#define INIT_M  prhs[32]
#define INIT_I_GABAA  prhs[33]
#define INIT_I_SOMA  prhs[34]
#define INIT_I_PROX  prhs[35]
#define INIT_I_AMPA  prhs[36]
#define INIT_I_NMDA  prhs[37]

/* Output Arguments */
#define OUT_N       plhs[0]
#define OUT_T       plhs[1]
#define N_OUT       plhs[2]
#define N_EVENTS    plhs[3]
#define TRACE_VALS  plhs[4]
#define MEM_POT     plhs[5]
#define I_GABAA     plhs[6]
#define I_PROX      plhs[7]
#define I_SOMA      plhs[8]
#define I_AMPA      plhs[9]
#define I_NMDA      plhs[10]

/*-------------------- RANDOM NUMBER GENERATOR -----------------------------------*/
#define SHR3 (jz=jsr, jsr^=(jsr<<13), jsr^=(jsr>>17), jsr^=(jsr<<5),jz+jsr)
#define UNI (.5 + (signed) SHR3 * .2328306e-9)
#define RNOR (hz=SHR3, iz=hz&127, (abs(hz)<kn[iz])? hz*wn[iz] : nfix())
#define REXP (jz=SHR3, iz=jz&255, ( jz <ke[iz])? jz*we[iz] : efix())

static unsigned int iz,jz,jsr=123456789,kn[128],ke[256];
static int hz; static float wn[128],fn[128], we[256],fe[256];

float nfix(void) { /*provides RNOR if #define cannot */
    const float r = 3.442620f; static float x, y;
    for(;;){ x=hz*wn[iz];
        if(iz==0){ do{x=-log(UNI)*0.2904764; y=-log(UNI);} while(y+y<x*x);
        return (hz>0)? r+x : -r-x;
        }
        if( fn[iz]+UNI*(fn[iz-1]-fn[iz]) < exp(-.5*x*x) ) return x;
        hz=SHR3; iz=hz&127;if(abs(hz)<kn[iz]) return (hz*wn[iz]);
        } 
    }
    
float efix(void) { /*provides REXP if #define cannot */
    float x; 
    for(;;){
        if(iz==0) return (7.69711-log(UNI));
        x=jz*we[iz];
        if( fe[iz]+UNI*(fe[iz-1]-fe[iz]) < exp(-x) ) return (x);
        jz=SHR3; iz=(jz&255);
        if(jz<ke[iz]) return (jz*we[iz]);
    } 
}

/*--------This procedure sets the seed and creates the tables------*/
void zigset(unsigned int jsrseed) {
    const double m1 = 2147483648.0, m2 = 4294967296.;
    double dn=3.442619855899,tn=dn,vn=9.91256303526217e-3, q;
    double de=7.697117470131487, te=de, ve=3.949659822581572e-3;
    int i; jsr=jsrseed;

    /* Tables for RNOR: */ 
    q=vn/exp(-.5*dn*dn);
    kn[0]=(dn/q)*m1; kn[1]=0;
    wn[0]=q/m1; wn[127]=dn/m1;
    fn[0]=1.; fn[127]=exp(-.5*dn*dn);
    for(i=126;i>=1;i--) {
        dn=sqrt(-2.*log(vn/dn+exp(-.5*dn*dn)));
        kn[i+1]=(dn/tn)*m1; tn=dn;
        fn[i]=exp(-.5*dn*dn); wn[i]=dn/m1; 
    }
    
    /* Tables for REXP */ 
    q = ve/exp(-de);
    ke[0]=(de/q)*m2; ke[1]=0;
    we[0]=q/m2; we[255]=de/m2;
    fe[0]=1.; fe[255]=exp(-de);
    for(i=254;i>=1;i--) {
        de=-log(ve/de+exp(-de));
        ke[i+1]= (de/te)*m2; te=de;
        fe[i]=exp(-de); we[i]=de/m2; 
    }
}


/*------------------------------- SIMULATION CODE STARTS HERE -----------------------------------*/
static void build_ext( /* //doesn't look optimal - 04/02/04 */
                unsigned int    *extn,
                unsigned int    *extp,
                unsigned int    *n_ext,
                unsigned int    *n_paths,
                unsigned int    *delays,
                unsigned int    n_sources,
                unsigned int    p,
                unsigned int    tc,
                unsigned int    h,
                unsigned int    MaxEx
              )
{
    int path, newt;
    for(path = p; path < (n_paths[p]*n_sources+p); path+=n_sources){
      newt = (delays[path] + tc)%h; /* //circular addressing */
      extn[n_ext[newt]*h+newt] = p;      /* //store source */
      extp[n_ext[newt]*h+newt] = path;  /* //store path */
        n_ext[newt]++;
        if (n_ext[newt] >= MaxEx) mexErrMsgTxt("Bounding error. Increase MaxEx");
    }
}

static void IntFire3_GHS( 
                unsigned int    *out_n,
                unsigned int    *out_t,
                unsigned int    *n_out,
                unsigned int    *n_events,
                double          *trace_vals,

		double          *m,           /* initialises neuron state variables     */
		double          *gabaa,       /* and returns last values                */
		double          *ip,          /* so that network can be started          */
		double          *is,          /* using some snapshot                     */ 
		double          *ampa,        /* (excluding both refractory period and   */
		double          *nmda,        /* burst firing status)                    */

                unsigned int    *in_n,
                unsigned int    *in_t,
                mxArray         *LN,
                mxArray         *L_AMPA,
                mxArray         *L_NMDA,
                mxArray         *L_GABAA,
                mxArray         *LT,
                double          *max_prox,
                double          *max_soma,    
                double          *i_gate,
                unsigned int    *delays,
                unsigned int    *bounds,
                unsigned int    *n_paths,
                double          *spon,
                unsigned int    *t_onset,
                unsigned int    *t_off,
                unsigned int    *p_cells,
                unsigned int    *burst_t1, 
                unsigned int    *burst_t2, 
                double          *alphaCA, 
                double          *thetaCA,
                unsigned int    *ca_cells,
                double          *theta,
                double          *mlimit,
                double          *tau_AMPA,
                double          *tau_NMDA,
                double          *tau_GABAa,
                double          *tau_m,
                double          *R,
                double          *Wgts,
                double          *dps,
                unsigned int    *ips,
		double          *init_mem_pot, 
		double          *init_i_gabaa, 
		double          *init_i_soma, 
		double          *init_i_prox, 
		double          *init_i_ampa, 
		double          *init_i_nmda
              )
{
    /* unbundle parameters */
    double sigma = dps[0];     /* //Noise std */
    double ref = dps[1];       /* //Refractory membrane potential reset value */
    double dt = dps[2];
    double step_size = dps[3];
    double PSP_sigma = dps[4];
    double M_hat = dps[5];
    
    unsigned int h = ips[0];
    unsigned int MaxEx = ips[1];
    unsigned int MaxOut = ips[3];
    unsigned int n_neurons = ips[4];
    unsigned int n_sources = ips[5];
    unsigned int n_in = ips[6];
    unsigned int time_steps = ips[7];
    unsigned int ref_period = ips[8];
    unsigned int trace = ips[9];
    
    unsigned int n_ext_total = 0;
    unsigned int n_out_val = 0;
    unsigned int tc = 0;  
    unsigned int in = 0;
    unsigned int t,p,loop,loop2,source,first,last;
    int idum[] = {-1};
    double w, psp_noise, pre_m, soma_gate, prox_gate, total_gate;
    
    unsigned int p_index = 0;
    
    /* mxArray *AMPA, *NMDA, *GABAA, *IP, *IS, *M, */
    mxArray  *I_EXT, *I_CA, *B_END, *T_CA, *T1, *E_AMPA, *E_NMDA, *E_GABAA, *EM, *AS, *EXTN, *EXTP, *N_EXT;
    /* double *ampa, *nmda, *gabaa, *m, *ip, *is; */
    double *i_ext, *i_ca, *as, *e_ampa, *e_nmda, *e_gabaa, *em, *lw_ampa, *lw_nmda, *lw_gabaa;
    unsigned int *extn, *extp, *n_ext, *ln, *lt;
    unsigned int type;
    int *t1, *t_ca, *b_end;
    
    int t1_dims[] = {1,0};
    int en_dims[] = {1,0};
    int ne_dims[] = {1,0};

    t1_dims[1] = n_neurons;
    en_dims[0] = h; en_dims[1] = MaxEx;
    ne_dims[1] = h;
    
    /*Create state variable arrays*/
    /* AMPA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    NMDA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    GABAA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    IP = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    IS = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    M = mxCreateDoubleMatrix(1,n_neurons,mxREAL); */

    I_EXT = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    I_CA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    AS = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    E_AMPA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    E_NMDA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    E_GABAA = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    EM = mxCreateDoubleMatrix(1,n_neurons,mxREAL);
    B_END = mxCreateNumericArray(2,t1_dims,mxINT32_CLASS,mxREAL);
    T_CA = mxCreateNumericArray(2,t1_dims,mxINT32_CLASS,mxREAL);
    T1 = mxCreateNumericArray(2,t1_dims,mxINT32_CLASS,mxREAL);
    
    /*assign pointers to state variables*/
    /* ampa = mxGetPr(AMPA);
    nmda = mxGetPr(NMDA);
    gabaa = mxGetPr(GABAA);
    ip = mxGetPr(IP);
    is = mxGetPr(IS);
    m = mxGetPr(M); */
    i_ext = mxGetPr(I_EXT);
    i_ca = mxGetPr(I_CA);
    as = mxGetPr(AS);
    e_ampa = mxGetPr(E_AMPA);
    e_nmda = mxGetPr(E_NMDA);
    e_gabaa = mxGetPr(E_GABAA);
    em = mxGetPr(EM);
    b_end = (int*)mxGetData(B_END);
    t_ca = (int*)mxGetData(T_CA);
    t1 = (int*)mxGetData(T1);
    
    /*Create queue variables*/
    EXTN = mxCreateNumericArray(2, en_dims, mxUINT32_CLASS, mxREAL);
    EXTP = mxCreateNumericArray(2, en_dims, mxUINT32_CLASS, mxREAL);
    N_EXT = mxCreateNumericArray(2, ne_dims, mxUINT32_CLASS, mxREAL);
    
    /*assign pointers to queue variables*/
    extn = mxGetData(EXTN);  
    extp = mxGetData(EXTP);
    n_ext = mxGetData(N_EXT);
    
     /*----- set random number seed ---------*/
    zigset(1);
    
    /* //Initialise t1, ee, ei, em, as */
    for(p = 0; p < n_neurons; p++){
        t_ca[p] = -1;
        t1[p] = -1;
        e_ampa[p] = exp(-dt/tau_AMPA[p]);
        e_nmda[p] = exp(-dt/tau_NMDA[p]);
        e_gabaa[p] = exp(-dt/tau_GABAa[p]);
        em[p] = exp(-dt/tau_m[p]);
        as[p] = R[p]*(1-em[p]);
        i_ca[p] = 0;
        i_ext[p] = 0;
        b_end[p] = burst_t1[p] + burst_t2[p]; /* // sum these here for efficiency */
        /* // initialise membrane and synaptic variables (history) */
        m[p] = init_mem_pot[p]; /* RNOR * 0.002 + 0.00; */
        ip[p] = init_i_prox[p];
        is[p] = init_i_soma[p];
        gabaa[p] = init_i_gabaa[p]; /* RNOR * 0.0001 - 0.0001; */
        ampa[p] = init_i_ampa[p];  /* RNOR * 0.0001 + 0.0001; */
        nmda[p] = init_i_nmda[p];  /* RNOR * 0.0001 + 0.0001; */
        /* //printf("Initial is %f\n",m[p]); */	
     }    
    
    
    /*Begin main loop over time steps*/
    
    printf("Tracing neuron... %d \n",trace);
    for(t = 0; t < time_steps; t++){
        /* //first trace update 'trace' is theneuron index you specified in matlab */
        /* // *************** customise detailed traces here **************** */
        trace_vals[t] = gabaa[trace]; /* distal inhibtory */
        trace_vals[time_steps + t] = ampa[trace] + nmda[trace]; /* total excitatory current */
        trace_vals[time_steps*2 + t] = m[trace];                /* membrane potential */
        /* trace_vals[time_steps + t] = gabaa[trace]; ** distal inhibitory */
        /* trace_vals[time_steps + t] = ip[trace]; ** proximal inhibitory */
        /* trace_vals[time_steps + t] = is[trace]; ** somatic inhibitory */
        /* trace_vals[time_steps + t] = i_ca[trace]; ** Ca++ inhibitory */
        
      for(p = 0; p < n_neurons; p++){ /* //loop over all neurons */
	
	    pre_m = m[p];	/* // membrane potential on previous time-step was */
        	
            if ((t-t1[p]) > ref_period) { 
	      /* //printf("Somatic gate %f \n",1 - (is[p] / max_soma[p])); */
                
                /* Option 1: scale gate size by max for that neuron */
                /* soma_gate = 1 - (is[p] / max_soma[p]);
                prox_gate = 1 - (ip[p] / max_prox[p]);  */
                
                /* Option 2: scale gate size by max for the network */
                soma_gate = 1 - (is[p] / M_hat);
                prox_gate = 1 - (ip[p] / M_hat);
                total_gate = 1 - (prox_gate + soma_gate) / 2; 

                /* // hard-limit gates */
                if (soma_gate < 0) soma_gate = 0;
                if (prox_gate < 0) prox_gate = 0;
                if (total_gate < 0) total_gate = 0;
                
                /* trace gate value */
        if (p == trace) trace_vals[time_steps + t] = total_gate;  /* // trace gating variables */
                
                
		/* // ********** The membrane equation ****************** */
        
                /* // gate inputs and spontaneous current */
                /* //m[p] = m[p]*em[p] + ((((e[p] + i[p]) * prox_gate)  + spon[p] + i_ca[p]) * soma_gate) + i_ext[p]; // gate all psps */
                
                /* // no gating, just increase size of clamping current  */
                /* m[p] = m[p]*em[p] + ampa[p] + nmda[p] + gabaa[p] + spon[p] + i_ca[p] + i_ext[p] + i_gate[p] * total_gate;
                
                /* gating and clamp current combined - shunt I_ca */
               /* m[p] = m[p]*em[p] + (((ampa[p] + nmda[p] + gabaa[p]) * prox_gate)  + i_ca[p]) * soma_gate + i_ext[p] + spon[p] + i_gate[p] * total_gate; */

	       /* gating and clamp current combined - don't shunt shunt I_ca */
               m[p] = m[p]*em[p] + (((ampa[p] + nmda[p] + gabaa[p]) * prox_gate)) * soma_gate + i_ext[p] + spon[p] + i_ca[p] + i_gate[p] * total_gate;
                
                /* // add noise - do this here because now have two thresholds etc */
                m[p] = m[p] + RNOR * sigma;
                
                /* // hard-limit to replicate reversal potential  */
                if(m[p] < mlimit[p]) m[p] = mlimit[p];
            }
            
            /* // check for external current input */
            if (p_cells[p]==1 & t > t_onset[p_index] & t < t_off[p_index]) {
	      /* // printf("p_cell %f\n",p_cells[p]); */
                i_ext[p] = step_size * as[p];
            }
            else i_ext[p] = 0;
            
            /* // do Ca current checks - is burst cell, not bursting, and rebounds */
            if (ca_cells[p]==1 & i_ca[p]==0 & pre_m < thetaCA[p] & m[p] > thetaCA[p]){    /* // rebound bursting */
	      /* //printf("burst triggered at %d by pre m %f and m %f \n",t,pre_m,m[p]); */
                  t_ca[p] = 1;
                  i_ca[p] = alphaCA[p];
            }
            else{
                if (ca_cells[p]==1 & i_ca[p] > 0){
                      t_ca[p]++;
                      /* //printf("t CA %d\n",t_ca[p]); */
                      if (t_ca[p] >= burst_t1[p] & t_ca[p] <= b_end[p]){
			/* // calculation for ramp - move constant calculation to preamble for efficiency */
                          i_ca[p] = alphaCA[p] - (alphaCA[p] / burst_t2[p])*(t_ca[p] - burst_t1[p]);
                      }
                      else{
                          if(t_ca[p] > b_end[p]){  
                             t_ca[p] = 0;
                             i_ca[p] = 0;
                          }
                      }
                }
            }        
            /* // update PSP totals */
            gabaa[p] *= e_gabaa[p]; 
            ip[p] *= e_gabaa[p];
            is[p] *= e_gabaa[p];
            ampa[p] *= e_ampa[p];
            nmda[p] *= e_nmda[p];

            
            /* // ******************** end membrane potential ****************  */
        }
        
      
            
      /* //Process external events */
        for(loop = tc; loop < (n_ext[tc]*h+tc); loop+=h){
            source = extn[loop];
            ln = mxGetData(mxGetCell(LN,source));   /* //get pointer to target neurons */
            lw_ampa = mxGetPr(mxGetCell(L_AMPA,source));     /* //get pointer to connection weights */            
            lw_nmda = mxGetPr(mxGetCell(L_NMDA,source));     /* //get pointer to connection weights */            
            lw_gabaa = mxGetPr(mxGetCell(L_GABAA,source));     /* //get pointer to connection weights */            
            lt = mxGetData(mxGetCell(LT,source));     /* //get pointer to connection type */
            first = bounds[extp[loop]];
            last = bounds[extp[loop]+n_sources];
            for(loop2 = first; loop2 < last; loop2++){
                n_ext_total++;
                p = ln[loop2];      /* //Neuron */
                if (p != 0) {       /* // then is connnected */
			        type = lt[loop2];   /* //Type                 */
			        psp_noise =  RNOR * PSP_sigma * w;	/* // scale PSP noise by weight */
			        /* //psp_noise = 0; */
              
					  switch(type){
					    case 0: ampa[p] += lw_ampa[loop2] + psp_noise;          /* //excitatory event */
					            nmda[p] += lw_nmda[loop2] + psp_noise;
                                            break;
					    case 1: gabaa[p] += lw_gabaa[loop2] + psp_noise;          /* //inhibitory event - distal */
                                            break;
					    case 2: ip[p] += lw_gabaa[loop2] + psp_noise;         /* //inhibitory event - proximal */
                                            break;
					    case 3: is[p] += lw_gabaa[loop2] + psp_noise;         /* //inhibitory event - soma */
                                            break;
                        default: printf("Shouldn't be here \n");
                      }
                }                
            } /* //loop2     */
        } /* //loop - external events */
        
        /* // check if current has been injected - if so, increment time array */
        if (t>0 & fmod(t,t_off[p_index]) == 0){
	  /* //printf("On loop %d t_onset is %f, t_off %f\n",t,t_onset[p_index],t_off[p_index]); */
            p_index++;
        }
        
        n_ext[tc] = 0; /* //reset number of external events */
        /* //printf("PSP noise.. %f \n",psp_noise); */
       
        for(p = 0; p < n_neurons; p++){ /* //fire loop */
	  /* //if (m[p] > (theta[p] + (RNOR * sigma ))){ //fire if above threshold + noise */
	  if (m[p] > theta[p]) { /* //fire if above threshold - noise now added to membrane potential */
                if (p == trace)
		  trace_vals[time_steps*2 + t] = m[trace]+0.04; /* //fake spikes */
                m[p] = ref;
                t1[p] = t;
                out_n[n_out_val] = p+1;
                out_t[n_out_val] = t+1;
                n_out_val++;
                if (n_out_val > MaxOut) mexErrMsgTxt("Bounding error. Increase MaxOut");
                build_ext(extn,extp,n_ext,n_paths,delays,n_sources,p,tc,h,MaxEx);
            }   
        } /* //fire loop */

        /* //Add extrinsic events to queue */
        while ((in_t[in]) == t){
            build_ext(extn,extp,n_ext,n_paths,delays,n_sources,in_n[in],tc,h,MaxEx);
            if (++in >= n_in) break;
        }
        
        tc++;
        if (tc == h) tc = 0; /* //circular addressing */
    }
    n_out[0] = n_out_val;
    n_events[0] = n_ext_total;
    printf("number of spike sources: \t%d\n",n_sources);
    printf("number of events: \t\t\t%d\n",n_ext_total);
    printf("number of spikes: \t\t\t%d\n",n_out_val);
}

void mexFunction( int nlhs, mxArray *plhs[], 
          int nrhs, const mxArray *prhs[] )
{
    double *dps, *spon, *alphaCA, *thetaCA, *theta, *mlimit, *tau_AMPA, *tau_NMDA, *tau_GABAa, *tau_m, *R, *trace_vals, *max_prox, *max_soma, *i_gate, *Wgts, *mem_pot, *i_gabaa, *i_soma, *i_prox, *i_ampa, *i_nmda, *init_mem_pot, *init_i_gabaa, *init_i_soma, *init_i_prox, *init_i_ampa, *init_i_nmda;
    unsigned int *ips, *in_n, *in_t, *out_n, *out_t, *n_out, *n_events, *delays, *bounds, *n_paths,*p_cells,*t_onset, *t_off, *burst_t1, *burst_t2, *ca_cells;
    int dims[] = {1,0};
    int dims2[] = {1,1};
    
    /* Check for proper number of arguments */
    /*if (nrhs != 31) { 
        mexErrMsgTxt("Thirty-one input arguments required."); 
    } else if (nlhs > 5) {
        mexErrMsgTxt("Too many output arguments."); 
    } */

    /* Assign input pointers - pass cell arrays unprocessed */
    in_n = mxGetData(IN_N);
    in_t = mxGetData(IN_T);
    max_prox = mxGetPr(MAX_P);
    max_soma = mxGetPr(MAX_S);
    i_gate = mxGetPr(IG);
    delays = mxGetData(DELAYS);
    bounds = mxGetData(BOUNDS);
    n_paths = mxGetData(N_PATHS);
    spon = mxGetPr(SPON);
    t_onset = mxGetData(T_ONSET);
    t_off = mxGetData(T_OFF);
    p_cells = mxGetData(P_CELLS);
    burst_t1 = mxGetData(B_T1);
    burst_t2 = mxGetData(B_T2);
    alphaCA = mxGetPr(A_CA);
    thetaCA = mxGetPr(THETA_CA);
    ca_cells = mxGetData(C_CELLS);
    theta = mxGetPr(THETA);
    mlimit = mxGetPr(MLIMIT);
    tau_AMPA = mxGetPr(TAU_AMPA);
    tau_NMDA = mxGetPr(TAU_NMDA);
    tau_GABAa = mxGetPr(TAU_GABAA);
    tau_m = mxGetPr(TAU_M);
    R = mxGetPr(R_IN);
    Wgts = mxGetPr(SCALAR_WS);
    dps = mxGetPr(DPS);
    ips = mxGetData(IPS);
    init_mem_pot = mxGetData(INIT_M);
    init_i_gabaa = mxGetData(INIT_I_GABAA);
    init_i_soma = mxGetData(INIT_I_SOMA); 
    init_i_prox = mxGetData(INIT_I_PROX);  
    init_i_ampa = mxGetData(INIT_I_AMPA); 
    init_i_nmda = mxGetData(INIT_I_NMDA);


    /* Create matrices for the return argument */ 
    dims[1] = ips[3]; /* //MaxOut */
    OUT_N = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
    OUT_T = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
    N_OUT = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
    N_EVENTS = mxCreateNumericArray(2, dims2, mxUINT32_CLASS, mxREAL);
    TRACE_VALS = mxCreateDoubleMatrix(ips[7],3,mxREAL);
    MEM_POT = mxCreateDoubleMatrix(ips[4],1,mxREAL);    
    I_GABAA =  mxCreateDoubleMatrix(ips[4],1,mxREAL); 
    I_PROX =  mxCreateDoubleMatrix(ips[4],1,mxREAL); 
    I_SOMA =  mxCreateDoubleMatrix(ips[4],1,mxREAL); 
    I_AMPA =  mxCreateDoubleMatrix(ips[4],1,mxREAL); 
    I_NMDA =  mxCreateDoubleMatrix(ips[4],1,mxREAL); 

    /* Assign output pointers */ 
    out_n = mxGetData(OUT_N);
    out_t = mxGetData(OUT_T);
    n_out = mxGetData(N_OUT);
    n_events = mxGetData(N_EVENTS);
    trace_vals = mxGetPr(TRACE_VALS);
    mem_pot = mxGetData(MEM_POT);
    i_gabaa = mxGetData(I_GABAA);
    i_prox = mxGetData(I_PROX);
    i_soma = mxGetData(I_SOMA);
    i_ampa = mxGetData(I_AMPA);    
    i_nmda = mxGetData(I_NMDA);    


    

    /* Do the actual computations in a subroutine */
    IntFire3_GHS(out_n,out_t,n_out,n_events,trace_vals,mem_pot,i_gabaa,i_soma,i_prox,i_ampa,i_nmda,
		 in_n,in_t,LN_IN,AMPA_IN,NMDA_IN,GABAA_IN,LT_IN,max_prox,max_soma,i_gate,delays,
		 bounds,n_paths,spon,t_onset,t_off,p_cells,burst_t1,burst_t2,alphaCA,thetaCA,
		 ca_cells,theta,mlimit,tau_AMPA,tau_NMDA,tau_GABAa,tau_m,R,Wgts,dps,ips,
		 init_mem_pot, init_i_gabaa, init_i_prox, init_i_soma, init_i_ampa, init_i_nmda); 
}
