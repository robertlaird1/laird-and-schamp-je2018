//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                          //
// intransitivity - Models intransitive competition with competitive outcomes matrices                                      //
// mex function code for MATLAB (type "[mat, world, tsabund, tsrich] = intransitivity(L, s, gmax, n, compmat_tr, frameint)" //
// Pseudorandom numbers generated with Mersenne Twister, seeded from clock and process ID                                   //
// Robert A. Laird, University of Lethbridge, 2016, robert.laird@uleth.ca, except where noted otherwise                     //
//                                                                                                                          //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <mex.h>

using namespace std;

// FUNCTION PROTOTYPES

int rich(const int* abund, const int n); // determines the species richness from an n-member vector of abundances
double twist(void); // gives a random double on [0, 1) using Mersenne Twister

// MAIN FUNCTION

void sim(
        // output (lhs) variables
        double mat_lhs[],           // competition matrix
        double worldmov_lhs[],      // state of the world 
        double tsabund_lhs[],       // time series: abundance values
        double tsrich_lhs[],        // time series: species richness
        
        // input (rhs) variables
        double L_rhs[],             // square-root of population size
        double s_rhs[],             // spatial (1) or aspatial (0)
        double gmax_rhs[],          // maximum number of generations
        double n_rhs[],             // number of species
        double compmat_tr_rhs[],    // top-right of competition matrix
        double frameint_rhs[]       // interval of frames for state of the world
        )
{
    
    // DECLARE USER-GENERATED PARAMETERS AND ANCILLARY VARIABLES 
    
    const int L = static_cast<int>(L_rhs[0]);                   // square root of population size
    const int s = static_cast<int>(s_rhs[0]);                   // interaction type (0 = well mixed, 1 = lattice)
    const int gmax = static_cast<int>(gmax_rhs[0]);             // max number of generations
    const int n = static_cast<int>(n_rhs[0]);                   // number of species
    bool* compmat_tr = new bool[n*(n - 1)/2];                   // top-right of competition matrix in row format 
    for (int i = 0; i < n*(n - 1)/2; i++) {compmat_tr[i] = (compmat_tr_rhs[i] != 0);}
    const int frameint = static_cast<int>(frameint_rhs[0]);     // interval between frames
    
    // DECLARE / INITIALIZE PARAMETERS AND VARIABLES BASED ON USER'S INPUT
    
    // Define and initialize compmat, the full competition matrix (1 means row species wins, 0 means col species wins (or row = col))
    bool* compmat = new bool[n*n]; 
    int ii = 0; // a counter used multiple times in this code
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                compmat[i*n + j] = 0;
            }
            else {
                compmat[i*n + j] = (compmat_tr[ii] != 0);
                compmat[j*n + i] = ((1 - compmat_tr[ii]) != 0);
                ii++;
            }
        }
    }
    
    // Define and intialize world, the competitive arena
    // Define and intialize abund, a vector of the abundance of species 0 to n - 1
    int* world = new int[L*L];
    int* abund = new int[n];
    for (int j = 0; j < n; j++) {abund[j] = 0;}
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            world[i*L + j] = static_cast<int>(twist()*n); // random integer between 0 and n - 1
            abund[world[i*L + j]]++; // initial number of 0s, 1s, 2s, ..., (n - 1)
        }
    }
    
    // Define and initialize worldmov, which records movie frames of the state of the world
    int nframes = static_cast<int>(gmax/frameint + 1); // number of frames to record (will truncate/round down to int)
    int frameID = 0; // gives the number of the frame we are on
    int* worldmov = new int[nframes*L*L]; // movie of the state of the world (equivalent to [nframes, L*L])
    for (int j = 0; j < L*L; j++) {worldmov[frameID*L*L + j] = world[j];} // record initial state of the world
    frameID++; // increment the frame number
    
    // Define and intialize tsabund, the time series for abundance
    int* tsabund = new int[(gmax + 1)*n]; // equivalent to [gmax + 1, n]
    for (int j = 0; j < n; j++) {tsabund[0*n + j] = abund[j];} // intitial abundance conditions
    
    // Define and initialize tsrich, the time series for results
    double* tsrich = new double[(gmax + 1)]; // equivalent to [gmax + 1, 1]
    tsrich[0] = rich(abund, n); // initial species richness
        
    // SIMULATION
    int xr; // first focal individual row
    int xc; // first focal individual col
    int yr; // second focal individual row
    int yc; // second focal individual col
    int yd; // direction of second focal individual (spatial version only)    
    bool monoculture = 0; // set to 1 when/if a monoculture is reached within gmax generations
    if (rich(abund, n) == 1) {monoculture = 1;} // i.e., if there is currently one species
    int g = 0; // index for generation loop    
    while ((monoculture == 0) && (g < gmax)) { // keep looping until monoculture is reached to a maximum of gmax generations
        for (int t = 0; t < L*L; t++) { // time loop
            
            // first focal individual row and col (random integers between 0 and L - 1)
            xr = static_cast<int>(twist()*L);
            xc = static_cast<int>(twist()*L);  
            
            // second focal individual row and col (random integers between 0 and L - 1)
            if (s == 0) { // aspatial
                do {
                    yr = static_cast<int>(twist()*L);
                    yc = static_cast<int>(twist()*L); 
                } while ((yr == xr) && (yc == xc)); // ensures no self neighbours
            }
            else { // spatial (von Neumann neighbourhood; wrapping (toroidal) boundaries)
                yd = static_cast<int>(twist()*4); // random integer between 0 and 3
                if (yd == 0) { // north
                    if (xr == 0) {yr = L - 1;} else {yr = xr - 1;}
                    yc = xc;
                }
                else if (yd == 1) { // west
                    yr = xr;
                    if (xc == 0) {yc = L - 1;} else {yc = xc - 1;}
                }
                else if (yd == 2) { // east
                    yr = xr; 
                    if (xc == L - 1) {yc = 0;} else {yc = xc + 1;}
                }
                else { // south
                    if (xr == L - 1) {yr = 0;} else {yr = xr + 1;}
                    yc = xc;
                }    
            }            
               
            // replacement
            if (compmat[world[xr*L + xc]*n + world[yr*L + yc]] == 1) { // x beats y
                abund[world[xr*L + xc]]++; // increment x's count
                abund[world[yr*L + yc]]--; // decrement y's count
                world[yr*L + yc] = world[xr*L + xc]; // y replaced by clone of x
            }
            else { // y beats x
                abund[world[yr*L + yc]]++; // increment y's count
                abund[world[xr*L + xc]]--; // decrement x's count
                world[xr*L + xc] = world[yr*L + yc]; // x replaced by clone of y
            }
            
        } // end of t loop
        
        // record state of the world (if on a suitable generation)  
        if ((g + 1)%frameint == 0) {
            for (int j = 0; j < L*L; j++) {worldmov[frameID*L*L + j] = world[j];} // record current state of the world
            frameID++; // increment the frame number
        }
        
        // record tsabund data
        for (int j = 0; j < n; j++) {tsabund[(g + 1)*n + j] = abund[j];} // enumerate number of 0s, 1s, ..., (n - 1)'s
        
        // record tsrich data
        tsrich[(g + 1)] = rich(abund, n); // current species richness
                        
        // determine if the population has reached a monoculture
        if (rich(abund, n) == 1) {monoculture = 1;} // i.e., if there is currently one species
        
        g++; // increment generation 
        
    } // end of g loop; end of simulation
    
    // fill in state of the world as a static monoculture for remaining frames
    // fill in rest of both time series
    for (int i = g; i < gmax; i++) {
        if ((i + 1)%frameint == 0) {
            for (int j = 0; j < L*L; j++) {worldmov[frameID*L*L + j] = world[j];} // record current (monoculture) state of the world
            frameID++; // increment the frame number
        }
        for (int j = 0; j < n; j++) {tsabund[(i + 1)*n + j] = abund[j];}
        tsrich[i + 1] = rich(abund, n);    
    }
    
    // output competition matrix
    ii = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            mat_lhs[ii] = compmat[i*n + j]; // convert to column-major order 
            ii++;
        }
    }
    
    //  output worldmov (convert to 1-based indexing)
    ii = 0;
    for (int j = 0; j < L*L; j++) {
        for (int i = 0; i < nframes; i++) {
            worldmov_lhs[ii] = worldmov[i*L*L + j] + 1; // convert to column-major order (+1 converts to 1-based indexing)
            ii++;
        }
    }
    
    // output tsabund time series
    ii = 0;
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < gmax + 1; i++) {
            tsabund_lhs[ii] = tsabund[i*n + j]; // convert to column-major order
            ii++;
        }
    }
    
    // output tsrich time series
    for (int i = 0; i < gmax + 1; i++) {tsrich_lhs[i] = tsrich[i];}
    
    
    // delete pointers
    delete[] compmat_tr; compmat_tr = NULL;
    delete[] compmat; compmat = NULL;
    delete[] world; world = NULL;
    delete[] abund; abund = NULL;
    delete[] worldmov; worldmov = NULL;
    delete[] tsabund; tsabund = NULL;
    delete[] tsrich; tsrich = NULL;
            
}

// SUPPLEMENTARY FUNCTIONS

// rich
int rich(const int* abund, const int n) // determines the species richness from an ns-member vector of abundances
{
    int richness = 0;
    for (int j = 0; j < n; j++) {
        if (abund[j] > 0) richness++;
    }    
    return richness;    
}

// twist
double genrand_real2(void); // function prototype for Mersenne Twister code
double twist(void) // gives a random double on [0, 1) using Mersenne Twister
// Mersenne Twister code given in Appendix 
{
    return genrand_real2();
}

// MEX FUNCTION
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    // arguments
    double* mat_lhs;           // 0. competition matrix
    double* worldmov_lhs;      // 1. state of the world 
    double* tsabund_lhs;       // 2. time series: abundance values
    double* tsrich_lhs;        // 3. time series: species richness
    double* L_rhs;             // 0. square-root of population size
    double* s_rhs;             // 1. spatial (1) or aspatial (0)
    double* gmax_rhs;          // 2. maximum number of generations
    double* n_rhs;             // 3. number of species
    double* compmat_rhs;       // 4. top-right of competition matrix
    double* frameint_rhs;      // 5. interval of frames for state of the world
        
    // simple error checking
    if (nrhs != 6) {mexErrMsgTxt("Function INTRANSITIVITY requires exactly 6 input arguments");}
    if (nlhs > 4) {mexErrMsgTxt("Function INTRANSITIVITY requires 0-4 output arguments");}
    
    // get some useful quantities
    const int n_in = static_cast<int>(*mxGetPr(prhs[3]));
    const int L_in = static_cast<int>(*mxGetPr(prhs[0]));
    const int gmax_in = static_cast<int>(*mxGetPr(prhs[2]));
    const int frameint_in = static_cast<int>(*mxGetPr(prhs[5]));
    const int nframes_in = static_cast<int>(gmax_in/frameint_in + 1);
    
    // more error checking
    if (n_in >= 14) {mexErrMsgTxt("Function INTRANSITIVITY requires n < 13");}
    
    // initiate mrows and ncols, used to size output matrices
    int mrows;
    int ncols;
    
    // prepare zeroth output matrix (mat_lhs)
    mrows = n_in;
    ncols = n_in;
    plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL); 
    
    // prepare first output matrix (worldmov_lhs)
    mrows = nframes_in; 
    ncols = L_in*L_in;
    plhs[1] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
    
    // prepare second output matrix (tsabund_lhs)
    mrows = gmax_in + 1;
    ncols = n_in;
    plhs[2] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
    
    // prepare third output matrix (tsrich_lhs)
    int nresultscols_in = 9;
    mrows = gmax_in + 1;
    ncols = 1;
    plhs[3] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
    
    // define lhs variables
    mat_lhs = mxGetPr(plhs[0]);
    worldmov_lhs = mxGetPr(plhs[1]);
    tsabund_lhs = mxGetPr(plhs[2]);
    tsrich_lhs = mxGetPr(plhs[3]);
    
    // define rhs variables
    L_rhs = mxGetPr(prhs[0]);
    s_rhs = mxGetPr(prhs[1]);
    gmax_rhs = mxGetPr(prhs[2]);
    n_rhs = mxGetPr(prhs[3]);
    compmat_rhs = mxGetPr(prhs[4]);
    frameint_rhs = mxGetPr(prhs[5]);
    
    // run simulation
    sim(mat_lhs, worldmov_lhs, tsabund_lhs, tsrich_lhs,
            L_rhs, s_rhs, gmax_rhs, n_rhs, compmat_rhs, frameint_rhs);
    
}

// APPENDIX - Mersenne Twister code

// Original code by Takuji Nishimura and Makoto Matsumoto (2002)
// Included and modified as allowed in the material below:

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>

// MODIFIED: added to allow for seeding from clock and process ID
#include <ctime>
#include <process.h>
#define GETPID _getpid

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

// OMITTED: "void init_by_array(unsigned long init_key[], int key_length)" 

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(getpid()*time(NULL)); //MODIFIED: Original: init_genrand(5489UL); /* a default initial seed is used */    
            
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

// OMITTED: "long genrand_int31(void)"
// OMITTED: "genrand_real1(void)"

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

// OMITTED: "genrand_real3(void)"
// OMITTED: "genrand_res53(void)"
// OMITTED: "main(void)"