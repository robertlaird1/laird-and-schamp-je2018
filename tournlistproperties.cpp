// When compliled, creates a Matlab function called 'tournlistproperties'
// Call "U = tournlistproperties(tournlist, n);"
// tournlist is a list of tournaments
// n is the number of nodes (i.e., species) in the tournament
//
// Robert Laird, University of Lethbridge, October 2016 (robert.laird@uleth.ca)

#include <mex.h>

using namespace std;

// SUBROUTINES AND FUNCTION PROTOTYPES

void slaters_ij(int* ij, const bool* c, const int n); // finds the distance between a tournament and the nearest hierarchy (Slater's i) and the number of such hierarchies (Slater's j)
void bezembinder(int* bez, const bool* c, const int n); // finds Bezembinders rho and delta_prime
int intrans_triads(const bool* c, const int n); // finds the number of intransitive triads in tournament c   
bool unbeatable(const bool* c, const int n); // determines whether there is an unbeatable species
bool alwaysbeatable(const bool* c, const int n); // determines whether there is a species beaten by every species
long long int fact(const int n); // gives the factorial (works for n <= 20)
void nextperm(int* currperm, const int n); // gives the next permutation of a sequence of integers

void main(double properties_lhs[], double tournlist_rhs[], double n_rhs[])
{
    
    // number of nodes 
    const int n = static_cast<int>(n_rhs[0]);
    
    // number of non-isomorphic tournaments (from Sloane's A000568)
    // when n is...          0, 1, 2, 3, 4,  5,  6,   7,    8,      9,      10
    const int A000568[11] = {1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056};
    
    // get list of tournaments from Matlab
    bool* tournlist = new bool[A000568[n]*n*(n - 1)/2];
    int ii = 0;
    for (int j = 0; j < n*(n - 1)/2; j++) {
        for (int i = 0; i < A000568[n]; i++) {
            tournlist[i*n*(n - 1)/2 + j] = (tournlist_rhs[ii] != 0);
            ii++;
        }
    }
    
    // determine properties for every tournament in the list 
    int numprops = 7; // number of properties to examine
    bool* tourn = new bool[n*(n - 1)/2];
    int ij[2];
    int bez[2];
    int* properties_list = new int[A000568[n]*numprops];
    for (int i = 0; i < A000568[n]; i++) {
        for (int j = 0; j < n*(n - 1)/2; j++) {
            tourn[j] = tournlist[i*n*(n - 1)/2 + j];
        }
        if (n > 2) {
            slaters_ij(ij, tourn, n);
            bezembinder(bez, tourn, n);
            properties_list[i*numprops + 0] = ij[0];
            properties_list[i*numprops + 1] = ij[1];
            properties_list[i*numprops + 2] = intrans_triads(tourn, n);
            properties_list[i*numprops + 3] = static_cast<int>(unbeatable(tourn, n));
            properties_list[i*numprops + 4] = static_cast<int>(alwaysbeatable(tourn, n));
            properties_list[i*numprops + 5] = bez[0];
            properties_list[i*numprops + 6] = bez[1];
        }
        else { // if n <= 2, must be a hierarchy
            properties_list[i*numprops + 0] = 0; // no reversals necessary
            properties_list[i*numprops + 2] = 0; // no intransitive triads 
            if (n == 2) {
                properties_list[i*numprops + 1] = 2; // reversal of only edge begets another hierarchy!
                properties_list[i*numprops + 3] = 1; // one of the two always wins
                properties_list[i*numprops + 4] = 1; // one of the two always loses
            }
            else { // n == 1
                properties_list[i*numprops + 1] = 1; // only one species, only one hierarchy
                properties_list[i*numprops + 3] = 1; // it doesn't have anyone to beat, but define as 1 
                properties_list[i*numprops + 4] = 0; // it doesn't have anyone to lose to, but define as 0
            }
            properties_list[i*numprops + 5] = 0;
            properties_list[i*numprops + 6] = 0;
        }
    }
    
    // output properties_list to Matlab
    ii = 0;
    for (int j = 0; j < numprops; j++) {
        for (int i = 0; i < A000568[n]; i++) {
            properties_lhs[ii] = properties_list[i*numprops + j];
            ii++;
        }
    }
            
    // delete pointers
    delete[] tournlist; tournlist = NULL;
    delete[] tourn; tourn = NULL;
    delete[] properties_list; properties_list = NULL;
    
}

void slaters_ij(int* ij, const bool* c, const int n)
// finds the distance between a tournament and the nearest hierarchy
{    

    // create a hierarchy h (i.e., all losses in upper right)
    bool* h = new bool[n*(n - 1)/2];
    for (int i = 0; i < n*(n - 1)/2; i++) {h[i] = 0;}
    
    // turn c and h into full n*n competitive outcomes matrices
    bool* C = new bool[n*n];
    bool* H = new bool[n*n];
    int ii = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                C[i*n + j] = 0;
                H[i*n + j] = 0;
            }
            else {
                C[i*n + j] = (c[ii] != 0);
                C[j*n + i] = ((1 - c[ii]) != 0);
                H[i*n + j] = (h[ii] != 0);
                H[j*n + i] = ((1 - h[ii]) != 0);
                ii++;
            }
        }
    }
    
    // determine number of reversals to convert C into every perumutation of H
    // the minimum is the distance to nearest hierarchy
    int* ord = new int[n]; // order of the rows and columns in permuted H
    for (int i = 0; i < n; i++) {ord[i] = i;} // original order of rows and columns
    int cand_s; // number of reversals to turn C to current permutation hierarchy ('candidate s')
    bool* Htemp = new bool[n*n]; // temp matrix for swapping rows
    bool* Htest = new bool[n*n]; // permuted H matrix 
    int x; // row or column for H permutation
    int s; // distance (num reversals) to *nearest* hierarchy
    int* distcounts = new int[n*(n - 1)/2]; // counts of distances to hierarchies
    for (int i = 0; i < n*(n - 1)/2; i++) {distcounts[i] = 0;} // initialize distcounts
    for (long long int k = 0; k < fact(n); k++) { // have to check every permutation
        
        // re-order H according to ord to get Htest
        for (int i = 0; i < n; i++) {
            x = ord[i]; // re-order rows
            for (int j = 0; j < n; j++) {
                Htemp[i*n + j] = H[x*n + j];
            }
        }
        for (int j = 0; j < n; j++) {
            x = ord[j]; // re-order columns
            for (int i = 0; i < n; i++) {
                Htest[i*n + j] = Htemp[i*n + x];
            }
        }
        
        // determine the distance (cand_s or 'candidate s') between C and Htest (only need to check top right)
        cand_s = 0; // reset sum to 0
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (C[i*n + j] != Htest[i*n + j]) {cand_s++;}
            }
        }
        
        // check to see whether we are closer than a previous permutation of the hierarchy
        if (k == 0) { // in the first round, we have nothing to compare with
            s = cand_s;
        }
        else {
            if (cand_s < s) {s = cand_s;}
        }
        
        // update distcounts
        distcounts[cand_s]++;
        
        // go on to the next permutation of H
        nextperm(ord, n);
        
    }
    
    // delete pointers
    delete[] h; h = NULL;
    delete[] C; C = NULL;
    delete[] H; H = NULL;
    delete[] ord; ord = NULL;
    delete[] Htemp; Htemp = NULL;
    delete[] Htest; Htest = NULL;
    
    // update ij
    ij[0] = s; // s is equivalent to Slater's i
    ij[1] = distcounts[s]; // distcounts[s] is equivalent to Slater's j 
    
}

void bezembinder(int* bez, const bool* c, const int n) 
// finds Bezembinders rho and delta_prime
{

    // convert c into full matrix form
    bool* C = new bool[n*n];
    int k = 0;
    for (int row = 0; row < n; row++) {
        for (int col = row; col < n; col++) {
            if (row == col) {
                C[row*n + col] = 0;
            }
            else {
                C[row*n + col] = (c[k] != 0);
                C[col*n + row] = ((1 - c[k]) != 0);
                k++;
            }
        }
    }
    
    // get score sequence
    int* a = new int[n];
    bool* a_used = new bool[n];
    for (int row = 0; row < n; row++) {
        a[row] = 0;
        for (int col = 0; col < n; col++) {
            a[row] += C[row*n + col];
        }
        a_used[row] = 0;
    }
    
    // sort score sequence
    int* a_sort = new int[n];
    int minscore;
    int x;
    for (int ii = 0; ii < n; ii++) {
        minscore = n - 1;
        for (int jj = 0; jj < n; jj++) {
            if ((a[jj] <= minscore) && (a_used[jj] == 0)) {
                minscore = a[jj];
                x = jj;
            }
        }
        a_sort[ii] = minscore;
        a_used[x] = 1;
    }
    
    // cumulative sorted score sequence
    int* b = new int[n];
    b[0] = a_sort[0];
    for (int ii = 1; ii < n; ii++) {
        b[ii] = b[ii - 1] + a_sort[ii];
    }
    
    // i-choose-2 sequence
    int* iC2 = new int[n];
    int i;
    for (int ii = 0; ii < n; ii++) {
        i = ii + 1;
        iC2[ii] = i*(i - 1)/2;
    }
    
    // find nu
    int nu = 0;
    for (int ii = 0; ii < n; ii++) {
        if (b[ii] == iC2[ii]) {nu++;}
    }
    
    // find size vector eta and content vector theta
    int* eta = new int[nu];
    int* theta = new int[nu];
    int i_last; // index the last time b[ii] == iC2[ii]
    k = 0;
    for (int ii = 0; ii < n; ii++) {
        i = ii + 1;
        if (b[ii] == iC2[ii]) {
            if (k == 0) {
                eta[k] = i;
            }
            else {
                eta[k] = i - i_last;
            }
            theta[k] = eta[k]*(eta[k] - 1)/2;
            i_last = i;
            k++;
        }
    }
    
    // find rho
    int rho = nu - n;
    for (int ii = 0; ii < nu; ii++) {rho += theta[ii];}
    
    // find delta_prime (doesn't divide by nC2)
    int delta_prime = (rho - nu + n);
    
    // delete pointers
    delete[] C; C = NULL;
    delete[] a; a = NULL;
    delete[] a_sort; a_sort = NULL;
    delete[] a_used; a_used = NULL;
    delete[] b; b = NULL;
    delete[] iC2; iC2 = NULL;
    delete[] eta; eta = NULL;
    delete[] theta; theta = NULL;
    
    // update bez
    bez[0] = rho;
    bez[1] = delta_prime;
}

int intrans_triads(const bool* c, const int n) 
// finds the number of intransitive triads in tournament c
{
    // turn c into full n*n competitive outcomes matrix
    bool* C = new bool[n*n];
    int ii = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                C[i*n + j] = 0;
            }
            else {
                C[i*n + j] = (c[ii] != 0);
                C[j*n + i] = ((1 - c[ii]) != 0);
                ii++;
            }
        }
    }
    
    bool t[3]; // test triad (upper right)
    int numintrans = 0; // number of intransitive triads
    for (int i = 0; i < n - 2; i++) {
        for (int j = i + 1; j < n - 1; j++) {
            for (int k = j + 1; k < n; k++) {
                t[0] = C[i*n + j]; 
                t[1] = C[i*n + k];
                t[2] = C[j*n + k];
                if (((t[0] == 0) && (t[1] == 1) && (t[2] == 0)) || ((t[0] == 1) && (t[1] == 0) && (t[2] == 1))) {
                    numintrans++;
                }
            }
        }
    }
    
    // delete pointers
    delete[] C; C = NULL;
    
    // return
    return numintrans;
    
}

bool unbeatable(const bool* c, const int n) 
// determines whether there is an unbeatable species
{
    // turn c into full n*n competitive outcomes matrix
    bool* C = new bool[n*n];
    int ii = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                C[i*n + j] = 0;
            }
            else {
                C[i*n + j] = (c[ii] != 0);
                C[j*n + i] = ((1 - c[ii]) != 0);
                ii++;
            }
        }
    }
    
    int maxwins = 0;
    int wins;
    for (int i = 0; i < n; i++) {
        wins = 0;
        for (int j = 0; j < n; j++) {
            wins += C[i*n + j];
        }
        if (wins > maxwins) {maxwins = wins;}
    }
    
    delete[] C; C = NULL;
    
    if (maxwins == n - 1) {
        return 1;
    }
    else {
        return 0;
    }
    
}

bool alwaysbeatable(const bool* c, const int n)
// determines whether there is a species beaten by every species
{
    
    // turn c into full n*n competitive outcomes matrix
    bool* C = new bool[n*n];
    int ii = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                C[i*n + j] = 0;
            }
            else {
                C[i*n + j] = (c[ii] != 0);
                C[j*n + i] = ((1 - c[ii]) != 0);
                ii++;
            }
        }
    }
    
    int minwins = n - 1;
    int wins;
    for (int i = 0; i < n; i++) {
        wins = 0;
        for (int j = 0; j < n; j++) {
            wins += C[i*n + j];
        }
        if (wins < minwins) {minwins = wins;}
    }
    
    delete[] C; C = NULL;
    
    if (minwins == 0) {
        return 1;
    }
    else {
        return 0;
    }
    
}

long long int fact(const int n)
// gives the factorial (works for n <= 20)
{
    long long int a = 1;
    for (int i = 1; i <= n; i++) {a *= i;}
    return a;
}

void nextperm(int* currperm, const int n)
// gives the next permutation of a sequence of integers
{    

    // find the tail
    bool tailfound = 0; // switch to see if tail is found
    int i = n - 1; // start at the back
    int tailstart = 0; // default if whole thing is tail
    while ((tailfound == 0) && (i >= 1)) {
        if (currperm[i - 1] < currperm[i]) {
            tailstart = i;
            tailfound = 1;
        }        
        i--;
    }
    
    // make a swap between the last element of the head, and the farthest-back element of the tail
    // that is larger than the last element of the head
    bool swapdone = 0; // determines whether the swap has been done or not (default is zero)
    i = n - 1; // start at the back of the tail
    int temp; // temp integer for swapping
    if (tailfound == 1) { // don't bother doing if whole thing is tail; just reverse the whole thing
        while ((i >= tailstart) && (swapdone == 0)) { // go as far as the start of the tail
            if (currperm[i] > currperm[tailstart - 1]) {
                temp = currperm[i];
                currperm[i] = currperm[tailstart - 1];
                currperm[tailstart - 1] = temp;
                swapdone = 1;
            }
            i--;
        }        
    }
    
    // reverse the tail
    int* temptail = new int[n - tailstart]; // temporary tail for swapping
    for (int j = 0; j < n - tailstart; j++) {temptail[j] = currperm[n - j - 1];}
    for (int j = 0; j < n - tailstart; j++) {currperm[tailstart + j] = temptail[j];}
    
    // delete pointers
    delete[] temptail; temptail = NULL;   
    
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *properties_lhs, *tournlist_rhs, *n_rhs;
    
    const int A000568_in[11] = {1, 1, 1, 2, 4, 12, 56, 456, 6880, 191536, 9733056};
    const int n_in = static_cast<int>(*mxGetPr(prhs[1]));
    
    // error check
    // code only works for up to n = 10. Beyond that, the number of entries in tournlist exceeds 2^31 - 1, the maximum size of data type 'int'
    if (n_in > 10) {
        mexErrMsgTxt("Function tournlistproperties does not apply to n > 10 (the number of entries in tournlist exceeds 2^31 - 1, the maximum size of data type 'int')");
    }

    // prepare properties_lhs
    int mrows = A000568_in[n_in]; 
    int ncols = 7;  
    plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, (mwSize)ncols, mxREAL);
    
    properties_lhs = mxGetPr(plhs[0]);
    
    tournlist_rhs = mxGetPr(prhs[0]);
    n_rhs = mxGetPr(prhs[1]);
    
    main(properties_lhs, tournlist_rhs, n_rhs);
    
}