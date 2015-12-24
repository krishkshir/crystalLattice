// Program to solve the Closest Lattice Vector Problem (CVP) for a general crystal lattice (with basis) in 3 dimensions
/* The basic idea is that given a lattice basis A (columns are basis vectors) and a arbitrary point in 3D r, to solve the linear system Am = r to get fractional (real) 
   coordinates m and round it to get integer coordinates (of the unit cell) n. For small displacements about lattice sites, this method will work most of the times
   and this can be checked by ensuring the 1-norm distance between m and n < 1/2.
   If this is not the case, then construct parallelopipeds around n of increasing size (credits to Prof. Christos Papadimitrou for suggesting this technique) until the
   sphere with radius equal to the current lowest distance (in Euclidean norm in 3D) between r and the current closest lattice vector is completely within the parallelopiped.
   For a lattice with basis, repeat this procedure with r-b where b is the basis atom position for each basis atom.
 */
#include <iostream>
#include <cmath>
// END include libraries
using namespace std;
// BEGIN GLOBAL PARAMETERS
const int g_dim = 3; // no. of dimensions
const double g_a = 1.0; // lattice constant
const double g_A[g_dim][g_dim] = {0, 0.5, 0.5,  0.5, 0, 0.5,  0.5, 0.5, 0}; // columns denote basis vectors of lattice
const double g_Ainv[g_dim][g_dim] = {-1, 1, 1,  1, -1, 1,  1, 1, -1}; // inverse of A
const int g_nbasis = 1; // no of basis atoms
const double g_b[g_nbasis][g_dim] = {0, 0, 0}; // basis atom positions
const double g_initializermax = 100*g_a; // used to initialize the shortest distance calculated so far in CVP problem 
double distCalc(const double a_currPos[g_dim],const int a_unitCellID[g_dim],double a_eqbPos[g_dim])
{ // calculate distance between given arbitrary point currPos and a lattice site specified with unitCellID while also storing lattice position vector in eqbPos
    for (int i=0; i<g_dim; ++i)
        a_eqbPos[i] = g_a*(g_A[i][0]*a_unitCellID[0]+g_A[i][1]*a_unitCellID[1]+g_A[i][2]*a_unitCellID[2]); // R = A*n
    double dist = 0.0; // store distance
    for (int i = 0; i < 3; ++i)
        dist += ( a_currPos[i] - a_eqbPos[i] ) * ( a_currPos[i] - a_eqbPos[i] ) ; // Euclieadn formula for squared distance
    return sqrt(dist);
}
int possUnitCells(const int a_ni[g_dim],const int a_K,int a_np[][g_dim])
{ // elaborate all the lattice vectors on this parallelopiped
    int count = 0; // running count of lattice vectors
    for (int corn=1; corn > -2; corn -=2){
        // corn=1 => vertex with greatest unit cell coordinates in this parallelopiped, corn=-1 => vertex with the smallest
        int nc[g_dim];
        for (int i=0; i < g_dim; ++i)
            nc[i] = a_ni[i]+a_K*corn; // set corner coordinates
        for (int face = abs(corn); face < abs(corn) + g_dim; ++face){
            // face = 1 => face with normal along +/- x-axis etc.
            int bsd = face % g_dim; // increment along this unit cell axis for big step
            int ssd = (face+1) % g_dim; // increment along this unit cell axis for small step
            for (int bst = 0; bst < 2*a_K+1; ++bst){ // take big step
                int ns[g_dim]; // store unit cell coordinates of starting position
                for (int i=0; i < g_dim; ++i)
                    ns[i] = nc[i]; // initialize start position to corner
                ns[bsd] -= corn*bst; // proceed along vector a[bsd] bst steps
                for (int sst = 0; sst < 2*a_K+1; ++sst){ // take small step
                    int nt[g_dim]; // store unit cell coordinates of trial position
                    for (int i=0; i < g_dim; ++i)
                        nt[i] = ns[i]; // initialize trial position to start
                    nt[ssd] -= corn*sst; // proceed along vector a[ssd] sst steps
                    for (int i=0; i < g_dim; ++i)
                        a_np[count][i] = nt[i];
                    ++count;
                }
            }
        }
    }
    return count;
}
double CVP3D(const double a_r[g_dim],int a_n[g_dim],double a_R[g_dim])
{ // finds the closest lattice vector to a_r and stores it in a_n, a_R and returns the distance between a_r and a_R
    // initial approximation is to solve linear system of equation R' = Ainv*r
    double m[g_dim];
    for (int i=0; i<g_dim; ++i)
        m[i] = (g_Ainv[i][0]*a_r[0]+g_Ainv[i][1]*a_r[1]+g_Ainv[i][2]*a_r[2])/g_a; // approx. R = Ainv*n
    for (int i=0; i<g_dim; ++i)
        a_n[i] = round(m[i]);
    double d = distCalc(a_r,a_n,a_R); // stores current min. distance from a_r to closest lattice vector
    if (abs(m[0]-a_n[0]) + abs(m[1]-a_n[1]) + abs(m[2]-a_n[2]) < 0.5) return d; // r is closer than half nearest neighbor distance from R
    // else we're forced to do the method of expanding parallelopipeds
    // cout<<"Rounding didn't give a lattice site closer than half nearest neighbor distance, d = "<<d<<"...\n";
    int K = 0; // denotes size of parallelopiped about a_n
    while (1){
        ++K;
        // cout<<"K = "<<K;
        int np[2*g_dim*(2*K+1)*(2*K+1)][g_dim]; // total possible no. of vertices
        int poss = possUnitCells(a_n,K,np); // no. of possible unit cells returned by possUnitCells
        // cout<<", np =  "<<poss; // check to make sure this is equal to 6(2K+1)^2
        double dk = g_initializermax; // stores min. distance between r and parallelopiped, initialized to a very large value
        int nk[g_dim]; // store unit cell coordinates of closest lattice vector of this parallelopiped
        double Rk[g_dim]; // store closest lattice vector position of this parallelopiped
        for(int nposs=0; nposs < poss; ++nposs){
            // compute distance of current position from this lattice site
            double Rt[g_dim];
            double dt = distCalc(a_r,np[nposs],Rt);
            if (dt < dk) { // candidate for min. over current parallelopiped
                dk = dt; 
                for (int i=0; i < g_dim; ++i){
                    nk[i] = np[nposs][i]; Rk[i] = Rt[i];
                }
            }
        }
        // cout<<", dk = "<<dk<<endl;
        if (dk < d) { // found a closer vertex to r in this parallelopiped
            d = dk;
            for (int i=0; i < g_dim; ++i){
                a_n[i] = nk[i]; a_R[i] = Rk[i];
            }
        }
        else return d; // sphere of radius equal to  distance from closest lattice vector so far lies completely within parallelopiped
    }
}
double latticeSiteCalc(const double a_r[g_dim],int a_n[g_dim+1],double a_R[g_dim])
{ // find closest crystal lattice point (incl. basis atoms) to current position
    double d = g_initializermax; // stores min. distance obtained so far
    for (int i=0; i < g_nbasis; ++i) { // solve CVP3D problem for each basis atom
        int nb[g_dim]; double Rb[g_dim]; double rb[g_dim];
        for (int j=0; j < g_dim; ++j)
            rb[j] = a_r[j] - g_b[i][j]; // move origin by basis atom position 
        double db = CVP3D(rb,nb,Rb); // and then solve CVP3D problem
        if (db < d) { // this basis atom is closest to current position so far
            d = db;
            for (int j=0; j< g_dim; ++j) { // update closest crystal lattice position
                a_n[j] = nb[j]; a_R[j] = Rb[j] + g_b[i][j]; // add back basis atom position to R to get its coordinates w.r.t. original origin
            }
            a_n[g_dim] = i; // the last coordinate gives the basis atom corresponding to this lattice site
        }
    }
    return d;
}
int main()
{
    cout<<"Solving closest vector problem for 3D crystal lattice\n";
    double r[g_dim]; int n[g_dim+1]; double R[g_dim];
    cout<<"Enter coordinates of arbitrary point: ";
    cin>>r[0]>>r[1]>>r[2];
    double d = latticeSiteCalc(r,n,R);
    cout<<"Closest lattice vector has coordinates: ("<<R[0]<<", "<<R[1]<<", "<<R[2]<<") at a distance of "<<d<<endl;
}
