// Author: Shrikant Kshirsagar
// Generate a supercell of a given size of a given crystal structure
// and compute the Force Constants between the central atom and all
// other atoms in the supercell
#include <iostream>
#include <iomanip>
#include <cmath>
#include "MyAtom.cpp"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define DELTA(a,b) (((a)==(b))?(1):(0)) // Kroenecker delta
using namespace std;

//---------------------- User Data ----------------------------------------------------
int g_supercellSize[3]={6,6,6}; // define size of supercell in terms of 
                                // no. of conventional unit cells along each dimension
                                // each of the 3 nos. must be >=1
double g_latticeConst = 5.2672; // lattice constant (a) in Angstrom
double g_A = 99275, g_B = 64.264; // LJ factors for Ar (Energy in eV, distance in A)

int main()
{
    MyCrystalStructure Ar_FCC; // We'll build FCC Ar crystal lattice
    // define primitive lattice vectors for FCC lattice
    double a1[3] = {0,0.5,0.5};
    double a2[3] = {0.5,0,0.5};
    double a3[3] = {0.5,0.5,0};
    Ar_FCC.m_noBasis = 1; // Ar has one basis atom per lattice site
    double basisPos[3] = {0.0,0.0,0.0}; // the atom is taken to be at it's lattice site
    // for additional basis atoms add more loops
    // REMEMBER to the basis position by a
//---------------------- User Data ----------------------------------------------------
    equateTriplet<double*>(Ar_FCC.m_a[0],a1); // set primitive lattice vectors
    equateTriplet<double*>(Ar_FCC.m_a[1],a2);
    equateTriplet<double*>(Ar_FCC.m_a[2],a3);
    equateTriplet<double*>(Ar_FCC.m_posBasis[0],basisPos); // set basis atom position
    int Nx = g_supercellSize[0], Ny = g_supercellSize[1], Nz = g_supercellSize[2];
    int noAtoms = 4*Nx*Ny*Nz + 2*(Nx*Ny + Ny*Nz + Nz*Nx) + Nx+Ny+Nz + 1;
    // total no. of atoms in the supercell of the given size
    // Set no. of conventional unit cells along each direction
    int plusXCells = MAX(Nx/2,1); // Allocate half to +X direction, min. 1
                               // NB: Nx/2 performs integer division
    int minusXCells = MAX(Nx - plusXCells,0); // Allocate remaining cells to -X
    int plusYCells = MAX(Ny/2,1);
    int minusYCells = MAX(Ny - plusYCells,0);
    int plusZCells = MAX(Nz/2,1);
    int minusZCells = MAX(Nz - plusZCells,0);
    double xmax = plusXCells * g_latticeConst; // max x-coordinate in supercell
    double xmin = -minusXCells * g_latticeConst;
    double ymax = plusYCells * g_latticeConst; 
    double ymin = -minusYCells * g_latticeConst;
    double zmax = plusZCells * g_latticeConst; 
    double zmin = -minusZCells * g_latticeConst;
    // -------------- Output --------------------------------------------------------------------------------------------------
    cout<<setw(120)<<setfill('-')<<"\n";
    cout<<"Input parameters:\n";
    cout<<"Supercell size: "<<g_supercellSize[0]<<"-by-"<<g_supercellSize[1]<<"-by-"<<g_supercellSize[2]<<" (conventional unit cells) which corresponds to (distances in A):\n";
    cout<<xmin<<" <= x <= "<<xmax<<",  "<<ymin<<" <= y <=  "<<ymax<<",  "<<zmin<<" <= z <= "<<zmax<<"\n";    
    // -------------- Output --------------------------------------------------------------------------------------------------
    int atomID  = 0; // stores current atom ID
    double coord[3]; // stores current x,y,z coordinates
    int unitCellID[3];// current unit cell ID
    MyAtom* atoms = new  MyAtom[noAtoms]; // array to hold all atoms in the supercell
    int nlim = (Nx + Ny + Nz)/2; // limits on n1, n2, n3 see 4/7/14
    // Create the supercell and ID all atoms in it
    for (int n1dir = 0; n1dir < 2; n1dir++)
    {// loop over +n1 then -n1 values 
        int n1inc = 1 - 2*n1dir;
       for (int n2dir = 0; n2dir < 2; n2dir++)
       { 
           int n2inc = 1 - 2*n2dir;
           for (int n3dir = 0; n3dir < 2; n3dir++)
           { 
               int n3inc = 1 - 2*n3dir;
               for (int n1 = -n1dir; abs(n1) <= nlim; n1 += n1inc)
               {
                   for (int n2 = -n2dir; abs(n2) <= nlim; n2 += n2inc)
                   {
                       for (int n3 = -n3dir; abs(n3) <= nlim; n3 += n3inc)
                       {
                           if ( (abs(n2 + n3) <= Nx) && (abs(n3 + n1) <= Ny) && (abs(n1 + n2) <= Nz) )
                           {
                               unitCellID[0] = n1;
                               unitCellID[1] = n2;
                               unitCellID[2] = n3;
                               coord[0] = (n2 + n3)*g_latticeConst/2;
                               coord[1] = (n3 + n1)*g_latticeConst/2;
                               coord[2] = (n1 + n2)*g_latticeConst/2;
                               atoms[atomID].m_setAtom(atomID,coord,unitCellID,0); // for more basis atoms need more statements for each one of them
                               ++atomID; // increment atom ID
                           }
                       }
                   }
               }
           }
       }
    }
    // -------------- Output --------------------------------------------------------------------------------------------------
    cout<<setw(120)<<setfill('-')<<"\n";
    cout<<"Atom IDs, positions, unit cell nos. and basis atom :\n"<<setfill(' ');
    cout<<setw(5)<<"ID"<<setw(12)<<"x"<<setw(12)<<"y"<<setw(12)<<"z"<<setw(4)<<"n1"<<setw(4)<<"n2"<<setw(4)<<"n3"<<setw(4)<<"b\n";
    cout<<setw(120)<<setfill('-')<<"\n"<<setfill(' ');
    for (int i=0; i<noAtoms; ++i)
    {
        cout<<setw(5)<<i;
        for (int j=0; j<3; ++j)
            cout<<setprecision(8)<<setw(12)<<atoms[i].m_position[j];
        for (int j=0; j<3; ++j)
            cout<<setw(4)<<atoms[i].m_unitCell[j];
        cout<<setw(4)<<atoms[i].m_basis<<"\n";;
    }
    // -------------- Output --------------------------------------------------------------------------------------------------    
    double Rij[3],Pi[noAtoms][3],Phi[noAtoms][3][3],Psi[noAtoms][3][3][3]; // stores vector joining atom i to atom j and FCs (for cubic FCs store only Psi_iij)
    MyAtom atomj; // store data for each atom in the supercell
    //double C1sh,C2sh[3][3],C2sc[3],C3sc[3][3][3]; // accumalators to compute self-interaction terms
    // initialize the accumalators
    //C1sh = 0.0;
    /*for(int a=0; a<3; ++a)
    {
        //C2sc[a] = 0.0;
        Pi[a] = 0.0;
        for(int b=0; b<3; ++b)
        {
            //C2sh[a][b] = 0.0;
            for(int c=0; c<3; ++c)
                C3sc[a][b][c] = 0.0;
        }
    }*/
    for (int j=1; j<noAtoms; ++j)
    { // loop over all neighbor atoms of central atom (atom 0)
        atomj = atoms[j];
        equateTriplet<double*>(Rij,atomj.m_position); // because atom i is at origin, it's coordinates are (0,0,0)
        double RijSq = 0.0; // calculate sum of distance squared
        for (int k=0; k<3; ++k)
            RijSq += Rij[k] * Rij[k]; // Pythagorean distance formula
        double C1ij = -6*g_A/pow(RijSq,7) + 3*g_B/pow(RijSq,4);
        double C2ij = 84*g_A/pow(RijSq,8) - 24*g_B/pow(RijSq,5); // compute C_2(i,j) and C_3(i,j)
        double C3ij = -448*g_A/pow(RijSq,9) + 80*g_B/pow(RijSq,6);
        //C1sh += C1ij; // for computing self-interaction term
        for (int a=0; a<3; ++a)
        { // loop over alpha Cartesian directions
            Pi[j][a] = 2*C1ij*Rij[a];
            //C2sc[a] += C2ij*Rij[a]; // for computing self-interaction term
            for (int b=0; b<3 ;++b)
            { // loop over beta Cartesian directions
                Phi[j][a][b] = -2*(C1ij*DELTA(a,b) + C2ij*Rij[a]*Rij[b]); // Phi_{i,i}^{\alpha,\beta}(i,j)
                //C2sh[a][b] += C2ij*Rij[a]*Rij[b]; // for computing self-interaction term
                for (int c=0; c<3; ++c)  // loop over gamma Cartesian directions
                {
                    //Psi[j][a][b][c] = 2*(C2ij*(Rij[a]*DELTA(b,c) + Rij[b]*DELTA(a,c) + Rij[c]*DELTA(a,b)) + 3*C3ij*Rij[a]*Rij[b]*Rij[c]); // Psi_iij
                    Psi[j][a][b][c] = 0;
                    //C3sc[a][b][c] += C3ij*Rij[a]*Rij[b]*Rij[c]; // for computing self-interaction term
                }
            }
        }
    }
    // Computing the self interaction harmonic and cubic FCs using the Acoustic Sum Rule
    for (int a=0; a<3; ++a)
    {
        Pi[0][a]=0;
        for (int b=0; b<3; ++b)
        {
            Phi[0][a][b] = 0; //2*(C1sh*DELTA(a,b) + C2sh[a][b]);
            for (int c=0; c<3; ++c)
                {
                    Psi[0][a][b][c] = 0;//-6*(C2sc[b]*DELTA(a,c) + C3sc[a][b][c]);
                    for (int j=1; j<noAtoms; ++j)
                    {
                        Pi[0][a] -= Pi[j][a];
                        Phi[0][a][b] -= Phi[j][a][b];
                        Psi[0][a][b][c] -= Psi[j][a][b][c];
                    }
                }

        }
    }
    // -------------- Output --------------------------------------------------------------------------------------------------    
    cout<<setw(120)<<setfill('-')<<"\n";
    cout<<"Table of interatomic force constants between atom 0 and other atoms\n";
    cout<<"(Here, atom i = 0 & neighbor atoms are denoted by k & for cubic terms it is understood that j=i; a = alpha, b = beta, c = gamma)\n";
    cout<<"NB: Energy in eV, distances in A\n";
    cout<<setw(120)<<setfill('-')<<"\n"<<setfill(' ');
    cout<<setw(5)<<"jID"<<setw(5)<<"a"<<setw(5)<<"b"<<setw(5)<<"c"<<setw(15)<<"Pi"<<setw(15)<<"Phi"<<setw(15)<<"Psi\n"; // Only Psi_iij as that satisfies ASR in the regular sense
    for (int j=0; j<noAtoms; ++j)
    {
        for (int a=0; a<3; ++a)
        {
            for (int b=0; b<3; ++b)
            {
                for (int c=0; c<3; ++c)
                    cout<<setw(5)<<j<<setw(5)<<(a+1)<<setw(5)<<(b+1)<<setw(5)<<(c+1)<<setprecision(5)<<setw(15)<<Pi[j][a]<<setw(15)<<Phi[j][a][b]<<setw(15)<<Psi[j][a][b][c]<<endl;
            }
        }
    }
    // -------------- Output --------------------------------------------------------------------------------------------------    
    /*
    double cterm,hterm,RjiSq,C2ij,C3ij,Phiii[3][3],Phiij[3][3],Phijj[3][3],Psiiii[3][3][3],Psiiij[3][3][3],Psiijj[3][3][3],Psijjj[3][3][3];
    // for storage of FCs and intermediate values for their calculation
    MyAtom atomj; // store data for each atom in the supercell
    for (int j=1; j<noAtoms; ++j)
    { // loop over all neighbor atoms of central atom (atom 0)
        atomj = atoms[j];
        equateTriplet<double*>(Rji,atomj.m_position); // because atom i is at origin, it's coordinates are (0,0,0)
        RjiSq = 0.0; // calculate sum of distance squared
        for (int k=0; k<3; ++k)
            RjiSq += Rji[k] * Rji[k]; // Pythagorean distance formula
        C2ij = 84*g_A/pow(RjiSq,8) - 24*g_B/pow(RjiSq,5); // compute C_2(i,j) and C_3(i,j)
        C3ij = -448*g_A/pow(RjiSq,9) + 80*g_B/pow(RjiSq,6);
        for (int a=0; a<3; ++a)
        { // loop over alpha Cartesian directions
            for (int b=0; b<3 ;++b)
            { // loop over beta Cartesian directions
                Phiii[a][b] = C2ij*Rji[a]*Rji[b]; // Phi_{i,i}^{\alpha,\beta}(i,j)
                Phiij[a][b] = -2*Phiii[a][b]; // Phi_{i,j}^{\alpha,\beta}(i,j)
                Phijj[a][b] = Phiii[a][b]; // Phi_{j,j}^{\alpha,\beta}(i,j)
                for (int c=0; c<3; ++c)
                { // loop over gamma Cartesian directions
                    cout<<setw(5)<<j<<setw(3)<<(a+1)<<setw(3)<<(b+1)<<setw(3)<<(c+1);
                    cterm = C3ij*Rji[a]*Rji[b]*Rji[c]; // first term in Psi
                    hterm = C2ij*Rji[b]*DELTA(c,a) ; // second term in Psi
                    Psiiii[a][b][c] = -hterm - cterm ; // Psi_{i,i,i}^{\alpha,\beta,gamma}(i,j)
                    Psijjj[a][b][c] = -Psiiii[a][b][c]; // Psi_{j,j,j}^{\alpha,\beta,gamma}(i,j)
                    Psiiij[a][b][c] = 2*hterm + C2ij*Rji[c]*DELTA(a,b) + 3*cterm; // Psi_{i,i,j}^{\alpha,\beta,gamma}(i,j)
                    Psiijj[a][b][c] = -2*hterm - C2ij*Rji[a]*DELTA(b,c) - 3*cterm; // Psi_{i,j,j}^{\alpha,\beta,gamma}(i,j)
                    Psiiii[a][b][c]=Psiiij[a][b][c]=Psiijj[a][b][c]=Psijjj[a][b][c]=0; // for getting pure harmonic
                cout<<setprecision(5)<<setw(14)<<Phiii[a][b]<<setw(14)<<Phiij[a][b]<<setw(14)<<Phijj[a][b]<<setw(14)<<Psiiii[a][b][c]<<setw(14)<<Psiiij[a][b][c]<<setw(14)<<Psiijj[a][b][c]<<setw(14)<<Psijjj[a][b][c]<<endl;
                }
            }
        }
    }
    */
    delete[] atoms;
    return 0;
}

