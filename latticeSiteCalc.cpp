#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#define LATTICECONST 5.2672
#define DELTAn 1.5 // Limits on max \Delta n_i
using namespace std;

double distCalc(const double a_currPos[3],const int a_unitCellID[3],double a_eqbPos[3])
{
  int n1 = a_unitCellID[0], n2 = a_unitCellID[1], n3 = a_unitCellID[2];
  a_eqbPos[0] = 0.5*LATTICECONST*(n2 + n3);
  a_eqbPos[1] = 0.5*LATTICECONST*(n1 + n3);
  a_eqbPos[2] = 0.5*LATTICECONST*(n1 + n2);
  double dist = 0.0;
  for (int i = 0; i < 3; ++i)
      dist += ( a_currPos[i] - a_eqbPos[i] ) * ( a_currPos[i] - a_eqbPos[i] ) ;
  dist = sqrt(dist);
  cout<<"Distance from lattice site with unit cell ID: "<<n1<<"   "<<n2<<"   "<<n3<<" is: "<<dist<<"\n";
  return dist;
}

int possUnitCells(const double a_m,int a_possUnitCellID[])
{
    int i;
    int s = ceil(a_m - DELTAn);
    int L = floor(a_m + DELTAn);
    for (i = 0; i <= (L-s); ++i)
        a_possUnitCellID[i] = s+i;
    return i;
}

int main()
{
  double checkSqDist = 0.5*LATTICECONST/sqrt(2); // nearest lattice site max. distance squared
  cout<<"LATTICE CONSTANT = "<<LATTICECONST<<", therefore ";
  cout<<"checkDist= "<<checkSqDist<<"\n";
  double r[3],R[3]; // store coordinates of current and equilibrium lattice positions
  int unitCellID[3]; // store unit cell ID integers
  cout<<"Enter current positions: ";
  cin>>r[0]>>r[1]>>r[2];
  double m1, m2, m3; // multiples of lattice primitive vectors contained in currPos
  m1 = ((-r[0]+r[1]+r[2])/LATTICECONST);
  m2 = ((r[0]-r[1]+r[2])/LATTICECONST);
  m3 = ((r[0]+r[1]-r[2])/LATTICECONST);
  cout<<"Current position ("<<r[0]<<", "<<r[1]<<", "<<r[2]<<") corresponds to the following multiple of lattice basis vectors: ";
  cout<<m1<<"   "<<m2<<"   "<<m3<<"\n";
  unitCellID[0] = round(m1); unitCellID[1] = round(m2); unitCellID[2] = round(m3);
  double sqDist = distCalc(r,unitCellID,R); 
  if ( sqDist > checkSqDist )
  {
      cout<<"Rounding unit cell IDs didn't work!\n";
      int maxnoPossni = floor(2 * DELTAn + 1); // max.no Possni
      int noPossn1,noPossn2,noPossn3,possn1[maxnoPossni],possn2[maxnoPossni],possn3[maxnoPossni],neighCellID[3];
      double R_t[3];
      noPossn1 = possUnitCells(m1,possn1);
      noPossn2 = possUnitCells(m2,possn2);
      noPossn3 = possUnitCells(m3,possn3);
      bool gotit = 0;
      for (int i = 0; (i < noPossn1) && (!gotit); ++i)
      {
          for (int j = 0; (j < noPossn2) && (!gotit); ++j)
          {
              for (int k = 0; (k < noPossn3) && (!gotit); ++k)
              {
                  neighCellID[0] = possn1[i]; neighCellID[1] = possn2[j]; neighCellID[2] = possn3[k];
                  double sqDist_t = distCalc(r,neighCellID,R_t);
                  if ( sqDist_t < sqDist )
                  { // try to find the minimum distance between current position and lattice site
                      sqDist = sqDist_t;
                      unitCellID[0] = neighCellID[0]; unitCellID[1] = neighCellID[1]; unitCellID[2] = neighCellID[2];
                      R[0] = R_t[0]; R[1] = R_t[1]; R[2] = R_t[2];
                  }
                  //if ( sqDist < checkSqDist )
                      //gotit = 1;
              }
          }
      }
  }
  cout<<"Nearest unit cell ID: "<<"   "<<unitCellID[0]<<"   "<<unitCellID[1]<<"   "<<unitCellID[2];
  cout<<" corresponding to equilibrium lattice position ("<<R[0]<<", "<<R[1]<<", "<<R[2]<<")\n"; 
  cout<<"Displacement from lattice site = "<<50*sqDist/checkSqDist<<"% of nearest neighbor distance.\n"; // nearest neighbor dist = a/sqrt(2) = 2*checkSqDist
  return 0;
}
