// Author: Shrikant Kshirsagar
#include "MyAtom.h"
#include "MyCrystalStructure.h"

template <typename T>
void equateTriplet(T a_lhs,T a_rhs)
{ // generic function template to equate triplets of any type
    for (int i=0; i<3; ++i)
        a_lhs[i] = a_rhs[i];
}

void MyAtom::m_setAtom(int a_ID, double a_position[3], int a_unitCell[3], int a_basis)
{
    m_ID = a_ID;
    equateTriplet<double*>(m_position,a_position);
    equateTriplet<int*>(m_unitCell,a_unitCell);
    m_basis = a_basis;
}

/*double calculateDistanceSquared(MyAtom a_atomi, MyAtom a_atomj)
{// function to compute square of distance between two atoms
    double result = 0.0,Ri[3],Rj[3]; // Ri & Rj store the lattice positions of the 2 atoms
    equateTriplet<double*>(Ri,a_atomi.m_position); // store atom i's position
    equateTriplet<double*>(Rj,a_atomj.m_position); //store atom j's position
    for (int k=0; k<3; ++k)
        result += (Rj[k] - Ri[k])*(Rj[k] - Ri[k]); // Pythagorean distance formula
    return result;
}*/
