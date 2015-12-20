// Author: Shrikant Kshirsagar
// Atom class - store atom attributes
#ifndef MYATOM_H
#define MYATOM_H
#include "MyCrystalStructure.h"

class MyAtom
{
    public:
        int m_ID; // atom ID
        double m_position[3]; // atom position in Cartesian coordinates
        int m_unitCell[3]; // triplet identifying the unit cell of the atom
        int m_basis; // identify the basis no. of atom 
        void m_setAtom(int a_ID, double a_position[3], int a_unitCell[3], int a_basis);
        // set atom data
};
#endif
