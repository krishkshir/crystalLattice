// Author: Shrikant Kshirsagar
// Defines crystal structure
#ifndef MYCRYSTALSTRUCTURE_H
#define MYCRYSTALSTRUCTURE_H

class MyCrystalStructure
{
    public:
        double m_a[3][3]; // primitive lattice vectors for crystal structure
                         // each row gives the components of one vector
        int m_noBasis; // no. of basis atoms in given crystal structure
        double m_posBasis[1][3];  // store positions of basis atoms in a unit
                                 // cell; first index must equal m_noBasis
};
#endif

