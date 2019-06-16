/*---------------------------------------------------------------------------*\
|                plicVofSolver | Copyright (C) 2019 Dezhi Dai                 |
-------------------------------------------------------------------------------
License
    This file is part of plicVofSolver which is an extension to OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::plicInterfaceField

Description
    A field of plicInterface.

    Reference:
        \verbatim
            Dai, Dezhi and Tong, Albert Y. (2019).
            Analytical interface reconstruction algorithms in the PLIC‐VOF
            method for 3D polyhedral unstructured meshes
            International Journal for Numerical Methods in Fluids
            doi 10.1002/fld.4750
            url https://doi.org/10.1002/fld.4750
        \endverbatim

SourceFiles
    plicInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef plicInterfaceField_H
#define plicInterfaceField_H

#include "plicInterface.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plicInterface Declaration
\*---------------------------------------------------------------------------*/
class plicInterfaceField
{
private:

    // Private data

        //- Number of elements in field
        label size_;

        //- Vector of plicInterfaces
        plicInterface* plicInterfaces_;


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct from VOF (alpha) field
        plicInterfaceField(volScalarField& alpha1);


    // Member operators

        //- Return element of plicInterfaceField
        plicInterface& operator[](const label i);

        //- Return element of constant plicInterfaceField
        const plicInterface& operator[](const label i) const;

    // Member functions

        //- Return element of plicInterfaceField
        plicInterface& interface(const label i);

        //- Return element of constant plicInterfaceField
        const plicInterface& interface(const label i) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // End namespace Foam

// ************************************************************************* //