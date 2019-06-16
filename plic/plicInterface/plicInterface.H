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
    Foam::plicInterface

Description
    plicInterface is a 3D plan, which is constructured from its unit normal vector $vec{n}$ and signed distance $D_0$.

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
    plicInterface.C

\*---------------------------------------------------------------------------*/

#ifndef plicInterface_H
#define plicInterface_H

#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plicInterface Declaration
\*---------------------------------------------------------------------------*/
class plicInterface
{
private:

    // Private data

        //- unit normal vector $\vec{n}$
        vector n_;

        //- interface center
        vector X_;

        //- signed distance $D_0$
        scalar D_;


public:

    // Static data members

        static const char* const typeName;


        // Position of a point relative to the Interface
        enum pointSide
        {
            PNORMAL,    //point is on the Interface (interface)
            PABOVE,     //point is above the Interface (inside gas phase)
            PBELOW      //point is below the Interface (inside fluid phase)
        };


    // Constructors

        //- Null constructor
        plicInterface();

        //- Construct from normal vector through the origin
        plicInterface(const vector& normalVector);

        //- Construct from normal vector and signed distance
        plicInterface(const vector& normalVector, const scalar signedDistance);

        //- Construct from normal vector and point of an interface
        plicInterface(const vector& normalVector, const point& basePoint);


    // Member Functions

        //- Return unit normal vector
        vector n() const;

        //- Return unit normal vector
        vector& n();

        //- Return interface center
        point X() const;

        //- Return interface center
        point& X();

        //- Return signed distance
        scalar D() const;

        //- Return signed distance
        scalar& D();

        //- Return the signed distance form one point to the interface
        scalar signedDistance(const point&) const;

        //- Return the side of the interface that the point is on
        pointSide sideOfPoint(const point&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //