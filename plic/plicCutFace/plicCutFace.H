/*---------------------------------------------------------------------------*\
|                plicVofSolver | Copyright (C) 2019 Dezhi Dai                 |
-------------------------------------------------------------------------------
                   isoAdvector | Copyright (C) 2016-2017 DHI
                 Modified work | Copyright (C) 2016-2017 OpenCFD Ltd.
                 Modified work | Copyright (C) 2018 Johan Roenby
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
    Foam::plicCutFace

Description
    Class for cutting a face, faceI, of an fvMesh, mesh_, at its intersection
    with a plicInterface.

    Reference:
        \verbatim
            Dai, Dezhi and Tong, Albert Y. (2019).
            Analytical interface reconstruction algorithms in the PLIC‐VOF
            method for 3D polyhedral unstructured meshes
            International Journal for Numerical Methods in Fluids
            doi 10.1002/fld.4750
            url https://doi.org/10.1002/fld.4750
        \endverbatim

    plicCutFace is modified from isoCutFace, whihc is part of isoAdvector.

    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

SourceFiles
    plicCutFace.C

\*---------------------------------------------------------------------------*/

#ifndef plicCutFace_H
#define plicCutFace_H

#include "fvMesh.H"
#include "plicInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plicCutFace Declaration
\*---------------------------------------------------------------------------*/

class plicCutFace
{
private:

    // Private data

        //- Ref to the fvMesh whose faces to be cut with a plicInterface
        const fvMesh& mesh_;

        //- Point along first cut edge where plicInterface cuts edge
        scalar firstEdgeCut_;

        //- Point along last cut edge where plicInterface cuts edge
        scalar lastEdgeCut_;

        //- Index in mesh_.faces()[faceI_] of first fully submerged face point
        label firstFullySubmergedPoint_;

        //- Index in mesh_.faces()[faceI_] of last fully submerged face point
        label nFullySubmergedPoints_;

        //- Storage for centre of subface
        point subFaceCentre_;

        //- Storage for area vector of subface
        vector subFaceArea_;

        //- Storage for subFacePoints
        DynamicList<point> subFacePoints_;

        //- Storage for surfacePoints (two endpoints of the intersected edge)
        DynamicList<point> surfacePoints_;

        //- Boolean telling if subface centre and area have been calculated
        bool subFaceCentreAndAreaIsCalculated_;


    // Private Member Functions

        //- Calculate centre and area vector of subface
        void calcSubFaceCentreAndArea();

        //- Calculate cut points along edges of face
        //  Returns the face status, where:
        //  -1: face is fully below the interface
        //   0: face is cut
        //  +1: face is fully above the interface
        label calcSubFace
        (
            const plicInterface& interface,
            const pointField& points,
            const labelList& pLabels
        );

        //- Calculate subface and surface points
        void subFacePoints
        (
            const pointField& points,
            const labelList& pLabels
        );

        //- Calculate surface points
        void surfacePoints
        (
            const pointField& points,
            const labelList& pLabels
        );


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct from fvMesh
        plicCutFace(const fvMesh& mesh);


    // Member functions

        //- Calculate cut points along edges of face with given label faceI
        label calcSubFace
        (
            const label faceI,
            const plicInterface& interface
        );

        //- Calculate cut points along edges of face with given point field
        label calcSubFace
        (
            const pointField& points,
            const plicInterface& interface
        );

        //- Return centre of subface
        const point& subFaceCentre();

        //- Return area vector of subface
        const vector& subFaceArea();

        //- Return subface points
        const DynamicList<point>& subFacePoints() const;

        //- Return surface points
        const DynamicList<point>& surfacePoints() const;

        //- Initialize all storage
        void clearStorage();

        //- Calculate volumetric face transport during dt with given
        //  plicInterface and its velocity
        scalar timeIntegratedFaceFlux
        (
            const label faceI,
            const plicInterface& interface,
            const scalar Un0,
            const scalar dt,
            const scalar phi,
            const scalar magSf
        );

        //- Calculate time integrated area for a face during dt
        scalar timeIntegratedArea
        (
            const pointField& fPts,
            const scalarField& pTimes,
            const scalar dt,
            const scalar magSf,
            const scalar Un0,
            const plicInterface& interface
        );

        // Calculate two endpoints of the face-interface intersection edge
        void cutPoints
        (
            const pointField& pts,
            const plicInterface& interface,
            DynamicList<point>& cutPoints
        );

        //- Calculate quadrilateral area coefficients
        void quadAreaCoeffs
        (
            const DynamicList<point>& pf0,
            const DynamicList<point>& pf1,
            scalar& alpha,
            scalar& beta
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //