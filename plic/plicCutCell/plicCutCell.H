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
    Foam::plicCutCell

Description
    Class for cutting a cell, cellI, of an fvMesh, mesh_, at its intersection
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

    plicCutCell is modified from isoCutCell, whihc is part of isoAdvector.

    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

SourceFiles
    plicCutCell.C

\*---------------------------------------------------------------------------*/

#ifndef plicCutCell_H
#define plicCutCell_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "plicCutFace.H"
#include "plicInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plicCutCell Declaration
\*---------------------------------------------------------------------------*/

class plicCutCell
{
private:

    // Private data

        //- Ref to the fvMesh whose cells to be cut with a plicInterface
        const fvMesh& mesh_;

        //- Cell to cut
        label cellI_;

        //- plicCell field
        plicInterfaceField& plicInterfaceField_;

        //- A plicCutFace object to reach its face cutting functionality
        plicCutFace plicCutFace_;

        //- List of face labels of plicCutFaces
        DynamicList<label> plicCutFaces_;

        //- List of point lists each defining a plicCutFace
        DynamicList<DynamicList<point>> plicCutFacePoints_;

        //- List of face centres for plicCutFaces
        DynamicList<point> plicCutFaceCentres_;

        //- List of face area vectors for plicCutFaces
        DynamicList<vector> plicCutFaceAreas_;

        //- Storage for subFace edges belonging to interface
        DynamicList<DynamicList<point>> plicFaceEdges_;

        //- Points constituting the cell-interface intersection
        DynamicList<point> plicFacePoints_;

        //- Face centre of the interface
        point plicFaceCentre_;

        //- Face normal of the interface (from liquid to gas)
        vector plicFaceArea_;

        //- Cell centre of the subcell of cellI which is "fully submerged"
        point subCellCentre_;

        //- Volume of fully submerged subcell
        scalar subCellVolume_;

        //- Volume of Fluid for cellI (subCellVolume_ / mesh_.V()[cellI])
        scalar VOF_;

        //- List of fully submerged faces
        DynamicList<label> fullySubFaces_;

        //- A cell status label taking one of the values:
        //  -1: cell is fully below the interface
        //   0: cell is cut
        //  +1: cell is fully above the interface
        label cellStatus_;

        //- Boolean telling if subcell centre and volume have been calculated
        bool subCellCentreAndVolumeCalculated_;

        //- Boolean telling if interface centre and area have been calculated
        bool plicFaceCentreAndAreaCalculated_;


    // Private Member Functions

        //- Calculate centre and volume of subcell
        void calcSubCellCentreAndVolume();

        //- Calculate centre and area vector of plicInterface
        void calcPlicFaceCentreAndArea();

        //- Reconstruct plicInterface points from face-interface intersection
        //  edges to ensure a counter-clockwise or clockwise point label
        //  arrangement
        void calcPlicFacePointsFromEdges();


public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct from fvMesh and plicInterfaceField
        plicCutCell(const fvMesh&, plicInterfaceField&);


    // Member functions

        //- Calculate subcell
        label calcSubCell(const label cellI, const plicInterface& interface);

        //- Return subcell centre
        const point& subCellCentre();

        //- Return subcell volume
        scalar subCellVolume();

        //- Return list of point lists each defining a plicCutFace
        const DynamicList<DynamicList<point>>& plicCutFacePoints()
        {
            return plicCutFacePoints_;
        }

        //- Return label list of fully submerged faces
        const DynamicList<label>& fullySubFaces()
        {
            return fullySubFaces_;
        }

        //- Return plicInterface points
        const DynamicList<point>& plicFacePoints();

        //- Reruen plicInterface centre
        const point& plicFaceCentre();

        //- Return plicInterface area vector
        const vector& plicFaceArea();

        //- Return fraction value
        scalar volumeOfFluid();

        //- Initialize all storage
        void clearStorage();

        //- Find signed distance with given fraction value and normal vector
        label findSignedDistance
        (
            const label cellI,
            const scalar alpha1
        );

        void volumeOfFluid(volScalarField& alpha1, const plicInterface& interface);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
