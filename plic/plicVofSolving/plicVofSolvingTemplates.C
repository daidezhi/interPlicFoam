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

\*---------------------------------------------------------------------------*/

#include "plicVofSolving.H"

// ************************************************************************* //

template<typename Type>
Type Foam::plicVofSolving::faceValue
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label faceI
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        return f.primitiveField()[faceI];
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[faceI - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << faceI
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return pTraits<Type>::zero;
        }

        const label patchFacei = pp.whichFace(faceI);
        return f.boundaryField()[patchi][patchFacei];
    }
}


template<typename Type>
void Foam::plicVofSolving::setFaceValue
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& f,
    const label faceI,
    const Type& value
) const
{
    if (mesh_.isInternalFace(faceI))
    {
        f.primitiveFieldRef()[faceI] = value;
    }
    else
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

        // Boundary face. Find out which face of which patch
        const label patchi = pbm.patchID()[faceI - mesh_.nInternalFaces()];

        if (patchi < 0 || patchi >= pbm.size())
        {
            FatalErrorInFunction
                << "Cannot find patch for face " << faceI
                << abort(FatalError);
        }

        // Handle empty patches
        const polyPatch& pp = pbm[patchi];
        if (isA<emptyPolyPatch>(pp) || pp.empty())
        {
            return;
        }

        const label patchFacei = pp.whichFace(faceI);

        f.boundaryFieldRef()[patchi][patchFacei] = value;
    }
}


// ************************************************************************* //
