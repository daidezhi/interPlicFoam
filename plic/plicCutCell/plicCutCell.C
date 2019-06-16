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

#include "plicCutCell.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::plicCutCell::typeName = "plicCutCell";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicCutCell::plicCutCell(const fvMesh& mesh, plicInterfaceField& pif)
:
    mesh_(mesh),
    cellI_(-1),
    plicInterfaceField_(pif),
    plicCutFace_(plicCutFace(mesh_)),
    plicCutFaces_(10),
    plicCutFacePoints_(10),
    plicCutFaceCentres_(10),
    plicCutFaceAreas_(10),
    plicFaceEdges_(10),
    plicFacePoints_(10),
    plicFaceCentre_(vector::zero),
    plicFaceArea_(vector::zero),
    subCellCentre_(vector::zero),
    subCellVolume_(-10),
    VOF_(-10),
    fullySubFaces_(10),
    cellStatus_(-1),
    subCellCentreAndVolumeCalculated_(false),
    plicFaceCentreAndAreaCalculated_(false)
{
    clearStorage();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plicCutCell::calcSubCellCentreAndVolume()
{
    if (cellStatus_ == 0)   // cell is cut
    {
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0.0;

        // Estimate the approximate cell centre as the average of face centres
        label nCellFaces(1+plicCutFaceCentres_.size()+fullySubFaces_.size());
        vector cEst = plicFaceCentre_ + sum(plicCutFaceCentres_);
        forAll(fullySubFaces_, facei)
        {
            cEst += mesh_.faceCentres()[fullySubFaces_[facei]];
        }
        cEst /= scalar(nCellFaces);


        // Contribution to subcell centre and volume from plicface
        const scalar pyr3Vol0 =
            max(mag(plicFaceArea_ & (plicFaceCentre_ - cEst)), VSMALL);

        // Calculate face-pyramid centre
        const vector pc0 = 0.75*plicFaceCentre_ + 0.25*cEst;

        // Accumulate volume-weighted face-pyramid centre
        subCellCentre_ += pyr3Vol0*pc0;

        // Accumulate face-pyramid volume
        subCellVolume_ += pyr3Vol0;

        // Contribution to subcell centre and volume from cut faces
        forAll(plicCutFaceCentres_, facei)
        {
            // Calculate 3*face-pyramid volume
            scalar pyr3Vol =
                max
                (
                    mag
                    (
                        plicCutFaceAreas_[facei]
                      & (plicCutFaceCentres_[facei] - cEst)
                    ),
                    VSMALL
                );

            // Calculate face-pyramid centre
            vector pc = 0.75*plicCutFaceCentres_[facei] + 0.25*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        // Contribution to subcell centre and volume from fully submerged faces
        forAll(fullySubFaces_, i)
        {
            const label facei = fullySubFaces_[i];
            const point& fCentre = mesh_.faceCentres()[facei];
            const vector& fArea = mesh_.faceAreas()[facei];

            // Calculate 3*face-pyramid volume
            scalar pyr3Vol = max(mag(fArea & (fCentre - cEst)), VSMALL);

            // Calculate face-pyramid centre
            vector pc = 0.75*fCentre + 0.25*cEst;

            // Accumulate volume-weighted face-pyramid centre
            subCellCentre_ += pyr3Vol*pc;

            // Accumulate face-pyramid volume
            subCellVolume_ += pyr3Vol;
        }

        subCellCentre_ /= subCellVolume_;
        subCellVolume_ /= scalar(3);
        VOF_ = subCellVolume_ / mesh_.cellVolumes()[cellI_];

        subCellCentreAndVolumeCalculated_ = true;
    }
    else if (cellStatus_ == 1)
    {
        // Cell fully above interface
        subCellCentre_ = vector::zero;
        subCellVolume_ = 0.0;
        VOF_ = 0.0;
    }
    else if (cellStatus_ == -1)
    {
        // Cell fully below interface
        subCellCentre_ = mesh_.cellCentres()[cellI_];
        subCellVolume_ = mesh_.cellVolumes()[cellI_];
        VOF_ = 1.0;
    }
}


void Foam::plicCutCell::calcPlicFaceCentreAndArea()
{
    // Initial guess of face centre from edge points
    point fCentre = vector::zero;
    label nEdgePoints = 0;
    forAll(plicFaceEdges_, ei)
    {
        DynamicList<point>& edgePoints = plicFaceEdges_[ei];
        forAll(edgePoints, pi)
        {
            fCentre += edgePoints[pi];
            nEdgePoints++;
        }
    }

    fCentre /= nEdgePoints;

    vector sumN = vector::zero;
    scalar sumA = 0.0;
    vector sumAc = vector::zero;

    forAll(plicFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = plicFaceEdges_[ei];
        const label nPoints = edgePoints.size();
        for (label pi = 0; pi < nPoints-1; pi++)
        {
            const point& nextPoint = edgePoints[pi + 1];

            vector c = edgePoints[pi] + nextPoint + fCentre;
            vector n = (nextPoint - edgePoints[pi])^(fCentre - edgePoints[pi]);
            scalar a = mag(n);

            // Edges may have different orientation
            sumN += Foam::sign(n & sumN)*n;
            sumA += a;
            sumAc += a*c;
        }
    }

    // This is to deal with zero-area faces. Mark very small faces
    // to be detected in e.g., processorPolyPatch.
    if (sumA < ROOTVSMALL)
    {
        plicFaceCentre_ = fCentre;
        plicFaceArea_ = vector::zero;
    }
    else
    {
        plicFaceCentre_ = sumAc / sumA / scalar(3);
        plicFaceArea_ = 0.5 * sumN;
    }


    // Check plicFaceArea_ direction and change if not pointing out of subcell
    if ((plicFaceArea_ & (plicFaceCentre_ - subCellCentre())) < 0.0)
    {
        plicFaceArea_ *= (-1.0);
    }

    plicFaceCentreAndAreaCalculated_ = true;
}


void Foam::plicCutCell::calcPlicFacePointsFromEdges()
{
    // Defining local coordinates with zhat along plicface normal and xhat from
    // plicface centre to first point in plicFaceEdges_
    const vector zhat = plicFaceArea_ / mag(plicFaceArea_);
    vector xhat = plicFaceEdges_[0][0] - plicFaceCentre_;
    xhat = (xhat - (xhat & zhat)*zhat);
    xhat /= mag(xhat);
    vector yhat = zhat ^ xhat;
    yhat /= mag(yhat);

    // Calculating plicface point angles in local coordinates
    DynamicList<point> unsortedPlicFacePoints(3*plicFaceEdges_.size());
    DynamicList<scalar> unsortedPlicFacePointAngles(3*plicFaceEdges_.size());
    forAll(plicFaceEdges_, ei)
    {
        const DynamicList<point>& edgePoints = plicFaceEdges_[ei];
        forAll(edgePoints, pi)
        {
            const point& p = edgePoints[pi];
            unsortedPlicFacePoints.append(p);
            unsortedPlicFacePointAngles.append
            (
                Foam::atan2
                (
                    ((p - plicFaceCentre_) & yhat),
                    ((p - plicFaceCentre_) & xhat)
                )
            );
        }
    }

    // Sorting plicface points by angle and inserting into plicFacePoints_
    labelList order(unsortedPlicFacePointAngles.size());
    Foam::sortedOrder(unsortedPlicFacePointAngles, order);
    plicFacePoints_.append(unsortedPlicFacePoints[order[0]]);
    for (label pi = 1; pi < order.size(); pi++)
    {
        if
        (
            mag
            (
                unsortedPlicFacePointAngles[order[pi]]
               -unsortedPlicFacePointAngles[order[pi-1]]
            ) > 1e-8
        )
        {
            plicFacePoints_.append(unsortedPlicFacePoints[order[pi]]);
        }
    }
}


Foam::label Foam::plicCutCell::calcSubCell
(
    const label cellI,
    const plicInterface& interface
)
{
    // Populate plicCutFaces_, plicCutFacePoints_, fullySubFaces_,
    // plicFaceCentres_ and plicFaceArea_.

    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    clearStorage();
    cellI_ = cellI;
    const cell& c = mesh_.cells()[cellI];

    forAll(c, fi)
    {
        const label faceI = c[fi];

        const label faceStatus = plicCutFace_.calcSubFace(faceI, interface);

        if (faceStatus == 0)
        {
            // Face is cut
            plicCutFacePoints_.append(plicCutFace_.subFacePoints());
            plicCutFaceCentres_.append(plicCutFace_.subFaceCentre());
            plicCutFaceAreas_.append(plicCutFace_.subFaceArea());
            plicFaceEdges_.append(plicCutFace_.surfacePoints());
        }
        else if (faceStatus == -1)
        {
            // Face fully below
            fullySubFaces_.append(faceI);
        }
    }

    if (plicCutFacePoints_.size())
    {
        // Cell cut at least at one face
        cellStatus_ = 0;
        calcPlicFaceCentreAndArea();

        // In the rare but occuring cases where a cell is only touched at a
        // point or a line the isoFaceArea_ will have zero length and here the
        // cell should be treated as either completely empty or full.
        if (mag(plicFaceArea_) < TSMALL)
        {
            if (fullySubFaces_.empty())
            {
                // Cell fully above interface
                cellStatus_ = 1;
            }
            else
            {
                // Cell fully below interface
                cellStatus_ = -1;
            }
        }
    }
    else if (fullySubFaces_.empty())
    {
        // Cell fully above interface
        cellStatus_ = 1;
    }
    else
    {
        // Cell fully below interface
        cellStatus_ = -1;
    }

    return cellStatus_;
}


const Foam::point& Foam::plicCutCell::subCellCentre()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return subCellCentre_;
}


Foam::scalar Foam::plicCutCell::subCellVolume()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return subCellVolume_;
}


const Foam::DynamicList<Foam::point>& Foam::plicCutCell::plicFacePoints()
{
    if (cellStatus_ == 0 && plicFacePoints_.size() == 0)
    {
        calcPlicFacePointsFromEdges();
    }

    return plicFacePoints_;
}


const Foam::point& Foam::plicCutCell::plicFaceCentre()
{
    if (!plicFaceCentreAndAreaCalculated_)
    {
        calcPlicFaceCentreAndArea();
    }

    return plicFaceCentre_;
}


const Foam::vector& Foam::plicCutCell::plicFaceArea()
{
    if (!plicFaceCentreAndAreaCalculated_)
    {
        calcPlicFaceCentreAndArea();
    }

    return plicFaceArea_;
}


Foam::scalar Foam::plicCutCell::volumeOfFluid()
{
    if (!subCellCentreAndVolumeCalculated_)
    {
        calcSubCellCentreAndVolume();
    }

    return VOF_;
}


void Foam::plicCutCell::clearStorage()
{
    cellI_ = -1;
    plicCutFace_.clearStorage();
    plicCutFaces_.clear();
    plicCutFacePoints_.clear();
    plicCutFaceCentres_.clear();
    plicCutFaceAreas_.clear();
    plicFaceEdges_.clear();
    plicFacePoints_.clear();
    plicFaceCentre_ = vector::zero;
    plicFaceArea_ = vector::zero;
    subCellCentre_ = vector::zero;
    subCellVolume_ = -10;
    VOF_ = -10;
    fullySubFaces_.clear();
    cellStatus_ = -1;
    subCellCentreAndVolumeCalculated_ = false;
    plicFaceCentreAndAreaCalculated_ = false;
}


Foam::label Foam::plicCutCell::findSignedDistance
(
    const label cellI,
    const scalar alpha1
)
{
    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    // Get unit normal vector of interface inside cellI
    const vector interNormal(plicInterfaceField_.interface(cellI).n());

    // Finding cell vertex extrema values
    const labelList& pLabels = mesh_.cellPoints(cellI);
    scalarField Dvert(pLabels.size());
    forAll(pLabels, pi)
    {
        Dvert[pi] = -(interNormal & mesh_.points()[pLabels[pi]]);
    }
    labelList order(Dvert.size());
    sortedOrder(Dvert, order, typename UList<scalar>::greater(Dvert));

    /*
    Binary Bracketing (BB) procedure
    Details can be found in:
        \verbatim
            L{\'o}pez, Joaqu{\'\i}n and Hern{\'a}ndez, J. (2008).
            Analytical and geometrical tools for 3D volume of fluid
            methods in general grids
            Journal of Computational Physics
            doi 10.1016/j.jcp.2008.03.010
            url https://doi.org/10.1016/j.jcp.2008.03.010
        \endverbatim
    */
    scalar DLow = Dvert[order.first()];
    scalar DUp  = Dvert[order.last()];
    label pLabelLow = 0;
    label pLabelUp  = Dvert.size() - 1;
    scalar alphaLow = 0.0;
    scalar alphaUp  = 1.0;
    scalar pLabelTmp, DTmp, alphaTmp;

    while ((pLabelUp - pLabelLow) > 1)
    {
        pLabelTmp = round(0.5 * (pLabelUp+pLabelLow));
        DTmp = Dvert[order[pLabelTmp]];
        calcSubCell(cellI, plicInterface(interNormal, DTmp));
        alphaTmp = volumeOfFluid();

        if(mag(alphaTmp-alpha1) < TSMALL)
        {
            plicInterfaceField_.interface(cellI).D() = DTmp;
            plicInterfaceField_.interface(cellI).X() = plicFaceCentre_;

            return cellStatus_;
        }

        if (alphaTmp > alpha1)
        {
            pLabelUp = pLabelTmp;
            DUp = DTmp;
            alphaUp = alphaTmp;
        }
        else
        {
            pLabelLow = pLabelTmp;
            DLow = DTmp;
            alphaLow = alphaTmp;
        }
    }


    if(mag(DLow - DUp) < TSMALL)
    {
        plicInterfaceField_.interface(cellI).D() = 0.5 * (DLow+DUp);
        calcSubCell(cellI, plicInterfaceField_.interface(cellI));
        plicInterfaceField_.interface(cellI).X() = plicFaceCentre_;

        return cellStatus_;
    }

    // Finding 2 additional points
    scalar DOneOfThree = DLow + (DUp - DLow) / scalar(3);
    calcSubCell(cellI, plicInterface(interNormal, DOneOfThree));
    scalar alphaOneOfThree = volumeOfFluid() - alphaLow;

    scalar DTwoOfThree = DLow + (DUp - DLow) * (scalar(2) / scalar(3));
    calcSubCell(cellI, plicInterface(interNormal, DTwoOfThree));
    scalar alphaTwoOfThree = volumeOfFluid() - alphaLow;

    // Calculate coefficients a, b and c by using Eq. (20)
    scalar a, b, c, d;

    scalar alphaTrap(alphaUp - alphaLow);

    a = scalar(13.5)  * alphaOneOfThree +
        scalar(-13.5) * alphaTwoOfThree +
        scalar(4.5)   * alphaTrap;

    b = scalar(-22.5) * alphaOneOfThree +
        scalar(18)    * alphaTwoOfThree +
        scalar(-4.5)  * alphaTrap;

    c = scalar(9)     * alphaOneOfThree +
        scalar(-4.5)  * alphaTwoOfThree +
        scalar(1)     * alphaTrap;

    d = alphaLow - alpha1;

    // Find the root by using Newton method
    scalar lambda = 0.5;    // Initial guess is set to 0.5
    label nIter = 0;
    while (nIter < 100)
    {
        scalar f(a*pow3(lambda) + b*sqr(lambda) + c*lambda + d);
        scalar fPrime(scalar(3)*a*sqr(lambda) + scalar(2)*b*lambda + c);
        scalar lambdaNew(lambda - (f / fPrime));

        // Convergence tolerance is 1e-14
        if (mag(lambdaNew - lambda) < TSMALL)
        {
            break;
        }
        lambda = lambdaNew;
        nIter++;
    }

    // Calculate $D_0$
    scalar D0(DLow + (DUp - DLow) * lambda);

    // Update subcell with $D_0$
    plicInterfaceField_.interface(cellI).D() = D0;
    calcSubCell(cellI, plicInterfaceField_.interface(cellI));
    plicInterfaceField_.interface(cellI).X() = plicFaceCentre_;

    return cellStatus_;
}


void Foam::plicCutCell::volumeOfFluid
(
    volScalarField& alpha1,
    const plicInterface& interface
)
{
    // Setting internal field
    scalarField& alphaIn = alpha1;
    forAll(alphaIn, celli)
    {
        const label cellStatus = calcSubCell(celli, interface);
        if (cellStatus != 1)
        {
            // If cell not entirely above plicsurface
            alphaIn[celli] = volumeOfFluid();
        }
    }

    // Setting boundary alpha1 values
    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].size() > 0)
        {
            const label start = mesh_.boundary()[patchi].patch().start();
            scalarField& alphap = alpha1.boundaryFieldRef()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

            forAll(alphap, patchFacei)
            {
                const label facei = patchFacei + start;
                const label faceStatus = plicCutFace_.calcSubFace(facei, interface);

                if (faceStatus != 1)
                {
                    // Face not entirely above plicsurface
                    alphap[patchFacei] =
                        mag(plicCutFace_.subFaceArea())/magSfp[patchFacei];
                }
            }
        }
    }
}


// ************************************************************************* //