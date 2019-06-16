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

#include "plicCutFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::plicCutFace::typeName = "plicCutFace";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicCutFace::plicCutFace
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    firstEdgeCut_(-1),
    lastEdgeCut_(-1),
    firstFullySubmergedPoint_(-1),
    nFullySubmergedPoints_(0),
    subFaceCentre_(vector::zero),
    subFaceArea_(vector::zero),
    subFacePoints_(10),
    surfacePoints_(4),
    subFaceCentreAndAreaIsCalculated_(false)
{
    clearStorage();
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::plicCutFace::calcSubFaceCentreAndArea()
{
    const label nPoints = subFacePoints_.size();

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        subFaceCentre_ = sum(subFacePoints_) / scalar(3);
        subFaceArea_ =
            0.5
           *(
                (subFacePoints_[1] - subFacePoints_[0])
               ^(subFacePoints_[2] - subFacePoints_[0])
            );
    }
    else if (nPoints > 0)
    {
        vector sumN = vector::zero;
        scalar sumA = 0.0;
        vector sumAc = vector::zero;
        const point fCentre = sum(subFacePoints_) / scalar(nPoints);

        for (label pi = 0; pi < nPoints; pi++)
        {
            const point& nextPoint = subFacePoints_[subFacePoints_.fcIndex(pi)];

            vector c = subFacePoints_[pi] + nextPoint + fCentre;
            vector n = (nextPoint - subFacePoints_[pi])
                      ^(fCentre - subFacePoints_[pi]);
            scalar a = magSqr(n);

            sumN += n;
            sumA += a;
            sumAc += a * c;
        }

        // This is to deal with zero-area faces. Mark very small faces
        // to be detected in e.g., processorPolyPatch.
        if (sumA < ROOTVSMALL)
        {
            subFaceCentre_ = fCentre;
            subFaceArea_ = vector::zero;
        }
        else
        {
            subFaceCentre_ = (1.0/3.0) * sumAc / sumA;
            subFaceArea_ = 0.5 * sumN;
        }
    }

    subFaceCentreAndAreaIsCalculated_ = true;
}


Foam::label Foam::plicCutFace::calcSubFace
(
    const plicInterface& interface,
    const pointField& points,
    const labelList& pLabels
)
{
    const scalar TSMALL(10.0*SMALL);

    // Face status set to one of the values:
    //  -1: face is fully below the interface
    //   0: face is cut
    //  +1: face is fully above the interface
    label faceStatus;

    const label nPoints = pLabels.size();
    scalarList r_(nPoints, 0.0);

    // Compute the signed distances from face points to the interface
    forAll(pLabels, pi)
    {
        r_[pi] = interface.signedDistance(points[pLabels[pi]]);
    }

    label pl1 = 0;

    // If vertex is very close to interface, lift them slightly to avoid
    // dealing with the many special cases of a face being touched either
    // at a single point, along an edge, or the entire face being on the
    // interface.
    if(mag(r_[pl1]) < TSMALL)
    {
        r_[pl1] += sign(r_[pl1]) * TSMALL;
    }

    forAll(pLabels, pi)
    {
        label pl2 = (pi + 1) % nPoints;

        if(mag(r_[pl2]) < TSMALL)
        {
            r_[pl2] += sign(r_[pl2]) * TSMALL;
        }

        if(r_[pl1] < 0.0)
        {
            nFullySubmergedPoints_ += 1;

            if(r_[pl2] > 0.0)
            {
                lastEdgeCut_ = r_[pl1] / (r_[pl1] - r_[pl2]);
            }
        }
        else if(r_[pl1] > 0.0 && r_[pl2] < 0.0)
        {
            if(firstFullySubmergedPoint_ == -1)
            {
                firstFullySubmergedPoint_ = pl2;

                firstEdgeCut_ = r_[pl1] / (r_[pl1] - r_[pl2]);
            }
        }
        pl1 = pl2;

    }

    if(firstFullySubmergedPoint_ != -1)
    {
        // Face is cut
        faceStatus = 0;
        subFacePoints(points, pLabels);
    }
    else if(r_[pl1] > 0.0)
    {
        // Face entirely above interface
        faceStatus = 1;
    }
    else
    {
        // Face entirely below interface
        faceStatus = -1;
    }

    return faceStatus;
}


void Foam::plicCutFace::subFacePoints
(
    const pointField& points,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();

    surfacePoints(points, pLabels);

    forAll(surfacePoints_, pi)
    {
        subFacePoints_.append(surfacePoints_[pi]);
    }

    for (label pi = 0; pi < nFullySubmergedPoints_; pi++)
    {
        subFacePoints_.append
        (
            points[pLabels[(firstFullySubmergedPoint_ + pi) % nPoints]]
        );
    }
}


void Foam::plicCutFace::surfacePoints
(
    const pointField& points,
    const labelList& pLabels
)
{
    const label nPoints = pLabels.size();

    const label n = firstFullySubmergedPoint_ + nFullySubmergedPoints_;

    label pl1 = pLabels[(n - 1) % nPoints];

    label pl2 = pLabels[n % nPoints];

    surfacePoints_.append
    (
        points[pl1] + lastEdgeCut_ * (points[pl2] - points[pl1])
    );

    pl1 = pLabels[(firstFullySubmergedPoint_ - 1 + nPoints) % nPoints];
    pl2 = pLabels[firstFullySubmergedPoint_];

    surfacePoints_.append
    (
        points[pl1] + firstEdgeCut_ * (points[pl2] - points[pl1])
    );
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::plicCutFace::calcSubFace
(
    const label faceI,
    const plicInterface& interface
)
{
    clearStorage();
    const labelList& pLabels = mesh_.faces()[faceI];
    const pointField& points = mesh_.points();
    return calcSubFace(interface, points, pLabels);
}


Foam::label Foam::plicCutFace::calcSubFace
(
    const pointField& points,
    const plicInterface& interface
)
{
    clearStorage();
    const labelList pLabels(identity(points.size()));
    return calcSubFace(interface, points, pLabels);
}


const Foam::point& Foam::plicCutFace::subFaceCentre()
{
    if (!subFaceCentreAndAreaIsCalculated_)
    {
        calcSubFaceCentreAndArea();
    }
    return subFaceCentre_;
}


const Foam::vector& Foam::plicCutFace::subFaceArea()
{
    if (!subFaceCentreAndAreaIsCalculated_)
    {
        calcSubFaceCentreAndArea();
    }
    return subFaceArea_;
}


const Foam::DynamicList<Foam::point>& Foam::plicCutFace::subFacePoints() const
{
    return subFacePoints_;
}


const Foam::DynamicList<Foam::point>& Foam::plicCutFace::surfacePoints() const
{
    return surfacePoints_;
}


void Foam::plicCutFace::clearStorage()
{
    firstEdgeCut_ = -1;
    lastEdgeCut_ = -1;
    firstFullySubmergedPoint_ = -1;
    nFullySubmergedPoints_ = 0;
    subFaceCentre_ = vector::zero;
    subFaceArea_ = vector::zero;
    subFacePoints_.clear();
    surfacePoints_.clear();
    subFaceCentreAndAreaIsCalculated_ = false;
}


Foam::scalar Foam::plicCutFace::timeIntegratedFaceFlux
(
    const label faceI,
    const plicInterface& interface,
    const scalar Un0,
    const scalar dt,
    const scalar phi,
    const scalar magSf
)
{
    const scalar TSMALL(10.0*SMALL);

    // Find sorted list of times where the interface will arrive at face
    // points with given initial position x0 and normal velocity Un0.

    // Get points for this face
    const face& f = mesh_.faces()[faceI];
    const pointField fPts(f.points(mesh_.points()));
    const label nPoints = fPts.size();

    // Get initial position x0 and unit normal vector n0 of interface
    const point& x0 = interface.X();
    const vector& n0 = interface.n();

    scalarField pTimes(fPts.size());
    if(mag(Un0) > TSMALL)
    {
        // Here we estimate time of arrival to the face points from their
        // normal distance to the initial interface and the interface normal
        // velocity.

        pTimes = ((fPts - x0) & n0) / Un0;

        scalar dVf = 0.0;

        // Check if pTimes changes direction more than twice when looping face
        // this may happen when the face holds warpage
        label nShifts = 0;
        forAll(pTimes, pi)
        {
            const label oldEdgeSign =
                sign(pTimes[(pi + 1) % nPoints] - pTimes[pi]);
            const label newEdgeSign =
                sign(pTimes[(pi + 2) % nPoints] - pTimes[(pi + 1) % nPoints]);

            if (newEdgeSign != oldEdgeSign)
            {
                nShifts++;
            }
        }
        if (nShifts == 2)
        {
            dVf = phi / magSf * timeIntegratedArea
                                (
                                    fPts,
                                    pTimes,
                                    dt,
                                    magSf,
                                    Un0,
                                    interface
                                );
        }
        else if (nShifts > 2)
        {
            // Triangle decompose the face
            pointField fPts_tri(3);
            scalarField pTimes_tri(3);
            fPts_tri[0] = mesh_.faceCentres()[faceI];
            pTimes_tri[0] = ((fPts_tri[0] - x0) & n0) / Un0;
            for (label pi = 0; pi < nPoints; pi++)
            {
                fPts_tri[1] = fPts[pi];
                pTimes_tri[1] = pTimes[pi];
                fPts_tri[2] = fPts[(pi + 1) % nPoints];
                pTimes_tri[2] = pTimes[(pi + 1) % nPoints];
                const scalar magSf_tri =
                    mag
                    (
                        0.5
                       *(fPts_tri[2] - fPts_tri[0])
                       ^(fPts_tri[1] - fPts_tri[0])
                    );
                const scalar phi_tri = phi*magSf_tri/magSf;
                dVf += phi_tri / magSf_tri *timeIntegratedArea
                                            (
                                                fPts_tri,
                                                pTimes_tri,
                                                dt,
                                                magSf_tri,
                                                Un0,
                                                interface
                                            );
            }
        }

        return dVf;
    }
    else
    {
        // Un0 is almost zero and interface is treated as stationary
        calcSubFace(faceI, interface);
        const scalar alphaf = mag(subFaceArea() / magSf);

        return phi*dt*alphaf;
    }
}


Foam::scalar Foam::plicCutFace::timeIntegratedArea
(
    const pointField& fPts,
    const scalarField& pTimes,
    const scalar dt,
    const scalar magSf,
    const scalar Un0,
    const plicInterface& interface
)
{
    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    // Initialise time integrated area returned by this function
    scalar tIntArea = 0.0;

    // Finding ordering of vertex points
    labelList order(pTimes.size());
    sortedOrder(pTimes, order);
    const scalar firstTime = pTimes[order.first()];
    const scalar lastTime = pTimes[order.last()];

    // Dealing with case where face is not cut by interface during time
    // interval [0, dt] because face was already passed by interface
    if (lastTime <= 0.0)
    {
        // If all face cuttings were in the past and cell is filling up (Un0>0)
        // then face must be full during whole time interval
        tIntArea = magSf * dt * scalar(pos0(Un0));
        return tIntArea;
    }

    // Dealing with case where face is not cut by interface during time
    // interval [0, dt] because dt is too small for interface to reach closest
    // face point
    if (firstTime >= dt)
    {
        // If all cuttings are in the future but non of them within [0, dt]
        // then if cell is filling up (Un0 > 0) face must be empty during
        // whole time interval
        tIntArea = magSf * dt * (1.0 - scalar(pos0(Un0)));
        return tIntArea;
    }

    // If we reach this point in the code some part of the face will be swept
    // during [tSmall, dt-tSmall]. However, it may be the case that there are
    // no vertex times within the interval. This will happen sometimes for
    // small time steps where both the initial and the final face-interface
    // intersection line (FIIL) will be along the same two edges.

    // Face-interface intersection line (FIIL) to be swept across face
    DynamicList<point> FIIL(3);
    // Submerged area at beginning of each sub time interval time
    scalar initialArea = 0.0;
    //Running time keeper variable for the integration process
    scalar time = 0.0;

    // Special treatment of first sub time interval
    if (firstTime > 0)
    {
        // If firstTime > 0 the face is uncut in the time interval
        // [0, firstTime] and hence fully submerged in fluid A or B.
        // If Un0 > 0 cell is filling up and it must initially be empty.
        // If Un0 < 0 cell must initially be fully immersed in fluid A.
        time = firstTime;
        initialArea = magSf * (1.0 - scalar(pos0(Un0)));
        tIntArea = initialArea*time;
        cutPoints
        (
            fPts,
            plicInterface(interface.n(), interface.D() - time*Un0),
            FIIL
        );
    }
    else
    {
        // If firstTime <= 0 then face is initially cut and we must
        // calculate the initial submerged area and FIIL:
        time = 0.0;
        // Note: calcSubFace assumes well-defined 2-point FIIL!!!!
        calcSubFace(fPts, interface);
        initialArea = mag(subFaceArea());
        cutPoints
        (
            fPts,
            plicInterface(interface.n(), interface.D() - time*Un0),
            FIIL
        );
    }

    // Making sorted array of all vertex times that are between
    // max(0,firstTime) and dt and further than tSmall from the previous time.
    DynamicList<scalar> sortedTimes(pTimes.size());
    {
        scalar prevTime = time;
        const scalar tSmall = max(1e-6*dt, TSMALL);
        forAll(order, ti)
        {
            const scalar timeI = pTimes[order[ti]];
            if ( timeI > prevTime + tSmall && timeI <= dt)
            {
                sortedTimes.append(timeI);
                prevTime = timeI;
            }
        }
    }

    // Sweeping all quadrilaterals corresponding to the intervals defined above
    forAll(sortedTimes, ti)
    {
        const scalar newTime = sortedTimes[ti];
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints
        (
            fPts,
            plicInterface(interface.n(), interface.D() - newTime*Un0),
            newFIIL
        );

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea += (newTime - time) *
            (initialArea + sign(Un0) * (alpha/3.0 + 0.5*beta));
        // Adding quad area to submerged area
        initialArea += sign(Un0) * (alpha + beta);

        FIIL = newFIIL;
        time = newTime;
    }

    // Special treatment of last time interval
    if (lastTime > dt)
    {
        // FIIL will end up cutting the face at dt
        // New face-interface intersection line
        DynamicList<point> newFIIL(3);
        cutPoints
        (
            fPts,
            plicInterface(interface.n(), interface.D() - dt*Un0),
            newFIIL
        );

        // quadrilateral area coefficients
        scalar alpha = 0, beta = 0;
        quadAreaCoeffs(FIIL, newFIIL, alpha, beta);
        // Integration of area(t) = A*t^2+B*t from t = 0 to 1
        tIntArea += (dt - time) *
            (initialArea + sign(Un0) * (alpha/3.0 + 0.5*beta));
    }
    else
    {
        // FIIL will leave the face at lastTime and face will be fully in fluid
        // A or fluid B in the time interval from lastTime to dt.
        tIntArea += magSf * (dt - lastTime) * scalar(pos0(Un0));
    }

    return tIntArea;
}


void Foam::plicCutFace::cutPoints
(
    const pointField& pts,
    const plicInterface& interface,
    DynamicList<point>& cutPoints
)
{
    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    const label nPoints = pts.size();
    scalarList r_(nPoints, 0.0);

    forAll(pts, pi)
    {
        r_[pi] = interface.signedDistance(pts[pi]);
    }

    scalar r1(r_[0]);

    // Snapping vertex value to f0 if very close (needed for 2D cases)
    if(mag(r1) < TSMALL)
    {
        r1 = 0.0;
    }

    forAll(pts, pi)
    {
        label pi2 = (pi + 1) % nPoints;
        scalar r2 = r_[pi2];

        // Snapping vertex value
        if(mag(r2) < TSMALL)
        {
            r2 = 0.0;
        }

        if((r1 > 0.0 && r2 < 0.0) || (r1 < 0.0 && r2 > 0.0))
        {
            const scalar s = r1 / (r1 - r2);
            cutPoints.append(pts[pi] + s * (pts[pi2] - pts[pi]));
        }
        else if(r1 == 0.0)
        {
            cutPoints.append(pts[pi]);
        }
        r1 = r2;
    }
}


void Foam::plicCutFace::quadAreaCoeffs
(
    const DynamicList<point>& pf0,
    const DynamicList<point>& pf1,
    scalar& alpha,
    scalar& beta
) const
{
    // Tolerance
    const scalar TSMALL(10.0*SMALL);

    // Number of points in provided face-interface intersection lines
    const label np0 = pf0.size();
    const label np1 = pf1.size();

    // quad area coeffs such that area(t) = alpha*t^2 + beta*t.
    // With time interval normalised, we have full quadArea = alpha + beta
    // and time integrated quad area = alpha/3 + beta/2;
    alpha = 0.0;
    beta = 0.0;

    if (np0 && np1)
    {
        // Initialising quadrilateral vertices A, B, C and D
        vector A(pf0[0]);
        vector C(pf1[0]);
        vector B(pf0[0]);
        vector D(pf1[0]);

        if (np0 == 2)
        {
            B = pf0[1];
        }
        else if (np0 > 2)
        {
            WarningInFunction
                << "Vertex face was cut at pf0 = " << pf0 << endl;
        }

        if (np1 == 2)
        {
            D = pf1[1];
        }
        else if (np1 > 2)
        {
            WarningInFunction
                << "Vertex face was cut at pf1 = " << pf1 << endl;
        }

        // Swapping pf1 points if pf0 and pf1 point in same general direction
        // (because we want a quadrilateral ABCD where pf0 = AB and pf1 = CD)
        if (((B - A) & (D - C)) > 0)
        {
            vector tmp = D;
            D = C;
            C = tmp;
        }

        // Defining local coordinates (xhat, yhat) for area integration of
        // swept quadrilateral ABCD such that A = (0,0), B = (Bx,0),
        // C = (Cx,Cy) and D = (Dx,Dy) with Cy = 0 and Dy > 0.

        const scalar Bx = mag(B - A);

        vector xhat(vector::zero);
        if (Bx > TSMALL)
        {
            // If |AB| > 0 ABCD we use AB to define xhat
            xhat = (B - A) / mag(B - A);
        }
        else if (mag(C - D) > TSMALL)
        {
            // If |AB| ~ 0 ABCD is a triangle ACD and we use CD for xhat
            xhat = (C - D) / mag(C - D);
        }
        else
        {
            return;
        }

        // Defining vertical axis in local coordinates
        vector yhat = D - A;
        yhat -= ((yhat & xhat) * xhat);

        if (mag(yhat) > TSMALL)
        {
            yhat /= mag(yhat);

            const scalar Cx = (C - A) & xhat;
            const scalar Cy = mag((C - A) & yhat);
            const scalar Dx = (D - A) & xhat;
            const scalar Dy = mag((D - A) & yhat);

            // area = ((Cx - Bx)*Dy - Dx*Cy)/6.0 + 0.25*Bx*(Dy + Cy);
            alpha = 0.5 * ((Cx - Bx) * Dy - Dx * Cy);
            beta = 0.5 * Bx * (Dy + Cy);
        }
    }
    else
    {
        WarningInFunction
            << "Vertex face was cut at " << pf0 << " and at " << pf1 << endl;
    }
}


// ************************************************************************* //