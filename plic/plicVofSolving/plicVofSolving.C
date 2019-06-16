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
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "interpolationCellPointFace.H"
#include "volPointInterpolation.H"
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
#include "upwind.H"
#include "cellSet.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::plicVofSolving::typeName = "plicVofSolving";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicVofSolving::plicVofSolving
(
    volScalarField& alpha1,
    const surfaceScalarField& phi,
    const volVectorField& U
)
:
    // General data
    mesh_(alpha1.mesh()),
    dict_(mesh_.solverDict(alpha1.name())),
    alpha1_(alpha1),
    alpha1In_(alpha1.ref()),
    plicInterfaceField_(alpha1),
    phi_(phi),
    U_(U),
    dVf_
    (
        IOobject
        (
            "dVf_",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVol, 0.0)
    ),

    // CPU time
    orientationTime_(0.0),
    reconstructionTime_(0.0),
    advectionTime_(0.0),

    // Mass error
    massTotalIni_(gSum(alpha1_.primitiveField() * mesh_.V())),
    massConservationError_(0.0),

    // Tolerances and solution controls
    nAlphaBounds_(dict_.lookupOrDefault<label>("nAlphaBounds", 3)),
    surfCellTol_(dict_.lookupOrDefault<scalar>("surfCellTol", 1e-8)),
    smoothedAlphaGrad_
    (
        dict_.lookupOrDefault<bool>("smoothedAlphaGrad", false)
    ),
    writePlicFacesToFile_
    (
        dict_.lookupOrDefault<bool>("writePlicFaces", false)
    ),

    // Cell cutting data
    mixedCells_(label(0.2*mesh_.nCells())),
    cellStatus_(label(0.2*mesh_.nCells())),
    plicCutCell_(mesh_, plicInterfaceField_),
    plicCutFace_(mesh_),
    cellIsBounded_(mesh_.nCells(), false),
    checkBounding_(mesh_.nCells(), false),
    bsFaces_(label(0.2*(mesh_.nFaces() - mesh_.nInternalFaces()))),
    bsUn0_(bsFaces_.size()),
    bsInterface0_(bsFaces_.size()),

    // Parallel run data
    procPatchLabels_(mesh_.boundary().size()),
    surfaceCellFacesOnProcPatches_(0)
{
    // Prepare lists used in parallel runs
    if(Pstream::parRun())
    {
        // Force calculation of required demand driven data (else parallel
        // communication may crash)
        mesh_.cellCentres();
        mesh_.cellVolumes();
        mesh_.faceCentres();
        mesh_.faceAreas();
        mesh_.magSf();
        mesh_.boundaryMesh().patchID();
        mesh_.cellPoints();
        mesh_.cellCells();
        mesh_.cells();

        // Get boundary mesh and resize the list for parallel comms
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        surfaceCellFacesOnProcPatches_.resize(patches.size());

        // Append all processor patch labels to the list
        forAll(patches, patchi)
        {
            if
            (
                isA<processorPolyPatch>(patches[patchi]) &&
                patches[patchi].size() > 0
            )
            {
                procPatchLabels_.append(patchi);
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::plicVofSolving::timeIntegratedFlux()
{
    // Get time step
    const scalar dt = mesh_.time().deltaTValue();

    // Create object for interpolating velocity to plicface centres
    interpolationCellPoint<vector> UInterp(U_);
    //interpolationCellPointFace<vector> UInterp(U_);

    // Clear out the data for re-use and reset list containing information
    // whether cells could possibly need bounding
    checkBounding_ = false;

    // Get necessary references
    const scalarField& phiIn = phi_.primitiveField();
    const scalarField& magSfIn = mesh_.magSf().primitiveField();
    scalarField& dVfIn = dVf_.primitiveFieldRef();

    // Get necessary mesh data
    const labelListList& cellCells = mesh_.cellCells();
    const cellList& cellFaces = mesh_.cells();
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();

    // Loop through all mixed cells
    forAll(mixedCells_, cellI)
    {
        checkBounding_[mixedCells_[cellI]] = true;

        if(cellStatus_[cellI] != 0) continue;

        const plicInterface& interface0 = plicInterfaceField_.interface
                                        (
                                            mixedCells_[cellI]
                                        );
        const point& x0 = interface0.X();
        const vector& n0 = interface0.n();

        // Get the speed of the plicInterface by interpolating velocity and
        // dotting its normal vector
        const scalar Un0 = UInterp.interpolate(x0, mixedCells_[cellI]) & n0;

        // Estimate time integrated flux through each downwind face
        // Note: looping over all cell faces - in reduced-D, some of
        //       these faces will be on empty patches
        const cell& celliFaces = cellFaces[mixedCells_[cellI]];
        forAll(celliFaces, fi)
        {
            const label facei = celliFaces[fi];

            if(mesh_.isInternalFace(facei))
            {
                bool isDownwindFace = false;
                label otherCell = -1;

                if (mixedCells_[cellI] == own[facei])
                {
                    if(phiIn[facei] > 10*SMALL)
                    {
                        isDownwindFace = true;
                    }

                    otherCell = nei[facei];
                }
                else
                {
                    if(phiIn[facei] < -10*SMALL)
                    {
                        isDownwindFace = true;
                    }

                    otherCell = own[facei];
                }

                if (isDownwindFace)
                {
                    dVfIn[facei] = plicCutFace_.timeIntegratedFaceFlux
                    (
                        facei,
                        interface0,
                        Un0,
                        dt,
                        phiIn[facei],
                        magSfIn[facei]
                    );
                }

                // We want to check bounding of neighbour cells to
                // surface cells as well:
                checkBounding_[otherCell] = true;

                // Also check neighbours of neighbours.
                // Note: consider making it a run time selectable
                // extension level (easily done with recursion):
                // 0 - only neighbours
                // 1 - neighbours of neighbours
                // 2 - ...
                const labelList& nNeighbourCells = cellCells[otherCell];
                forAll(nNeighbourCells, ni)
                {
                    checkBounding_[nNeighbourCells[ni]] = true;
                }
            }
            else
            {
                bsFaces_.append(facei);
                bsUn0_.append(Un0);
                bsInterface0_.append(interface0);

                // Note: we must not check if the face is on the
                // processor patch here.
            }
        }
    }

    // Get references to boundary fields
    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const surfaceScalarField::Boundary& phib = phi_.boundaryField();
    const surfaceScalarField::Boundary& magSfb = mesh_.magSf().boundaryField();
    surfaceScalarField::Boundary& dVfb = dVf_.boundaryFieldRef();
    const label nInternalFaces = mesh_.nInternalFaces();

    // Loop through boundary surface faces
    forAll(bsFaces_, i)
    {
        // Get boundary face index (in the global list)
        const label facei = bsFaces_[i];
        const label patchi = boundaryMesh.patchID()[facei - nInternalFaces];
        const label start = boundaryMesh[patchi].start();

        if(phib[patchi].size())
        {
            const label patchFacei = facei - start;
            const scalar phiP = phib[patchi][patchFacei];

            if (phiP > 10*SMALL)
            {
                const scalar magSf = magSfb[patchi][patchFacei];

                dVfb[patchi][patchFacei] = plicCutFace_.timeIntegratedFaceFlux
                (
                    facei,
                    bsInterface0_[i],
                    bsUn0_[i],
                    dt,
                    phiP,
                    magSf
                );

                // Check if the face is on processor patch and append it to
                // the list if necessary
                checkIfOnProcPatch(facei);
            }
        }
    }

    // Synchronize processor patches
    syncProcPatches(dVf_, phi_);
}


void Foam::plicVofSolving::normaliseAndSmooth
(
    volVectorField& cellN
)
{
    const labelListList& cellPoints = mesh_.cellPoints();
    const vectorField& cellCentres = mesh_.cellCentres();
    const pointField& points = mesh_.points();

    vectorField& cellNIn = cellN.primitiveFieldRef();
    cellNIn /= (mag(cellNIn) + SMALL);

    if (smoothedAlphaGrad_)
    {
        vectorField vertexN(mesh_.nPoints(), vector::zero);
        vertexN = volPointInterpolation::New(mesh_).interpolate(cellN);
        vertexN /= (mag(vertexN) + SMALL);

        // Interpolate vertex normals back to cells
        forAll(cellNIn, celli)
        {
            const labelList& cp = cellPoints[celli];
            vector cellNi = vector::zero;
            const point& cellCentre = cellCentres[celli];
            forAll(cp, pointI)
            {
                point vertex = points[cp[pointI]];
                scalar w = 1.0/mag(vertex - cellCentre);
                cellNi += w*vertexN[cp[pointI]];
            }
            cellNIn[celli] = cellNi/(mag(cellNi) + SMALL);
        }
    }
}


void Foam::plicVofSolving::setDownwindFaces
(
    const label cellI,
    DynamicLabelList& downwindFaces
) const
{
    // Get necessary mesh data and cell information
    const labelList& own = mesh_.faceOwner();
    const cellList& cells = mesh_.cells();
    const cell& c = cells[cellI];

    downwindFaces.clear();

    // Check all faces of the cell
    forAll(c, fi)
    {
        // Get face and corresponding flux
        const label faceI = c[fi];
        const scalar phi = faceValue(phi_, faceI);

        if (own[faceI] == cellI)
        {
            if (phi > 10*SMALL)
            {
                downwindFaces.append(faceI);
            }
        }
        else if (phi < -10*SMALL)
        {
            downwindFaces.append(faceI);
        }
    }

    downwindFaces.shrink();
}


void Foam::plicVofSolving::limitFluxes()
{
    // Get time step value
    const scalar dt = mesh_.time().deltaT().value();

    volScalarField alphaNew(alpha1_ - fvc::surfaceIntegrate(dVf_));
    scalar maxAlphaMinus1 = gMax(alphaNew) - 1; // max(alphaNew - 1);
    scalar minAlpha = gMin(alphaNew);       // min(alphaNew);
    const scalar aTol = 1.0e-12;            // Note: tolerances
    const label nUndershoots = 20;          // sum(neg(alphaNew + aTol));
    const label nOvershoots = 20;           // sum(pos(alphaNew - 1 - aTol));
    cellIsBounded_ = false;

    Info<< "plicVofSolving: Before conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    // Loop number of bounding steps
    for(label n = 0; n < nAlphaBounds_; n++)
    {
        if (maxAlphaMinus1 > aTol) // Note: tolerances
        {
            surfaceScalarField dVfcorrected("dVfcorrected", dVf_);
            DynamicList<label> correctedFaces(3*nOvershoots);
            boundFromAbove(alpha1In_, dVfcorrected, correctedFaces);

            forAll(correctedFaces, fi)
            {
                label faceI = correctedFaces[fi];

                // Change to treat boundaries consistently
                setFaceValue(dVf_, faceI, faceValue(dVfcorrected, faceI));
            }

            syncProcPatches(dVf_, phi_);
        }

        if (minAlpha < -aTol) // Note: tolerances
        {
            scalarField alpha2(1.0 - alpha1In_);
            surfaceScalarField dVfcorrected
            (
                "dVfcorrected",
                phi_*dimensionedScalar("dt", dimTime, dt) - dVf_
            );

            // phi_ and dVf_ have same sign and dVf_ is the portion of
            // phi_*dt that is water.
            // dVfcorrected -= dVf_;
            // If phi_ > 0 then dVf_ > 0 and mag(phi_*dt-dVf_) < mag(phi_*dt)
            // as it should.
            // If phi_ < 0 then dVf_ < 0 and mag(phi_*dt-dVf_) < mag(phi_*dt)
            // as it should.
            DynamicList<label> correctedFaces(3*nUndershoots);
            boundFromAbove(alpha2, dVfcorrected, correctedFaces);
            forAll(correctedFaces, fi)
            {
                label faceI = correctedFaces[fi];

                // Change to treat boundaries consistently
                scalar phi = faceValue(phi_, faceI);
                scalar dVcorr = faceValue(dVfcorrected, faceI);
                setFaceValue(dVf_, faceI, phi*dt - dVcorr);
            }

            syncProcPatches(dVf_, phi_);
        }
    }
}


void Foam::plicVofSolving::boundFromAbove
(
    const scalarField& alpha1,
    surfaceScalarField& dVf,
    DynamicList<label>& correctedFaces
)
{
    correctedFaces.clear();
    scalar aTol = 10*SMALL; // Note: tolerances

    const scalarField& meshV = mesh_.cellVolumes();
    const scalar dt = mesh_.time().deltaTValue();

    DynamicList<label> downwindFaces(10);
    DynamicList<label> facesToPassFluidThrough(downwindFaces.size());
    DynamicList<scalar> dVfmax(downwindFaces.size());
    DynamicList<scalar> phi(downwindFaces.size());

    // Loop through alpha cell centred field
    forAll(alpha1, cellI)
    {
        if(checkBounding_[cellI])
        {
            const scalar Vi = meshV[cellI];
            scalar alpha1New = alpha1[cellI] - netFlux(dVf, cellI)/Vi;
            scalar alphaOvershoot = alpha1New - 1.0;
            scalar fluidToPassOn = alphaOvershoot*Vi;
            label nFacesToPassFluidThrough = 1;

            bool firstLoop = true;

            // First try to pass surplus fluid on to neighbour cells that are
            // not filled and to which dVf < phi*dt
            while (alphaOvershoot > aTol && nFacesToPassFluidThrough > 0)
            {
                facesToPassFluidThrough.clear();
                dVfmax.clear();
                phi.clear();

                cellIsBounded_[cellI] = true;

                // Find potential neighbour cells to pass surplus phase to
                setDownwindFaces(cellI, downwindFaces);

                scalar dVftot = 0;
                nFacesToPassFluidThrough = 0;

                forAll(downwindFaces, fi)
                {
                    const label facei = downwindFaces[fi];
                    const scalar phif = faceValue(phi_, facei);
                    const scalar dVff = faceValue(dVf, facei);
                    const scalar maxExtraFaceFluidTrans = mag(phif*dt - dVff);

                    // dVf has same sign as phi and so if phi>0 we have
                    // mag(phi_[facei]*dt) - mag(dVf[facei]) = phi_[facei]*dt
                    // - dVf[facei]
                    // If phi < 0 we have mag(phi_[facei]*dt) -
                    // mag(dVf[facei]) = -phi_[facei]*dt - (-dVf[facei]) > 0
                    // since mag(dVf) < phi*dt

                    if (maxExtraFaceFluidTrans/Vi > aTol)
                    {
                        // Last condition may be important because without
                        // this we will flux through uncut downwind faces
                        //if (maxExtraFaceFluidTrans/Vi > aTol &&
                        //mag(dVfIn[facei])/Vi > aTol)

                        facesToPassFluidThrough.append(facei);
                        phi.append(phif);
                        dVfmax.append(maxExtraFaceFluidTrans);
                        dVftot += mag(phif*dt);
                    }
                }

                forAll(facesToPassFluidThrough, fi)
                {
                    const label faceI = facesToPassFluidThrough[fi];
                    scalar fluidToPassThroughFace =
                        fluidToPassOn*mag(phi[fi]*dt)/dVftot;

                    nFacesToPassFluidThrough +=
                        pos(dVfmax[fi] - fluidToPassThroughFace);

                    fluidToPassThroughFace =
                        min(fluidToPassThroughFace, dVfmax[fi]);

                    scalar dVff = faceValue(dVf, faceI);
                    dVff += sign(phi[fi])*fluidToPassThroughFace;
                    setFaceValue(dVf, faceI, dVff);

                    if(firstLoop)
                    {
                        checkIfOnProcPatch(faceI);
                        correctedFaces.append(faceI);
                    }
                }

                firstLoop = false;
                alpha1New = alpha1[cellI] - netFlux(dVf, cellI)/Vi;
                alphaOvershoot = alpha1New - 1.0;
                fluidToPassOn = alphaOvershoot*Vi;
            }
        }
    }
}


Foam::scalar Foam::plicVofSolving::netFlux
(
    const surfaceScalarField& dVf,
    const label cellI
) const
{
    scalar dV = 0;

    // Get face indices
    const cell& c = mesh_.cells()[cellI];

    // Get mesh data
    const labelList& own = mesh_.faceOwner();

    forAll(c, fi)
    {
        const label facei = c[fi];
        const scalar dVff = faceValue(dVf, facei);

        if (own[facei] == cellI)
        {
            dV += dVff;
        }
        else
        {
            dV -= dVff;
        }
    }

    return dV;
}


void Foam::plicVofSolving::syncProcPatches
(
    surfaceScalarField& dVf,
    const surfaceScalarField& phi
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    if(Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        // Send
        forAll(procPatchLabels_, i)
        {
            const label patchi = procPatchLabels_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UOPstream toNbr(procPatch.neighbProcNo(), pBufs);
            const scalarField& pFlux = dVf.boundaryField()[patchi];

            const List<label>& surfCellFacesOnProcPatch =
                surfaceCellFacesOnProcPatches_[patchi];

            const UIndirectList<scalar> dVfPatch
            (
                pFlux,
                surfCellFacesOnProcPatch
            );

            toNbr << surfCellFacesOnProcPatch << dVfPatch;
        }

        pBufs.finishedSends();


        // Receive and combine
        forAll(procPatchLabels_, patchLabeli)
        {
            const label patchi = procPatchLabels_[patchLabeli];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchi]);

            UIPstream fromNeighb(procPatch.neighbProcNo(), pBufs);
            List<label> faceIDs;
            List<scalar> nbrdVfs;

            fromNeighb >> faceIDs >> nbrdVfs;

            // Combine fluxes
            scalarField& localFlux = dVf.boundaryFieldRef()[patchi];

            forAll(faceIDs, i)
            {
                const label facei = faceIDs[i];
                localFlux[facei] = - nbrdVfs[i];
            }
        }

        // Reinitialising list used for minimal parallel communication
        forAll(surfaceCellFacesOnProcPatches_, patchi)
        {
            surfaceCellFacesOnProcPatches_[patchi].clear();
        }
    }
}


void Foam::plicVofSolving::checkIfOnProcPatch(const label faceI)
{
    if(!mesh_.isInternalFace(faceI))
    {
        const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
        const label patchi = pbm.patchID()[faceI - mesh_.nInternalFaces()];

        if (isA<processorPolyPatch>(pbm[patchi]) && pbm[patchi].size())
        {
            const label patchFacei = pbm[patchi].whichFace(faceI);
            surfaceCellFacesOnProcPatches_[patchi].append(patchFacei);
        }
    }
}


void Foam::plicVofSolving::getMixedCellList()
{
    // Loop through all cells
    forAll(alpha1In_, cellI)
    {
        if(isAMixedCell(cellI))
        {
            mixedCells_.append(cellI);
            cellStatus_.append(-100);
        }
    }
}


void Foam::plicVofSolving::preProcess()
{
    // Clear out the data for re-use
    clearPlicInterfaceData();

    getMixedCellList();

    Info<< "plicVofSolving: Number of mixed cells = "
        << returnReduce(mixedCells_.size(), sumOp<label>()) << endl;
}


void Foam::plicVofSolving::orientation()
{
    scalar startTime = mesh_.time().elapsedCpuTime();

    volVectorField cellNormals("gradAlpha", fvc::grad(alpha1_));

    normaliseAndSmooth(cellNormals);

    // Loop through mixed cells
    forAll(mixedCells_, cellI)
    {
        plicInterfaceField_.interface(mixedCells_[cellI]).n() =
        -cellNormals.operator[](mixedCells_[cellI]);
    }


    // Method used in the paper
    /*
    const cellList& cells(mesh_.cells());
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const vectorField& Sf(mesh_.Sf());
    const scalarField& vol(mesh_.V());

    surfaceScalarField alphaFace("interpolateAlpha", fvc::interpolate(alpha1_));

    const scalar deltaN = (1e-8/pow(average(alpha1_.mesh().V()), 1.0/3.0)).value();

    // Loop through mixed cells
    forAll(mixedCells_, cellI)
    {
        const cell& c(cells[mixedCells_[cellI]]);

        vector gradACellI(vector(0, 0, 0));

        forAll(c, faceI)
        {
            scalar alphaFace(0);

            if(mesh_.isInternalFace(c[faceI]))
            {
                scalar weight(vol[nei[c[faceI]]] / (vol[own[c[faceI]]] +
                            vol[nei[c[faceI]]]));

                alphaFace = weight * alpha1In_[own[c[faceI]]] +
                            (1-weight) * alpha1In_[nei[c[faceI]]];
            }
            else
            {
                alphaFace = alpha1In_[mixedCells_[cellI]];
            }

            if(own[c[faceI]] == mixedCells_[cellI])
            {
                gradACellI -= alphaFace * Sf[c[faceI]];
            }
            else
            {
                gradACellI += alphaFace * Sf[c[faceI]];
            }
        }

        scalar gradMag(mag(gradACellI));

        if(gradMag < 10 * SMALL)
        {
            gradACellI = vector(1, 0, 0);
        }
        else
        {
            gradACellI /= (gradMag);
        }

        plicInterfaceField_.interface(mixedCells_[cellI]).n() = gradACellI;

    }
    */

    orientationTime_ += (mesh_.time().elapsedCpuTime() - startTime);
}


void Foam::plicVofSolving::reconstruction()
{
    scalar startTime = mesh_.time().elapsedCpuTime();

    // Storage for plicInterface points. Only used if writePlicFacesToFile_
    DynamicList<List<point> > plicFacePts;

    forAll(mixedCells_, cellI)
    {
        cellStatus_[cellI] = plicCutCell_.findSignedDistance
        (
            mixedCells_[cellI],
            alpha1In_[mixedCells_[cellI]]
        );

        if (writePlicFacesToFile_ && mesh_.time().writeTime())
        {
            plicFacePts.append(plicCutCell_.plicFacePoints());
        }
    }

    if (writePlicFacesToFile_ && mesh_.time().writeTime())
    {
        writePlicFaces(plicFacePts);
    }

    reconstructionTime_ += (mesh_.time().elapsedCpuTime() - startTime);
}


void Foam::plicVofSolving::advection()
{
    scalar startTime = mesh_.time().elapsedCpuTime();

    // Initialising dVf with upwind values
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_) * mesh_.time().deltaT();

    // Calculate volumetric face transport during dt
    timeIntegratedFlux();

    // Adjust dVf for unbounded cells
    limitFluxes();

    // Advect the free surface
    alpha1_ -= fvc::surfaceIntegrate(dVf_);
    alpha1_.correctBoundaryConditions();

    scalar maxAlphaMinus1 = gMax(alpha1In_) - 1;
    scalar minAlpha = gMin(alpha1In_);
    Info<< "plicVofSolving: After  conservative bounding: min(alpha) = "
        << minAlpha << ", max(alpha) = 1 + " << maxAlphaMinus1 << endl;

    applyBruteForceBounding();

    massConservationError_ = (gSum(alpha1_.primitiveField() * mesh_.V()) - massTotalIni_) / massTotalIni_;

    advectionTime_ += (mesh_.time().elapsedCpuTime() - startTime);
}


Foam::surfaceScalarField Foam::plicVofSolving::alphaPhi()
{
    // Initialising dVf with upwind values
    dVf_ = upwind<scalar>(mesh_, phi_).flux(alpha1_) * mesh_.time().deltaT();

    timeIntegratedFlux();

    // Adjust dVf for unbounded cells
    limitFluxes();

    return dVf_/mesh_.time().deltaT();
}


void Foam::plicVofSolving::applyBruteForceBounding()
{
    bool alpha1Changed = false;

    scalar snapAlphaTol = dict_.lookupOrDefault<scalar>("snapTol", 0.0);
    if (snapAlphaTol > 0)
    {
        alpha1_ = alpha1_
                * pos0(alpha1_ - snapAlphaTol)
                * neg0(alpha1_ - (1.0 - snapAlphaTol))
                + pos0(alpha1_ - (1.0 - snapAlphaTol));

        alpha1Changed = true;
    }

    bool clip = dict_.lookupOrDefault<bool>("clip", true);
    if(clip)
    {
        alpha1_ = min(scalar(1.0), max(scalar(0.0), alpha1_));
        alpha1Changed = true;
    }

    if(alpha1Changed)
    {
        alpha1_.correctBoundaryConditions();
    }
}


void Foam::plicVofSolving::writePlicFaces
(
    const DynamicList<List<point>>& plicFacePts
) const
{
    // Writing PLIC reconstructured poly-faces to .obj file for inspection,
    // e.g. in paraview
    // It should be noted that the Tecplot doesn't support poly obj file

    const fileName dirName
    (
        Pstream::parRun() ?
            mesh_.time().path()/".."/"plicFaces"/mesh_.time().timeName()
          : mesh_.time().path()/"plicFaces"/mesh_.time().timeName()
    );
    const word fName
    (
        "plicFaces.obj"
    );

    if (Pstream::parRun())
    {
        // Collect points from all the processors
        List<DynamicList<List<point>>> allProcFaces(Pstream::nProcs());
        allProcFaces[Pstream::myProcNo()] = plicFacePts;
        Pstream::gatherList(allProcFaces);

        if (Pstream::master())
        {
            mkDir(dirName);
            OFstream os(dirName/fName);

            if (!os.good())
            {
                FatalErrorIn
                (
                    "Foam::plicVofSolving::writePlicFaces (const DynamicList<List<point>>& plicFacePts)"
                )
                    << "Cannot open file for writing " << fName
                    << exit(FatalError);
            }

            Info<< nl << "plicVofSolving: writing PLIC faces to file: "
                << os.name() << nl << endl;

            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFacePts =
                    allProcFaces[proci];

                forAll(procFacePts, i)
                {
                    const List<point>& facePts = procFacePts[i];

                    forAll(facePts, pointI)
                    {
                        point pt(facePts[pointI]);
                        os  << "v "
                            << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
                    }
                }
            }

            os << nl;

            label iStart = 1;
            forAll(allProcFaces, proci)
            {
                const DynamicList<List<point>>& procFacePts =
                    allProcFaces[proci];

                forAll(procFacePts, i)
                {
                    const List<point>& facePts = procFacePts[i];

                    if (facePts.size() > 0)
                    {
                        os << "f";
                        forAll(facePts, pointI)
                        {
                            os << " " << iStart;
                            iStart++;
                        }

                        os << nl;
                    }
                }
            }
        }
    }
    else
    {
        mkDir(dirName);
        OFstream os(dirName/fName);

        if (!os.good())
        {
            FatalErrorIn
            (
                "Foam::plicVofSolving::writePlicFaces (const DynamicList<List<point>>& plicFacePts)"
            )
                << "Cannot open file for writing " << fName
                << exit(FatalError);
        }

        Info<< nl << "plicVofSolving: writing PLIC faces to file: "
            << os.name() << nl << endl;

        forAll(plicFacePts, faceI)
        {
            const List<point>& facePts = plicFacePts[faceI];

            forAll(facePts, pointI)
            {
                point pt(facePts[pointI]);
                os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
            }
        }

        os << nl;

        label iStart = 1;
        forAll(plicFacePts, faceI)
        {
            const List<point>& facePts = plicFacePts[faceI];

            if (facePts.size() > 0)
            {
                os << "f";
                forAll(facePts, pointI)
                {
                    os << " " << iStart;
                    iStart++;
                }

                os << nl;
            }
        }
    }
}


// ************************************************************************* //