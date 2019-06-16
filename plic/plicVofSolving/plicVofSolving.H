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
    Foam::plicVofSolving

Description
    Calculates the new VOF (alpha) field after time step dt with given initial
    VOF field, a velocity field U and face fluxes phi. The fluid transport
    calculation is based on PLIC and advecting the interface across the mesh
    faces with the velocity field interpolated to the interfaces.

    Reference:
        \verbatim
            Dai, Dezhi and Tong, Albert Y. (2019).
            Analytical interface reconstruction algorithms in the PLIC‐VOF
            method for 3D polyhedral unstructured meshes
            International Journal for Numerical Methods in Fluids
            doi 10.1002/fld.4750
            url https://doi.org/10.1002/fld.4750
        \endverbatim

    plicVofSolving is modified from isoAdvection, whihc is part of isoAdvector.

    Reference:
        \verbatim
            Roenby, J., Bredmose, H. and Jasak, H. (2016).
            A computational method for sharp interface advection
            Royal Society Open Science, 3
            doi 10.1098/rsos.160405
        \endverbatim

SourceFiles
    plicVofSolving.C
    plicVofSolvingTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef plicVofSolving_H
#define plicVofSolving_H

#include "fvMesh.H"
#include "volFieldsFwd.H"
#include "surfaceFields.H"
#include "className.H"
#include "plicCutCell.H"
#include "plicCutFace.H"
#include "plicInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class plicVofSolving Declaration
\*---------------------------------------------------------------------------*/

class plicVofSolving
{
private:

protected:

    // Private data types

        typedef DynamicList<label>          DynamicLabelList;
        typedef DynamicList<scalar>         DynamicScalarList;
        typedef DynamicList<vector>         DynamicVectorList;
        typedef DynamicList<point>          DynamicPointList;
        typedef DynamicList<plicInterface>  DynamicPlicInterfaceList;


    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Dictionary for plicVofSolving controls
        const dictionary dict_;

        //- Reference to alpha1 (VOF) field
        volScalarField& alpha1_;

        //- Often used reference to alpha1 internal field
        scalarField& alpha1In_;

        //- plicInterface field
        plicInterfaceField plicInterfaceField_;

        //- Reference to flux field
        const surfaceScalarField& phi_;

        //- Reference to velocity field
        const volVectorField& U_;

        //- Face volumetric liquid transport
        surfaceScalarField dVf_;

        //- Time spent calculating the interface orientations
        scalar orientationTime_;

        //- Time spent reconstructing the interfaces
        scalar reconstructionTime_;

        //- Time spent advecting the fraction field
        scalar advectionTime_;

        //- Total mass at initial time
        scalar massTotalIni_;

        //- Mass conservation error
        scalar massConservationError_;


        // Switches and tolerances

            //- Number of alpha bounding steps
            label nAlphaBounds_;

            //- Tolerance for marking surface (mixed) cells:
            //  Those with surfCellTol_ < alpha1 < 1 - surfCellTol_
            scalar surfCellTol_;

            //- Switch controlling whether to use a smoothed alpha gradient
            bool smoothedAlphaGrad_;

            //- Print plicfaces in a <case>/plicFaces/time/plicFaces.obj file.
            //  Intended for post-process
            bool writePlicFacesToFile_;


        // Cell and face cutting

            //- List of surface cell labels
            DynamicLabelList mixedCells_;

            //- List of surface cell status
            DynamicLabelList cellStatus_;

            //- Cell cutting object
            plicCutCell plicCutCell_;

            //- Face cutting object
            plicCutFace plicCutFace_;

            //- Bool list for cells that have been touched by bounding step
            boolList cellIsBounded_;

            //- True for all surface cells and their neighbours
            boolList checkBounding_;

            //- Storage for boundary faces downwind to a surface cell
            DynamicLabelList bsFaces_;

            //- Storage for boundary surface plicInterface speed
            DynamicScalarList bsUn0_;

            //- Storage for boundary surface plicInterface
            DynamicPlicInterfaceList bsInterface0_;


        // Additional data for parallel runs

            //- List of processor patch labels
            DynamicLabelList procPatchLabels_;

            //- For each patch if it is a processor patch this is a list of the
            //  face labels on this patch that are downwind to a surface cell.
            //  For non-processor patches the list will be empty.
            List<DynamicLabelList> surfaceCellFacesOnProcPatches_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        plicVofSolving(const plicVofSolving&) = delete;

        //- Disallow default bitwise copy assignment
        void operator=(const plicVofSolving&) = delete;


        // VofSolving functions

            //- For each face calculate volumetric face transport during dt
            void timeIntegratedFlux();

            //- Function used to normalise and smoothen grad(alpha)
            void normaliseAndSmooth
            (
                volVectorField& cellN
            );

            //- For a given cell return labels of faces fluxing out of this
            //  cell (based on sign of phi)
            void setDownwindFaces
            (
                const label cellI,
                DynamicLabelList& downwindFaces
            ) const;

            // Limit fluxes
            void limitFluxes();

            // Bound fluxes
            void boundFromAbove
            (
                const scalarField& alpha1,
                surfaceScalarField& dVfcorrected,
                DynamicLabelList& correctedFaces
            );

            //- Given the face volume transport dVf calculates the total volume
            //  leaving a given cell. Note: cannot use dVf member because
            //  netFlux is called also for corrected dVf
            scalar netFlux
            (
                const surfaceScalarField& dVf,
                const label cellI
            ) const;

            //- Determine if a cell is a surface (mixed) cell
            bool isAMixedCell(const label cellI) const
            {
                return
                (
                    surfCellTol_ < alpha1In_[cellI] &&
                    alpha1In_[cellI] < 1.0 - surfCellTol_
                );
            }

            //- Determine if a cell is an empty cell
            bool isAnEmptyCell(const label cellI) const
            {
                return
                (
                    alpha1In_[cellI] <= surfCellTol_
                );
            }

            //- Determine if a cell is a full cell
            bool isAFullCell(const label cellI) const
            {
                return
                (
                    alpha1In_[cellI] >= 1.0 - surfCellTol_
                );
            }

            //- Clear out plicInterface data
            void clearPlicInterfaceData()
            {
                mixedCells_.clear();
                cellStatus_.clear();
                bsFaces_.clear();
                bsUn0_.clear();
                bsInterface0_.clear();

                if (mesh_.topoChanging())
                {
                    // Introduced resizing to cope with changing meshes
                    checkBounding_.resize(mesh_.nCells());
                    cellIsBounded_.resize(mesh_.nCells());
                }

                checkBounding_ = false;
                cellIsBounded_ = false;
            }


        // Face value functions needed for random face access where the face
        // can be either internal or boundary face

            //- Return face value for a given Geometric surface field
            template<typename Type>
            Type faceValue
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label faceI
            ) const;

            //- Set face value for a given Geometric surface field
            template<typename Type>
            void setFaceValue
            (
                GeometricField<Type, fvsPatchField, surfaceMesh>& f,
                const label faceI,
                const Type& value
            ) const;


        // Parallel run handling functions

            //- Synchronize dVf across processor boundaries using upwind value
            void syncProcPatches
            (
                surfaceScalarField& dVf,
                const surfaceScalarField& phi
            );

            //- Check if the face is on processor patch and append it to the
            //  list of surface cell faces on processor patches
            void checkIfOnProcPatch(const label faceI);


public:

    // Static data members

        static const char* const typeName;


    //- Constructors

        //- Construct given alpha, phi and velocity field.
        //  Note: phi should be divergence free up to a sufficient tolerance
        plicVofSolving
        (
            volScalarField& alpha1,
            const surfaceScalarField& phi,
            const volVectorField& U
        );


    //- Destructor
    virtual ~plicVofSolving()
    {}


    // Member functions

        //- Get label list of mixed cells
        void getMixedCellList();

        //- Initialization
        void preProcess();


        // Three fundamental steps in PLIC-VOF method

            //- Calculate interface unit normal vectors
            void orientation();

            //- Reconstruct the interfaces
            void reconstruction();

            //- Advect the free surface. Updates alpha field, taking into
            //  account multiple calls within a single time step.
            void advection();

        //- Calculate alpha flux
        surfaceScalarField alphaPhi();


        //- Apply the bounding based on user inputs
        void applyBruteForceBounding();


        // Access functions

            //- Return alpha field
            const volScalarField& alpha() const
            {
                return alpha1_;
            }

            //- Return the controls dictionary
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return mass flux
            tmp<surfaceScalarField> getRhoPhi
            (
                const dimensionedScalar rho1,
                const dimensionedScalar rho2
            ) const
            {
                return tmp<surfaceScalarField>
                (
                    new surfaceScalarField
                    (
                        "rhoPhi",
                        (rho1 - rho2)*dVf_/mesh_.time().deltaT() + rho2*phi_
                    )
                );
            }


            //- Get VOF solving time
            scalar VofSolvingTime() const
            {
                return (orientationTime_ +
                        reconstructionTime_ +
                        advectionTime_);
            }

            //- Get time in orientation step
            scalar orientationTime() const
            {
                return orientationTime_;
            }

            //- Get time in reconstruction step
            scalar reconstructionTime() const
            {
                return reconstructionTime_;
            }

            //- Get time in advection step
            scalar advectionTime() const
            {
                return advectionTime_;
            }

            //- Get mass conservation error
            scalar massConservationError() const
            {
                return massConservationError_;
            }

            //- Write plicface points to .obj file
            void writePlicFaces
            (
                const DynamicList<List<point>>& plicFacePts
            ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "plicVofSolvingTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //