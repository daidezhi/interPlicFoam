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

\*---------------------------------------------------------------------------*/

#include "plicInterfaceField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::plicInterfaceField::typeName = "plicInterfaceField";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicInterfaceField::plicInterfaceField(volScalarField& alpha1)
:
    size_(alpha1.mesh().nCells())
{
    this->plicInterfaces_ = new plicInterface[size_];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::plicInterface& Foam::plicInterfaceField::operator[](const label i)
{

    return plicInterfaces_[i];
}


const Foam::plicInterface& Foam::plicInterfaceField::operator[](const label i) const
{

    return plicInterfaces_[i];
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::plicInterface& Foam::plicInterfaceField::interface(const label i)
{

    return plicInterfaces_[i];
}


const Foam::plicInterface& Foam::plicInterfaceField::interface(const label i) const
{

    return plicInterfaces_[i];
}


// ************************************************************************* //