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

#include "plicInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::plicInterface::typeName = "plicInterface";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::plicInterface::plicInterface()
:
    n_(vector::one),
    X_(vector::zero),
    D_(0.0)
{}


Foam::plicInterface::plicInterface(const vector& normalVector)
:
    n_(normalVector),
    X_(vector::zero),
    D_(0.0)
{}


Foam::plicInterface::plicInterface
(
    const vector& normalVector,
    const scalar signedDistance
)
:
    n_(normalVector),
    X_(vector::zero),
    D_(signedDistance)
{}


Foam::plicInterface::plicInterface
(
    const vector& normalVector,
    const point& basePoint
)
:
    n_(normalVector),
    X_(basePoint),
    D_(-(normalVector & basePoint))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::plicInterface::n() const
{
    return n_;
}


Foam::vector& Foam::plicInterface::n()
{
    return n_;
}


Foam::point Foam::plicInterface::X() const
{
    return X_;
}


Foam::point& Foam::plicInterface::X()
{
    return X_;
}


Foam::scalar Foam::plicInterface::D() const
{
    return D_;
}


Foam::scalar& Foam::plicInterface::D()
{
    return D_;
}


Foam::scalar Foam::plicInterface::signedDistance(const point& p) const
{
    return (p & n()) + D();
}


Foam::plicInterface::pointSide Foam::plicInterface::sideOfPoint(const point& p) const
{
    const scalar dist(signedDistance(p));

    if(mag(dist) <= SMALL)
        return PNORMAL;
    else
        return (dist > 0.0 ? PABOVE : PBELOW);
}


// ************************************************************************* //