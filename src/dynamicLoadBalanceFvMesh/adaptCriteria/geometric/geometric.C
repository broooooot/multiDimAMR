/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Henning Scheufler, DLR  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "geometric.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "volPointInterpolation.H"
#include "labelIOField.H"
#include "syncTools.H"
#include "volumeType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(geometric, 0);
addToRunTimeSelectionTable
(
    adaptCriteria,
    geometric,
    dictionary
);

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geometric::geometric
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    adaptCriteria(mesh, dict),
    surftype_(dict.get<word>("surftype")),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.get<word>("surftype"),
            IOobject
            (
                "surfaces",
                mesh.time().constant(),
                "triSurface",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    )
{

}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::geometric::~geometric()
{}

// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //
void Foam::geometric::reReadDictionary
(
    const dictionary& reReadDict
)
{
    if (reReadDict.get<word>("surftype") != surftype_)
    {
        FatalErrorInFunction
            << "The surface type can not be changed during runtime"
            << exit(FatalError);
    }
}

Foam::bitSet
Foam::geometric::refinementCellCandidates() const
{
    Info
        << "Selecting refinement candidate cells based on "
        << surftype_<< endl;

    bitSet refineCells(mesh().nCells(),false);

    List<volumeType> insideGeom;
    surfPtr_->getVolumeType(mesh().cellCentres(), insideGeom);

    // get cellLevels
    const labelIOList& curCellLevel =
        mesh().lookupObject<labelIOList>("cellLevel");

    forAll(insideGeom, cellI)
    {
        if
        (
            insideGeom[cellI] == volumeType::INSIDE
            && curCellLevel[cellI] < maxCellLevel_
        )
        {
            // Cell value is within the bounds, append cell for potential
            // refinement
            refineCells.set(cellI);
        }
    }

    const label nLocalCandidates = refineCells.count();
    const label nCandidates = returnReduce(nLocalCandidates, sumOp<label>());

    Info
        << "-> "<< nCandidates<< " cells selected"<< endl;

    // Print out some information
    // Info<< "Selection algorithm " << type() << " selected "
    //     << returnReduce(refineCells.toc().size(), sumOp<label>())
    //     << " cells as refinement candidates."
    //     << endl;

    // Return the list in the Xfer container to prevent copying
    if(negate_)
    {
        refineCells = ~refineCells;
    }
    return refineCells;
}


Foam::bitSet
Foam::geometric::unrefinementPointCandidates() const
{
    Info
        << "Selecting unrefinement candidate points based on "
        << surftype_<< endl;
    
    bitSet unrefinePoints(mesh().nPoints(),false);

    List<volumeType> insideGeom;
    surfPtr_->getVolumeType(mesh().cellCentres(), insideGeom);

    // get cellLevels
    const labelIOList& curCellLevel =
        mesh().lookupObject<labelIOList>("cellLevel");

    forAll (mesh().points(), pI)
    {
        // Get point value
        for (const auto celli : mesh().pointCells()[pI])
        {
            bool unRefPoint = true;
            // all cells need to fullfill this criteria
            if
            (
                !(insideGeom[celli] == volumeType::INSIDE
                || curCellLevel[celli] > maxCellLevel_)
                || curCellLevel[celli] <= minCellLevel_
            )
            {
                unRefPoint = false;
                break;
            }
            
            if(unRefPoint)
            {
                unrefinePoints.set(pI);
            }
        }
    }

    syncTools::syncPointList(mesh(), unrefinePoints, minEqOp<unsigned int>(), 0);

    const label nCandidates = returnReduce(unrefinePoints.count(), sumOp<label>());

    Info
        << "-> "<< nCandidates << " split points selected"<< endl;

    // Print out some information
    // Info<< "Selection algorithm " << type() << " selected "
    //     << returnReduce(unrefinePoints.toc().size(), sumOp<label>())
    //     << " points as unrefinement candidates."
    //     << endl;

    if(negate_)
    {
        unrefinePoints = ~unrefinePoints;
    }

    return unrefinePoints;
}


// ************************************************************************* //
