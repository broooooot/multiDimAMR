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
    Vuko Vukcevic, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "composedAdaptCriteria.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(composedAdaptCriteria, 0);
addToRunTimeSelectionTable
(
    adaptCriteria,
    composedAdaptCriteria,
    dictionary
);

}

const Foam::Enum
<
    Foam::composedAdaptCriteria::operationType
>
Foam::composedAdaptCriteria::operationTypeNames_
({

    { operationType::opOr, "or" },
    { operationType::opAnd, "and" },

});

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::composedAdaptCriteria::composedAdaptCriteria
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    adaptCriteria(mesh, dict),
    operation_(operationTypeNames_.get("operation", dict)),
    criteriaSelections_()
{
    // Read basic refinement selections
    PtrList<entry> criteriaEntries
    (
        dict.lookup("criteria")
    );

    Info << "Using multiple criteria for mesh refinement"<< nl<< "( ";
    forAll(criteriaEntries, Index)
    {
        Info << criteriaEntries[Index].name() << " ";
    }
    Info<< ")"<< endl;

    criteriaSelections_.setSize(criteriaEntries.size());

    forAll (criteriaSelections_, brsI)
    {
        criteriaSelections_.set
        (
            brsI,
            adaptCriteria::New
            (
                mesh,
                criteriaEntries[brsI].dict()
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::composedAdaptCriteria::~composedAdaptCriteria()
{}
            
// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //
void Foam::composedAdaptCriteria::reReadDictionary
(
    const dictionary& reReadDict
)
{       
    PtrList<entry> criteriaEntries
    (
        reReadDict.lookup("criteria")
    );

    forAll(criteriaSelections_, brsI)
    {
        criteriaSelections_[brsI].reReadDictionary(criteriaEntries[brsI].dict());
    }
}    

Foam::bitSet
Foam::composedAdaptCriteria::refinementCellCandidates() const
{
    bitSet refineCells(mesh().nCells(),false);

    // Loop through all base refinement selections
    forAll(criteriaSelections_, brsI)
    {
        // Get refinement candidates from this base selection algorithm. Note:
        // list is transferred
        const bitSet curRefCandidates
        (
            criteriaSelections_[brsI].refinementCellCandidates()
        );

        if(brsI == 0)
        {
            refineCells = curRefCandidates;
        }
        else
        {
            switch (operation_)
            {
                // 11001 | 01101 = 11101
                case opOr:
                {
                    refineCells = refineCells | curRefCandidates;
                    break;
                }
                // 11001 & 01101 = 01001
                case opAnd:
                {
                    refineCells = refineCells & curRefCandidates;
                    break;
                }
                default:
                {
                    FatalErrorInFunction << "operation not available" << endl;
                    break;
                }
            }
        }        
    }

    Info
        << "Selected " << returnReduce(refineCells.count() , sumOp<label>())
        << " candidate cells which";

        switch (operation_)
        {
            case opOr:
            {
                Info<< " meet at least one criterion."<< endl;
                break;
            }
            case opAnd:
            {
                Info<< " meet all criteria."<< endl;
                break;
            }
        }

    // // Create storage for collection of final cell candidates. Assume that
    // // one fifth of the cells will be marked to prevent excessive resizing
    // DynamicList<label> refinementCandidates(mesh().nCells()/5);

    // // Get number of active selection algorithms
    // const label nBaseSelections = criteriaSelections_.size();

    // // Loop through all cells and collect final refinement candidates
    // forAll (nSelections, cellI)
    // {
    //     if (nSelections[cellI] > 0) //nBaseSelections)
    //     {
    //         // Cell has been marked by all selection algorithms, append it
    //         refinementCandidates.append(cellI);
    //     }
    // }

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
Foam::composedAdaptCriteria::unrefinementPointCandidates() const
{
    bitSet unrefinePoints(mesh().nPoints(),true);

    // Loop through all base refinement selections
    forAll (criteriaSelections_, brsI)
    {
        // Get unrefinement candidates from this base selection algorithm. Note:
        // list is transferred
        const bitSet curUnRefCandidates
        (
            criteriaSelections_[brsI].unrefinementPointCandidates()
        );

        // Increment the number of selections for selected cells
        // forAll (curRefCandidates, i)
        // {
        //     ++nSelections[curRefCandidates[i]];
        // }
        if(brsI == 0)
        {
            unrefinePoints = curUnRefCandidates;
        }
        else
        {
            switch (operation_)
            {
                // refine cell if it matches at least ONE criteria
                // so unrefine only if it matches ALL unrefinement criteria
                case opOr:
                {
                    unrefinePoints = unrefinePoints & curUnRefCandidates;
                    break;
                }
                // refine cell if it matches ALL criteria
                // so unrefine if it matches even ONE unrefinement criteria
                case opAnd:
                {
                    unrefinePoints = unrefinePoints | curUnRefCandidates;
                    break;
                }
                default:
                {
                    FatalErrorInFunction << "operation not available" << endl;
                    break;
                }
            }
        }
    }

    Info
        << "Selected " << returnReduce(unrefinePoints.count() , sumOp<label>())
        << " candidate points which";

        switch (operation_)
        {
            case opOr:
            {
                Info<< " meet all criteria."<< endl;
                break;
            }
            case opAnd:
            {
                Info<< " meet at least one criterion."<< endl;
                break;
            }
        }

    // // Create storage for collection of final point candidates. Assume that one
    // // tenth of the points will be marked to prevent excessive resizing
    // DynamicList<label> unrefinementCandidates(mesh().nPoints()/10);

    // // Get number of active selection algorithms
    // const label nBaseSelections = criteriaSelections_.size();

    // // Loop through all points and collect final unrefinement candidates
    // forAll (nSelections, pointI)
    // {
    //     if (nSelections[pointI] == nBaseSelections)
    //     {
    //         // Point has been marked by all selection algorithms, append it
    //         unrefinementCandidates.append(pointI);
    //     }
    // }

    // Print out some information
    // Info<< "Selection algorithm " << type() << " selected "
    //     << returnReduce(unrefinePoints.toc().size(), sumOp<label>())
    //     << " points as unrefinement candidates."
    //     << endl;

    // Return the list in the Xfer container to prevent copying
    if(negate_)
    {
        unrefinePoints = ~unrefinePoints;
    }

    return unrefinePoints;
}


// ************************************************************************* //
