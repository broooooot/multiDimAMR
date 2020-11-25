# dynamicAMR

## Descrition

dynamic meshing with load balancing for hexahedral meshes in 3D and 2D

Library is based on:

Rettenmaier, Daniel, et al. "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM." SoftwareX 10 (2019): 100317.

link:
https://www.sciencedirect.com/science/article/pii/S2352711018301699

port to the OpenFOAM+ version v2006

refinement selection algoritm is based on foam extended 4.1

## Changes to Henning Scheuflers Version

- adapted the code to allow compilation with OpenFOAM v2006
- replaced deprecated operators

## Getting Started

### Prerequisites

Install OpenFOAM v2006:

```
https://www.openfoam.com/download/release-history.php
```

### Installing
clone and compile this library
```
 git clone https://github.com/broooooot/multiDimAMR
 cd multiDimAMR
 ./Allwmake
```
### Usage

add to controlDict:
```
libs
(
   "libdynamicLoadBalanceFvMesh.so"	
);
```

system/decomposeParDict
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2006                                   |
|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;

method          scotch;

constraints
{
    refinementHistoryMultiDim
    {
        //- Decompose cells such that all cell originating from single cell
        //  end up on same processor
        type    refinementHistoryMultiDim;
    }
}
```

system/balanceParDict:
```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2006                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      balanceParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains 4;

method          ptscotch; //clsutered //scotch

constraints
{
    refinementHistoryMultiDim
    {
        //- Decompose cells such that all cell originating from single cell
        //  end up on same processor
        type    refinementHistoryMultiDim;
    }
}
// ************************************************************************* //
```

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments



