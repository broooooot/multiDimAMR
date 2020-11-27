# dynamicAMR

## Description

dynamic meshing with load balancing for hexahedral meshes in 3D and 2D

Library is based on:

Rettenmaier, Daniel, et al. "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM." SoftwareX 10 (2019): 100317.

link:
https://www.sciencedirect.com/science/article/pii/S2352711018301699

port to the OpenFOAM+ version v2006

refinement selection algoritm is based on foam extended 4.1

## Changes wrt Henning Scheuflers Version

### general
- adapted the code to allow compilation with OpenFOAM v2006
- replaced deprecated operators
- load balancing does not work at the moment (don't know why)

### new features
- allow refinement based on gradient or curl of a field via extra keyword in dynamicMeshDict

### removed 
- folder and files which are not needed, focus is on adaptive mesh refinement (tutorials might follow at a later point)
- removed possibility for XOR combination of refinement criteria (as I see no point in having that)

## Getting Started

### Prerequisites

Install OpenFOAM v2006:

```
https://www.openfoam.com/download/release-history.php
```

### Installing
Clone and compile this library
```
 git clone https://github.com/broooooot/multiDimAMR
 cd multiDimAMR
 ./Allwmake
```
### Usage

Add this line to the controlDict
```
libs
(
   "libdynamicLoadBalanceFvMesh.so"	
);
```

and add dynamicMeshDict inside the constant/ folder. You'll find
a template in the dict/ folder. 

The balancing process does not work at the moment so you don't need
the balanceParDict.

## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments



