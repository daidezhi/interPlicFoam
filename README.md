# interPlicFoam

---

:warning: **IMPORTANT** :warning: The latest `interPlicFoam` is available at [geometricVofExt](https://github.com/daidezhi/geometricVofExt).

---

## Introduction

```interPlicFoam``` is a two-phase flow solver derived from interFoam for two incompressible, isothermal and immiscible fluids using the PLIC-VOF method. The fraction advecting follows the algorithms developed in [isoAdvector](https://github.com/isoAdvector/isoAdvector). The interface inside a mixed polygonal/polyhedral cell is approximated with an orientated plane which defined as:

![](http://latex.codecogs.com/gif.latex?\\vec{n}\cdot\vec{X}+D_0=0),

where ![](http://latex.codecogs.com/gif.latex?\\vec{n}) is the unit orientation vector, ![](http://latex.codecogs.com/gif.latex?\\vec{X}) arbitrary point on the plane and ![](http://latex.codecogs.com/gif.latex?D_0) the signed distance.

In the interface reconstruction step, the unit orientation vector ![](http://latex.codecogs.com/gif.latex?\\vec{n}) is evaluated by using the gradient of ```alpha``` field. Then the signed distance ![](http://latex.codecogs.com/gif.latex?D_0) is computed by employing the Standard algorithm developed in:

[Dezhi Dai and Albert Y. Tong (2019). "Analytical interface reconstruction algorithms in the PLIC‐VOF method for 3D polyhedral unstructured meshes." International Journal for Numerical Methods in Fluids.](https://doi.org/10.1002/fld.4750)


## Compatibility

The source code of ```interPlicFoam``` is developed and maintained for the latest OpenFOAM-plus release (currently ```OpenFOAM-v1812```). A script for generating code for other releases and old versions will be provided in the near future.

*Please notify me via the email address below if you found any errors or bugs, and I will fix them as soon as possible.*


## Installation

**Before compiling ```interPlicFoam```, make sure that the OpenFOAM environment has been set properly.**

1. Download ```interPlicFoam``` from this page

2. Build ```libplicVofSolving.so```
```bash
cd plic
wmake libso
```

3. Build ```interPlicFoam```
```bash
cd ..
wmake
```

*All of the compiling commands above have been integrated into ```Allwmake``` script.*


## Usage

The usage of ```interPlicFoam``` is similar with ```interFoam```, however, two extra settings should be considered for ```interPlicFoam```:

* ```gradSchemes``` of ```alpha.water``` (```alpha1```) field (*used for computing interface orientation vectors*), e.g.,
```c++
gradSchemes
{
    default             Gauss linear;
    gradAlpha           Gauss pointLinear;
}
```

* ```solvers``` of ```alpha.water``` (```alpha1```) field, i.e.,
```c++
"alpha.water.*"
{
    surfCellTol         1e-8;   // Tolerance for marking mixed cells
    nAlphaBounds        3;      // Number of alpha bounding steps
    snapTol             1e-8;   // Tolerance of fraction value snapping
    clip                true;   // Switch of fraction value clipping
    smoothedAlphaGrad   false;  // Switch of smoothed alpha gradient

    writePlicFaces      true;   // Switch of reconstructed interface outputting

    nAlphaSubCycles     1;      // Number of alpha sub-cycles

    // Note: cAlpha is not used by interPlicFoam but must
    // be specified because interfacePropertes object
    // reads it during construction.
    cAlpha              1;
}
```


## Demos

Two dam-breaking tutorial cases are available in http://dx.doi.org/10.17632/wm5w5g3kzt.1 (```damBreak.tar.gz``` and ```damBreakKleefsman.tar.gz```).

>**It should be noted that ```damBreakKleefsmanFull.tar.gz``` is the numerical simulation results associated with the paper https://doi.org/10.1002/fld.4750.**

>**The commands are integrated into the ```Allrun``` and ```Allrun-parallel``` scripts for each case.**

### ```damBreak``` (2D)

1. Geometry and mesh

<img src="https://i.imgur.com/FQTp96O.png" width="400"><img src="https://i.imgur.com/jkn2VvK.png" width="400">

2. Time history of the interface profiles

<img src="https://imgur.com/lKXr6KS.gif" width="400">

3. Reconstructed interface output

<img src="https://imgur.com/b5sXIPJ.png" width="800">

<img src="https://imgur.com/NShuSop.png" width="800">


### ```damBreakKleefsman``` (3D)

1. Geometry and mesh

<img src="https://imgur.com/ljBgzWJ.png" width="500">

<img src="https://imgur.com/haWjXMc.png" width="500">

2. Time history of the interface profiles

<img src="https://i.imgur.com/8ibnl2u.gif">

3. Reconstructed interface output

<img src="https://imgur.com/6z9RrUB.png">


## Change Log

### 06/16/2019

* Initial release


## Contributors

* Dezhi Dai, UT Arlington, dezhi.dai@mavs.uta.edu (Developer)
