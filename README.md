# interPlicFoam

## Introduction

```interPlicFoam``` is a two-phase flow solver derived from interFoam for two incompressible, isothermal and immiscible fluids using the PLIC-VOF method. The fraction advecting follows the algorithms developed in [isoAdvector](https://github.com/isoAdvector/isoAdvector). The interface inside a mixed polygonal/polyhedral cell is approximated with an orientated plane which defined as:

![](http://latex.codecogs.com/gif.latex?\\vec{n}\cdot\vec{X}+D_0=0),

where ![](http://latex.codecogs.com/gif.latex?\\vec{n}) is the unit orientation vector, ![](http://latex.codecogs.com/gif.latex?\\vec{X}) abritry point on the plane and ![](http://latex.codecogs.com/gif.latex?D_0) the signed distance.

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


## Demos

Two dam-breaking tutorial cases are available in http://dx.doi.org/10.17632/wm5w5g3kzt.1 (```damBreak.tar.gz``` and ```damBreakKleefsman.tar.gz```).

>**It should be noted that ```damBreakKleefsmanFull.tar.gz``` is the numerical simulation results associated with the paper https://doi.org/10.1002/fld.4750.**

>**The commands are integrated into the ```Allrun``` and ```Allrun-parallel``` script for each case.**

### ```damBreak``` (2D)



### ```damBreakKleefsman``` (3D)

#[Imgur](https://i.imgur.com/Z8kGlzk.gif)

<img src="https://i.imgur.com/Z8kGlzk.gif">

## Change Log

### 06/16/2019

* Initial release


## Contributors

* Dezhi Dai, UT Arlington, dezhi.dai@mavs.uta.edu (Developer)