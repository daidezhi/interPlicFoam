# interPlicFoam

## Introduction

```interPlicFoam``` is a two-phase flow solver derived from interFoam for two incompressible, isothermal and immiscible fluids using the PLIC-VOF method. The fraction advecting follows the algorithms developed in [isoAdvector](https://github.com/isoAdvector/isoAdvector). The interface inside a mixed polygonal/polyhedral cell is approximated with an orientated plane which defined as:

![](http://latex.codecogs.com/gif.latex?\\vec{n}\cdot\vec{X}+D_0=0),

where ![](http://latex.codecogs.com/gif.latex?\\vec{n}) is the unit orientation vector, ![](http://latex.codecogs.com/gif.latex?\\vec{X}) abritry point on the plane and ![](http://latex.codecogs.com/gif.latex?D_0) the signed distance.

In the interface reconstruction step, the unit orientation vector ![](http://latex.codecogs.com/gif.latex?\\vec{n}) is evaluated by using the gradient of ```alpha``` field. Then the signed distance ![](http://latex.codecogs.com/gif.latex?D_0) is computed by employing the Standard algorithm developed in:

[Dezhi Dai and Albert Y. Tong (2019). "Analytical interface reconstruction algorithms in the PLIC‐VOF method for 3D polyhedral unstructured meshes." International Journal for Numerical Methods in Fluids.](https://doi.org/10.1002/fld.4750)
