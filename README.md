# interPlicFoam

## Introduction

```interPlicFoam``` is a two-phase flow solver derived from interFoam for two incompressible, isothermal and immiscible fluids using the PLIC-VOF method. The fraction advecting follows the algorithms developed in [isoAdvector](https://github.com/isoAdvector/isoAdvector). The interface inside a mixed polygonal/polyhedral cell is approximated with an orientated plane which defined as:

![](http://latex.codecogs.com/gif.latex?\\vec{n}\cdot\vec{X}+D_0=0),

where ![](http://latex.codecogs.com/gif.latex?\\vec{n}) is the unit orientation vector, ![](http://latex.codecogs.com/gif.latex?\\vec{X}) abritry point on the plane and ![](http://latex.codecogs.com/gif.latex?\\D_0) the signed distance.
