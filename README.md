# interPlicFoam

## Introduction

```interPlicFoam``` is a two-phase flow solver derived from interFoam for two incompressible, isothermal and immiscible fluids using the PLIC-VOF method. The fraction advecting follows the algorithms developed in [isoAdvector](https://github.com/isoAdvector/isoAdvector). The interface inside a polygonal/polyhedral cell is approximated as an orientated plane which defined as:

![](http://latex.codecogs.com/gif.latex?\\\vec{n}\cdot\vec{X}+D_0=0)
