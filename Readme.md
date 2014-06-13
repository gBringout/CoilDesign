# Coil Design #
The repository contains many Matlab script used to design coils.
The available technique is:
* [BEM using stream function.](#bem-stream-function-technique)

## BEM stream function technique ##
Here, BEM formulation is used along with stream function to formulate the coil design as an optimization. 3 Basic examples are presented to make coils for [Magnetic Particle Imaging](http://en.wikipedia.org/wiki/Magnetic_particle_imaging) scanner according to an upcoming publication. The examples are:
+ A circular quadrupole,
+ A circular drive coil,
+ A planar drive coil.

### Instalation ###
To work, this script required 3 externals scripts:
* [SphericalHamornics](https://github.com/gBringout/SphericalHarmonics)
* [regu](http://www.imm.dtu.dk/~pcha/Regutools/) package from C Hansen for the Tikhonov implementation
* [OPTI TOOLBOX](http://www.i2c2.aut.ac.nz/Wiki/OPTI/)  from J. Currie for the QP implementation 

To start, you have to adapt streamBEM/makeCoil6.m to your installation, and then run it.
The definition of the wanted coils is done in a separate files. One files per coil.

# Contributing #
You are welcome to contribute to this repository. Please read the associated license.
