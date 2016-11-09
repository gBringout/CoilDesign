# Coil Design #
The repository contains many Matlab script used to design coils.
The available techniques are:
* [Wirepath,](#wirepath)
* [BEM using stream function,](#bem-stream-function-technique)
* [Iron core shaping using Multipole Expansion.](#iron-core-shaping)

## Installation ##
To work, this script required 3 externals scripts:
* [SphericalHamornics](https://github.com/gBringout/SphericalHarmonics)
* [regu](http://www.imm.dtu.dk/~pcha/Regutools/) package from C Hansen for the Tikhonov implementation
* [YALMIP TOOLBOX](http://users.isy.liu.se/johanl/yalmip/)  from Johan Löfberg for the QP implementation 

To start, you have to adapt streamBEM/makeCoil6.m to your installation, and then run it.
The definition of the wanted coils is done in a separate files. One files per coil.

# Implemented design techniques #
## Wirepath ##
The easiest approach to approximate coil centroids: simply defining loops of current at given positions, to make solenoids or gradient coils for example.

## BEM stream function technique ##
<div align="center">
        <img width="45%" src="/streamBEM/examples/DriveCurrentAndWire.jpg" alt="stream function & wire" title="Stream function & wire"</img>
        <img height="0" width="8px">
        <img width="35%" src="/streamBEM/examples/Proto.png" alt="Protype" title="Protype"></img>
</div>

Here, BEM formulation is used along with stream function to formulate the coil design as an optimization. 3 Basic examples are presented to make coils for [Magnetic Particle Imaging](http://en.wikipedia.org/wiki/Magnetic_particle_imaging) scanner according to [this publication](http://gael-bringout.com/public/Bringout%202014%20-%20Coil%20Design%20for%20Magnetic%20Particle%20Imaging%20Application%20for%20a%20Preclinical%20Scanner.pdf). The examples are:
+ A circular quadrupole,
+ A circular drive coil,
+ A planar drive coil,
+ A y-gradient coil for MRI,
+ A drive coil using a low-density mesh.

This has been extended to evaluate the effect of induced current in surfaces close to the coil on the coil efficiency and field topology. An example for a 
+ shielded y-drive coil for a MPI scanner

is provided.

Please refer yourself to [this publication](http://www.gael-bringout.com/public/Bringout%202014%20-%20Performance%20of%20shielded%20electromagnet%20-%20evaluation%20under%20low.pdf) and [my thesis](http://gael-bringout.com/public/DissGaelBringout_Online_Compressed.pdf) for more details.

## Iron core shaping ##
Using a multipole expansion techniques, the shape of iron core magnet can be designed. This has been used in [this conference contribution](http://www.gael-bringout.com/public/Bringout%202015%20-%20Performance%20and%20safety%20evaluation%20of%20a%20human%20sized%20FFL%20imager%20concept.pdf) to design a human sized MPI FFL scanner.
Please refer yourself to [my thesis](http://gael-bringout.com/public/DissGaelBringout_Online_Compressed.pdf) for more details.

# Contributing #
You are welcome to contribute to this repository. Please read the associated license.
To thank you for your (coming) help, please have this pig coil:

![Alt text](/streamBEM/examples/smallPigCoilStream.gif?raw=true "Stream function on a pig")
![Alt text](/streamBEM/examples/smallPigCoilWire.gif?raw=true "Stream function on a pig")
