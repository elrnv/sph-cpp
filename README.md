Multi Material SPH Framework
============================

This framework is developed to promote research of Smoothed Particle Hydrodynamic
(SPH) flows. It was originally a course project demonstrating two SPH
implementations from:

  * "Particle-Based Fluid Simulation for Interactive Applications" by M. Muller, D. Charypar and M. Gross
  * "Weakly Compressible SPH for Free Surface Flows" by M. Becker, M. Teschner

In addition, the "Implicit Incompressible SPH" method presented by M. Ihmsen et.
al. is partially implemented.

The goal of this project is to establish a framework which is easily extendable
to develop multi-material and multi-phase SPH schemes. Currently, it is possible
to have SPH fluids from both implemented methods to interact.

Throughout the code, the method by Muller et al is referenced by *MCG03*, method
by Beckner and Teschner is referenced by *BT07*, and the IISPH method is
referenced by *ICS13*. Note, however, these methods aren't implemented exactly
as they were originally concieved.


Overview
========

This project uses *boost*, *Qt 5.3*, *Assimp 3.1.1*, and *libconfig 1.4.9*
libraries. It is by no means complete, but I believe it is a good starting point
for someone who is studying SPH.

The application reads config files stored in the ``data/`` directory named
``scene#.cfg``, which are loaded with the keyboard shortcuts ``0-9``.


Setup
=====

A small environment setup script is provided in ``.setup.sh`` for unix-like
enviroments with a bash shell. Run 
````
$ source .setup.sh
````
to get a shortcut to built executable binaries as well as the ``bin/`` directory
in your ``PATH``, and the environment variable ``SPH`` which tells the
executable where to search for input and output files.
If you don't have bash, you must define the SPH environment variable to
point to the root project directory manually in order to run the program.


Build
=====

This project relies on ``qmake`` to generate unix-style Makefiles. The standard
way to build and run the project is to run
````
$ qmake               # creates Makefiles for debug and release builds
$ make                # compiles a default release build
$ source .setup.sh    # sets the enviroment variable SPH
$ sph                 # runs the binary
````


Supported Platforms
===================

I only know the code builds on Mac OS X, probably Mavericks and Yosemite. The
project uses Qt so it shouldn't be difficult to port onto other platforms. Let
me know if you succeed.


Supplementary Material
======================

A demo of what this application should look like can be seen on
[youtube](http://youtu.be/zRr84SF6FZw).


Since this project was part of a course, I had to write a SIGGRAPH-style
[report](https://elrnv.com/static/cs888/cs888proj.pdf),
which may be helpful in guiding the interested reader through the fundamentals
of SPH and the necessary background needed to understand how the code works.

