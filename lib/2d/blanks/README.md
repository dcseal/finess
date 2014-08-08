# Default (blank) templates for 2D FINESS

This is a collection of potential templates in the [2d](../README.md) branch
of [FINESS](../../../README.md) that can be used for each application.  It is
advisable that the user copies an existing application in place of copying or
modifying these.

The code assumes that each of these function have been defined in order to
compile.  For example, please visit the application directory.

A typical application is required to define the following

* AuxFunc.cpp
* FluxFunc.cpp
* ProjectLeftEig.cpp
* ProjectRightEig.cpp
* QinitFunc.cpp
* SetBndValues.cpp
* SetWaveSpd.cpp
* SourceTermFunc.cpp

Please visit the source code for a description of these functions.
