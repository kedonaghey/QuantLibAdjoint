
# QuantLibAdjoint

QuantLibAdjoint implements adjoint algorithmic differentiation (AAD)
in QuantLib using
[TapeScript](http://tapescript.org),  an open source (Apache license) 
C++ library with support for scalar and vector AAD.

QuantLibAdjoint repository is a fork of the
master QuantLib repository located at
[github.com/lballabio/QuantLib](http://github.com/lballabio/QuantLib).

QuantLib, QuantLibAdjoint, and TapeScript are Non-Copylefted Free Software and OSI Certified Open Source Software.

The QuantLib project (<http://quantlib.org>) is aimed at providing a
comprehensive software framework for quantitative finance. QuantLib is
a free/open-source library for modeling, trading, and risk management
in real-life.

Bugs can be reported  at
<https://github.com/compatibl/QuantLibAdjoint/issues>; if you have a patch
for a bugfix or new feature, you can open a pull request instead.

## TapeScript

[TapeScript](http://tapescript.org) is an open source library for adjoint algorithmic differentiation
(AAD) developed and maintained by CompatibL. It can be downloaded from
[github.com/compatibl/tapescript](http://github.com/compatibl/tapescript) and used free of charge in academic or commercial applications.

TapeScript supports vector AAD (tape compression), an approach in which
each slot of the calculation record (AAD tape) can store not only a
single double number, but also an entire array of values. Vector AAD
can lead to performance gain of several orders of magnitude due to the reduction
of tape size.

### TapeScript features:

* Scalar AAD
* Vector AAD (tape compression)
* APIs for C++, C#, and Java
* Complex numbers
* Works with Boost and QuantLib
* Multithreading support

## About CompatibL

CompatibL offers turnkey solutions for XVA and regulatory capital
as well as custom development, integration, and consultancy.

Check out TapeLib, CompatibL's product suite for AAD that incorporates:

* A C++ library extending TapeScript with features specific to quantitative finance
* An application platform for interactive AAD

### TapeLib features:

* Document database preserving AAD data
* Tape database
* Finance-specific atomics including adjointable AMC
* User defined atomics
* Tape cutting and splicing
* Parallel tape execution
* Specialized gate checking API
* Excel addin, desktop and client
