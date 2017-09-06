.. # --------------------------------------------------------------------------
   #
   # Copyright
   #   Markus Wittmann, 2016-2017
   #   RRZE, University of Erlangen-Nuremberg, Germany
   #   markus.wittmann -at- fau.de or hpc -at- rrze.fau.de
   #
   #   Viktor Haag, 2016
   #   LSS, University of Erlangen-Nuremberg, Germany
   #
   #  This file is part of the Lattice Boltzmann Benchmark Kernels (LbmBenchKernels).
   #
   #  LbmBenchKernels is free software: you can redistribute it and/or modify
   #  it under the terms of the GNU General Public License as published by
   #  the Free Software Foundation, either version 3 of the License, or
   #  (at your option) any later version.
   #
   #  LbmBenchKernels is distributed in the hope that it will be useful,
   #  but WITHOUT ANY WARRANTY; without even the implied warranty of
   #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   #  GNU General Public License for more details.
   #
   #  You should have received a copy of the GNU General Public License
   #  along with LbmBenchKernels.  If not, see <http://www.gnu.org/licenses/>.
   #
   # --------------------------------------------------------------------------

.. title:: LBM Benchmark Kernels Documentation 


===================================
LBM Benchmark Kernels Documentation
===================================

.. sectnum::
.. contents::

Compilation
===========

The benchmark framework currently supports only Linux systems and the GCC and 
Intel compilers. Every other configuration probably requires adjustment inside
the code and the makefiles. Further some code might be platform or at least
POSIX specific.

The benchmark can be build via ``make`` from the ``src`` subdirectory. This will
generate one binary which hosts all implemented benchmark kernels. 

Binaries are located under the ``bin`` subdirectory and will have different names
depending on compiler and build configuration.

Debug and Verification
----------------------

:: 

  make

Running ``make`` without any arguments builds the debug version (BUILD=debug) of
the benchmark kernels, where no optimizations are performed,  line numbers and
debug symbols are included as well as ``DEBUG`` will be defined.  The resulting
binary will be found in the ``bin`` subdirectory and named
``lbmbenchk-linux-<compiler>-debug``.
 
Without any further specification the binary includes verification
(``VERIFICATION=on``), statistics (``STATISTICS``), and VTK output
(``VTK_OUTPUT=on``) enabled. 

Please note that the generated binary will therefore
exhibit a poor performance.

Benchmarking
------------

To generate a binary for benchmarking run make with ::

  make BENCHMARK=on BUILD=release

Here BUILD=release turns optimizations on and BENCHMARK=on disables
verfification, statistics, and VTK output.

Release and Verification
------------------------

Verification with the debug builds can be extremely slow. Hence verification
capabilities can be build with release builds: ::

  make BUILD=release

Compilers
---------

Currently only the GCC and Intel compiler under Linux are supported. Between
both configuration can be chosen via ``CONFIG=linux-gcc`` or
``CONFIG=linux-intel``.

Options Summary
---------------

Options that can be specified when building the framework with make:

============= ======================= ============ ==========================================================
name          values                  default      description
------------- ----------------------- ------------ ----------------------------------------------------------
TARCH         --                      --           Via TARCH the architecture the compiler generates code for can be overriden. The value depends on the chose compiler.
BENCHMARK     on, off                 off          If enabled, disables VERIFICATION, STATISTICS, VTK_OUTPUT.
BUILD         debug, release          debug        No optimization, debug symbols, DEBUG defined.
CONFIG        linux-gcc, linux-intel  linux-intel  Select GCC or Intel compiler. 
ISA           avx, sse                avx          Determines which ISA extension is used for macro definitions. This is *not* the architecture the compiler generates code for.
OPENMP        on, off                 on           OpenMP, i.\,e.\. threading support.
STATISTICS    on, off                 off          View statistics, like density etc, during simulation. 
VERIFICATION  on, off                 off          Turn verification on/off.
VTK_OUTPUT    on, off                 off          Enable/Disable VTK file output.
============= ======================= ============ ==========================================================

Invocation
==========

Running the binary will print among the GPL licence header a line like the following:
 
  LBM Benchmark Kernels 0.1, compiled Jul  5 2017 21:59:22, type: verification

if verfication was enabled during compilation or

  LBM Benchmark Kernels 0.1, compiled Jul  5 2017 21:59:22, type: benchmark

if verfication was disabled during compilation.

Command Line Parameters
-----------------------

Running the binary with ``-h`` list all available parameters: ::

  Usage:
  ./lbmbenchk -list
  ./lbmbenchk
      [-dims XxYyZ] [-geometry box|channel|pipe|blocks[-<block size>]] [-iterations <iterations>] [-lattice-dump-ascii]
      [-rho-in <density>] [-rho-out <density] [-omega <omega>] [-kernel <kernel>]
      [-periodic-x]
      [-t <number of threads>]
      [-pin core{,core}*]
      [-verify]
      -- <kernel specific parameters>

  -list           List available kernels.

  -dims XxYxZ     Specify geometry dimensions.

  -geometry blocks-<block size>
                  Geometetry with blocks of size <block size> regularily layout out.


If an option is specified multiple times the last one overrides previous ones.
This holds also true for ``-verify`` which sets geometry dimensions,
iterations, etc, which can afterward be override, e.g.: ::

  $ bin/lbmbenchk-linux-intel-release -verfiy -dims 32x32x32

Kernel specific parameters can be opatained via selecting the specific kernel
and passing ``-h`` as parameter: ::

  $ bin/lbmbenchk-linux-intel-release -kernel -- -h
  ...
  Kernel parameters:
  [-blk <n>] [-blk-[xyz] <n>]

  
A list of all available kernels can be obtained via ``-list``: ::

  $ ../bin/lbmbenchk-linux-gcc-debug -list
  Lattice Boltzmann Benchmark Kernels (LbmBenchKernels) Copyright (C) 2016, 2017 LSS, RRZE
  This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.
  This is free software, and you are welcome to redistribute it under certain conditions.

  LBM Benchmark Kernels 0.1, compiled Jul  5 2017 21:59:22, type: verification
  Available kernels to benchmark:
     list-aa-pv-soa
     list-aa-ria-soa
     list-aa-soa
     list-aa-aos
     list-pull-split-nt-1s-soa
     list-pull-split-nt-2s-soa
     list-push-soa
     list-push-aos
     list-pull-soa
     list-pull-aos
     push-soa
     push-aos
     pull-soa
     pull-aos
     blk-push-soa
     blk-push-aos
     blk-pull-soa
     blk-pull-aos


Benchmarking
============

Correct benchmarking is a nontrivial task. Whenever benchmark results should be
created make sure the binary was compiled with:  

- ``BENCHMARK=on`` and
- ``BUILD=release`` and 
- the correct ISA for macros is used, selected via ``ISA`` and
- use ``TARCH`` to specify the architecture the compiler generates code for.

During benchmarking pinning should be used via the ``-pin`` parameter. Running
a benchmark with 10 threads an pin them to the first 10 cores works like ::

  $ bin/lbmbenchk-linux-intel-release ... -t 10 -pin $(seq -s , 0 9)

Things the binary does nor check or controll:

- transparent huge pages: when allocating memory small 4 KiB pages might be
  replaced with larger ones. This is in general a good thing, but if this is
  really the case, depends on the system settings.

- CPU/core frequency: For reproducible results the frequency of all cores
  should be fixed.

- NUMA placement policy: The benchmark assumes a first touch policy, which
  means the memory will be placed at the NUMA domain the touching core is
  associated with. If a different policy is in place or the NUMA domain to be
  used is already full memory might be allocated in a remote domain. Accesses
  to remote domains typically have a higher latency and lower bandwidth.

- System load: interference with other application, espcially on desktop
  systems should be avoided.

- Padding: most kernels do not care about padding against cache or TLB
  thrashing. Even if the number of (fluid) nodes suggest everything is fine,
  through parallelization still problems might occur.

- CPU dispatcher function: the compiler might add different versions of a
  function for different ISA extensions. Make sure the code you might think is
  executed is actually the code which is executed.


Acknowledgements
================

This work was funded by BMBF, grant no. 01IH15003A (project SKAMPY).

This work was funded by KONWHIR project OMI4PAPS.



.. |datetime| date:: %Y-%m-%d %H:%M

Document was generated at |datetime|.

