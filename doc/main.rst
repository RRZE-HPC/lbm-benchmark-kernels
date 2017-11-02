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

  make BUILD=debug BENCHMARK=off

Running ``make`` with ``BUILD=debug`` builds the debug version of
the benchmark kernels, where no optimizations are performed,  line numbers and
debug symbols are included as well as ``DEBUG`` will be defined.  The resulting
binary will be found in the ``bin`` subdirectory and named
``lbmbenchk-linux-<compiler>-debug``.
 
Specifying ``BENCHMARK=off`` turns on verification
(``VERIFICATION=on``), statistics (``STATISTICS=on``), and VTK output
(``VTK_OUTPUT=on``) enabled. 

Please note that the generated binary will therefore
exhibit a poor performance.

Benchmarking
------------

To generate a binary for benchmarking run make with ::

  make 

As default ``BENCHMARK=on`` and ``BUILD=release`` is set, where
BUILD=release turns optimizations on and ``BENCHMARK=on`` disables
verfification, statistics, and VTK output.

Release and Verification
------------------------

Verification with the debug builds can be extremely slow. Hence verification
capabilities can be build with release builds: ::

  make BENCHMARK=off 

Compilers
---------

Currently only the GCC and Intel compiler under Linux are supported. Between
both configuration can be chosen via ``CONFIG=linux-gcc`` or
``CONFIG=linux-intel``.


Cleaning
--------

For each configuration and build (debug/release) a subdirectory under the
``src/obj`` directory is created where the dependency and object files are
stored.
With ::

  make CONFIG=... BUILD=... clean

a specific combination is select and cleaned, whereas with ::

  make clean-all

all object and dependency files are deleted.


Options Summary
---------------

Options that can be specified when building the framework with make:

============= ======================= ============ ==========================================================
name          values                  default      description
------------- ----------------------- ------------ ----------------------------------------------------------
BENCHMARK     on, off                 on           If enabled, disables VERIFICATION, STATISTICS, VTK_OUTPUT. If disabled enables the three former options.
BUILD         debug, release          release      No optimization, debug symbols, DEBUG defined.
CONFIG        linux-gcc, linux-intel  linux-intel  Select GCC or Intel compiler. 
ISA           avx, sse                avx          Determines which ISA extension is used for macro definitions. This is *not* the architecture the compiler generates code for.
OPENMP        on, off                 on           OpenMP, i.\,e.\. threading support.
STATISTICS    on, off                 off          View statistics, like density etc, during simulation. 
TARCH         --                      --           Via TARCH the architecture the compiler generates code for can be overridden. The value depends on the chosen compiler.
VERIFICATION  on, off                 off          Turn verification on/off.
VTK_OUTPUT    on, off                 off          Enable/Disable VTK file output.
============= ======================= ============ ==========================================================

Invocation
==========

Running the binary will print among the GPL licence header a line like the following: ::
 
  LBM Benchmark Kernels 0.1, compiled Jul  5 2017 21:59:22, type: verification

if verfication was enabled during compilation or ::

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

  $ bin/lbmbenchk-linux-intel-release -kernel kernel-name -- -h
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

Kernels
-------

The following list shortly describes available kernels:

- push-soa/push-aos/pull-soa/pull-aos:
  Unoptimized kernels (but stream/collide are already fused) using two grids as
  source and destination. Implement push/pull semantics as well structure of
  arrays (soa) or array of structures (aos) layout.

- blk-push-soa/blk-push-aos/blk-pull-soa/blk-pull-aos:
  The same as the unoptimized kernels without the blk prefix, except that they support
  spatial blocking, i.e. loop blocking of the three loops used to iterate over
  the lattice. Here manual work sharing for OpenMP is used.

- list-push-soa/list-push-aos/list-pull-soa/list-pull-aos:
  The same as the unoptimized kernels without the list prefix, but for indirect addressing.
  Here only a 1D vector of is used to store the fluid nodes, omitting the
  obstacles. An adjacency list is used to recover the neighborhood associations.

- list-pull-split-nt-1s-soa/list-pull-split-nt-2s-soa:
  Optimized variant of list-pull-soa. Chunks of the lattice are processed as
  once. Postcollision values are written back via nontemporal stores in 18 (1s)
  or 9 (2s) loops.

- list-aa-aos/list-aa-soa:
  Unoptimized implementation of the AA pattern for the 1D vector with adjacency
  list. Supported are array of structures (aos) and structure of arrays (soa)
  data layout is supported.

- list-aa-ria-soa:
  Implementation of AA pattern with intrinsics for the 1D vector with adjacency
  list. Furthermore it contains a vectorized even time step and run length
  coding to reduce the loop balance of the odd time step.

- list-aa-pv-soa:
  All optimizations of list-aa-ria-soa. Additional with partial vectorization
  of the odd time step.


Note that all array of structures (aos) kernels might require blocking
(depending on the domain size) to reach the performance of their structure of
arrays (soa) counter parts.

The following table summarizes the properties of the kernels. Here **D** means
direct addressing, i.e. full array, **I** means indirect addressing, i.e. 1D
vector with adjacency list, **x** means supported, whereas **--** means unsupported.
The loop balance B_l is computed for D3Q19 model with double precision floating
point for PDFs (8 byte) and 4 byte integers for the index (adjacency list).
As list-aa-ria-soa and list-aa-pv-soa support run length coding their effective
loop balance depends on the geometry. The effective loop balance is printed
during each run.


====================== =========== =========== ===== ======== ======== ============
kernel name            prop. step  data layout addr. parallel blocking B_l [B/FLUP]
====================== =========== =========== ===== ======== ======== ============
push-soa               OS          SoA         D     x         --      456
push-aos               OS          AoS         D     x         --      456
pull-soa               OS          SoA         D     x         --      456
pull-aos               OS          AoS         D     x         --      456
blk-push-soa           OS          SoA         D     x         x       456
blk-push-aos           OS          AoS         D     x         x       456
blk-pull-soa           OS          SoA         D     x         x       456
blk-pull-aos           OS          AoS         D     x         x       456
list-push-soa          OS          SoA         I     x         x       528
list-push-aos          OS          AoS         I     x         x       528
list-pull-soa          OS          SoA         I     x         x       528
list-pull-aos          OS          AoS         I     x         x       528
list-pull-split-nt-1s  OS          SoA         I     x         x       376
list-pull-split-nt-2s  OS          SoA         I     x         x       376
list-aa-soa            AA          SoA         I     x         x       340
list-aa-aos            AA          AoS         I     x         x       340
list-aa-ria-soa        AA          SoA         I     x         x       304-342
list-aa-pv-soa         AA          SoA         I     x         x       304-342
====================== =========== =========== ===== ======== ======== ============

Benchmarking
============

Correct benchmarking is a nontrivial task. Whenever benchmark results should be
created make sure the binary was compiled with:  

- ``BENCHMARK=on`` (default if not overriden) and
- ``BUILD=release`` (default if not overriden) and 
- the correct ISA for macros is used, selected via ``ISA`` and
- use ``TARCH`` to specify the architecture the compiler generates code for.

During benchmarking pinning should be used via the ``-pin`` parameter. Running
a benchmark with 10 threads an pin them to the first 10 cores works like ::

  $ bin/lbmbenchk-linux-intel-release ... -t 10 -pin $(seq -s , 0 9)

Things the binary does nor check or controll:

- transparent huge pages: when allocating memory small 4 KiB pages might be
  replaced with larger ones. This is in general a good thing, but if this is
  really the case, depends on the system settings (check e.g. the status of
  ``/sys/kernel/mm/transparent_hugepage/enabled``).
  Currently ``madvise(MADV_HUGEPAGE)`` is used for allocations which are aligned to
  a 4 KiB page, which should be the case for the lattices. 
  This should result in huge pages except THP is disabled on the machine.
  (NOTE: madvise() is used if ``HAVE_HUGE_PAGES`` is defined, which is currently
  hard coded defined in ``Memory.c``).

- CPU/core frequency: For reproducible results the frequency of all cores
  should be fixed.

- NUMA placement policy: The benchmark assumes a first touch policy, which
  means the memory will be placed at the NUMA domain the touching core is
  associated with. If a different policy is in place or the NUMA domain to be
  used is already full memory might be allocated in a remote domain. Accesses
  to remote domains typically have a higher latency and lower bandwidth.

- System load: interference with other application, espcially on desktop
  systems should be avoided.

- Padding: For SoA based kernels the number of (fluid) nodes is automatically
  adjusted so that no cache or TLB thrashing should occur. The parameters are
  optimized for current Intel based systems. For more details look into the
  padding section.

- CPU dispatcher function: the compiler might add different versions of a
  function for different ISA extensions. Make sure the code you might think is
  executed is actually the code which is executed.

Padding
-------

With correct padding cache and TLB thrashing can be avoided. Therefore the
number of (fluid) nodes used in the data layout is artificially increased.

Currently automatic padding is active for kernels which support it. It can be
controlled via the kernel parameter (i.e. parameter after the ``--``)
``-pad``. Supported values are ``auto`` (default), ``no`` (to disable padding),
or a manual padding.

Automatic padding tries to avoid cache and TLB thrashing and pads for a 32
entry (huge pages) TLB with 8 sets and a 512 set (L2) cache. This reflects the
parameters of current Intel based processors.

Manual padding is done via a padding string and has the format
``mod_1+offset_1(,mod_n+offset_n)``, which specifies numbers of bytes.
SoA data layouts can exhibit TLB thrashing. Therefore we want to distribute the
19 pages with one lattice (36 with two lattices) we are concurrently accessing
over as much sets in the TLB as possible.
This is controlled by the distance between the accessed pages, which is the
number of (fluid) nodes in between them and can be adjusted by adding further
(fluid) nodes.
We want the distance d (in bytes) between two accessed pages to be e.g. 
**d % (PAGE_SIZE * TLB_SETS) = PAGE_SIZE**. 
This would distribute the pages evenly over the sets. Hereby **PAGE_SIZE * TLB_SETS**
would be our ``mod_1`` and **PAGE_SIZE** (after the =) our ``offset_1``.
Measurements show that with only a quarter of half of a page size as offset
higher performance is achieved, which is done by automatic padding.
On top of this padding more paddings can be added. They are just added to the
padding string and are separated by commas.

A zero modulus in the padding string has a special meaning. Here the
corresponding offset is just added to the number of nodes. A padding string
like ``-pad 0+16`` would at a static padding of two nodes (one node = 8 b).


Geometries
==========

TODO: supported geometries: channel, pipe, blocks


Results
=======

TODO


Licence
=======

The Lattice Boltzmann Benchmark Kernels are licensed under GPLv3.


Acknowledgements
================

This work was funded by BMBF, grant no. 01IH15003A (project SKAMPY).

This work was funded by KONWHIR project OMI4PAPS.



.. |datetime| date:: %Y-%m-%d %H:%M

Document was generated at |datetime|.

