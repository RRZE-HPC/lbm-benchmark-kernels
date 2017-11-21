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

Introduction
============

The lattice Boltzmann (LBM) benchmark kernels are a collection of LBM kernel
implementations.

**AS SUCH THE LBM BENCHMARK KERNELS ARE NO FULLY EQUIPPED CFD SOLVER AND SOLELY
SERVES THE PURPOSE OF STUDYING POSSIBLE PERFORMANCE OPTIMIZATIONS AND/OR
EXPERIMENTS.**

Currently all kernels utilize a D3Q19 discretization and the
two-relaxation-time (TRT) collision operator [ginzburg-2008]_.
All operations are carried out in double precision arithmetic.

Compilation
===========

The benchmark framework currently supports only Linux systems and the GCC and 
Intel compilers. Every other configuration probably requires adjustment inside
the code and the makefiles. Furthermore some code might be platform or at least
POSIX specific.

The benchmark can be build via ``make`` from the ``src`` subdirectory. This will
generate one binary which hosts all implemented benchmark kernels. 

Binaries are located under the ``bin`` subdirectory and will have different names
depending on compiler and build configuration.

Compilation can target debug or release builds. Combined with both build types
verification can be enabled, which increases the runtime and hence is not
suited for benchmarking.


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


Release and Verification
------------------------

Verification with the debug builds can be extremely slow. Hence verification
capabilities can be build with release builds: ::

  make BENCHMARK=off 


Benchmarking
------------

To generate a binary for benchmarking run make with ::

  make 

As default ``BENCHMARK=on`` and ``BUILD=release`` is set, where
``BUILD=release`` turns optimizations on and ``BENCHMARK=on`` disables
verfification, statistics, and VTK output.

See Options Summary below for further description of options which can be
applied, e.g. TARCH as well as the Benchmarking section.

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

Options that can be specified when building the suite with make:

============= ======================= ============ ==========================================================
name          values                  default      description
============= ======================= ============ ==========================================================
BENCHMARK     on, off                 on           If enabled, disables VERIFICATION, STATISTICS, VTK_OUTPUT. If disabled enables the three former options.
BUILD         debug, release          release      debug: no optimization, debug symbols, DEBUG defined. release: optimizations enabled.
CONFIG        linux-gcc, linux-intel  linux-intel  Select GCC or Intel compiler. 
ISA           avx, sse                avx          Determines which ISA extension is used for macro definitions of the intrinsics. This is *not* the architecture the compiler generates code for.
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

Kernel specific parameters can be obtained via selecting the specific kernel
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
  
Intel Compiler
--------------

For the Intel compiler one can specify depending on the target ISA extension:

- AVX:          ``TARCH=-xAVX``
- AVX2 and FMA: ``TARCH=-xCORE-AVX2,-fma``
- AVX512:       ``TARCH=-xCORE-AVX512``
- KNL:          ``TARCH=-xMIC-AVX512``

Compiling for an architecture supporting AVX (Sandy Bridge, Ivy Bridge): ::

  make ISA=avx TARCH=-xAVX


Compiling for an architecture supporting AVX2 (Haswell, Broadwell): ::

  make ISA=avx TARCH=-xCORE-AVX2,-fma

WARNING: ISA is here still set to ``avx`` as currently we have the FMA intrinsics not
implemented. This might change in the future.


Compiling for an architecture supporting AVX-512 (Skylake): ::

  make ISA=avx TARCH=-xCORE-AVX512

WARNING: ISA is here still set to ``avx`` as currently we have no implementation for the
AVX512 intrinsics. This might change in the future.


Pinning
-------

During benchmarking pinning should be used via the ``-pin`` parameter. Running
a benchmark with 10 threads and pin them to the first 10 cores works like ::

  $ bin/lbmbenchk-linux-intel-release ... -t 10 -pin $(seq -s , 0 9)


General Remarks
---------------

Things the binary does nor check or control:

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

- System load: interference with other application, especially on desktop
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

TODO: supported geometries: channel, pipe, blocks, fluid


Performance Results
===================

The sections lists performance values measured on several machines for
different kernels and geometries.
The **RFM** column denotes the expected performance as predicted by the
Roofline performance model [williams-2008]_. 
For performance prediction of each kernel a memory bandwidth benchmark is used
which mimics the kernels memory access pattern and the kernel's loop balance
(see [kernels]_ for details).

Haswell, Intel Xeon E5-2695 v3
------------------------------

- Haswell architecture, AVX2, FMA
- 14 cores, 2,3 GHz
- 2 x 7 cores in cluster-on-die (CoD) mode enabled
- SMT enabled

memory bandwidth:

- copy-19              47.3 GB/s
- copy-19-nt-sl        47.1 GB/s
- update-19            44.0 GB/s

geometry dimensions:  500x100x100

=========================    =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =====
kernel                            pipe   blocks-2   blocks-4   blocks-6   blocks-8  blocks-10  blocks-15  blocks-16  blocks-20  blocks-25  blocks-32  RFM
=========================    =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =====
blk-push-aos                     58.82      49.85      57.34      59.90      61.37      62.17      65.30      64.00      67.54      64.46      69.69   104 
blk-push-soa                     32.32      33.46      34.02      34.64      35.06      35.04      36.31      35.44      37.20      35.14      37.95   104
blk-pull-aos                     56.97      51.41      56.09      57.92      59.98      59.83      63.37      61.55      65.50      63.11      67.02   104
blk-pull-soa                     49.29      46.23      47.50      51.97      51.27      49.52      55.23      53.13      54.50      49.79      57.90   104
aa-aos                           91.35      66.14      76.80      84.76      83.63      91.36      93.46      92.62      93.91      92.25      92.93   145
aa-soa                           75.51      65.68      70.94      71.36      73.83      75.46      74.84      79.48      83.28      77.70      82.72   145
aa-vec-soa                       93.85      83.44      91.58      93.96      94.35      96.62     101.76      96.72     106.37     102.60     110.28   145
list-push-aos                    80.29      80.97      80.95      81.10      81.37      82.44      81.77      81.49      80.72      81.93      80.93   83
list-push-soa                    47.52      42.65      45.28      46.64      43.46      40.59      44.94      46.55      41.53      45.98      44.86   83
list-pull-aos                    85.30      82.97      86.43      83.42      86.33      83.70      86.43      83.77      83.10      85.89      84.44   83
list-pull-soa                    62.12      63.61      63.28      61.32      66.72      62.65      64.82      60.49      58.01      64.46      62.52   83
list-pull-split-nt-1s-soa       121.35     113.77     115.29     113.54     117.00     116.46     114.78     114.54     110.83     112.67     117.85   125
list-pull-split-nt-2s-soa       118.09     110.48     112.55     113.18     113.44     111.85     109.27     114.41     110.28     111.78     113.74   125
list-aa-aos                     121.28     118.63     119.00     118.50     121.99     119.11     118.83     121.47     121.62     126.18     120.12   129
list-aa-soa                     126.34     116.90     129.45     127.12     129.41     121.42     126.19     126.76     126.70     124.40     125.22   129
list-aa-ria-soa                 133.68     121.82     126.04     128.46     131.15     132.25     128.78     133.50     126.69     124.40     130.37   145
list-aa-pv-soa                  146.22     124.39     130.73     136.29     137.61     131.21     138.65     138.78     127.02     132.40     138.37   145
=========================    =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =====


Broadwell, Intel Xeon E5-2630 v4
--------------------------------

- Broadwell architecture, AVX2, FMA
- 10 cores, 2.2 GHz
- SMT disabled

memory bandwidth:

- copy-19              48.0 GB/s
- copy-nt-sl-19        48.2 GB/s
- update-19            51.1 GB/s

geometry dimensions:  500x100x100

=========================   =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =======
kernel                           pipe   blocks-2   blocks-4   blocks-6   blocks-8  blocks-10  blocks-15  blocks-16  blocks-20  blocks-25  blocks-32  RFM
=========================   =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =======
blk-push-aos                    55.75      47.62      54.57      57.10      58.49      59.00      61.72      60.56      64.05      61.10      66.03  105
blk-push-soa                    30.06      31.09      32.13      32.54      32.74      32.72      33.81      33.19      34.90      33.21      35.75  105
blk-pull-aos                    53.80      48.61      53.08      54.99      56.08      56.68      59.20      58.12      61.49      58.71      63.45  105
blk-pull-soa                    46.96      46.61      48.84      49.70      50.33      50.46      52.36      51.39      54.20      51.61      55.71  105
aa-aos                          91.40      66.99      78.47      83.38      86.62      88.62      92.98      91.54      97.08      94.93      98.90  168
aa-soa                          83.01      69.96      75.85      77.72      79.01      79.29      82.38      80.11      85.70      83.91      87.69  168
aa-vec-soa                     112.03      96.52     105.32     109.76     112.55     113.82     120.55     118.37     126.30     121.37     131.94  168
list-push-aos                   75.13      74.18      75.20      75.42      75.24      75.99      75.80      75.80      75.54      76.22      76.21   97
list-push-soa                   40.99      38.14      39.00      38.89      38.89      39.67      39.87      39.28      39.35      40.08      40.13   97
list-pull-aos                   82.07      82.88      83.29      83.09      83.32      83.49      82.82      82.88      83.32      82.60      82.93   97
list-pull-soa                   62.07      60.40      61.89      61.39      62.43      60.90      60.48      62.80      62.50      61.10      60.38   97
list-pull-split-nt-1s-soa      125.81     120.60     121.96     122.34     122.86     123.53     123.64     123.67     125.94     124.09     123.69  128
list-pull-split-nt-2s-soa      122.79     117.16     118.86     119.16     119.56     119.99     120.01     120.03     122.64     120.57     120.39  128
list-aa-aos                    128.13     127.41     129.31     129.07     129.79     129.63     129.67     129.94     129.12     128.41     129.72  150
list-aa-soa                    141.60     139.78     141.58     142.16     141.94     141.31     142.37     142.25     142.43     141.40     142.26  150
list-aa-ria-soa                141.82     134.88     140.15     140.72     141.67     140.51     141.18     141.29     142.97     141.94     143.25  168
list-aa-pv-soa                 164.79     140.95     159.24     161.78     162.40     163.04     164.69     164.38     165.11     165.75     166.09  168
=========================   =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =======


Skylake, Intel Xeon Gold 6148
-----------------------------

- Skylake architecture, AVX2, FMA, AVX512
- 20 cores, 2.4 GHz
- SMT enabled

memory bandwidth:

- copy-19                  89.7 GB/s
- copy-19-nt-sl            92.4 GB/s
- update-19                93.6 GB/s

geometry dimensions:  500x100x100


=========================    =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  ===
kernel                            pipe   blocks-2   blocks-4   blocks-6   blocks-8  blocks-10  blocks-15  blocks-16  blocks-20  blocks-25  blocks-32  RFM
=========================    =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  ===
blk-push-aos                    113.01      93.99     108.98     114.65     117.87     119.47     124.95     122.46     129.29     123.87     133.01  197
blk-push-soa                    100.21      98.87     103.63     105.56     107.02     107.27     111.61     109.83     116.16     110.51     110.29  197
blk-pull-aos                    118.45     102.54     114.12     117.82     122.69     124.31     130.58     127.85     135.72     129.65     139.94  197
blk-pull-soa                     82.60      83.36      87.13      88.39      88.84      88.96      92.48      90.93      95.79      91.92      98.64  197
aa-aos                          171.32     125.43     147.73     157.70     163.35     167.25     175.39     174.20     182.54     173.67     187.76  308
aa-soa                          180.85     152.39     165.84     152.59     171.90     175.76     184.94     182.34     189.43     180.30     193.54  308
aa-vec-soa                      208.03     181.51     195.86     203.41     209.08     212.34     224.05     219.49     234.31     225.92     245.22  308
list-push-aos                   158.81     164.67     162.93     163.05     165.22     164.31     164.66     160.78     164.07     165.19     164.06  177
list-push-soa                   134.60     110.44     110.17     132.01     132.95     133.46     134.37     134.33     135.12     134.91     137.87  177
list-pull-aos                   169.61     170.03     170.89     170.90     171.20     171.60     172.09     171.95     169.48     172.08     171.02  177
list-pull-soa                   120.50     116.73     118.62     118.00     120.99     118.15     117.17     121.41     120.83     120.00     118.74  177
list-pull-split-nt-1s-soa       225.59     224.18     225.10     226.34     226.01     230.37     227.50     228.42     227.39     231.65     227.35  246
list-pull-split-nt-2s-soa       219.20     214.63     217.61     218.13     219.07     221.01     219.88     220.09     220.62     221.68     220.58  246
list-aa-aos                     241.39     239.27     239.53     242.56     242.46     243.00     242.91     242.46     241.24     242.96     241.52  275
list-aa-soa                     273.73     268.49     268.48     271.79     275.29     274.56     277.18     272.67     274.21     275.24     278.21  275
list-aa-ria-soa                 288.42     261.89     273.26     284.84     283.88     288.29     290.72     289.81     293.36     290.75     292.93  308
list-aa-pv-soa                  303.35     267.21     289.18     294.96     294.36     298.16     300.45     301.71     302.37     302.88     304.46  308
=========================    =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  =========  ===

Licence
=======

The Lattice Boltzmann Benchmark Kernels are licensed under GPLv3.


Acknowledgements
================

This work was funded by BMBF, grant no. 01IH15003A (project SKAMPY).

This work was funded by KONWHIR project OMI4PAPS.


Bibliography
============

.. [ginzburg-2008]
 I. Ginzburg, F. Verhaeghe, and D. d'Humières. 
 Two-relaxation-time lattice Boltzmann scheme: About parametrization, velocity, pressure and mixed boundary conditions. 
 Commun. Comput. Phys., 3(2):427-478, 2008.

.. [williams-2008]
 S. Williams, A. Waterman, and D. Patterson. 
 Roofline: an insightful visual performance model for multicore architectures. 
 Commun. ACM, 52(4):65-76, Apr 2009. doi:10.1145/1498765.1498785


.. |datetime| date:: %Y-%m-%d %H:%M

Document was generated at |datetime|.

