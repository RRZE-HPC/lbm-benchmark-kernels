#!/bin/bash -l
# --------------------------------------------------------------------------
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
set -u
set -e

XTag="-test"

Build=release



Config=linux-intel

module purge
module add intel64

make clean-all
make -j CONFIG=$Config TAG=$XTag-debug
make -j CONFIG=$Config BUILD=$Build TAG=$XTag-v
make -j CONFIG=$Config BUILD=$Build TAG=$XTag-b BENCHMARK=on

BinaryV="../bin/lbmbenchk-$Config-$Build$XTag-v"
BinaryB="../bin/lbmbenchk-$Config-$Build$XTag-b"

./test-verification.sh "$BinaryV"



Config=linux-gcc

module purge
module add gcc

make clean-all

make -j CONFIG=$Config TAG=$XTag-debug
make -j CONFIG=$Config BUILD=$Build TAG=$XTag-v
make -j CONFIG=$Config BUILD=$Build TAG=$XTag-b BENCHMARK=on

BinaryV="../bin/lbmbenchk-$Config-$Build$XTag-v"
BinaryB="../bin/lbmbenchk-$Config-$Build$XTag-b"

./test-verification.sh "$BinaryV"

