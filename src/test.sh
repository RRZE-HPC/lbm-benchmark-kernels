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

# How many parallel processes during make.
NProc="10"

Build=release

if [ "$#" -lt 1 ]; then
  echo "Compiles and runs several test cases."
  echo ""
  echo "Usage: test.sh <config>"
  echo ""
  echo "Select a configuration via config: linux-gcc or linux-intel."
  exit 1
fi

if [ "$1" == "-h" -o "$1" == "-help" -o "$1" == "--help" ]; then
  echo "Compiles and runs several test cases."
  echo ""
  echo "Usage: test.sh <config>"
  echo ""
  echo "Select a configuration via config: linux-gcc or linux-intel."
  exit 1
fi

Config="$1"

make clean-all

make -j $NProc PRECISION=dp CONFIG=$Config TAG=$XTag-debug
make -j $NProc PRECISION=dp CONFIG=$Config BUILD=$Build TAG=$XTag-v VERIFICATION=on
make -j $NProc PRECISION=dp CONFIG=$Config BUILD=$Build TAG=$XTag-b BENCHMARK=on

BinaryVDp="../bin/lbmbenchk-$Config-$Build-dp$XTag-v"
BinaryBDp="../bin/lbmbenchk-$Config-$Build-dp$XTag-b"


make -j $NProc PRECISION=sp CONFIG=$Config TAG=$XTag-debug
make -j $NProc PRECISION=sp CONFIG=$Config BUILD=$Build TAG=$XTag-v VERIFICATION=on
make -j $NProc PRECISION=sp CONFIG=$Config BUILD=$Build TAG=$XTag-b BENCHMARK=on

BinaryVSp="../bin/lbmbenchk-$Config-$Build-sp$XTag-v"
BinaryBSp="../bin/lbmbenchk-$Config-$Build-sp$XTag-b"


echo "#"
echo "# [test.sh] ./test-verification.sh \"$BinaryVDp\""
echo "#"

./test-verification.sh "$BinaryVDp"

ExitCodeDp="$?"

echo "#"
echo "# [test.sh] ./test-verification.sh \"$BinaryVSp\""
echo "#"

./test-verification.sh "$BinaryVSp"

ExitCodeSp="$?"

ResultDp="errors occurred"
ResultSp="errors occurred"

if [ "$ExitCodeDp" == "0" ]; then ResultDp="OK"; fi
if [ "$ExitCodeSp" == "0" ]; then ResultSp="OK"; fi

echo "#"
echo "# [test.sh] test   double precision: $ResultDp   single precision: $ResultSp"
echo "#"

ExitCode="0"

if [ "$ExitCodeDp" != 0 -o "$ExitCodeSp" != 0 ]; then
  ExitCode="1"
fi

exit "$ExitCode"

