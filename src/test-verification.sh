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

Tmp="delme.test.sh.$(hostname).$$.tmp"
Binary="../bin/lbmbenchk-linux-intel-release"
NThreads="5"

TestsTotal="0"
TestsFailed="0"
TestsSucceeded="0"

if [ "$#" -ge "1" ]; then
  Binary="$1"
fi


function on_exit
{
  if [ -e "$Tmp" ]; then
    rm -f "$Tmp" 2>&1 || true
  fi
}

trap "on_exit" EXIT

function run_kernel
{
  local Binary="$1"
  local K="$2"      # Kernel name
  local T="$3"      # Number of threads
  local KernelArgs="${4:-""}"
  local BinaryArgs="${5:-""}"

  ((++TestsTotal))

  echo -n "$Binary -verify -kernel $K -t $T -pin $(seq -s , 0 $((T-1))) $BinaryArgs ${KernelArgs:+"-- "}$KernelArgs "

  if [ "$KernelArgs" == "" ]; then
    $Binary -verify -kernel $K -t $T -pin $(seq -s , 0 $((T-1))) $BinaryArgs > "$Tmp" 2>&1
  else
    $Binary -verify -kernel $K -t $T -pin $(seq -s , 0 $((T-1))) $BinaryArgs -- $KernelArgs > "$Tmp" 2>&1
  fi

  local ExitCode="$?"

  if [ "$ExitCode" != "0" ]; then
    echo ""
    cat "$Tmp"
    echo "$Binary -verify  -kernel $K"
    echo "Verification failed. Exit code = $ExitCode."
    ((++TestsFailed))
  else
    echo "OK"
    ((++TestsSucceeded))
  fi

}

for K in $("$Binary" -list | tail -n +7); do

  for T in $(seq 1 $NThreads); do

    run_kernel "$Binary" "$K" "$T"
    # run_kernel "$Binary" "$K" "$T" "" "-dims 17x17x17"

    # Check in the usage string, if the kernel accepts parameters for blocking.

    $Binary -kernel $K -- -h > "$Tmp" 2>&1

    LineParameterStart="$(grep -n "^Kernel parameters:" "$Tmp" | sed -e 's/:.*//')"

    if [ "$LineParameterStart" == "" ]; then
      continue
    fi

    tail -n +$LineParameterStart "$Tmp" | grep -q -- "-blk"
    ExitCode="$?"

    if [ "$ExitCode" == "0" ]; then
      # Kernel supports blocking
      run_kernel "$Binary" "$K" "$T" "-blk 7"
    fi

  done

done


echo "# Tests toal: $TestsTotal  succeeded: $TestsSucceeded  failed: $TestsFailed"
