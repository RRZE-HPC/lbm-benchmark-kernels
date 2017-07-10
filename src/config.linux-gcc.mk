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


# ------------------------------------------------------------------------
# C ompiler/linker to use.
# Flags are specified at the end of the file.
# ------------------------------------------------------------------------
CC              = gcc
LD              = gcc

# ------------------------------------------------------------------------

# Preprocessing variables.
D          = -D
I          = -I
PP_FLAGS  +=

# Architecture to optimize for.
TARCH	   ?= -mavx

# Generated dependencies, can be left empty.
MAKE_DEPEND = $(CC) -MM -MQ'$(OBJECT_DIR)/$(<:%.c=%.o)' -MF'$(DEP_DIR)/$(<:%.c=%.d)' $(PP_FLAGS) $< > /dev/null

# Generates dependencies, can be left empty.
# $(call make_depend,<source-file>,<source-file-to-generate-dependency-for)
define make_depend
	$(CC) -MM -MQ'$(OBJECT_DIR)/$(2:%.c=%.o)' -MF'$(DEP_DIR)/$(2:%.c=%.d)' $(PP_FLAGS) $1 > /dev/null
endef


ifeq (on,$(OPENMP))
  OPENMP_C_FLAGS  += -fopenmp
  OPENMP_LD_FLAGS += -fopenmp
endif


ifeq (release,$(BUILD))

  C_FLAGS      += -O3 $(TARCH)
  LD_FLAGS     += -O3 $(TARCH)

  PP_FLAGS     +=

else
ifeq (debug,$(BUILD))

  C_FLAGS      += -O0 $(TARCH) -g -ggdb
  LD_FLAGS     += -O0 $(TARCH) -g -ggdb

  PP_FLAGS     += $(D)DEBUG

else
  $(error unknown BUILD=$(BUILD), specify release or debug)
endif
endif

ifeq (on,$(ADDRESS_SANITIZER))
  # see https://github.com/google/sanitizers/wiki/AddressSanitizerFlags for details
  C_FLAGS      += -fsanitize=address
  LD_FLAGS     += -fsanitize=address
endif

ifeq (on,$(DEBUG_SYMBOLS))
  C_FLAGS      += -g -ggdb
  LD_FLAGS     += -g -ggdb
endif

C_FLAGS        += -Wall -Wuninitialized -Wunused-variable -Wshadow -fargument-noalias -fargument-noalias-anything -std=c99 \
 -MT $@ -MF $(patsubst $(OBJECT_DIR)/%.o,$(DEP_DIR)/%.d,$@) -MMD $(OPENMP_C_FLAGS)
LD_FLAGS       += -Wall $(OPENMP_LD_FLAGS)
LD_LIBS        += -lm

