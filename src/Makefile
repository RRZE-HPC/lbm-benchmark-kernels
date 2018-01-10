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

# Default configuration to chose.
# CONFIG ?= linux-gcc
CONFIG ?= linux-intel

BUILD ?= release

BENCHMARK ?= on

# If we compile for BENCHMARK all verifcation and statistics are disabled,
# if not specified otherwise.
ifeq (on,$(BENCHMARK))
  VERIFICATION ?= off
  VTK_OUTPUT   ?= off
  STATISTICS   ?= off
else
  VERIFICATION ?= on
  STATISTICS   ?= on
  VTK_OUTPUT   ?= on
endif


COLOR ?= on

OPENMP ?= on

VTK_OUTPUT_ASCII ?= off

LID_DRIVEN_CAVITY ?= off

ISA ?= avx

LIKWID ?= off

# Which floating point precision to use: dp (double precision) or sp (single preicision)
PRECISION ?= dp

# Global settings for the Makefile
SHELL      = sh

ifdef V
  Q        =
else
  Q        = @
endif

AR         = ar
CMP        = cmp
CHMOD      = chmod
CPP        = cpp
DIFF       = diff
ECHO       = echo
ECHO_E     = echo -e
ECHO_N     = echo -n
GIT        = git
GREP       = grep
MKDIR      = mkdir
MKDIR_P    = mkdir -p
NM         = nm
RM         = rm
RM_F       = rm -f
RM_RF      = rm -rf
TAR        = tar
TAR_CZF    = tar czf
TOUCH      = touch
SED        = sed



# Where to store objects and dependency files.
OBJECT_DIR      = obj/$(CONFIG)-$(BUILD)$(PREC)$(TAG)
DEP_DIR         = obj/$(CONFIG)-$(BUILD)$(PREC)$(TAG)-dep

# Sources to consider.
SOURCES_C = Main.c Memory.c Geometry.c Kernel.c \
            Vtk.c Pinning.c Padding.c Lattice.c

# ------------------------------------------------------------------------
# NO CHANGE BELOW SHOULD BE NEEDED
# ------------------------------------------------------------------------

OBJ_C = $(foreach SOURCE,$(SOURCES_C),$(OBJECT_DIR)/$(SOURCE:%.c=%.o)) \
		$(OBJECT_DIR)/BenchKernelD3Q19Common_PullSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19Common_PullAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19Common_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19Common_PushAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19_PullSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19_PullAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19_PushAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19List_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListCommon_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19List_PushAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListCommon_PushAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19List_PullSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListCommon_PullSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19List_PullAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListCommon_PullAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAa_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAaCommon_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAa_PushAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAaCommon_PushAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAaRia_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAaRiaCommon_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAaPv_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListAaPvCommon_PushSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListPullSplitNt_PullSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19ListPullSplitNtCommon_PullSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19Aa_AaAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19AaCommon_AaAoS.o \
		$(OBJECT_DIR)/BenchKernelD3Q19Aa_AaSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19AaCommon_AaSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19AaVec_AaSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19AaVecCommon_AaSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19AaVecSl_AaSoA.o \
		$(OBJECT_DIR)/BenchKernelD3Q19AaVecSlCommon_AaSoA.o

OBJ = $(OBJ_C)

include ./config.$(CONFIG).mk


# Colorize command lines.
ifeq (on,$(COLOR))
  COLOR_GRAY        =\033[1;30m
  COLOR_LIGHT_GRAY  =\033[0;37m
  COLOR_CYAN        =\033[0;36m
  COLOR_LIGHT_CYAN  =\033[1;36m
  COLOR_MAGENTA     =\033[0;35m
  COLOR_DARK_MAGENTA=\033[1;35m
  COLOR_NO          =\033[0m
  COLOR_RED         =\033[0;31m
  COLOR_DARK_RED    =\033[1;31m
  COLOR_GREEN       =\033[0;32m
  COLOR_DARK_GREEEN =\033[1;32m
else
  COLOR_GRAY        =
  COLOR_LIGHT_GRAY  =
  COLOR_CYAN        =
  COLOR_LIGHT_CYAN  =
  COLOR_MAGENTA     =
  COLOR_DARK_MAGENTA=
  COLOR_NO          =
  COLOR_RED         =
  COLOR_DARK_RED    =
  COLOR_GREEN       =
  COLOR_DARK_GREEEN =
endif

# Dependencies which are not source files, that require
# rebuilding the application.
REBUILD_DEPS = $(filter-out %.d,$(MAKEFILE_LIST)) $(OBJECT_DIR)/.target $(DEP_DIR)/.target

BUILD_CONFIG=

# Verify computed flow field
ifeq (on,$(VERIFICATION))
  PP_FLAGS	+= $(D)VERIFICATION
endif

ifeq (on,$(VTK_OUTPUT))
  PP_FLAGS += $(D)VTK_OUTPUT
endif

ifeq (on,$(STATISTICS))
  PP_FLAGS += $(D)STATISTICS
endif

ifeq (on,$(VTK_OUTPUT_ASCII))
  PP_FLAGS += $(D)VTK_OUTPUT_ASCII
endif

ifeq (on,$(LID_DRIVEN_CAVITY))
  PP_FLAGS += $(D)LID_DRIVEN_CAVITY
endif

ifeq (avx,$(ISA))
  PP_FLAGS += $(D)VECTOR_AVX
endif
ifeq (sse,$(ISA))
  PP_FLAGS += $(D)VECTOR_SSE
endif

# Required for posix_memalign
PP_FLAGS += $(D)_POSIX_C_SOURCE=200112L

ifeq (on,$(LIKWID))
  PP_FLAGS += $(LIKWID_INC) $(D)HAVE_LIKWID
  LD_LIBS  += $(LIKWID_LIB) -llikwid
endif


ifeq (dp,$(PRECISION))
  PP_FLAGS += $(D)PRECISION_DP
  PREC=-dp
else
ifeq (sp,$(PRECISION))
  PP_FLAGS += $(D)PRECISION_SP
  PREC=-sp
else
  $(error PRECISION is only be allowed to be sp (single precision) or dp (doble precision))
endif
endif

# ARCH can only be assigned a string without a space. The space is escaped as
# a comma which we have to replace here.

ifdef TARCH
  ARCH_W_COMMA := $(TARCH)

  COMMA_ := ,
  override TARCH := $(subst $(COMMA_), ,$(TARCH))

  $(info $(shell echo -e "\n$(COLOR_GREEN)INFO: Automatically expanding comma in TARCH variable:  TARCH=$(TARCH).$(COLOR_NO)\n"))

endif

.phony: all clean clean-all

#$(info $(shell $(ECHO_E) "# Configuration:   CONFIG=$(COLOR_CYAN)$(CONFIG)$(COLOR_NO)   BUILD=$(COLOR_CYAN)$(BUILD)$(COLOR_NO)   VERIFICATION=$(COLOR_CYAN)$(VERIFICATION)$(COLOR_NO)   STATISTICS=$(COLOR_CYAN)$(STATISTICS)$(COLOR_NO)   VTK_OUTPUT=$(COLOR_CYAN)$(VTK_OUTPUT)$(COLOR_NO)  OPENMP=$(COLOR_CYAN)$(OPENMP)$(COLOR_NO)   ISA=$(COLOR_CYAN)$(ISA)$(COLOR_NO)   LIKWID=$(COLOR_CYAN)$(LIKWID)$(COLOR_NO)   TARCH=$(COLOR_CYAN)$(TARCH)$(COLOR_NO)   building $(.DEFAULT_GOAL)..."))

$(info $(shell $(ECHO_E) "# Configuration: CONFIG=$(COLOR_CYAN)$(CONFIG)$(COLOR_NO)   BUILD=$(COLOR_CYAN)$(BUILD)$(COLOR_NO)  PRECISION=$(COLOR_CYAN)$(PRECISION)$(COLOR_NO)"))
$(info $(shell $(ECHO_E) "#                OPENMP=$(COLOR_CYAN)$(OPENMP)$(COLOR_NO)   ISA=$(COLOR_CYAN)$(ISA)$(COLOR_NO)   LIKWID=$(COLOR_CYAN)$(LIKWID)$(COLOR_NO)"))
$(info $(shell $(ECHO_E) "#                VERIFICATION=$(COLOR_CYAN)$(VERIFICATION)$(COLOR_NO)   STATISTICS=$(COLOR_CYAN)$(STATISTICS)$(COLOR_NO)   VTK_OUTPUT=$(COLOR_CYAN)$(VTK_OUTPUT)$(COLOR_NO)"))
$(info $(shell $(ECHO_E) "#                      target=$(.DEFAULT_GOAL)"))

$(info # Object dir:      $(OBJECT_DIR))
$(info # Dependency dir:  $(DEP_DIR))
# $(info # objects:        $(OBJ))

BIN_DIR=../bin

all: $(BIN_DIR)/lbmbenchk-$(CONFIG)-$(BUILD)$(BUILD_CONFIG)$(PREC)$(TAG)


# ------------------------------------------------------------------------
# TODO: Include dependencies, if they exist and no clean target is specified.

#ifeq (,$(SKIP_DEP_INCLUDE))
#  ifneq (,$(OBJ))
#    TMP_OBJ__ = $(patsubst $(OBJECT_DIR)/%,$(DEP_DIR)/%,$(OBJ))
#    $(warning TMP_OBJ__: $(TMP_OBJ__))
#    $(warning OBJ: $(OBJ))
#    $(warning DEPS:      $(patsubst %.o,%.d,$(TMP_OBJ__)))
#    $(warning $(patsubst $(OBJECT_DIR)/%.o,$(DEP_DIR)/%.d,$(OBJ)))
    -include $(patsubst $(OBJECT_DIR)/%.o,$(DEP_DIR)/%.d,$(OBJ))
#  endif
#endif


$(BIN_DIR):
	[ -d "$@" ] || mkdir -p "$@"

$(BIN_DIR)/lbmbenchk-$(CONFIG)-$(BUILD)$(BUILD_CONFIG)$(PREC)$(TAG): $(OBJ) $(REBUILD_DEPS) $(DEP_DIR)/.target | $(BIN_DIR)
	@$(ECHO_E) "linking: $(COLOR_CYAN)$@$(COLOR_NO)"
	$(LD) $(LD_FLAGS) -o $@ $(filter-out $(REBUILD_DEPS),$^) $(LD_LIBS)
	@$(ECHO_E) "# Builded binary:      $(COLOR_CYAN)$@$(COLOR_NO)"
	@$(ECHO_E) "# Configuration was:   CONFIG=$(COLOR_CYAN)$(CONFIG)$(COLOR_NO)   BUILD=$(COLOR_CYAN)$(BUILD)$(COLOR_NO)  PRECISION=$(COLOR_CYAN)$(PRECISION)$(COLOR_NO)"
	@$(ECHO_E) "#                      OPENMP=$(COLOR_CYAN)$(OPENMP)$(COLOR_NO)   ISA=$(COLOR_CYAN)$(ISA)$(COLOR_NO)   LIKWID=$(COLOR_CYAN)$(LIKWID)$(COLOR_NO)"
	@$(ECHO_E) "#                      VERIFICATION=$(COLOR_CYAN)$(VERIFICATION)$(COLOR_NO)   STATISTICS=$(COLOR_CYAN)$(STATISTICS)$(COLOR_NO)   VTK_OUTPUT=$(COLOR_CYAN)$(VTK_OUTPUT)$(COLOR_NO)"
	@$(ECHO_E) "#                      target=$(.DEFAULT_GOAL)"

$(OBJECT_DIR)/%_SoA.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_SOA$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_SOA -c $< -o $@

$(OBJECT_DIR)/%_PushSoA.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_SOA$(COLOR_NO) $(COLOR_MAGENTA)PROP_MODEL_PUSH$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_SOA $(D)PROP_MODEL_PUSH -c $< -o $@

$(OBJECT_DIR)/%_PushAoS.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_AOS$(COLOR_NO) $(COLOR_MAGENTA)PROP_MODEL_PUSH$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_AOS $(D)PROP_MODEL_PUSH -c $< -o $@

$(OBJECT_DIR)/%_PullSoA.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_SOA$(COLOR_NO) $(COLOR_MAGENTA)PROP_MODEL_PULL$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_SOA $(D)PROP_MODEL_PULL -c $< -o $@

$(OBJECT_DIR)/%_PullAoS.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_AOS$(COLOR_NO) $(COLOR_MAGENTA)PROP_MODEL_PULL$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_AOS $(D)PROP_MODEL_PULL -c $< -o $@

$(OBJECT_DIR)/%_AoS.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_AOS$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_AOS -c $< -o $@

$(OBJECT_DIR)/%_AaAoS.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_AOS$(COLOR_NO) $(COLOR_MAGENTA)PROP_MODEL_AA$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_AOS $(D)PROP_MODEL_AA -c $< -o $@

$(OBJECT_DIR)/%_AaSoA.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)    $(COLOR_MAGENTA)DATA_LAYOUT_SOA$(COLOR_NO) $(COLOR_MAGENTA)PROP_MODEL_AA$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) $(D)DATA_LAYOUT_SOA $(D)PROP_MODEL_AA -c $< -o $@

$(OBJECT_DIR)/%.o: %.c $(REBUILD_DEPS)
	@$(ECHO_E) "compiling: $(COLOR_CYAN)$@$(COLOR_NO)"
	$(CC) $(strip $(C_FLAGS)) $(strip $(PP_FLAGS)) -c $< -o $@


# ------------------------------------------------------------------------
# Following targets in princibal make only sure, that the directories
# are created if the do not exist already.

$(DEP_DIR)/.target:
	$(Q)$(MKDIR_P) $(DEP_DIR)
	$(Q)$(TOUCH) $@

# ------------------------------------------------------------------------
# Current configuration.
MAKE_CFG = SYSTEM=$(SYSTEM) // BUILD=$(BUILD) // MAKEOVERRIDES=\"$(strip $(MAKEOVERRIDES))\" // VERIFICATION=$(VERIFICATION) // STATISTICS=$(STATISTICS) // VTK_OUTPUT=$(VTK_OUTPUT) // VTK_OUTPUT_ASCII=$(VTK_OUTPUT_ASCII) // LID_DRIVEN_CAVITY=$(LID_DRIVEN_CAVITY) // ISA=$(ISA) // LIKWID=$(LIKWID) // PRECISION=$(PRECISION)

# Compare current configuration to the last one so we know when to
# rebuild this system/target despite when sources have not changed.
$(OBJECT_DIR)/.target::
	$(Q)$(MKDIR_P) $(OBJECT_DIR)
	@ ( [ -f $(OBJECT_DIR)/.target ] && $(ECHO) $(MAKE_CFG) | $(CMP) -s $(OBJECT_DIR)/.target - )           \
       || ( $(ECHO) $(MAKE_CFG) > $(OBJECT_DIR)/.target; $(ECHO) "# configuration variables changed, rebuilding ..." )


clean:
	-rm $(OBJECT_DIR)/*.o $(DEP_DIR)/*.d

clean-all:
	-[ -d obj ] && rm -r obj
