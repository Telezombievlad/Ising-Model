#==================================================================================================
# COMPILER FLAGS
#==================================================================================================

CCFLAGS += -std=c++17 -Werror -Wall -O3 -mavx -fno-stack-protector -pthread

#==================================================================================================
# INSTALLATION
#==================================================================================================

install : 
	# cnpy
	git clone https://github.com/rogersce/cnpy.git model/vendor/cnpy
	mkdir -p model/vendor model/vendor/cnpy/build
	cd model/vendor/cnpy/build && cmake ..
	cd model/vendor/cnpy/build && sudo make && sudo make install

#==================================================================================================
# CNPY STUFF 
#==================================================================================================

MKFILE_PATH := $(abspath $(lastword $(MAKEFILE_LIST)))
CUR_DIR     := $(dir $(MKFILE_PATH))

SRC_ABS = ${CUR_DIR}model

LINK_TO_CNPY_FLAGS = -L/usr/local -lcnpy -lz

#==================================================================================================
# MODEL
#==================================================================================================

MODEL_SRC = model/model.cpp
MODEL_EXE = model/model

RESULT_FILE = res/THM.npy

model_compile : ${MODEL_SRC}
	g++ ${CCFLAGS} ${MODEL_SRC} -o ${MODEL_EXE} ${LINK_TO_CNPY_FLAGS}

model_compile_rendering : ${MODEL_SRC}
	g++ -DRENDERING ${CCFLAGS} ${MODEL_SRC} -o ${MODEL_EXE} ${LINK_TO_CNPY_FLAGS}

run :
	${MODEL_EXE} 0.0001 0.0 ${RESULT_FILE}

#==================================================================================================
# EXPERIMENTS
#==================================================================================================

MODEL_ASM = optimization/listing.asm

compile_debug : ${MODEL_SRC}
	g++ -S ${CCFLAGS} -g ${MODEL_SRC} -o ${MODEL_ASM}
	g++    ${CCFLAGS} -g ${MODEL_SRC} -o ${MODEL_EXE} ${LINK_TO_CNPY_FLAGS}

profile : compile_debug
	valgrind --tool=callgrind --dump-instr=yes --collect-jumps=yes ${MODEL_EXE}



