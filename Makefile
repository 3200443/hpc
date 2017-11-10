FILE = main.c $(normal) $(SSE2) $(ROC)

normal = morpho.c mouvement.c nrutil.c test_mouvement.c vnrutil.c 
SSE2 = test_mouvement_SSE2.c mouvement_SSE2.c morpho_simd.c
ROC = matric_roc.c

# -- Paths ----------
SRC_PATH = src
OBJ_PATH = obj

INC_PATH = include

# -- Flags ----------
C_DEBUG_FLAGS = -g -O3
C_ARCH_FLAGS = -msse4.2
C_INC_FLAGS = -I$(INC_PATH)


CFLAGS =  $(C_DEBUG_FLAGS) $(C_ARCH_FLAGS) $(C_INC_FLAGS)
LDFLAGS = $(C_ARCH_FLAGS) $(C_INC_FLAGS)



#SOURCES=$(addprefix $(SRC_PATH)/, $(wildcard *.c))
#OBJECTS=$(addprefix ${OBJ_PATH}/, $(addsuffix .o, $(basename $(wildcard *.c))))
SRC = $(addprefix ${SRC_PATH}/, $(FILE))
OBJ = $(addprefix ${OBJ_PATH}/, $(addsuffix .o, $(basename $(FILE))))

#Macro
CC=gcc

EXEC=hpc_NO.exe


# -- Base rules ----------
$(OBJ_PATH)/%.o : $(SRC_PATH)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

#-----Main rule ----------
$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) -lm


clean:
	rm -f $(OBJ) $(EXEC)

test:
	./hpc_NO.exe
