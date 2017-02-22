CC        = mpiicc
GCFLAGS    = ./COMBINATORIAL_STRATEGIES/hsl_mc73/hsl_mc73d.o ./COMBINATORIAL_STRATEGIES/hsl_mc73/libhsl_mc73.a \
	     ./COMBINATORIAL_STRATEGIES/hsl_mc64/mc64d.o ./COMBINATORIAL_STRATEGIES/hsl_mc64/libmc64.a \
	     -L./PARALLEL/pardiso -lpardiso500-INTEL1301-X86-64 \
	     -limf -lblas -llapack -lm -fopenmp -lifcore -Wall -Ofast -march=native
SOURCES    = ./COMBINATORIAL_STRATEGIES/min_max.c \
	     ./COMBINATORIAL_STRATEGIES/matching.c\
             ./COMBINATORIAL_STRATEGIES/spectral.c\
             ./COMBINATORIAL_STRATEGIES/de_min.c  \
             ./COMBINATORIAL_STRATEGIES/fib_heap.c\
             ./COMBINATORIAL_STRATEGIES/list.c    \
	     ./COMMON_FILES/common.c              \
	     ./COMMON_FILES/matrix.c              \
	     ./PARALLEL/parallel.c                \
	     ./PARALLEL/pardiso.c                 \
	     program.c
OBJECTS    = $(SOURCES:.c=.o)
EXECUTABLE = program 

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(GCFLAGS) -o $@

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)
	rm -f data*

