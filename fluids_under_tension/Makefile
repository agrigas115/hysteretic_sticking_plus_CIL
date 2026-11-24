CFLAGS = -std=gnu99 -O3 -Wall -Wextra

all: cil_bd attractive_cil

# Hysteretic sticking model, both activities (--model C or --model P)
cil_bd: cil_bd.c multitau_isf.h multitau_scalar.h
	$(CC) $(CFLAGS) -o $@ $< -lm

# Attractive particle model with CIL-C (the absorbing state)
attractive_cil: attractive_cil.c
	$(CC) $(CFLAGS) -o $@ $< -lm

# Tuned to the build machine's CPU (faster; results can differ in the last bits from the generic build)
native: CFLAGS += -march=native
native: cil_bd attractive_cil

clean:
	rm -f cil_bd attractive_cil
