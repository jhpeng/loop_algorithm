#ifndef chain_h
#define chain_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#define NDEBUG
#include <assert.h>

typedef struct kink{
    double tau;
    int state[2];
    int spin_id;
    uint64_t key;
} kink;

typedef struct chain{
    int flag;
    int size;
    int n;
    kink* node;
    int state;
} chain;

/**
 * @brief Allocates memory for a new chain of a given size and initializes its properties.
 *
 * @param size The size of the chain to be allocated.
 * @return A pointer to the newly allocated chain.
 */
chain* chain_alloc(int size);

/**
 * @brief Frees the memory allocated for a chain.
 *
 * @param c Pointer to the chain to be freed.
 */
void chain_free(chain* c);

/**
 * @brief Reallocates the memory for a chain to a new size while preserving the chain's properties.
 *
 * @param c Pointer to the chain to be reallocated.
 * @param size The new size for the chain.
 */
void chain_realloc(chain* c, int size);

/**
 * @brief Inserts an array of kinks into the chain while maintaining the sorted order of kinks by tau values.
 *
 * @param c Pointer to the chain where the kinks will be inserted.
 * @param tau Pointer to the array of tau values for the kinks to be inserted.
 * @param key Pointer to the array of keys for the kinks to be inserted.
 * @param n Number of kinks to be inserted.
 * @param spin_id The spin ID for the kinks to be inserted.
 */
void chain_insert(
                chain* c, 
                double* tau, 
                uint64_t* key, 
                int n, 
                int spin_id);

/**
 * @brief Prints the current state of the chain, including the flag, size, number of kinks, and initial state.
 *
 * @param c Pointer to the chain to be printed.
 */
void chain_print_state(chain* c);
#endif
