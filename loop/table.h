#ifndef table_h
#define table_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#define NDEBUG
#include <assert.h>

#define NSPIN_MAX 4

extern uint64_t table_statistic_index;

typedef struct item{
    uint64_t key;
    int type;
    int nspin;
    int state[2*NSPIN_MAX];
    uint64_t link_key[2*NSPIN_MAX];
    int link_spin[2*NSPIN_MAX];
} item;

typedef struct table{
    int size;
    int n;
    uint64_t key;
    item* list;
} table;

/**
 * @brief Allocate memory for a new table with a given scale.
 * 
 * This function allocates memory for a new table of size 2^scale and initializes its elements.
 * Each element in the table is an item structure, which contains a key, type, nspin, and arrays for state and links.
 *
 * @param scale The power of 2 for determining the size of the table (size = 2^scale).
 * @return Pointer to the allocated table structure.
 */
table* table_alloc(int scale);

/**
 * @brief Free the memory allocated for a table structure.
 *
 * This function frees the memory allocated for the table structure and its list of items.
 *
 * @param t Pointer to the table structure to be freed.
 */
void table_free(table* t);

/**
 * @brief Reallocate memory for a table with double the original size.
 *
 * This function reallocates memory for a table with double the size of the original table and copies the elements.
 * The original table's list of items is freed after copying.
 *
 * @param t Pointer to the table structure to be reallocated.
 */
void table_realloc(table* t);

/**
 * @brief Calculate the hash value for a given key.
 *
 * This function calculates the hash value for a given key using the size of the table.
 *
 * @param t Pointer to the table structure.
 * @param key The key for which the hash value needs to be calculated.
 * @return The calculated hash value.
 */
int table_hash(table* t, uint64_t key);

/**
 * @brief Generate a unique key for a table.
 *
 * This function generates a unique key for a table by incrementing the table's key value
 * until an unused key is found.
 *
 * @param t Pointer to the table structure.
 * @return The generated unique key.
 */
uint64_t table_generate_key(table* t);

/**
 * @brief Search for an item in the table using a given key.
 *
 * This function searches for an item in the table using a given key by calculating the hash value
 * and returning the corresponding item in the table.
 *
 * @param t Pointer to the table structure.
 * @param key The key to be searched for in the table.
 * @return Pointer to the found item in the table.
 */
item* table_search_from_key(table* t, uint64_t key);

/**
 * @brief Print the current state of a table.
 *
 * This function prints the current state of a table, including its size, number of elements,
 * keys, types, nspins, and links.
 *
 * @param t Pointer to the table structure.
 */
void table_print_state(table* t);

#endif
