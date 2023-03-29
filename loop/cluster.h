#ifndef cluster_h
#define cluster_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NDEBUG
#include <assert.h>

#include "chain.h"
#include "table.h"


/**
 * @brief Updates the table of clusters.
 *
 * This function iterates through each item in the table and performs various
 * checks and updates to the item's state, link_key, and link_spin values.
 *
 * @param t Pointer to the table of clusters.
 */
void cluster_update_table(table* t);

/**
 * @brief Updates the chain of clusters.
 *
 * This function iterates through each node in the chain and updates its state
 * based on the state of the corresponding item in the table.
 *
 * @param c Pointer to the chain of clusters.
 * @param t Pointer to the table of clusters.
 */
void cluster_update_chain(chain* c, table* t);

/**
 * @brief Links vertices in a chain of clusters.
 *
 * This function iterates through each node in the chain and updates the link_key
 * and link_spin values of the corresponding items in the table.
 *
 * @param c Pointer to the chain of clusters.
 * @param t Pointer to the table of clusters.
 */
void cluster_link_vertex(chain* c, table* t);

/**
 * @brief Clusters the table based on a specific item and spin.
 *
 * This function modifies the table of clusters by linking spins based on their type.
 * It uses the provided item_id and spin_id as the starting point for clustering.
 *
 * @param t Pointer to the table of clusters.
 * @param item_id Index of the item in the table to start clustering from.
 * @param spin_id Index of the spin in the item to start clustering from.
 * @param rng Pointer to the GSL random number generator.
 */
static void cluster_clustering(table* t, int item_id, int spin_id, gsl_rng* rng);

/**
 * @brief Traverses the table and performs clustering on each item.
 *
 * This function iterates through the table and performs clustering on each item
 * by calling the cluster_clustering function.
 *
 * @param t Pointer to the table of clusters.
 * @param rng Pointer to the GSL random number generator.
 */
void cluster_traverse(table* t, gsl_rng* rng);

#endif
