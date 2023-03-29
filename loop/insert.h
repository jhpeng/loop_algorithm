#ifndef insert_h
#define insert_h

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NDEBUG
#include <assert.h>

#include "chain.h"
#include "table.h"

#define INSERT_MAX 100000

/**
 * @brief Inserts a horizontal graph between two chains.
 *
 * This function creates a horizontal graph between two chains (c1 and c2) and
 * updates the table of clusters accordingly. If the insert_tau array is NULL,
 * it allocates memory for the insert_tau, insert_size, and other related arrays.
 *
 * @param c1 Pointer to the first chain.
 * @param c2 Pointer to the second chain.
 * @param t Pointer to the table of clusters.
 * @param w The weight of the distribution.
 * @param beta The beta value used in the distribution.
 * @param rng Pointer to the GSL random number generator.
 */
void insert_horizontal_graph(
                chain* c1, 
                chain* c2, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

/**
 * @brief Inserts a triangular cut graph between three chains.
 *
 * This function inserts a triangular cut graph between the provided chains and
 * updates the table of clusters accordingly. If the insert_tau array is NULL,
 * it allocates memory for the insert_tau, chain_tau, and other related arrays.
 *
 * @param c Array of pointers to the three chains.
 * @param t Pointer to the table of clusters.
 * @param w The weight of the distribution.
 * @param beta The beta value used in the distribution.
 * @param rng Pointer to the GSL random number generator.
 */
void insert_triangular_cut_graph(
                chain** c, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

/**
 * @brief Inserts a triangular graph between three chains.
 *
 * This function inserts a triangular graph between the provided chains and
 * updates the table of clusters accordingly. If the insert_tau array is NULL,
 * it allocates memory for the insert_tau, chain_tau, and other related arrays.
 *
 * @param c Array of pointers to the three chains.
 * @param t Pointer to the table of clusters.
 * @param w The weight of the distribution.
 * @param beta The beta value used in the distribution.
 * @param rng Pointer to the GSL random number generator.
 */
void insert_triangular_graph(
                chain** c, 
                table* t, 
                double w, 
                double beta, 
                gsl_rng* rng);

#endif
