//
// Created by Claudio on 23/05/23.
//

#ifndef INTERPOLATIONFD_MAIN_AUX_H
#define INTERPOLATIONFD_MAIN_AUX_H

#include <stdio.h>
#include <stdlib.h>

unsigned int count_data_by_rows(int counter);

void read_two_col_data(double *x, double *y);

void double_ptr_alloc(double **ptr, unsigned int size);

void ptr_alloc_test(const double *pointer_to_allocate);

void test_file_open(const FILE *pointer_to_file);

void print_two_col_data(unsigned int n, const double *x, const double *y);

void interpolGSL(unsigned int n, const double *x, const double *y);

#endif //INTERPOLATIONFD_MAIN_AUX_H
