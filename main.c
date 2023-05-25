#include "main_aux.h"

int main() {
  unsigned int n;
  double *x;
  double *y;
  n = count_data_by_rows(0);
  double_ptr_alloc(&x, n);
  double_ptr_alloc(&y, n);
  read_two_col_data(x, y);
  print_two_col_data(n, x, y);
  interpolGSL(n, x, y);
  return 0;
}
