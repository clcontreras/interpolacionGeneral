//
// Created by Claudio on 23/05/23.
//

#include "main_aux.h"
#include <gsl/gsl_spline.h>

unsigned int count_data_by_rows(int counter) {
  counter = 0;
  FILE *data_input_file;
  data_input_file = fopen("../input_data", "r");
  test_file_open(data_input_file);
  fscanf(data_input_file, "%*[^\n]\n");
  while (!feof(data_input_file)) {
    fscanf(data_input_file, "%*[^\n]\n");
    counter++;
  }
  fclose(data_input_file);
  return counter;
}

void read_two_col_data(double *x, double *y) {
  FILE *data_input_file;
  data_input_file = fopen("../input_data", "r");
  test_file_open(data_input_file);
  fscanf(data_input_file, "%*[^\n]\n");
  double aux1;
  double aux2;
  double aux3;
  int i = 0;
  while (!feof(data_input_file)) {
    fscanf(data_input_file, "%lf %lf %lf", &aux1, &aux2, &aux3);
    x[i] = aux1;
    y[i] = aux2;
    i++;
  }
  fclose(data_input_file);
}

void ptr_alloc_test(const double *pointer_to_allocate) {
  if (pointer_to_allocate == NULL) {
    printf("Error al reservar memoria\n");
    exit(1);
  }
}

void test_file_open(const FILE *pointer_to_file) {
  if (pointer_to_file == NULL) {
    printf("No se pudo abrir el archivo\n");
    exit(1);
  }
}

void print_two_col_data(unsigned int n, const double *x, const double *y) {
  printf("Estos son los %d datos leidos:\n", n);
  for (int i = 0; i < n; i++) printf("%6.3lf  %6.3lf\n", x[i], y[i]);
}

void double_ptr_alloc(double **ptr, unsigned int size) {
  *ptr = (double *) malloc(size * sizeof(double));
  ptr_alloc_test(*ptr);
}

void interpolGSL(unsigned int n, const double *x, const double *y) {
  FILE *output_interp_GSL;
  gsl_interp_accel *acc;
  double xint;
  double fxint;
  double dfxint;
  double d2fxint;

  int N = 1000;
  double delta_x = (x[n - 1] - x[0]) / N;

  output_interp_GSL = fopen("../interpol_GSL_cubic", "w");
  test_file_open(output_interp_GSL);

  acc = gsl_interp_accel_alloc();
  gsl_spline *cubic_spline;
  cubic_spline = gsl_spline_alloc(gsl_interp_cspline, n);
  gsl_spline_init(cubic_spline, x, y, n);
  for (int s = 0; s <= N; s++) {
    xint = x[0] + s * delta_x;
    fxint = gsl_spline_eval(cubic_spline, xint, acc);
    dfxint = gsl_spline_eval_deriv(cubic_spline, xint, acc);
    d2fxint = gsl_spline_eval_deriv2(cubic_spline, xint, acc);
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint);
  }
  printf("La integral de la spline cubica es %.5f\n", gsl_spline_eval_integ(cubic_spline, x[0], x[n - 1], acc));
  gsl_spline_free(cubic_spline);
  fclose(output_interp_GSL);
  gsl_interp_accel_free(acc);

  output_interp_GSL = fopen("../interpol_GSL_Akima", "w");
  test_file_open(output_interp_GSL);

  acc = gsl_interp_accel_alloc();
  gsl_spline *akima_spline;
  akima_spline = gsl_spline_alloc(gsl_interp_akima, n);
  gsl_spline_init(akima_spline, x, y, n);
  for (int s = 0; s <= N; s++) {
    xint = x[0] + s * delta_x;
    fxint = gsl_spline_eval(akima_spline, xint, acc);
    dfxint = gsl_spline_eval_deriv(akima_spline, xint, acc);
    d2fxint = gsl_spline_eval_deriv2(akima_spline, xint, acc);
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint);
  }
  printf("La integral de la spline Akima es %.5f\n", gsl_spline_eval_integ(akima_spline, x[0], x[n - 1], acc));
  gsl_spline_free(akima_spline);
  fclose(output_interp_GSL);
  gsl_interp_accel_free(acc);

  output_interp_GSL = fopen("../interpol_GSL_linear", "w");
  test_file_open(output_interp_GSL);

  acc = gsl_interp_accel_alloc();
  gsl_spline *linear_spline;
  linear_spline = gsl_spline_alloc(gsl_interp_linear, n);
  gsl_spline_init(linear_spline, x, y, n);
  for (int s = 0; s <= N; s++) {
    xint = x[0] + s * delta_x;
    fxint = gsl_spline_eval(linear_spline, xint, acc);
    dfxint = gsl_spline_eval_deriv(linear_spline, xint, acc);
    d2fxint = gsl_spline_eval_deriv2(linear_spline, xint, acc);
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint);
  }
  printf("La integral de la spline lineal es %.5f\n", gsl_spline_eval_integ(linear_spline, x[0], x[n - 1], acc));
  gsl_spline_free(linear_spline);
  fclose(output_interp_GSL);
  gsl_interp_accel_free(acc);

  output_interp_GSL = fopen("../interpol_GSL_Steffen", "w");
  test_file_open(output_interp_GSL);

  acc = gsl_interp_accel_alloc();
  gsl_spline *steffen_spline;
  steffen_spline = gsl_spline_alloc(gsl_interp_steffen, n);
  gsl_spline_init(steffen_spline, x, y, n);
  for (int s = 0; s <= N; s++) {
    xint = x[0] + s * delta_x;
    fxint = gsl_spline_eval(steffen_spline, xint, acc);
    dfxint = gsl_spline_eval_deriv(steffen_spline, xint, acc);
    d2fxint = gsl_spline_eval_deriv2(steffen_spline, xint, acc);
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint);
  }
  printf("La integral de la spline Steffen es %.5f\n", gsl_spline_eval_integ(steffen_spline, x[0], x[n - 1], acc));
  gsl_spline_free(steffen_spline);
  fclose(output_interp_GSL);
  gsl_interp_accel_free(acc);
  
}
