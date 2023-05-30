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
  double xint; //valor de x donde queremos interpolar
  double fxint; //resultado de la interpolación
  double dfxint; //primera derivada
  double d2fxint; //valor de la segunda derivada

  int N = 1000; //interpolar a mil puntos
  //espaciamiento entre los puntos
  double delta_x = (x[n - 1] - x[0]) / N;

  //abrimos el archivo para almacenar datos a interpolar
  output_interp_GSL = fopen("../interpol_GSL_cubic", "w");
  //checamos si el archivo fué abierto correctamente
  test_file_open(output_interp_GSL);
//Asignamos la memoria al constructor de la tabla de búsqueda
  acc = gsl_interp_accel_alloc();
  //asignamos la memoria para interpolacion cubica
  gsl_spline *cubic_spline;
  //incializa y reserva memoria para la estructura de los datos
  cubic_spline = gsl_spline_alloc(gsl_interp_cspline, n);
  //inicializa la interpolación cúbica
  gsl_spline_init(cubic_spline, x, y, n);
  //interpolacion de datos y se almacena la interpolacion en un archivo
  for (int s = 0; s <= N; s++) {
      //punto donde vamos a interpolar (xint)
    xint = x[0] + s * delta_x;
    //valor de la interpolacion
    fxint = gsl_spline_eval(cubic_spline, xint, acc);
    //primera derivada de la interpolacion (valor interpolado)
    dfxint = gsl_spline_eval_deriv(cubic_spline, xint, acc);
    //segunda derivada de la interpolacion
    d2fxint = gsl_spline_eval_deriv2(cubic_spline, xint, acc);
    //almacenamos los datos en el archivo
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint);
  }
  //imprime el valor de la integral de la interpolacion, donde la funcion pide entrada la interpolacion cubica, el valor inicial, el valor final de la interpolacion y el operador
  printf("La integral de la spline cubica es %.5f\n", gsl_spline_eval_integ(cubic_spline, x[0], x[n - 1], acc));
  //liberamos la memoria que se reservó para la spline cubica
  gsl_spline_free(cubic_spline);
  //cerramos el archivo que abrimos
  fclose(output_interp_GSL);
  gsl_interp_accel_free(acc); //liberamos el espacio de memoria del operador acc tambien

  //abrimos el archivo para almacenar los datos a interpolar por el método de akima
  output_interp_GSL = fopen("../interpol_GSL_Akima", "w");
  //verificamos que el archivo fuera abierto correctamente
  test_file_open(output_interp_GSL);

  //Asignamos la memoria al constructor de la tabla de búsqueda
  acc = gsl_interp_accel_alloc();
  //asignamos la memoria a la interpolacion akima
  gsl_spline *akima_spline;
  //incializa y reserva memoria para la estructura de los datos
  akima_spline = gsl_spline_alloc(gsl_interp_akima, n);
  //inicializa la interpolacion akima
  gsl_spline_init(akima_spline, x, y, n);
  //interpolacion de datos, almacenamos la interpolacion en un archivo
  for (int s = 0; s <= N; s++) {
      //valor donde queremos interpolar para x
    xint = x[0] + s * delta_x; //valor de la interpolacion
    fxint = gsl_spline_eval(akima_spline, xint, acc); //valor de la primera derivada de la interpolacion
    dfxint = gsl_spline_eval_deriv(akima_spline, xint, acc); //valor de la segunda derivada de la interpolacion
    d2fxint = gsl_spline_eval_deriv2(akima_spline, xint, acc);
    //imprimimos los valores que acabamos de obtener, los almacenamos en un archivo
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint);
  }
  //imprime el valor de la integral de la interpolacion, donde la funcion pide de entrada la interp akima, el valor inicial, el valor final y el operador acc
  printf("La integral de la spline Akima es %.5f\n", gsl_spline_eval_integ(akima_spline, x[0], x[n - 1], acc));
  gsl_spline_free(akima_spline); //liberamos el espacio de memoria dedicado al archivo para la interpolacion akima
  fclose(output_interp_GSL); //cerramos el archivo
  gsl_interp_accel_free(acc); //libera el espacio de memoria para acc

  //abrimos el archivo para almacenar los datos de la interpolacion lineal
  output_interp_GSL = fopen("../interpol_GSL_linear", "w");
  test_file_open(output_interp_GSL); //checamos que el archivo fué abierto correctamente

  acc = gsl_interp_accel_alloc(); //asignamos la memoria al constructor de la tabla de busqueda
  gsl_spline *linear_spline; //asignamos el espacio de memoria a la interpolacion lineal
  linear_spline = gsl_spline_alloc(gsl_interp_linear, n); //incializa y reserva memoria para la estructura de los datos
  gsl_spline_init(linear_spline, x, y, n); //inicializa la interpolacion lineal
  for (int s = 0; s <= N; s++) {
    xint = x[0] + s * delta_x; //valor de x donde queremos interpolar
    fxint = gsl_spline_eval(linear_spline, xint, acc); //valor de la interpolacion
    dfxint = gsl_spline_eval_deriv(linear_spline, xint, acc); //valor de la primer derivada de la interpolacion
    d2fxint = gsl_spline_eval_deriv2(linear_spline, xint, acc); //valor de la segunda derivada de la interpolacion
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint); //imprimimos los valores obtenidos y los asignamos a un archivo
  }
  //imprimimos el valor de la integral de la interP lineal, donde nos pide el valor de la interpolacion, valor incial, valor final de la interpolacion y acc
  printf("La integral de la spline lineal es %.5f\n", gsl_spline_eval_integ(linear_spline, x[0], x[n - 1], acc));
  gsl_spline_free(linear_spline); //liberamos el espacio de memoria para el archivo de la interP lineal
  fclose(output_interp_GSL); //cerramos el archivo
  gsl_interp_accel_free(acc); //liberamos espacio de memoria

  //abrimos el archivo para almacenar los datos de la interP de Steffen
  output_interp_GSL = fopen("../interpol_GSL_Steffen", "w");
  test_file_open(output_interp_GSL); //checamos que el archivo fué abierto correctamente

  acc = gsl_interp_accel_alloc(); //asignamos la memoria al constructor de la tabla de busqueda
  gsl_spline *steffen_spline; //asignamos el espacio de memoria a la interpolacion de Steffen
  steffen_spline = gsl_spline_alloc(gsl_interp_steffen, n); //inicializa y reserva memoria para la estrcutura de los datos
  gsl_spline_init(steffen_spline, x, y, n); //inicializa la interpolacion steffen
  for (int s = 0; s <= N; s++) {
    xint = x[0] + s * delta_x; //valor de x donde queremos interpolar
    fxint = gsl_spline_eval(steffen_spline, xint, acc); //valor de la interpolacion
    dfxint = gsl_spline_eval_deriv(steffen_spline, xint, acc); //valor de la primer derivada de la interpolacion
    d2fxint = gsl_spline_eval_deriv2(steffen_spline, xint, acc); //valor de la segunda derivada de la interpolacion
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint); //imprime los valores obtenidos y los asigna a un archivo
  }
  //imprime el valor de la integral de la InterP steffen
  printf("La integral de la spline Steffen es %.5f\n", gsl_spline_eval_integ(steffen_spline, x[0], x[n - 1], acc));
  gsl_spline_free(steffen_spline); //libera el espacio de memoria para el archivo de la interP steffen
  fclose(output_interp_GSL); //cierra el archivo
  gsl_interp_accel_free(acc); //liberamos espacio de memoria de la tabla de busqueda
  
}
