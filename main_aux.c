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
  double xint; //Valor muda de x donde queremos interpolar
  double fxint; //Es nuestra función
  double dfxint; //Primera derivada de nuestra función
  double d2fxint; //Segunda derivada de nuestra función

  int N = 1000; //interpolación con mil puntos
  //Es el espaciamiento entre dos puntos
  double delta_x = (x[n - 1] - x[0]) / N;

  //se abre el archivo para almacenar datos que vamos a interpolar
  output_interp_GSL = fopen("../interpol_GSL_cubic", "w");
  //checar si el archivo se abrio correctamente
  test_file_open(output_interp_GSL);
  //Se le asigna la memoria al constructor de la tabla de la memoria

  acc = gsl_interp_accel_alloc();
  //asignamos la memoria para hacer la interpolación cubica
  gsl_spline *cubic_spline;
  //Inicia la estructura de los datos
  cubic_spline = gsl_spline_alloc(gsl_interp_cspline, n);
  //Y se iniciaria la interpolación cubica
  gsl_spline_init(cubic_spline, x, y, n);
  //interpolación de los datos y almacenación de estos en un archivo
  for (int s = 0; s <= N; s++) {
      //interpolación
    xint = x[0] + s * delta_x;
    //interpolar en X_i
    fxint = gsl_spline_eval(cubic_spline, xint, acc);
    //Interpolación evaluada (como una función de X_i)
    dfxint = gsl_spline_eval_deriv(cubic_spline, xint, acc);
    //primera derivada de nuestra función
    d2fxint = gsl_spline_eval_deriv2(cubic_spline, xint, acc);
    //Segunda derivada
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint); //Almacenamiento de los datos recabados en un archivo señalado (output_interp_GSL, que es un puntero)
  }

  printf("La integral de la spline cubica es %.5f\n", gsl_spline_eval_integ(cubic_spline, x[0], x[n - 1], acc)); //Impresión de nuestra integral
  //donde se utiliza la función de la interpolación cubica, señalando los valores iniciales, el valor de la ultima interpolación y añadiendo el acc
  gsl_spline_free(cubic_spline); //Esta función libera el espacio de memoria que se reservó de nuestra interpolación cubica
  fclose(output_interp_GSL); //Cerramos el archivo que abrimos
  gsl_interp_accel_free(acc); //Liberar el espacio de memoria reservado de nuestro acc

  output_interp_GSL = fopen("../interpol_GSL_Akima", "w"); //Es para abrir el archivo sobre nuestro metodo akima
  test_file_open(output_interp_GSL); //Se le asigna la memoria al constructor de la tabla de memoria

  acc = gsl_interp_accel_alloc(); //Asignamos la memoria para la interpolación cubica
  gsl_spline *akima_spline; //Es para estructurar los datos
  akima_spline = gsl_spline_alloc(gsl_interp_akima, n); //Es para iniciar la interpolación
  gsl_spline_init(akima_spline, x, y, n); //Aqui se interpolan los datos y se almacenan en el archivo
  for (int s = 0; s <= N; s++) {
      //aqui se hace la interpolación
    xint = x[0] + s * delta_x;
    //función para evaluar la interpolación (con X_i)
    fxint = gsl_spline_eval(akima_spline, xint, acc);
    dfxint = gsl_spline_eval_deriv(akima_spline, xint, acc); //Primera derivada
    d2fxint = gsl_spline_eval_deriv2(akima_spline, xint, acc); //Segunda derivada
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint); //Aqui se imprime los datos de la función y sus derivadas
  }
  printf("La integral de la spline Akima es %.5f\n", gsl_spline_eval_integ(akima_spline, x[0], x[n - 1], acc));
  //Impresión de la integral donde se utiliza la función de la interpolación de akima, el valor inicial, el ultimo valor y el acc
  gsl_spline_free(akima_spline); //Libera el espacio de memoria de nuestra función akima
  fclose(output_interp_GSL); //Cierra el archivo
  gsl_interp_accel_free(acc); //Libera el espacio de memoria de nuestro acc

  output_interp_GSL = fopen("../interpol_GSL_linear", "w"); //Para el archivo de nuestro metodo de imterpolación lineal
  test_file_open(output_interp_GSL); //Aqui se asigna la memoria al constructor de los datos de memoria

  acc = gsl_interp_accel_alloc(); //Asignamos la memoria para la interpolación lineal
  gsl_spline *linear_spline; //Se estructuran los datos
  linear_spline = gsl_spline_alloc(gsl_interp_linear, n); //Se inicia la interpolación
  gsl_spline_init(linear_spline, x, y, n); //Se interpolan los datos y se guardan en el archivo
  for (int s = 0; s <= N; s++) {
      //se hace la interpolación
    xint = x[0] + s * delta_x;
    fxint = gsl_spline_eval(linear_spline, xint, acc); //función que evalua la interpolación
    dfxint = gsl_spline_eval_deriv(linear_spline, xint, acc); //Primera derivada
    d2fxint = gsl_spline_eval_deriv2(linear_spline, xint, acc); //Segunda derivada
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint); //Se imprime los datos de la función y su derivada
  }
  printf("La integral de la spline lineal es %.5f\n", gsl_spline_eval_integ(linear_spline, x[0], x[n - 1], acc));
    //Impresión de la integral donde se utiliza la función de la interpolación lineal, el valor inicial, el ultimo valor y el acc
  gsl_spline_free(linear_spline); //Libera el espacio de memoria de nuestra función lineal
  fclose(output_interp_GSL); //Cierra el archivo
  gsl_interp_accel_free(acc); //Libera el espacio de memoria de nuestro acc

  output_interp_GSL = fopen("../interpol_GSL_Steffen", "w"); //Es para abrir el archivo de la interpolación de Steffen
  test_file_open(output_interp_GSL); //Aqui se asigna la memoria al constructor de los datos de memoria

  acc = gsl_interp_accel_alloc(); //Asignamos la memoria para la interpolación de Steffen
  gsl_spline *steffen_spline; //Estructuración de los datos
  steffen_spline = gsl_spline_alloc(gsl_interp_steffen, n); //Se inicia la interpolación
  gsl_spline_init(steffen_spline, x, y, n); //Se interpolan y se guardan en el archivo
  for (int s = 0; s <= N; s++) {
      //Se hace la interpolación
    xint = x[0] + s * delta_x;
    fxint = gsl_spline_eval(steffen_spline, xint, acc); //Función de la interpolación
    dfxint = gsl_spline_eval_deriv(steffen_spline, xint, acc); //Primera derivada
    d2fxint = gsl_spline_eval_deriv2(steffen_spline, xint, acc); //Segunda derivada
    fprintf(output_interp_GSL, "%.5f %.5f %.5f %.5f\n", xint, fxint, dfxint, d2fxint); //Se imprimen los datos de la función y las derivadas
  }
  printf("La integral de la spline Steffen es %.5f\n", gsl_spline_eval_integ(steffen_spline, x[0], x[n - 1], acc));
    //Impresión de la integral donde se utiliza la función de la interpolación lineal, el valor inicial, el ultimo valor y el acc
  gsl_spline_free(steffen_spline); //Libera el espacio de memoria de nuestra función
  fclose(output_interp_GSL); //Cierra el archivo
  gsl_interp_accel_free(acc); //Libera el espacio de memoria de nuestro acc
  
}
