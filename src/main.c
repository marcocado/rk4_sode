#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>

#include "test_function.h"
#include "runge_kutta_4.h"

int main(void){
    double start_value = 0.0, end_value = 1.0;    // Start and end values of the range of the calculation
    unsigned int n = 20;  // Number of calculation steps
    double step_width;
    unsigned int size_array = n+1;
    bool is_memory_allocated = false;

    // Calculation of the step width
    step_width = (end_value - start_value)/ (double) n;

    // Boundaries conditions for the second order differential equation
    double start_y0 = 0.0;
    double start_yd0 = 4.0;

    // Initialize arrays, pointers to arrays and assign pointers to the arrays
    double (*x_array)[size_array] = malloc(sizeof(double)*size_array);
    if(x_array != NULL){
        is_memory_allocated = true;
    }
    else{
        is_memory_allocated = false;
    }

    double (*y_array)[size_array] = malloc(sizeof(double)*size_array);
    if(y_array != NULL){
        is_memory_allocated = true;
    }
    else{
        is_memory_allocated = false;
    }

    double (*yd_array)[size_array] = malloc(sizeof(double)*size_array);
    if(yd_array != NULL){
        is_memory_allocated = true;
    }
    else{
        is_memory_allocated = false;
    }

    double (*ydd_array)[size_array] = malloc(sizeof(double)*size_array);
    if(ydd_array != NULL){
        is_memory_allocated = true;
    }
    else{
        is_memory_allocated = false;
    }

    if(is_memory_allocated == true){
        // Assign the start conditions
        (*x_array)[0] = start_value;
        (*y_array)[0] = start_y0;
        (*yd_array)[0] = start_yd0;
        (*ydd_array)[0] = function(start_value, start_y0, start_yd0);

        // Start the iterative calculation process
        RungeKutta4(size_array, *x_array, *y_array, *yd_array, *ydd_array, step_width);
    }
    else{
        printf("Due to lacking memory allocation, the numerical process can't be started\n");
    }

    return 0;
}
