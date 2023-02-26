#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>


double function(double x, double y, double yd){
    /*
    function:   y'' = f(x, y, y')
    */
    return -2.0*yd + 3.0*y + 2.0*x;
}

void calculate_step(double (*xArray)[], double (*yArray)[], double (*ydArray)[], double (*yddArray)[], double h, int current_step){
    /*
    Calculation of the classical Runge-Kutta 4th order method for 2nd order ordinary differential equations for one step
    */
    double k1, k2, k3, k4;
    double m1, m2, m3, m4;
    int before_step = current_step - 1;

    // Calculation of the temporary helper variables
    k1 = h*(*ydArray)[before_step];
    m1 = h*function((*xArray)[before_step], (*yArray)[before_step], (*ydArray)[before_step]);
    k2 = h*((*ydArray)[before_step] + m1/2.0);
    m2 = h*function((*xArray)[before_step]+h/2.0, (*yArray)[before_step]+k1/2.0, (*ydArray)[before_step]+m1/2.0);
    k3 = h*((*ydArray)[before_step] + m2/2.0);
    m3 = h*function((*xArray)[before_step]+h/2.0, (*yArray)[before_step]+k2/2.0, (*ydArray)[before_step]+m2/2.0);
    k4 = h*((*ydArray)[before_step] + m3);
    m4 = h*function((*xArray)[before_step]+h, (*yArray)[before_step]+k3, (*ydArray)[before_step]+m3);

    // Calculation and assignment of the variables with pointers
    (*xArray)[current_step] = (*xArray)[before_step] + h;
    (*yArray)[current_step] = (*yArray)[before_step] + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
    (*ydArray)[current_step] = (*ydArray)[before_step] + (1.0/6.0)*(m1 + 2.0*m2 + 2.0*m3 + m4);
    (*yddArray)[current_step] = function((*xArray)[before_step], (*yArray)[before_step], (*ydArray)[before_step]);
}

void print_step(double (*xArray)[], double (*yArray)[], double (*ydArray)[], double (*yddArray)[], int current_step){
    /*
    Print the output for each calculation step
    */
    printf("----------------------------\n");
    printf("x:      %.3f\n", (*xArray)[current_step]);
    if((*yArray)[current_step]>10000){
        printf("y:      %.4e\n", (*yArray)[current_step]);
    }
    else{
        printf("y:      %.3f\n", (*yArray)[current_step]);
    }
    if((*ydArray)[current_step]>10000){
        printf("yd:     %.4e\n", (*ydArray)[current_step]);
    }
    else{
        printf("yd:     %.3f\n", (*ydArray)[current_step]);
    }
    if((*yddArray)[current_step]>10000){
        printf("ydd:    %.4e\n", (*yddArray)[current_step]);
    }
    else{
        printf("ydd:    %.3f\n", (*yddArray)[current_step]);
    }
    printf("----------------------------\n\n");
}

void runga_kutta_4(int sizeArray, double (*xArray)[], double (*yArray)[], double (*ydArray)[], double (*yddArray)[], double h){
    /*
    Iterate through all points and calculate the x, y, y' and y'' values for the step
    */
    int i;
    print_step(*xArray, *yArray, *ydArray, *yddArray, 0);
    for(i=1;i<sizeArray;i++){
        calculate_step(*xArray, *yArray, *ydArray, *yddArray, h, i);
        print_step(*xArray, *yArray, *ydArray, *yddArray, i);
    }
}

int main(void){
    unsigned int i;
    double a = 0.0, b = 1.0;    // Start and end values of the range of the calculation
    unsigned int n = 100;  // Number of calculation steps
    double h;
    unsigned int sizeArray = n+1;
    bool memoryAllocation = true;

    // Calculation of the step width
    h = (b - a)/ (double) n;

    // Boundaries conditions for the second order differential equation
    double y0 = 0.0;
    double yd0 = 4.0;

    // Initialize arrays, pointers to arrays and assign pointers to the arrays
    double (*ptrxArray)[sizeArray] = malloc(sizeof(double*)*sizeArray);
    double *xArray = malloc(sizeof(double)*sizeArray);
    if(xArray != NULL && ptrxArray != NULL){
        ptrxArray = xArray;
    }
    else{
        memoryAllocation = false;
    }

    double (*ptryArray)[sizeArray] = malloc(sizeof(double*)*sizeArray);
    double *yArray = malloc(sizeof(double)*sizeArray);
    if(yArray != NULL && ptryArray != NULL){
        ptryArray = yArray;
    }
    else{
        memoryAllocation = false;
    }

    double (*ptrydArray)[sizeArray] = malloc(sizeof(double*)*sizeArray);
    double *ydArray = malloc(sizeof(double)*sizeArray);
    if(ydArray != NULL && ptrydArray != NULL){
        ptrydArray = ydArray;
    }
    else{
        memoryAllocation = false;
    }

    double (*ptryddArray)[sizeArray] = malloc(sizeof(double*)*sizeArray);
    double *yddArray = malloc(sizeof(double)*sizeArray);
    if(yddArray != NULL && ptryddArray != NULL){
        ptryddArray = yddArray;
    }
    else{
        memoryAllocation = false;
    }

    if(memoryAllocation == true){
        // Assign the start conditions
        (*ptrxArray)[0] = a;
        (*ptryArray)[0] = y0;
        (*ptrydArray)[0] = yd0;
        (*ptryddArray)[0] = function(a, y0, yd0);

        // Start the iterative calculation process
        runga_kutta_4(sizeArray, *ptrxArray, *ptryArray, *ptrydArray, *ptryddArray, h);
    }
    else{
        printf("Due to lacking memory allocation, the numerical process can't be started\n");
    }

    return 0;
}
