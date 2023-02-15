#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>


double function(double x, double y, double yd){
    /*
    function:   y'' = f(x, y, y')
    */
    return -2*yd + 3*y + 2*x;
}

void stepCalculation(double (*xArray)[], double (*yArray)[], double (*ydArray)[], double (*yddArray)[], double h, int currStep){
    /*
    Calculation of the classical Runge-Kutta 4th order method for 2nd order ordinary differential equations for one step
    */
    double k1, k2, k3, k4;
    double m1, m2, m3, m4;
    int beforeStep = currStep - 1;

    // Calculation of the temporary helper variables
    k1 = h*(*ydArray)[beforeStep];
    m1 = h*function((*xArray)[beforeStep], (*yArray)[beforeStep], (*ydArray)[beforeStep]);
    k2 = h*((*ydArray)[beforeStep] + m1/2.0);
    m2 = h*function((*xArray)[beforeStep]+h/2.0, (*yArray)[beforeStep]+k1/2.0, (*ydArray)[beforeStep]+m1/2.0);
    k3 = h*((*ydArray)[beforeStep] + m2/2);
    m3 = h*function((*xArray)[beforeStep]+h/2.0, (*yArray)[beforeStep]+k2/2.0, (*ydArray)[beforeStep]+m2/2.0);
    k4 = h*((*ydArray)[beforeStep] + m3);
    m4 = h*function((*xArray)[beforeStep]+h, (*yArray)[beforeStep]+k3, (*ydArray)[beforeStep]+m3);

    // Calculation and assignment of the variables with pointers
    (*xArray)[currStep] = (*xArray)[beforeStep] + h;
    (*yArray)[currStep] = (*yArray)[beforeStep] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    (*ydArray)[currStep] = (*ydArray)[beforeStep] + (1.0/6.0)*(m1 + 2*m2 + 2*m3 + m4);
    (*yddArray)[currStep] = function((*xArray)[beforeStep], (*yArray)[beforeStep], (*ydArray)[beforeStep]);
}

void printCalculation(double (*xArray)[], double (*yArray)[], double (*ydArray)[], double (*yddArray)[], int currStep){
    /*
    Print the output for each calculation step
    */
    printf("----------------------------\n");
    printf("x:      %.3f\n", (*xArray)[currStep]);
    if((*yArray)[currStep]>10000){
        printf("y:      %.4e\n", (*yArray)[currStep]);
    }
    else{
        printf("y:      %.3f\n", (*yArray)[currStep]);
    }
    if((*ydArray)[currStep]>10000){
        printf("yd:     %.4e\n", (*ydArray)[currStep]);
    }
    else{
        printf("yd:     %.3f\n", (*ydArray)[currStep]);
    }
    if((*yddArray)[currStep]>10000){
        printf("ydd:    %.4e\n", (*yddArray)[currStep]);
    }
    else{
        printf("ydd:    %.3f\n", (*yddArray)[currStep]);
    }
    printf("----------------------------\n\n");
}

void RungeKutta4(int sizeArray, double (*xArray)[], double (*yArray)[], double (*ydArray)[], double (*yddArray)[], double h){
    /*
    Iterate through all points and calculate the x, y, y' and y'' values for the step
    */
    int i;
    for(i=1;i<sizeArray;i++){
        stepCalculation(*xArray, *yArray, *ydArray, *yddArray, h, i);
        printCalculation(*xArray, *yArray, *ydArray, *yddArray, i);
    }
}

int main(void){
    unsigned int i;
    double a = 0, b = 1;    // Start and end values of the range of the calculation
    unsigned int n = 1000000000;  // Number of calculation steps
    double h;
    unsigned int sizeArray = n+1;
    bool memoryAllocation = true;

    // Calculation of the step width
    h = (b - a)/ (double) n;

    // Boundaries conditions for the second order differential equation
    double y0 = 0;
    double yd0 = 4;

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
        RungeKutta4(sizeArray, *ptrxArray, *ptryArray, *ptrydArray, *ptryddArray, h);
    }
    else{
        printf("Due to lacking memory allocation, the numerical process can't be started\n");
    }

    return 0;
}
