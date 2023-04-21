void calculate_step(double (*x_array)[], double (*y_array)[], double (*yd_array)[], double (*ydd_array)[], double h, int current_step){
    /*
    Calculation of the classical Runge-Kutta 4th order method for 2nd order ordinary differential equations for one step
    */
    double k1, k2, k3, k4;
    double m1, m2, m3, m4;
    int before_step = current_step - 1;

    // Calculation of the temporary helper variables
    k1 = h*(*yd_array)[before_step];
    m1 = h*function((*x_array)[before_step], (*y_array)[before_step], (*yd_array)[before_step]);
    k2 = h*((*yd_array)[before_step] + m1/2.0);
    m2 = h*function((*x_array)[before_step]+h/2.0, (*y_array)[before_step]+k1/2.0, (*yd_array)[before_step]+m1/2.0);
    k3 = h*((*yd_array)[before_step] + m2/2.0);
    m3 = h*function((*x_array)[before_step]+h/2.0, (*y_array)[before_step]+k2/2.0, (*yd_array)[before_step]+m2/2.0);
    k4 = h*((*yd_array)[before_step] + m3);
    m4 = h*function((*x_array)[before_step]+h, (*y_array)[before_step]+k3, (*yd_array)[before_step]+m3);

    // Calculation and assignment of the variables with pointers
    (*x_array)[current_step] = (*x_array)[before_step] + h;
    (*y_array)[current_step] = (*y_array)[before_step] + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
    (*yd_array)[current_step] = (*yd_array)[before_step] + (1.0/6.0)*(m1 + 2.0*m2 + 2.0*m3 + m4);
    (*ydd_array)[current_step] = function((*x_array)[before_step], (*y_array)[before_step], (*yd_array)[before_step]);
}

void print_current_step(double (*x_array)[], double (*y_array)[], double (*yd_array)[], double (*ydd_array)[], int current_step){
    /*
    Print the output for each calculation step
    */
    printf("----------------------------\n");
    printf("x:      %.3f\n", (*x_array)[current_step]);
    if((*y_array)[current_step]>10000){
        printf("y:      %.4e\n", (*y_array)[current_step]);
    }
    else{
        printf("y:      %.3f\n", (*y_array)[current_step]);
    }
    if((*yd_array)[current_step]>10000){
        printf("yd:     %.4e\n", (*yd_array)[current_step]);
    }
    else{
        printf("yd:     %.3f\n", (*yd_array)[current_step]);
    }
    if((*ydd_array)[current_step]>10000){
        printf("ydd:    %.4e\n", (*ydd_array)[current_step]);
    }
    else{
        printf("ydd:    %.3f\n", (*ydd_array)[current_step]);
    }
    printf("----------------------------\n\n");
}

void RungeKutta4(int sizeArray, double (*x_array)[], double (*y_array)[], double (*yd_array)[], double (*ydd_array)[], double h){
    /*
    Iterate through all points and calculate the x, y, y' and y'' values for the step
    */
    int i;
    print_current_step(*x_array, *y_array, *yd_array, *ydd_array, 0);
    for(i=1;i<sizeArray;i++){
        calculate_step(*x_array, *y_array, *yd_array, *ydd_array, h, i);
        print_current_step(*x_array, *y_array, *yd_array, *ydd_array, i);
    }
}