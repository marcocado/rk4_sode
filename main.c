#include<stdio.h>
#include<stdbool.h>
#include<math.h>

struct mValues{
    double m1;
    double m2;
    double m3;
    double m4;
};

struct kValues{
    double k1;
    double k2;
    double k3;
    double k4;
};

struct funcValues{
    double xn;
    double yn;
    double ydn;
};

double function(double x, double y, double yd){
    /*
    Here you can assert the equation for y'' = f(x, y, y'), you can also use math.h functions
    */
   double f = x + yd - y;
   return f;
}

double approximation(struct funcValues *currP, struct kValues *k, struct mValues *m, double h){
    /*
    First of all, the helper variables k1, k2, k3, k4 and m1, m2, m3, m4 are calculated. For better 
    readability and memory management the helper variables and function values are passed in with pointers
    to the structures.
    Afterwards y and y' of the current step are calculated.
    */
    k->k1 = h*currP->ydn;
    m->m1 = h*function(currP->xn, currP->yn, currP->ydn);
    k->k2 = h*(currP->ydn + m->m1/2);
    m->m2 = h*function(currP->xn + h/2, currP->yn + k->k1/2, currP->ydn + m->m1/2);
    k->k3 = h*(currP->ydn + m->m2/2);
    m->m3 = h*function(currP->xn + h/2, currP->yn + k->k2/2, currP->ydn + m->m2/2);
    k->k4 = h*(currP->ydn + m->m3);
    m->m4 = h*function(currP->xn + h, currP->yn + k->k3, currP->ydn + m->m3);

    currP->yn = currP->yn + 0.16666*(k->k1 + 2*k->k2 + 2*k->k3 + k->k4);
    currP->ydn = currP->yn + 0.16666*(m->m1 + 2*m->m2 + 2*m->m3 + m->m4);
}

void main(){
    // Declare the variables & structs
    int i;
    int n = 20;
    double h;
    struct funcValues startValues = {1, 2, 1};

    // Calculate the step increment
    h = (6 - startValues.xn)/n;
    printf("The step increment is: %.2f\n", h);

    // Initialize the struct and struct pointers
    struct kValues k;
    struct kValues *kPointer = &k;
    struct mValues m;
    struct mValues *mPointer = &m;
    struct funcValues f = startValues;
    struct funcValues *fPointer = &f;

    // Iterate through all step increments and calculate y and y'
    for(i=0; i<n; i++){
        approximation(fPointer, kPointer, mPointer, h);
        fPointer->xn = fPointer->xn + h;
        printf("\nFor the current step %.2f\n", fPointer->xn);
        printf("y:  %.3f\n", fPointer->yn);
        printf("y': %.3f\n", fPointer->ydn);
    }
}
