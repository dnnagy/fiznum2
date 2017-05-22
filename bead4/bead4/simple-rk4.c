//
//  main.c
//  bead4
//
//  Created by Nagy Daniel on 2017. 05. 05..
//  Copyright © 2017. Nagy Daniel. All rights reserved.
//
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#define T long double
#define EPSILON 1e-9 // Ennel kisebb elemeket nem engedunk a foatloba !!
#define nduint unsigned long int
#define SWAPD(U,X,Y) {U aux = (X); (X) = (Y); (Y) = aux;}
#define ABS(X) ((X) > 0 ? (X) : (-(X)))

#ifndef DEBUG_LN
#define DEBUG_LN(X, ...) printf((X), ##__VA_ARGS__);printf("\n");
#endif


typedef T (*diff_ptr)(T, T); // dy/dx = diff_ptr(x,y), y = f(x).

T rk4_step(T xn, T yn,T h, diff_ptr F){
    T k1 = h*F(xn, yn);
    T k2 = h*F(xn + h/2.0, yn + k1/2.0);
    T k3 = h*F(xn + h/2.0, yn + k2/2.0);
    T k4 = h*F(xn + h, yn + k3);
    
    T ynp1 = yn + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    return ynp1;
};

void solve_RK4(diff_ptr F, T x0, T xn, T y0, nduint steps, const char* outputFileName){
    
    T xk = x0;
    T yk = y0;
    T h = (xn - x0)/steps;
    
    FILE* optr = fopen(outputFileName, "w");
    
    nduint i=0;
    while (i<steps) {
        yk = rk4_step(xk, yk, h, F);
        xk = xk+h;
        DEBUG_LN("step %lu: %Lf %Lf", i, xk, yk);
        fprintf(optr, "%Lf %Lf\n", xk, yk);
        i++;
    }
    fclose(optr);
    return;
};


T func(T x, T y){
    return 1 - sin(x)*x*x;
}

int main(int argc, const char * argv[]) {
    diff_ptr eq = &func;
    
    solve_RK4(eq, 0, 4, 0, 5000, "output.dat");
    return 0;
}
