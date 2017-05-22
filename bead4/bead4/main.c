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


//Differential equation: dy/dx = diff_ptr(x,y), y = f(x).
typedef T (*diff_ptr)(T, T);

//Differential equation system: dx1/dt = diff_sys_ptr[0](t, x1, x2, ...); dx2/dt = diff_sys_ptr[1](t, x1, x2, ...) ...
//If we have a system of N equations, each function should have N+2 params. The N+2nd is N.
typedef T (*diff_sys_eq)(T*, nduint);
typedef diff_sys_eq* diff_sys_ptr;

// Calculate RK4 step of equation
T rk4_step(T xn, T yn,T h, diff_ptr F){
    T k1 = h*F(xn, yn);
    T k2 = h*F(xn + h/2.0, yn + k1/2.0);
    T k3 = h*F(xn + h/2.0, yn + k2/2.0);
    T k4 = h*F(xn + h, yn + k3);
    
    T ynp1 = yn + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);
    return ynp1;
};

//Solve a simple equation
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

//Copy vector
void copy_v(T* dst, T* src, nduint N){
    for (nduint i = 0; i<N; i++) {
        dst[i] = src[i];
    }
}

void print_v_f(FILE* o, T* v, nduint N){
    for (nduint i = 0; i<N; i++) {
        fprintf(o, "%Lf ", v[i]);
    }
    fprintf(o, "\n");
}
// Calculate RK4 step of the system
void rk4_sys_step(T* xn,
                  T** xnp1,
                  T** x_buf,
                  T* k1,
                  T* k2,
                  T* k3,
                  T* k4,
                  T h,
                  diff_sys_ptr system,
                  nduint nEq,
                  nduint nParam){
    
    // Use x_buf to calculate k2,k3,k4 from xn and k1
    for (int i=0; i<nEq; i++) {
        k1[i] = h*system[i](xn, nParam);
    }
    
    //(*x_buf)[0] = xn[0] + h/2.0;
    for (int j=1; j<nParam; j++) {
        (*x_buf)[j] = xn[j] + k1[j-1]/2.0; // j-1 !!!
    }
    
    for (int i=0; i<nEq; i++) {
        k2[i] = h*system[i]((*x_buf), nParam);
    }
    
    //(*x_buf)[0] = xn[0] + h/2.0;
    for (int j=1; j<nParam; j++) {
        (*x_buf)[j] = xn[j] + k2[j-1]/2.0;
    }
    
    for (int i=0; i<nEq; i++) {
        k3[i] = h*system[i]((*x_buf), nParam);
    }
    
    //(*x_buf)[0] = xn[0] + h;
    for (int j=1; j<nParam; j++) {
        (*x_buf)[j] = xn[j] + k3[j-1];
    }
    
    for (int i=0; i<nEq; i++) {
        k4[i] = h*system[i]((*x_buf), nParam);
    }
    
    for (int i=1; i<nParam; i++) {
        (*xnp1)[i] = xn[i] + (1.0/6.0)*(k1[i-1] + 2*k2[i-1] + 2*k3[i-1] + k4[i-1]);
    }
    
    (*xnp1)[0] = xn[0] + h;
}

void solve_rk4_sys(diff_sys_ptr system,
                   nduint nEq,
                   T* x0,
                   T tn,
                   nduint steps){
    
    nduint nParam = nEq+1;
    T* xnp1 = (T*)malloc(nEq*sizeof(T));
    
    T* k1 = (T*)malloc(nEq*sizeof(T));
    T* k2 = (T*)malloc(nEq*sizeof(T));
    T* k3 = (T*)malloc(nEq*sizeof(T));
    T* k4 = (T*)malloc(nEq*sizeof(T));
    T* xn = (T*)malloc(nParam*sizeof(T));
    T* x_buf = (T*)malloc(nParam*sizeof(T));
    
    T h = (tn-x0[0])/steps;
    
    copy_v(xn, x0, nParam);
    
    FILE* out = fopen("results.dat","w");
    while (xn[0] <= tn) {
        print_v_f(out, xn, nParam);
        rk4_sys_step(xn, &xnp1, &x_buf, k1, k2, k3, k4, h, system, nEq, nParam);
        copy_v(xn, xnp1, nParam);
    }
    fclose(out);
    
}

T func(T x, T y){
    return 1 - sin(x)*x*x;
}


/// Csillapitott oszcill.
T dydx(T* x, nduint nP){
    return x[2];
}

T dzdx(T* x, nduint nP){
    return -3.0*x[1] - 0.05*x[2];
}
void testOscill(){
    diff_sys_eq sys[2] = {&dydx, &dzdx};
    T x0[3] = {0,3.0,0};
    
    solve_rk4_sys(sys, 2, x0, 90, 90000);
}

///Lorenz attraktor
T dxdt(T* x, nduint n){
    T sigma = 10.0;
    return sigma*(x[2]-x[1]);
    
}
T dydt(T* x, nduint n){
    T rho = 28.0;
    return x[1]*(rho - x[3]) - x[2];
    
}
T dzdt(T* x, nduint n){
    T beta = 8.0/3.0;
    return x[1]*x[2] - beta*x[3];
}

void testLorenz(){
    diff_sys_eq sys[3] = {&dxdt, &dydt, &dzdt};
    T x0[4] = {0.0, 1.0,42.0,1.0};
    
    solve_rk4_sys(sys, 3, x0, 10000, 100000);
}


/// Fold-hold rendszer
T MFold = 5.9736e+24; //kg
T MHold = 7.349e+22; //kg
T Rapog = 405500000; //m
T Rperig = 363300000; //m
T Vapog = 964; // m/s
T Vperig = 1076; // m/s
T G = 6.67384e-11;

// i  0  1   2   3    4
// X: t, xh, yh, vhx, vhy
T dxhdt(T* x, nduint np){
    return x[3];
};
T dyhdt(T* x, nduint np){
    return x[4];
};
T dvhxdt(T* x, nduint np){
    T r = sqrtl(x[1]*x[1] + x[2]*x[2]);
    return -G*MFold*x[1]/(r*r*r);
}
T dvhydt(T* x, nduint np){
    T r = sqrtl(x[1]*x[1] + x[2]*x[2]);
    return -G*MFold*x[2]/(r*r*r);
}

void testEarthMoon(){
    
    diff_sys_eq sys[4] = {&dxhdt, &dyhdt, &dvhxdt, &dvhydt};
    T x0[5] = {0.0, Rapog/2, 0.0, 0.0, Vapog};
    
    
    solve_rk4_sys(sys, 4, x0, 500000000, 1000000);
}

int main(int argc, const char * argv[]) {
    
    //Simple RK4
    //diff_ptr eq = &func;
    //solve_RK4(eq, 0, 4, 0, 5000, "output.dat");
    
    testEarthMoon();
    //testOscill();
    //testLorenz();
    return 0;
}
