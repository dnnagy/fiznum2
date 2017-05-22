//
//  main.c
//  bead3
//
//  Created by Nagy Daniel on 2017. 04. 27..
//  Copyright © 2017. Nagy Daniel. All rights reserved.
//

#include <stdio.h>
#include "gauss.h"

nduint countLines(const char* file){
    FILE* fp = fopen(file, "r");
    char ch;
    nduint lines = 0;
    while(!feof(fp))
    {
        ch = fgetc(fp);
        if(ch == '\n')
        {
            lines++;
        }
    }
    return lines;
}

void read_input(struct ndMatrix* const mt, nduint cols,const char* filename){
    mt->R = countLines(filename);
    mt->data = (T*)malloc(mt->C*mt->R*sizeof(T));
    
    FILE* in = fopen(filename, "r");
    nduint i = 0;
    nduint ok = 1;
    
    while (ok) {
        if ( fscanf(in, "%lf", (mt->data) + i*mt->C) == 1 ) {
            for (nduint j=1; j<mt->C; j++) {
                if( fscanf(in, "%lf", (mt->data) + i*mt->C + j ) != 1){
                    ok = 0;
                }
            }
            i++;
        }else{
            ok = 0;
        }
    }
}


void print_results(const char* yfilename,
                   const char* paramsFilename,
                   const struct ndMatrix* const input,
                   const struct ndMatrix* const plan,
                   nduint N,
                   nduint rend,
                   T* const a){
    /// Kiiratni az eredmenyeket
    FILE* ki = fopen(yfilename, "w");
    
    
    for (nduint i = 0; i<plan->R; i++) {
        T y=0;
        T sigma = E(input, i, N+1);
        for (nduint k=0; k<plan->C; k++) {
            y+= *(a+k)*E(plan, i, k);
        }
        fprintf(ki, "%lf %lf\n", E(input, i, N), y*sigma);
    }
    fclose(ki);
    
    FILE* ki2 = fopen(paramsFilename, "w");
    for(int k=0; k<N*rend+1; k++) {
        fprintf(ki2, "%lf\n", *(a+k));
    }
    fclose(ki2);

};

//Megcsinalni a terv matrixot az inputbol
struct ndMatrix* make_plan_matrix(struct ndMatrix* const input, nduint N, nduint rend){
    struct ndMatrix* plan = (struct ndMatrix*)malloc(sizeof(struct ndMatrix));
    plan->R = input->R;
    plan->C = 1 + rend*N;
    
    plan->data = (T*)malloc(sizeof(T)*plan->C*plan->R);
    
    
        
    for (nduint j=0; j<input->R; j++) {
        nduint col = 1;
        for (nduint k=0; k<N; k++) { // Minden parameterre.
            
            T param = E(input, j, k);
            E(plan, j, 0) = 1.00;
            
            for (nduint r=1; r<=rend; r++) {
                T e = param;
                for (nduint h=1; h<r; h++) {
                    e*=param;
                }
                E(plan, j, col) = e;
                col++;
            }
        }
    }
    
    //Tervmatrix sorait elosztani a szigmakkal
    for (nduint i=0; i<plan->R; i++) {
        T sigma = E(input, i, N+1);
        for (nduint j=0; j<plan->C; j++) {
            E(plan, i, j) = E(plan, i, j)/sigma;
        }
    }
    return plan;
}

void solve(struct ndMatrix* const input, nduint N, nduint rend){
    struct ndMatrix* plan = make_plan_matrix(input, N, rend);
    
    NDDEBUGF("\n Tervmatrix:", *plan)
    
    struct ndMatrix* planT = transpose(plan);
    
    NDDEBUGF("\n Tervmatrix transzponalt:", *planT)
    
    struct ndMatrix* planTplan = product(planT, plan);
    
    NDDEBUGF("\n XTX:", *planTplan);
    
    //free(plan);
    
    //print_to_file(plan, "terv.dat");
    //print_to_file(planTplan, "invertalni.dat");
    
    struct ndMatrix* inverse = (struct ndMatrix*)malloc(sizeof(struct ndMatrix));
    nduint* colflips = (nduint*)calloc(planTplan->C, sizeof(nduint));
    for (int k=0; k<planTplan->C; k++) {
        colflips[k]=k;
    }
    
    inverse = solveInverse(planTplan, &colflips);
    free(planTplan);
    
    NDDEBUGF("\n Az inverz: \n", *inverse);
    struct ndMatrix* result = product(inverse, planT);
    
    free(planT);
    free(inverse);
    
    NDDEBUGF("\n Eredmeny matrix: \n", *result);
    
    //Most az y vektort meg kell szorozni a result-al, ez lesz a parameter tomb.
    T* y = (T*)malloc(sizeof(T)*result->C);
    for (int k=0; k<result->C; k++) {
        T sigma = E(input, k, N+1);
        *(y+k) = E(input, k, N)/sigma; // El kell osztani a sigmaval
    }
    
    T* a = dot_product(result, y);
    
    printf("\n Eredmeny: \n a= \n");
    for(int k=0; k<N*rend+1; k++) {
        printf("%lf\n", *(a+k));
    }
    
    print_results("results.dat", "params.dat",  input, plan, N, rend, a);
    
    free(plan);
    free(result);
}

int main(int argc, const char * argv[]) {
    
    if (argc == 4) {
        nduint N = atoi(argv[1]);
        nduint rend = atoi(argv[2]);
        
        struct ndMatrix* input = (struct ndMatrix*)malloc(sizeof(struct ndMatrix));
        input->C = N + 2; //Ennyi oszlopa van a bemenetnek
        read_input(input, N ,argv[3]);
        
        NDDEBUGF("\nBemeneti mátrix: \n", *input)
        
        solve(input, N, rend);
    }else{
        printf("Kell 2 parameter: N, polinom rendje, filename");
    }
    return 0;
}
