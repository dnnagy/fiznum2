//
//  main.c
//  bead2
//
//  Created by Nagy Daniel on 28/02/2017.
//  Copyright © 2017 Nagy Daniel. All rights reserved.
//
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#define T double
#define EPSILON 1e-9 // Ennel kisebb elemeket nem engedunk a foatloba !!
#define nduint unsigned int
#define SWAPD(U,X,Y) {U aux = (X); (X) = (Y); (Y) = aux;}
#define ABS(X) ((X) > 0 ? (X) : (-(X)))

typedef nduint* ndArray; // olyan vektor, amiben az oszlopok vagy sorok poziciojat taroljuk

struct ndMatrix{
    T* data;
    nduint R, C;
};


#define NDDEBUG(Y) printf("\nLine: %d:", __LINE__); print_mat((Y));
#define E(M, I, J)  *((M)->data + (I)*(M)->C + (J))  //Dereferalunk egy elemet a matrixbol

// felcs. az i-edik es a j-edik sort
void swap_rows(struct ndMatrix* const mat, nduint i, nduint j){
    for (int k = 0; k < mat->C; k++) {
        SWAPD(T, E(mat, i, k) , E(mat, j, k) );
    }
}

// d-edik sorhoz hozzaadjuk az s-edik sor a-szorosat.
void add_mult_row(struct ndMatrix* const mat, nduint s, nduint d, T a){
    for (nduint i = 0; i < mat->C; i++) {
        E(mat, d, i) += (E(mat, s, i))*a;
    }
}

//Kivonja az i-edik sorbol a j-ediket
void diffRow(struct ndMatrix* const mat, nduint i, nduint j){
    add_mult_row(mat, j, i, -1);
}

// felcs. az i-edik es a j-edik oszlopot , recordFlips=true => feljegyezzük a sorcserét
void swap_cols(struct ndMatrix* const mat, nduint i, nduint j, ndArray cf, unsigned int recordFlips){
    //for (int k = 0; k < mat->R*mat->C; k+=mat->C) {
    //    SWAPD(T, *(mat->data + k + i), *(mat->data + k + j));
    //}
    
    for (int r = 0; r<mat->R; r++) {
        SWAPD(T, E(mat, r, i), E(mat, r, j));
    }
    
    if (recordFlips != 0) {
        SWAPD(nduint, cf[i], cf[j]);
    }
}


//Kibovitjuk a colflips tombot egy elemmel
void expand(unsigned int** a, int N){
    unsigned int* tmp = malloc((N+1)*sizeof(unsigned int));
    
    for (int k = 0; k<N; k++) {
        *(tmp + k) = *( (*a) + k );
    }
    *(tmp + N) = N;
    
    free(*a);
    *a = tmp;
    return;
};


//Hozzáadunk a mátrixhoz egy oszlopot, a jobb oldalra.
void addCol(struct ndMatrix* const mt, T* col, ndArray* ptrToCf){
    
    expand(ptrToCf, mt->C);
    
    T* newData = (T*)malloc(sizeof(T)*(mt->C+1)*(mt->R));
    
    int k=0;
    for (nduint i=0; i < mt->R; i++) {
        for (nduint j=0; j < mt->C; j++) {
           *(newData + i*mt->C + j + k) = *(mt->data + i*mt->C + j);
        }
        *(newData + (i+1)*mt->C + k) = *(col + i);
        k++;
    }
    
    mt->C++;
    free(mt->data);
    mt->data = newData;
}

//hozzaadunk egy egysegmatrixot. ez felteszi, hogy csak NxN-es a matrix
void addIdentity(struct ndMatrix* const mt, ndArray* cfArray){
    T* col = (T*)calloc(mt->R, sizeof(T));
    
    for (int i=0; i<mt->R; i++) {
        *(col+i)=1;
        addCol(mt, col, cfArray);
        *(col+i) = 0; //vissza kell allitani
    }
}

//megnezzuk, hogy jol vannak-e az oszlopok?
nduint checkColOrder(ndArray cf, nduint C){
    int k=0;
    while (k<C && cf[k] == k) {
        k++;
    }
    return k==C;
}

//oszlopok visszacserelese az eredeti sorrendben
void reorderCols(struct ndMatrix* const mat, ndArray cf){
    
    while (checkColOrder(cf, mat->C) == 0) {
        for(int k=0; k<mat->C; k++){
            swap_cols(mat, cf[cf[k]], cf[k], cf, 1);
        }
    }
    
}

// az r-edik sort beszorozza a-val
void mult_row(struct ndMatrix* const mat, nduint r, T a){
    for (nduint i = mat->C*r; i < mat->C*r + mat->C; i++) {
        *(mat->data + i) *= a;
    }
}

//megkeresi az row-adik sorban az abszolut ertekben legnagyobb elemet, ez a retval-adik oszlopban lesz, es visszaadja
T find_pivot_in_row(struct ndMatrix mat, nduint row, nduint* const retval){
    int k = 0;
    T max = ABS(E(&mat, row, 0)); //ABS(*(mat.data + row*mat.C));
    T elem = E(&mat, row, k);
    while ( k < mat.C ){
        elem = E(&mat, row, k);
        if ( ABS(elem) >= max) {
            max = ABS(elem);
            *retval = k;
        }
        k++;
    }
    return E(&mat, row, *retval);
}

//Legnagyobb abszolut erteku elem megtalalasa oszlopban, foatlo alatt
T find_abs_max_col(struct ndMatrix mat, nduint col, nduint* const retval){
    nduint k=col;
    T max = E(&mat, k, col);//ABS(*(mat.data + mat.C*col + col));
    T elem = E(&mat, k, col);
    while (k < mat.R) {
        elem = E(&mat, k, col);
        if( ABS(elem) > max ){
            max = ABS(elem);
            *retval = k;
        }
        k++;
    }
    return E(&mat, *retval, col);
}

//Matrix beolvasas
void read_matrix(struct ndMatrix* const mt, const char* filename){
    mt->data = (T*)malloc(mt->C*mt->R*sizeof(T));
    
    FILE* in = fopen(filename, "r");
    for (nduint i=0; i< mt->R; i++) {
        for (nduint j=0; j< mt->C; j++) {
          fscanf(in, "%lf", (mt->data) + i*mt->C + j);
        }
    }
}

//print matrix for debugging
void print_mat(struct ndMatrix m){
    printf("\n");
    for (int i=0; i < m.R ; i++) {
        for (int j=0; j < m.C; j++) {
            printf("%lf ", *(m.data + i*m.C + j));
        }
        printf("\n");
    }
}

///Fo algoritmus matrix inverz meghatarozasara
struct ndMatrix* solveInverse(struct ndMatrix* const mt, ndArray* cf){
    
    assert(mt->C == mt->R);
    assert(mt->data != NULL);
    nduint N=mt->R; //eredetileg NxN-es matrixunk volt
    
    //Eloszor kiegeszitjuk egy egysegmatrixszal
    addIdentity(mt, cf);
    
    // NDDEBUG(*mt)
    
    //Elozetes pivotolas
    nduint i = 0;
    while (i != N) {
        //Pivotalas
        unsigned int k = 0;
        find_pivot_in_row(*mt, i, &k);
        mult_row(mt, i, 1.0/E(mt, i, k));// NDDEBUG(*mt)
        i++;
    }
    
    //Vegigmegyunk a foatlon
    nduint idx = 0;
    while (idx != N) {
        
        //printf("\nFoatlo beli elem M[%d][%d] = %lf\n", idx, idx, E(mt, idx,idx));
        
        // Ha a foatlo beli elem < Epsilon, akkor meg kell keresni az oszlopban a legnagyobb elemet, es fel kell cserelni a 2 sort
        if ( ABS(E(mt, idx, idx)) < EPSILON ) {
            
            nduint k = 0;
            //Meg kell nezni, hogy a foatlo alatt legyen legalabb 1 nem nulla elem
            if ( find_abs_max_col(*mt, idx, &k) < EPSILON ) {
                exit(1);
            }else{
                swap_rows(mt, idx, k);// printf("Megcsereltuk a %d es %d sorokat:", idx, k); NDDEBUG(*mt)
            }
        }
        
        //Kinullazzuk az adott oszlopot a foatlo alatt
        nduint row = idx+1; //kovetkezo sortol indulunk
        T szorzo;
        while ( row < mt->R ) {
            szorzo = -1*(E(mt, row, idx))/(E(mt, idx, idx));
            add_mult_row(mt, idx, row, szorzo );// printf("Hozzaadjuk a %d. sorhoz a %d. sort szorozva %lf", row, idx, szorzo); NDDEBUG(*mt)
            row++;
        }
        idx++;
    }
    
    //Ha a foatlo alatt kinullaztunk mindent, akkor elkezdjuk folotte is kinullazni
    //Ehhez vegigmegyunk a foatlon a masodik elemtol kezdve
    idx = 1;
    while (idx != N) {
        
        //Kinullazzuk az adott oszlopot a foatlo felett
        nduint row = idx; //ettol a sortol indulunk
        T szorzo;
        while ( row > 0 ) {
            szorzo = -1*(E(mt, row-1, idx))/(E(mt, idx, idx));
            add_mult_row(mt, idx, row-1, szorzo );// NDDEBUG(*mt)
            row--;
        }
        idx++; //NDDEBUG(*mt)
    }

    //Minden sort beallitunk ugy, hogy a bal oldalon 1 legyen
    idx=0;
    while ( idx < N ) {
        T szorzo = 1.00/E(mt, idx, idx);
        mult_row(mt, idx, szorzo);
        idx++;
    }
    
    //Visszaadjuk az inverz matrixot
    struct ndMatrix* inv = (struct ndMatrix*)malloc(sizeof(struct ndMatrix));
    inv->C = N;
    inv->R = N;
    
    inv->data = (T*)malloc(sizeof(T)*N*N);
    
    //Bemasoljuk az inverz matrixot
    for (int i=0; i < inv->R ; i++) {
        for (int j=0; j < inv->C; j++) {
           E(inv, i, j) = E(mt, i, j + N);
        }
    }
    
    return inv;
}

////////////////////
///Teszt fuggvenyek
////////////////////
void testSwapRows(struct ndMatrix* const m){
    printf("Swapping rows %d and %d\n", 0, 1);
    swap_rows(m, 0, 1);
    printf("Swapping rows %d and %d\n", 2, 3);
    swap_rows(m, 2, 3);
    print_mat(*m);
}

//Oszlop csere teszt
void testSwapcCols(struct ndMatrix* const m, ndArray* cf){
    swap_cols(m, 0, 2, *cf, 1);
    printf("\nOsszecserelt oszlopok: \n" );
    
    print_mat(*m);
    printf("\nViszzarendezzuk: \n");
    
    reorderCols(m, *cf);
    print_mat(*m);
}

//oszlop hozzaadas teszt
void testAddCol(struct ndMatrix* const m, ndArray* cf){
    T col[4] = {55.0, 55.0, 55.0, 55.0};
    
    printf("oszlop hozzaadasa...\n");
    addCol(m, col, cf);
    print_mat(*m);
}

void testFindMax(struct ndMatrix m){
    nduint k = 0;
    printf("Legnagyobb elem: %lf\n", find_pivot_in_row(m, 1, &k));
}


int main(int argc, const char * argv[]) {

    if (argc == 3) {
        
        struct ndMatrix mtrx; //ez a matrix
        
        ndArray colflips = NULL; //Itt taroljuk, hogy melyik oszlopokat csereltuk fel
        

        nduint N=atoi(argv[1]);
        
        mtrx.R = N;
        mtrx.C = N;
        read_matrix(&mtrx, argv[2]);

        //Ezekben taroljuk a sor es az oszlop csereket
        colflips = (nduint*)calloc(mtrx.C, sizeof(nduint));
        for (int k=0; k<mtrx.C; k++) {
            colflips[k]=k;
        }
        
        printf("Eredeti mátrix: \n");
        print_mat(mtrx);
        printf("\n");
       
        //testSwapRows(&mtrx);
        //testAddCol(&mtrx, &colflips);
        //testSwapcCols(&mtrx, &colflips);
        
        //Debugoljuk a derefet
        //printf("M(0,0): %f\n", DEREF(&mtrx,0,0));
        //printf("M(1,1): %f\n", DEREF(&mtrx,1,1));
        //printf("M(3,2): %f\n", DEREF(&mtrx,3,2));
        
        
        struct ndMatrix* inverz = solveInverse(&mtrx, &colflips);
        printf("Inverz matrix (sorokat cserelve): \n");
        print_mat(*inverz);
    
    }else{
        printf("HIBA! Ket argumentum kell: N, filename");
    }

    return 0;
}
