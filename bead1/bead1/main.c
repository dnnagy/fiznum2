//
//  main.c
//  bead1
//
//  Created by Nagy Daniel on 16/02/2017.
//  Copyright © 2017 Nagy Daniel. All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include "string.h"
#include <assert.h>

// Linuxon 0-ra kell allitani
#define FIRST_ARG_IDX 1 //Mac-en az argv[0] a fullpath, linuxon az elso option

struct vector{
    float num;
    struct vector* prev;
    struct vector* next;
};

struct vector* vecRoot = NULL;

unsigned long int vsize(struct vector* root){
    unsigned long int vs = 0;
    struct vector* i = root;
    while (i != NULL) {
        i = i->next;
        vs++;
    }
    return vs;
}

float dotProduct(struct vector* v1, struct vector* v2){
    assert(vsize(v1) == vsize(v2));
    struct vector* i = v1;
    struct vector* j = v2;
    float dp = 0;
    while (i != NULL && j != NULL) {
        dp += i->num*j->num;
        i = i->next;
        j = j->next;
    }
    return dp;
};

void push_back(struct vector* root, struct vector* item){
    struct vector* current = root;
    while (current->next != NULL) {
        current = current->next;
    }
    current->next = item;
    item->prev = current;
}

void read_vec(const char* filename){
    vecRoot = (struct vector*) malloc(sizeof(struct vector));
    vecRoot->prev=NULL;
    vecRoot->next=NULL;
    
    FILE* in = fopen(filename, "r");
    float k;
    fscanf(in, "%f", &vecRoot->num);
    while (fscanf(in,"%f", &k) == 1) {
        struct vector* v = (struct vector*) malloc(sizeof(struct vector));
        v->num = k;
        v->next = NULL;
        push_back(vecRoot, v);
    }
}

struct matrix{
    struct matrix* next;
    struct vector* row;
};

struct matrix* matRoot = NULL;

void new_row(struct matrix* root, struct matrix* newItem){
    struct matrix* c = root;
    while (c->next != NULL && c->row != NULL) {
        c = c->next;
    }
    c->next = newItem;
}

struct vector* make_row(char* string){
    struct vector* root = (struct vector*)malloc(sizeof(struct vector));
    root->next=NULL;
    root->prev=NULL;
    char* pch;
    pch = strtok (string," \t\n");
    float f = 0;
    while (pch != NULL)
    {
        struct vector* v = (struct vector*)malloc(sizeof(struct vector));
        f = atof(pch);
        v->num = f;
        v->next = NULL;
        push_back(root, v);
        pch = strtok (NULL, " \t\n");
    }
    return root;
}

void read_matrix(const char* filename){
    matRoot = (struct matrix*)malloc(sizeof(struct matrix));
    matRoot->next = NULL;
    matRoot->row = (struct vector*)malloc(sizeof(struct vector));
    matRoot->row->next = NULL;
    matRoot->row->prev = NULL;
    matRoot->row->num = 0;

    FILE* in = fopen(filename, "r");
    
    char* line = NULL;
    size_t len = 0;
    ssize_t read;
    
    int count = 0;
    while ((read = getline(&line, &len, in)) != -1) {
        struct matrix* row = (struct matrix*)malloc(sizeof(struct matrix));
        row->next = NULL;
        row->row = make_row(line);//(struct vector*)malloc(sizeof(struct vector));
        new_row(matRoot, row);
        count++;
    }
    printf("number of lines read: %d", count);
}

unsigned long numofRows(struct matrix* root){
   unsigned long ret = 0;
    struct matrix* i = root;
    while (i->next != NULL) {
        i = i->next;
        ret++;
    }
    return ret;
}

//Csak debug
//void print_vector(struct vector* root){
//    struct vector* i = root;
//    while (i != NULL) {
//        printf("%f ", i->num);
//        i=i->next;
//    }
//}
//void print_matrix(struct matrix* root){
//    struct matrix* i = root->next;
//    while (i != NULL) {
//        print_vector(i->row->next);
//        printf("\n");
//        i = i->next;
//    }
//}


void solve(struct matrix* matrix, struct vector* vector){
    struct matrix* row = matrix;
    while (row != NULL) {
        printf("%f\n", dotProduct(row->row->next, vector));
        row = row->next;
    }
}

void argParse(int argc, const char * argv[]){
    if (argc == FIRST_ARG_IDX + 2) {
        read_matrix(argv[FIRST_ARG_IDX]);
        read_vec(argv[FIRST_ARG_IDX + 1]);
        
        solve(matRoot->next, vecRoot);
        
    }else{
        printf("Nem jo a bemenet!");
    }
}

int main(int argc, const char * argv[]) {
    argParse(argc, argv);
    return 0;
}
