 /*
C MATRIX

Author: Guilherme Arruda

Github: https://github.com/ohananoshi/C-Matrix

Created in: 09/05/23

Last updated: 14/06/23

*/


//=============================== HEADERS ==============================
#ifndef EMBEDDED_MODE
    #pragma once
    #include <stdio.h>
    #include <stdlib.h>
    #include <stdint.h>
    #include <stdarg.h>
    #include <limits.h>
    #include <math.h>
#endif
#ifdef EMBEDDED_MODE
    #ifndef _INC_STDIO
        #include <stdio.h>
    #endif
    #ifndef _INC_STDLIB
        #include <stdlib.h>
    #endif
    #ifndef _STDINT_H
        #include <stdint.h>
    #endif
    #ifndef _STDARG_H
        #include <stdarg.h>
    #endif
    #ifndef _GCC_LIMITS_H_
        #include <limits.h>
    #endif
    #ifndef _MATH_H_
        #include <math.h>
    #endif
#endif

#ifndef _C_MATRIX_HEADER_
    #define _C_MATRIX_HEADER_

//=================================== CONSTANTS ===================================

#define FILLED 1
#define NOT_FILLED 0 

//=================================== DATATYPES ===================================

enum iterate_modes{
    ALL = 0,
    BY_ROW,
    BY_COLUMN,
    BY_SLICE,
    CUSTOM
};

enum return_modes{
    INPLACE = -1,
    NEW = -2
};

enum operations{
    SUM = 0,
    SUB,
    MULT,
    DIVIDE,
    POWER
};

enum extract_mode{
    KEEP_ROWS_SIZE = 0,
    KEEP_COLUMNS_SIZE,
    KEEP_SLICES_SIZE,
    CUT
};

enum expand_mode{
    ADD_ROWS = 0,
    ADD_COLUMNS,
    ADD_SLICES
};

typedef struct{
    uint32_t rows;
    uint32_t columns;
    uint32_t slices;
    double ***data;
}MATRIX;

//================================== FUNCTIONS ====================================
//math functions
static inline double divide(double x, double y){
    if(y == 0.000000000000000){
        fprintf(stderr, "DIVISION BY 0");
        return 0;
    }
    return x/y;
}
static inline double sum(double x, double y){
    return x+y;
}
static inline double sub(double x, double y){
    return x-y;
}
static inline double mult(double x, double y){
    return x*y;
}

//bool functions
static inline double is_bigger(double x, double y){ return x > y ? 1.0:0.0;}
static inline double is_smaller(double x, double y){ return x < y ? 1.0:0.0;}
static inline double is_different(double x, double y){ return x != y ? 1.0:0.0;}
static inline double is_equal(double x, double y){ return x == y ? 1.0:0.0;}

//array functions
double* arr(uint32_t size, ...){
    double* buffer = (double*)calloc(size, sizeof(double));

    va_list num;
    va_start(num, size);
        
    for(uint32_t i = 0; i < size; i++){
        buffer[i] = va_arg(num,double);
    }

    va_end(num);

    return buffer;
}

int* int_arr(uint32_t size, ...){
    int* buffer = (int*)calloc(size, sizeof(int));

    va_list num;
    va_start(num, size);

    for(uint32_t i = 0; i < size; i++){
        buffer[i] = va_arg(num,int);
    }

    va_end(num);

    return buffer;
}

int32_t* int_a_b(int32_t start, int32_t end, int32_t step){
    if(step == 0){
        fprintf(stderr,"ERROR: DIVISION BY ZERO");
        return NULL;
    }
    uint32_t size = ((end-start+1)/step) + 1;
    int32_t *buffer = (int32_t*)calloc(size, sizeof(int32_t));

    for(int32_t i = 0; i < size; i++){
        buffer[i] = start + i*step;
    }

    return buffer;
}

double* double_linspace(double start, double end, uint32_t size){
    if(size == 0){
        fprintf(stderr,"ERROR: ARRAY LENGTH CAN'T BE ZERO");
        return NULL;
    }

    double *buffer = (double*)calloc(size, sizeof(double));
    if(size == 1){
        buffer[0] = start;
        return buffer;
    }

    double step = (end - start)/((double)(size - 1));

    printf("\n%d \n", size);
    for(int32_t i = 0; i < size; i++){
        buffer[i] = start + ((double)i)*step;
    }

    return buffer;
}

double array_dot(const double* a, const double* b, uint32_t size){
    if(size == 0){
        fprintf(stderr,"ERROR: ARRAY LENGTH CAN'T BE ZERO");
        return NAN;
    }
    double result = 0.0;

    for(uint32_t i = 0; i < size; i++){
        result += a[i]*b[i];
    }

    return result;
}

//matrix functions
void matrix_free(MATRIX *matrix){
    if(matrix == NULL) return;

    for(uint32_t i = 0; i < matrix->slices; i++){
        for(uint32_t j = 0; j < matrix->rows; j++){
            free(matrix->data[i][j]);
        }
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);

    matrix = NULL;
}

void matrix_vfree(uint32_t count, ...){
    va_list matrix;
    va_start(matrix, count);

    for(uint32_t i = 0; i < count; i++){
        matrix_free(va_arg(matrix,MATRIX*));
    }
}

void matrix_print(const MATRIX a){

    for(uint32_t i = 0; i < a.slices; i++){
        printf("slice: %d\n", i);
        for(uint32_t j = 0; j < a.rows; j++){
            //printf("linha %d\n",j);
            for(uint32_t k = 0; k < a.columns; k++){
                printf("%f ", a.data[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }

}

void matrix_info(const MATRIX matrix){
    printf("SLICES: %d\nROWS: %d\nCOLUMNS: %d\n", matrix.slices, matrix.rows, matrix.columns);
}

MATRIX* matrix_init(uint32_t rows, uint32_t columns, uint32_t slices, uint8_t filled, ...){
    MATRIX* matrix = (MATRIX*)calloc(1,sizeof(MATRIX));
    matrix->columns = columns;
    matrix->rows = rows;
    matrix->slices = slices;
    matrix->data = (double***)calloc(slices, sizeof(double**));
    
    if(filled){
        va_list array;
        va_start(array, filled);
        double* buffer = va_arg(array, double*);
        va_end(array);

        for(uint32_t i = 0, count = 0; i < slices; i++){
            matrix->data[i] = (double**)calloc(rows, sizeof(double*));
            for(uint32_t j = 0; j < rows; j++){
                matrix->data[i][j] = (double*)calloc(columns, sizeof(double));
                for(uint32_t k = 0; k < columns; k++, count++){
                    matrix->data[i][j][k] = buffer[count];
                }
            }
        }
        free(buffer);

        return matrix;
    }

    for(uint32_t i = 0; i < slices; i++){
        matrix->data[i] = (double**)calloc(rows, sizeof(double*));
        for(uint32_t j = 0; j < rows; j++){
            matrix->data[i][j] = (double*)calloc(columns, sizeof(double));
        }
    }

    return matrix;
}

MATRIX* matrix_ones(uint32_t rows, uint32_t columns, uint32_t slices){
    MATRIX *ones = matrix_init(rows, columns, slices, NOT_FILLED);

    for(uint32_t i = 0; i < slices; i++){
        for(uint32_t j = 0; j < rows; j++){
            for(uint32_t k = 0; k < columns; k++){
                ones->data[i][j][k] = 1.0;
            }
        }
    }

    return ones;
}

void matrix_transfer(const MATRIX a, MATRIX *b){

    for(uint32_t i = 0; i < a.slices; i++){
        for(uint32_t j = 0; j < a.rows; j++){
            for(uint32_t k = 0; k < a.columns; k++){
                b->data[i][j][k] = a.data[i][j][k];
            }
        }
    }
}

MATRIX* matrix_copy(const MATRIX input_matrix){
    MATRIX* copy = matrix_init(input_matrix.rows, input_matrix.columns, input_matrix.slices, NOT_FILLED);

    for(uint32_t i = 0; i < copy->slices; i++){
        for(uint32_t j = 0; j < copy->rows; j++){
            for(uint32_t k = 0; k < copy->columns; k++){
                copy->data[i][j][k] = input_matrix.data[i][j][k];
            }
        }
    }

    return copy;
}

MATRIX* matrix_transpose(MATRIX** input_matrix, int8_t return_mode, uint8_t iterate_mode,...){
    
    if(iterate_mode == ALL){
        MATRIX *transpose = matrix_init((*input_matrix)->columns, (*input_matrix)->rows, (*input_matrix)->slices, NOT_FILLED);

        for(uint32_t i = 0; i < (*input_matrix)->slices; i++){
            for(uint32_t j = 0; j < (*input_matrix)->rows; j++){
                for(uint32_t k = 0; k < (*input_matrix)->columns; k++){
                    transpose->data[i][k][j] = (*input_matrix)->data[i][j][k];
                }
            }
        }

        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return transpose;
    }
    if(iterate_mode == CUSTOM){
        va_list info;
        va_start(info, iterate_mode);

        int count = va_arg(info, int);
        if(count > (*input_matrix)->slices){
            fprintf(stderr, "ERROR: INCOMPATIBLE NUMBER OF SLICES.");
            return NULL;
        }

        int* slices = va_arg(info, int*);
        va_end(info);

        MATRIX *transpose = matrix_init((*input_matrix)->columns, (*input_matrix)->rows, count, NOT_FILLED);

        for(uint16_t i = 0; i < count; i++){
            for(uint32_t i = 0; i < (*input_matrix)->rows; i++){
                for(uint32_t j = 0; j < (*input_matrix)->columns; j++){
                    if(slices[i] > (*input_matrix)->slices){
                        fprintf(stderr, "ERROR: SLICE NUMBER DOESN'T EXIST.");
                        return NULL;
                    }
                    transpose->data[slices[i]][j][i] = (*input_matrix)->data[slices[i]][i][j];
                }
            }
        }

        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return transpose;
    }

    fprintf(stderr,"ERROR: ITERATE MODE DOESN'T EXIST.");
    return NULL;
}

MATRIX* matrix_reshape(MATRIX** input_matrix, uint32_t new_rows, uint32_t new_columns, uint32_t new_slices, int8_t return_mode){
    if((new_rows*new_columns*new_slices) != (((*input_matrix)->columns) * ((*input_matrix)->rows) * ((*input_matrix)->slices))){
        fprintf(stderr,"MATRIX RESHAPE DIMENSIONS NOT MATCH");
        return NULL;
    }

    MATRIX* reshaped = matrix_init(new_rows, new_columns, new_slices, NOT_FILLED);
    
    for (uint32_t k = 0; k < (*input_matrix)->slices; k++) {
        for (uint32_t i = 0; i < (*input_matrix)->columns; i++) {
            for (uint32_t j = 0; j < (*input_matrix)->rows; j++) {
                uint32_t index = k * ((*input_matrix)->rows * (*input_matrix)->columns) + i * (*input_matrix)->rows + j;
                uint32_t new_k = index / (new_rows * new_columns);
                uint32_t new_i = (index % (new_rows * new_columns)) / new_rows;
                uint32_t new_j = (index % (new_rows * new_columns)) % new_rows;

                reshaped->data[new_k][new_j][new_i] = (*input_matrix)->data[k][j][i];
            }
        }
    }

    if(return_mode == INPLACE) matrix_free(*input_matrix);

    return reshaped;
}

MATRIX* matrix_expand(MATRIX** input_matrix, uint8_t expand_mode, uint16_t count, int8_t return_mode){
    if(expand_mode == ADD_COLUMNS){
        if((count + (*input_matrix)->columns) > UINT32_MAX){
            fprintf(stderr, "EXPANSION EXCEED MATRIX COLUMNS LENGTH LIMIT.");
            return NULL;
        }

        MATRIX* expanded = matrix_init((*input_matrix)->rows, (*input_matrix)->columns + count, (*input_matrix)->slices, NOT_FILLED);
        matrix_transfer(*(*input_matrix), expanded);

        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return expanded;
    }
    if(expand_mode == ADD_ROWS){
        if((count + (*input_matrix)->rows) > UINT32_MAX){
            fprintf(stderr, "EXPANSION EXCEED MATRIX ROWS LENGTH LIMIT.");
            return NULL;
        }

        MATRIX* expanded = matrix_init((*input_matrix)->rows + count, (*input_matrix)->columns, (*input_matrix)->slices, NOT_FILLED);
        matrix_transfer(*(*input_matrix), expanded);

        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return expanded; 
    }
    if(expand_mode == ADD_SLICES){
        if((count + (*input_matrix)->slices) > UINT32_MAX){
            fprintf(stderr, "EXPANSION EXCEED MATRIX SLICES LENGTH LIMIT.");
            return NULL;
        }

        MATRIX* expanded = matrix_init((*input_matrix)->rows, (*input_matrix)->columns, (*input_matrix)->slices + count, NOT_FILLED);
        matrix_transfer(*(*input_matrix), expanded);

        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return expanded;
    }
    
    fprintf(stderr, "ERROR: EXPAND MODE DOESN'T EXIST.");
    return NULL;
}

MATRIX* matrix_pop(MATRIX** input_matrix, int8_t return_mode, uint8_t extract_mode, uint8_t iterate_mode, ...){

    if(iterate_mode == BY_ROW){
        va_list info;
        va_start(info, iterate_mode);

        int size = va_arg(info, int);
        int* rows = va_arg(info, int*);
        va_end(info);

        MATRIX *filtered = matrix_init(size, (*input_matrix)->columns, (*input_matrix)->slices, NOT_FILLED);

        for(uint32_t i = 0; i < filtered->slices; i++){
            for(uint32_t j = 0; j < size; j++){
                for(uint32_t k = 0; k < filtered->columns; k++){
                    filtered->data[i][j][k] = (*input_matrix)->data[i][rows[j]][k];
                }
            }
        }

        free(rows);
        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return filtered;
    }
    if(iterate_mode == BY_COLUMN){
        va_list info;
        va_start(info, iterate_mode);

        int size = va_arg(info, int);
        int* columns = va_arg(info, int*);
        va_end(info);

        MATRIX *filtered = matrix_init((*input_matrix)->rows, size, (*input_matrix)->slices, NOT_FILLED);

        for(uint32_t i = 0; i < filtered->slices; i++){
            for(uint32_t j = 0; j < filtered->rows; j++){
                for(uint32_t k = 0; k < size; k++){
                    filtered->data[i][j][k] = (*input_matrix)->data[i][j][columns[k]];
                }
            }
        }

        free(columns);
        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return filtered;
    }
    if(iterate_mode == BY_SLICE){
        va_list info;
        va_start(info, iterate_mode);

        int size = va_arg(info, int);
        int* slices = va_arg(info, int*);
        va_end(info);

        MATRIX *filtered = matrix_init((*input_matrix)->rows, (*input_matrix)->columns, size, NOT_FILLED);

        for(uint32_t i = 0; i < size; i++){
            for(uint32_t j = 0; j < filtered->rows; j++){
                for(uint32_t k = 0; k < filtered->columns; k++){
                    filtered->data[i][j][k] = (*input_matrix)->data[slices[i]][j][k];
                }
            }
        }

        free(slices);
        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return filtered;
    }
    if(iterate_mode == CUSTOM){
        va_list info;
        va_start(info, iterate_mode);

        int* sizes = va_arg(info, int*);
        int* rows = va_arg(info, int*);
        int* columns = va_arg(info, int*);
        int* slices = va_arg(info, int*);
        va_end(info);

        if(extract_mode == KEEP_COLUMNS_SIZE){
            MATRIX *filtered = matrix_init(sizes[0], (*input_matrix)->columns, sizes[2], NOT_FILLED);
            for(uint32_t i = 0; i < sizes[2]; i++){
                for(uint32_t j = 0; j < sizes[0]; j++){
                    for(uint32_t k = 0; k < sizes[1]; k++){
                        filtered->data[i][j][columns[k]] = (*input_matrix)->data[slices[i]][rows[j]][columns[k]];
                    }
                }
            }

            free(slices);
            free(rows);
            free(columns);
            if(return_mode == INPLACE) matrix_free(*input_matrix);
            return filtered;
        }
        if(extract_mode == CUT){
            MATRIX *filtered = matrix_init(sizes[0], sizes[1], sizes[2], NOT_FILLED);
            for(uint32_t i = 0; i < sizes[2]; i++){
                for(uint32_t j = 0; j < sizes[0]; j++){
                    for(uint32_t k = 0; k < sizes[1]; k++){
                        filtered->data[i][j][k] = (*input_matrix)->data[slices[i]][rows[j]][columns[k]];
                    }
                }
            }

            free(slices);
            free(rows);
            free(columns);
            if(return_mode == INPLACE) matrix_free(*input_matrix);
            return filtered;
        }


    }

    return NULL;
}

MATRIX* matrix_add(MATRIX** a, const MATRIX b, int8_t return_mode, uint8_t iterate_mode, ...){

    MATRIX *sum = matrix_init((*a)->rows, (*a)->columns, (*a)->slices, NOT_FILLED);

    if(iterate_mode == ALL){
        for(uint32_t i = 0; i < (*a)->slices; i++){
            for(uint32_t j = 0; j < (*a)->rows; j++){
                for(uint32_t k = 0; k < (*a)->columns; k++){
                    sum->data[i][j][k] = (*a)->data[i][j][k] + b.data[i][j][k];
                }
            }
        }
    }
    if(iterate_mode == BY_ROW){
        va_list info;
        va_start(info, iterate_mode);

        int* sizes = va_arg(info, int*);
        int* rows = va_arg(info, int*);

        if(sizes[0] == -1){
            va_end(info);
            for(uint32_t i = 0; i < sum->slices; i++){
                for(uint32_t j = 0; j < sizes[1]; j++){
                    for(uint32_t k = 0; k < sum->columns; k++){
                        sum->data[i][rows[j]][k] = (*a)->data[i][rows[j]][k] + b.data[i][rows[j]][k];
                    }
                }
            }
            free(rows);
        }
        else{
            int* slices = va_arg(info, int*);
            va_end(info);
            for(uint32_t i = 0; i < sizes[0]; i++){
                for(uint32_t j = 0; j < sizes[1]; j++){
                    for(uint32_t k = 0; k < sum->columns; k++){
                        sum->data[slices[i]][rows[j]][k] = (*a)->data[slices[i]][rows[j]][k] + b.data[slices[i]][rows[j]][k];
                        //printf("slice: %d linha: %d coluna: %d -> %f\n", slices[i], rows[j], k, (*a)->data[slices[i]][rows[j]][k]);
                    }
                }
            }
            free(rows);
            free(slices);
        }
    }
    if(iterate_mode == BY_COLUMN){
        va_list info;
        va_start(info, iterate_mode);

        int* sizes = va_arg(info, int*);
        int* columns = va_arg(info, int*);

        if(sizes[0] == -1){
            va_end(info);
            for(uint32_t i = 0; i < sum->slices; i++){
                for(uint32_t j = 0; j < sum->rows; j++){
                    for(uint32_t k = 0; k < sizes[1]; k++){
                        sum->data[i][j][columns[k]] = (*a)->data[i][j][columns[k]] + b.data[i][j][columns[k]];
                    }
                }
            }
            free(columns);
        }
        else{
            int* slices = va_arg(info, int*);
            va_end(info);
            for(uint32_t i = 0; i < sizes[0]; i++){
                for(uint32_t j = 0; j < sum->rows; j++){
                    for(uint32_t k = 0; k < sizes[1]; k++){
                        sum->data[slices[i]][j][columns[k]] = (*a)->data[slices[i]][j][columns[k]] + b.data[slices[i]][j][columns[k]];
                    }
                }
            }
            free(columns);
            free(slices);
        }
    }
    if(iterate_mode == CUSTOM){
        va_list info;
        va_start(info, iterate_mode);
        
        int* sizes = va_arg(info, int*);
        int* rows = va_arg(info, int*);
        int* columns = va_arg(info, int*);
        int* slices = va_arg(info, int*);
        va_end(info);

        //printf("slice : %d linha: %d coluna: %d\n", sizes[0], sizes[1], sizes[2]);
        for(uint32_t i = 0; i < sizes[0]; i++){
            for(uint32_t j = 0; j < sizes[1]; j++){
                for(uint32_t k = 0; k < sizes[2]; k++){
                    sum->data[slices[i]][rows[j]][columns[k]] = (*a)->data[slices[i]][rows[j]][columns[k]] + b.data[slices[i]][rows[j]][columns[k]];
                    //printf("slice : %d linha: %d coluna: %d %d valor: %f \n",slices[i], rows[j], columns[k] ,k,sum->data[slices[i]][rows[j]][columns[k]]);
                }
            }
        }
        free(slices);
        free(sizes);
        free(rows);
        free(columns);
    }
    //else{
    //    fprintf(stderr, "ITERATE MODE DOESN'T EXIST");
    //    return NULL;
    //}
    if(return_mode == INPLACE) matrix_free(*a);
    //printf("-->aqui\n");
    //matrix_print(*sum);
    return sum;
}

MATRIX* matrix_product(MATRIX** a, const MATRIX b, int8_t return_mode, uint8_t iterate_mode, ...){
    if((*a)->columns != b.rows){
        fprintf(stderr, "ERROR: MATRIX PRODUCT DIMENSIONS NOT MATCH");
        return NULL;
    }

    if(iterate_mode == ALL){
        if(b.slices != (*a)->slices){
            fprintf(stderr, "ERROR: DIFFERENT NUMBER OF SLICES.");
            return NULL;
        }

        MATRIX *product = matrix_init((*a)->rows, b.columns, b.slices, NOT_FILLED);
        double *buffer = (double*)calloc((*a)->columns, sizeof(double));

        for(uint32_t k = 0; k < b.slices; k++){
            for(uint32_t i = 0; i < b.columns; i++){
                for(uint32_t j = 0; j < b.rows; j++){
                    buffer[j] = b.data[k][j][i];
                }
                for(uint32_t j = 0; j < (*a)->rows; j++){
                    product->data[k][j][i] = array_dot(buffer, (*a)->data[k][j], (*a)->columns);
                }
            }
        }
        free(buffer);

        if(return_mode == INPLACE) matrix_free(*a);
        return product;
    }
    if(iterate_mode == CUSTOM){
        if(b.slices != (*a)->slices){
            fprintf(stderr, "ERROR: DIFFERENT NUMBER OF SLICES.");
            return NULL;
        }

        va_list info;
        va_start(info, iterate_mode);

        int slice_count = va_arg(info, int);
        int* a_slices = va_arg(info, int*);
        int* b_slices = va_arg(info, int*);
        va_end(info);

        MATRIX *product = matrix_init((*a)->rows, b.columns, slice_count, NOT_FILLED);
        double *buffer = (double*)calloc((*a)->columns, sizeof(double));

        for(uint32_t k = 0; k < slice_count; k++){
            for(uint32_t i = 0; i < b.columns; i++){
                for(uint32_t j = 0; j < b.rows; j++){
                    buffer[j] = b.data[b_slices[k]][j][i];
                }
                for(uint32_t j = 0; j < (*a)->rows; j++){
                    product->data[k][j][i] = array_dot(buffer, (*a)->data[a_slices[k]][j], (*a)->columns);
                }
            }
        }
        free(buffer);

        if(return_mode == INPLACE) matrix_free(*a);
        return product;
    }

    fprintf(stderr, "ERROR: ITERATE MODE DOESN'T EXIST.");
    return NULL;
}

MATRIX* matrix_sum(MATRIX** input_matrix, int8_t iterate_mode, int8_t return_mode){
    
    if(iterate_mode == ALL){
        MATRIX* sum = matrix_init(1,1,1, NOT_FILLED);

        for(uint32_t i = 0; i < (*input_matrix)->slices; i++){
            for(uint32_t j = 0; j < (*input_matrix)->rows; j++){
                for(uint32_t k = 0; k < (*input_matrix)->columns; k++){
                    sum->data[0][0][0] += (*input_matrix)->data[i][j][k];
                }
            }
        }

        if(return_mode == INPLACE) matrix_free(*input_matrix);

        return sum;
    }
    if(iterate_mode == BY_COLUMN){
        MATRIX* sum = matrix_init(1,(*input_matrix)->columns,(*input_matrix)->slices, NOT_FILLED);

        for(uint32_t i = 0; i < (*input_matrix)->slices; i++){
            for(uint32_t j = 0; j < (*input_matrix)->columns; j++){
                for(uint32_t k = 0; k < (*input_matrix)->rows; k++){
                    sum->data[i][0][j] += (*input_matrix)->data[i][k][j];
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return sum;
    }
    if(iterate_mode == BY_ROW){
        MATRIX* sum = matrix_init((*input_matrix)->rows,1,(*input_matrix)->slices, NOT_FILLED);

        for(uint32_t i = 0; i < (*input_matrix)->slices; i++){
            for(uint32_t j = 0; j < (*input_matrix)->rows; j++){
                for(uint32_t k = 0; k < (*input_matrix)->columns; k++){
                    sum->data[i][j][0] += (*input_matrix)->data[i][j][k];
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*input_matrix);
        return sum;
    }
    if(iterate_mode == BY_SLICE){
        MATRIX* sum = matrix_init((*input_matrix)->rows,(*input_matrix)->columns,1, NOT_FILLED);

        for(uint32_t i = 0; i < (*input_matrix)->rows; i++){
            for(uint32_t j = 0; j < (*input_matrix)->columns; j++){
                for(uint32_t k = 0; k < (*input_matrix)->slices; k++){
                    sum->data[0][i][j] += (*input_matrix)->data[k][i][j];
                    //printf("slice : %d linha: %d coluna: %d %d valor: %f\n",i,j,k ,sum->data[0][i][j]);
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*input_matrix);
        //matrix_info(*sum);
        //matrix_print(*sum);
        return sum;
    }
    
    fprintf(stderr, "ITERATE MODE DOESN'T EXIST\n");
    return NULL;

}

MATRIX* matrix_1by1(MATRIX** a, const MATRIX b, uint8_t operation, int return_mode){
    if(((*a)->rows != b.rows) || ((*a)->columns != b.columns) || ((*a)->slices) != b.slices){
        fprintf(stderr, "MATRIX 1 BY 1 OPERATION ERROR: DIMENSIONS NOT MATCH\n");
        return NULL;
    };
    
    MATRIX* operand = matrix_init((*a)->rows, (*a)->columns, (*a)->slices, NOT_FILLED);

    if(operation == SUM){
        for(uint32_t i = 0; i < operand->slices; i++){
            for(uint32_t j = 0; j < operand->rows; j++){
                for(uint32_t k = 0; k < operand->columns; k++){
                    operand->data[i][j][k] = (*a)->data[i][j][k] + b.data[i][j][k];
                }
            }
        }

        if(return_mode == INPLACE) matrix_free(*a);
        return operand;
    }
    if(operation == SUB){
        for(uint32_t i = 0; i < operand->slices; i++){
            for(uint32_t j = 0; j < operand->rows; j++){
                for(uint32_t k = 0; k < operand->columns; k++){
                    operand->data[i][j][k] = (*a)->data[i][j][k] - b.data[i][j][k];
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*a);
        return operand;
    }
    if(operation == MULT){
        for(uint32_t i = 0; i < operand->slices; i++){
            for(uint32_t j = 0; j < operand->rows; j++){
                for(uint32_t k = 0; k < operand->columns; k++){
                    operand->data[i][j][k] = (*a)->data[i][j][k] * b.data[i][j][k];
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*a);
    }
    if(operation == DIVIDE){
        for(uint32_t i = 0; i < operand->slices; i++){
            for(uint32_t j = 0; j < operand->rows; j++){
                for(uint32_t k = 0; k < operand->columns; k++){
                    if(b.data[i][j][k] != 0) operand->data[i][j][k] = (*a)->data[i][j][k] / b.data[i][j][k];
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*a);
    }
    if(operation == POWER){
        for(uint32_t i = 0; i < operand->slices; i++){
            for(uint32_t j = 0; j < operand->rows; j++){
                for(uint32_t k = 0; k < operand->columns; k++){
                    operand->data[i][j][k] = pow((*a)->data[i][j][k], b.data[i][j][k]);
                }
            }
        }
        if(return_mode == INPLACE) matrix_free(*a);
    }

    return operand;
}

MATRIX* matrix_iterate(MATRIX** input_matrix, double (*action)(double, double), double y, int8_t return_mode, uint8_t iterate_mode, ...){

    MATRIX *output = matrix_copy(*(*input_matrix));

    if(iterate_mode == ALL){
        for(uint32_t i = 0; i < output->slices; i++){
            for(uint32_t j = 0; j < output->rows; j++){
                for(uint32_t k = 0; k < output->columns; k++){
                    output->data[i][j][k] = action((*input_matrix)->data[i][j][k], y);
                }
            }
        }
    }
    if(iterate_mode == BY_ROW){
        va_list info;
        va_start(info, iterate_mode);

        int* sizes = va_arg(info, int*);
        int* rows = va_arg(info, int*);
        int* slices = va_arg(info, int*);

        va_end(info);

        if(sizes[0] == -1){
            for(uint32_t i = 0; i < output->slices; i++){
                for(uint32_t j = 0; j < sizes[1]; j++){
                    for(uint32_t k = 0; k < output->columns; k++){
                        output->data[i][rows[j]][k] = action((*input_matrix)->data[i][rows[j]][k], y);
                    }
                }
            }
        }
        else{
            for(uint32_t i = 0; i < sizes[0]; i++){
                for(uint32_t j = 0; j < sizes[1]; j++){
                    for(uint32_t k = 0; k < output->columns; k++){
                        output->data[slices[i]][rows[j]][k] = action((*input_matrix)->data[slices[i]][rows[j]][k], y);
                    }
                }
            }
        }

        free(sizes);
        free(slices);
        free(rows);
        free(action);

    }
    if(iterate_mode == BY_COLUMN){
        va_list info;
        va_start(info, iterate_mode);

        int* sizes = va_arg(info, int*);
        int* columns = va_arg(info,int*);

        if(sizes[0] == -1){
            va_end(info);
            for(uint32_t i = 0; i < output->slices; i++){
                for(uint32_t j = 0; j < output->rows; j++){
                    for(uint32_t k = 0; k < sizes[1]; k++){
                        output->data[i][j][columns[k]] = action((*input_matrix)->data[i][j][columns[k]], y);
                    }
                }
            }
        }
        else{
            int* slices = va_arg(info, int*);
            va_end(info);
            
            for(uint32_t i = 0; i < sizes[0]; i++){
                for(uint32_t j = 0; j < output->rows; j++){
                    for(uint32_t k = 0; k < sizes[1]; k++){
                        output->data[slices[i]][j][columns[k]] = action((*input_matrix)->data[slices[i]][j][columns[k]], y);
                    }
                }
            }
            free(slices);
        }

        free(sizes);
        free(columns);
    }
    if(iterate_mode == CUSTOM){
        va_list info;
        va_start(info, iterate_mode);
        
        int* sizes = va_arg(info, int*);
        int* rows = va_arg(info, int*);
        int* columns = va_arg(info, int*);

        //printf("slice : %d linha: %d coluna: %d\n", sizes[0], sizes[1], sizes[2]);
        
        if(sizes[0] == -1){
            //va_end(info);
            for(uint32_t i = 0; i < output->slices; i++){
                for(uint32_t j = 0; j < sizes[1]; j++){
                    for(uint32_t k = 0; k < sizes[2]; k++){
                        output->data[i][rows[j]][columns[k]] = action((*input_matrix)->data[i][rows[j]][columns[k]], y);
                    }
                }
            }
            free(rows);
            free(columns);
        }
        else{
            int* slices = va_arg(info, int*);
            va_end(info);
            for(uint32_t i = 0; i < sizes[0]; i++){
                for(uint32_t j = 0; j < sizes[1]; j++){
                    for(uint32_t k = 0; k < sizes[2]; k++){
                        output->data[slices[i]][rows[j]][columns[k]] =  action((*input_matrix)->data[slices[i]][rows[j]][columns[k]], y);
                        //printf("slice : %d linha: %d coluna: %d %d valor: %f \n",slices[i], rows[j], columns[k] ,k,output->data[slices[i]][rows[j]][columns[k]]);
                    }
                }
            }
            free(slices);
        }
        free(sizes);
        free(rows);
        free(columns);
    }
    //else{
    //    fprintf(stderr, "ITERATE MODE DOESN'T EXIST");
    //    return NULL;
    //}

    if(return_mode == INPLACE) free(*input_matrix);

    return output;
}

MATRIX* array_to_matrix(uint32_t input_rows, uint32_t input_columns, uint32_t input_slices, ...){
    if(input_rows == 0 || input_columns == 0){
        fprintf(stderr, "\nERROR: DIMENSION 0 ARRAY IS NOT ALLOWED.\n");
        return NULL;
    }

    MATRIX* matrix;
    va_list array;
    va_start(array, input_slices);

    if(input_slices == 1 && input_rows == 1 && input_columns >= 1) matrix = matrix_init(1, input_columns, 1, FILLED, va_arg(array, double*));
    if(input_slices == 1 && input_rows > 1 && input_columns >= 1){
        matrix = matrix_init(input_rows, input_columns, 1, NOT_FILLED);
        double** arr = va_arg(array, double**);
        
        for(uint32_t i = 0; i < input_slices; i++){
            for(uint32_t j = 0; j < input_rows; j++){
                for(uint32_t k = 0; k < input_columns; k++){
                    matrix->data[i][j][k] = arr[j][k];
                }
            }
        }
    }
    if(input_slices > 1 && input_rows >= 1 && input_columns >= 1){
        matrix = matrix_init(input_rows, input_columns, input_slices, NOT_FILLED);
        double*** arr = va_arg(array, double***);

        for(uint32_t i = 0; i < input_slices; i++){
            for(uint32_t j = 0; j < input_rows; j++){
                for(uint32_t k = 0; k < input_columns; k++){
                    matrix->data[i][j][k] = arr[i][j][k];
                }
            }
        }
    }

    va_end(array);

    return matrix;
}

double*** matrix_to_array(MATRIX** a, uint8_t erase){
    double*** arr = (double***)calloc((uint32_t)(*a)->slices, sizeof(double**));
    
    for(uint32_t i = 0; i < (uint32_t)(*a)->slices; i++){
        arr[0] = (double**)calloc((uint32_t)(*a)->rows, sizeof(double*));

        for(uint32_t j = 0; j < (uint32_t)(*a)->rows; j++){
            arr[i][j] = (double*)calloc((uint32_t)(*a)->columns, sizeof(double));

            for(uint32_t k = 0; k < (uint32_t)(*a)->columns; k++){
                arr[i][j][k] = (*a)->data[i][j][k];
            }
        }
    }

    if(erase) matrix_free(*a);

    return arr;
}

#endif //END OF C MATRIX HEADER