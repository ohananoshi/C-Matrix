 /*
C MATRIX

Author: Guilherme Arruda

Github: https://github.com/ohananoshi/C-Matrix

Created in: 09/05/23

Last updated: 08/nov/24

*/


//=============================== HEADERS ==============================

#ifndef _C_MATRIX_HEADER_
    #define _C_MATRIX_HEADER_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdarg.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

//=================================== CONSTANTS ===================================

enum error_flag{
    DIVIDE_BY_0,
    ARRAY_WITH_LENGTH_0,
    NULL_ARRAY,
    LENGTH_MISSMATCH,
    INVALID_PARAMETER
}error_flag;

enum AXIS{
    X,Y,Z
}AXIS;

enum FILL_MODE{
    BY_ROW,
    BY_COLUMN
}FILL_MODE;

//=================================== DATATYPES ===================================

typedef struct MATRIX_t{
    uint32_t rows;
    uint32_t columns;
    uint32_t layers;
    double ***data;
}MATRIX_t;

typedef struct slice_t{
    uint32_t start;
    uint32_t end;
    int32_t step;;
}slice_t;

//================================== COMMON FUNCTIONS ====================================================

int d_compare(const void *a, const void *b) {
    return (*((double*)a) > *((double*)b)) ? 1:-1;
}

//============================= FORMAT FUNCTIONS (FMT) ===============================================

slice_t* fmt_read(const char* fmt){
    slice_t *select = (slice_t*)calloc(1, sizeof(slice_t));
    sscanf(fmt,"[%u:%u:%d]%*s", &select->start, &select->end, &select->step);
    return select;
};

int8_t fmt_verify(const char* function_name , const uint32_t src_att, const slice_t* select, const char* att_name){

    if(select->start > src_att){
        fprintf(stderr, "ERROR:\n\tfunction: %s\n\tparameter: fmt\n\tmessage: START %s MUST BE <= SRC->%sS.\n", function_name, att_name, att_name);
        return LENGTH_MISSMATCH;
    }
    if(select->end > src_att){
        fprintf(stderr, "ERROR:\n\tfunction: %s\n\tparameter: fmt\n\tmessage: END %s MUST BE <= SRC->%sS.\n", function_name, att_name, att_name);
        return LENGTH_MISSMATCH;
    }
    if(select->step < 0) if((double)(-select->step) > (double)src_att){
        fprintf(stderr, "ERROR:\n\tfunction: %s\n\tparameter: fmt\n\tmessage: STEP MUST BE <= SRC->%sS.\n", function_name, att_name);
        return LENGTH_MISSMATCH;
    }else if((double)select->step > (double)src_att){
        fprintf(stderr, "ERROR:\n\tfunction: %s\n\tparameter: fmt\n\tmessage: STEP MUST BE <= SRC->%sS.\n", function_name, att_name);
        return LENGTH_MISSMATCH;
    }

    return 0;
}

//=============================== MATRIX FUNCTIONS =============================================

void matrix_init(MATRIX_t** dest, uint32_t rows, uint32_t columns, uint32_t layers){
    
    if(rows == 0){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_init\n\tparameter: rows\n\tmessage: rows cannot be 0.\n");
        exit(LENGTH_MISSMATCH);
    }
    if(columns == 0){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_init\n\tparameter: columns\n\tmessage: columns cannot be 0.\n");
        exit(LENGTH_MISSMATCH);
    }
    if(layers == 0){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_init\n\tparameter: layers\n\tmessage: layers cannot be 0.\n");
        exit(LENGTH_MISSMATCH);
    }
    
    (*dest) = (MATRIX_t*)calloc(1,sizeof(MATRIX_t));
    (*dest)->columns = columns;
    (*dest)->rows = rows;
    (*dest)->layers = layers;
    (*dest)->data = (double***)calloc(layers, sizeof(double**));

    for(uint32_t i = 0; i < layers; i++){
        (*dest)->data[i] = (double**)calloc(rows, sizeof(double*));
        for(uint32_t j = 0; j < rows; j++){
            (*dest)->data[i][j] = (double*)calloc(columns, sizeof(double*));
        }
    }
}

void matrix_free(MATRIX_t** matrix){
    if(matrix == NULL || *matrix == NULL) return;

    if((*matrix)->data != NULL){
        for(uint32_t i = 0; i < (*matrix)->layers; i++){
            for(uint32_t j = 0; j < (*matrix)->rows; j++){
                free((*matrix)->data[i][j]);
            }
            free((*matrix)->data[i]);
        }
        free((*matrix)->data);
    }
    
    free((*matrix));

    matrix = NULL;
}

void matrix_vfree(uint32_t count, ...){
    va_list arg;
    va_start(arg, count);

    for(uint32_t i = 0; i < count; i++){
        matrix_free(va_arg(arg, MATRIX_t**));
    }

    va_end(arg);
}

void matrix_print(const MATRIX_t* src){

    for(uint32_t i = 0; i < src->layers; i++){
        printf("layer: %d\n", i);
        for(uint32_t j = 0; j < src->rows; j++){
            for(uint32_t k = 0; k < src->columns; k++){
                printf("%lf ", src->data[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }

}

void matrix_info(const MATRIX_t* src){
    printf("LAYERS: %d\nROWS: %d\nCOLUMNS: %d\nSIZE: %lu bytes\n", src->layers, src->rows, src->columns, src->layers * src->rows * src->columns * sizeof(double));
}

void matrix_fill(MATRIX_t** dest, double x, const char* layer_fmt, const char* rows_fmt, const char* columns_fmt){
    
    if((*dest) == NULL || dest == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_fill\n\tparameter: dest\n\tmessage: dest CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    } 

    slice_t *l_slice = (slice_t*)calloc(1, sizeof(slice_t));
    slice_t *r_slice = (slice_t*)calloc(1, sizeof(slice_t));
    slice_t *c_slice = (slice_t*)calloc(1, sizeof(slice_t));

    if(layer_fmt == NULL){
        l_slice->start = 0;
        l_slice->end = (*dest)->layers;
        l_slice->step = 1;
    }else l_slice = fmt_read(layer_fmt);

    if(fmt_verify("matrix_fill", (*dest)->layers, l_slice, "LAYER")){
        free(l_slice);  
        free(r_slice);  
        free(c_slice); 
        exit(LENGTH_MISSMATCH);
    }

    if(rows_fmt == NULL){
        r_slice->start = 0;
        r_slice->end = (*dest)->rows;
        r_slice->step = 1;
    }else r_slice = fmt_read(rows_fmt);

    if(fmt_verify("matrix_fill", (*dest)->rows, r_slice, "ROW")){
        free(l_slice);  
        free(r_slice);  
        free(c_slice); 
        exit(LENGTH_MISSMATCH);
    }

    if(columns_fmt == NULL){
        c_slice->start = 0;
        c_slice->end = (*dest)->columns;
        c_slice->step = 1;
    }else c_slice = fmt_read(columns_fmt);

    if(fmt_verify("matrix_fill", (*dest)->columns, c_slice, "COLUMN")){
        free(l_slice);  
        free(r_slice);  
        free(c_slice); 
        exit(LENGTH_MISSMATCH);
    }

    for(uint32_t i = l_slice->start; i < l_slice->end; i += l_slice->step){
        for(uint32_t j = r_slice->start; j < r_slice->end; j += r_slice->step){
            for(uint32_t k = c_slice->start; k < c_slice->end; k += c_slice->step){
                (*dest)->data[i][j][k] = x;
            }
        }
    } 

    free(l_slice);  
    free(r_slice);  
    free(c_slice); 
}

void matrix_rotate(MATRIX_t** src, enum AXIS axis, const char* axis_fmt){
    if((*src) == NULL || src == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_rotate\n\tparameter: src\n\tmessage: SRC CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    } 

    slice_t *ax_slice = NULL;
    slice_t *ax_slice2 = NULL;

    double tmp;

    switch (axis)
    {
    case X:
    {   
        ax_slice = (slice_t*)calloc(1, sizeof(slice_t));

        if(axis_fmt == NULL){
            ax_slice->start = 0;
            ax_slice->end = (*src)->columns;
            ax_slice->step = 1;
        }else{
            ax_slice = fmt_read(axis_fmt);
            if(fmt_verify("matrix_verify", (*src)->columns, ax_slice, "COLUMN")) exit(LENGTH_MISSMATCH);
        }

        for(uint32_t i = 0; i < (*src)->layers; i++){
            for(uint32_t k = ax_slice->start; k < ax_slice->end; k += ax_slice->step){
                for(uint32_t j = 0; j < (*src)->rows/2; j++){
                    tmp = (*src)->data[i][j][k];
                    (*src)->data[i][j][k] = (*src)->data[i][(*src)->rows - j - 1][k];
                    (*src)->data[i][(*src)->rows - j - 1][k] = tmp;
                }
            }
        }

        free(ax_slice);
    }
    break;
    case Y:
    {   
        ax_slice = (slice_t*)calloc(1, sizeof(slice_t));

        if(axis_fmt == NULL){
            ax_slice->start = 0;
            ax_slice->end = (*src)->rows;
            ax_slice->step = 1;
        }else{
            ax_slice = fmt_read(axis_fmt);
            if(fmt_verify("matrix_verify", (*src)->columns, ax_slice, "ROW")) exit(LENGTH_MISSMATCH);
        }

        for(uint32_t i = 0; i < (*src)->layers; i++){
            for(uint32_t j = ax_slice->start; j < ax_slice->end; j += ax_slice->step){
                for(uint32_t k = 0; k < (*src)->columns/2; k++){
                    tmp = (*src)->data[i][j][k];
                    (*src)->data[i][j][k] = (*src)->data[i][j][(*src)->columns - k - 1];
                    (*src)->data[i][j][(*src)->columns - k - 1] = tmp;
                }
            }
        }

        free(ax_slice);
    }
    break;
    case Z:
    {   
        ax_slice = (slice_t*)calloc(1, sizeof(slice_t));
        ax_slice2 = (slice_t*)calloc(1, sizeof(slice_t));

        if(axis_fmt == NULL){
            ax_slice->start = 0;
            ax_slice->end = (*src)->rows;
            ax_slice->step = 1;

            ax_slice2->start = 0;
            ax_slice2->end = (*src)->columns;
            ax_slice2->step = 1;
        }else{
            ax_slice = fmt_read(axis_fmt);
            ax_slice2 = fmt_read(strchr(axis_fmt, ','));
            if(fmt_verify("matrix_verify", (*src)->columns, ax_slice2, "COLUMN")) exit(LENGTH_MISSMATCH);
            if(fmt_verify("matrix_verify", (*src)->rows, ax_slice, "ROW")) exit(LENGTH_MISSMATCH);
        }

        for(uint32_t i = 0; i < (*src)->layers/2; i++){
            for(uint32_t j = ax_slice->start; j < ax_slice->end; j += ax_slice->step){
                for(uint32_t k = ax_slice2->start; k < ax_slice2->end; k += ax_slice2->step){
                    tmp = (*src)->data[i][j][k];
                    (*src)->data[i][j][k] = (*src)->data[(*src)->layers - i - 1][j][k];
                    (*src)->data[(*src)->layers - i - 1][j][k] = tmp;
                }
            }
        }

        free(ax_slice);
        free(ax_slice2);
    }
    break;
    
    default:
        fprintf(stderr, "ERROR:\n\tfunction: matrix_rotate\n\tparameter: axis\n\tmessage: axis %hu is invalid.\n", axis);
        exit(INVALID_PARAMETER);
        break;
    }

}

void matrix_copy(const MATRIX_t* src, MATRIX_t** dest, const char* layer_fmt, const char* row_fmt, const char* column_fmt, bool same_shape){

    if(src == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_copy\n\tparameter: src\n\tmessage: SRC CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    }

    slice_t *l_slice = (slice_t*)calloc(1, sizeof(slice_t));
    slice_t *r_slice = (slice_t*)calloc(1, sizeof(slice_t));
    slice_t *c_slice = (slice_t*)calloc(1, sizeof(slice_t));
    
    if(layer_fmt == NULL){
        l_slice = (slice_t*)calloc(1, sizeof(slice_t));

        l_slice->start = 0;
        l_slice->end = src->layers;
        l_slice->step = 1;
    }else l_slice = fmt_read(layer_fmt);

    if(fmt_verify("matrix_fill", src->layers, l_slice, "LAYER")){
        free(l_slice);  
        free(r_slice); 
        free(c_slice);  
        exit(LENGTH_MISSMATCH);
    }

    if(row_fmt == NULL){
        r_slice = (slice_t*)calloc(1, sizeof(slice_t));

        r_slice->start = 0;
        r_slice->end = src->rows;
        r_slice->step = 1;
    }else r_slice = fmt_read(row_fmt);

    if(fmt_verify("matrix_fill", src->rows, r_slice, "ROW")){
        free(l_slice);  
        free(r_slice); 
        free(c_slice);   
        exit(LENGTH_MISSMATCH);
    }

    if(column_fmt == NULL){
        c_slice = (slice_t*)calloc(1, sizeof(slice_t));

        c_slice->start = 0;
        c_slice->end = src->columns;
        c_slice->step = 1;
    }else c_slice = fmt_read(column_fmt);

    if(fmt_verify("matrix_fill", src->columns, c_slice, "COLUMN")){
        free(l_slice);  
        free(r_slice); 
        free(c_slice); 
        exit(LENGTH_MISSMATCH);
    }

    if(same_shape){
        if(dest == NULL || (*dest) == NULL) matrix_init(dest, src->rows, src->columns, src->layers);

        for(uint32_t i = l_slice->start; i < l_slice->end; i += l_slice->step){
            for(uint32_t j = r_slice->start; j < r_slice->end; j += r_slice->step){
                for(uint32_t k = c_slice->start; k < c_slice->end; k += c_slice->step){
                    (*dest)->data[i][j][k] = src->data[i][j][k];
                }
            }
        }
    }else{
        if(dest == NULL || (*dest) == NULL){
            matrix_init(dest, ((r_slice->end - r_slice->start - 1)/r_slice->step)+1, ((c_slice->end - c_slice->start - 1)/c_slice->step)+1, ((l_slice->end - l_slice->start - 1)/l_slice->step)+1);
        }

        for(uint32_t i = l_slice->start, i2 = 0; i < l_slice->end; i += l_slice->step, i2++){
            for(uint32_t j = r_slice->start, j2 = 0; j < r_slice->end; j += r_slice->step, j2++){
                for(uint32_t k = c_slice->start, k2 = 0; k < c_slice->end; k += c_slice->step, k2++){
                    (*dest)->data[i2][j2][k2] = src->data[i][j][k];
                }
            }
        }
    }

    free(l_slice);  
    free(r_slice);  
    free(c_slice); 
}

void matrix_reshape(const MATRIX_t* src, MATRIX_t** dest,uint32_t new_rows, uint32_t new_columns, uint32_t new_layers, enum FILL_MODE fill_mode){
    if(src == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_reshape\n\tparameter: src\n\tmessage: SRC CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    }
    
    if((new_rows*new_columns*new_layers) != ((src->columns) * (src->rows) * (src->layers))){
        fprintf(stderr,"ERROR:\n\tfunction: matrix_reshape\n\tparamter: new_rows | new_columns | new_layers\n\tmessage: MATRIX RESHAPE DIMENSIONS NOT MATCH");
        exit(LENGTH_MISSMATCH);
    }

    matrix_init(dest, new_rows, new_columns, new_layers);
    uint32_t idx_destino = 0;
    
    switch(fill_mode){
    
        case BY_COLUMN:
            for (uint32_t k = 0; k < src->layers; k++) {
                for (uint32_t i = 0; i < src->columns; i++) {
                    for (uint32_t j = 0; j < src->rows; j++) {
                        uint32_t index = k * (src->rows * src->columns) + i * src->rows + j;
                        uint32_t new_k = index / (new_rows * new_columns);
                        uint32_t new_i = (index % (new_rows * new_columns)) / new_rows;
                        uint32_t new_j = (index % (new_rows * new_columns)) % new_rows;

                        (*dest)->data[new_k][new_j][new_i] = src->data[k][j][i];
                    }
                }
            }
            break;
        case BY_ROW:

            for (uint32_t i = 0; i < src->layers; i++) {
                for (uint32_t j = 0; j < src->rows; j++) {
                    for (uint32_t k = 0; k < src->columns; k++) {
                        uint32_t linha_destino = idx_destino / (new_rows * new_columns);
                        uint32_t coluna_destino = (idx_destino / new_columns) % new_rows;
                        uint32_t profundidade_destino = idx_destino % new_columns;

                        (*dest)->data[linha_destino][coluna_destino][profundidade_destino] = src->data[i][j][k];
                        idx_destino++;
                    }
                }
            }

            break;

        default:
            fprintf(stderr, "ERROR:\n\tfunction: matrix_reshape\n\tparameter: fill_mode\n\tmessage: fill_mode %hu does not exist.\n", fill_mode);
            exit(INVALID_PARAMETER);
    }
    
}

void matrix_transpose(const MATRIX_t* src, MATRIX_t** dest, enum AXIS axis){
    
    if(src == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_transpose\n\tparameter: src\n\tmessage: SRC CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    }

    if(axis == Y){
        if(src->columns != src->layers) matrix_init(dest, src->rows, src->layers, src->columns);
        else matrix_init(dest, src->rows, src->columns, src->layers);
        
        for(uint32_t i = 0; i < src->layers; i++){
            for(uint32_t j = 0; j < src->rows; j++){
                for(uint32_t k = 0; k < src->columns; k++){
                    (*dest)->data[k][j][i] = src->data[i][j][k];
                }
            }
        }
    }else if(axis == X){
        if(src->rows != src->layers) matrix_init(dest, src->layers, src->columns, src->rows);
        else matrix_init(dest, src->rows, src->columns, src->layers);
        
        for(uint32_t i = 0; i < src->layers; i++){
            for(uint32_t j = 0; j < src->rows; j++){
                for(uint32_t k = 0; k < src->columns; k++){
                    (*dest)->data[j][i][k] = src->data[i][j][k];
                }
            }
        }
    }else if(axis == Z){
        if(src->rows != src->columns) matrix_init(dest, src->columns, src->rows, src->layers);
        else matrix_init(dest, src->rows, src->columns, src->layers);
        
        for(uint32_t i = 0; i < src->layers; i++){
            for(uint32_t j = 0; j < src->rows; j++){
                for(uint32_t k = 0; k < src->columns; k++){
                    (*dest)->data[i][k][j] = src->data[i][j][k];
                }
            }
        }
    }

}

void matrix_func(const MATRIX_t* src, MATRIX_t** dest, double func(double), const char* layer_fmt, const char* row_fmt, const char* column_fmt, bool same_shape){
    if(src == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_func\n\tparameter: src\n\tmessage: SRC CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    }

    slice_t *l_slice, *r_slice, *c_slice;
    
    if(layer_fmt == NULL){
        l_slice->start = 0;
        l_slice->end = (*dest)->layers;
        l_slice->step = 1;
    }else l_slice = fmt_read(layer_fmt);

    if(fmt_verify("matrix_fill", (*dest)->layers, l_slice, "LAYER")){
        free(l_slice);  
        free(r_slice);  
        free(c_slice); 
        exit(LENGTH_MISSMATCH);
    }

    if(row_fmt == NULL){
        r_slice->start = 0;
        r_slice->end = (*dest)->rows;
        r_slice->step = 1;
    }else r_slice = fmt_read(row_fmt);

    if(fmt_verify("matrix_fill", (*dest)->rows, r_slice, "ROW")){
        free(l_slice);  
        free(r_slice);  
        free(c_slice); 
        exit(LENGTH_MISSMATCH);
    }

    if(column_fmt == NULL){
        c_slice->start = 0;
        c_slice->end = (*dest)->columns;
        c_slice->step = 1;
    }else c_slice = fmt_read(column_fmt);

    if(same_shape){
        if(dest == NULL || (*dest) == NULL) matrix_init(dest, src->rows, src->columns, src->layers);

        for(uint32_t i = l_slice->start; i < l_slice->end; i += l_slice->step){
            for(uint32_t j = r_slice->start; j < r_slice->end; j += r_slice->step){
                for(uint32_t k = c_slice->start; k < c_slice->end; k += c_slice->step){
                    (*dest)->data[i][j][k] = func(src->data[i][j][k]);
                }
            }
        }
    }else{
        if(dest == NULL || (*dest) == NULL){
            matrix_init(dest, ((r_slice->end - r_slice->start - 1)/r_slice->step)+1, ((c_slice->end - c_slice->start - 1)/c_slice->step)+1, ((l_slice->end - l_slice->start - 1)/l_slice->step)+1);
        }

        for(uint32_t i = l_slice->start, i2 = 0; i < l_slice->end; i += l_slice->step, i2++){
            for(uint32_t j = r_slice->start, j2 = 0; j < r_slice->end; j += r_slice->step, j2++){
                for(uint32_t k = c_slice->start, k2 = 0; k < c_slice->end; k += c_slice->step, k2++){
                    (*dest)->data[i2][j2][k2] = func(src->data[i][j][k]);
                }
            }
        }
    }

    free(l_slice);  
    free(r_slice);  
    free(c_slice);
}

int8_t matrix_compare(const void* a, const void* b){
    if(a == NULL || b == NULL){
        fprintf(stderr, "ERROR:\n\tfunction: matrix_compare\n\tparameter: a || b\n\tmessage: a || b CAN'T BE NULL\n");
        exit(NULL_ARRAY);
    }

    MATRIX_t *x = (MATRIX_t*)a, *y = (MATRIX_t*)b;

    if(x->layers != y->layers || x->rows != y->rows || x->columns != y->columns) return -1;
    for(uint32_t i = 0; i < x->layers; i++){
        for(uint32_t j = 0; j < x->rows; j++){
            for(uint32_t k = 0; k < x->columns; k++){
                if(x->data[i][j][k] != y->data[i][j][k]) return 1;
            }
        }
    }

    return 0;
}

void matrix_sort(const MATRIX_t* src, MATRIX_t** dest){
    MATRIX_t *tmp;
    double *dtmp = (double*)calloc(src->layers*src->rows*src->columns, sizeof(double));
    
    matrix_reshape(src, &tmp, 1, src->layers*src->rows*src->columns, 1, 1);

    for(uint64_t i = 0; i < src->layers*src->rows*src->columns; i++) dtmp[i] = tmp->data[0][0][i];

    qsort(dtmp, src->layers*src->rows*src->columns, sizeof(double), d_compare);
   
    for(uint64_t i = 0; i < src->layers*src->rows*src->columns; i++) tmp->data[0][0][i] = dtmp[i];

    matrix_reshape(tmp, dest, src->rows, src->columns, src->layers, 2);

    free(tmp);
    free(dtmp);
}

#endif //END OF C MATRIX HEADER