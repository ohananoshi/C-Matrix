### C Matrix

##### Abstraction to matrices

#### Headers
-  stdio 
-  stlib
-  stdint
-  stdarg
-  string
-  math

#### DATATYPES

```C
//Matrix representation itself

typedef struct MATRIX_t{
    uint32_t rows;
    uint32_t columns;
    uint32_t slices;
    double ***data;
}MATRIX_t;

//Python like slicing to facilitate matrix traversal
typedef struct slice_t{
    uint32_t start;
    uint32_t end;
    int32_t step;;
}slice_t;
```

#### Error Flags

```C
DIVIDE_BY_0
ARRAY_WITH_LENGTH_0
NULL_ARRAY
LENGTH_MISSMATCH
INVALID_PARAMETER
```

#### AXIS

```C
X,Y,Z
```

#### FILL_MODE

```C
BY_ROW
BY_COLUMN
```

#### Functions

- #### Common Functions

    ```C
    int d_compare(const void *a, const void *b)
    
    Used as function pointer to qsort
    
    ```

    ```C
    slice_t* fmt_read(const char* fmt) 
    
    Python like slicing pattern.

    a   -> iteration start
    b-1 -> iteration ends
    c   -> step

    Reads [a:b:c] and converts into slice_t pointer
    
    ```

    ```C
    int8_t fmt_verify(const char* function_name , const uint32_t src_att, const slice_t* select, const char* att_name) 
    
    Verify if fmt string was correctly passed
    
    ```
    
- #### Matrix functions
    - #### matrix_init
    ```C
    void matrix_init(MATRIX_t** dest, uint32_t rows, uint32_t columns, uint32_t layers);

    Create a MATRIX_t* pointer dynamically allocated.

    ERRORS:
        if rows == 0 or columns == 0 or layers == 0
        function will exit with LENGTH_MISSMATCH.

    ```

    - #### matrix_free
    ```C
    void matrix_free(MATRIX_t** matrix);

    Free memory of matrix pointer.
    ```

    - #### matrix_vfree
    ```C
    void matrix_vfree(uint32_t count, ...);

    Free memory of a sequence of matrices.

    Parameters
        count: number of matrices.
        ...  : MATRIX_t** values.

    ```

    - #### matrix_print
    ```C
    void matrix_print(const MATRIX_t* src);

    Print matrix layer by layer.
    ```

    - #### matrix_info
    ```C
    void matrix_info(const MATRIX_t* src);

    Print number of layers, rows columns and size in bytes of matrix.

    ```

    - #### matrix_fill
    ```C
    void matrix_fill(MATRIX_t** dest, double x, const char* layer_fmt, const char* rows_fmt, const char* columns_fmt);

    
    Fills the matrix with the value indicated in 'x'.

    Parameters:
        dest        : matrix passed by reference.
        x           : fill value
        layer_fmt   : layer filling pattern "[a:b:c]"
            if layer_fmt == NULL -> all layers will be filled.

        rows_fmt    : row filling pattern "[a:b:c]"
            if rows_fmt == NULL -> all rows will be filled.

        columns_fmt : columns filling pattern "[a:b:c]"
            if columns_fmt == NULL -> all columns will be filled.

        
    ERRORS:
        exit with LENGTH_MISSMATCH if any fmt be bigger than dest dimensions.
        exit with NULL_ARRAY if dest == NULL

    ```
  
    - #### matrix_copy
    ```C
    void matrix_copy(const MATRIX_t* src, MATRIX_t** dest, const char* layer_fmt, const char* row_fmt, const char* column_fmt, bool same_shape);

    Copy matrices. Completely or by parts.

    Parameters:
        src : source matrix
        dest : dest matrix
        layer_fmt   : layer copy pattern "[a:b:c]"
            if layer_fmt == NULL -> all layers will be copied.

        rows_fmt    : row copy pattern "[a:b:c]"
            if rows_fmt == NULL -> all rows will be copied.

        columns_fmt : columns copy pattern "[a:b:c]"
            if columns_fmt == NULL -> all columns will be copied.

        same_shape  :
            if true : dest will receive data according to the fmt pattern
            if false: dest will have just enough len in each dimension to receive data

    ERRORS:
        LENGTH_MISSMATCH
        NULL_ARRAY
    ```
    
    - #### matrix_rotate

    ```C
    void matrix_rotate(MATRIX_t** src, enum AXIS axis, const char* axis_fmt);

    Flip values around axis without change dimensions.

    Parameters:
        src: source/dest matrix
        axis: rotation axis (X,Y,Z)
        axis_fmt: layers/rows/columns that will be rotated

    ERRORS:
        LENGTH_MISSMATCH
        NULL_ARRAY
    ```
  

    - #### matrix_reshape
    ```C
    void matrix_reshape(const MATRIX_t* src, MATRIX_t** dest,uint32_t new_rows, uint32_t new_columns, uint32_t new_layers, uint8_t fill_mode);

    Reshape the matrix while maintaining its volume.

    Parameters:
        src         : source matrix
        dest        : dest matrix
        new_rows    : new rows count
        new_columns : new columns count
        new_layers  : new layers count
        fill_mode   :
            BY_ROWS    : fills the new matrix row by row
            BY_COLUMNS : fills the new matrix column by column

    ERRORS:
        NULL_ARRAY
        LENGTH_MISSMATCH
        INVALID_PARAMETER
    ```

    - #### matrix_tranpose
    See [Matrix Transpose](https://en.wikipedia.org/wiki/Transpose)

    ```C
    void matrix_transpose(const MATRIX_t* src, MATRIX_t** dest, enum AXIS axis);

    Transpose the matrix.

    Parameters:
        src  : source matrix
        dest : dest matrix
        axis : rotation axis (X,Y,Z)

    ERRORS:
        NULL_ARRAY
    ```
    
    - #### matrix_func
    ```C
    void matrix_func(const MATRIX_t* src, MATRIX_t** dest, double func(double), const char* layer_fmt, const char* row_fmt, const char* column_fmt, bool same_shape)

    Apply a function to each cell in the matrix.

    Parameters:
        src         : source matrix
        dest        : dest matrix
        func        : function
        layer_fmt   : layer pattern "[a:b:c]"
            if layer_fmt == NULL -> all layers will be used.

        rows_fmt    : row pattern "[a:b:c]"
            if rows_fmt == NULL -> all rows will be used.

        columns_fmt : columns pattern "[a:b:c]"
            if columns_fmt == NULL -> all columns will be used.

        same_shape  :
            if true : dest will receive data according to the fmt pattern
            if false: dest will have just enough len in each dimension to receive data

    ERRORS:
        NULL_ARRAY
        LENGTH_MISSMATCH
    ```
    
    - #### matrix_compare
    ```C
    int8_t matrix_compare(const void* a, const void* b)

    Compares two matrices.

    return:
        dim(a) != dim(b) -> -1
        cell_ijk(a) != cell_ijk(b) -> 1
        dim(a) == dim(b) and cell_ijk(a) == cell_ijk(b) -> 0

    ERRORS:
        NULL_ARRAY
    ```
   
    - #### matrix_sort
    ```C
    void matrix_sort(const MATRIX_t* src, MATRIX_t** dest);

    Sort the matrix.
    ```