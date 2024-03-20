# Sharp and Blunt matrices are consistent

    Code
      waldo::compare(first_100_elements_sharp, sharp_sparse_matrix)
    Output
          is(old)         | is(new)              
      [1] "dgCMatrix"     - "dsCMatrix"       [1]
      [2] "CsparseMatrix" | "CsparseMatrix"   [2]
      [3] "dsparseMatrix" | "dsparseMatrix"   [3]
      [4] "generalMatrix" - "symmetricMatrix" [4]
      [5] "AnyMatrix"     -                      
      [6] "V3Matrix"      -                      
      [7] "dMatrix"       | "dMatrix"         [5]
      [8] "sparseMatrix"  | "sparseMatrix"    [6]
      [9] "compMatrix"    | "compMatrix"      [7]
      
      `old@i[1:45]`: 0 3 5 7 9 1 2 9 1 2 and 35 more...
      `new@i[1:24]`: 0         1     1 2            ...
      
      `old@p`: 0 5 8 13 21 26 32 37 42 44 and 1 more...
      `new@p`: 0 1 2  4  7  9 13 17 21 23           ...
      
      `old@Dimnames[[1]]` is a character vector ('AL669831.5', 'NOC2L', 'HES4', 'ISG15', 'SDF4', ...)
      `new@Dimnames[[1]]` is NULL
      
      `old@x[1:29]`: 5 1 1 1 1 5 1 2 1 5 and 19 more...
      `new@x[1:12]`: 5         5 1     5            ...
      
      `old@x[31:45]`: 1 1 1 1 1 5 1 1 1 1 and 5 more...
      `new@x[14:24]`: 1 1 1     5 1 1 1             ...
      
      `old@uplo` is absent
      `new@uplo` is a character vector ('U')

