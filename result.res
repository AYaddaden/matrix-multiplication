General infos:
Block size: 32
Matrix size: 1024 x 1024
_______________________________________________________
Multiply without optimization
Operation took: 8.309747 s
Multiply with temporary var optimization
Operation took: 5.413826 s
The multiplication gave the same result for the two matrix
Multiply with OpenMP optimization
Operation took: 5.593702 s
The multiplication gave the same result for the two matrix
Multiply with matrix transpose optimization
Operation took: 2.852428 s
The multiplication gave the same result for the two matrix
Multiply with matrix transpose and omp optimization
Operation took: 2.988738 s
The multiplication gave the same result for the two matrix
Multiply with 2D tiling optimization
Operation took: 4.061197 s
The multiplication gave the same result for the two matrix
Multiply with 3D tiling optimization
Operation took: 4.765365 s
The multiplication gave the same result for the two matrix
Multiply with 2D tiling and matrix transposition optimization
Operation took: 2.749773 s
The multiplication gave the same result for the two matrix
Multiply with 3D tiling and matrix transposition optimization
Operation took: 2.743929 s
The multiplication gave the same result for the two matrix
Multiply with 2D tiling and matrix transposition optimization with OpenMP
Operation took: 2.747584 s
The multiplication gave the same result for the two matrix
-----------------------------------------------------------------------
|                                   Method | Taken time | Improvement |
-----------------------------------------------------------------------
|                            Non-optimized |   8.309747 |       -3.87 |
-----------------------------------------------------------------------
|                         tmp optimization |   5.413826 |       32.33 |
-----------------------------------------------------------------------
|                      OpenMP optimization |   5.593702 |       30.08 |
-----------------------------------------------------------------------
|                   transpose optimization |   2.852428 |       64.34 |
-----------------------------------------------------------------------
|            transpose_OpenMP optimization |   2.988738 |       62.64 |
-----------------------------------------------------------------------
|                   2D_tiling optimization |   4.061197 |       49.24 |
-----------------------------------------------------------------------
|                   3D_tiling optimization |   4.765365 |       40.43 |
-----------------------------------------------------------------------
|         2D_tiling_transpose optimization |   2.749773 |       65.63 |
-----------------------------------------------------------------------
|         3D_tiling_transpose optimization |   2.743929 |       65.70 |
-----------------------------------------------------------------------
|  2D_tiling_transpose_OpenMP optimization |   2.747584 |       65.66 |
-----------------------------------------------------------------------
