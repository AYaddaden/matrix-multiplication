# matrix-multiplication

Implementation of matrix multiplication algorithme combining different techniques:
* Naive implementation
* Improve with a temporary variable
* OpenMP parallelization
* Loop tiling
* Matrix transposition

Those techniques are also combined together to measure computation time.

Here's a preview of an execution on a Intel&#174; i7 7500U CPU 

                                   Method | Taken time | Improvement 
---------------------------------------------------------------------
                            Non-optimized |   8.309747 |      0.00%
                         tmp optimization |   5.413826 |     32.33% 
                      OpenMP optimization |   5.593702 |     30.08% 
                   transpose optimization |   2.852428 |     64.34% 
            transpose_OpenMP optimization |   2.988738 |     62.64% 
                   2D_tiling optimization |   4.061197 |     49.24% 
                   3D_tiling optimization |   4.765365 |     40.43% 
         2D_tiling_transpose optimization |   2.749773 |     65.63% 
         3D_tiling_transpose optimization |   2.743929 |     65.70% 
  2D_tiling_transpose_OpenMP optimization |   2.747584 |     65.66% 
