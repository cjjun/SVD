# SVD implementation
## Summary
This is a naive implementation of SVD with $O(n^3)$ following two steps:
-  Diagonizing symmetric matrix by $A^TA$ or $AA^T$
- Calculating $U,S,V$

## Diagonize symmetric matrix
The basic ideas for diagonization of symmetric matrix are:

- Use the iteration $x_{n+1} = ||Ax_{n}||$ , $x = x / ||x||$ to find x with respect to eigenvalue with largest absolute value. **THIS  WORKABLE AND STABLE ONLY FOR SYMMETRIC MATRIX**

- Every iteration amplifies the components in different eigenvalues subspaces with amplitude of norm of respecting eigenvalues. Using $Gram-smith$ normalization to eliminate projection on known eigenvector space. This has to be done onece every several iterations since even projection on previous eigenvectors with $10^{-16}$ precision error can be amplified only after several iterations and  take dominance.

- For zeros vectors, we use $B = A + I$ instead to iterate. The previous eigenvalue - eigenvetor pairs are not changed with just increase 1 in each eigenvalues, while for remaining zero egienvalues would be changed to 1. We only do so to obtain orthogonal basis for null space.
  
## SVD decomposition  
Given target matrix $[A]_{n\times m} = USV^T$, suppose $m\ge n$, process for calculating SVD decomposition:

   (1) $U\Sigma U^T = AA^T$, diagonize $AA^T$ to get $U$ and $\Sigma = SS^T$. Extract singular value from $\Sigma$ by square root and yield $S$.

   (2) $SV^T = U^TA$, where $S$, $U$, $A$ are known. Solve the linear matrix equation and yield part of $V^T$. However, we only get part of vectors and still need to find the complimentary vectors that orthogonal to known vectors. This can be found with random vector generation and Gram-Smith.

## Complexity
   Step (1) takes $O(m^2n)$ to take matrix multiplication and $O(m^3)$ for diagonal decomposition. 

   Step (2) takes $O(n^3)$. In total it takes $O(n^3)$.