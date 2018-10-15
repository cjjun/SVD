This is a ugly implementation of SVD with O(n^3)

The basic ideas for diagonization of symmetric matrix are:

  (1) Use the iteration x_{n+1} = ||Ax_{n}|| , x = x / ||x|| to find x with respect to eigenvalue with largest norm. THIS WORKABLE AND STABLE ONLY FOR SYMMETRIC MATRIX
  
  (2) Every iteration amplifies the components in different eigenvalues subspaces with amplitude of norm of respecting eigenvalues. Using Gram-smith normalization to eliminate projection on previous eigenvector space. This has to be done every iteration since even projection on previous eigenvectors with 1e-16 precision error can be amplified only after several iterations and  takE dominance.
  
  (3) For zeros vectors, we use B = A + I to iterate. The previous eigenvalue - eigenvetor pairs are not changed with just increase 1 in each eigenvalues, while for remaining zero egienvalues would be changed to 1. We only do so to obtain orthogonal basis for null space.
  
 The process for calculating SVD decomposition:
 
  (1) MV = A^T * A , MV = A * A^T. Seperately diagonize A MU and MV to obtain U,S,V.
  
  (2) The U and V might not be paired. That's to say qi * ui * vi^T might differ from the final decompoition with just a minus signal. So we use AV = US, AV and S to slightly adjust U. This process continues untill qi is no longer greater than 0. Since it makes no difference for qi = 0.
