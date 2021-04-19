# RLCEtoSAT
  1.fieldmath provides a method for mathematical calculations in the polynomial domain, 
    and we use it to perform Gaussian elimination on the matrix.

  2.reedsolo provides a method for constructing ReedSolomon codes, 
    we use it to obtain primitive polynomials and generator polynomials
    At the same time, various calculation methods in the polynomial domain are provided.
    The most important thing is to simplify the calculation, construct the exponential logarithm table, 
    so that the calculation in the polynomial domain can be completed by looking up the table
  
  
  3.The public key of the RLCE scheme is calculated and generated in srlce.
    You can choose to input parameters n, k, t, m
    According to these parameters will generate：
    Generating Matrix of Generalized RS Code-Gs(g,n）
    Sparse matrix A-getMatrixA(k,w)
    Non-singular matrix-Snonsingula(n,k)
    Non-standard form of public key-Gpub(g,n,k,w)
  
  4.A random error vector is generated in randomerror

  5.In H_G, the check matrix is obtained by generating the matrix.

  6.Tocnf output the cnf file

#RUN 
1. In srlce input n,k,t,m
2. run  Tocnf and you will get a "out.cnf" file

