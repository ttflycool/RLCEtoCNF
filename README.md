# RLCETOSAT
# fieldmath provides a method for mathematical calculations in the polynomial domain, 
  and we use it to perform Gaussian elimination on the matrix.

# reedsolo provides a method for constructing ReedSolomon codes, 
  we use it to obtain primitive polynomials and generator polynomials
  
# The public key of the RLCE scheme is calculated and generated in srlce
  You can choose to input parameters n, k, t, m
  
# A random error vector is generated in randomerror

# In H_G, the check matrix is obtained by generating the matrix.

# Tocnf outputs the processed cnf file

RUN 
1. In srlce input n,k,t,m
2. run  Tocnf and you will get a "out.cnf" file

