# Smith-normal-form-of-q-Varchenko-matrices
This repository contains code relevant to calculating the Smith normal form of the $q$-Varchenko matrices corresponding to the hyperplane arrangements that form $n$-cubes and $n$-simplices in $\mathbb{R}^n$. 

populate_varchenko_hypercube.m and populate_varchenko_simplex.m generate the $q$-Varchenko matrices corresponding to the hyperplane arrangements that form $n$-cubes and $n$-simplices in $\mathbb{R}^n$. Each function takes a single parameter, which is the desired dimension. 

smith_normalize_hypercube.m calculates the Smith normal form of the $q$-Varchenko matrix corresponding to the hyperplane arrangement forming an $n$-cube in $\mathbb{R}^n$ by algorithmically generating all of the required transition and permutation matrices. This function calls populate_varchenko_hypercube in the process, so it only needs a single parameter: the desired dimension. However, giving the function any secondary parameter will cause it to print out the transition and permutation matrices it uses along the way. 

smith_normalize_simplex.m calculates the Smith normal form of the $q$-Varchenko matrix corresponding to the hyperplane arrangement forming an $n$-simplex in $\mathbb{R}^n$ by algorithmically generating all of the required transition and permutation matrices. This function calls populate_varchenko_simplex in the process, so it only needs a single parameter: the desired dimension. However, giving the function any secondary parameter will cause it to print out the transition and permutation matrices it uses along the way. 
