
source("Gaussian Elimination.R")

### Examples ###
A = matrix(c(0,4,6,8,10,4,5,2,10,12,17,18), 
           byrow = TRUE, nrow = 3, ncol = 4)

B = matrix(c(0,4,6,8,5,4,5,2,10,12,17,18), 
           byrow = TRUE, nrow = 3, ncol = 4)

C = matrix(c(10,4,6,8,10,4,5,2,0,12,17,18), 
           byrow = TRUE, nrow = 3, ncol = 4)

D = matrix(c(6,-1,2,1,5,-3,3,4,3,-2,1,14), 
           byrow = TRUE, nrow = 3, ncol = 4)

E = matrix(c(3,2,0,-1,4,1,0,-2,6,4,0,3), 
           byrow = TRUE, nrow = 3, ncol = 4)

Fm = matrix(c(2,-1,3,1,1,3,-2,1,3,-2,5,1), 
            byrow = TRUE, nrow = 3, ncol = 4)

G = matrix(c(1,2,-1,2,1,1,2,9,2,3,-3,-1), 
            byrow = TRUE, nrow = 3, ncol = 4)

H = matrix(c(3,-4,0,-26,2,3,0,28,0,0,0,0), 
           byrow = TRUE, nrow = 3, ncol = 4)

I = matrix(c(3,4,0,-1,2,5,0,-3,0,0,0,0), 
               byrow = TRUE, nrow = 3, ncol = 4)

J = matrix(c(1,-2,1,-1,-2,1,2,-5,3,-1,2,3,1,-3,8,-9), 
            byrow = TRUE, nrow = 4, ncol = 4)

K = matrix(c(3,4,-5,6,39,6,5,-6,5,43,9,-4,2,3,6,0,2,-3,1,13), 
           byrow = TRUE, nrow = 4, ncol = 5)

L = matrix(c( c(81,27,9,3,1,2),c(1,1,1,1,1,0),c(-1,-1,-1,-1,1,-1),
              c(108,27,6,1,0,0),c(-4,3,-2,1,0,0)), 
           byrow = TRUE, nrow = 5, ncol = 6)

Ah = matrix(c(1,2,0,0,1,1,1,0,2,3,1,0), 
            byrow = TRUE, nrow = 3, ncol = 4)

Bh = matrix(c(1,-2,3,0,-1,2,-3,0,2,-4,6,0), 
            byrow = TRUE, nrow = 3, ncol = 4)

Ch = matrix(c(6,-1,1,-12,-5,6,-2,2,-8,8,3,0,2,-4,5,3,-1,0,-4,9), 
            byrow = TRUE, nrow = 4, ncol = 5)

Dh = matrix(c(2,4,2,1,2,3,-3,1,0,-6,1,1,-2,0,8,4,-2,3,0,-10), 
            byrow = TRUE, nrow = 4, ncol = 5)

Eh = matrix(c(23,7,5,9,23,7,4,3,0,13,16,19), 
           byrow = TRUE, nrow = 3, ncol = 4)

bigM = Eh = matrix(sample(1:100, size = 110, replace = TRUE), 
                   byrow = TRUE, nrow = 10, ncol = 11)

bigMax = Eh = matrix(sample(1:100, size = 10100, replace = TRUE), 
                   byrow = TRUE, nrow = 100, ncol = 101)

### Test ###
# 3x4
gaussJordanElmination(A) ## L = {-3/5+1/10x3;2-3/2x3;x3}  
gaussJordanElmination(B) # L = {0;-7;6}     
gaussJordanElmination(C) # L = {0;-7;6}      

gaussJordanElmination(Eh) # L = {0.8896321;-5.9230769;6} similar to C

gaussJordanElmination(D) # L = {2.75;-9;-12.25}
gaussJordanElmination(E) # nicht Lösbar
gaussJordanElmination(Fm) # nicht Lösbar
gaussJordanElmination(G) # L = {1;2;3}
gaussJordanElmination(Bh) # L = {2x2;-3x3;x3}
# 2x3
gaussJordanElmination(H) # L = {2;8}
gaussJordanElmination(I) # L = {1;-1}
# 4x4
gaussJordanElmination(J) # L = {2;1;-1}
# 4x5
gaussJordanElmination(K) # L = {1;2;-2;3}
gaussJordanElmination(Ch) # L = {13;6;-5;6} 
gaussJordanElmination(Dh) # L = {0.5;1.5;-3;1}
# 5x6
gaussJordanElmination(L) # L = {-0.04861111;0.11805556;0.32638889;0.10416667;-0.5}
gaussJordanElmination(Ah) # L = {-2x3;x2;x3}

# 10x11
gaussJordanElmination(bigM) 
# L = { -0.96339722; 0.52102096; -0.39412268; -0.14762585; 0.22688343; 
# 0.31465101; -0.16597164; 0.44578498; -0.02456013; 1.38072059 }

gaussJordanElmination(bigMax) 

