### Gaussian Elimination
# Argument: matrix
# return: matrix in reduced row echelon form
gaussJordanElmination = function(matrixA){
  
  # Stops if Argument is not a matrix
  stopifnot(is.matrix(matrixA))
  
  # transpose matrix 
  matrixT = t(matrixA)
  
  # first coefficient != 0
  if(matrixT[1,1] == 0){
    # swap with one colum below
    for(i in 2:ncol(matrixT)){
      if(matrixT[1,i] != 0){
        matrixT = swapCol(matrixT,1,i)
        break
      }
    }
  }
  
  # zeros in row 1
  for(i in 2:ncol(matrixT)){
    if(matrixT[1,i] != 0){
      matrixT = calcZero(matrixT, 1, i, comp = 1)
    }
  }
  # check placement of zeros and swap columns if necessary
  if(matrixT[1,2] == 0 & matrixT[2,2] == 0){
    matrixT = swapCol(matrixT,2,3)
  }
  
  # outer loop: number of rows where zeros are created
  for(j in 2:(nrow(matrixT)-1)){
    for(i in ncol(matrixT):1){  
      
      zerosPrev = length(which(matrixT == 0))
      
      if(matrixT[j,i] != 0 & matrixT[j,j] != 0){  
        if(i != j){
          
          matrixT = calcZero(matrixT, j, i, comp = j) 
         
          # spaltenvertauschung wenn mehrere nullen entstanden sind
          zeros = length(which(matrixT == 0))
          if( (zerosPrev + 1) != zeros & i < ncol(matrixT)){ # gelöscht:matrixT[j,j] == 0 | matrixT[j-1,i] == 0) &
            
            for(k in ncol(matrixT):i){
              if( i > 1 & (length(which(matrixT[,k] == 0)) < length(which(matrixT[,i] == 0))) ){
                matrixT = swapCol(matrixT,i,k)
              }
            }
          }
        }
      }
    }
  }
  
  # row echelon form
  # print(matrixT)
  
  # change to reduced row echelon form
  for(i in 1:ncol(matrixT)){
    if(matrixT[i,i] != 0){
      matrixT[,i] = round(matrixT[,i] * (1/matrixT[i,i]), digits = 8)
    }
  }
  
  return(t(matrixT))
}

### Swaps colA with colB ###
# Argument: transposed matrix
# return: matrix with swapped columns
swapCol = function(matrixT,colA, colB){
  temp = matrixT[,colB]
  matrixT[,colB] = matrixT[,colA]
  matrixT[,colA] = temp
  return(matrixT)
}

#### Calculates zeros ###
# Arguments: transposed matrix; c = colum; r = row, comp = coefficent to compare 
# to
# return: matrix with min. one more zero
calcZero = function(matrixT,r,c,comp){
  #check if zero can be produced without changing col1
  if(matrixT[comp,comp] == matrixT[r,c]){
    # subtract c1 from c2                     
    matrixT[,c] = matrixT[,c] - matrixT[,r]
  }else{
    # change c1: divide by itself and multiply by head col3
    temp = matrixT[,r] / matrixT[comp,comp]
    temp = temp * matrixT[r,c]
    # subtract c1 from c2
    matrixT[,c] = matrixT[,c] - temp
  }
  
  return(matrixT)
}
