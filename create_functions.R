### default matrix model (no conservation efforts)

# function to create M matrix
createM <- function(){
  ### dispersal matrix elements
  ## early stage movement 
  # between patches
  l21 <- l12 <- .12
  l31 <- l13 <- .03
  l42 <- l24 <- .1
  l53 <- l35 <- .06
  l54 <- l45 <- .1
  # stay in patch
  l11 <- .85
  l22 <- .78
  l33 <- .91
  l44 <- .8
  l55 <- .84
  # no movement
  l41 <- l14 <- 0
  l51 <- l15 <- 0
  l52 <- l25 <- 0
  l32 <- l23 <- 0
  l34 <- l43 <- 0

  L <- matrix(
        c(l11,l12,l13,l14,l15,
          l21,l22,l23,l24,l25,
          l31,l32,l33,l34,l35,
          l41,l42,l43,l44,l45,
          l51,l52,l53,l54,l55),
        nrow=5,byrow=TRUE)

  ## late stage movement 
  # between patches
  v21 <- v12 <- .25
  v31 <- v13 <- .05
  v42 <- v24 <- .25
  v53 <- v35 <- .1
  v54 <- v45 <- .25
  # stay in patch
  v11 <- .7
  v22 <- .5
  v33 <- .85
  v44 <- .5
  v55 <- .65
  # no movement
  v41 <- v14 <- 0
  v51 <- v15 <- 0
  v52 <- v25 <- 0
  v32 <- v23 <- 0
  v34 <- v43 <- 0

  V <- matrix(
    c(v11,v12,v13,v14,v15,
      v21,v22,v23,v24,v25,
      v31,v32,v33,v34,v35,
      v41,v42,v43,v44,v45,
      v51,v52,v53,v54,v55),
    nrow=5,byrow=TRUE)

  ## combined transition matrix
  Z = matrix(rep(0,25), nrow=5, byrow=TRUE)
  Mp1 = cbind(L,Z)
  Mp2 = cbind(Z,V)
  M = rbind(Mp1,Mp2)
  M
}

# function to create B matrix
createB <- function(){
  ## patch 1 demography
  a11_1 <- 2/3
  a12_1 <- 2
  a21_1 <- 1/6
  a22_1 <- 4/5
  ## patch 2 demography
  a11_2 <- 1/3
  a12_2 <- 1
  a21_2 <- 1/6
  a22_2 <- 1/2
  ## patch 3 demography
  a11_3 <- 2/3
  a12_3 <- 3
  a21_3 <- 1/6
  a22_3 <- 79/80
  ## patch 4 demography
  a11_4 <- 1/3
  a12_4 <- 1
  a21_4 <- 1/6
  a22_4 <- 1/2
  ## patch 5 demography
  a11_5 <- 2/3
  a12_5 <- 3/2
  a21_5 <- 1/6
  a22_5 <- 4/5

  ## B11 matrix
  B11 = matrix(
    c(a11_1,0,0,0,0,
      0,a11_2,0,0,0,
      0,0,a11_3,0,0,
      0,0,0,a11_4,0,
      0,0,0,0,a11_5),
    nrow=5,
    byrow=TRUE)

  ## B12 matrix
  B12 = matrix(
    c(a12_1,0,0,0,0,
      0,a12_2,0,0,0,
      0,0,a12_3,0,0,
      0,0,0,a12_4,0,
      0,0,0,0,a12_5),
    nrow=5,
    byrow=TRUE)

  ## B21 matrix
  B21 = matrix(
    c(a21_1,0,0,0,0,
      0,a21_2,0,0,0,
      0,0,a21_3,0,0,
      0,0,0,a21_4,0,
      0,0,0,0,a21_5),
    nrow=5,
    byrow=TRUE)

  ## B22 matrix
  B22 = matrix(
    c(a22_1,0,0,0,0,
      0,a22_2,0,0,0,
      0,0,a22_3,0,0,
      0,0,0,a22_4,0,
      0,0,0,0,a22_5),
    nrow=5,
    byrow=TRUE)

  # create B matrix
  Bp1 = cbind(B11, B12)
  Bp2 = cbind(B21,B22)
  B = rbind(Bp1,Bp2)
  B
}

# function to calculate heuristic
heuristic <- function(n, operation){
  
  # weights
  w <- c(1, 1, 1, 1, 2)
  # juvenile weight
  j <- .85
  # EDM cost
  edm <- 0.5
  # Puppies cost
  puppies <- 50

  value <- (n[1] * w[1] + n[2] * w[2] + n[3] * w[3] + 
    n[4] * w[4] + n[5] * w[5]) * j

  value <- value + (n[6] * w[1] + n[7] * w[2] + n[8] * 
    w[3] + n[9] * w[4] + n[10] * w[5])

  if(operation > 0 && operation < 6) {
    value <- value + edm
  } else if(operation >= 6) {
    value <- value + puppies
  } # else nop

  value
}
