### EDM Concert Effects
# assume that we have a matrix M that has been created

# A function to augment M, based on EDM effects
# Input: M = transition matrix
#        n = assignment number (E1->21,E2->31,E3->42,E4->53,E5->54)
#        r = reduction coefficient
#        d = killing coefficient (1-d is the rate at which the blocked
#                                   migrants survive, default is 1)
# Output: the properly augmented matrix
augmentM <- function(M,n,r=.01,d=1){
  
  if (n==1) {
    # juvenile
    ov_j <- M[2,1]
    M[1,1] <- M[1,1] + (1-r)*ov_j*d
    M[2,2] <- M[2,2] + (1-r)*ov_j*d
    M[2,1] <- r*ov_j
    M[1,2] <- r*ov_j
    # adult
    ov_a <- M[7,6]
    M[6,6] <- M[6,6] + (1-r)*ov_a*d
    M[7,7] <- M[7,7] + (1-r)*ov_a*d
    M[7,6] <- r*ov_a
    M[6,7] <- r*ov_a
  } else if (n==2) {
    # juvenile
    ov_j <- M[3,1]
    M[1,1] <- M[1,1] + (1-r)*ov_j*d
    M[3,3] <- M[3,3] + (1-r)*ov_j*d
    M[3,1] <- r*ov_j
    M[1,3] <- r*ov_j
    # adult
    ov_a <- M[8,6]
    M[6,6] <- M[6,6] + (1-r)*ov_a*d
    M[8,8] <- M[8,8] + (1-r)*ov_a*d
    M[8,6] <- r*ov_a
    M[6,8] <- r*ov_a
  } else if (n==3) {
    # juvenile
    ov_j <- M[4,2]
    M[2,2] <- M[2,2] + (1-r)*ov_j*d
    M[4,4] <- M[4,4] + (1-r)*ov_j*d
    M[4,2] <- r*ov_j
    M[2,4] <- r*ov_j
    # adult
    ov_a <- M[9,7]
    M[7,7] <- M[7,7] + (1-r)*ov_a*d
    M[9,9] <- M[9,9] + (1-r)*ov_a*d
    M[9,7] <- r*ov_a
    M[7,9] <- r*ov_a
  } else if (n==4) {
    # juvenile
    ov_j <- M[5,3]
    M[3,3] <- M[3,3] + (1-r)*ov_j*d
    M[5,5] <- M[5,5] + (1-r)*ov_j*d
    M[5,3] <- r*ov_j
    M[3,5] <- r*ov_j
    # adult
    ov_a <- M[10,8]
    M[8,8] <- M[8,8] + (1-r)*ov_a*d
    M[10,10] <- M[10,10] + (1-r)*ov_a*d
    M[10,8] <- r*ov_a
    M[8,10] <- r*ov_a
  } else {
    # juvenile
    ov_j <- M[5,4]
    M[4,4] <- M[4,4] + (1-r)*ov_j*d
    M[5,5] <- M[5,5] + (1-r)*ov_j*d
    M[5,4] <- r*ov_j
    M[4,5] <- r*ov_j
    # adult
    ov_a <- M[10,9]
    M[9,9] <- M[9,9] + (1-r)*ov_a*d
    M[10,10] <- M[10,10] + (1-r)*ov_a*d
    M[10,9] <- r*ov_a
    M[9,10] <- r*ov_a
  }
  M
}



