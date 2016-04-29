
# A function to create matrix M.
# Default matrix model (no conservation efforts)
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

### Puppy Squadron Effects
# A function to augment the population vector, based
# on the puppy squadron effects
# Input: n = current pop. vector (assume it's not a matrix...can be changed)
#        op = operation number (6->1, 7->2, 8->3, 9->4)
#        ju = juvenile reduction number (default=30)
#        ad = adult reduction number (default=15)
# Output: the properly augmented pop. vector
augmentN <- function(n, op, ju=30, ad=15){
  if (op==6) {
    n <- n - c(ju,0,0,0,0,ad,0,0,0,0)
  } else if (op==7) {
    n <- n - c(0,ju,0,0,0,0,ad,0,0,0)
  } else if (op==8) {
    n <- n - c(0,0,ju,0,0,0,0,ad,0,0)
  } else {
    n <- n - c(0,0,0,ju,0,0,0,0,ad,0)
  }
  n[n<0]=0
  n
}

# A function to calculate heuristic of a node
# Input: n = current node
#        operation = conservation action taken
# Output: numerical heuristic value of a node
heuristic <- function(n, operation){
  
  # weights
  w <- c(1, 1, 1, 1, 5)
  # juvenile weight
  j <- .75
  # EDM cost
  edm <- 5
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

# Prints out ASCII art of a model at time t
# Input: sm = solution_model
#         t = time segment of model to print out
# Output: screen output of model art
print_model <- function(sm = solution_model, t = 15){
  
  cat("The model over time is: \n\n")
  
  colnames(sm) <- c("n0", "n1", "n2", "n3", "n4", "n5", "n6", 
                  "n7", "n8", "n9", "n10", "n11", "n12", 
                  "n13", "n14", "n15") 
  rownames(sm) <- c("j1", "j2", "j3", "j4", "j5", "a1", "a2", "a3", "a4", "a5")
  print(sm)
  Sys.sleep(.5)
  cat("\n\n")
  
  
  
  
  cat(sprintf("The model at time t = %d is:  \n\n", t))
  cat(sprintf("     --------         --------\n"))
  cat(sprintf("     |%06s|         |%06s|\n", strtrim(sm[4, t+1],6), strtrim(sm[5, t+1],6)))
  cat(sprintf("     |------|=========|------|\n"))
  cat(sprintf("     |%06s|         |%06s|\n", strtrim(sm[9, t+1],6), strtrim(sm[10, t+1],6)))
  cat(sprintf("     --------         --------\n"))
  cat(sprintf("       //                ||   \n"))
  cat(sprintf("      //          ^  ^   ||   \n"))
  cat(sprintf("     //         ^   ^    ||   \n"))
  cat(sprintf("    //        ^   ^  ^   ||   \n"))
  cat(sprintf("   //          ^  ^   ^  ||   \n"))
  cat(sprintf("--------     ^    ^  ^   ||   \n"))
  cat(sprintf("|%06s|       ^  ^   ^  ||   \n", strtrim(sm[2, t+1],6)))
  cat(sprintf("|------|    ^  ^   ^     ||   \n"))
  cat(sprintf("|%06s|     ^   ^  ^    ||   \n", strtrim(sm[7, t+1],6)))
  cat(sprintf("--------   ^   ^   ^     ||   \n"))
  cat(sprintf("   ||         ^  ^       ||   \n"))
  cat(sprintf("   ||       ^   ^        ||   \n"))
  cat(sprintf("   ||         ^          ||   \n"))
  cat(sprintf("   ||      ^   ^         ||   \n"))
  cat(sprintf("   ||         ^          ||   \n"))
  cat(sprintf("--------              --------\n"))
  cat(sprintf("|%06s|              |%06s|\n", strtrim(sm[1, t+1],6), strtrim(sm[3, t+1],6)))
  cat(sprintf("|------|==============|------|\n"))
  cat(sprintf("|%06s|              |%06s|\n", strtrim(sm[6, t+1],6), strtrim(sm[7, t+1],6)))
  cat(sprintf("--------              --------\n\n\n"))
}

# Function to print table of operations
# Input: 
print_ops <- function(ops = solution_operations){
  solutions <- data.frame(stringsAsFactors = TRUE)
  for(i in 1:(length(ops)-1))
  {
    if (ops[i+1]==0) {
      solutions[i,1] <- "No operation"
    } else if (ops[i+1]==1) {
      solutions[i,1] <- "EDM between 1 and 2"
    } else if (ops[i+1]==2) {
      solutions[i,1] <- "EDM between 1 and 3"
    } else if (ops[i+1]==3) {
      solutions[i,1] <- "EDM between 2 and 4"
    } else if (ops[i+1]==4) {
      solutions[i,1] <- "EDM between 3 and 5"
    } else if (ops[i+1]==5) {
      solutions[i,1] <- "EDM between 4 and 5"
    } else if (ops[i+1]==6) {
      solutions[i,1] <- "Puppy Squad in 1"
    } else if (ops[i+1]==7) {
      solutions[i,1] <- "Puppy Squad in 2"
    } else if (ops[i+1]==8) {
      solutions[i,1] <- "Puppy Squad in 3"
    } else {
      solutions[i,1] <- "Puppy Squad in 4"
    }
  }
  colnames(solutions) <- c("Operation")
  rows <- c(1:(length(ops)-1))
  years <- rep("Year",(length(ops)-1))
  colons <- rep(":",(length(ops)-1))
  names <- paste(years,rows,sep=" ")
  names <- paste(names,colons, sep="")
  rownames(solutions)<-names
  solutions
}

# Main method of simulation
# Input: j = initial number of juveniles
#        a = initial number of adults
#        t = time of simulation
# Output: solution model, ops, heuristics (global vars)
main <- function(j = 100, a = 80, t = 15, fun = 10){
  
  if(fun > 0) {
    cat("Certainly! I'll start that simulation up for you.\n")
    cat("\n")
    cat("\n")
    cat("\n")
    Sys.sleep(1)
    cat("BEGINNING ANTI AQUATIC WOOLLY MAMMOTH EXTERMINATION SIMULATION PROGRAM\n")
    for(lines in 0:fun)
    {
      marks = "!"
      for(numMarks in 0:floor(runif(1, 0, 8))){
        marks <- paste(marks, "!", sep = "")
      }
      Sys.sleep(runif(1, 0, 0.2))
      cat(sprintf("EXTERMINATE%s\n", marks))
    }
    Sys.sleep(.2)
    cat("EXTERMINATE!!!!!!!!!!!!!!!\n")
    cat("\n")
    cat("Your simulation has been completed! Anything else?\n")
  }
  
  
  solution_model <<- matrix(0, 10, t+1)
  solution_model[,1] <<- c(j,0,0,0,0,a,0,0,0,0) #n0
  # solution_model[,1] <<- c(j,10,10,0,0,a,5,5,0,0) #n0 
  solution_operations <<- rep(0, t+1)
  solution_heuristics <<- rep(0, t+1)
  solution_heuristics[1] <<- heuristic(solution_model[,1], 0)
  
  B <- createB()
  
  for(year in 1:t) {
  
    # default operation: do nothing
    current_node <- B%*%createM()%*%solution_model[,year]
    current_operation <- 0
    current_heuristic <- heuristic(current_node, current_operation)
  
    for(operation in 1:9) {
      M <- createM()
      n <- solution_model[,year]
  
      if (operation <= 5) {
        M <- augmentM(M, operation)
      } else {
        n <- augmentN(n, operation)
      }
  
      node <- B%*%M%*%n
      heur <- heuristic(node, operation)
  
      if(heur < current_heuristic) {
        current_node <- node
        current_operation <- operation
        current_heuristic <- heur
      }
  
    }
  
    solution_model[,year+1] <<- current_node
    solution_operations[year+1] <<- current_operation
    solution_heuristics[year+1] <<- current_heuristic
  }
  
  
  # print out solutions
  # print(solution_operations)
  # print_model(solution_model, t)

    
}

mammoths <- function(){
  cat("Sure thing, Chris!\n\n")
  Sys.sleep(.5)
  cat("Let's use a starting population of 100 juveniles and 80 adults.\n\n")
  Sys.sleep(.5)
  cat("This is what happens to the AWH population after 15 years, with no intervention.\n\n")
  n0 = c(100,0,0,0,0,80,0,0,0,0)
  M <- createM()
  B <- createB()
  # number of time units - 1
  t = 16
  n=matrix(0,10,t)
  n[,1]= n0
  for (i in 2:t)
  {
    n[,i] = B%*%M%*%n[,i-1]
  }
  
  print_model(n)
  
  Sys.sleep(2)
  cat("\n\nAnything else?\n")
}

print_sim <- function()
{
  print_model()
  print_ops()
}

final <- function()
{
  main(150, 100, 15, 0)
  print_model()
  print_ops()
}

startup <- function(){
  
  cat("\014") 
  cat("Hello, Chris!  What can I do for you today?\n") 

}

cleanup <- function() {
  ENV <- globalenv()
  ll <- ls(envir = ENV)
  ll <- ll[ll != "cleanup"]
  rm(list = ll, envir = ENV)
  cat("\014") 
}

startup()
