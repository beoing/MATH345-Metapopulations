
solution_model <- matrix(0, 10, t+1)
solution_model[,1] = c(60,0,0,0,0,50,0,0,0,0) #n0
solution_operations <- rep(0, t+1)
solution_heuristics <- rep(0, t+1)
solution_heuristics[1] = heuristic(solution_model[,1], 0)

B <- createB()
t <- 15

for(year in 1:t) {

  # default operation: do nothing
  current_node <- B%*%createM()%*%solution_model[,year]
  current_operation <- 0
  current_heuristic <- heuristic(current_node, current_solution)

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

  solution_model[,year+1] <- current_node
  solution_operations[year+1] <- current_operation
  solution_heuristics[year+1] <- current_heuristic
}

# print out solutions

