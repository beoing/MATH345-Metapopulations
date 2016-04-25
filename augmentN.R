### Puppy Squadron Effects

# A function to augment the population vector, based
# on the puppy squadron effects
# Input: n = current pop. vector (assume it's not a matrix...can be changed)
#        op = operation number (6->1, 7->2, 8->3, 9->4)
#        ju = juvenile reduction number (default=30)
#        ad = adult reduction number (default=20)
# Output: the properly augmented pop. vector
augmentN <- function(n, op, ju=30, ad=20){
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
