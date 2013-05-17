div_data <- function(data){
  dead = which(data$V45 == 1)
  health = which(data$V45 == 0)
  
  
  deadlen = length(dead)
  healthlen = length(health)
  
  d = sample(1:deadlen,ceiling(deadlen/2))
  h = sample(1:healthlen,ceiling(healthlen/2))
  
  c(health[h],dead[d])
}