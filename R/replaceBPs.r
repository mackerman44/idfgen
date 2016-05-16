
replaceBPs.Population <- function(x, bps, reps) {
  
   for (i in 1:length(bps)) {
        x$Scores[x$Scores==bps[i]] <- reps[i]
   }
   x
}