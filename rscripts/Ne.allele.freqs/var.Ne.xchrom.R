# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# Calculating the variance effective population size for 
# haplodiploid model using the effective population size
# formula for X-linked loci. Use for X-chromosomes too.

var.ne <- data.frame(matrix(0,8, 9))
osr <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)
colnames(var.ne) <- c("N", "common.sex", osr)
commsex <- c(1000, 500, 100, 50,1000, 500, 100, 50)
var.ne$N <- commsex
var.ne$common.sex <- rep(c("male", "female"), each=4)

GetVarNe <- function(female, male){
  x <- (9 * female * male) / (2 * female + 4 * male)
  return(x)
}

for(i in 1:nrow(var.ne)){ 
  for(j in 3:ncol(var.ne)){
    if(var.ne$common.sex[i] == "male"){
      females <- var.ne$N[i] * as.numeric(colnames(var.ne)[j])
      males <- var.ne$N[i]
    }
    if(var.ne$common.sex[i] == "female"){
      males <- var.ne$N[i] * as.numeric(colnames(var.ne)[j])
      females <- var.ne$N[i]
    }
    var.ne[i,j] <- round(GetVarNe(female = females, male = males))
  }
}

var.ne.xchrom <- var.ne[-(1:4),]
write.csv(var.ne, file= "var.ne.xchrom")


                     
