# Julio Rincones-Gamboa
# jgamboa@bio.tamu.edu
# Calculating the effective population size for the first part of the OSR/SA model

var.ne <- data.frame(matrix(0, 8, 9))
osr <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05)
colnames(var.ne) <- c("N", "common.sex", osr)
commsex <- c(1000, 500, 100, 50,1000, 500, 100, 50)
var.ne$N <- commsex
var.ne$common.sex <- rep(c("male", "female"), each=4)

GetVarNe <- function(female, male){
  x <- (4 * female * male) / (female + male)
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

var.ne.autosomes <- var.ne[-(1:4),]
write.csv(var.ne.autosomes , file ="../../figures/var.ne.autosomes.csv")
