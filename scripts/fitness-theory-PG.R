mal <- read.csv("rare.male.250.iter.csv", as.is=T)
bar <- mal[43085,2:10]
bar <- mal[98046, 2:10]
bar

s <- bar$s
h <- bar$h

# p^2 + 2pq + q^2

#Sex chromosomes

h = 1
#X is beneficial

bar$X^2 * 1 +
  2*bar$X*(1-bar$X) * (1+h*s) +
  (1-bar$X)^2 * (1 + s)

# XY
h =0

#Y is beneficial
bar$X*bar$Y                 * (1 + s) +
  bar$X* (1-bar$Y)          * (1+h*s) +
  (1-bar$X) * bar$Y         * (1+h*s) +
  (1-bar$Y) * (1-bar$X)     * (1 )

#.       male     female
# 00     1         1
# 01     1+hs      1/1+hs  
# 11     1+s       1/1+s

#.       male     female
# 00     1         1+s
# 01     1+hs      1+hs  
# 11     1+s       1

h=1

s <- seq(from=.01, to=.99, length.out=100)

wfem <- 1/(1+s)
wmal <- 1+s
plot(1/wfem~s, ylim=c(1,2))
points(wmal/1~s,col="green")



temp.cur<-dat[dat$gens > 1000,]
