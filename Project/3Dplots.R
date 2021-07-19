# Example file to show 2D functions in 3D plots
# Assume here that dim = 2 (dim is the number of variables)
# R. Le Riche
rm(list=ls()) # clear environment

library("rgl") # library for plots
source("./test_functions.R")

# function dimension
dim <- 2
# choose a function from test_functions.R file
fun <- quadratic
# upper and lower bounds
LB<-c(-5,5) 
UB<-c(5,5)

# start drawing the function (necessarily dim=2)
no.grid <- 300
x1 <- seq(LB[1], UB[1], length.out=no.grid)
x2 <- seq(LB[2], UB[2], length.out=no.grid)
x.grid <- expand.grid(x1, x2)
z <- apply(x.grid, 1, fun)
z.grid <- matrix(z, no.grid)

### 2D contour plot
# png(filename="./contour.png") # save the contour in the current directory
contour(x1, x2, z.grid, nlevels=20, xlab="x1", ylab="x2")
# dev.off()

#### 3D rgl plot
open3d()
surface3d(x1, x2, z.grid, col= "lightblue")
title3d("quadratic", col="blue", font=4)
decorate3d()
aspect3d(1, 1, 1)
rgl.snapshot("./fileofplot.png", fmt="png", top=T)
# 



