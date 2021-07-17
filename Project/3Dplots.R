# Example file to show 2D functions in 3D plots
# Assume here that dim = 2 (dim is the number of variables)
# R. Le Riche
rm(list=ls())

library("rgl") # library for plots
source("./test_functions.R")

# function dimension
dim <- 2
# choose a function from test_functions.R file
fun <- quadratic
# upper and lower bounds
LB<--5
UB<-5

# start drawing the function (necessarily dim=2)
no.grid <- 300
x <- seq(LB, UB, ,no.grid)
x.grid <- expand.grid(x, x)
z <- apply(x.grid, 1, fun)
z.grid <- matrix(z, no.grid)

#### 3D rgl plot
open3d()
surface3d(x, x, z.grid, col= "lightblue")
title3d("quadratic", col="blue", font=4)
decorate3d()
aspect3d(1, 1, 1)
rgl.snapshot("./fileofplot.png", fmt="png", top=T)
# 

### 2D contour plot
# png(filename="./contour.png") # save the contour in the current directory
 contour(x, x, z.grid, nlevels=20, xlab="x", ylab="y")
# dev.off()

