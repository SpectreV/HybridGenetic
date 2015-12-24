library(plot3D)

placement <- read.csv(file="placement.hybrid.csv", head=FALSE,sep=" ")
index = placement[,1]
x = placement[,2]
y = placement[,3]
w = placement[,4]
l = placement[,5]
x1 = x + w
y1 = y + l
d = placement[,6]
rect2D(x0=x, y0=y, x1=x1, y1=y1,colvar = 1:10, alpha = 0.3, lwd = 2, main="rect2D")