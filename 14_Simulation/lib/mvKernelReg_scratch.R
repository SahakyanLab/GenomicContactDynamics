https://bookdown.org/egarpor/NP-UC3M/kre-ii-multmix.html

library(np)
library(rgl)

n <- 300
set.seed(123456)
X <- mvtnorm::rmvnorm(n = n, mean = c(0, 0),
                      sigma = matrix(c(2, 0.5, 0.5, 1.5), nrow = 2, ncol = 2))

m <- function(x) 0.5 * (x[, 1]^2 + x[, 2]^2)
epsilon <- rnorm(n)
Y <- m(x = X) + epsilon

# Plot sample and regression function
rgl::plot3d(x = X[, 1], y = X[, 2], z = Y, xlim = c(-3, 3), ylim = c(-3, 3),
            zlim = c(-4, 10), xlab = "X1", ylab = "X2", zlab = "Y")
lx <- ly <- 50
x_grid <- seq(-3, 3, l = lx)
y_grid <- seq(-3, 3, l = ly)
xy_grid <- as.matrix(expand.grid(x_grid, y_grid))
rgl::surface3d(x = x_grid, y = y_grid,
               z = matrix(m(xy_grid), nrow = lx, ncol = ly),
               col = "lightblue", alpha = 1, lit = FALSE)

# Local constant fit

# An alternative for calling np::npregbw without formula
bw0 <- np::npregbw(xdat = X, ydat = Y, regtype = "lc")
kre0 <- np::npreg(bws = bw0, exdat = xy_grid) # Evaluation grid is now a matrix
plot(kre0, plot.errors.method="asymptotic", plot.errors.style="band")

fig <- plot_ly(x = x_grid, y = y_grid, z = matrix(kre0$mean, nrow = lx, ncol = ly))
fig <- add_surface(p=fig)
fig <- add_surface(p=fig, z = matrix(kre0$mean, nrow = lx, ncol = ly)+100, opacity = 0.98)
rgl::surface3d(x = x_grid, y = y_grid,
               z = matrix(kre0$mean, nrow = lx, ncol = ly),
               col = "red", alpha = 0.25, lit = FALSE)
