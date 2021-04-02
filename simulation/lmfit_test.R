n <- 15000
p <- 10

X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- matrix(rnorm(n * 100), n, 100) # no intercept!

system.time(lm.fit(x = X, y = y[, 1]))
system.time(replicate(100, lm.fit(x = X, y = y[, 1])))
