n <- 15000
p <- 10

X <- matrix(rnorm(n * p), n, p) # no intercept!
y <- matrix(rnorm(n * 100), n, 100) # no intercept!

system.time(lm.fit(x = X, y = y[, 1]))
system.time(replicate(100, lm.fit(x = X, y = y[, 1])))

y <- rep(Y, each=length(beta)) - D%*%beta
x1 <- W
x2 <- cbind(G,W)

p1 <- ncol(x1)
p2 <- ncol(x2)
n <- nrow(y)

rss1 <- colSums(lm.fit(x=x1,y=y)$residuals^2)
rss2 <- colSums(lm.fit(x=x2,y=y)$residuals^2)

Fstat <- ((rss1 - rss2)/(p2 - p1))/(rss2/(n - p2))
