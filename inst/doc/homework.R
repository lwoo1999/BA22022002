## -----------------------------------------------------------------------------
m1 <- matrix(1, nr = 2, nc = 2)
m2 <- matrix(data = c(1, 2, 3, 4), nr = 2, nc = 2)

## -----------------------------------------------------------------------------
m1

## -----------------------------------------------------------------------------
m2

## -----------------------------------------------------------------------------
m1 %*% m2

## -----------------------------------------------------------------------------
knitr::kable(head(iris))

## -----------------------------------------------------------------------------
plot(iris$Sepal.Length, iris$Petal.Length, col=iris$Species)

## -----------------------------------------------------------------------------
cor.test(iris$Sepal.Length, iris$Petal.Length, method="pearson")

## -----------------------------------------------------------------------------
cor.test(iris$Sepal.Length, iris$Petal.Length, method="spearman")

## -----------------------------------------------------------------------------
cor.test(iris$Sepal.Length, iris$Petal.Length, method="kendall")

## -----------------------------------------------------------------------------
my.sample <- function(x, size, prob = NULL) {
    n <- length(x)

    if (is.null(prob)) {
        prob <- rep(1/n, n)
    }

    cp <- cumsum(prob)

    u <- runif(size)

    ret <- x[findInterval(u, cp) + 1]

    return(ret)
}

## -----------------------------------------------------------------------------
my.sample(1:3, 20)

## -----------------------------------------------------------------------------
x <- my.sample(1:5, 1000, prob=c(1,2,3,4,5) / 15)
hist(x, prob = TRUE, breaks = 0:5)

## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
x <- numeric(length(u))
x[u < 0.5] <- log(2*u[u < 0.5])
x[u > 0.5] <- -log(2*(1-u[u > 0.5]))

hist(x, prob = TRUE, main = expression(f(x)==0.5*exp(-abs(y))))
y <- seq(-5, 5, .1)
lines(y, 0.5*exp(-abs(y)))

## -----------------------------------------------------------------------------
my.beta <- function(a, b, x) {
    return(gamma(a + b) / gamma(a) / gamma(b) * x^(a-1) * (1-x)^(b - 1))
}

my.genbeta <- function(a, b, n) {
    k <- 0
    y <- numeric(n)
    max.beta <- ((a - 1) / (a + b - 2)) ^ (a - 1) * ((b - 1) / (a + b - 2)) ^ (b - 1) * gamma(a + b) / gamma(a) / gamma(b)


    while (k < n) {
        u <- runif(1)
        x <- runif(1)
        if (my.beta(a, b, x) / max.beta > u) {
            k <- k + 1
            y[k] <- x
        }
    }

    return(y)
}

## -----------------------------------------------------------------------------
x = my.genbeta(3, 2, 1000)
hist(x, prob = TRUE, main = expression(f(x)==Beta(3,2)))
y <- seq(0, 1, .01)
lines(y, my.beta(3, 2, y))

## -----------------------------------------------------------------------------
genfe <- function(n) {
    u1 <- runif(n) * 2 - 1
    u2 <- runif(n) * 2 - 1
    u3 <- runif(n) * 2 - 1

    y <- u3
    cond <- (abs(u3) >= abs(u2)) & (abs(u3) >= abs(u1))
    y[cond] <- u2[cond]

    return(y) 
}

## -----------------------------------------------------------------------------
x <- genfe(1000)
hist(x, prob = TRUE, main = "f_e")
y <- seq(-1, 1, .01)
lines(y, 0.75*(1-y^2))

## -----------------------------------------------------------------------------
buffon <- function(rho) {
    n <- 1e6
    l <- 1
    d <- l/rho

    X <- runif(n,0,d/2)
    Y <- runif(n,0,pi/2)
    pihat <- 2*l/d/mean(l/2*sin(Y)>X)
    return(pihat)
}

## -----------------------------------------------------------------------------
pihat1 <- replicate(100, buffon(1/3))
pihat2 <- replicate(100, buffon(2/3))
pihat3 <- replicate(100, buffon(1))
c(sd(pihat1), sd(pihat2), sd(pihat3))

## -----------------------------------------------------------------------------
u <- runif(1000)
u1 <- exp(u)
u2 <- exp(1-u)

cov.u1.u2 = cov(u1, u2)
var.u1.add.u2 = var(u1 + u2)

c(cov.u1.u2, var.u1.add.u2)

## -----------------------------------------------------------------------------
var.u1 = 0.5 * var.u1.add.u2 - cov.u1.u2

## -----------------------------------------------------------------------------
0.5 - 0.5 * cov.u1.u2 / var.u1

## -----------------------------------------------------------------------------
# antithetic variate approach
u <- runif(1000)
u1 <- exp(u)
u2 <- exp(1-u)

theta1 <- (u1 + u2) / 2

# simple Monte Carlo method
u <- runif(1000)
theta2 <- exp(u)

c(mean(theta1), mean(theta2))

## -----------------------------------------------------------------------------
c(var(theta1), var(theta2), 1 - var(theta1)/var(theta2))

## -----------------------------------------------------------------------------
u <- runif(10000)
x <- 1/(1-u)
theta <- x^2 /sqrt(2*pi) * exp(-x^2/2) / (1/x^2)

mean(theta)

## -----------------------------------------------------------------------------
M <- 10000; k <- 5
r <- M/k
theta <- numeric(M)
g <- function(x) exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
f <- function(x) exp(-x) / (1 - exp(-1))

for (j in 1:k) {
    u <- runif(r, (j-1)/k, j/k)
    x <- - log(1 - u * (1 - exp(-1)))
    theta[(r*(j-1)+1):(r*j)] <- g(x)/f(x)
}

mean(theta)

## -----------------------------------------------------------------------------
M <- 10000
N <- 0

n <- 20
alpha <- 0.05

for (i in 1:M) {
    sample <- rchisq(n, df=2)
    upper <- mean(sample) + sd(sample)/sqrt(n)*qt(1-alpha/2, df=n-1)
    lower <- mean(sample) - sd(sample)/sqrt(n)*qt(1-alpha/2, df=n-1)

    if (2 < upper && 2 > lower) N = N + 1
}

N/M

## -----------------------------------------------------------------------------
M <- 10000
n <- 20

p.val <- numeric(M)

for (i in 1:M) {
    sample <- rchisq(n, df=1)
    mu <- 1
    t <- (mean(sample) - mu) / (sd(sample) / sqrt(n))
    p.val[i] <- 2*(1-pt(abs(t), n-1))
}

mean(p.val<=0.05)

## -----------------------------------------------------------------------------
M <- 10000
n <- 20

p.val <- numeric(M)

for (i in 1:M) {
    sample <- runif(n, 0, 2)
    mu <- 1
    t <- (mean(sample) - mu) / (sd(sample) / sqrt(n))
    p.val[i] <- 2*(1-pt(abs(t), n-1))
}

mean(p.val<=0.05)

## -----------------------------------------------------------------------------
M <- 10000
n <- 20

p.val <- numeric(M)

for (i in 1:M) {
    sample <- rexp(n, 1)
    mu <- 1
    t <- (mean(sample) - mu) / (sd(sample) / sqrt(n))
    p.val[i] <- 2*(1-pt(abs(t), n-1))
}

mean(p.val<=0.05)

## -----------------------------------------------------------------------------
M <- 1000
m <- 1000

alpha <- 0.1

FWER.bonferroni <- numeric(M)
FDR.bonferroni <- numeric(M)
TPR.bonferroni <- numeric(M)

FWER.BH <- numeric(M)
FDR.BH <- numeric(M)
TPR.BH <- numeric(M)


for (i in 1:M) {
    p <- numeric(m)
    null.is.true <- logical(m)
    p[1:(0.95*m)] <- runif(0.95*m)
    null.is.true[1:(0.95*m)] <- TRUE
    p[(0.95*m+1):m] <- rbeta(0.05*m, 0.1, 1)

    p.bonferroni <- p.adjust(p, method = "bonferroni")
    p.BH <- p.adjust(p, method = "BH")

    FWER.bonferroni[i] <- any(p.bonferroni < alpha & null.is.true)
    FDR.bonferroni[i] <- sum(p.bonferroni < alpha & null.is.true) / sum(p.bonferroni < alpha)
    TPR.bonferroni[i] <- sum(p.bonferroni < alpha & !null.is.true) / sum(!null.is.true)

    FWER.BH[i] <- any(p.BH < alpha & null.is.true)
    FDR.BH[i] <- sum(p.BH < alpha & null.is.true) / sum(p.BH < alpha)
    TPR.BH[i] <- sum(p.BH < alpha & !null.is.true) / sum(!null.is.true)
}

tb <- matrix(c(mean(FWER.bonferroni),mean(FDR.bonferroni),mean(TPR.bonferroni), 
               mean(FWER.BH),mean(FDR.BH),mean(TPR.BH)), 
             nrow = 2, ncol = 3, byrow = TRUE,
             dimnames = list(c("bonferroni", "BH"),
                             c("FWER", "FDR", "TPR")))

knitr::kable(tb)

## -----------------------------------------------------------------------------
library(boot)
B <- 1000
m <- 1000

b.lambda <- function (x, i) 1 / mean(x[i])

tb <- matrix(numeric(12), 
             nrow = 3, ncol = 4, byrow = TRUE,
             dimnames = list(c("n=5", "n=10", "n=20"),
                             c("bias", "se", "theoretical bias", "theoretical se")))

for (n in c(5, 10, 20)) {
    row <- stringr::str_c("n=", n)
    tb[row, "theoretical bias"] <- 2 / (n-1)
    tb[row, "theoretical se"] <- 2*n / (n-1) / sqrt(n-2)

    for (i in 1:m) {
        s <- rexp(n, 2)
        result <- boot(s, b.lambda, B)
        tb[row, "bias"] = tb[row, "bias"] + mean(result$t) - result$t0
        tb[row, "se"] = tb[row, "se"] + sd(result$t)
    }

    tb[row, "bias"] = tb[row, "bias"] / m
    tb[row, "se"] = tb[row, "se"] / m
}


knitr::kable(tb)

## -----------------------------------------------------------------------------
law <- cbind(
    LSAT = c(576, 635, 558, 578, 666, 580, 555, 661, 651, 605, 653, 575, 545, 572, 594),
    GPA = c(339, 330, 281, 303, 344, 307, 300, 343, 336, 313, 312, 274, 276, 288, 296)
)

r_ <- function(x, i) cor(x[i,1], x[i,2])
r <- function(x, i) {
    b <- boot(x[i,], r_, 100)
    return(c(cor(x[i,1], x[i,2]), var(b$t)))
}

obj <- boot(law, r, 2000)

boot.ci(obj, type="stud")

## -----------------------------------------------------------------------------
library(boot)

b.mean <- function (x, i) return(mean(x[i]))
obj <- boot(aircondit$hours, b.mean, 2000)
boot.ci(obj)

## -----------------------------------------------------------------------------
library(bootstrap)

b.theta <- function(x, i) {
    lambdas <- eigen(cor(x[i,]))$values
    return(lambdas[1] / sum(lambdas))
}

n <- 88

theta.hat <- b.theta(scor, 1:n)
theta.jack <- numeric(n)

for(i in 1:n){
    theta.jack[i] <- b.theta(scor, (1:n)[-i])
}

bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-theta.hat)^2))

c(bias.jack=bias.jack, se.jack=se.jack)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)

n <- length(magnetic) # in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)) # n choose 2 error pair, 2*C_n^2


# fit models on leave-two-out samples
i <- 1 # indice of e_i
for (k in 2:n) {
    for (l in 1:(k-1)) {
        y <- magnetic[-k][-l]
        x <- chemical[-k][-l]

        J1 <- lm(y ~ x)
        yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
        e1[i] <- magnetic[k] - yhat1
        yhat1 <- J1$coef[1] + J1$coef[2] * chemical[l]
        e1[i+1] <- magnetic[l] - yhat1

        J2 <- lm(y ~ x + I(x^2))
        yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
        e2[i] <- magnetic[k] - yhat2
        yhat2 <- J2$coef[1] + J2$coef[2] * chemical[l] + J2$coef[3] * chemical[l]^2
        e2[i+1] <- magnetic[l] - yhat2

        J3 <- lm(log(y) ~ x)
        logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
        yhat3 <- exp(logyhat3)
        e3[i] <- magnetic[k] - yhat3
        logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[l]
        yhat3 <- exp(logyhat3)
        e3[i+1] <- magnetic[l] - yhat3

        J4 <- lm(log(y) ~ log(x))
        logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
        yhat4 <- exp(logyhat4)
        e4[i] <- magnetic[k] - yhat4
        logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[l])
        yhat4 <- exp(logyhat4)
        e4[i+1] <- magnetic[l] - yhat4

        i = i + 2
    }
}

c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
library(boot)

boot.cvm <- function(z, ix, sizes) {
    n <- sizes[1]; m <- sizes[2]; nm <- n + m
    z <- z[ix]

    xs <- z[1:n]; ys <- z[(n+1):nm]
    F <- ecdf(xs)
    G <- ecdf(ys)

    m*n/(m+n)^2 * ( sum((F(xs) - G(xs))^2) + sum((F(ys) - G(ys))^2) )
}

perm.cvm <- function(xs, ys) {
    boot.obj <- boot(data = c(xs, ys), statistic = boot.cvm, R = 999,
                     sim = "permutation", sizes = c(length(xs), length(ys)))

    ts <- c(boot.obj$t0,boot.obj$t)
    mean(ts>=ts[1])
}

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

perm.cvm(x, y)

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
    X <- x - mean(x)
    Y <- y - mean(y)
    outx <- sum(X > max(Y)) + sum(X < min(Y))
    outy <- sum(Y > max(X)) + sum(Y < min(X))

    return(max(c(outx, outy)))
}

boot.c5 <- function(z, ix, sizes) {
    n <- sizes[1]; m <- sizes[2]; nm <- n + m
    z = z[ix]

    count5test(z[1:n], z[(n+1):nm])
}

perm.c5 <- function(xs, ys) {
    

    boot.obj <- boot(data = c(xs, ys), statistic = boot.c5, R = 999,
                     sim = "permutation", sizes = c(length(xs), length(ys)))

    ts <- c(boot.obj$t0,boot.obj$t)
    mean(ts>=ts[1])
}

## -----------------------------------------------------------------------------
x <- rnorm(100, 0, 1)
y <- rnorm(20, 0, 1)
perm.c5(x, y)

## -----------------------------------------------------------------------------
x <- rnorm(30, 0, 1)
y <- rnorm(70, 0, 2)
perm.c5(x, y)

## -----------------------------------------------------------------------------
produce.a <- function(N, b1, b2, b3, f0) {
    X1 <- rpois(N, 1)
    X2 <- rexp(N, 1)
    X3 <- rbinom(N, 1, 0.5)

    g <- function(a) {
        tmp <- exp(-a - b1*X1 - b2*X2 - b3*X3)

        mean(1/(1+tmp)) - f0
    }

    res <- uniroot(g, c(-100, 0))
    res$root
}

## -----------------------------------------------------------------------------
f0s = c(1e-1, 1e-2, 1e-3, 1e-4)
as = sapply(f0s, function (f0) produce.a(1e6, 0, 1, -1, f0))
as

## -----------------------------------------------------------------------------
plot(-log10(f0s), as)
lines(-log10(f0s), as)

## -----------------------------------------------------------------------------
rwMC <- function (N, f, rwsd, x0) {
    # proposal distribution for random walk
    rg <- function(xt) rnorm(1, xt, rwsd)
    dg <- function(xt, x) dnorm(x, xt, rwsd)

    # chain
    xs <- numeric(N)
    xs[1] <- x0

    # acceptance rate
    ac <- numeric(N)
    ac[1] <- 1

    for (i in 1:(N-1)) {
        Y <- rg(xs[i])
        U <- runif(1)

        if (U <= f(Y)*dg(Y, xs[i])/f(xs[i])/dg(xs[i], Y)) {
            xs[i+1] = Y
            ac[i+1] = 1
        } else {
            xs[i+1] = xs[i]
        }
    }

    list(chain=xs, ac=mean(ac))
}

## -----------------------------------------------------------------------------
f <- function (x) 1/2*exp(-abs(x))
N <- 1e3

res <- rwMC(N, f, 1, 5)
plot(res$chain, type="l", main=paste("sd = 1\nacceptance rate = ", res$ac))

## -----------------------------------------------------------------------------
res <- rwMC(N, f, 0.5, 5)
plot(res$chain, type="l", main=paste("sd = 0.5\nacceptance rate = ", res$ac))

## -----------------------------------------------------------------------------
res <- rwMC(N, f, 0.1, 5)
plot(res$chain, type="l", main=paste("sd = 0.1\nacceptance rate = ", res$ac))

## -----------------------------------------------------------------------------
mc.binorm <- function(N) {
    xs <- numeric(N)
    ys <- numeric(N)

    for (i in 1:(N-1)) {
        xs[i+1] <- rnorm(1, 0.9*ys[i], 1-0.9^2)
        ys[i+1] <- rnorm(1, 0.9*xs[i+1], 1-0.9^2)
    }

    list(x = xs, y = ys)
}

res <- mc.binorm(1e4)
# discard burn-in
res$x <- res$x[1e3:1e4]
res$y <- res$y[1e3:1e4]

## -----------------------------------------------------------------------------
plot(res$x, res$y)

## -----------------------------------------------------------------------------
fit <- lm(y ~ x, res)
summary(fit)

## -----------------------------------------------------------------------------
qqnorm(fit$residuals)
qqline(fit$residuals)

## -----------------------------------------------------------------------------
plot(fitted(fit), fit$residuals)

## -----------------------------------------------------------------------------
R.hat <- function (chains) {
    Wn <- mean(apply(chains, 2, var))
    var.hat <- var(c(chains))

    var.hat / Wn
}

## -----------------------------------------------------------------------------
f <- function(x, sigma) {
    if (any(x < 0)) return (0)
    stopifnot(sigma > 0)
    return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}


m <- 1e4
sigma <- 4
chains <- matrix(0, ncol=10, nrow=m)
chains[1,] <- rchisq(10, df=1)
k <- 0
u <- matrix(runif(m*10), ncol=10, nrow=m)

for (i in 2:m) {
    for (j in 1:10) {
        xt <- chains[i-1, j]
        y <- rchisq(1, df = xt)
        num <- f(y, sigma) * dchisq(xt, df = y)
        den <- f(xt, sigma) * dchisq(y, df = xt)
        if (u[i, j] <= num/den) chains[i, j] <- y else {
            chains[i, j] <- xt
            k <- k+1 #y is rejected
        }
    }

    if (R.hat(chains[1:i,]) < 1.2) {
        chains = chains[1:i,]
        break
    }
}

## -----------------------------------------------------------------------------
n <- 2:nrow(chains)
R.hats <- lapply(n, function (n) R.hat(chains[1:n,]))

plot(n, R.hats, type="l", ylim=c(1,2))

## -----------------------------------------------------------------------------
library(coda)

chains.coda <- do.call(mcmc.list, lapply(1:10, function (i) as.mcmc(chains[,i])))

gelman.diag(chains.coda)
gelman.plot(chains.coda)

## -----------------------------------------------------------------------------
us <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
vs <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)

estimate.em <- function (us, vs, max_iter) {
    l <- 1
    n <- length(us)
    
    for (i in 1:max_iter) {
        ex <- numeric(n)
        for (j in 1:n) {
            ex[j] <- integrate(function(x) x*dexp(x, l), us[j], vs[j])$value / 
                    (pexp(vs[j], l) - pexp(us[j], l))
        }
        tmp <- n / sum(ex)

        if (abs(tmp - l) < 1e-6) {
            return(tmp)
        }

        l <- tmp
    }

    l
}

estimate.obs <- function (us, vs) {
    n <- length(us)

    logl <- function(l) {
        ret <- 0
        for (i in 1:n) {
            ret = ret + log(pexp(vs[i], l) - pexp(us[i], l))
        }
        ret
    }

    optimise(logl, c(1e-6, 1), maximum=TRUE)$maximum
}

## -----------------------------------------------------------------------------
estimate.em(us, vs, 1e6)
estimate.obs(us, vs)

## -----------------------------------------------------------------------------
solve.game <- function(A) {
    #solve the two player zero-sum game by simplex method
    #optimize for player 1, then player 2
    #maximize v subject to ...
    #let x strategies 1:m, and put v as extra variable
    #A1, the <= constraints
    #
    min.A <- min(A)
    A <- A - min.A #so that v >= 0
    max.A <- max(A)
    A <- A / max(A)
    m <- nrow(A)
    n <- ncol(A)
    it <- n^3
    a <- c(rep(0, m), 1) #objective function
    A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
    b1 <- rep(0, n)
    A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
    b3 <- 1
    sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
    maxi=TRUE, n.iter=it)
    #the ’solution’ is [x1,x2,...,xm | value of game]
    #
    #minimize v subject to ...
    #let y strategies 1:n, with v as extra variable
    a <- c(rep(0, n), 1) #objective function
    A1 <- cbind(A, rep(-1, m)) #constraints <=
    b1 <- rep(0, m)
    A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
    b3 <- 1
    sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
                  maxi=FALSE, n.iter=it)
    soln <- list("A" = A * max.A + min.A,
        "x" = sx$soln[1:m],
        "y" = sy$soln[1:n],
        "v" = sx$soln[m+1] * max.A + min.A)
    soln
}

## -----------------------------------------------------------------------------
A <- matrix(c(  0,-2,-2,3,0,0,4,0,0,
                2,0,0,0,-3,-3,4,0,0,
                2,0,0,3,0,0,0,-4,-4,
                -3,0,-3,0,4,0,0,5,0,
                0,3,0,-4,0,-4,0,5,0,
                0,3,0,0,4,0,-5,0,-5,
                -4,-4,0,0,0,5,0,0,6,
                0,0,4,-5,-5,0,0,0,6,
                0,0,4,0,0,5,-6,-6,0), 9, 9)

library(boot) #needed for simplex function

solve.game(A)

## -----------------------------------------------------------------------------
B <- A + 2

solve.game(B)

## -----------------------------------------------------------------------------
dim(c(1,2,3))

## -----------------------------------------------------------------------------
m <- matrix(1:9, ncol=3, nrow=3)
is.matrix(m)
is.array(m)

## -----------------------------------------------------------------------------
dfm <- data.frame(x = 1:3, y = c("1", "2", "3"))

as.matrix(dfm)

## -----------------------------------------------------------------------------
dfm <- data.frame()
ncol(dfm)
nrow(dfm)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
dfm <- data.frame(x = c(1, 2, 3), y = c(0, 2, 8))
lapply(dfm, scale01)

## -----------------------------------------------------------------------------
dfm <- data.frame(x = c(1, 2, 3), y = c("0", "2", "8"))

idx <- sapply(dfm, class) == "numeric"
lapply(dfm[idx], scale01)

## -----------------------------------------------------------------------------
dfm <- data.frame(x = c(1, 2, 3), y = c(0, 2, 8))

vapply(dfm, sd, numeric(1))

## -----------------------------------------------------------------------------
dfm <- data.frame(x = c(1, 2, 3), y = c("0", "2", "8"))

idx <- sapply(dfm, class) == "numeric"
vapply(dfm[idx], sd, numeric(1))

## -----------------------------------------------------------------------------
mc.r <- function(N, a, b, n) {
    xs <- numeric(N)
    ys <- numeric(N)

    for (i in 1:(N-1)) {
        xs[i+1] <- rbinom(1, n, ys[i])
        ys[i+1] <- rbeta(1, xs[i+1] + a, n - xs[i+1] + b)
    }

    list(x = xs, y = ys)
}

## -----------------------------------------------------------------------------
library(Rcpp)

mc.rcpp <- cppFunction("
List mc(int N, double a, double b, int n) {
    NumericVector xs(N);
    NumericVector ys(N);
    
    xs[0] = 0;
    ys[0] = 0;
    
    
    for (int i = 0; i < N-1; i++) {
        xs[i+1] = R::rbinom(n, ys[i]);
        ys[i+1] = R::rbeta(xs[i+1] + a, n - xs[i+1] + b);
    }
    
    return List::create(Named(\"x\") = xs, Named(\"y\") = ys);
}
")

## -----------------------------------------------------------------------------
library(microbenchmark)
microbenchmark(mc.r(1e4, 3, 3, 4), mc.rcpp(1e4, 3, 3, 4))

