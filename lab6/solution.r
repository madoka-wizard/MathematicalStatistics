# Ионин Василий Андреевич (вариант № 14)
# Лабораторная работа № 6

spearman_criterion <- function(alpha, X, Y) {
  stopifnot(length(X) == length(Y))
  n <- length(X)
  z_alpha <- qnorm(1 - alpha / 2)
  X_rank <- rank(X)
  Y_rank <- rank(Y)
  rho <- 1 - (6 / (n * (n^2 - 1))) * sum((X_rank - Y_rank)^2)
  return(sqrt(n) * abs(rho) > z_alpha)
}

wilcoxon_criterion <- function(alpha, X, Y) {
  n <- length(X)
  m <- length(Y)
  z_alpha <- qnorm(1 - alpha / 2)
  Z_rank <- rank(c(X, Y))
  T_nm <- sum(Z_rank[1:n])
  return(abs((T_nm - (n / 2) * (n + m + 1)) / (sqrt(m * n * (m + n + 1) / 12))) > z_alpha)
}


N <- 1000
M <- 800
eps <- 0.05

# Независимые величины
X <- rnorm(N)
Y <- rnorm(M)
spearman_criterion(eps, X, Y)
wilcoxon_criterion(eps, X, Y)

# Алгебраическая зависимость
X <- rnorm(N)
Y <- sin(X) + 2 * X
spearman_criterion(eps, X, Y)

# Величина и смесь
X <- rnorm(N)
p <- 0.2
p_1 <- p
p_2 <- 1 - p
scheme <- sample(x = 1:2, M, replace = TRUE, prob = c(p_1, p_2))
unif_samples <- runif(M, 4, 3.5)
norm_samples <- rnorm(M)
for (i in 1:M) {
  Y[i] <- switch(scheme[i], unif_samples[i], norm_samples[i])
}
wilcoxon_criterion(eps, X, Y)
