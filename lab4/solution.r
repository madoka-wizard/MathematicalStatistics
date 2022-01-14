# Ионин Василий Андреевич (вариант № 14)
# Лабораторная работа № 4

library(CPAT)
library(Pareto)

pearson_criterion <- function(alpha, intervals_count, real_distribution, samples) {
  n <- length(samples)
  mask <- numeric(intervals_count)
  step <- ceiling(n / intervals_count)
  var <- sort(samples)
  mask[1] <- var[1]
  stepcount <- 1
  for (i in 2:(intervals_count - 1)) {
    if (step * (i - 1) <= n) {
      mask[i] <- var[step * (i - 1)]
      continue
    }
    mask[i] <- var[n] + stepcount
    stepcount <- stepcount + 1
  }
  mask[intervals_count] <- 1

  p <- numeric(intervals_count)
  p[1] <- real_distribution(mask[1])
  for (i in 2:intervals_count) {
    p[i] <- real_distribution(mask[i]) - real_distribution(mask[i - 1])
  }

  v <- numeric(intervals_count)
  stepcount <- 1
  for (i in 1:n) {
    if (stepcount == intervals_count) {
      v[stepcount] <- v[stepcount] + 1
    } else if (var[i] <= mask[stepcount]) {
      v[stepcount] <- v[stepcount] + 1
    } else {
      stepcount <- stepcount + 1
    }
  }
  return(sum((v - n * p)^2 / (n * p), na.rm = TRUE) <= qchisq(1 - alpha, intervals_count - 1))
}

kolmogorov_criterion <- function(alpha, real_distribution, samples) {
  z_alpha <- CPAT:::qkolmogorov(1 - alpha)
  var <- sort(samples)
  n <- length(samples)
  D_1 <- 0
  D_2 <- 0
  for (i in 1:n) {
    D_1 <- max(D_1, abs(i / n - real_distribution(var[i])))
    D_2 <- max(D_2, abs(real_distribution(var[i]) - (i - 1) / n))
  }
  return(max(D_1, D_2) <= z_alpha / sqrt(n))
}

Z <- rbinom(300, 1, 0.5)
Z_mean <- sum(Z) / 300

F_b <- function(p) {
  return(pbinom(p, 1, 0.5))
}

pearson_criterion(0.05, 10, F_b, Z)

F_p <- function(p) {
  return(ppois(p, Z_mean))
}

pearson_criterion(0.05, 10, F_p, Z)

N <- 1000

# Параметры смеси
p <- 0.2    # Коэффициент смешивания
p_1 <- 0.2
p_2 <- 1 - p_1

# Параметры для распределения Парето
a <- 4
x_0 <- 0.5  # Коэффициент масштаба

# Параметры для нормального распределения
l <- -1     # Левая граница интервала
r <- 1      # Правая граница интервала

X <- sample(x = 1:2, N, replace = TRUE, prob = c(p_1, p_2))
Y <- numeric(N)
F_1_samples <- rPareto(N, a, x_0)
F_2_samples <- runif(N, l, r)
for (i in 1:N)
  Y[i] <- switch(X[i], F_1_samples[i], F_2_samples[i])

F_1 <- function(x) {
  return(pPareto(x, a, x_0))
}

F_2 <- function(x) {
  return(punif(x, min = l, max = r))
}

F_m <- function(p) {
  return(p_1 * F_1(x) + p_2 * F_2(x))
}

F_n <- function(p) {
  return(pnorm(p))
}

pearson_criterion(0.05, 100, F_m, Y)
kolmogorov_criterion(0.05, F_m, Y)

pearson_criterion(0.05, 100, F_n, Y)
kolmogorov_criterion(0.05, F_n, Y)