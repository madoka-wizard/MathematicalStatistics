# Ионин Василий Андреевич (вариант № 14)
# Лабораторная работа № 2

library(Pareto)

############## Подзадача 1 ##############
# Задать объём выборки N
N <- 1000

# Задать параметры соответствующего распределения

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

# Найти аналитически функцию распределения F моделируемого закона и плотность распределения
# Вычислить математическое ожидание и дисперсию этого распределения

F_1 <- function(x) {
  return(pPareto(x, a, x_0))
}

F_2 <- function(x) {
  return(punif(x, min = l, max = r))
}

F <- function(x) {                                            # Функция распределения
  return(p_1 * F_1(x) + p_2 * F_2(x))
}

f_1 <- function(x) {
  return(dPareto(x, a, x_0))
}

f_2 <- function(x) {
  return(dunif(x, min = l, max = r))
}

f <- function(x) {                                            # Плотность распределения
  return(p_1 * f_1(x) + p_2 * f_2(x))
}

# Смоделировать с помощью генераторов случайных чисел выборку из заданного распределения
X <- sample(x = 1:2, N, replace = TRUE, prob = c(p_1, p_2))
Y <- numeric(N)
F_1_samples <- rPareto(N, a, x_0)
F_2_samples <- runif(N, l, r)
for (i in 1:N)
  Y[i] <- switch(X[i], F_1_samples[i], F_2_samples[i])

Y <- sort(Y, decreasing = FALSE)

############## Подзадача 2 ##############

# Построение графиков функции распределения F и эмпирической функции распределения F_n
from <- as.integer(Y[1]) - 1
to <- as.integer(Y[N]) + 1
t <- seq(from, to, by = 0.1)

counter <- numeric(N)
m <- round(1.72 * N^(1 / 3))
frequences <- numeric(m + 1)
h <- (to - from) / m
limits <- seq(from, to, by = h)
for (i in 2:(m - 1)) {
  j <- counter[i] + 1
  while (j <= N & Y[j] < limits[i])
    j <- j + 1
  counter[i] <- j - 1
  frequences[i - 1] <- counter[i] - counter[i - 1]
}
frequences <- frequences / (N * h)
hist(Y, freq = FALSE, breaks = m, main = "Плотность распределения", ylab = "",
     xlab = "x")
lines(limits, frequences, type = 'l', col = "blue")
lines(t, f(t), type = 'l', col = "green")

step <- N^-0.2
xs <- seq(from, to, by = 0.001)
l <- length(xs)

solver <- function(kerne) {
  f <- numeric(l)
  for (i in 1:l) {
    f[i] <- sum(kernel((Y - xs[i]) / step)) / (N * step)
  }
  plot(xs, f1)
  lines(xs, f1, type = 'l', col = "red")
}

kernel_epanechnikov <- function(x) {
  return(ifelse(abs(x) <= 1, (3 / 4) * (1 - x^2), 0))
}

solver(kernel_epanechnikv)

kernel_triangular <- function(x) {
  return(ifelse(abs(x) <= 1, 1 - abs(x), 0))
}

solver(kernel_triangular)

kernel_squared <- function(x) {
  return(ifelse(abs(x) <= 1, (15 / 16) * (1 - x^2), 0))
}

solver(kernel_squared)

kernel_normal <- function(x) {
  return(exp(-x^2 / 2) / sqrt(2 * pi))
}

solver(kernel_normal)

kernel_rectangular <- function(x) {
  return(ifelse(abs(x) <= 1, 1 / 2, 0))
}

solver(kernel_rectangular)

kernel_dexp <- function(x) {
  return(exp(-abs(x)) / 2)
}

solver(kernel_dexp)