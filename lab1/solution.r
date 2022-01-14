# Ионин Василий Андреевич (вариант № 14)
# Лабораторная работа № 1

library(Pareto)

############## Подзадача 1 ##############
# Задать объём выборки N
N <- 1000

############## Подзадача 2 ##############
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

############## Подзадача 3 ##############
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

F_1_mean <- a * x_0 / (a - 1)
F_2_mean <- (l + r) / 2

real_mean <- p_1 * F_1_mean + p_2 * F_2_mean                  # Теоретическое среднее

F_1_variance <- (x_0 / (a - 1))^2 * a / (a - 2)
F_2_variance <- (r - l)^2 / 12

real_variance <- p_1^2 * F_1_variance + p_2^2 * F_2_variance  # Теоретическая дисперсия

############## Подзадача 4 ##############
# Смоделировать с помощью генераторов случайных чисел выборку из заданного распределения
X <- sample(x = 1:2, N, replace = TRUE, prob = c(p_1, p_2))
Y <- numeric(N)
F_1_samples <- rPareto(N, a, x_0)
F_2_samples <- runif(N, l, r)
for (i in 1:N)
  Y[i] <- switch(X[i], F_1_samples[i], F_2_samples[i])

############## Подзадача 5 ##############
# Построить графики F и эмпирической функции распределения F_n
F_n <- ecdf(Y)

from <- as.integer(min(Y))
to <- as.integer(max(Y))
t <- seq(from, to, by = 0.1)

plot(t, F(t), main = "F - зеленый, F_n - синий", ylab = "", xlab = "x", type = 'l', col = "green")
lines(t, F_n(t), type = 'l', col = "blue")

############## Подзадача 6 ##############
# Построить график плотности распределения и гистограмму относительных частот
m <- round(1.72 * N^(1 / 3))
counter <- numeric(m)
frequences <- numeric(m + 1)
h <- (to - from) / m
limits <- seq(from, to, by = h)
Y <- sort(Y)

for (i in 2:m - 1) {
  j <- counter[i] + 1
  while (j <= N & Y[j] < limits[i])
    j <- j + 1
  counter[i] <- j - 1
  frequences[i - 1] <- counter[i] - counter[i - 1]
}

frequences <- frequences / (N * h)
hist(Y, freq = FALSE, breaks = m, main = "Гистограмма и плотность", ylab = "", xlab = "x", col = "blue")
lines(x, f(x), type = 'l', col = "green")