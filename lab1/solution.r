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

############## Подзадача 7 ##############
# Вычислить значения оценок мат. ожидания и дисперсии по полученной выборке
emperical_mean <- mean(Y)
emperical_variance <- var(Y)

############## Подзадача 8 ##############
# Построить асимптотические доверительные интервалы для "неизвестных" мат. ожидания и дисперсии
confidence_mean <- function(alpha) {
  l <- emperical_mean - qnorm(1 - alpha / 2) * sqrt(emperical_variance / N)
  r <- emperical_mean + qnorm(1 - alpha / 2) * sqrt(emperical_variance / N)
  return(c(l, r))
}

confidence_variance <- function(alpha) {
  M_4 <- sum((Y - emperical_mean)^4) / N
  l <- emperical_variance - qnorm(1 - alpha / 2) * sqrt(M_4 - emperical_variance^2) / sqrt(N)
  r <- emperical_variance + qnorm(1 - alpha / 2) * sqrt(M_4 - emperical_variance^2) / sqrt(N)
  return(c(l, r))
}

alpha_1 <- 0.1
mean_interval_1 <- confidence_mean(alpha_1)
variance_interval_1 <- confidence_variance(alpha_1)

alpha_2 <- 0.05
mean_interval_2 <- confidence_mean(alpha_2)
variance_interval_2 <- confidence_variance(alpha_2)

alpha_3 <- 0.0026
mean_interval_3 <- confidence_mean(alpha_3)
variance_interval_3 <- confidence_variance(alpha_3)

############## Подзадача 9 ##############
# Статистические оценки парметров асимметрии и эксцесса
B_N <- mean((Y - emperical_mean)^3) / emperical_variance^(3 / 2)
E_N <- mean((Y - emperical_mean)^4) / emperical_variance^2 - 3