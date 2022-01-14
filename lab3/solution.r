# Ионин Василий Андреевич (вариант № 14)
# Лабораторная работа № 3

N <- 18
Y <- c(180, 158, 190, 180, 170, 175, 174, 167, 199, 190, 165, 170, 174, 182, 183, 167, 172, 175)
x_0 <- rep(1, N)
x_1 <- c(1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1)
x_2 <- c(170, 170, 175, 180, 170, 178, 175, 174, 193, 170, 178, 165, 198, 183, 178, 176, 181, 176)
x_3 <- c(170, 155, 195, 170, 165, 179, 163, 167, 173, 166, 156, 160, 176, 167, 172, 164, 187, 164)
x_4 <- c(72, 65, 80, 70, 60, 53, 65, 56, 85, 58, 55, 65, 50, 69, 94, 68, 65, 85)

############## Подзадача 1 ##############
# Метод наименьших квадратов
X <- matrix(c(x_0, x_1, x_2, x_3, x_4), N, 5)
B <- t(X) %*% X
theta <- solve(B) %*% t(X) %*% Y

sq_optimizer <- function(x_1, x_2, x_3, x_4, y) {
  return(sum((y - (theta[1] +
    x_1 * theta[2] +
    x_2 * theta[3] +
    x_3 * theta[4] +
    x_4 * theta[5]))^2))
}

############## Подзадача 2 ##############
# Остаточная дисперсия
rss <- sq_optimizer(x_1, x_2, x_3, x_4, Y)
vr <- rss / (N - 5)

############## Подзадача 3 ##############
make_model <- function(x) {
  return(x %*% theta)
}

make_model(c(1, 1, 181, 169, 69))

############## Подзадача 4 ##############
# Ищем доверительные интервалы
p <- 0.95
q <- 1 - p

compute_intervals <- function(i) {
  return(c(theta[i] - qt(p = 1 - q / 2, df = N - 5) *
    vr *
    sqrt(solve(B)[i, i]), theta[i] + qt(p = 1 - q / 2, df = N - 5) *
    vr *
    sqrt(solve(B)[i, i])))
}

for (i in 1:5) {
  compute_intervals(i)
}

intervals_sigma <- function() {
  return(c(rss / qchisq(p = 1 - q / 2, df = N - 5),
           rss / qchisq(p = q / 2, df = N - 5)))
}

intervals_sigma()

intervals_estimate <- function() {
  return(c(ex %*% theta - qt(p = 1 - q / 2, df = N - 5) *
    vr *
    sqrt(ex %*% solve(B) %*% t(t(ex))),
           ex %*% theta + qt(p = 1 - q / 2, df = N - 5) *
             vr *
             sqrt(ex %*% solve(B) %*% t(t(ex)))))
}

intervals_estimate()

############## Подзадача 5 ##############
compute_importance <- function(i) {
  return(abs(theta[i] / (vr * sqrt(solve(B)[i, i]))) <= qt(p = 1 - q / 2, df = N - 5))
}

for (i in 1:5) {
  compute_importance(i)
}

############## Подзадача 6 ##############
# Вычисляем коэффициент детерминации
Y_mean <- mean(Y)
d <- 1 - vr / (sum((Y - Y_mean)^2))