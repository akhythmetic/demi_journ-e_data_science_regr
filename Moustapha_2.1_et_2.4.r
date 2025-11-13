---
  title: "Demi-journée Data Science : Evaluation des méthodes de sélection de variable"
author: "M1 MIASHS 2025"
date: "2025"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## 2.1 Modèle linéaire (cas simple)

```{r}
library(MASS)

generate.lm <- function(n, p, p0, sigma2, seed = 123) {
  set.seed(seed)
  beta <- rep(0, p)
  S_star <- sample(1:p, p0)
  beta[S_star] <- runif(p0, 1, 2) * sample(c(-1, 1), p0, replace = TRUE)
  sigmaMatrix <- diag(p)
  X <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigmaMatrix)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2))
  y <- X %*% beta + epsilon
  list(y = y, X = X, beta = beta, sigmaMatrix = sigmaMatrix)
}

data <- generate.lm(50, 10, 5, 60)
data
```

## 2.2 Structure de dépendance des prédicteurs

```{r}
# install.packages("mvtnorm") # si pas déjà installé
library(mvtnorm)

# Dépendance longitudinale (AR(1)) : Sigma[i,j] = rho^|i-j|
generate.lm.long <- function(n.train, p, p0, sigma2, rho, n.test = 10 * n.train) {
  stopifnot(n.train > 0, p > 0, p0 >= 0, p0 <= p, sigma2 >= 0)
  stopifnot(is.numeric(rho), abs(rho) < 1)
  
  n <- n.train + n.test
  train <- 1:n.train
  test  <- (n.train + 1):n
  
  idx <- 0:(p - 1)
  sigmaMatrix <- outer(idx, idx, function(i, j) rho^abs(i - j))
  
  beta <- numeric(p)
  if (p0 > 0) {
    S.star <- sample.int(p, p0, replace = FALSE)
    beta[S.star] <- runif(p0, 1, 2) * sample(c(-1, 1), p0, replace = TRUE)
  }
  
  X <- rmvnorm(n, mean = rep(0, p), sigma = sigmaMatrix)
  noise <- rnorm(n, sd = sqrt(sigma2))
  y <- as.vector(X %*% beta + noise)
  
  list(y = y, X = X, beta = beta, sigma2 = sigma2, sigmaMatrix = sigmaMatrix, train = train, test = test)
}

set.seed(1)
dataset <- generate.lm.long(n.train = 100, p = 20, p0 = 5, sigma2 = 1, rho = 0.6, n.test = 200)

dim(dataset$X)
sum(dataset$beta != 0)
round(dataset$sigmaMatrix[1:4,1:4], 3)

df_tr <- data.frame(y = dataset$y[dataset$train], dataset$X[dataset$train, ])
colnames(df_tr) <- c("y", paste0("x", 1:ncol(dataset$X)))
fit <- lm(y ~ . - 1, data = df_tr)
summary(fit)

Xte <- dataset$X[dataset$test, , drop = FALSE]
yte <- dataset$y[dataset$test]
yhat <- as.vector(Xte %*% coef(fit))
mean((yte - yhat)^2)
```




## 2.4 Implémentation des méthodes de sélection de modèle

# Chargement des packages
library(MASS)
library(leaps)
library(ggplot2)
library(gridExtra)

# Fonctions de sélection : Stepwise AIC, Stepwise BIC, Best Subset
getStepAIC <- function(X, y) {
  df <- data.frame(y = y, X)
  mod <- lm(y ~ ., data = df)
  step_aic <- step(mod, direction = "both", trace = FALSE)
  coef(step_aic)
}

getStepBIC <- function(X, y) {
  df <- data.frame(y = y, X)
  mod <- lm(y ~ ., data = df)
  step_bic <- step(mod, direction = "both", k = log(nrow(X)), trace = FALSE)
  coef(step_bic)
}

getBestSubset <- function(X, y, crit = "bic") {
  df <- data.frame(y = y, X)
  regfit <- regsubsets(y ~ ., data = df, nvmax = ncol(X), really.big = TRUE)
  regsum <- summary(regfit)
  
  if (crit == "bic") best <- which.min(regsum$bic)
  else if (crit == "cp") best <- which.min(regsum$cp)
  else if (crit == "adjr2") best <- which.max(regsum$adjr2)
  else stop("Critère inconnu.")
  
  coef(regfit, best)
}

# Génération des données
data <- generate.lm(100, 10, 4, 10)
X <- data$X
y <- data$y
df <- data.frame(y, X)

# Best Subset : exploration des modèles
regfit.full <- regsubsets(y ~ ., data = df, nvmax = ncol(X))
summary.full <- summary(regfit.full)
names(summary.full)

# Graphiques des critères de sélection
par(mfrow = c(2, 2))

plot(summary.full$cp, xlab = "Nombre de variables", ylab = "Cp",
     type = "b", pch = 19, col = "darkred")
points(which.min(summary.full$cp), min(summary.full$cp), col = "blue", pch = 19)

plot(summary.full$bic, xlab = "Nombre de variables", ylab = "BIC",
     type = "b", pch = 19, col = "darkgreen")
points(which.min(summary.full$bic), min(summary.full$bic), col = "blue", pch = 19)

plot(summary.full$adjr2, xlab = "Nombre de variables", ylab = "R² ajusté",
     type = "b", pch = 19, col = "purple")
points(which.max(summary.full$adjr2), max(summary.full$adjr2), col = "blue", pch = 19)

title("Critères de sélection – Best Subset")

# Best Subset : modèle optimal selon BIC
best.bic <- which.min(summary.full$bic)
coef.best.bic <- coef(regfit.full, best.bic)

# Stepwise AIC et BIC
mod_full <- lm(y ~ ., data = df)
mod_AIC <- step(mod_full, direction = "both", trace = FALSE)
mod_BIC <- step(mod_full, direction = "both", k = log(nrow(df)), trace = FALSE)

# Comparaison des coefficients estimés
beta_AIC <- coef(mod_AIC)
beta_BIC <- coef(mod_BIC)
beta_BEST <- coef.best.bic

comp <- data.frame(
  beta_vrai = round(data$beta, 2),
  beta_AIC  = 0,
  beta_BIC  = 0,
  beta_BEST = 0
)

names(comp) <- c("beta_vrai", "beta_AIC", "beta_BIC", "beta_BEST")

for (j in names(beta_AIC)[-1])  comp[j, "beta_AIC"]  <- beta_AIC[j]
for (j in names(beta_BIC)[-1])  comp[j, "beta_BIC"]  <- beta_BIC[j]
for (j in names(beta_BEST)[-1]) comp[j, "beta_BEST"] <- beta_BEST[j]

comp

# Diagnostics du modèle final
par(mfrow = c(2, 2))
plot(mod_BIC, which = 1)



#3.1 Évaluation des performances


perf <- function(X_test, y_test, beta, beta.star) {
  
  nzero <- which(beta != 0)
  zero  <- which(beta == 0)
  
  true.nzero <- which(beta.star != 0)
  true.zero  <- which(beta.star == 0)
  
  TP <- sum(nzero %in% true.nzero)
  TN <- sum(zero %in%  true.zero)
  FP <- sum(nzero %in% true.zero)
  FN <- sum(zero %in%  true.nzero)
  
  recall    <- TP/(TP + FN) ## also recall and sensitivity
  specificity   <- TN/(FP + TN) ## specificity
  precision <- TP/(TP + FP) ## also PPR
  recall[TP + FN == 0] <- NA
  specificity[TN + FP == 0] <- NA
  precision[TP + FP == 0] <- NA
  
  rmse <- sqrt(mean((beta - beta.star)^2, na.rm = TRUE))
  rerr <- sqrt(mean((y_test - X_test %*% beta)^2))
  res <-  round(c(precision,recall,specificity, rmse, rerr),4)
  res[is.nan(res)] <- 0
  names(res) <- c("precision","recall","specificity","rmse", "prediction") 
  res
}

#3.2 Planning de simulations

library(tibble)

res <- tribble(
  ~method,       ~mse, ~err,  ~acc, ~sen, ~spe, ~n.p, ~sigma2, ~simu,
  "bestsubset",   0.3, -0.25, 0.92, 0.9,  0.75, 0.5,   0.75,     1,
  "stepwiseAIC",  1.65, -0.17, 0.92, 0.9,  0.75, 0.5,   0.75,     1,
  "stepwiseBIC",  0.51,  0.07, 0.92, 0.9,  0.75, 0.5,   0.75,     1
)

res

library(tidyverse) # pour l'utilisation du pipe %>%
getOneSimu <- function(i) {
  data <- generate.lm(n.train=100, p=50, p0=10, sigma2=1, n.test=10*n.train100)
  
  beta_AIC <- getStepAIC(data$X[data$train, ], data$y[data$train])
  beta_BIC <- getStepBIC(data$X[data$train, ], data$y[data$train])
  
  res <-
    data.frame(
      rbind(
        perf(data$X[data$test, ], data$y[data$test], beta_AIC, data$beta),
        perf(data$X[data$test, ], data$y[data$test], beta_BIC, data$beta)
      )
    ) %>% 
    add_column(method = c("stepAIC", "stepBIC")) %>% 
    add_column(simu_label = i)
  res
}

library(tidyverse) # pour la fonction Reduce()
library(pbmcapply)
n_simu <- 10
res <- Reduce("rbind", pbmclapply(1:n_simu, getOneSimu, mc.cores = 2))

getOneSimu <- function(i) {
  data <- generate.lm(n = 100, p = 50, p0 = 10, sigma2 = 1)
  
  beta_AIC  <- getStepAIC(data$X, data$y)
  beta_BIC  <- getStepBIC(data$X, data$y)
  beta_BEST <- getBestSubset(data$X, data$y, crit = "bic")
  
  res <-
    data.frame(
      rbind(
        perf(data$X, data$y, beta_AIC, data$beta),
        perf(data$X, data$y, beta_BIC, data$beta),
        perf(data$X, data$y, beta_BEST, data$beta)
      )
    ) %>%
    add_column(method = c("stepAIC", "stepBIC", "bestSubset")) %>%
    add_column(simu_label = i)
  res
}

getOneSimu(1)

