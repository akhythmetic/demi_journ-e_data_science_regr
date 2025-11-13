## =========================
## 0) Packages
## =========================
# install.packages("MASS") # au besoin
library(MASS)

## =========================
## 1) Génération (libre : X ~ N(0, I_p)) + split train/test
## =========================
generate.lm.libre <- function(n.train, n.test, p, p0, sigma2, seed = 123) {
  stopifnot(n.train > 0, n.test > 0, p > 0, p0 >= 0, p0 <= p, sigma2 >= 0)
  set.seed(seed)
  
  n <- n.train + n.test
  train <- 1:n.train
  test  <- (n.train + 1):n
  
  # Vrai beta sparse
  beta <- rep(0, p)
  S_star <- if (p0 > 0) sample.int(p, p0, replace = FALSE) else integer(0)
  if (length(S_star)) {
    beta[S_star] <- runif(p0, 1, 2) * sample(c(-1, 1), p0, replace = TRUE)
  }
  
  # X ~ N(0, I_p)
  sigmaMatrix <- diag(p)
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = sigmaMatrix)
  
  # Bruit et réponse
  eps <- rnorm(n, sd = sqrt(sigma2))
  y <- as.vector(X %*% beta + eps)
  
  list(y = y, X = X, beta = beta, sigmaMatrix = sigmaMatrix,
       S_star = S_star, train = train, test = test)
}

## =========================
## 2) Méthodes
## =========================

# 2.1 MCO (OLS, sans intercept car X centrées)
getOLS <- function(X, y) {
  df <- data.frame(y = as.vector(y), X)
  colnames(df) <- c("y", paste0("x", seq_len(ncol(X))))
  fit <- lm(y ~ . - 1, data = df)
  b <- coef(fit); b[is.na(b)] <- 0
  p <- ncol(X)
  beta_hat <- numeric(p); names(beta_hat) <- paste0("x", 1:p)
  beta_hat[match(names(b), names(beta_hat))] <- as.numeric(b)
  list(beta_hat = beta_hat, fit = fit, support = which(beta_hat != 0))
}

# 2.2 MCO + test de Student : sélection par p-values < alpha puis refit
getOLS_pselect <- function(X, y, alpha = 0.05) {
  df <- data.frame(y = as.vector(y), X); colnames(df) <- c("y", paste0("x", 1:ncol(X)))
  full <- lm(y ~ . - 1, data = df)
  pvals <- summary(full)$coefficients[, "Pr(>|t|)"]
  keep <- which(pvals < alpha)
  p <- ncol(X)
  if (length(keep) == 0) {
    return(list(beta_hat = rep(0, p), keep = integer(0), pvals = pvals))
  }
  Xs <- X[, keep, drop = FALSE]
  b_partial <- coef(lm(as.vector(y) ~ Xs - 1))
  beta_hat <- rep(0, p); beta_hat[keep] <- as.numeric(b_partial)
  list(beta_hat = beta_hat, keep = keep, pvals = pvals)
}

# 2.3 Stepwise BIC (backward) + barplot des coefficients retenus
getStepBIC_report <- function(X, y, plot = TRUE, main = "Stepwise BIC — Coefficients retenus") {
  df <- data.frame(y = as.vector(y), X); colnames(df) <- c("y", paste0("x", 1:ncol(X)))
  full <- lm(y ~ . - 1, data = df)
  k_bic <- log(nrow(X))
  st <- MASS::stepAIC(full, direction = "backward", k = k_bic, trace = FALSE)
  bic_val <- as.numeric(BIC(st))
  coefs <- coef(st)
  
  # Graphe barplot des coefficients retenus
  if (plot) {
    if (is.null(coefs)) {
      plot.new(); title(main = "Stepwise BIC — aucun prédicteur retenu")
    } else {
      ord <- order(coefs)
      op <- par(mar = c(7,4,4,1)); on.exit(par(op), add = TRUE)
      barplot(coefs[ord], las = 2, col = "gray85", border = NA,
              main = sprintf("%s\nBIC = %.2f (k = log(n) = %.2f)", main, bic_val, k_bic),
              ylab = "Coefficient")
      abline(h = 0, col = "gray40", lty = 2); box()
    }
  }
  
  # Remettre les coefficients en plein format (longueur p)
  p <- ncol(X)
  beta_hat <- rep(0, p); names(beta_hat) <- paste0("x", 1:p)
  if (!is.null(coefs)) {
    idx <- match(names(coefs), colnames(df)[-1])
    beta_hat[idx] <- as.numeric(coefs)
  }
  list(beta_hat = beta_hat, bic = bic_val, model = st)
}

## =========================
## 3) Métriques
## =========================
pred_mse <- function(X, y, beta_hat) {
  yhat <- as.vector(X %*% beta_hat)
  mean((as.vector(y) - yhat)^2)
}

## =========================
## 4) Démo complète (libre)
## =========================
set.seed(123)
  data_libre <- generate.lm.libre(n.train = 140, n.test = 60, p = 20, p0 = 5, sigma2 = 1)

# Découpage
Xtr <- data_libre$X[data_libre$train, , drop = FALSE]
ytr <- data_libre$y[data_libre$train]
Xte <- data_libre$X[data_libre$test,  , drop = FALSE]
yte <- data_libre$y[data_libre$test]

# --- MCO
ols <- getOLS(Xtr, ytr)

# --- Test de Student : p-values + refit
ols_p <- getOLS_pselect(Xtr, ytr, alpha = 0.05)

# Graphe des p-values (triées) + “volcano” simple
pvals <- ols_p$pvals
op <- par(mfrow = c(1,2), mar = c(7,4,3,1))
ord <- order(pvals)
plot(pvals[ord], pch = 16, xaxt = "n", ylab = "p-value", main = "Test de Student (p-values)")
axis(1, at = seq_along(ord), labels = paste0("x", ord), las = 2, cex.axis = 0.7)
abline(h = 0.05, col = 2, lty = 2)

plot(abs(ols$beta_hat), -log10(pvals), pch = 16,
     xlab = "|beta_hat (OLS)|", ylab = "-log10(p-value)",
     main = "Effet vs significativité")
abline(h = -log10(0.05), col = 2, lty = 2)
par(op)

# --- Stepwise BIC (trace un barplot des coefficients retenus)
sbic <- getStepBIC_report(Xtr, ytr, plot = TRUE)

# --- BIC de référence (modèle nul et modèle plein) pour situer
df_tr <- data.frame(y = ytr, Xtr); colnames(df_tr) <- c("y", paste0("x", 1:ncol(Xtr)))
mod_null <- lm(y ~ 0,     data = df_tr)     # sans intercept
mod_full <- lm(y ~ . - 1, data = df_tr)     # toutes variables
bic_null <- BIC(mod_null); bic_full <- BIC(mod_full)

print(c(BIC_null = bic_null, BIC_full = bic_full, BIC_step = sbic$bic))

# --- MSE test des trois approches
mse_ols  <- pred_mse(Xte, yte, ols$beta_hat)
mse_stud <- pred_mse(Xte, yte, ols_p$beta_hat)
mse_bic  <- pred_mse(Xte, yte, sbic$beta_hat)

print(c(MSE_OLS = mse_ols, MSE_Student = mse_stud, MSE_StepBIC = mse_bic))

# (option) nombre de variables retenues
cat("Vars retenues — Student:", length(ols_p$keep),
    " | StepBIC:", sum(sbic$beta_hat != 0), "\n")
