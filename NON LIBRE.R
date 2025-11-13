# install.packages("mvtnorm") # si pas déjà installé
library(mvtnorm)
library(MASS)

# Dépendance longitudinale (AR(1)) : Sigma[i,j] = rho^|i-j|
generate.lm.long <- function(n.train, p, p0, sigma2, rho, n.test = 10 * n.train) {
  stopifnot(n.train > 0, p > 0, p0 >= 0, p0 <= p, sigma2 >= 0)
  stopifnot(is.numeric(rho), abs(rho) < 1)  # pour garantir Σ définie positive
  
  n <- n.train + n.test
  train <- 1:n.train
  test  <- (n.train + 1):n
  
  # 1) Matrice de covariance AR(1)
  idx <- 0:(p - 1)
  sigmaMatrix <- outer(idx, idx, function(i, j) rho^abs(i - j))
  
  # 2) Vrai support et coefficients (|β| ~ U[1,2], signe ±)
  beta <- numeric(p)
  if (p0 > 0) {
    S.star <- sample.int(p, p0, replace = FALSE)
    beta[S.star] <- runif(p0, 1, 2) * sample(c(-1, 1), p0, replace = TRUE)
  }
  
  # 3) Génération de X ~ N(0, Σ) et du bruit ε ~ N(0, sigma2)
  X <- rmvnorm(n, mean = rep(0, p), sigma = sigmaMatrix)
  noise <- rnorm(n, sd = sqrt(sigma2))
  
  # 4) Réponse
  y <- as.vector(X %*% beta + noise)
  
  list(
    y = y,
    X = X,
    beta = beta,
    sigma2 = sigma2,
    sigmaMatrix = sigmaMatrix,
    train = train,
    test = test
  )
}
set.seed(123)
dataset <- generate.lm.long(n.train = 140, p = 20, p0 = 5, sigma2 = 1, rho = 0.6, n.test = 60)
dataset

# Vérifs
dim(dataset$X)                 # 300 x 20
sum(dataset$beta != 0)         # 5
round(dataset$sigmaMatrix[1:4,1:4], 3)  # vérifie la structure AR(1)
# ---- DIAGNOSTICS RAPIDES ----
Xtr <- dataset$X[dataset$train, , drop = FALSE]
ytr <- dataset$y[dataset$train]
Xte <- dataset$X[dataset$test,  , drop = FALSE]
yte <- dataset$y[dataset$test]
beta.star <- dataset$beta

# Refit des méthodes (assure-toi d'utiliser les fonctions corrigées)
ols   <- getOLS(Xtr, ytr)
ols_p <- getOLS_pselect(Xtr, ytr, alpha = 0.05)
sbic  <- getStepBIC_report(Xtr, ytr, plot = FALSE)  # plot = FALSE juste pour debug

# 1) Affiche les vrais coefficients et indices
cat("VRAI beta (non nuls) indices:", which(abs(beta.star) > 1e-6), "\n")
print(round(beta.star[which(abs(beta.star) > 1e-6)], 4))

# 2) Supports estimés
sup_ols    <- which(abs(ols$beta_hat) > 1e-6)
sup_stud   <- which(abs(ols_p$beta_hat) > 1e-6)
sup_sbic   <- which(abs(sbic$beta_hat) > 1e-6)

cat("Support OLS   :", sup_ols, "\n")
cat("Support Stud. :", sup_stud, "\n")
cat("Support StepB :", sup_sbic, "\n")

# 3) Compare vecteurs béta estimés (are they identical?)
cat("ols_p == sbic (all equal?) :", all.equal(as.numeric(ols_p$beta_hat), as.numeric(sbic$beta_hat)), "\n")
cat("ols == ols_p (all equal?)  :", all.equal(as.numeric(ols$beta_hat), as.numeric(ols_p$beta_hat)), "\n")

# 4) Affiche coefficients non nuls pour chaque méthode (value + index)
cat("\nCoefficients non nuls (OLS):\n"); print(round(ols$beta_hat[sup_ols], 4))
cat("\nCoefficients non nuls (Student):\n"); print(round(ols_p$beta_hat[sup_stud], 4))
cat("\nCoefficients non nuls (StepBIC):\n"); print(round(sbic$beta_hat[sup_sbic], 4))

# 5) Affiche p-values du modèle full (pour inspection)
df_full <- data.frame(y = ytr, Xtr); colnames(df_full) <- c("y", paste0("x", 1:ncol(Xtr)))
full_fit <- lm(y ~ . - 1, data = df_full)
pv <- summary(full_fit)$coefficients[, "Pr(>|t|)"]
cat("\nP-values (full model) :\n"); print(round(pv, 6))

# 6) Si ols_p and sbic are identical, show which elements equal
if (isTRUE(all.equal(as.numeric(ols_p$beta_hat), as.numeric(sbic$beta_hat)))) {
  cat("\nols_p and sbic beta_hat are exactly equal. Showing their vector:\n")
  print(round(ols_p$beta_hat, 6))
}

# 7) Quick re-run with smaller train to see variability (optional)
cat("\n--- Optional quick rerun with n.train = 120 (more noise) ---\n")
set.seed(999)
tmp <- generate.lm.long(n.train = 120, p = 20, p0 = 6, sigma2 = 1, rho = 0.6, n.test = 60)
Xtr2 <- tmp$X[tmp$train, , drop = FALSE]; ytr2 <- tmp$y[tmp$train]
Xte2 <- tmp$X[tmp$test, , drop = FALSE];  yte2 <- tmp$y[tmp$test]
o2   <- getOLS(Xtr2, ytr2)
s2   <- getOLS_pselect(Xtr2, ytr2)
b2   <- getStepBIC_report(Xtr2, ytr2, plot = FALSE)
out2 <- rbind(
  OLS = perf(Xte2, yte2, o2$beta_hat, tmp$beta),
  Student = perf(Xte2, yte2, s2$beta_hat, tmp$beta),
  StepBIC = perf(Xte2, yte2, b2$beta_hat, tmp$beta)
)
print(out2)

# Ajustement (sans intercept, X centré autour de 0)
df_tr <- data.frame(y = dataset$y[dataset$train], dataset$X[dataset$train, ])
colnames(df_tr) <- c("y", paste0("x", 1:ncol(dataset$X)))
fit <- lm(y ~ . - 1, data = df_tr)
summary(fit)

# MSE test
Xte <- dataset$X[dataset$test, , drop = FALSE]
yte <- dataset$y[dataset$test]
yhat <- as.vector(Xte %*% coef(fit))
mean((yte - yhat)^2)

## =========================
## 2.4 — MÉTHODES : MCO & ORACLE
## =========================

## ---- MCO (moindres carrés ordinaires, sans intercept) ----
getOLS <- function(X, y) {
  df <- data.frame(y = y, X)
  colnames(df) <- c("y", paste0("x", seq_len(ncol(X))))
  fit <- lm(y ~ . - 1, data = df)               # -1 : pas d'intercept, cohérent avec X centré
  b <- coef(fit); b[is.na(b)] <- 0
  
  p <- ncol(X)
  beta_hat <- numeric(p)
  names(beta_hat) <- paste0("x", seq_len(p))
  idx <- match(names(b), names(beta_hat))
  beta_hat[idx] <- as.numeric(b)
  
  list(beta_hat = beta_hat, support = which(abs(beta_hat) > 1e-8), fit = fit)
}

getStepBIC_report <- function(X, y, plot = TRUE, main = "Stepwise (BIC) — Coefficients retenus") {
  # 1) Modèle plein (sans intercept)
  df <- data.frame(y = y, X)
  colnames(df) <- c("y", paste0("x", seq_len(ncol(X))))
  full <- lm(y ~ . - 1, data = df)
  
  # 2) Stepwise backward avec pénalité BIC
  k_bic <- log(nrow(X))
  st <- stepAIC(full, direction = "backward", k = k_bic, trace = FALSE)
  
  # 3) BIC du modèle final
  bic_val <- as.numeric(BIC(st))
  
  # 4) Récupérer les coefficients sélectionnés (nommés)
  coefs <- coef(st)
  # coefs peut être NULL si aucun prédicteur retenu (rare mais possible)
  if (is.null(coefs)) {
    if (plot) {
      plot.new()
      title(main = "Stepwise (BIC) — Aucun prédicteur retenu")
    }
    return(bic_val)
  }
  
  # 5) Graphe : barplot des coefficients sélectionnés (avec signe)
  if (plot) {
    op <- par(mar = c(6, 4, 4, 2) + 0.1)  # marges pour étiquettes
    on.exit(par(op), add = TRUE)
    
    ord <- order(coefs)  # tri pour une lecture plus claire
    bp <- barplot(coefs[ord],
                  las = 2, # étiquettes verticales
                  main = sprintf("%s\nBIC = %.3f | k = log(n) = %.2f", main, bic_val, k_bic),
                  ylab = "Coefficient", col = "gray85", border = NA)
    abline(h = 0, lty = 2, col = "gray40")
    box()
  }
  
  # 6) Retourne le BIC (numérique)
  bic_val
}
Xtr <- dataset$X[dataset$train, , drop = FALSE]
ytr <- dataset$y[dataset$train]

bic_num <- getStepBIC_report(Xtr, ytr)  # trace le graphe et renvoie le BIC numérique
bic_num

perf <- function(X_test, y_test, beta, beta.star, tol = 1e-6) {
  nzero      <- which(abs(beta)      > tol)
  zero       <- which(abs(beta)      <= tol)
  true.nzero <- which(abs(beta.star) > tol)
  true.zero  <- which(abs(beta.star) <= tol)
  
  TP <- sum(nzero %in% true.nzero)
  TN <- sum(zero  %in% true.zero)
  FP <- sum(nzero %in% true.zero)
  FN <- sum(zero  %in% true.nzero)
  
  recall      <- ifelse(TP + FN == 0, NA, TP/(TP + FN))
  specificity <- ifelse(TN + FP == 0, NA, TN/(TN + FP))
  precision   <- ifelse(TP + FP == 0, NA, TP/(TP + FP))
  
  rmse <- sqrt(mean((beta - beta.star)^2))
  rerr <- sqrt(mean((y_test - X_test %*% beta)^2))
  
  res <- round(c(precision, recall, specificity, rmse, rerr), 4)
  res[is.nan(res)] <- 0
  names(res) <- c("precision","recall","specificity","rmse","prediction")
  res
}
getStepBIC_report <- function(X, y, plot = TRUE, main = "Stepwise (BIC) — Coefficients retenus") {
  df <- data.frame(y = y, X)
  colnames(df) <- c("y", paste0("x", seq_len(ncol(X))))
  full <- lm(y ~ . - 1, data = df)
  
  k_bic <- log(nrow(X))
  st <- stepAIC(full, direction = "backward", k = k_bic, trace = FALSE)
  
  bic_val <- as.numeric(BIC(st))
  coefs <- coef(st)  # peut être NULL
  
  # Remise en format complet (longueur p)
  p <- ncol(X)
  beta_hat <- rep(0, p); names(beta_hat) <- paste0("x", 1:p)
  if (!is.null(coefs)) {
    # harmonise les noms si besoin
    nm <- names(coefs)
    nm <- gsub(".*\\$", "", nm)   # enlève "df$" ou "dataset$" si jamais
    nm <- gsub("^X", "x", nm)     # X1 -> x1
    idx <- match(nm, names(beta_hat))
    beta_hat[idx] <- as.numeric(coefs)
  }
  
  if (plot) {
    if (is.null(coefs)) {
      plot.new(); title(main = "Stepwise (BIC) — Aucun prédicteur retenu")
    } else {
      ord <- order(coefs); op <- par(mar = c(6,4,4,2) + 0.1); on.exit(par(op), add = TRUE)
      barplot(coefs[ord], las = 2,
              main = sprintf("%s\nBIC = %.3f | k = log(n) = %.2f", main, bic_val, k_bic),
              ylab = "Coefficient", col = "gray85", border = NA)
      abline(h = 0, lty = 2, col = "gray40"); box()
    }
  }
  
  list(beta_hat = beta_hat, bic = bic_val, model = st)
}
# Split
Xtr <- dataset$X[dataset$train, , drop = FALSE]
ytr <- dataset$y[dataset$train]
Xte <- dataset$X[dataset$test,  , drop = FALSE]
yte <- dataset$y[dataset$test]
beta.star <- dataset$beta

# Fits
ols   <- getOLS(Xtr, ytr)
ols_p <- getOLS_pselect(Xtr, ytr, alpha = 0.05)
sbic  <- getStepBIC_report(Xtr, ytr, plot = TRUE)

# Perfs (⚠️ StepBIC doit utiliser sbic$beta_hat, pas ols_p$beta_hat)
perf_OLS     <- perf(Xte, yte, ols$beta_hat,   beta.star)
perf_Student <- perf(Xte, yte, ols_p$beta_hat, beta.star)
perf_BIC     <- perf(Xte, yte, sbic$beta_hat,  beta.star)

res <- rbind(OLS = perf_OLS,
             Student = perf_Student,
             StepBIC = perf_BIC)
print(res)
cat("BIC_step =", round(sbic$bic, 3), "\n")
