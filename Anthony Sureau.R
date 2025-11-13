train <- as.data.frame(cbind(dataset$y, dataset$X))

#Modèle linéaire

library(leaps)
out <- regsubsets(data$y ~ . , data=train1, 
                  nbest=1, nvmax=10, really.big=FALSE) 
bss <- summary(out) 
bss.size <- as.numeric(rownames(bss$which))
intercept <- lm(data$y~ 1, data=train1) 
bss.best.rss <-
  c(sum(resid(intercept)^2), tapply(bss$rss, bss.size, min)) 
plot(0:10, bss.best.rss, ylim=c(30, 135), type="b", xlab="subset size", 
     ylab="RSS", col="red2" ) 
points(bss.size, bss$rss, pch=20, col="gray", cex=0.7) 

bss.best.cp <- tapply(bss$cp , bss.size, min) 
plot(1:10, bss.best.cp, type="b", xlab="subset size", ylab="Cp", col="red2" ) 
points(bss.size, bss$cp, pch=20, col="gray", cex=0.7)

bss.best.bic <- tapply(bss$bic , bss.size,min) 
plot(1:10, bss.best.bic, type="b", xlab="subset size", ylab="BIC",
     col="red2" )
points(bss.size, bss$bic, pch=20, col="gray", cex=0.7)

#Modèle non linéaire

library(leaps)
out <- regsubsets(dataset$y ~ . , data=train, 
                  nbest=1, nvmax=10, really.big=FALSE) 
bss <- summary(out) 
bss.size <- as.numeric(rownames(bss$which))
intercept <- lm(dataset$y~ 1, data=train) 
bss.best.rss <-
  c(sum(resid(intercept)^2), tapply(bss$rss, bss.size, min)) 
plot(0:10, bss.best.rss, ylim=c(30, 135), type="b", xlab="subset size", 
     ylab="RSS", col="red2" ) 
points(bss.size, bss$rss, pch=20, col="gray", cex=0.7) 

#RSS

library(leaps) 
out <- regsubsets(dataset$y ~ . , data=train, nbest=100, really.big=TRUE) 
bss <- summary(out)
bss.size <- as.numeric(rownames(bss$which))
intercept <- lm(dataset$y ~ 1, data=train)
bss.best.rss <- c(sum(resid(intercept)^2), tapply(bss$rss , bss.size,
                                                  min))
plot(0:8, bss.best.rss, type="b", 
     xlab="subset size", ylab="RSS", col="red2" ) 
points(bss.size, bss$rss, pch=20, col="gray", cex=0.7)

#C_p

bss.best.cp <- tapply(bss$cp , bss.size, min) 
plot(1:8, bss.best.cp, type="b", xlab="subset size", ylab="Cp", col="red2" ) 
points(bss.size, bss$cp, pch=20, col="gray", cex=0.7)

#BIC

bss.best.bic <- tapply(bss$bic , bss.size,min) 
plot(1:8, bss.best.bic, type="b", xlab="subset size", ylab="BIC",
     col="red2" )
points(bss.size, bss$bic, pch=20, col="gray", cex=0.7)


library(leaps)

set.seed(123)      # permet la reproductibilité
n.sim <- 30        # nombre de simulations
nvmax  <- 10       # nombre max de variables dans les sous-modèles

# vecteurs pour stocker les tailles optimales
best_cp_sizes  <- numeric(n.sim)
best_bic_sizes <- numeric(n.sim)

for (i in 1:n.sim) {
  #Génération données simulées
  data_i <- generate.lm(n.train = 100, p = 20, p0 = 5, sigma2 = 1)
  
  # jeu d'apprentissage
  train_i <- data.frame(
    y = as.numeric(data_i$y[data_i$train]),
    data_i$X[data_i$train, , drop = FALSE]
  )
  
  out_i <- regsubsets(y ~ ., data = train_i, nbest = 1, nvmax = nvmax, really.big = FALSE)
  bss_i <- summary(out_i)
  
  #Taille du meilleur modèle pour chaque critère
  best_cp_sizes[i]  <- which.min(bss_i$cp)
  best_bic_sizes[i] <- which.min(bss_i$bic)
}

#Résumé des résultats
cat("Moyenne des tailles optimales selon Cp :", mean(best_cp_sizes), "\n")
cat("Moyenne des tailles optimales selon BIC :", mean(best_bic_sizes), "\n")

#Visualisation
boxplot(best_cp_sizes, best_bic_sizes, names = c("Cp", "BIC"),
        col = c("skyblue", "tomato"),
        main = paste(n.sim, "simulations - Taille du meilleur modèle"),
        ylab = "Nombre de variables sélectionnées")