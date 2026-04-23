#install.packages("igraph") 
#install.packages('blockmodeling')
library(igraph)
library(igraphdata)

nodes <- read.csv("political-books-nodes.csv")
edges <- read.csv("political-books-edges.csv")
# проверим
head(nodes)
head(edges)

G <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
length(V(G)$political_ideology)
vcount(G)


##########################################
############ пункт 4 #####################
##########################################


# ассортативность

# Ассортативность по степени вершины
r_degree <- assortativity_degree(G, directed = FALSE)
cat("Ассортативность по степени:", r_degree, "\n")

# Категориальная ассортативность (по идеологии)
r_type <- round(assortativity_nominal(G, as.numeric(as.factor(V(G)$political_ideology))), 2)
cat("Ассортативность по типу (идеология):", r_type, "\n")



# Гомофилия

V(G)$name <- seq_len(vcount(G))
edgeList <- as_edgelist(G)
faction <- matrix(c(V(G)$name, as.numeric(factor(V(G)$political_ideology))),
                  nrow = vcount(G), ncol = 2)

FromLabel <- as.numeric(faction[match(as.numeric(edgeList[,1]), as.numeric(faction[,1])), 2])
ToLabel   <- as.numeric(faction[match(as.numeric(edgeList[,2]), as.numeric(faction[,1])), 2])
edgeType  <- FromLabel + ToLabel

m <- ecount(G)
p <- 2 * m / (vcount(G) * (vcount(G) - 1))
group_sizes <- table(as.numeric(faction[,2]))

exp_dyad <- group_sizes[1] * (group_sizes[1] - 1) * p / 2
obs_dyad <- sum(FromLabel == 1 & ToLabel == 1)
dyadicity <- obs_dyad / exp_dyad

exp_hetero <- group_sizes[1] * group_sizes[2] * p
obs_hetero <- sum(FromLabel != ToLabel)
heterophilicity <- obs_hetero / exp_hetero

cat("Dyadicity:", round(dyadicity,3),
    "Heterophilicity:", round(heterophilicity,3), "\n")



##########################################
############ пункт 5 #####################
##########################################


# Компоненты связности
cl <- components(G)
cat("Число компонент связности:", cl$no, "\n") # граф полностью связный 

par(mfrow=c(1,2), mar=c(0,0,1,0))

ceb <- cluster_edge_betweenness(G)
fg  <- cluster_fast_greedy(G)

L <- layout_with_fr(G, niter=1000, area=vcount(G)^2.3)

plot(ceb, G, layout=L, vertex.size=5, vertex.label=NA, edge.arrow.size=0.2, main="Edge Betweenness")
plot(fg,  G, layout=L, vertex.size=5, vertex.label=NA, edge.arrow.size=0.2, main="Fast-Greedy")

cat("Компонент связности:", components(G)$no, "\n")
cat("Модульность (Betweenness):", modularity(ceb), "\n")
cat("Модульность (Fast-Greedy):", modularity(fg), "\n")


# интерпретация кластеров

table(membership(ceb), V(G)$political_ideology)
table(membership(fg),  V(G)$political_ideology)



# block modeling

# бинарная модель

library(igraph)
library(blockmodeling)

# Матрица смежности (0/1)
A <- as.matrix(as_adjacency_matrix(G))

# Binary blockmodeling (структурная эквивалентность)
res_bin <- optRandomParC(
  M = A,
  k = 3,                # три блока (идеологии)
  approaches = "bin",   # бинарный блокмоделинг
  blocks = c("nul", "com"),
  rep = 100,
  nCores = 0
)

# Визуализация
plotMat(M = A, clu = clu(res_bin), main = "Binary Blockmodeling")
round(funByBlocks(res_bin), 2) # средние плотности внутри блоков


# логистическая регрессия 

library(blockmodeling)

# Если веса в графе отсутствуют — используем матрицу 0/1
A <- as.matrix(as_adjacency_matrix(G))

# Выбираем порог (например, медиану ненулевых связей)
m <- median(A[A > 0])

# Valued blockmodeling
res_val <- optRandomParC(
  M = A,
  k = 3,
  rep = 100,
  preSpecM = m,         # аналог "порогового" разделения
  approach = "val",     # valued (логистический)
  blocks = c("nul", "com"),
  nCores = 0
)

# Визуализация и блоковые плотности
plotMat(M = funByBlocks(res_val), main = "Block densities (valued)")
round(funByBlocks(res_val), 2)




##########################################
############ пункт 6 #####################
##########################################


# install.packages(c("igraph", "Matrix", "RSSL", "caret", "dplyr"))

library(RSSL)
library(igraph)
library(Matrix)
library(caret)       # для confusionMatrix и метрик
library(dplyr)



# Категории
classes <- c("liberal", "neutral", "conservative")
y <- factor(V(G)$political_ideology, levels = classes)

n <- vcount(G)
A <- igraph::as_adjacency_matrix(G, sparse = TRUE)  # разреженная матрица смежности

# Проверим баланс классов
table(y)

# создаем функцию вручную

label_propagation <- function(A, y_train, alpha = 0.9, maxit = 1000, tol = 1e-6) {
  n <- nrow(A)
  classes <- levels(y_train)
  K <- length(classes)
  
  # нормировка: S = D^{-1}A
  deg <- Matrix::rowSums(A)
  Dinv <- Diagonal(n, 1 / pmax(deg, 1))
  S <- Dinv %*% A
  
  # one-hot матрица меток
  Y <- Matrix(0, n, K, sparse = TRUE)
  idx <- which(!is.na(y_train))
  Y[cbind(idx, as.integer(y_train[idx]))] <- 1
  
  # Инициализация распределения F
  F <- Y
  
  # Список помеченных узлов
  L <- which(!is.na(y_train))
  
  for (t in seq_len(maxit)) {
    F_new <- alpha * (S %*% F) + (1 - alpha) * Y
    # клампим помеченные узлы
    F_new[L, ] <- Y[L, ]
    if (max(abs(F_new - F)) < tol) break
    F <- F_new
  }
  
  pred_idx <- max.col(as.matrix(F), ties.method = "first")
  pred <- factor(classes[pred_idx], levels = classes)
  return(pred)
}

# Пропорциональное стратифицированное деление 
train_idx <- createDataPartition(y, p = 0.2, list = FALSE) # оставляем только 20% меток. 
# функция createDataPartition пропорционально убирает метки (то есть в оставшихся данных такая же пропорция классов)
test_idx <- setdiff(seq_len(n), train_idx)

# Размеченные метки только для train
y_train <- y
y_train[test_idx] <- NA  # скрываем метки test-узлов


# финальное обучение с лучшим α на train
final_pred <- label_propagation(A, y_train, alpha = 0.7)

# Оценка на test
cm_final <- confusionMatrix(final_pred[test_idx], y[test_idx])
cm_final

cat("Macro-F1 =", mean(cm_final$byClass[, "F1"], na.rm = TRUE), "\n")





