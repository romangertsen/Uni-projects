install.packages("igraph")
library(igraph)
install.packages("igraphdata")
library(igraphdata)
nodes <- read.csv("political-books-nodes.csv", header=T) 
links <- read.csv("political-books-edges.csv", header=T)
g <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
g

vcount
ecount
edge_density(g)
is.connected(g)
mean(degree(g))
diameter(g)
average.path.length(g)
transitivity(g, type = "global")
vcount(g)
ecount(g)
edge_density(g)
plot(g, vertex.label=NA, layout = layout_with_fr)
plot(g, vertex.label=NA, layout = layout_with_kk)
plot(g, vertex.label=NA, layout = layout_nicely)

ideology_colors <- c("conservative" = "red", "liberal" = "blue", "neutral" = "green")
V(g)$color <- ideology_colors[V(g)$political_ideology]
V(g)$size <- 4 + sqrt(degree(g)) * 1.5
V(g)$frame.color <- "white"
V(g)$frame.width <- 0.5
E(g)$width <- 0.8
E(g)$color <- "gray85"
plot(g,
     vertex.label = NA,
     layout = layout_with_fr,
     edge.arrow.size = 0,
     margin = c(0, 0, 0, 0))




# модель Эрдеша-Реньи
p <- edge_density(g)
g.random <- sample_gnp(n = vcount(g), p = p)

# есть ли визуальные отличия?
par(mfrow = c(1, 2))
plot(g, vertex.label = NA, 
     vertex.color = c("conservative" = "red", "liberal" = "blue", "neutral" = "gray")[V(g)$political_ideology],
     main = "Реальная сеть\nполитических книг")
plot(g.random, vertex.label = NA, 
     main = "Случайная сеть\nЭрдёша-Реньи")
par(mfrow = c(1, 1))

log(vcount(g))/vcount(g) # как соотносится p с пороговыми значениями?

# распределение степеней вершины
N <- round(1 + log(vcount(g), 2)) + 1 # правило Стерджеса
par(mfrow = c(1, 2))
hist(degree(g), main = "Распределение степеней\nреальной сети", breaks = N, prob = TRUE, xlab = "Степень")
hist(degree(g.random), main = "Распределение степеней\nслучайной сети", breaks = N, prob = TRUE, xlab = "Степень")
par(mfrow = c(1, 1))

# построение "доверительного интервала" для средней длины кратчайшего пути
gl <- vector('list', 1000)
m <- rep(NA, 1000)
for(i in 1:1000){
  gl[[i]] <- sample_smallworld(dim = 1, size = vcount(karate), nei = 3, p = 0.05)
  m[i] <- diameter(gl[[i]])
}
hist(m, main = "Распределение средней длины пути\nв случайных сетях", xlab = "Средняя длина пути")
abline(v = diameter(g), col = "red", lty = 3, lwd = 2)
legend("topright", legend = "Реальная сеть", col = "red", lty = 3, lwd = 2)

# построение "доверительного интервала" для предпочтительного присоединения
gl <- vector('list', 1000)
m <- rep(NA, 1000)
for(i in 1:1000){
  gl[[i]] <- sample_gnp(n = vcount(g), p = p)
  m[i] <- assortativity_degree(gl[[i]])
}
hist(m, xlim = range(c(-1, 1)), 
     main = "Распределение ассортативности по степени\nв случайных сетях", xlab = "Ассортативность по степени")
abline(v = assortativity_degree(g), col = "red", lty = 3, lwd = 2)
legend("topright", legend = "Реальная сеть", col = "red", lty = 3, lwd = 2)

# "рассыпем" на вершины политическую идеологию случайно в заданной пропорции
gl <- vector('list', 1000)
m <- rep(NA, 1000)
for(i in 1:1000){
  gl[[i]] <- sample_gnp(n = vcount(g), p = p)
  # Случайно присваиваем идеологии в тех же пропорциях
  random_ideology <- sample(V(g)$political_ideology)
  m[i] <- assortativity(gl[[i]], types1 = as.factor(random_ideology))
}
hist(m, xlim = range(c(-1, 1)), 
     main = "Распределение ассортативности по идеологии\nв случайных сетях", xlab = "Ассортативность по идеологии")
abline(v = assortativity(g, types1 = as.factor(V(g)$political_ideology)), 
       col = "red", lty = 3, lwd = 2)
legend("topright", legend = "Реальная сеть", col = "red", lty = 3, lwd = 2)

# модель Watts-Strogatz
g.WS <- sample_smallworld(dim = 1, size = vcount(g), nei = 4, p = 0.05)
par(mfrow = c(1, 2))
plot(g.WS, vertex.label = NA, layout = layout_in_circle, vertex.size = 5,
     main = "Модель малого мира\nWatts-Strogatz (круг)")
plot(g.WS, vertex.label = NA, layout = layout_nicely, vertex.size = 5,
     main = "Модель малого мира\nWatts-Strogatz (авто)")
par(mfrow = c(1, 1))

gl_global_ws <- vector('list', 1000)
m_global_ws <- rep(NA, 1000)

for(i in 1:1000){
  gl_global_ws[[i]] <- sample_smallworld(dim = 1, size = vcount(g), nei = 4, p = 0.05)
  m_global_ws[i] <- transitivity(gl_global_ws[[i]], type = "global")
}


# Доверительный интервал для СРЕДНЕЙ ДЛИНЫ ПУТИ в малом мире
gl_path <- vector('list', 1000)
m_path <- rep(NA, 1000)

for(i in 1:1000){
  gl_path[[i]] <- sample_smallworld(dim = 1, size = vcount(g), nei = 4, p = 0.05)
  m_path[i] <- mean_distance(gl_path[[i]])
}

hist(m_path, xlim = range(c(0, 10)))
abline(v = mean_distance(g), col = "red", lty = 3, lwd = 2)

# Доверительный интервал для АССОРТАТИВНОСТИ ПО СТЕПЕНИ в малом мире
gl_assort <- vector('list', 1000)
m_assort <- rep(NA, 1000)

for(i in 1:1000){
  gl_assort[[i]] <- sample_smallworld(dim = 1, size = vcount(g), nei = 4, p = 0.05)
  m_assort[i] <- assortativity_degree(gl_assort[[i]])
}

hist(m_assort, xlim = range(c(-1, 1)))
abline(v = assortativity_degree(g), col = "red", lty = 3, lwd = 2)


# можно убедиться, что в реальной сети локальная плотность отличается от малого мира
gl <- vector('list', 1000)
m <- rep(NA, 1000)
for(i in 1:1000){
  gl[[i]] <- sample_smallworld(dim = 1, size = vcount(g), nei = 4, p = 0.05)
  m[i] <- transitivity(gl[[i]], type = "global")
}
hist(m, main = "Распределение кластерного коэффициента\nв сетях малого мира", xlab = "Кластерный коэффициент")
abline(v = transitivity(g, type = "global"), col = "red", lty = 3, lwd = 2)
legend("topright", legend = "Реальная сеть", col = "red", lty = 3, lwd = 2)
transitivity(g, type = "global")
# модель Барабаши-Альберта
g.BA <- sample_pa(n = vcount(g), power = 1, directed = FALSE)
par(mfrow = c(1, 2))
hist(degree(g), main = "Распределение степеней\nреальной сети", breaks = N, prob = TRUE, xlab = "Степень")
hist(degree(g.BA), main = "Распределение степеней\nмодели Барабаши-Альберта", breaks = N, prob = TRUE, xlab = "Степень")
par(mfrow = c(1, 1))

# у политических книг 105 вершин, посмотрим на степенной закон
plot(degree.distribution(g), main = "Распределение степеней\nреальной сети", 
     xlab = "Степень", ylab = "Вероятность")
plot(degree.distribution(g), log = 'xy', main = "Распределение степеней\nв логарифмическом масштабе",
     xlab = "Степень (log)", ylab = "Вероятность (log)")

# сравним с моделью Барабаши-Альберта для большего графа
g.BA.large <- sample_pa(500, power = 1, directed = FALSE)
par(mfrow = c(1, 2))
plot(degree.distribution(g.BA.large), main = "Распределение степеней\nмодели БА (500 вершин)")
plot(degree.distribution(g.BA.large), log = 'xy', main = "Распределение степеней\nмодели БА (логарифмы)")
par(mfrow = c(1, 1))




p <- edge_density(g)

# Доверительный интервал для LOCAL transitivity
gl_local <- vector('list', 1000)
m_local <- rep(NA, 1000)

for(i in 1:1000){
  gl_local[[i]] <- sample_gnp(n = vcount(g), p = p)
  local_trans <- transitivity(gl_local[[i]], type = "local")
  m_local[i] <- mean(local_trans, na.rm = TRUE)
}

hist(m_local, xlim = range(c(0, 1)), 
     main = "Доверительный интервал для LOCAL транзитивности\n(среднее локальных коэффициентов)",
     xlab = "Mean Local Transitivity")
real_local <- mean(transitivity(g, type = "local"), na.rm = TRUE)
abline(v = real_local, col = "red", lty = 3, lwd = 2)
legend("topright", legend = "Реальная сеть", col = "red", lty = 3, lwd = 2)

# Доверительный интервал для GLOBAL transitivity
gl_global <- vector('list', 1000)
m_global <- rep(NA, 1000)

for(i in 1:1000){
  gl_global[[i]] <- sample_gnp(n = vcount(g), p = p)
  m_global[i] <- transitivity(gl_global[[i]], type = "global")
}

hist(m_global, xlim = range(c(0, 1)), 
     main = "Доверительный интервал для GLOBAL транзитивности\n(кластерный коэффициент графа)",
     xlab = "Global Transitivity")
real_global <- transitivity(g, type = "global")
abline(v = real_global, col = "red", lty = 3, lwd = 2)
legend("topright", legend = "Реальная сеть", col = "red", lty = 3, lwd = 2)



# 1. Расчет всех мер центральности
centrality_measures <- data.frame(
  Book = V(g)$Label,
  Ideology = V(g)$political_ideology,
  
  # Центральность по степени (Degree Centrality)
  Degree = degree(g),
  
  # Центральность по близости (Closeness Centrality)
  Closeness = closeness(g),
  
  # Центральность по посредничеству (Betweenness Centrality)
  Betweenness = betweenness(g),
  
  # Центральность по собственному вектору (Eigenvector Centrality)
  Eigenvector = eigen_centrality(g)$vector,
  
  # PageRank
  PageRank = page_rank(g)$vector
)

# 2. Нормализация мер центральности (кроме Betweenness, которая уже нормализована)
centrality_measures$Degree_norm <- centrality_measures$Degree / max(centrality_measures$Degree)
centrality_measures$Closeness_norm <- centrality_measures$Closeness / max(centrality_measures$Closeness)
centrality_measures$Eigenvector_norm <- centrality_measures$Eigenvector / max(centrality_measures$Eigenvector)
centrality_measures$PageRank_norm <- centrality_measures$PageRank / max(centrality_measures$PageRank)

# 3. Топ-10 книг по каждой мере центральности

cat("\n=== ЦЕНТРАЛЬНОСТЬ ПО СТЕПЕНИ (популярность) ===\n")
top_degree <- centrality_measures[order(-centrality_measures$Degree), ][1:10, c("Book", "Ideology", "Degree")]
print(top_degree)

cat("\n=== ЦЕНТРАЛЬНОСТЬ ПО БЛИЗОСТИ (доступность) ===\n")
top_closeness <- centrality_measures[order(-centrality_measures$Closeness), ][1:10, c("Book", "Ideology", "Closeness")]
print(top_closeness)

cat("\n=== ЦЕНТРАЛЬНОСТЬ ПО ПОСРЕДНИЧЕСТВУ (мосты) ===\n")
top_betweenness <- centrality_measures[order(-centrality_measures$Betweenness), ][1:10, c("Book", "Ideology", "Betweenness")]
print(top_betweenness)

cat("\n=== СОБСТВЕННЫЙ ВЕКТОР (влияние) ===\n")
top_eigenvector <- centrality_measures[order(-centrality_measures$Eigenvector), ][1:10, c("Book", "Ideology", "Eigenvector")]
print(top_eigenvector)

cat("\n=== PAGERANK (важность) ===\n")
top_pagerank <- centrality_measures[order(-centrality_measures$PageRank), ][1:10, c("Book", "Ideology", "PageRank")]
print(top_pagerank)

library(dplyr)
# 4. Сравнительный анализ по идеологиям
cat("\n=== СРАВНЕНИЕ ПО ИДЕОЛОГИЯМ ===\n")
ideology_analysis <- centrality_measures %>%
  group_by(Ideology) %>%
  summarise(
    Count = n(),
    Avg_Degree = mean(Degree),
    Avg_Betweenness = mean(Betweenness),
    Avg_Closeness = mean(Closeness),
    Avg_Eigenvector = mean(Eigenvector),
    Avg_PageRank = mean(PageRank)
  )
print(ideology_analysis)

# 5. Визуализация распределения центральностей по идеологиям
library(ggplot2)

# График центральности по степени
ggplot(centrality_measures, aes(x = Ideology, y = Degree, fill = Ideology)) +
  geom_boxplot() +
  labs(title = "Распределение центральности по степени по идеологиям",
       y = "Degree Centrality") +
  scale_fill_manual(values = c("conservative" = "red", "liberal" = "blue", "neutral" = "gray"))

# График центральности по посредничеству
ggplot(centrality_measures, aes(x = Ideology, y = Betweenness, fill = Ideology)) +
  geom_boxplot() +
  labs(title = "Распределение центральности по посредничеству по идеологиям",
       y = "Betweenness Centrality") +
  scale_fill_manual(values = c("conservative" = "red", "liberal" = "blue", "neutral" = "gray"))

# 6. Корреляция между мерами центральности
correlation_matrix <- cor(centrality_measures[, c("Degree", "Closeness", "Betweenness", "Eigenvector", "PageRank")])
cat("\n=== КОРРЕЛЯЦИЯ МЕЖДУ МЕРАМИ ЦЕНТРАЛЬНОСТИ ===\n")
print(round(correlation_matrix, 3))

# 7. Анализ "самых влиятельных" книг
# Книги, входящие в топ-10 по всем мерам центральности
top_books_all_measures <- centrality_measures %>%
  mutate(
    Rank_Degree = rank(-Degree),
    Rank_Betweenness = rank(-Betweenness),
    Rank_Eigenvector = rank(-Eigenvector),
    Total_Rank = Rank_Degree + Rank_Betweenness + Rank_Eigenvector
  ) %>%
  arrange(Total_Rank) %>%
  head(10) %>%
  select(Book, Ideology, Degree, Betweenness, Eigenvector, Total_Rank)