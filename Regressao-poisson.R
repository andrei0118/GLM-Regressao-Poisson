# Base de dados:

#Dados de interesse: número de prêmios conquistados pelos estudantes de um escola secundária.
# Variaveis:
# - número de prêmios
# - tipo de programa que estudante foi matriculado
# - pontuação no exame final de matemática



# Leitura dos dados_____________________________________________________________

dados <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")

# Recodificação de variáveis explicativa em níves de fatores
dados <- within(dados, {
  prog <- factor(prog, levels = 1:3, 
                 labels = c("General", "Academic","Vocational"))
  id <- factor(id)
})

head(dados)

## Análise Descritiva___________________________________________________________

par(mfrow = c(1,3))

# Gráfico histograma Número de prêmios
hist(dados$num_awards, xlab = "Número de prêmios", ylab = "Frequência",
     col = "lightblue", main = "Análise Descritiva Gráfica ", cex.axis = 0.8, cex.lab = 0.8)
# Gráfico de dispersão Nota de Matemática x Número de prêmios
plot(dados$math, dados$num_awards, xlab = "Nota em matemática", ylab = "Número de prêmios",
     cex.axis = 0.8, cex.lab = 0.8)
# Gráfico boxplot Tipo de Programa x Número de prêmios
boxplot(num_awards ~ prog, data = dados, xlab = " ", ylab = "Número de prêmios", 
        cex.axis = 0.8, cex.lab = 0.8, las = 2)


## Ajuste do modelo linear generalizado com distribuição de poisson_____________

modelo = glm(formula = num_awards ~ prog + math, family = "poisson", data = dados)
summary(modelo)

# OBS: p-valor < 0.05 o número esperado de prêmios de um estudante matriculado no
# programa vocacional não difere significativamente do programa geral(intercepct fator referencia),
# (valor - p não significativo)

## Exponencial do Intervalo de confiança para os parãmetros estimados___________

exp(confint(modelo))

## Intervalo de Predição pelo tipo de programa academic_________________________

dados.academic = subset(dados, prog == "Academic")
head(dados.academic)
# Separação do banco de dados e criação de novo banco
novos.dados = data.frame(math = seq(min(dados.academic$math),
                                    max(dados.academic$math), by = 0.1), prog = "Academic")
novos.dados$fit = predict(modelo, newdata = novos.dados, type = "response")

#valores esperados (fit) para número de premios conforme sua nota em matemática:
head(novos.dados)

#Cáculo dos limites superior e inferior:
# calulo de percentil 2.5% para uma distribuição de poisson com valores esperados
li.inf = qpois(0.025, novos.dados$fit)

# calulo de percentil 97.5% para uma distribuição de poisson com valores esperados
li.sup = qpois(0.975, novos.dados$fit)

# Gráfico
par(mfrow = c(1, 1))

# Gráfico de dispersão Nota da prova de matemática x Número de prêmios
plot(num_awards ~ math, color = prog, data = dados.academic, 
     ylab = "Número de prêmios", xlab = "Nota na Prova de Matemática", main = "Programa Acadêmico")
grid()
lines(fit ~ math, data = novos.dados, lwd = 2)
lines(novos.dados$math, li.inf, lty = 2, col = 2, lwd = 2)
lines(novos.dados$math, li.sup, lty = 2, col = 2, lwd = 2)
# obs: 3 pontos fora do intervalo

# Deviance - Análise do ajuste do modelo________________________________________

deviance = with(modelo, cbind(Deviance = deviance, 
                              "Graus de Liberdade" = df.residual,
                              "P-valor" = pchisq(deviance, df.residual, 
                                                 lower.tail = FALSE)))
deviance
# H0: O modelo está bem ajustado
# H1: O modelo não está bem ajustado
# Obs: Com o p -valor obtido, não reijeitamos a hipotese de que o modelo está bem ajustado

# Análise de Resíduos do modelo_________________________________________________

par(mfrow = c(2, 2))
res <- residuals(modelo, type = "deviance");

# Gráfico de dispersão Valores ajustados x Desvio Residual
plot(modelo$fitted, res, xlab = "Valores ajustados", ylab = "Desvio Residual",
     pch = 19, col = "black", cex.axis = 0.6, cex.lab = 0.6, cex = 0.6, ylim = c(-3, 3))
abline(h = c(-3, 0, 3), col = 'red', lty = 2)

# Gráfico dispersão Ordem x Desvio Residual
plot(res, ylab = "Desvio Residual", xlab = "Ordem", 
     cex = 0.6, pch=19,col="black", cex.axis = 0.6, cex.lab = 0.6, ylim = c(-3,3)); 
abline(h=c(-3,0,3), col = 'red', lty =2)

# Gráfico de Desvio Residual
hist(res, main = " ", ylab = "Frequência", xlab = "Desvio Residual", 
     cex.axis = 0.6, cex.lab = 0.6)
qqnorm(res, cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.6, cex = 0.6)

# Gráfico de envelope___________________________________________________________

# Cálculo da matriz H
retorna.H <- function(modelo){
  X = model.matrix(modelo)
  W = diag(modelo$weights)
  M = solve(t(X)%*%W%*%X)
  H = sqrt(W)%*%X%*%M%*%t(X)%*%sqrt(W)
  
  # Extrai os elementos da diagonal da Matriz H
  h = diag(H)
  return(h)}

h = retorna.H(modelo)
residuos.modelo = resid(modelo, type = "deviance")/sqrt((1 - h))
n <- nrow(dados)
m = 100
matriz.residuos <- matrix(NA, n, m)

for(i in 1:m){
  nresp <- rpois(n, fitted(modelo))
  fit <- glm(nresp ~  prog + math, data = dados, family = poisson(link = "log"))
  h = retorna.H(fit)
  matriz.residuos[,i] = sort(resid(fit, type = "deviance")/sqrt(1-h))
}

# Cálculo dos percentis
intervalos = apply(matriz.residuos, 1, function(x) 
  quantile(x, c(0.025, 0.975)))

# Cálculo da média
med <- apply(matriz.residuos, 1, mean)

# Amplitude de variação dos resíduos
faixa = range(residuos.modelo, intervalos)


# Gráfico

qqnorm(residuos.modelo, xlab = "Percentil da Normal Padrão", ylab = "Desvio residual padronizado", 
       ylim = faixa, pch = 16, cex.lab = 0.7, cex.axis = 0.7, main = " ")

# Adicionando o envelope simulado ao gráfico
par(new = T)
qqnorm(intervalos[1,], axes = F, xlab = " ", ylab = " ", type = "l", 
       ylim = faixa, lty = 1, cex = 0.6, main = " ")
par(new = T)
qqnorm(intervalos[2,], axes = F, xlab = " ", ylab = " ", type = "l",
       ylim = faixa, lty = 1, cex = 0.6, main = " ")
par(new = T)
qqnorm(med, axes = F, xlab = " ", ylab = " ", type = "l",
       ylim = faixa, lty = 2, cex = 0.6, main = " ")
# Obs: Apesar da analise de residuos não terem demostrado graficos satisfatorios,
# o gráfico de envelope nos mostra que sim o modelo está bem ajustado

# Alterando e comparando as funções de ligação__________________________________

y = c(2, 3, 6, 7, 8, 9, 10, 12, 15)
x = c(-1, -1, 0, 0, 0, 0, 1, 1, 1)
dados = data.frame(x = x, y = y)

modelo1 <- glm(y ~  x, data = dados, family = poisson(link = "log"))
modelo2 <- glm(y ~ x, data = dados, family = poisson(link = "identity"))
modelo3 <- glm(y ~  x, data = dados, family = poisson(link = "sqrt"))


deviance.m = c(modelo1$deviance, modelo2$deviance, modelo3$deviance)
aic.m = c(modelo1$aic, modelo2$aic, modelo3$aic)
nome = c("Modelo 1", "Modelo 2", "Modelo 3")
frame = data.frame(Modelo = nome, Deviance = deviance.m, AIC = aic.m)
head(frame)
