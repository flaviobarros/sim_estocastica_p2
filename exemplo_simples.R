## Exemplo de código 
## Carregando pacotes
library(arules)
library(arulesViz)

## Rodando exemplo
## generate data using R package 'MultiOrd'
set.seed(200)
library(MultiOrd)

l <- 5
n <- 1000
mp <- rep(0.5, l-1)
bcor <- diag( x=1, nrow = l-1, ncol = l-1 )
bcor[1, l-1] <- 0.8
bcor[l-1, 1] <- 0.8
bcor[2, l-1] <- 0
bcor[l-1, 2] <- 0
bcor[3, l-1] <- 0.2
bcor[l-1, 3] <- 0.2
validation.CorrMat( mp, bcor)
dd <- generate.binary( n, mp, bcor)
data <- cbind(dd, 1- dd[, l-1])
colnames(data) <- c( paste( "I", 1:(l-2), sep = ""), "r1", "r2")

## Salvar como csv
write.csv(x = data, file = "simulado.csv")

## Compilado daqui em diantes
library(compiler)
enableJIT(3)

## Response being the last second item
rhs <- dim(data)[2]-1 # the last second item to be in the rhs
lhs_offset <- c( dim(data)[2]) # column numbers that are not contained in the lhs
M <- 1000 # number of arules need to be sampled. M = 1000 in the reference paper.
ig <- 10 # the value for the tuning parameter 3, 6, 10
result <- RSarules( data = data, rhs = rhs, M = M , ig = ig, lhs_offset = lhs_offset )
regras_salvas1 <- inspect(result$sampled_rules)
cbind(regras_salvas1, result$measures)

## Response being the last second item
rhs2 <- dim(data)[2] # the last second item to be in the rhs
lhs_offset2 <- c( dim(data)[2]-1) # column numbers that are not contained in the lhs
M <- 1000 # number of arules need to be sampled. M = 1000 in the reference paper.
ig <- 10 # the value for the tuning parameter 3, 6, 10
result2 <- RSarules( data = data, rhs = rhs2, M = M , ig = ig, lhs_offset = lhs_offset2 )
cbind(inspect(result2$sampled_rules), result2$measures)

## Comparando com o apriori
library(arules)

# rm.duplicates is to remove duplicates because the Apriori algorithm cannot have duplicates
data_transactions = as(data, "transactions")
inspect(data_transactions)

## Obtendo regras
rules <- apriori(data_transactions, parameter = list(supp = 0.001, conf = 0.4),
                 appearance = list(rhs = c("r2")))

## Ordenando regras por confiança
rules_conf <- sort (rules, by="confidence", decreasing=TRUE)
inspect(rules_conf)

########################## DIAGNÓSTICO ########################################
####### trace plot
plot(result$sampled_p[2000:3000], type='l', ylim = c(0,1.1))

## Usando o coda
library(coda)
p <- mcmc(result$sampled_p[200:3000])

## Avaliando a autocorrelação e plot
autocorr(p)
autocorr.plot(p)

## Avaliando a estatítica Geweke
geweke.diag(p) ## está entre -1,96 e 1,96 então está ok

## Tamanho efeivo da amostra
effectiveSize(p) ## autocorrelações negativas produzem isso

## Raftery
raftery.diag(p)
