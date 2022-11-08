## Carregando pacotes
library(arules)
library(stringr)
library(progress)

## Carrega a função RSarules
source("RSarules.R")

## Carrega o dados Groceries
data(Groceries)

## Preparando o Groceries para fazer o exemplo
## Converte transações para matrix
mm <- t(as(Groceries,"ngCMatrix"))
groceries <- as.data.frame(as.matrix(mm))

## Transforma TRUE e FALSE em 0 e 1
# Convert all to numeric
cols <- sapply(groceries, is.logical)
groceries[,cols] <- lapply(groceries[,cols], as.numeric)


## Compilado daqui em diantes
library(compiler)
enableJIT(3)

## Response being the last second item
M <- 500 # number of arules need to be sampled. M = 1000 in the reference paper.
ig <- 10 # the value for the tuning parameter 3, 6, 10
rhs <- 25 # whole milk
result <- RSarules( data = as.matrix(groceries), rhs = rhs, M = M , ig = ig)
result
inspect(result$sampled_rules)

## Colocando os labels corretos
num = as.numeric(str_extract(string = result$sampled_rules@itemInfo$labels, pattern = "\\d+"))
result$sampled_rules@itemInfo$labels = colnames(groceries)[num]
inspect(result$sampled_rules)

## Salvando o resultado como matriz
resultado = t(as(result$sampled_rules,"ngCMatrix"))
resultado <- as.data.frame(as.matrix(resultado))
