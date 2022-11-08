## Regras de associação
library(arules)

## Carregando dados
data(Groceries)  
head(Groceries)
inspect(head(Groceries))

## Obtendo número de itens por transação
size(head(Groceries)) 

## Converte transações para uma lista
LIST(head(Groceries, 3))

## Obtendo regras
rules <- apriori (Groceries, parameter = list(supp = 0.001, conf = 0.5, maxlen=2))

## Ordenando regras por confiança
rules_conf <- sort (rules, by="confidence", decreasing=TRUE)
inspect(rules_conf)

## Converte transações para matrix
mm <- t(as(Groceries,"ngCMatrix"))
groceries <- as.data.frame(as.matrix(mm))

## Transforma TRUE e FALSE em 0 e 1
# Convert all to numeric
cols <- sapply(groceries, is.logical)
groceries[,cols] <- lapply(groceries[,cols], as.numeric)
