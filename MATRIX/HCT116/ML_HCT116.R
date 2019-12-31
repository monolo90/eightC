suppressMessages(library(caret))
# leemos matriz
tab <- read.csv("YY1_HCT116_MATRIX.csv", h = T)
# PRUEBA PRIMERO CON 2000 AL AZAR Y LUEGO REPITES CON TODO
#tab <- tab[sample(1:nrow(tab), 2000),]
# dividimos en testing y training (20% vs 80%)
training.idx <- createDataPartition(tab$LOOP, p = .8, list = FALSE)
# conuntos test y train 
training <- tab[training.idx,]
testing  <- tab[-training.idx,]
# parámetros de validación cruzada, fold = 5
fit_params <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
# entrenamos modelo con random forests (rf)
set.seed(123)
rfFit <- train(LOOP ~ ., data = training, 
	method = "rf", 
        trControl = fit_params,
        verbose = FALSE)
# representar el accuracy (sobre el bagging de random forest)
pdf("1_accuracy.pdf");plot(rfFit);dev.off()
# importancias de los predictores del  modelo
importances <- varImp(rfFit)
pdf("2_importances.pdf");plot(importances);dev.off()
# matriz de confusión sobre el training data
conf.mat <- rfFit$finalModel
conf.mat
# aplicamos el modelo al conjunto test
pred <- predict(rfFit, testing)
# matriz de confusión sobre el training data
confusionMatrix(pred, testing$LOOP)
# data frame con las probabilidades de cada loop en test
probs <- predict(rfFit, testing, type="prob")
head(probs)
# representamos curva ROC
suppressMessages(library(pROC))
ROC <- roc(predictor=probs$`YES`,
               response=testing$LOOP,
               levels=rev(levels(testing$LOOP)))
ROC$auc
pdf("3_ROC.pdf");plot(ROC,main="ROC");dev.off()
# selección de atributos
set.seed(763)
# configuramos los parámetros del entrenamiento
control <- rfeControl(functions=rfFuncs, method="cv", number=5)
# entrenamos
rfe_res <- rfe(tab[,1:(ncol(tab) - 1)], tab$LOOP, sizes=1:(ncol(tab) - 1), rfeControl=control)
# imprimimos el resultado
print(rfe_res)
predictors(rfe_res)
# graficamos el accuracy para cada combinación de atributos
pdf("4_accuracy_feature_selection.pdf");plot(rfe_res, type=c("g", "o"));dev.off()
