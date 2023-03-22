rm(list=ls())
setwd("C:/Users/thais/Documents/PhD/FINISHED/Sneaker project/Boostedregressions/")

#load data
final_LHS <- read.csv("CDMetaPop_finaldata.csv")

#### Data exploration ####

sum(final_LHS$propsneaker < 1)#Number of simulation where sneaker where there are territorial males left at the end. 
sum(final_LHS$term1_Tobs < 1)#NUmber of simulation where sneaker were fixed.
mean(final_LHS$TObs_first_AA[final_LHS$term1_Tobs < 1])#mean number of years to sneaker fixation
final_LHS$TObs_first_AA[final_LHS$term1_Tobs < 1] #These are the years where sneaker became fixed.
mean(final_LHS$Sobs_ter[final_LHS$term1_Tobs==1])#mean number of territorial males in fixed sim. 
final_LHS$Sobs_ter[final_LHS$term1_Tobs==1]#There are the  number of territorial males left at the end in sim that where not fixed. 

hist(final_LHS$Sneakmetric,breaks=20) 
mean(final_LHS$Sneakmetric)

#distribution of the explanatory variables (for supmat)
test 
pdf("explanatoryvar.pdf")
par(mfrow=c(2,2))
dotchart(final_LHS$AssortativeMate_Factor,main="Assortative Mating factor",ylab="LHS #")
dotchart(final_LHS$mature_eqn_slope,main="Sneaker maturation slope")
dotchart(final_LHS$N0,main="Sneaker (%) at initiation",xlab="value",ylab="LHS #")
dotchart(final_LHS$Prop_sizemort,main="Sneaker mortality",xlab="value")
dev.off()

#### GLM don't do a great job at modelling the data ####

#try with binomial first.
model1 <- glm(final_LHS$Sneakmetric~AssortativeMate_Factor+mature_eqn_slope+N0+Prop_sizemort,
              family = binomial(link="logit"),
              data=final_LHS)
car::vif(model1) # as expected there is no problem of multicollinearity in model1.
summary(model1) 
#model selection
drop1(model1, test = "Chi") #assortative mate factor is significant based on LRT...
#explanatory power
pseudoR2 <- (model1$null.deviance - model1$deviance) / model1$null.deviance
pseudoR2 # this model explain ~58% of the variation
#check assumption of the model
library(DHARMa)
simResids1 <- simulateResiduals(model1)
plot(simResids1)#significant deviation
# KS test for correct distribution of residuals
testUniformity(simResids1)#not unif disribued
testDispersion(simResids1)
testOutliers(simResids1, type = "bootstrap")
testZeroInflation(simResids1)
#unfortunately this doesnt look so good
# the model is violating teh assumptions of residual distribution and equality of variance. 


#lets try with quasibinomial.  
model2 <- glm(final_LHS$Sneakmetric~AssortativeMate_Factor+mature_eqn_slope+N0+Prop_sizemort,
              family = quasibinomial,
              data=final_LHS)
summary(model2) 
drop1(model2, test = "Chi") 
#this says all but mature eqn are significant. N0, followed by Assortative mate, explain most of the variance. 
#Explanatory power is good
pseudoR2 <- (model2$null.deviance - model2$deviance) / model2$null.deviance
pseudoR2 # this model explain ~58% of the variation
#hard to check the assumptions of quasiB. 
scatter.smooth(fitted(model2), resid(model2))
abline(h=0, lty=2)


#lets try with observation-level random effect in binomial.  
final_LHS$obs<-seq(nrow(final_LHS)) #create observation-level random
library(lme4)
model1_OLRE <- glmer(final_LHS$Sneakmetric~AssortativeMate_Factor+mature_eqn_slope+N0+Prop_sizemort+(1|obs),
                   family = binomial(link="logit"),
                   data=final_LHS) 
simResids1OLRE <- simulateResiduals(model1_OLRE)
plot(simResids1OLRE) 
#violates assumption of residual distribution, equality of variance, and dispersion....

model3 <- glm(final_LHS$Sneakmetric~AssortativeMate_Factor+mature_eqn_slope+N0+Prop_sizemort,
              family = poisson,
              data=final_LHS)
simresid2 <- simulateResiduals(model3)
plot(simresid2)#significant deviation



#Look at plots of relationships for model 2.
pdf("GLMquasibin.pdf")
par(mfrow=c(2,2))
plot(SysFin_sneakProp_Final~N0,data=final_LHS,pch=21,bg="grey",cex=1.5,xlab="N0",ylab="Sneaker (%) in metapop")
plot(SysFin_sneakProp_Final~AssortativeMate_Factor,data=final_LHS,pch=21,bg="grey",cex=1.5,xlab="Assortative mating factor",ylab="")
plot(SysFin_sneakProp_Final~Prop_sizemort,data=final_LHS,pch=21,bg="grey",cex=1.5,xlab="Proportion size mortality",ylab="Sneaker (%) in metapop")
plot(SysFin_sneakProp_Final~mature_eqn_slope,data=final_LHS,pch=21,cex=1.5,xlab="mature eqn slope",ylab="")
dev.off()
# Clearly, those models do not do a good job at describing the data (we cannot rpedict from them)

#### More data exploration for usefulness of interactions####
library("plot3D")
scatter3D(final_LHS$AssortativeMate_Factor,final_LHS$N0,final_LHS$Sneakmetric,pch=16,
                xlab="AssortativeMate_Factor",ylab="N0",zlab="Relative sneaker success")
scatter3D(final_LHS$AssortativeMate_Factor,final_LHS$Prop_sizemort,final_LHS$Sneakmetric,pch=16,
          xlab="AssortativeMate_Factor",ylab="Prop_sizemort",zlab="Relative sneaker success")


library(plotly)
p <- plot_ly(x=final_LHS$AssortativeMate_Factor,y=final_LHS$N0,z=final_LHS$Sneakmetric,color=final_LHS$Sneakmetric)
#create axis titles (unecessarily complicated with 3d plots in plotly...)
axx <- list(title = "AssortativeMate_Factor")
axy <- list(title = "N0")
axz <- list(title = "SysFin_sneakProp_Final")
#plot
p %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))


#### Fit a single prunned regression tree####
#The perks: regression trees are nonparametric
#works by partitioning the data into smaller groups (nodes) of similar response values, and a constant within each subgroup (the mean)
#Root node= the one at the top; terminal nodes=leaves; branches connect the nodes. 

#Useful refs
#https://bradleyboehmke.github.io/HOML/DT.html
#https://uc-r.github.io/regression_trees

library(rpart)  
#create data frame with only variables of interest
df <- data.frame(final_LHS$Sneakmetric,final_LHS$AssortativeMate_Factor,final_LHS$mature_eqn_slope,final_LHS$N0,final_LHS$Prop_sizemort)

# Create training (70%) and test (30%) sets for the AmesHousing::make_ames() data.
# Use set.seed for reproducibility
set.seed(123)
library(rsample)
ames_split <- initial_split(df, prop = .7)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

#Produce plot for supmat (this time with different colours for training vs testing set)
ames_train2 <- ames_train
ames_train2$set <- rep("training",70) 
ames_test2 <- ames_test
ames_test2$set <- rep("testing",30)  
df2 <- rbind(ames_train2,ames_test2)
df2 <- df2[ order(as.numeric(row.names(df2))), ] #reorder by rownames
my_cols <- c("#E69F00","#999999")
pdf("explanatoryvar.pdf")
par(mfrow=c(2,2))
dotchart(df2$final_LHS.AssortativeMate_Factor,main="Assortative Mating factor",ylab="LHS #",groups=as.factor(df2$set),gcolor=my_cols,color=my_cols[as.factor(df2$set)],pch=19)
dotchart(df2$final_LHS.mature_eqn_slope,main="Sneaker maturation slope",groups=as.factor(df2$set),gcolor=my_cols,color=my_cols[as.factor(df2$set)],pch=19)
dotchart(df2$final_LHS.N0,main="Sneaker (%) at initiation",xlab="value",ylab="LHS #",groups=as.factor(df2$set),gcolor=my_cols,color=my_cols[as.factor(df2$set)],pch=19)
dotchart(df2$final_LHS.Prop_sizemort,main="Sneaker mortality",xlab="value",groups=as.factor(df2$set),gcolor=my_cols,color=my_cols[as.factor(df2$set)],pch=19)
dev.off()


#fit basic tree
m1 <- rpart(
  formula = final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
  data    = ames_train,
  method  = "anova"
)
m1#based on 70 observations, we can see that the first split occurs on assortative mate factor. 
#The number of simulatios with c<2.04 is listed (35), with average sneaker propo (0.53) and SSE (0.58)  
library(rpart.plot)
rpart.plot(m1,yesno=2,type=5,extra=100,under=TRUE) 
#nodes show the variables with values on the branches 
#terminal leaves show the predicted value, and under it the percentage of observation in the node. 
#Behind the scenes rpart is automatically applying a range of cost complexity values to prune the tree. 
#To compare the error for each value, rpart performs a 10-fold cross validation so that the error associated with a given 
# value is computed on the hold-out validation data
plotcp(m1)

#The important decisions relate to the complexity of the tree - tradeoff of overly complex tree is overfitting/lack of generalizations. 
#the main two paraneters affectibg it are the number of observations/node (minsplit), and the depth (maxdepth). 
#By default, minsplit=20, maxdepth=30
#lets see what happens when we generate a fully grown tree (no pruning). 
m2 <- rpart(
  formula = final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
  data    = ames_train,
  method  = "anova", 
  control = list(cp = 0, xval = 10)
)
rpart.plot(m2,yesno=2,type=5,extra=100,under=TRUE) 
plotcp(m2) #this shows that we do not gain a lot, compared to m1, in terms of error reduction. 
abline(v = 2, lty = "dashed")

#lets try to change the minimum splitting size and the max depth of the tree. 
#create a grid with a combination of parameter to automate the search
hyper_grid <- expand.grid(
  minsplit = seq(5, 20, 1),
  maxdepth = seq(2, 15, 1)
)
head(hyper_grid)
#iterate through each combination and save models
models <- list()
for (i in 1:nrow(hyper_grid)) {
    # get minsplit, maxdepth values at row i
  minsplit <- hyper_grid$minsplit[i]
  maxdepth <- hyper_grid$maxdepth[i]
  # train a model and store in the list
  models[[i]] <- rpart(
    formula = final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
    data    = ames_train,
    method  = "anova",
    control = list(minsplit = minsplit, maxdepth = maxdepth)
  )
}
# function to get optimal cost complexity
get_cp <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  cp <- x$cptable[min, "CP"] 
}
# function to get minimum error
get_min_error <- function(x) {
  min    <- which.min(x$cptable[, "xerror"])
  xerror <- x$cptable[min, "xerror"] 
}
library(dplyr) 
#identify the best model 
hyper_grid %>%
  mutate(
    cp    = purrr::map_dbl(models, get_cp),
    error = purrr::map_dbl(models, get_min_error)
  ) %>%
  arrange(error) %>%
  top_n(-5, wt = error)
#suggests that the best model has a minsplit of 20 and a max depth of 10 (xerror of 0.45)
#run and plot this optimal model
optimal_tree <- rpart(
  formula = final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
  data    = ames_train,
  method  = "anova",
  control = list(minsplit = 20, maxdepth = 10, cp = 0.0527)
)
#plot
pdf("optimal_tree_sneakerprop.pdf")
rpart.plot(optimal_tree,type=5,extra=100,under=TRUE) 
dev.off()
#see how well it predicts our data
pred <- predict(optimal_tree, newdata = ames_test)
library(Metrics)
rmse_mod3 <- Metrics::rmse(actual=ames_test$final_LHS.Sneakmetric, #Actual values
                  predicted = pred ) # This suggests that, on average, our predicted values are 16% off. 
#here we end up with a single prunned tree. 
#lets see if we can improve predictions using 100 bagged unpruned trees!

#### Use bagging to improve predictive performance#### 
#bagging (bootstrap aggregating) is a general method for fitting multiple versions of a prediction model and 
# then combining them into an aggregated prediction. 
# make bootstrapping reproducible
set.seed(123)
# train bagged model
library(ipred)  
bagged_m1 <- bagging(
  formula = final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
  data    = ames_train,
  nbagg=100, #100 bootstrao rep
  coob    = TRUE,
  control = rpart.control(minsplit = 2, cp = 0) #specify no prunning and a minsplit of 2
)
bagged_m1 #out of bag error = 0.135

# assess 10-100 bagged trees
ntree <- 10:100
rmse <- vector(mode = "numeric", length = length(ntree))
for (i in seq_along(ntree)) {
  # reproducibility
  set.seed(123)
  
  # perform bagged model
  model <- bagging(
    formula = final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
    data    = ames_train,
    coob    = TRUE,
    nbagg   = ntree[i]
  )
  # get OOB error
  rmse[i] <- model$err
}
plot(ntree, rmse, type = 'l', lwd = 2) #you can tell that the predictive power doesnt change much with increasing number of trees(rmse varies by ~1%). 

library(vip)
vip(bagged_m1, num_features = 4, bar = FALSE)

library(caret)
# Specify 10-fold cross validation
ctrl <- trainControl(method = "cv",  number = 10) 

# CV bagged model
bagged_cv <- train(
  final_LHS.Sneakmetric~final_LHS.AssortativeMate_Factor+final_LHS.mature_eqn_slope+final_LHS.N0+final_LHS.Prop_sizemort,
  data = ames_train,
  method = "treebag",
  trControl = ctrl,
  nbagg = 200,
  importance = TRUE
)

# assess results
bagged_cv
p0 <- vip::vip(bagged_cv, num_features = 4,geom="point") #from this figure we can tell which variable influence the output the most. 
pdf("baggedtree_sneakerprop_variables.pdf")
p0
dev.off()
# Construct partial dependence plots
p1 <- pdp::partial(
  bagged_cv, 
  pred.var = "final_LHS.AssortativeMate_Factor",
  grid.resolution = 20
) %>% 
  autoplot()

p2 <- pdp::partial(
  bagged_cv, 
  pred.var = "final_LHS.N0", 
  grid.resolution = 20
) %>% 
  autoplot()

p3 <- pdp::partial(
  bagged_cv, 
  pred.var = "final_LHS.Prop_sizemort", 
  grid.resolution = 20
) %>% 
  autoplot()
pdf("baggedtree_sneakerprop_variableshape.pdf")
gridExtra::grid.arrange(p1, p2,p3, nrow = 2, ncol=2)
dev.off()
  
pred <- predict(bagged_cv, ames_test)
rmse_bagged_cv <- Metrics::rmse(actual=ames_test$final_LHS.Sneakmetric, #Actual values
                  predicted = pred ) 
rmse_bagged_cv
#RMSE is ~15.7%

#https://uc-r.github.io/regression_trees

#### Use boosting to see if performance can be furter improved####

#these models are appealing because they can identify and model the most relevant interactions...Plus, boosting reduces biases.
library(dismo)

#create data frame with only variables of interest
df <- data.frame(final_LHS$Sneakmetric,final_LHS$AssortativeMate_Factor,final_LHS$mature_eqn_slope,final_LHS$N0,final_LHS$Prop_sizemort)

# Create training (70%) and test (30%) sets.
# Use set.seed for reproducibility
set.seed(123)
library(rsample)
ames_split <- initial_split(df, prop = .8)
ames_train <- training(ames_split)
ames_test  <- testing(ames_split)

#run simple model. 
ex1 <- gbm.step(data=df, gbm.x = 2:5, #predictors
         gbm.y = 1, #response
         family = "gaussian",#my data countains no zero, it is not really a binary outcome. 
         tree.complexity = 3, #tree complexity; best to stick with 2 or 3 since the data set is relatively small (<500)
         learning.rate = 0.01,#tradeoff with tc; try to get a combination that gives >1,000 trees.  
         bag.fraction = 0.7)

#To identify the optimal parameter combination, I am going to run models for a range of parameter combinations. 

#create a grid of the parameters I would like to explore
hyper_grid <- expand.grid(
  lr = c(0.01,0.005,0.001,0.0005),#learning rate: cannot fit model for lr>0.01 with my data.  
  tc = seq(2, 4, 1), #tree complexity (number of nodes): allowing 2 to 4-way interactions.
  bf=seq(0.5,0.70,0.10) #random sample of the data selected
)

###Run for a combination of parameter (above) and store models in a list. 
models <- list()
for (i in 1:nrow(hyper_grid)) {
  # get minsplit, maxdepth values at row i
  treecomp <- hyper_grid$tc[i]
  learnrat <- hyper_grid$lr[i]
  bagfract <- hyper_grid$bf[i]
    # train a model and store in the list
    models[[i]] <- gbm.step(data=ames_train, gbm.x = 2:5, #predictors
           gbm.y = 1, #response
           family = "gaussian",#my data countains no zero, it is not really a binary outcome. 
           tree.complexity = treecomp, #tree complexity; best to stick with 2 or 3 since the data set is relatively small (<500)
           learning.rate = learnrat,#tradeoff with tc; try to get a combination that gives >1,000 trees.  
           bag.fraction = bagfract)
}

#Function to get metrics
#get_metrics <- function(x){
#pd <- x$cv.statistics$deviance.mean #predictive deviance
#nt <- tail(x$trees.fitted,1)
#metrics <- cbind(pd,nt)
#return(metrics)#number of trees
#}

get_MSE <- function(x){
pred = predict(x,ames_test)
MAE <- mean(abs(ames_test$final_LHS.Sneakmetric -pred)) #mean abso error
nt <- tail(x$trees.fitted,1) #number of trees
metrics <- cbind(MAE,nt)
return(metrics)
} 

out <- t(as.data.frame(sapply(models,get_MSE)))
colnames(out) <- c("MAE","Ntrees")
final <- cbind(hyper_grid,out)
format(final[order(final$MAE),],scientific=FALSE)
write.csv(format(final[order(final$MAE),],scientific=FALSE),"paramsel.csv") #save it
#When ran on the whole dataset, the best model has a lr=0.005, tc=2, and bf=0.5. NTrees=1,400. 

angaus.fin <- gbm.step(data=ames_train, gbm.x = 2:5, #predictors
                            gbm.y = 1, #response
                            family = "gaussian",#my data countains no zero, it is not really a binary outcome. 
                            tree.complexity = 2, #tree complexity; best to stick with 2 or 3 since the data set is relatively small (<500)
                            learning.rate = 0.005,#tradeoff with tc; try to get a combination that gives >1,000 trees.  
                            bag.fraction = 0.5) #proportion of obs used in selecting variables. influences stochasticity. This is the default. 
angaus.fin$self.statistics
Rsq <- (angaus.fin$self.statistics$mean.null-angaus.fin$self.statistics$mean.resid)/angaus.fin$self.statistics$mean.null
Rsq #(mean tot deviance - mean resid dev)/mean tot deviance

saveRDS(angaus.fin, file = "Fin_GBM.rds")

summary(angaus.fin) # look at the influence of the predictors. 
#measures are based on the number of times a variable is selected for splitting, weighted by the squared improvement to the model as a result of each split, and aver- aged over all trees (Friedman & Meulman 2003). The relative influence (or contribution) of each variable is scaled so that the sum adds to 100, with higher numbers indicating stronger influence on the response

#lets simplify the model, which is generally recommended for small datasets. 
angaus.fin_simp <- gbm.simplify(angaus.fin)

#vertical red line and model output suggest dropping one predictor improves model perf slightly. 
angaus.tc2.lr01.simp_sub <- gbm.step(data=ames_train,
                                  gbm.x=angaus.fin_simp$pred.list[[1]], #This line is what drops the extra variable. 
                                 gbm.y=1,
                                  family="gaussian", tree.complexity=2, learning.rate=0.01)

saveRDS(angaus.tc2.lr01.simp_sub, file = "Fin_GBM_sim.rds")
angaus.tc2.lr01.simp_sub <- readRDS("Fin_GBM_sim.rds")
#now lets interpret this model!
var <- as.data.frame(summary(angaus.tc2.lr01.simp_sub)) #effect of the variables on the response
var$var <- c("Assortative Mating Factor","Sneaker (%) at initialization","Sneaker mortality")

pdf("Variables.pdf")
barplot(var$rel.inf,main="Relative influence of predictor variables",cex.main=1.3,
        xaxt="n", yaxt="n")
axis(2, tck=1, col.ticks="light gray")
barplot(var$rel.inf,main="Relative influence of predictor variables",
        names.arg=c("Assortative Mate Factor","Sneaker at initialization","Sneaker mortality"),
        add=TRUE)
dev.off()        

#partial dependence plots; shows the effect of a variable after accountinf for the effect of the other
pdf("BRT_variableshape_smoothed.pdf")
gbm.plot(angaus.tc2.lr01.simp_sub, n.plots=4, plot.layout=c(2, 2), write.title = FALSE)
dev.off()

pdf("BRT_variableshape_notsmoothed.pdf")
gbm.plot.fits(angaus.tc2.lr001.simp)
dev.off()

#
library(ggBRT)
list_pred <- plot.gbm.4list(angaus.tc2.lr01.simp_sub) #produces a list of predictors. 
testboot<-gbm.bootstrap.functions(angaus.tc2.lr01.simp_sub,list_pred,n.reps=1000)#make 1,000 bootstrap samples
saveRDS(testboot, file = "testboot.rds")

testboot <- readRDS("testboot.rds")
B1 <- ggPD_boot(angaus.tc2.lr01.simp_sub,predictor=1,list.4.preds=list_pred,booted.preds=testboot$function.preds,cex.line=1,col.line="black",col.ci="red")
B2 <- ggPD_boot(angaus.tc2.lr01.simp_sub,predictor=2,list.4.preds=list_pred,booted.preds=testboot$function.preds,cex.line=1,col.line="black",col.ci="red")
B3 <- ggPD_boot(angaus.tc2.lr01.simp_sub,predictor=3,list.4.preds=list_pred,booted.preds=testboot$function.preds,cex.line=1,col.line="black",col.ci="red")

B1b <- B1+xlab("Assortative Mate Factor (71.6%)")+theme(panel.grid.major=element_line(colour="grey"))
B2b <- B2+xlab("Sneaker (%) at initialization (20.6%)")+ylim(-0.17,0.1)
B3b <- B3+xlab("Sneaker mortality (7.8%)")+ylim(-0.1,0.1)

pdf("BRTpartialplots2.pdf")
gridExtra::grid.arrange(B1b,B2b,B3b)
dev.off()
#Evaluate model
pred = predict(angaus.tc2.lr01.simp_sub,ames_test)
cor(pred,ames_test$final_LHS.Sneakmetric,method="spearman") #correlation is 0.79
#Mean absolute error
mean(abs(ames_test$final_LHS.Sneakmetric -pred)) #on avergae, values are 9.7% off
