##############################################################
## Simulation - differences between having random
##              intercepts or slopes
##              does it inflates 1st or 2nd type errors?
## code by Scandola M., idea by Tidoni E.
##############################################################

library(mvtnorm)
library(clusterGeneration)
library(lmerTest)
library(car)

## is it better y ~ Group * NWCond + (NWCond | ID)
##           or y ~ Group * NWCond + (1 | ID/NWCond) ?

## NGroups = number of simulated groups (between-subjects)
## NTrials = number of trials each condition
## NWCond  = number of levels for within-subjects conds
## NSubj   = total number of subjects
NGroups <- 3
NWCond  <- 3
NSubj   <- 30
NTrials <- 30

## coefficients for fixed effects
betasH0 <- c( 0 , 0 , 0 , 0 , 0 , 0 , 0 ,  0 , 0)
betasH1 <- c( 0 , 0 , 0 , 0 , 0 , 0,0.4 ,0.4, 0)

## output list
output <- list()

for(iteration in 1:2000){
  ##################
  ## random generation of random effects for each participant
  ##################
  stddev <- rep(0.5, 9)
  print( iteration )
  corMat <- matrix(
    c(
      1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
      0.2, 1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
      0.2,0.2,  1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
      0.2,0.2,0.2,   1, 0.2, 0.2, 0.2, 0.2, 0.2,
      0.2,0.2,0.2, 0.2,   1, 0.2, 0.2, 0.2, 0.2,
      0.2,0.2,0.2, 0.2, 0.2,   1, 0.2, 0.2, 0.2,
      0.2,0.2,0.2, 0.2, 0.2, 0.2,   1, 0.2, 0.2,
      0.2,0.2,0.2, 0.2, 0.2, 0.2, 0.2,   1, 0.2,
      0.2,0.2,0.2, 0.2, 0.2, 0.2, 0.2, 0.2,   1
    ),
    ncol=9
  )
  
  covMat <- stddev %*% t(stddev) * corMat
  random.effects <- rmvnorm(NSubj, mean = c( 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ), sigma = covMat)
  
  ##################
  ## adding random effects to fixed effects
  ##################
  coefsH0 <- coefsH1 <- matrix( NA , ncol = ncol(random.effects) , nrow = nrow(random.effects))
  for( i in 1:ncol(random.effects) ){
    coefsH0[ , i ] <- random.effects[ , i] + betasH0[ i ]
    coefsH1[ , i ] <- random.effects[ , i] + betasH1[ i ]
  }
  
  ##################
  ## matrix design
  ##################
  data.sim <- expand.grid(
    trial      = 1:NTrials,
    ID         = factor(1:NSubj),
    WCond      = factor(1:NWCond)
  )
  
  ## group is between-subjects, we must be sure that subjects are not in more than one group
  data.sim$Group <- NA
  
  gg <- rep(1:3,each=NSubj/NGroups)
  
  for( i in 1:length(levels(data.sim$ID))){
    data.sim$Group[data.sim$ID==as.character(i)] <- gg[i]
  }
  data.sim$Group <- factor( data.sim$Group )
  
  contrasts( data.sim$Group ) <- contr.sum( n = 3 )
  contrasts( data.sim$WCond ) <- contr.sum( n = 3 )
  
  ##################
  ## d.v. generation
  ##################
  yH0 <- yH1 <- rep( times = nrow(data.sim) , NA )
  
  for( i in 1:length( levels( data.sim$ID ) ) ){
    
    sel <- which( data.sim$ID == as.character(i) )
    
    mm  <- model.matrix(~ Group * WCond , data = data.sim[ sel , ] )
    
    yH0[sel] <- mm %*% as.matrix(coefsH0[ i , ]) + rnorm( n = NTrials * NWCond )
    
    yH1[sel] <- mm %*% as.matrix(coefsH1[ i , ]) + rnorm( n = NTrials * NWCond )
      
  }
  
  ##################
  ## models fitting
  ##################
  start.time         <- Sys.time()
  slope.H0           <- suppressMessages( lmer( yH0 ~ Group * WCond + (WCond | ID) , data = data.sim ) )
  elapsed.time.sH0   <- Sys.time() - start.time
  
  start.time         <- Sys.time()
  slope.H1 <- suppressMessages( lmer( yH1 ~ Group * WCond + (WCond | ID) , data = data.sim ) )
  elapsed.time.sH1   <- Sys.time() - start.time
  
  start.time         <- Sys.time()
  int.H0 <- suppressMessages( lmer( yH0 ~ Group * WCond + ( 1 | ID/WCond) , data = data.sim ) )
  elapsed.time.iH0   <- Sys.time() - start.time
  
  start.time         <- Sys.time()
  int.H1 <- suppressMessages( lmer( yH1 ~ Group * WCond + ( 1 | ID/WCond) , data = data.sim ) )
  elapsed.time.iH1   <- Sys.time() - start.time
  
  ##################
  ## models analysis
  ##################
  start.time         <- Sys.time()
  AsH1 <- Anova(slope.H1 , type = 3 , test = "F")
  elapsed.time.AsH1   <- Sys.time() - start.time
  
  start.time         <- Sys.time()
  AiH1 <- Anova(int.H1   , type = 3 , test = "F")
  elapsed.time.AiH1   <- Sys.time() - start.time
  
  start.time         <- Sys.time()
  AsH0 <- Anova(slope.H0 , type = 3 , test = "F")
  elapsed.time.AsH0   <- Sys.time() - start.time
  
  start.time         <- Sys.time()
  AiH0 <- Anova(int.H0   , type = 3 , test = "F")
  elapsed.time.AiH0   <- Sys.time() - start.time
  
  output[[iteration]] <- data.frame(
    cbind(
      coefs     = rownames(AsH1),
      iteration = iteration,
      AsH1      = AsH1$`Pr(>F)`,
      AsH0      = AsH0$`Pr(>F)`,
      AiH1      = AiH1$`Pr(>F)`,
      AiH0      = AiH0$`Pr(>F)`,
      elapsed.time.sH0,
      elapsed.time.sH1,
      elapsed.time.iH0,
      elapsed.time.iH1,
      elapsed.time.AsH1,
      elapsed.time.AsH0,
      elapsed.time.AiH1,
      elapsed.time.AiH0
    )
  )
  
  write.csv2(do.call("rbind", output), file = "output-simulation-random-effects.csv")
}

sessionInfo()
#########################################