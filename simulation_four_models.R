##############################################################
## Simulation - differences between having random
##              intercepts or slopes
##              does it inflates 1st or 2nd type errors?
## code by Scandola M., idea by Tidoni E.
##############################################################

library(mvtnorm)
library(clusterGeneration)
library(lme4)
library(car)

## is it better y ~ Group * NWCond + (NWCond | ID)           code: slope / s
##           or y ~ Group * NWCond + (1 | ID/NWCond)         code: intC  / iC  (intercept complete)
##           or y ~ Group * NWCond + (1 | ID:NWCond)         code: intP  / iP  (intercept partial)
##           or y ~ Group * NWCond + (1 | ID)?               code: int0  / i0  (intercept minimal)

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
  stddev <- runif(9, min = 0.3, max = 1)
  print( iteration )
  corMat <- rcorrmatrix(9)
  
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
  elapsed.time.sH0   <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  slope.H1 <- suppressMessages( lmer( yH1 ~ Group * WCond + (WCond | ID) , data = data.sim ) )
  elapsed.time.sH1   <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  intC.H0 <- suppressMessages( lmer( yH0 ~ Group * WCond + ( 1 | ID/WCond) , data = data.sim ) )
  elapsed.time.iCH0  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  intC.H1 <- suppressMessages( lmer( yH1 ~ Group * WCond + ( 1 | ID/WCond) , data = data.sim ) )
  elapsed.time.iCH1  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  intP.H0 <- suppressMessages( lmer( yH0 ~ Group * WCond + ( 1 | ID:WCond) , data = data.sim ) )
  elapsed.time.iPH0  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  intP.H1 <- suppressMessages( lmer( yH1 ~ Group * WCond + ( 1 | ID:WCond) , data = data.sim ) )
  elapsed.time.iPH1  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  int0.H0 <- suppressMessages( lmer( yH0 ~ Group * WCond + ( 1 | ID) , data = data.sim ) )
  elapsed.time.i0H0  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  int0.H1 <- suppressMessages( lmer( yH1 ~ Group * WCond + ( 1 | ID) , data = data.sim ) )
  elapsed.time.i0H1  <- difftime( Sys.time(), start.time , units = "secs" )
  
  ##################
  ## models analysis
  ##################
  start.time         <- Sys.time()
  AsH1 <- Anova(slope.H1 , type = 3 , test = "F")
  elapsed.time.AsH1  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  AsH0 <- Anova(slope.H0 , type = 3 , test = "F")
  elapsed.time.AsH0  <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  AiCH1 <- Anova(intC.H1   , type = 3 , test = "F")
  elapsed.time.AiCH1 <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  AiCH0 <- Anova(intC.H0   , type = 3 , test = "F")
  elapsed.time.AiCH0 <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  AiPH1 <- Anova(intP.H1   , type = 3 , test = "F")
  elapsed.time.AiPH1 <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  AiPH0 <- Anova(intP.H0   , type = 3 , test = "F")
  elapsed.time.AiPH0 <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  Ai0H1 <- Anova(int0.H1   , type = 3 , test = "F")
  elapsed.time.Ai0H1 <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  Ai0H0 <- Anova(int0.H0   , type = 3 , test = "F")
  elapsed.time.Ai0H0 <- difftime( Sys.time(), start.time , units = "secs" )
  
  ##################
  ## models checks
  ##################
  SsH1  <- isSingular(slope.H1)
  SsH0  <- isSingular(slope.H0)
  SiCH1 <- isSingular(intC.H1)
  SiCH0 <- isSingular(intC.H0)
  SiPH1 <- isSingular(intP.H1)
  SiPH0 <- isSingular(intP.H0)
  Si0H1 <- isSingular(int0.H1)
  Si0H0 <- isSingular(int0.H0)
  
  MsH1  <- ifelse(is.null(slope.H1@optinfo$conv$lme4$messages),"",slope.H1@optinfo$conv$lme4$messages)
  MsH0  <- ifelse(is.null(slope.H0@optinfo$conv$lme4$messages),"",slope.H0@optinfo$conv$lme4$messages)
  MiCH1 <- ifelse(is.null(intC.H1@optinfo$conv$lme4$messages),"",intC.H1@optinfo$conv$lme4$messages)
  MiCH0 <- ifelse(is.null(intC.H0@optinfo$conv$lme4$messages),"",intC.H0@optinfo$conv$lme4$messages)
  MiPH1 <- ifelse(is.null(intP.H1@optinfo$conv$lme4$messages),"",intP.H1@optinfo$conv$lme4$messages)
  MiPH0 <- ifelse(is.null(intP.H0@optinfo$conv$lme4$messages),"",intP.H0@optinfo$conv$lme4$messages)
  Mi0H1 <- ifelse(is.null(int0.H1@optinfo$conv$lme4$messages),"",int0.H1@optinfo$conv$lme4$messages)
  Mi0H0 <- ifelse(is.null(int0.H0@optinfo$conv$lme4$messages),"",int0.H0@optinfo$conv$lme4$messages)
  
  ##################
  ## output
  ##################
  output[[iteration]] <- data.frame(
    cbind(
      coefs      = rownames(AsH1),
      iteration  = iteration,
      ## anova tables
      AsH1       = AsH1$`Pr(>F)`,
      AsH0       = AsH0$`Pr(>F)`,
      AiCH1      = AiCH1$`Pr(>F)`,
      AiCH0      = AiCH0$`Pr(>F)`,
      AiPH1      = AiPH1$`Pr(>F)`,
      AiPH0      = AiPH0$`Pr(>F)`,
      Ai0H1      = Ai0H1$`Pr(>F)`,
      Ai0H0      = Ai0H0$`Pr(>F)`,
      ## elapsed times for models
      elapsed.time.sH0,
      elapsed.time.sH1,
      elapsed.time.iCH0,
      elapsed.time.iCH1,
      elapsed.time.iPH0,
      elapsed.time.iPH1,
      elapsed.time.i0H0,
      elapsed.time.i0H1,
      ## elapsed times for anovas
      elapsed.time.AsH1,
      elapsed.time.AsH0,
      elapsed.time.AiCH1,
      elapsed.time.AiCH0,
      elapsed.time.AiPH1,
      elapsed.time.AiPH0,
      elapsed.time.Ai0H1,
      elapsed.time.Ai0H0,
      ## models singularity
      SsH1,
      SsH0,
      SiCH1,
      SiCH0,
      SiPH1,
      SiPH0,
      Si0H1,
      Si0H0,
      ## models convergence messages
      MsH1,
      MsH0,
      MiCH1,
      MiCH0,
      MiPH1,
      MiPH0,
      Mi0H1,
      Mi0H0
    )
  )
  
  write.csv2(do.call("rbind", output), file = "output-simulation-random-effects.csv")
}

sessionInfo()
#########################################