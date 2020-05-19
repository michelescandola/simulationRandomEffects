##############################################################
## Simulation - differences between having random
##              intercepts or slopes
##              does it impact convergence and singularity?
## code by Scandola M., idea by Tidoni E.
##############################################################

library(mvtnorm)
library(clusterGeneration)
library(lme4)
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

## parameters for the random effects
stddev <- rep(0.5, 9)

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

## output list
output <- list()

for(iteration in 1:2000){
  ##################
  ## random generation of random effects for each participant
  ##################

  print( iteration )

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
  ## models checks
  ##################
  SsH1 <- isSingular(slope.H1)
  SiH1 <- isSingular(int.H1)
  SsH0 <- isSingular(slope.H0)
  SiH0 <- isSingular(int.H0)
  
  OsH1 <- length(slope.H1@optinfo$warnings)==0
  OiH1 <- length(int.H1@optinfo$warnings)==0
  OsH0 <- length(slope.H0@optinfo$warnings)==0
  OiH0 <- length(int.H0@optinfo$warnings)==0
  
  MsH1 <- ifelse(is.null(slope.H1@optinfo$conv$lme4$messages),"",slope.H1@optinfo$conv$lme4$messages)
  MiH1 <- ifelse(is.null(int.H1@optinfo$conv$lme4$messages),"",int.H1@optinfo$conv$lme4$messages)
  MsH0 <- ifelse(is.null(slope.H0@optinfo$conv$lme4$messages),"",slope.H0@optinfo$conv$lme4$messages)
  MiH0 <- ifelse(is.null(int.H0@optinfo$conv$lme4$messages),"",int.H0@optinfo$conv$lme4$messages)
  
  output[[iteration]] <- data.frame(
    cbind(
      iteration = iteration,
      SsH1      = SsH1,
      SsH0      = SsH0,
      SiH1      = SiH1,
      SiH0      = SiH0,
      OsH1,
      OsH0,
      OiH1,
      OiH0,
      MsH1,
      MiH1,
      MsH0,
      MiH0
    )
  )
  
  write.csv2(do.call("rbind", output), file = "output-simulation-checks.csv")
}

sessionInfo()
#########################################