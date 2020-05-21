##############################################################
## Simulation - differences between having random
##              intercepts or slopes
##              does it inflates 1st or 2nd type errors?
## code by Scandola M., idea by Tidoni E.
##
##  ALL WITHIN-SUBJECTS FACTORS
##  ADDED RANDOM INTERCEPT FOR THE STIMULUS
##
##############################################################

library(mvtnorm)
library(clusterGeneration)
library(lme4)
library(car)

## is it better y ~ Cond1 * Cond2 + (Cond1 * Cond2 | ID)           code: full.slope / fs
##           or y ~ Cond1 * Cond2 + (1 | ID:Cond1:Cond2)           code: full.int   / fi
##           or y ~ Cond1 * Cond2 + (1 | ID:Cond1)                 code: part.int   / pi
##           or y ~ Cond1 * Cond2 + (Cond2 | ID:Cond1)             code: full.mix   / fm
##           or y ~ Cond1 * Cond2 + (Cond1 | ID)                   code: part.slope / ps
##           or y ~ Cond1 * Cond2 + (1 | ID)                       code: minimal    / mn

## NStimuli= number of stimuli
## NCond1  = number of levels for within-subjects cond1
## NTrials = number of trials each condition
## NCond2  = number of levels for within-subjects cond2
## NSubj   = total number of subjects
NCond1  <- 3
NCond2  <- 3
NSubj   <- 30
NTrials <- 5
NStimuli<- 6

## coefficients for fixed effects
betasH0 <- c( 0 , 0 , 0 , 0 , 0 , 0 , 0 ,  0 , 0)
betasH1 <- c( 0 , 0 , 0 , 0 , 0 , 0,0.4 ,0.4, 0)

## generation of random effects for stimuli
ran.stim <- rnorm( NStimuli )
ran.stim <- rep( ran.stim, each = NTrials)

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
    Cond1      = factor(1:NCond1),
    Cond2      = factor(1:NCond2)
  )
  
  contrasts( data.sim$Cond1 ) <- contr.sum( n = 3 )
  contrasts( data.sim$Cond2 ) <- contr.sum( n = 3 )
  
  ##################
  ## d.v. generation
  ##################
  yH0 <- yH1 <- rep( times = nrow(data.sim) , NA )
  
  for( i in 1:length( levels( data.sim$ID ) ) ){
    
    sel <- which( data.sim$ID == as.character(i) )
    
    mm  <- model.matrix(~ Cond1 * Cond2 , data = data.sim[ sel , ] )
    
    yH0[sel] <- mm %*% as.matrix(coefsH0[ i , ]) + rnorm( n = NTrials * NCond1 * NCond2 )
    
    yH1[sel] <- mm %*% as.matrix(coefsH1[ i , ]) + rnorm( n = NTrials * NCond1 * NCond2 )
      
  }
  
  ##################
  ## models fitting
  ##################
  start.time         <- Sys.time()
  fs.H0              <- suppressMessages( lmer( yH0 ~ Cond1 * Cond2 + (Cond1 * Cond2 | ID) + (1 | Stimulus)  , data = data.sim ) )
  et.fs.H0           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  fs.H1              <- suppressMessages( lmer( yH1 ~ Cond1 * Cond2 + (Cond1 * Cond2 | ID) + (1 | Stimulus)  , data = data.sim ) )
  et.fs.H1           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  ps.H0              <- suppressMessages( lmer( yH0 ~ Cond1 * Cond2 + (Cond1 | ID) + (1 | Stimulus)  , data = data.sim ) )
  et.ps.H0           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  ps.H1              <- suppressMessages( lmer( yH1 ~ Cond1 * Cond2 + (Cond1 | ID) + (1 | Stimulus)  , data = data.sim ) )
  et.ps.H1           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  fi.H0              <- suppressMessages( lmer( yH0 ~ Cond1 * Cond2 + ( 1 | ID:Cond1:Cond2) + (1 | Stimulus)  , data = data.sim ) )
  et.fi.H0           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  fi.H1              <- suppressMessages( lmer( yH1 ~ Cond1 * Cond2 + ( 1 | ID:Cond1:Cond2) + (1 | Stimulus)  , data = data.sim ) )
  et.fi.H1           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  pi.H0              <- suppressMessages( lmer( yH0 ~ Cond1 * Cond2 + ( 1 | ID:Cond1) + (1 | Stimulus)  , data = data.sim ) )
  et.pi.H0           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  pi.H1              <- suppressMessages( lmer( yH1 ~ Cond1 * Cond2 + ( 1 | ID:Cond1) , data = data.sim ) )
  et.pi.H1           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  fm.H0              <- suppressMessages( lmer( yH0 ~ Cond1 * Cond2 + (Cond1 | ID:Cond2) + (1 | Stimulus)  , data = data.sim ) )
  et.fm.H0           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  fm.H1              <- suppressMessages( lmer( yH1 ~ Cond1 * Cond2 + (Cond1 | ID:Cond2) + (1 | Stimulus)  , data = data.sim ) )
  et.fm.H1           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  mn.H0              <- suppressMessages( lmer( yH0 ~ Cond1 * Cond2 + (1 | ID) + (1 | Stimulus)  , data = data.sim ) )
  et.mn.H0           <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  mn.H1              <- suppressMessages( lmer( yH1 ~ Cond1 * Cond2 + (1 | ID) + (1 | Stimulus)  , data = data.sim ) )
  et.mn.H1           <- difftime( Sys.time(), start.time , units = "secs" )
  
  ##################
  ## models analysis
  ##################
  start.time         <- Sys.time()
  A.fs.H0            <- Anova(fs.H0 , type = 3 , test = "F")
  et.A.fs.H0         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.fs.H1            <- Anova(fs.H1 , type = 3 , test = "F")
  et.A.fs.H1         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.ps.H0            <- Anova(ps.H0 , type = 3 , test = "F")
  et.A.ps.H0         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.ps.H1            <- Anova(ps.H1 , type = 3 , test = "F")
  et.A.ps.H1         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.fi.H0            <- Anova(fi.H0 , type = 3 , test = "F")
  et.A.fi.H0         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.fi.H1            <- Anova(fi.H1 , type = 3 , test = "F")
  et.A.fi.H1         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.pi.H0            <- Anova(pi.H0 , type = 3 , test = "F")
  et.A.pi.H0         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.pi.H1            <- Anova(pi.H1 , type = 3 , test = "F")
  et.A.pi.H1         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.fm.H0            <- Anova(fm.H0 , type = 3 , test = "F")
  et.A.fm.H0         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.fm.H1            <- Anova(fm.H1 , type = 3 , test = "F")
  et.A.fm.H1         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.mn.H0            <- Anova(mn.H0 , type = 3 , test = "F")
  et.A.mn.H0         <- difftime( Sys.time(), start.time , units = "secs" )
  
  start.time         <- Sys.time()
  A.mn.H1            <- Anova(mn.H1 , type = 3 , test = "F")
  et.A.mn.H1         <- difftime( Sys.time(), start.time , units = "secs" )
  
  ##################
  ## models checks
  ##################
  S.fs.H0  <- isSingular(fs.H0)
  S.fs.H1  <- isSingular(fs.H1)
  S.ps.H0  <- isSingular(ps.H0)
  S.ps.H1  <- isSingular(ps.H1)
  S.fi.H0  <- isSingular(fi.H0)
  S.fi.H1  <- isSingular(fi.H1)
  S.pi.H0  <- isSingular(pi.H0)
  S.pi.H1  <- isSingular(pi.H1)
  S.fm.H0  <- isSingular(fm.H0)
  S.fm.H1  <- isSingular(fm.H1)
  S.mn.H0  <- isSingular(mn.H0)
  S.mn.H1  <- isSingular(mn.H1)
  
  M.fs.H0  <- ifelse(is.null(fs.H0@optinfo$conv$lme4$messages),"",fs.H0@optinfo$conv$lme4$messages)
  M.fs.H1  <- ifelse(is.null(fs.H1@optinfo$conv$lme4$messages),"",fs.H1@optinfo$conv$lme4$messages)
  M.ps.H0  <- ifelse(is.null(ps.H0@optinfo$conv$lme4$messages),"",ps.H0@optinfo$conv$lme4$messages)
  M.ps.H1  <- ifelse(is.null(ps.H1@optinfo$conv$lme4$messages),"",ps.H1@optinfo$conv$lme4$messages)
  M.fi.H0  <- ifelse(is.null(fi.H0@optinfo$conv$lme4$messages),"",fi.H0@optinfo$conv$lme4$messages)
  M.fi.H1  <- ifelse(is.null(fi.H1@optinfo$conv$lme4$messages),"",fi.H1@optinfo$conv$lme4$messages)
  M.pi.H0  <- ifelse(is.null(pi.H0@optinfo$conv$lme4$messages),"",pi.H0@optinfo$conv$lme4$messages)
  M.pi.H1  <- ifelse(is.null(pi.H1@optinfo$conv$lme4$messages),"",pi.H1@optinfo$conv$lme4$messages)
  M.fm.H0  <- ifelse(is.null(fm.H0@optinfo$conv$lme4$messages),"",fm.H0@optinfo$conv$lme4$messages)
  M.fm.H1  <- ifelse(is.null(fm.H1@optinfo$conv$lme4$messages),"",fm.H1@optinfo$conv$lme4$messages)
  M.mn.H0  <- ifelse(is.null(mn.H0@optinfo$conv$lme4$messages),"",mn.H0@optinfo$conv$lme4$messages)
  M.mn.H1  <- ifelse(is.null(mn.H1@optinfo$conv$lme4$messages),"",mn.H1@optinfo$conv$lme4$messages)
  
  ##################
  ## output
  ##################
  output[[iteration]] <- data.frame(
    cbind(
      coefs      = rownames(S.fs.H0),
      iteration  = iteration,
      
      ## anova tables
      A.fs.H0       = A.fs.H0$`Pr(>F)`,
      A.fs.H1       = A.fs.H1$`Pr(>F)`,
      A.ps.H0       = A.ps.H0$`Pr(>F)`,
      A.ps.H1       = A.ps.H1$`Pr(>F)`,
      A.fi.H0       = A.fi.H0$`Pr(>F)`,
      A.fi.H1       = A.fi.H1$`Pr(>F)`,
      A.pi.H0       = A.pi.H0$`Pr(>F)`,
      A.pi.H1       = A.pi.H1$`Pr(>F)`,
      A.fm.H0       = A.fm.H0$`Pr(>F)`,
      A.fm.H1       = A.fm.H1$`Pr(>F)`,
      A.mn.H0       = A.mn.H0$`Pr(>F)`,
      A.mn.H1       = A.mn.H1$`Pr(>F)`,
      
      
      ## elapsed times for models
      et.fs.H0,
      et.fs.H1,
      et.ps.H0,
      et.ps.H1,
      et.fi.H0,
      et.fi.H1,
      et.pi.H0,
      et.pi.H1,
      et.fm.H0,
      et.fm.H1,
      et.mn.H0,
      et.mn.H1,
      
      
      ## elapsed times for anovas
      et.A.fs.H0,
      et.A.fs.H1,
      et.A.ps.H0,
      et.A.ps.H1,
      et.A.fi.H0,
      et.A.fi.H1,
      et.A.pi.H0,
      et.A.pi.H1,
      et.A.fm.H0,
      et.A.fm.H1,
      et.A.mn.H0,
      et.A.mn.H1,
      
      
      ## models singularity
      S.fs.H0,
      S.fs.H1,
      S.ps.H0,
      S.ps.H1,
      S.fi.H0,
      S.fi.H1,
      S.pi.H0,
      S.pi.H1,
      S.fm.H0,
      S.fm.H1,
      S.mn.H0,
      S.mn.H1,
      
      
      ## models convergence messages
      M.fM.H0,
      M.fM.H1,
      M.pM.H0,
      M.pM.H1,
      M.fi.H0,
      M.fi.H1,
      M.pi.H0,
      M.pi.H1,
      M.fm.H0,
      M.fm.H1,
      M.mn.H0,
      M.mn.H1
    )
  )
  
  write.csv2(do.call("rbind", output), file = "output-all-within-stimuli.csv")
}

sessionInfo()
#########################################