## Code description ####
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Author: Jesus C. Compaire
## Institution: Centro para el Estudio de Sistemas Marinos (CESIMAR)
## Position: Postdoctoral researcher
## Contact details: jesus.canocompaire@uca.es
## Date created: Dec-2023
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Code to replicate the statistical analyses and generate figures 
## performed in the manuscript:
## Compaire, J.C., Simionato C.G.  Moreira, D. & Acha, E.M. (2024).
## Modeling environmental effects on fishery landings: A case study of 
## whitemouth croaker (Micropogonias furnieri) in the Southwestern
## Atlantic shelf, 108806, https://doi.org/10.1016/j.ecss.2024.108806
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
## === FUNCTIONS TO PERFORM Generalized Linear Models (GLMs) ANALYSES ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to get all combinations among explanatory variables - ADDITIVE 
## -- -- -- -- -- -- -- -- -- -- -- -- --
comb_glm <- function(xnam, ynam){ 
  lstls <- list()
  cnum <- length(xnam) # length explanatory variables
  # List comprehension. Getting all combinations among explanatory variables
  lstls <- to_list(for(i in 1:cnum) 
    c(combn(xnam, i, simplify = FALSE)))
  for(i in 1:cnum){
    k <- length(lstls[[i]])
    for(j in 1:k){
      lstls[[i]][[j]] <- (as.formula(paste(
        paste(ynam), " ~ ", paste(lstls[[i]][[j]], collapse = "+"))))
    }
  }
  lst <- matrix(unlist(lstls), ncol = 1, byrow = TRUE)
  lst <- c(as.formula(paste(paste(ynam), " ~ 1")), lst) # Adding NULL model
}
# 
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to get all combinations among explanatory variables - INTERACCTION
## -- -- -- -- -- -- -- -- -- -- -- -- --
# considering the interaction of "x1" and "x2" with "x3"
comb_glmi <- function(xnam, ynam, inter){
  lstls <- list()
  cnum <- length(xnam) # length explanatory variables
  # List comprehension. Getting all combinations of explanatory variables
  lstls <- to_list(for(i in 1:cnum) 
    c(combn(xnam, i, simplify = FALSE)))
  for(i in 1:cnum){
    k <- length(lstls[[i]])
    for(j in 1:k){
      covariates <- lstls[[i]][[j]] # covariates names
      if (
        inter[1] %in% covariates &
        inter[2] %in% covariates & !(inter[3] %in% covariates)
      ){
        inter_t <- paste(inter[c(1,2)], collapse = "*")
      } else if (
        inter[1] %in% covariates &
        !inter[2] %in% covariates & (inter[3] %in% covariates)
      ){
        inter_t <- paste(inter[c(1,3)], collapse = "*")
      } else if (
        inter[1] %in% covariates &
        inter[2] %in% covariates & (inter[3] %in% covariates)
      ){
        inter_t <- paste(
          paste(inter[c(1,2)], collapse = "*"), "+",
          paste(inter[c(1,3)], collapse = "*"))
      } else {
        inter_t <- NULL
      }
      if (
        !is.null(inter_t)
      ){
        lstls[[i]][[j]] <- (as.formula(paste(
          paste(ynam), " ~ ", paste(covariates, collapse = "+"),
          " + ", inter_t)))
      } else {
        lstls[[i]][[j]] <- (as.formula(paste(
          paste(ynam), " ~ ", paste(covariates, collapse = "+"))))
      }
    }
  }
  lst <- matrix(unlist(lstls), ncol = 1, byrow = TRUE)
  lst <- c(as.formula(paste(paste(ynam), " ~ 1")), lst) # Adding NULL model
}
# 
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to perform GLM analyses
## -- -- -- -- -- -- -- -- -- -- -- -- --
try_glm <- function(eq, dataset, i, j, f, th){
  model <- tryCatch(
    glm(eq, family = f, data = dataset,
        na.action = na.omit),
    sink(stdout(), type = "message"),
    error = function(e){
      sink(message(paste0(
        "Analysis not performed, an error occurred for equation ", j, ":\n",
        eq, ".\n"), e), type = 'message')
    }
  )
  if (!is.null(model)){
    DE <- round(
      (1-model[["deviance"]]/model[["null.deviance"]])*100, digits = 2)
    K = length(model$coefficients)
    n = length(model$linear.predictors)
    AICc <- model$aic + ((2 * K * (K + 1) / (n - K -1)))
    pv <- summary(model)$coef[, "Pr(>|t|)"]
    covs <- attr(model$terms, "term.labels")
    cov_i <- unique(unlist(strsplit((covs[grepl(":", covs, fixed = T)]),':')))
    cov_s <- covs[-grep(pattern = ':', covs)]
    covs_d <- setdiff(cov_s, cov_i)
    thr <- 0.05
    if (
      (all(pv < thr) & length(pv) > 1) |
      (any(pv[2:length(pv)][grepl(":", covs, fixed = T)] < thr) &
       all(pv[covs_d] < thr))
    ){
      cat(paste0('percentile: ',
                 paste0(th),
                 '; family:',
                 paste0(f[["family"]],'(link=',f[["link"]],')', sep = ""),
                 '; lag:',
                 sprintf("%02d", i),
                 '; eq:',
                 sprintf("%02d", j),
                 '; DE:',
                 round(DE, digits = 2),
                 '; AICc:',
                 round(AICc, digits = 2), "\n"))
      message(as.character(eq))
      print(summary(model))
      print(check_collinearity(model))
      cat('---------- ######### ---------- ',"\n")
      output <- list('lag' = i, 'eq' = j,
                     'DE' =  DE, 'AICc' = AICc)
    } else {
      DE <- NA; AICc <- NA
      output <- list('lag' = i, 'eq' = j,
                     'DE' =  DE, 'AICc' = AICc)
    }
  }
  return(output)
}
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to check some GLMs manually
## -- -- -- -- -- -- -- -- -- -- -- -- --
handyf_glm <- function(k, ix, l, j){
  print(th[k])
  # List comprehension to get index (columns number) of variables
  vq <- c(th[k], 'total')
  idx_var <- to_vec(for(i in 1:length(vq))
    str_which(names_col, vq[i], negate = FALSE))
  names_var <- names_col[idx_var] # variable names complete
  # Getting climatic indices by pairs
  ci <- c('soi','sam','enso')
  cixs <- to_list(for(i in 1:length(ci)) 
    c(combn(ci, i, simplify = FALSE)))
  cixss <- list()
  cixss[[1]] <- cixs[[2]] # pairs of climatic indexes
  cixss[[1]][[4]] <- unlist(cixs[[3]]) # all climatic indexes
  names_var_ix <- NA # variable names short for each climatic index
  # Getting two of the climatic indices and 'th' landings to omit them
  cith <- paste(paste(cixss[[1]][[ix]], collapse = '|'), 'th', sep = '|')
  names_var_ix = names_var[
    -grep(pattern = cith, names_var)]
  # List comprehension to get index (columns number) of variables
  idx_ss <- to_vec(for(i in 1:length(names_var_ix))
    str_which(names_col, names_var_ix[i], negate = FALSE))
  # Removing tag "th"
  tag <- c(th[k])
  for(i in 1:length(tag)){
    names_var_ix <- gsub(tag[i],'', names_var_ix)
  }
  print(names_var_ix[1])
  # Covariates and response variable names
  xnames <- head(names_var_ix, -1); ynames <- tail(names_var_ix, 1)
  # GLMs formulas ADDITIVE glms_f <- comb_glm(xnames, ynames) 
  # GLMs formulas INTERACTION
  inter <- c("river", "kd", "sss") # interaction terms
  glms_f <- comb_glmi(xnames, ynames, inter) # getting equations
  if(
    ix < length(cixss[[1]])
  ){
    glms_f <- glms_f[grep(
      pattern = paste(cixss[[1]][[4]], collapse = '|'), (glms_f))]
  } else {
    glms_f <- glms_f
  }
  spq <- t(sp[c(idx_ss), , l]) # dataset to each lag
  spq <- as.data.frame(spq)
  # spq$V1 <- as.factor(spq$V1) # climatic index as qualitative variable
  colnames(spq) <- names_var_ix 
  eq <- glms_f[[j]]
  my_list <- list(spq, eq)
  return(my_list)
}
# 
## -- -- -- -- -- -- -- -- -- -- -- -- --
## === FUNCTIONS TO PERFORM Generalized Additive Models (GAMs) ANALYSES ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
# f(x) to get all combinations among explanatory variables - ADDITIVE
## -- -- -- -- -- -- -- -- -- -- -- -- --
comb_gam <- function(xnam, ynam, kdf){ 
  lstls <- list()
  cnum <- length(xnam) # length explanatory variables
  # List comprehension. Getting all combinations among explanatory variables
  lstls <- to_list(for(i in 1:cnum) 
    c(combn(xnam, i, simplify = FALSE)))
  for(i in 1:cnum){
    k <- length(lstls[[i]])
    for(j in 1:k){
      lstls[[i]][[j]] <- (as.formula(paste(
        paste(ynam), " ~ ", 
        paste("s(", lstls[[i]][[j]], 
              ", k = ", kdf, ", bs = 'cr')", collapse = "+"))))
    }
  }
  lst <- matrix(unlist(lstls), ncol = 1, byrow = TRUE)
  lst <- c(as.formula(paste(paste(ynam), " ~ 1")), lst) # Adding NULL model
}
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to get all combinations among explanatory variables - INTERACCTION
## -- -- -- -- -- -- -- -- -- -- -- -- --
# considering the interaction of "x1" and "x2" with "x3"
comb_gami <- function(xnam, ynam, kdf, inter){
  lstls <- list()
  cnum <- length(xnam) # length explanatory variables
  # List comprehension. Getting all combinations of explanatory variables
  lstls <- to_list(for(i in 1:cnum) 
    c(combn(xnam, i, simplify = FALSE)))
  for(i in 1:cnum){
    k <- length(lstls[[i]])
    for(j in 1:k){
      covariates <- lstls[[i]][[j]] # covariates names
      if (
        inter[1] %in% covariates &
        inter[2] %in% covariates & !(inter[3] %in% covariates)
      ){
        inter_t <- paste(
          "ti(", paste(inter[c(1,2)],
                       collapse = ", "), ", k = ", kdf, ", bs = 'cr')")
      } else if (
        inter[1] %in% covariates &
        !inter[2] %in% covariates & (inter[3] %in% covariates)
      ){
        inter_t <- paste(
          "ti(", paste(inter[c(1,3)],
                       collapse = ", "), ", k = ", kdf, ", bs = 'cr')")
      } else if (
        inter[1] %in% covariates &
        inter[2] %in% covariates & (inter[3] %in% covariates)
      ){
        inter_t <- paste(
          paste("ti(", paste(inter[c(1,2)],
                             collapse = ", "), ", k = ", kdf, ", bs = 'cr')"),
          "+",
          paste("ti(", paste(inter[c(1,3)],
                             collapse = ", "), ", k = ", kdf, ", bs = 'cr')"))
      } else {
        inter_t <- NULL
      }
      if (
        !is.null(inter_t)
      ){
        lstls[[i]][[j]] <- (as.formula(paste(
          paste(ynam), " ~ ",
          paste("s(", lstls[[i]][[j]], 
                ", k = ", kdf, ", bs = 'cr')", collapse = "+"),
          "+", inter_t)))
      } else {
        lstls[[i]][[j]] <- (as.formula(paste(
          paste(ynam), " ~ ",
          paste("s(", lstls[[i]][[j]], 
                ", k = ", kdf, ", bs = 'cr')", collapse = "+"))))
      }
    }
  }
  lst <- matrix(unlist(lstls), ncol = 1, byrow = TRUE)
  lst <- c(as.formula(paste(paste(ynam), " ~ 1")), lst) # Adding NULL model
}
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to perform GAM analyses
## -- -- -- -- -- -- -- -- -- -- -- -- --
try_gam <- function(eq, dataset, i, j, f, th){
  model <- tryCatch(
    gam(eq, 
        family = f,
        data = dataset,
        na.action = na.omit, method = "REML"),
    error = function(e){
    }
    # To print error message replace previous two lines with the next ones
    # sink(stdout(), type = "message"),
    # error = function(e){
    # message(paste0(
    # # sink(message(paste0(
    #   "Analysis not performed, an error occurred for equation ", j, ":\n",
    #   eq, ".\n"), e, type = 'message')
    # }
  )
  if (!is.null(model)){
    DE <- round(
      (1-model[["deviance"]]/model[["null.deviance"]])*100, digits = 2)
    # K = length(model$coefficients)
    # n = length(model$linear.predictors)
    AICc <- MuMIn::AICc(model)
    pvi <- summary.gam(model)[["p.pv"]]
    pvc <- summary(model)$s.table[, "p-value"]
    pv <- c(pvi, pvc)
    names(pv) <- c("(Intercept)", rownames(summary(model)$s.table))
    if (length(pv) <= 1)
    {
      DE <- NA; AICc <- NA
      output <- list('lag' = i, 'eq' = j,
                     'DE' =  DE, 'AICc' = AICc)
    } else {
      covs <- rownames(summary(model)$s.table)
      cov_i <- unique(unlist(strsplit((covs[grepl("ti", covs, fixed = T)]),',')))
      cov_i <- gsub('ti','', cov_i)
      cov_i <- gsub('[^[:alnum:] ]','', cov_i) 
      cov_s <- covs[startsWith(covs, 's(') ] 
      cov_s <- gsub("[\\(\\)]", "",
                    regmatches(cov_s, gregexpr("\\(.*?\\)", cov_s)))
      covs_d <- setdiff(cov_s, cov_i)
      if (length(covs_d) > 0)
      {
        cnum <- length(covs_d)
        covs_dl <- to_list(for(i in 1:cnum)
          covs[grepl(covs_d[i], covs, fixed = T)]
        )
        covs_dpv <- unlist(covs_dl)
      } else {
        covs_dpv <- covs_d
      }
      thr <- 0.05
      if (
        (all(pv < thr) & length(pv) > 1) |
        (any(pv[2:length(pv)][covs[grepl("ti", covs, fixed = T)]] < thr) &
         all(pv[covs_dpv] < thr))
      ){
        cat(paste0('percentile:',
                   paste0(th),
                   '; family:',
                   paste0(f[["family"]],'(link=',f[["link"]],')', sep = ""),
                   '; lag:',
                   sprintf("%02d", i),
                   '; eq:',
                   sprintf("%02d", j),
                   '; DE:',
                   round(DE, digits = 2),
                   '; AICc:',
                   round(AICc, digits = 2), "\n"))
        # print(as.character(eq))
        print(summary(model))
        print(concurvity(model, full = T))
        cat('---------- ######### ---------- ',"\n")
        output <- list('lag' = i, 'eq' = j,
                       'DE' =  DE, 'AICc' = AICc)
      } else {
        DE <- NA; AICc <- NA
        output <- list('lag' = i, 'eq' = j,
                       'DE' =  DE, 'AICc' = AICc)
      }
    }
  }
  return(output)
}
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to check some GAMs manually
## -- -- -- -- -- -- -- -- -- -- -- -- --
handyf_gam <- function(k, ix, l, j){
  print(th[k])
  # List comprehension to get index (columns number) of variables
  vq <- c(th[k], 'total')
  idx_var <- to_vec(for(i in 1:length(vq))
    str_which(names_col, vq[i], negate = FALSE))
  names_var <- names_col[idx_var] # variable names complete
  # Getting climatic indices by pairs
  ci <- c('soi','sam','enso')
  cixs <- to_list(for(i in 1:length(ci)) 
    c(combn(ci, i, simplify = FALSE)))
  cixss <- list()
  cixss[[1]] <- cixs[[2]] # pairs of climatic indexes
  cixss[[1]][[4]] <- unlist(cixs[[3]]) # all climatic indexes
  names_var_ix <- NA # variable names short for each climatic index
  # Getting two of the climatic indices and 'th' landings to omit them
  cith <- paste(paste(cixss[[1]][[ix]], collapse = '|'), 'th', sep = '|')
  names_var_ix = names_var[
    -grep(pattern = cith, names_var)]
  # List comprehension to get index (columns number) of variables
  idx_ss <- to_vec(for(i in 1:length(names_var_ix))
    str_which(names_col, names_var_ix[i], negate = FALSE))
  # Removing tag "th"
  tag <- c(th[k])
  for(i in 1:length(tag)){
    names_var_ix <- gsub(tag[i],'', names_var_ix)
  }
  print(names_var_ix[1])
  # Covariates and response variable names
  xnames <- head(names_var_ix, -1); ynames <- tail(names_var_ix, 1)
  # GAMs formulas ADDITIVE gams_f <- comb_gam(xnames, ynames, kdf) 
  # GAMs formulas INTERACTION
  inter <- c("river", "kd", "sss") # interaction terms
  gams_f <- comb_gami(xnames, ynames, "4", inter) # getting equations
  if(
    ix < length(cixss[[1]])
  ){
    gams_f <- gams_f[grep(
      pattern = paste(cixss[[1]][[4]], collapse = '|'), (gams_f))]
  } else {
    gams_f <- gams_f
  }
  spq <- t(sp[c(idx_ss), , l]) # dataset to each lag
  spq <- as.data.frame(spq)
  # spq$V1 <- as.factor(spq$V1) # climatic index as qualitative variable
  colnames(spq) <- names_var_ix 
  eq <- gams_f[[j]]
  my_list <- list(spq, eq)
  return(my_list)
  
}
# 
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to get x-axis values with specific precision and number of intervals
## -- -- -- -- -- -- -- -- -- -- -- -- --
xaxs <- function(xval, idx = i, precision, nintervals){
  xv <- xval[[i]][['x']]
  mn <- round_any(min(xv), precision, f = floor)
  mx <- round_any(max(xv), precision, f = ceiling)
  xvalues <- list(c(mn, mx, (mx-mn)/nintervals))
}
# 
## === FUNCTIONS TO SUMMARISE OUTPUTS FROM GLMs & GAMs ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to visualize number of significant model per index and their DE values
## -- -- -- -- -- -- -- -- -- -- -- -- --
vis_n_de <- function(model_summary, fmly, th, ix){
  n  <- to_vec(for(i in 1:length(ix))
    sapply(model_summary[[fmly]][[th]][[ix[i]]][1], NROW))
  de <- (to_vec(for(i in 1:length(ix))
    model_summary[[fmly]][[th]][[ix[i]]][['DE']]))
  return(list(n, de))
}
# 
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to get a vector containing lag or DE according to family and index
## -- -- -- -- -- -- -- -- -- -- -- -- --
vector_de_lag <- function(lst, fmly, th, ix, vrbl){
  output <- to_vec(for(i in 1:length(ix))
    (lst[[fmly]][[th]][[ix[i]]][[vrbl]]))
  return(output)
}
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to evaluate deviance explained (max, mean, sd) by each lag 
## -- -- -- -- -- -- -- -- -- -- -- -- --
table_de_lag <- function(dev, lagv, tag){
  df <- data.frame(lagv, dev)
  names(df) <- c('lag','DE')
  dfmx <- aggregate(df$DE, list(df$lag), FUN=max) 
  dfav <- aggregate(df$DE, list(df$lag), FUN=mean) 
  dfsd <- aggregate(df$DE, list(df$lag), FUN=sd)
  dfn  <- aggregate(df$DE, list(df$lag), FUN=length)
  df2  <- data.frame(dfmx, dfav$x, dfsd$x, dfn$x)
  names(df2) <- c('lag', 'max','avg','sd', 'n')
  df2$error <- tag
  # Add NA elements where there is no significant models
  dfna = data.frame(matrix(NA, nrow = 11, ncol = 5))
  dfna[,1] <- c(0:10) # numer of lags
  names(dfna) <- c('lag', 'max','avg','sd', 'n')
  dfna$error <- tag
  ldif <- dfna[,1][!dfna[,1] %in% df2[,1]]
  df3 <- rbind(df2,dfna[ldif+1,])
  output <- df3[order(df3$lag,decreasing=F),]
  rownames(output) <- NULL 
  return(output)
}
#
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to evaluate the contribution (deviance explained) of each variable
## -- -- -- -- -- -- -- -- -- -- -- -- --
de_ev <- function(model_complete, model_simplified){
  round(
    (model_simplified$deviance - model_complete$deviance) /
      model_complete$null.deviance * 100, digits = 2)
}
## === FIGURES ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## Customize Figures
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
## Fig. 1 - MAP ####
map_aucfz <- function(df, plygn = NULL, plygn2 = NULL,
                      scale = NULL, clabs = NULL, ...){
  if (!is.data.frame(df))
    stop("'df' must be a dataframe containing longitude, latitude and the
         elevation/depth values")
  if (is.null(plygn)) {
    ppolygon <- geom_polygon()
  } else {
    if (!is.data.frame(plygn)) {
      stop("'polygon' must be a dataframe containing longitude and latitude
           values")
    } else {
      ppolygon <- geom_polygon(data = plygn, aes(x = plygn[,1], y = plygn[,2]),
                               colour = "grey25", fill = NA, 
                               linetype = 5, linewidth = 0.5)
    }
  }
  if (is.null(plygn2)) {
    ppolygon2 <- geom_polygon()
  } else {
    if (!is.data.frame(plygn2)) {
      stop("'polygon' must be a dataframe containing longitude and latitude
           values")
    } else {
      ppolygon2 <- geom_polygon(data = plygn2, 
                                aes(x = plygn2[,1], y = plygn2[,2]),
                                colour = "grey25", fill = NA, 
                                linetype = 1, linewidth = 0.5)
    }
  }
  if (is.null(scale) || scale == "discrete") {
    if (is.null(clabs)) {
      stop("Default 'scale_fill' is discrete, so 'clabs' must be defined as a
    character vector that contains depth scale values. If you want to draw the
    figure using a continuous scale, please choose 'scale = continuous'")
    } else {
      # Discrete scale
      clabs <- clabs
      depth <- df[,4]
      scl <- scale_fill_manual(values =
                                 rev(c(blues9[2:length(unique(df[,4]))])),
                               na.value =  "gray75",
                               limits = clabs, labels = rev(clabs),
                               guide = guide_legend(reverse = TRUE))
    }
  } else if (!is.null(scale) & scale == "continuous") {
    if(!is.null(clabs)){
      warning("omiting unused argument 'clabs'")
    } else {
    }
    # Continuous scale
    depth <- df[,3]
    scl <- scale_fill_distiller(palette="Blues", na.value="gray75",
                                limits = c(min(df[,3]), 0))
  } else if (!is.null(scale) || scale != "discrete" || scale != "continuous"){
    stop("The argument scale must be 'discrete' or 'continuous'")
  }
  # Nested function of the ggsn package modified to use "Times" as font family
  rulermap <- function(data = NULL, location = "bottomright", dist = NULL,
                       dd2km = NULL, model = NULL, height = 0.02, 
                       st.dist = 0.02, dist_unit = "km", st.bottom = TRUE,
                       st.size = 5, st.color = "black", st.family = "Times",
                       box.fill = c("black", "white"), 
                       box.color = "black", border.size = 1, 
                       x.min = NULL, x.max = NULL, y.min = NULL, y.max = NULL,
                       anchor = NULL, facet.var = NULL, facet.lev = NULL, ...){
    if (is.null(data)) {
      if (is.null(x.min) | is.null(x.max) |
          is.null(y.min) | is.null(y.max) ) {
        stop('If data is not defined, x.min, x.max, y.min and y.max must be.')
      }
      data <- data.frame(long = c(x.min, x.max), lat = c(y.min, y.max))
    }
    if (is.null(dd2km)) {
      stop("dd2km should be logical.")
    }
    if (any(class(data) %in% "sf")) {
      xmin <- sf::st_bbox(data)["xmin"]
      xmax <- sf::st_bbox(data)["xmax"]
      ymin <- sf::st_bbox(data)["ymin"]
      ymax <- sf::st_bbox(data)["ymax"]
    } else {
      xmin <- min(data$long)
      xmax <- max(data$long)
      ymin <- min(data$lat)
      ymax <- max(data$lat)
    }
    if (location == 'bottomleft') {
      if (is.null(anchor)) {
        x <- xmin
        y <- ymin
      } else {
        x <- as.numeric(anchor['x'])
        y <- as.numeric(anchor['y'])
      }
      direction <- 1
      
    }
    if (location == 'bottomright') {
      if (is.null(anchor)) {
        x <- xmax
        y <- ymin
      } else {
        x <- as.numeric(anchor['x'])
        y <- as.numeric(anchor['y'])
      }
      direction <- -1
      
    }
    if (location == 'topleft') {
      if (is.null(anchor)) {
        x <- xmin
        y <- ymax
      } else {
        x <- as.numeric(anchor['x'])
        y <- as.numeric(anchor['y'])
      }
      direction <- 1
      
    }
    if (location == 'topright') {
      if (is.null(anchor)) {
        x <- xmax
        y <- ymax
      } else {
        x <- as.numeric(anchor['x'])
        y <- as.numeric(anchor['y'])
      }
      direction <- -1
      
    }
    if (!st.bottom) {
      st.dist <-
        y + (ymax - ymin) * (height + st.dist)
    } else {
      st.dist <- y - (ymax - ymin) * st.dist
    }
    height <- y + (ymax - ymin) * height
    
    if (dd2km) {
      if (dist_unit == "m") {
        dist <- dist / 1e3
      }
      break1 <- maptools::gcDestination(lon = x, lat = y,
                                        bearing = 90 * direction,
                                        dist = dist, dist.units = 'km',
                                        model = model)[1, 1]
      break2 <- maptools::gcDestination(lon = x, lat = y,
                                        bearing = 90 * direction,
                                        dist = dist*2, dist.units = 'km',
                                        model = model)[1, 1]
    } else {
      if (location == 'bottomleft' | location == 'topleft') {
        break1 <- x + dist * 1e3
        break2 <- x + dist * 2e3
      } else {
        break1 <- x - dist * 1e3
        break2 <- x - dist * 2e3
      }
      
    }
    box1 <- data.frame(x = c(x, x, rep(break1, 2), x),
                       y = c(y, height, height, y, y), group = 1)
    box2 <- data.frame(x = c(rep(break1, 2), rep(break2, 2), break1),
                       y=c(y, rep(height, 2), y, y), group = 1)
    if (!is.null(facet.var) & !is.null(facet.lev)) {
      for (i in 1:length(facet.var)){
        if (any(class(data) == "sf")) {
          if (!is.factor(data[ , facet.var[i]][[1]])) {
            data[ , facet.var[i]] <- factor(data[ , facet.var[i]][[1]])
          }
          box1[ , facet.var[i]] <- factor(facet.lev[i],
                                          levels(data[ , facet.var[i]][[1]]))
          box2[ , facet.var[i]] <- factor(facet.lev[i],
                                          levels(data[ , facet.var[i]][[1]]))
        } else {
          if (!is.factor(data[ , facet.var[i]])) {
            data[ , facet.var[i]] <- factor(data[ , facet.var[i]])
          }
          box1[ , facet.var[i]] <- factor(facet.lev[i],
                                          levels(data[ , facet.var[i]]))
          box2[ , facet.var[i]] <- factor(facet.lev[i],
                                          levels(data[ , facet.var[i]]))
        }
        
      }
    }
    if (dist_unit == "km") {
      legend <- cbind(text = c(0, dist, dist * 2), row.names = NULL)
    }
    if (dist_unit == "m") {
      legend <- cbind(text = c(0, dist * 1e3, dist * 2e3), row.names = NULL)
    }
    
    gg.box1 <- geom_polygon(data = box1, aes(x, y),
                            fill = utils::tail(box.fill, 1),
                            color = utils::tail(box.color, 1),
                            size = border.size)
    gg.box2 <- geom_polygon(data = box2, aes(x, y), fill = box.fill[1],
                            color = box.color[1],
                            size = border.size)
    x.st.pos <- c(box1[c(1, 3), 1], box2[3, 1])
    if (location == 'bottomright' | location == 'topright') {
      x.st.pos <- rev(x.st.pos)
    }
    if (dist_unit == "km") {
      legend2 <- cbind(data[1:3, ], x = unname(x.st.pos), y = unname(st.dist),
                       label = paste0(legend[, "text"], c("", "", "km")))
    }
    if (dist_unit == "m") {
      legend2 <- cbind(data[1:3, ], x = unname(x.st.pos), y = unname(st.dist),
                       label = paste0(legend[, "text"], c("", "", "m")))
    }
    if (!is.null(facet.var) & !is.null(facet.lev)) {
      for (i in 1:length(facet.var)){
        if (any(class(data) == "sf")) {
          legend2[ , facet.var[i]] <- factor(facet.lev[i],
                                             levels(data[ , facet.var[i]][[1]]))
        } else {
          legend2[ , facet.var[i]] <- factor(facet.lev[i],
                                             levels(data[ , facet.var[i]]))
        }
      }
    } else if (!is.null(facet.var) & is.null(facet.lev)) {
      facet.levels0 <- unique(as.data.frame(data)[, facet.var])
      facet.levels <- unlist(unique(as.data.frame(data)[, facet.var]))
      legend2 <- do.call("rbind", replicate(length(facet.levels),
                                            legend2, simplify = FALSE))
      if (length(facet.var) > 1) {
        facet.levels0 <- expand.grid(facet.levels0)
        legend2[, facet.var] <-
          facet.levels0[rep(row.names(facet.levels0), each = 3), ]
      } else {
        legend2[, facet.var] <- rep(facet.levels0, each = 3)
      }
    }
    gg.legend <- geom_text(data = legend2, aes(x, y, label = label),
                           size = st.size, color = st.color,
                           family = st.family, ...)
    return(list(gg.box1, gg.box2, gg.legend))
  }
  # Coordinate numbers to plot on x- and y-axes
  xnum <- rev(seq(round(-max(df[,1])),
                  round(-min(df[,1])), by = 2))
  ynum <- rev(seq(round_any(-max(df[,2]), 1, ceiling),
                  round_any(-min(df[,2]), 1, ceiling), by = 2))
  ggplot() +
    ## Bathymetry
    geom_raster(aes(df[,1], df[,2], fill = depth), data = df) +
    coord_equal() +
    scl +
    labs(fill = "Depth (m)") +
    ## Polygon
    ppolygon +
    ## Polygon2
    ppolygon2 +
    ## Map ruler
    rulermap(x.min = -59.2, x.max = -59.2, y.min = -33.6, y.max = -33.7,
             location = "topleft", dist = 50, dist_unit = "km",
             dd2km = TRUE, model = "WGS84", 
             height = 0.7, st.dist = 0.9, st.size = 3.5, st.family = "Times",
             box.fill = c("black", "white"),
             box.color = "black", border.size = 0.5) +
    ## North arrow
    annotation_north_arrow(location = "tl",
                           pad_x = unit(0.2, "cm"),
                           pad_y = unit(0.25, "cm"),
                           style = north_arrow_fancy_orienteering) +
    ## Country names
    annotate(
      geom = "text", x = -58.7, y = -36.2, label = "Argentina", 
      fontface = "plain", family = "Times", color = "grey22", size = 6) +
    annotate(
      geom = "text", x = -56, y = -33.8, label = "Uruguay", 
      fontface = "plain", family = "Times", color = "grey22", size = 6) +
    annotate(
      geom = "text", x = -53.2, y = -33.3, label = "Brazil", 
      fontface = "plain", family = "Times", color = "grey22", size = 3) +
    ## Rio de la Plata name
    annotate(
      geom = "text", x = -57.8, y = -34.7, label = "Rio de la Plata", 
      fontface = "plain", family = "Times", color = "grey22", size = 4.5,
      angle = -31) +
    ## Punta del Este name
    annotate(
      geom = "text", x = -54.95, y = -34.75, label = "PE", 
      fontface = "plain", family = "Times", color = "grey22", size = 4.5,
      angle = 0) +
    ## Punta Rasa del Cabo San Antonio Plata name
    annotate(
      geom = "text", x = -56.9, y = -36.5, label = "PR", 
      fontface = "plain", family = "Times", color = "grey22", size = 4.5,
      angle = 0) +
    ## Axes and theme
    scale_y_continuous(expand=c(0,0),
                       sec.axis = dup_axis(), # right axis labels
                       breaks = -ynum,
                       labels = to_vec(for(i in 1:length(ynum))
                         (paste0(ynum[i], "\u00B0S")))) +
    scale_x_continuous(expand=c(0,0),
                       sec.axis = dup_axis(), # top axis labels
                       breaks = -xnum,
                       labels = to_vec(for(i in 1:length(xnum))
                         (paste0(xnum[i], "\u00B0W")))) +
    theme_bw()+
    theme(
      panel.grid.major = element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(vjust = -0.5, color = "black"),
      axis.text.y = element_text(hjust = -0.5, color = "black"),
      axis.text.x.top = element_blank(), # do not show top / right axis labels
      axis.text.y.right = element_blank(),
      axis.ticks.length = unit(-2, "mm"),
      legend.text = element_text(size = 12),
      legend.key = element_rect(fill = "white", colour = "transparent"),
      legend.position = c(.983, .33),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(5, 5, 5, 5)) +
    theme(
      text = element_text(family = "Times"))
}
# 
## Figs. 2 & 3 - ERROR DISTRIBUTIONS COMPARISON ####
plot_de_lag <- function(df){
  if (!is.data.frame(df))
    stop("'df' must be a dataframe containing DE values at different time langs
         for both error distributions")
  ggplot(df, aes(x = factor(lag, levels = seq(0,10,1)),
                 y = max, fill = error)) +
    geom_bar(stat = "identity", color = "black",
             position = position_dodge()) +
    # geom_errorbar(aes(ymin  = avg - sd, ymax  = avg + sd),
    #               width = 0.5, size  = 0.5, position = position_dodge(.9)) +
    geom_point(aes(x = factor(lag, levels = seq(0,10,1)),
                   y = avg),
               position = position_dodge(0.9), size = 3,
               show.legend=FALSE) +
    scale_fill_manual(values = c("white", "grey50")) +
    coord_cartesian(ylim = c(0,100)) +
    scale_y_continuous(breaks=c(0, 20, 40, 60, 80, 100), expand = c(0, 0)) +
    # geom_text(aes(x = lag, y = 5,
    #          label= n,
    #          na.rm = T,
    #          group = error),
    #          position = position_dodge(0.9)) +
    theme_tufte() +
    theme(axis.title = element_text(face = "bold")) +
    theme(legend.title=element_blank(),
          legend.position= "none" ) +
    theme(axis.line = element_line(colour = "black"),
         axis.ticks = element_line(colour = "black") +
    theme(axis.text.y = element_text(color = "black", size = 12)) +
    theme(axis.text.x = element_text(color = "black", size = 12)) +
    labs(x = "", y = "") +
    theme(axis.text=element_text(size=12),
          axis.title.y =  element_blank(),
          axis.title= element_blank())
}
# 
## Fig. 4 - RESPONSE PLOT - GLM ####
glm_effect_single.plot <- function(var_x, model,
                                   xmin, xmax, interval,
                                   yymin, yymax,
                                   col_line){
  if (!is.character(var_x))
    stop("'var_x' must be a character string")
  if (!is.list(model))
    stop("'model' must be a list containing the 'glm' function output")
  if (!is.double(xmin) || !is.double(xmin))
    stop("'xmin and xmax' must be double values")
  if (!is.list(intervals))
    stop("'intervals' must be a list of numeric vectors")
  if (!is.character(col_line))
    stop("'col_line' must be a character string containing colour pattern")
  # 
  xl <- list(seq(xmin, xmax, interval))
  names(xl) <- var_x
  obj <- effect(var_x, model,
                xlevels = xl)
  df <- data.frame(
    x = obj[["x"]][[var_x]],
    y = 1/(obj[["fit"]]),
    ymin = 1/obj[["upper"]],
    ymax = 1/obj[["lower"]]
  )
  # 
  ggplot(data = df, aes(x = x, y = y)) +
    geom_line(col = col_line[1]) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax),
                alpha = 0.5, fill = col_line[1]) +
    geom_rug(
      data = obj[["data"]], aes(
        x = obj[["data"]][[var_x]],
        y = Inf),
      col = "black", sides = "b", size = 0.3) +
    scale_x_continuous(breaks = seq(xmin, xmax, by = interval*2)) +
    ylim(yymin, yymax) +  
    labs(y = "Landings (t)", x = var_x) +
    theme_tufte() +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    theme(
      legend.title = element_blank(),
      legend.justification = c(0, 1),
      legend.position = c(0.05, 0.90),
      legend.direction = "vertical",
      legend.spacing.x = unit(0.25, 'cm'),
      legend.spacing.y = unit(1, 'cm'),
      legend.key.size = unit(1, "cm"),
      legend.background = element_blank(),
      legend.text = element_text(
        colour = "black",
        size = 14,
        face = "italic"
      )
    )
}
# 
## Fig. 5 - RESPONSE PLOT - GAM ####
gam_effect_single.plot <- function(model, intercept, var_x,
                                   xmin, xmax, interval,
                                   yymin, yymax,
                                   col_line){
  if (!is.list(model))
    stop("'model' must be a list containing the 'gam' function output")
  if (!is.numeric(intercept))
    stop("'intercept' must be a numeric vector")
  if (!is.character(var_x))
    stop("'var_x' must be a character string")
  if (!is.double(xmin) || !is.double(xmin))
    stop("'xmin and xmax' must be double values")
  if (!is.list(intervals))
    stop("'intervals' must be a list of numeric vectors")
  if (!is.character(col_line))
    stop("'col_line' must be a character string containing colour pattern")
  # 
  plot_elements <- plot(model, select = var_x, xlim = c(xmin, xmax))
  idx <- 1:length(plot_elements)
  new_names <- to_list(for(i in 1:length(idx))
    plot_elements[[i]][["xlab"]])
  names(plot_elements) <- unlist(new_names)
  rug = plot_elements[[var_x]]$raw
  df <- data.frame(
    x = plot_elements[[var_x]]$x,
    y = 1/(intercept + plot_elements[[var_x]]$fit),
    ymin = 1/(intercept + (plot_elements[[var_x]]$fit - plot_elements[[var_x]]$se)),
    ymax = 1/(intercept + (plot_elements[[var_x]]$fit + plot_elements[[var_x]]$se)),
    rug = c(rug, rep(NA, length(plot_elements[[var_x]][["x"]]) - length(rug)))
  )
  # 
  ggplot(data = df, aes(x = x, y = y)) +
    geom_line(col = col_line[1]) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax),
                alpha = 0.5, fill = col_line[1]) +
    geom_rug(
      data = df,
      aes(
        x = rug,
        y = Inf),
      col = "black", sides = "b", size = 0.3, inherit.aes = FALSE) +
    # geom_hline(yintercept = 1/intercept, col = "black", alpha = 0.7, lwd = 0.7, lty = 2) +
    scale_x_continuous(breaks = seq(xmin, xmax, by = interval*2),
                       labels = label_number(accuracy = 0.01)) +
    ylim(yymin, yymax) +
    labs(y = "Landings (t)", x = var_x) +
    theme_tufte() +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    theme(
      legend.title = element_blank(),
      legend.justification = c(0, 1),
      legend.position = c(0.05, 0.90),
      legend.direction = "vertical",
      legend.spacing.x = unit(0.25, 'cm'),
      legend.spacing.y = unit(1, 'cm'),
      legend.key.size = unit(1, "cm"),
      legend.background = element_blank(),
      legend.text = element_text(
        colour = "black",
        size = 14,
        face = "italic"
      )
    )
}
# 
## Figs. S1A & S1B - OBSERVED vs FITTED VALUES ####
plot_obs_fit <- function(obsfit, model_pal, xp, yp, upperidx = NULL){
  if (!is.list(obsfit))
    stop("'obsfit' must be a list containing the 'appraise' function output")
  if (!is.list(model_pal))
    stop("'model_pal' must be a list containing the 'lm' function output")
  if (!is.double(xp) || !is.double(yp))
    stop("'xp and yp' must be double values")
  if (!is.null(upperidx) && !is.character(upperidx))
    upperidx <- as.character(upperidx)
  # 
  ggplot(as.data.frame(obsfit["data"]),
         aes(x = (obsfit[["data"]][[2]]), y = (obsfit[["data"]])[[1]])) + 
    geom_point() +
    # adjustment line
    stat_smooth(method = "lm", se = F, fullrange = T, color = "black") +
    # perfect agreement line
    geom_line(data = model_pal,
              aes(x = model_pal[["model"]][[2]],
                  y = model_pal[["model"]][[1]]),
              lty = "dotted", lwd = 1, col = "black") +
    # r and p-value 
    annotate(
      geom = "text", x = xp, y = yp, 
      label = paste0("r-pearson = ", r2, ";  p " , pv),
      fontface = "plain", family = "Times", color = "black", size = 5
    ) + 
    theme_tufte() +
    labs(x = "Fitted landings (t)",
         y = "Observed landings (t)") +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    ggtitle(upperidx) +
    theme(plot.title = element_text(size=23, face = "bold",
                                    hjust = 0.01, vjust = -4))
}
# 
