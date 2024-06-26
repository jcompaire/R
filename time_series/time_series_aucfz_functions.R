## Code description ####
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Author: Jesus C. Compaire
## Institution: Centro de Investigaciones del Mar y la Atmósfera (CIMA)
## Position: Postdoctoral researcher
## Contact details: jesus.canocompaire@uca.es
## Date created: Aug-2023
## -- -- -- -- -- -- -- -- -- -- -- -- --
##
## Code to replicate the statistical analyses and generate figures 
## performed in the manuscript:
## Compaire, J.C., Acha, E.M., Moreira, D. & Simionato C.G. (2024).
## Time series modeling of coastal fishery landings on the Southwestern
## Atlantic shelf: influence of environmental drivers
##
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
## === ALGORITHM TO COMPARE ARIMA MODELS ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## Running different ARIMA models based on all combinations between AR and
## MA parameters provided without exceeding of six parameters.
## The output shows AIC values of those models whose
## p-value is greater than 0.05 (which means that so residuals are independent)
## -- -- -- -- -- -- -- -- -- -- -- -- --
run_models <- function(y, pi, qi, d, Pi, Qi, D, per){
  p=NULL; q=NULL; P=NULL; Q=NULL # Loop elements
  for(p in 1:length(pi)){
    for(q in 1:length(qi)){
      for(P in 1:length(Pi)){
        for(Q in 1:length(Qi)){
          if(p+d+q+P+D+Q<=10){
            model<-arima(x=y,
                         order = c((p-1),d,(q-1)),
                         seasonal = list(order=c((P-1),D,(Q-1)),
                                         period=per))
            pval<-Box.test(model$residuals, lag=log(length(model$residuals)),
                           type = "Ljung-Box")
            nT <- model[["nobs"]]
            AICc <- model$aic + (2*(p+q+1)*(p+q+2))/(nT-p-q-2)
            if (pval$p.value >= 0.05) {
              cat(p-1,d,q-1,P-1,D,Q-1,per,
                  'AICc=', AICc,
                  ' p-VALUE=', pval$p.value,'\n')
            }
          }
        }
      }
    }
  }
}
# 
## === BOX-COX TRANSFORMATION ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## Run a Box Cox power transformation 
## -- -- -- -- -- -- -- -- -- -- -- -- --
tran_boxcox <- function(x,lambda){
  if (lambda == 0) {
    xt_box <- log(x)
  } else {
    xt_box <- (x^lambda - 1)/lambda
  }
}
#
## === WAVELET ANALYSIS ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## f(x) to perform wavelet coherency analysis between two time series
## -- -- -- -- -- -- -- -- -- -- -- -- --
coherency_analysis <- function(dataframe, vars_key){
  df <- dataframe
  wv_list <- list()
  for (j in 1:length(vars_key)){
    print(paste0(vars_key[j],'... ',j,'/',length(vars_key)))
    my.data <- as.data.frame(cbind(
      x = df[[vars_key[j]]],
      y = df$landings))
    my.data[["date"]] <- as.POSIXct(df$date)
    my.data <- na.omit(my.data)
    print(my.data$date[1])
    my.wc <- analyze.coherency(my.data, my.pair = c("x","y"),
                               loess.span = 0, dt = 1, dj = 1/100,
                               lowerPeriod = 2,
                               upperPeriod = 128,
                               window.type.t = 2, window.type.s = 2,
                               window.size.t = 2, window.size.s = 2,
                               make.pval = TRUE, n.sim = 1000,
                               verbose = F)
    wv_list[[j]] <-assign(paste0("wv_",vars_key[j]), my.wc)
  }
  return(wv_list)
}
#
## === FIGURES ====
## -- -- -- -- -- -- -- -- -- -- -- -- --
## Customize Figures
## -- -- -- -- -- -- -- -- -- -- -- -- --
#
## Fig. 1 - MAP ####
map_aucfz <- function(df, plygn = NULL, scale = NULL, clabs = NULL, ...){
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
      geom = "text", x = -57, y = -35.1, label = "Rio de la Plata", 
      fontface = "plain", family = "Times", color = "grey22", size = 4.5,
      angle = -40) +
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
      axis.text.x = element_text(vjust = -0.5),
      axis.text.y = element_text(hjust = -0.5),
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
## Fig. 2 - TIME SERIES: MONTHLY PEAKS ####
monthly_peaks.plot <- function(df){
  if (!is.data.frame(df))
    stop("'df' must be a dataframe containing landings and dates")
  startT <- df$date_num[1] 
  endT <- df$date_num[length(df$date_num)]
  start.end <- c(startT, endT)
  n_year <- length(unique(df$Year))
  months_n = rep(c("J", "F", "M", "A", "M", "J",
                   "J", "A", "S", "O", "N", "D"),
                 n_year)
  season_c = rep(c(1, 2, 3, 4), each = 3, n_year)
  p <- ggplot(df, aes(x = date_num, y = Landing)) +
    geom_line(alpha = 1) +
    scale_x_date(limits = c(startT,endT+48), date_breaks = "36 month",
      labels = date_format("%Y"), expand = c(0, 0)) +
    geom_point(shape = 15, fill = "green", color = "white", size = 4.5) +
    geom_point(pch = months_n, size = 3.5, color = season_c) +
    annotate("text",x = as.Date(c(
        "1999-01-01", "2005-01-01", "2011-01-01", "2017-01-01")),
        y = max(df$Landing) + max(df$Landing)/7,
        label = c("Summer", "Autumn", "Winter", "Spring"),
        color = c(1, 2, 3, 4), fontface = 2) +
    theme_tufte() +
    labs(y = "Landings (t)", x = "Dates") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain")) +
    scale_y_continuous(breaks = round(seq(
      0, round_any(max(df$Landing), 1000, f = ceiling),
      by = (round_any(max(df$Landing), 1000, f = ceiling))/10),1)) +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    theme(legend.title = element_blank(),
          legend.position = c(0.3, 0.97),
          legend.justification = c("left", "top"),
          legend.direction = "horizontal",
          legend.key = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, "cm"))
  # Vertical lines to highlight the separation between years
  vl <- list()
  januaries <- c(as.Date(df$date_num[seq(1, length(df$date_num), 12)]))
    for (k in 1:length(januaries)) {
      vl[[k]] <- annotate('segment', x = januaries[k], xend = januaries[k],
                          y = 0, yend = max(df$Landing) + max(df$Landing)/12,
                          linetype = "dotted", size = 0.75, alpha = 0.5,
                          color = "black")
      }
  pl <- p + vl
  pl
}
# 
## Fig. 3 & Figs. 4b, 5b, 6b - CORRELOGRAMS ####
correlograms.plot <- function(x, lags = NULL, main = NULL,
                              upperindex = NULL, ...){
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  if (is.null(lags)) {
    lags <- 36
  } else if ((lags %% 1 == 0) == FALSE || (lags < 0)) {
    stop('The argument "lags" must be a positive integer')
  }
  if (is.null(main) & is.null(upperindex)) {
    mytitle <- deparse(substitute(x))
    mysubt <- NULL
    sz <- 20
    fc <- 'plain'
    hz <- 0.5
    vt <- 1
  } else if (!is.null(main) & is.null(upperindex)) {
    if (!is.character(main))
      stop("'main' must be a character string")
    mytitle <- paste(main)
    mysubt <- NULL
    sz <- 20
    fc <- 'italic'
    hz <- 0.5
    vt <- 1
  } else if (is.null(main) & !is.null(upperindex)) {
    if (!is.character(upperindex))
      stop("'upperindex' must be a character string")
    mytitle <- NULL
    mysubt <- paste(upperindex)
    sz <- 20
    fc <- 'bold'
    hz <- 0.5
    vt <- -5
  } else if (!is.null(main) & !is.null(upperindex)) {
    if (!is.character(main))
      stop("'main' amd 'upperindex' must be a character string")
    mytitle <- paste(main)
    mysubt <- paste(upperindex)
    sz <- 20
    fc <- 'italic'
    hz <- 0.5
    vt <- -5
  }
  acf_v <- acf(x, lag.max = lags, type = "correlation", plot = FALSE)$acf
  p1 <- ggAcf(x, lag.max = lags) +
    theme_tufte() +
    scale_x_continuous(breaks = round(seq(0, lags, by = 2),1)) +
    scale_y_continuous(breaks = round(seq(
      round_any(min(acf_v), 1/10, f = floor),
      round_any(max(acf_v[2:length(acf_v)]), 1/10, f = ceiling),
      by = 0.1),1)) +
      theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    labs(title = mytitle, subtitle = mysubt) +
    theme(plot.title = element_text(
      size=sz, face = fc, hjust = hz, vjust = vt)) +
    theme(plot.subtitle = element_text(
      size=sz-4, face = "bold", hjust = hz-0.55, vjust = vt+10)) 
  pacf_v <- acf(x, lag.max = lags, type = "partial", plot = FALSE)$acf
  p2 <- ggAcf(x, lag.max = lags, type = 'partial') +
  labs(title = "") + theme_tufte() +
  scale_x_continuous(breaks = round(seq(0, lags, by = 2),1)) +
  scale_y_continuous(breaks = round(seq(
    round_any(min(pacf_v), 1/10, f = floor),
    round_any(max(pacf_v[2:length(pacf_v)]), 1/10, f = ceiling),
    by = 0.1),1)) +
  theme(axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour="black")) +
  theme(axis.text = element_text(color = "black", size = 12)) +
  labs(x="Lags") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="plain"))
pt <- p1 / p2
pt
}
#
## Figs. 4a, 5a, 6a - RESIDUALS_TEMPORAL SERIES ####
residuals_ts.plot <- function(x, upperindex = NULL, ...){
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(is.null(upperindex)){
    upperindex = NULL
    } else if (!is.character(upperindex))
    stop("'upperindex' must be a character string")
  startT = format(as.yearmon(time(x)[1]), format = "%Y")
  endT = format(as.yearmon(time(x)[length(x)]+1), format = "%Y")
  autoplot(ts(x, start = c(startT,1),
    frequency = 12),lty=1, lwd=.5) +
      labs(y = 'Residuals', x = 'Dates') + theme_tufte() +
      scale_x_continuous(breaks = round(seq(startT, endT, by = 2),1)) +
      scale_y_continuous(breaks = round(seq(
        round(min(stdres)), round(max(stdres)), by = 1),1)) +
      theme(axis.line = element_line(colour = "black"),
            axis.ticks = element_line(colour="black")) +
      theme(axis.text = element_text(color = "black", size = 12)) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="plain")) +
      ggtitle(upperindex) +
      theme(plot.title = element_text(size=16, face = "bold",
                                      hjust = -0.05, vjust = -4))
}
#
## Figs. 4c, 5c, 6c - RESIDUALS_HISTOGRAM ####
residuals_hist.plot <- function(x, pval = NULL, upperindex = NULL, ...){
  if (!is.numeric(x))
    stop("'x' must be a numeric vector")
  if(is.null(pval)){
    mytitle <- NULL
  } else if (is.numeric(pval)){
    mytitle <- paste("Ljung-Box p-value = ", round((pval), digits = 3))
  } else if (!is.numeric(pval)){
    stop("'p-value' from Ljung-Box test must be a numeric value")
  }
  if(is.null(upperindex)){
    upperindex = NULL
  } else if (!is.character(upperindex))
    stop("'upperindex' must be a character string")
  # Nested function of ggplot modified to plot the density curve in blue
  gghistogram_B <- function (x, add.normal = FALSE, add.kde = FALSE,
                             add.rug = TRUE, bins, boundary = 0){
    {
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")", 
             call. = FALSE)
      }
      else {
        if (missing(bins)) {
          bins <- min(500, grDevices::nclass.FD(na.exclude(x)))
        }
        data <- data.frame(x = as.numeric(c(x)))
        binwidth <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/bins
        p <- ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes(x), 
                                  data = data, color="gray23", fill = "gray23",
                                  binwidth = binwidth, boundary = boundary) + 
          ggplot2::xlab(deparse(substitute(x)))
        if (add.normal || add.kde) {
          xmin <- min(x, na.rm = TRUE)
          xmax <- max(x, na.rm = TRUE)
          if (add.kde) {
            h <- stats::bw.SJ(x)
            xmin <- xmin - 3 * h
            xmax <- xmax + 3 * h
          }
          if (add.normal) {
            xmean <- mean(x, na.rm = TRUE)
            xsd <- sd(x, na.rm = TRUE)
            xmin <- min(xmin, xmean - 3 * xsd)
            xmax <- max(xmax, xmean + 3 * xsd)
          }
          xgrid <- seq(xmin, xmax, length.out = 512)
          if (add.normal) {
            df <- data.frame(x = xgrid, y = length(x) * 
                               binwidth * stats::dnorm(xgrid, xmean, xsd))
            p <- p + ggplot2::geom_line(ggplot2::aes(df$x, 
                                                     df$y), 
                                        linetype="dashed", col = "blue")
          }
          if (add.kde) {
            kde <- stats::density(x, bw = h, from = xgrid[1], 
                                  to = xgrid[512], n = 512)
            p <- p + ggplot2::geom_line(ggplot2::aes(x = kde$x, 
                                                     y = length(x) * binwidth * kde$y), col = "#67a9ff")
          }
        }
        if (add.rug) {
          p <- p + ggplot2::geom_rug(ggplot2::aes(x))
        }
        return(p)
      }
    }
  }
  gghistogram_B(x, add.normal = TRUE, add.rug = TRUE) +
    theme_tufte() +
    scale_x_continuous(breaks = round(seq(
      round(min(x)), round(max(x)), by = 1),1)) +
    scale_y_continuous(breaks = round(seq(0, length(x)/6, by = 8),1)) +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    geom_text() +
    labs(x = 'Residuals', y = 'Frequency',
         title = mytitle,
         subtitle = upperindex) +
    theme(plot.title = element_text(hjust = 0.01, vjust = -15),
          plot.subtitle = element_text(size = 16, face = 'bold',
                                       hjust = -0.05, vjust = 1.8))
}
#
## Figs. 4d, 5d, 6d - RESIDUALS_p-VALUES ####
residuals_pv.plot <- function(df, upperindex = NULL,...){
  if (!is.data.frame(df))
    stop("'df' must be a dataframe got from LjungBoxTest")
  if (is.null(upperindex)) {
    mytitle = NULL
  } else if (!is.null(upperindex) & is.character(upperindex)) {
    mytitle <- upperindex
  } else if (!is.character(upperindex)){
  stop("'upperindex' must be a character string")
  }
  pv <- df[,3]
  ggplot(df, aes(x = m, y = pvalue)) +
    geom_point() +
    labs(title = "") + theme_tufte() +
    scale_x_continuous(breaks = round(seq(
      min(boxresult[,1]), max(boxresult[,1]), by = 2),1)) +
    scale_y_continuous(n.breaks = 5, limits = c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    labs(x="Lags", y = "p-value for Ljung-Box statistic") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    geom_hline(yintercept=0.05, linetype="dashed", color = "blue") +
    ggtitle(mytitle) +
    theme(plot.title = element_text(size=16, face = "bold",
                                    hjust = -0.05, vjust = 0.5))
}
## Figs. 7, 8, 9 - TIME SERIES: OBSERVED, FITTED AND FORECAST VALUES  ####
obs_fit_for.plot <- function(x, y, z, ufit, lfit, ufor, lfor,
                             main = NULL, line_col, shaded_col ){
  if (!is.numeric(x) || !is.numeric(y) || !is.numeric(z))
    stop("'x, y, z, ufit, lfit, ufor, lfor' must be a numeric vector")
  if (is.null(main)) {
    mytitle = NULL
  } else if (!is.character(main)){
    stop("'main' must be a character string")
  }
  startT = format(as.yearmon(time(x)[1]), format = "%Y")
  endT = format(as.yearmon(time(x)[length(x)]+1), format = "%Y")
  startTfor = format(as.yearmon(time(z)[1]), format = "%Y")
  autoplot(x, series=" Original", lty=1, lwd=1.3) +
    autolayer(y, lty = 1, lwd=0.8, series=paste0("Fitted")) +
    autolayer(ts(z, start = c(startTfor,1), frequency = 12),
              series="Forecast", lty = 1, lwd=.8) +
    geom_ribbon(data = y, aes(x = time(ufit),
                                      ymax = ufit, ymin = lfit),
                fill = shaded_col[1], alpha = 0.23) +
    geom_ribbon(data = z,
                aes(x = time(ufor),ymax = ufor, ymin = lfor),
                fill = shaded_col[2], alpha = 0.23) +
    scale_color_manual(values = line_col) +
    labs(y = 'Landings (t)', x = 'Dates') + theme_tufte() +
    scale_x_continuous(breaks = round(seq(startT, endT, by = 2),1)) +
    scale_y_continuous(
      breaks = round(seq(0,
                         round_any(max(up95fit), 1000, f = ceiling),
                         by = 1000),1)) +
    theme(axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour="black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    ggtitle(paste(main)) +
    theme(plot.title = element_text(size=20, face = "italic",
                                    hjust = 0.04, vjust = -4.5)) +
    theme(legend.title = element_blank(),
          legend.position = c(0.3, 0.97),
          legend.justification = c("left", "top"),
          legend.direction = "horizontal",
          legend.key = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1, "cm"))
}
## Fig. 10 - HEATMAP ====
heatmap.plot <- function(df, upperindex = NULL, ylabs = NULL, ...){
  if (!is.data.frame(df))
    stop("'df' must be a dataframe containing species and variable names,
         r-pearson values, and lags")
  if(is.null(upperindex)){
    upperindex = NULL
  } else if (!is.character(upperindex))
    stop("'upperindex' must be a character string")
  if (is.null(ylabs)){
    plt_ylabs = FALSE
    et = element_blank()
  } else if (!is.null(ylabs)){
    plt_ylabs = TRUE
    et = element_text(hjust = -0.4)
  }
  # Define the desired order of variables
  desired_order <- c("wmi", "chl", "sss", "kd", "sst", "river")
  # Convert variable to a factor with the desired order
  df$variable <- factor(df$variable, levels = desired_order)
  # 
  ggplot(data = df) +
    geom_tile(aes(x = lag, y = variable, fill = r)) +
    scale_fill_gradientn(
      colors = c("navyblue", "white", "white", "white", "red"),
      breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
      limits = c(-0.2, 0.2),
      na.value = "white"
    ) +
    labs(
      x = "Lag (months)", y = '',
      title = df$species[1],
      subtitle = upperindex
    ) +
    scale_x_continuous(
      breaks = seq(0, 15, by = 6),
      guide = "axis_minor",
      minor_breaks = seq(0, 15, by = 1),
      expand = c(0,0)
    ) +
    scale_y_discrete(
      expand = c(0,0),
      guide = "axis_minor"
    ) +
    guides(
      x.sec = "axis_minor"
    ) + # y.sec = "axis_minor"
    theme_tufte() +
    theme(
      panel.border = element_rect(fill = NA, color = "black", size = 1),
      panel.grid.major = element_blank()
    ) +
    theme(
      axis.ticks.length.y = unit(-0.25, "cm"),
      axis.ticks.length.x = unit(-0.25, "cm"),
      ggh4x.axis.ticks.length.minor = rel(0.5)
    ) +
    theme(
      axis.text = element_text(color = "black", size = 16),
      axis.text.y.left = et, #element_text(hjust = -0.4),
      axis.text.y.right = element_blank(),
      axis.text.x = element_text(vjust = -3),
      axis.text.x.top = element_blank(),
      axis.title.x = element_text(size = 15, margin = margin(t = 20))
    ) +
    theme(
      legend.position = "bottom", legend.key.width = unit(2.5, "cm"),
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    ) +
    theme(
      plot.title = element_text(size = 20, face = 'italic',
                                hjust = 0.5, vjust = -5),
      plot.subtitle = element_text(size = 23, face = 'bold',
                                   hjust = -0.01, vjust = 2)
    )
}
## Figs. 11, 12, 13 - WAVELET COHERENCY AND PHASE ANALYSES ====
getWavelets.plot = function(wvc, var_name, sp_name, ulab, cb = NULL) {
  if (!is.list(wvc))
    stop("'wvc' must be a list containing analysis coherence results")
  if (!is.character(var_name) || !is.character(sp_name) ||
      !is.character(ulab))
    stop("'var_name, sp_name, ulab' must be a character string")
  if (is.null(cb)){
    plt_lgnd = FALSE
  } else if (!is.null(cb)){
    plt_lgnd = TRUE
  }
  titlet <- bquote(bold(.(var_name))*bold(" and fishery landings of ")*
                     bolditalic(.(sp_name)))
  # wvc <-mysignal
  date_axis <- wvc[["series"]][["date"]]
  taxis <- seq(
    as.POSIXct(date_axis[1]),
    as.POSIXct(date_axis[length(date_axis)]), by = 'year')
  par(cex.lab = 1.5, cex.axis = 1.3)
  wc.image(
    wvc,
    which.image = "wc",
    n.levels = 250, color.key = "interval",
    maximum.level = 1, exponent = 1,
    siglvl.contour = 0.05, siglvl.arrow = 0.05, which.arrow.sig = "wc",
    main = title(titlet, cex.main = 1.7, line = 0.9),
    plot.legend = plt_lgnd,
    legend.params = list(lab = "Wavelet coherence levels",
                         label.digits = 1, mar = 5.5, cex.legend = 0.5),
    timelab = NULL, periodlab = "Period (months)",
    show.date = TRUE, date.format = "%F %T",
    spec.time.axis = list(at = taxis,
                          labels = year(as.Date(taxis, format = "%Y")),
                          las = 1),
    timetcl = -0.25,
    spec.period.axis = list(at = c(3, 6, 12, 24, 48, 84, 120)),
    periodtcl = -0.25,
    graphics.reset = F)
  mtext(paste0(ulab), side = 3, adj = -0.05, line = 0.5,
        cex = 2, font = 2)
  # return(wvc)
}
## Fig. 14 - AVERAGE WAVELET POWER SPECTRUM ENVIRONMENTAL ====
wavelet_power.plot <- function(wps_lines, color_patterns, lmts, lbls){
  if (!is.data.frame(wps_lines))
    stop("'wps_lines' must be a dataframe containing wavelet power spectrum
         results for each variable")
  ggplot() +
    geom_line(data = wps_lines, aes(x = period, y = wps, color = variable),
              linetype = 1,
              lwd = 1.5,
              alpha = 0.8) +
    scale_color_manual(values = colorBlindBlack8,
                       limits = lmts,
                       labels = lbls) +
    # coord_flip() +
    scale_x_continuous(name="Period (months)", 
                       breaks = c(0, 3, 6, 12, 24, 48, 84, 120),
                       labels = c(0, 3, 6, 12, 24, 48, 84, 120),
                       limits= c(2, 120)) +
    scale_y_continuous(name="Average wavelet power", 
                       breaks = c(seq(0, 1, 0.2)),
                       labels = c(seq(0, 1, 0.2)),
                       limits= c(0, 1)) +
    theme_tufte() +
    theme(
      legend.title = element_blank(),
      legend.justification = c(0, 1),
      legend.position = c(0.4, 0.9),
      legend.direction = "vertical",
      legend.spacing.x = unit(0.25, 'cm'),
      legend.spacing.y = unit(1, 'cm'),
      legend.key.size = unit(1, "cm"),
      legend.background = element_blank(),
      legend.text = element_text(
        colour = "black",
        size = 14,
        face = "italic")
    ) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"))
}
#                                       
## Fig. S1 - RELATIONSHIP CPUE - COMMERCIAL LANDINGS ####
cpue_landings.plot <- function(df){
  if (!is.data.frame(df))
    stop("'df' must be a dataframe containing species names, cpue and landings")
  sp <- unique(df$species)
  rs <- list(); ns <- list()
  for (i in 1:length(sp)) {
    subdf <- df[df$species == sp[i], ]
    ns[i] <- length(subdf$landings)
    mod_lm <- lm(cpue ~ landings, data = subdf)
    r <- sqrt(summary(mod_lm)[["r.squared"]])
    rs[i] <- round(r, 2)
    pvalue <- summary(mod_lm)$coefficients[, 4][2]
    if (pvalue < 0.001) {
      print(paste(sp[i], " -> p-value < 0.001", " / r = ", round(r, 2)))
    } else
      print(paste(sp[i], " -> p-value = ", round(pvalue, 3),
                  "/ r = ", round(r, 2)))
  }
 lm <-
    ggplot(df, aes(
      x = cpue,
      y = landings,
      color = species,
      fill = species
    )) +
    geom_point() +
    scale_x_continuous(breaks = round(seq(0, 2400, by = 200), 1)) +
    scale_y_continuous(breaks = round(seq(0, 110000, by = 10000), 1),
                       labels = comma) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.5) +
    scale_fill_brewer(palette = "Dark2") +
   annotate(
     geom = "text", x = 275, y = 1000,
     label = paste0("n= ", ns[1] ,"; r= ", rs[1]), 
     fontface = "plain", family = "Times", color = "black", size = 4) +
   annotate(
     geom = "text", x = 700, y = 13000,
     label = paste0("n= ", ns[2] ,"; r= ", rs[2]), 
     fontface = "plain", family = "Times", color = "black", size = 4) +
   annotate(
     geom = "text", x = 1500, y = 50000,
     label = paste0("n= ", ns[3] ,"; r= ", rs[3]), 
     fontface = "plain", family = "Times", color = "black", size = 4) +
    theme_tufte() +
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
    ) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    labs(y = "Landings (t)", x = "CPUE",
         caption = "Data source: DINARA's technical reports") +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "plain"))
  lm + scale_color_brewer(palette = "Dark2")
}
#
## Fig. S2 - BOX PLOTS - COMMERCIAL LANDINGS ####
bp_landings.plot <- function(df){
  out_mild <- {}; out_extr <- {}
  n <- {}
  i <- 0
  ux <- unique(df$Scientific_name)
  for(j in 1:length(ux)){
    i <- i+1
    idx <- which(toString(ux[j]) == df$Scientific_name)
    n[j] <- length(idx)
    sp <- df[idx,]
    # Get QUARTILES to detect outliers
    QT <- summary(sp$Landing)
    # Access elements in QUARTILE SUMMARY
    Min <- QT[[1]]; Q1 <- QT[[2]]; MD <- QT[[3]];
    AV <- QT[[4]]; Q3 <- QT[[5]]; Max <- QT[[6]];
    # Calculating IQR (interquartile range): Q3-Q1
    IQR = Q3-Q1
    # Minor & Extreme outliers fall outside (Q1 +- 1.5*IQR) & (Q3 +- 1.5*IQR)
    outliers_mild <- c(Q1 - 1.5*IQR , Q3 + 1.5*IQR)
    outliers_extr <- c(Q1 - 3*IQR , Q3 + 3*IQR)
    # Appending for each species
    out_mild <- rbind(out_mild, outliers_mild) 
    out_extr <- rbind(out_extr, outliers_extr)
    # Renaming rows according to each species common names
    rownames(out_mild)[i] <- ux[i]
    rownames(out_extr)[i] <- ux[i]
  }
  # Convert to DataFrame and renaming columns
  out_mild <- as.data.frame(out_mild); out_extr <- as.data.frame(out_extr); 
  colnames(out_mild) <- paste(c('lower','upper'))
  colnames(out_extr) <- paste(c('lower','upper'))
  upper_out.extr <- out_extr[2][[1]]
  # Kind of text
  bt <- element_text(face = "bold", color = "black") # Bold
  it <- element_text(face = "italic", color = "black") # Italic
  ggplot(df, aes(x= Scientific_name, y= Landing, fill = Scientific_name)) +
    geom_boxplot(outlier.colour="black", outlier.shape=16,
                 outlier.size=2, notch=FALSE, alpha = 1) +
    geom_segment(aes(y = upper_out.extr[1], yend = upper_out.extr[1],
                     x = 0.75, xend = 1.25),
                 linetype="dashed", colour = "black", size=.2) +
    annotate("text", x = 1, y = -250,
             label= paste("n = ",n[1]),
             family = "Times", size = 4.5, fontface = 1) +
    geom_segment(aes(y = upper_out.extr[2], yend = upper_out.extr[2],
                     x = 1.75, xend = 2.25),
                 linetype="dashed", colour = "black", size=.2) +
    annotate("text", x = 2, y = -250,
             label= paste("n = ",n[2]),
             family = "Times", size = 4.5, fontface = 1) +
    geom_segment(aes(y = upper_out.extr[3], yend = upper_out.extr[3],
                     x = 2.75, xend = 3.25),
                 linetype="dashed", colour = "black", size=.2) +
    annotate("text", x = 3, y = -250,
             label= paste("n = ",n[3]),
             family = "Times", size = 4.5, fontface = 1) +
    scale_fill_brewer(palette = 'Dark2') +
    theme(strip.background = element_blank(), strip.text.x = element_blank()) +
    theme(axis.text.x = element_blank()) +
    labs(x = "", y = "Landings (t)") +
    theme(axis.title = bt) +
    theme_tufte() +
    theme(legend.title=element_blank(),
          legend.justification=c(0,1), 
          legend.position= c(0.05, 0.95),
          legend.direction = "vertical",
          legend.spacing.x = unit(0.25, 'cm'),
          legend.spacing.y = unit(1, 'cm'),
          legend.key.size = unit(1, "cm"),
          legend.background = element_blank(),
          legend.text = element_text(colour = "black",
                                     size = 14,
                                     face = "italic")) +
    theme(axis.text.x=element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 12)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) 
}
#
## Figs. S3, S4 - DECOM.PLOT ====
decomp.plot <- function(x, main = NULL, type = NULL, ...){
  mylist <- c("x", "seasonal", "trend", "random", "figure", "type")
  if(!isTRUE(all(mylist %in% names(x), TRUE)))
    stop("x must be a list get from 'decompose' function")
  if(is.null(main))
    main <- paste("Decomposition of", x$type, "time series")
  startT = format(as.yearmon(time(x$x)[1]), format = "%Y")
  endT = format((as.yearmon(time(x$x)[length(x$x)])+1), format = "%Y")
  obs <- autoplot(x$x) +
    labs(title = main, y = 'Observed', x = '') +
    scale_x_continuous(breaks = round(seq(startT, endT, by = 1),1)) +
    scale_y_continuous(breaks = round(seq(
      0, round_any(max(na.omit(x$x)), 1000, f = ceiling), 
      by = (round_any(max(na.omit(x$x)), 1000, f = ceiling))/5),1),
      labels = comma) +
    theme_tufte() +
    theme(plot.title = element_text(size=20,
                                    hjust = 0.5)) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 10)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain"))
  tre <- autoplot(x$trend) +
    labs(y = 'Trend', x = '') +
    scale_x_continuous(breaks = round(seq(startT, endT, by = 1),1)) +
    scale_y_continuous(breaks = round(seq(
      0, round_any(max(na.omit(x$trend)), 1000, f = ceiling), 
      by = (round_any(max(na.omit(x$trend)), 1000, f = ceiling))/5),1),
      labels = comma) +
    theme_tufte() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 10)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "pt"))
  sea <- autoplot(x$seasonal) +
    labs(y = 'Seasonal', x = '') + theme_tufte() +
    scale_x_continuous(breaks = round(seq(startT, endT, by = 1),1)) +
    scale_y_continuous(breaks = round(seq(
      round_any(min(na.omit(x$seasonal)), 1000, f = floor),
      round_any(max(x$seasonal), 1000, f = ceiling),
      by = (round_any(max(na.omit(x$seasonal)), 1000, f = ceiling) -
              round_any(min(na.omit(x$seasonal)), 1000, f = floor))/5),1),
      labels = comma) +
    theme_tufte() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 10)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain"))
  ran <- autoplot(x$random) +
    labs(y = 'Random', x = 'Dates') +
    scale_x_continuous(breaks = round(seq(startT, endT, by = 1),1)) +
    scale_y_continuous(breaks = round(seq(
      round_any(min(na.omit(x$random)), 1000, f = floor),
      round_any(max(na.omit(x$random)), 1000, f = ceiling),
      by = (round_any(max(na.omit(x$random)), 1000, f = ceiling) -
              round_any(min(na.omit(x$random)), 1000, f = floor))/5),1),
      labels = comma) +
    theme_tufte() +
    theme(axis.line = element_line(colour = "black")) +
    theme(axis.text = element_text(color = "black", size = 10)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="plain")) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "pt"))
  par(mfrow = c(1, 1)) 
  p <- obs/tre/sea/ran 
  p
}
#
