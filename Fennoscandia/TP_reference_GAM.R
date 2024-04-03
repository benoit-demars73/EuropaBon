# GAM TP reference condition

# -------------------------------------------------------------------------
# Housekeeping

rm(list=ls()) # remove everything currently held in the R memory

graphics.off() # close all open graphic windows

# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(mgcv) 

setwd("~/Fennoscandia")

df <- read.table('Fennoscandia_1081lakes.txt',header=TRUE)
nrow(df) # 1081 lakes
head(df) 

str(df)
df$Country <- as.factor(df$Country)

df %>%
  filter(Ref_lake == 1) -> df.ref
nrow(df.ref) # 254
head(df.ref)

#########################################
####       Models and graphs         ####
#########################################

g.selec <- gam(log(TP) ~  log(Area_km2) +
                 s(log1p(Alk)) +
                 s(log(Colour)) +
       #          s(log(deg.sum)) +
                 s(X_3035,Y_3035) ,
               data = df.ref, 
               family = gaussian (link=identity), 
               method = "REML",
               trace = TRUE)
g.selec
summary(g.selec) # expl. deviance = 77.4%
AIC(g.selec)
# Choosing best model according to AIC and hand selection 
AIC (g0,g1,g2) 
# g0 full 270.9
# g1 full - x,y 324.8
# g2 full - deg-sum 280.0

# select g2 to avoid concurvity of deg.sum with x,y while keeping similar AIC

# model checking
coef(g.selec) # provide coefficients 
gam.check(g.selec) 
concurvity(g.selec, full = TRUE) # 0.56-0.65 -> OK
concurvity(g.selec, full = FALSE) # highest 0.59 between Alk and x,y -> OK


# gam residuals 
# also available: "pearson","scaled.pearson", "working" residuals
devi.res <- residuals(g.selec, type = "deviance")
boxplot(devi.res, horizontal = TRUE)
resp.res <- residuals(g.selec, type = "response")
boxplot(resp.res, horizontal = TRUE)


# Transform coordinate units to km
xy <- df[,5:6]
xy$X_km <- xy$X_3035/1000 # longitude in km
xy$Y_km <- xy$Y_3035/1000 # latitude in km


library(ncf)
sp.corr <- spline.correlog(
  xy$X_km, # units are km
  xy$Y_km, # units are km
  devi.res,
  w = NULL,
  df = NULL,
  type = "boot",
  resamp = 99,
  npoints = 50, # every 10 km with xmax = 500 km
  save = FALSE,
  filter = FALSE,
  fw = 0,
  max.it = 25,
  xmax = 500, # units are km
  latlon = FALSE,
  na.rm = FALSE,
  quiet = FALSE
)
sp.corr
summary(sp.corr)

autocorr.plot <- plot(sp.corr) # save as .tiff
# no spatial autocorrelation


#############################
#####    PARTIAL PLOT   #####
#############################

# Partial plots of model prediction 
library(psych) # help to calculate geometric.mean below
library(scales) # help with x axis log transformation
# custom margins
 par(mar=c(4,4,1,1))

# Alkalinity
sort(df.ref$Alk)
p.Alk <- seq(0.01, max(df.ref$Alk), by = 0.01) # minimum value set just above zero to avoid log issues

pgam <- predict(g.selec, newdata = data.frame(Alk = p.Alk,
                                              Area_km2 = geometric.mean(df.ref$Area_km2),
                                              Colour = geometric.mean(df.ref$Colour),
                                              X_3035 = mean(df.ref$X_3035),
                                              Y_3035 = mean(df.ref$Y_3035)),
                se.fit = T, type = "response")
p.Alk.plot   <- plot(pgam$fit ~ log(p.Alk), xaxt = "n", xlab = expression(paste("Alkalinity",
                                      " (mEq L"^-1, ")")), yaxt = "n", 
                                      ylab = "TP reference", ylim = log(c(1, 100)), type = "n") 
axis(1, at = log(c(0.01, 0.1, 1, 10)), 
     labels = c(0.01, 0.1, 1, 10), las = 2)
axis(2, at = log(c(1, 10, 100)), 
     labels = c(1, 10, 100), las = 2)
polygon(x= c(log(p.Alk), rev(log(p.Alk))), y = c(pgam$fit-2*pgam$se, rev(pgam$fit+2*pgam$se)),
        col = "grey", border = NA)
lines(pgam$fit ~ log(p.Alk), type = "l") 
#lines(pgam$fit-2*pgam$se ~ log(p.Alk), lty = 2)
#lines(pgam$fit+2*pgam$se ~ log(p.Alk), lty = 2)
Alk.rug <- df.ref[complete.cases(df.ref[ , c("Alk", "Area_km2", "Colour")]), "Alk"]
points(rep(0, length(Alk.rug)) ~ log(Alk.rug), pch = "+", col = "black")

p.Alk.plot   # save as tiff 400 x 300



# Colour
sort(df.ref$Colour)
p.Colour <- seq(min(df.ref$Colour), max(df.ref$Colour), by = 0.1) 

geometric.mean(df.ref$Alk+1)-1 # 0.23

pgam <- predict(g.selec, newdata = data.frame(Colour = p.Colour,
                                              Area_km2 = geometric.mean(df.ref$Area_km2),
                                              Alk = 0.23, # geometric mean
                                              X_3035 = mean(df.ref$X_3035),
                                              Y_3035 = mean(df.ref$Y_3035)),
                se.fit = T, type = "response")
p.Colour.plot   <- plot(pgam$fit ~ log(p.Colour), xaxt = "n", xlab = expression(paste("Colour",
                                                                                " (mg Pt L"^-1, ")")), yaxt = "n", 
                     ylab = "TP reference", ylim = log(c(1, 100)), type = "n") 
axis(1, at = log(c(0.1, 1, 10, 100, 1000)), 
     labels = c(0.1, 1, 10, 100, 1000), las = 2)
axis(2, at = log(c(1, 10, 100)), 
     labels = c(1, 10, 100), las = 2)
polygon(x= c(log(p.Colour), rev(log(p.Colour))), y = c(pgam$fit-2*pgam$se, rev(pgam$fit+2*pgam$se)),
        col = "grey", border = NA)
lines(pgam$fit ~ log(p.Colour), type = "l") 
#lines(pgam$fit-2*pgam$se ~ log(p.Colour), lty = 2)
#lines(pgam$fit+2*pgam$se ~ log(p.Colour), lty = 2)
Colour.rug <- df.ref[complete.cases(df.ref[ , c("Alk", "Area_km2", "Colour")]), "Colour"]
points(rep(0, length(Colour.rug)) ~ log(Colour.rug), pch = "+", col = "black")

p.Colour.plot   # save as tiff 400 x 300




# Area_km2
sort(df.ref$Area_km2)
p.Area_km2 <- seq(min(df.ref$Area_km2), max(df.ref$Area_km2), by = 0.1) 

geometric.mean(df.ref$Alk+1)-1 # 0.23

pgam <- predict(g.selec, newdata = data.frame(Area_km2 = p.Area_km2,
                                              Colour = geometric.mean(df.ref$Colour),
                                              Alk = 0.23, # geometric mean
                                              X_3035 = mean(df.ref$X_3035),
                                              Y_3035 = mean(df.ref$Y_3035)),
                se.fit = T, type = "response")
p.Area_km2.plot   <- plot(pgam$fit ~ log(p.Area_km2), xaxt = "n", xlab = expression(paste("Area_km2",
                                                                                      " (mg Pt L"^-1, ")")), yaxt = "n", 
                        ylab = "TP reference", ylim = log(c(1, 100)), type = "n") 
axis(1, at = log(c(0.1, 1, 10, 100, 1000)), 
     labels = c(0.1, 1, 10, 100, 1000), las = 2)
axis(2, at = log(c(1, 10, 100)), 
     labels = c(1, 10, 100), las = 2)
polygon(x= c(log(p.Area_km2), rev(log(p.Area_km2))), y = c(pgam$fit-2*pgam$se, rev(pgam$fit+2*pgam$se)),
        col = "grey", border = NA)
lines(pgam$fit ~ log(p.Area_km2), type = "l") 
#lines(pgam$fit-2*pgam$se ~ log(p.Area_km2), lty = 2)
#lines(pgam$fit+2*pgam$se ~ log(p.Area_km2), lty = 2)
Area_km2.rug <- df.ref[complete.cases(df.ref[ , c("Alk", "Area_km2", "Colour")]), "Area_km2"]
points(rep(0, length(Area_km2.rug)) ~ log(Area_km2.rug), pch = "+", col = "black")

p.Area_km2.plot   # save as tiff 400 x 300



# TP reference spatial patterns


#################################################
#### custom vis.gam to define colour palette ####
#################################################
# code from
# https://stackoverflow.com/questions/21750020/changing-the-colors-in-a-contour-plot-from-vis-gam-in-mgcv
library(grDevices)

jet.colors <-colorRampPalette(c("yellow", "red"))


#########################################
####  Modify the vis.gam function    ####
#########################################

# First modify vis.gam function to be able to backtransform the response variable
# Here we modify myvis.gam which inversed the colour palette
# https://stackoverflow.com/questions/15734789/changing-the-units-in-a-contour-plot-from-vis-gam-in-mgcv

# For usage just use color = 'jet'
# backtransformation of the response variable 

myvis.gam <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
                       col = NA, color = "heat", contour.col = NULL, se = -1, type = "link", 
                       plot.type = "persp", zlim = NULL, nCol = 50, ...) 
{
  fac.seq <- function(fac, n.grid) {
    fn <- length(levels(fac))
    gn <- n.grid
    if (fn > gn) 
      mf <- factor(levels(fac))[1:gn]
    else {
      ln <- floor(gn/fn)
      mf <- rep(levels(fac)[fn], gn)
      mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
      mf <- factor(mf, levels = levels(fac))
    }
    mf
  }
  dnm <- names(list(...))
  v.names <- names(x$var.summary)
  if (is.null(view)) {
    k <- 0
    view <- rep("", 2)
    for (i in 1:length(v.names)) {
      ok <- TRUE
      if (is.matrix(x$var.summary[[i]])) 
        ok <- FALSE
      else if (is.factor(x$var.summary[[i]])) {
        if (length(levels(x$var.summary[[i]])) <= 1) 
          ok <- FALSE
      }
      else {
        if (length(unique(x$var.summary[[i]])) == 1) 
          ok <- FALSE
      }
      if (ok) {
        k <- k + 1
        view[k] <- v.names[i]
      }
      if (k == 2) 
        break
    }
    if (k < 2) 
      stop("Model does not seem to have enough terms to do anything useful")
  }
  else {
    if (sum(view %in% v.names) != 2) 
      stop(paste(c("view variables must be one of", v.names), 
                 collapse = ", "))
    for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
                                 c("numeric", "factor"))) 
      stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
  }
  ok <- TRUE
  for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
    if (length(levels(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  else {
    if (length(unique(x$var.summary[[view[i]]])) <= 1) 
      ok <- FALSE
  }
  if (!ok) 
    stop(paste("View variables must contain more than one value. view = c(", 
               view[1], ",", view[2], ").", sep = ""))
  if (is.factor(x$var.summary[[view[1]]])) 
    m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
  else {
    r1 <- range(x$var.summary[[view[1]]])
    m1 <- seq(r1[1], r1[2], length = n.grid)
  }
  if (is.factor(x$var.summary[[view[2]]])) 
    m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
  else {
    r2 <- range(x$var.summary[[view[2]]])
    m2 <- seq(r2[1], r2[2], length = n.grid)
  }
  v1 <- rep(m1, n.grid)
  v2 <- rep(m2, rep(n.grid, n.grid))
  newd <- data.frame(matrix(0, n.grid * n.grid, 0))
  for (i in 1:length(x$var.summary)) {
    ma <- cond[[v.names[i]]]
    if (is.null(ma)) {
      ma <- x$var.summary[[i]]
      if (is.numeric(ma)) 
        ma <- ma[2]
    }
    if (is.matrix(x$var.summary[[i]])) 
      newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                          byrow = TRUE)
    else newd[[i]] <- rep(ma, n.grid * n.grid)
  }
  names(newd) <- v.names
  newd[[view[1]]] <- v1
  newd[[view[2]]] <- v2
  if (type == "link") 
    zlab <- paste("linear predictor")
  else if (type == "response") 
    zlab <- type
  else stop("type must be \"link\" or \"response\"")
  fv <- predict.gam(x, newdata = newd, se.fit = TRUE, type = type)
    fv$fit = exp(fv$fit)   #INSERT THIS LINE HERE WITH WHATEVER FUNCTION YOU WANT TO MODIFY THE FITTED VALUES BY
  z <- fv$fit
  if (too.far > 0) {
    ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
                             x$model[, view[2]], dist = too.far)
    fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
  }
  if (is.factor(m1)) {
    m1 <- as.numeric(m1)
    m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
  }
  if (is.factor(m2)) {
    m2 <- as.numeric(m2)
    m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
  }
  if (se <= 0) {
    old.warn <- options(warn = -1)
    av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
                 n.grid - 1)
    options(old.warn)
    max.z <- max(z, na.rm = TRUE)
    z[is.na(z)] <- max.z * 10000
    z <- matrix(z, n.grid, n.grid)
    surf.col <- t(av) %*% z %*% av
    surf.col[surf.col > max.z * 2] <- NA
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      min.z <- min(fv$fit, na.rm = TRUE)
      max.z <- max(fv$fit, na.rm = TRUE)
    }
    surf.col <- surf.col - min.z
    surf.col <- surf.col/(max.z - min.z)
    surf.col <- round(surf.col * nCol)
    con.col <- 1
    if (color == "heat") {
      pal <- heat.colors(nCol)
      con.col <- 3
    }
    else if (color == "topo") {
      pal <- topo.colors(nCol)
      con.col <- 2
    }
    else if (color == "cm") {
      pal <- cm.colors(nCol)
      con.col <- 1
    }
    else if (color == "terrain") {
      pal <- terrain.colors(nCol)
      con.col <- 2
    }
    else if (color == "gray" || color == "bw") {
      pal <- gray(seq(0.1, 0.9, length = nCol))
      con.col <- 1
    }
    ### colour customized here
    else if (color == 'jet') {
      pal <- jet.colors(nCol)
      con.col = 1
    }
    ####
    else stop("color scheme not recognised")
    if (is.null(contour.col)) 
      contour.col <- con.col
    surf.col[surf.col < 1] <- 1
    surf.col[surf.col > nCol] <- nCol
    if (is.na(col)) 
      col <- pal[as.array(surf.col)]
    z <- matrix(fv$fit, n.grid, n.grid)
    if (plot.type == "contour") {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                    sep = "")
      if (color != "bw") {
        txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
        txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                     ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                     sep = "")
        eval(parse(text = txt))
      }
      else {
        txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
    else {
      stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                    ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                    ifelse("main" %in% dnm, "", ",zlab=zlab"), ",...)", 
                    sep = "")
      if (color == "bw") {
        op <- par(bg = "white")
        txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                     stub, sep = "")
        eval(parse(text = txt))
        par(op)
      }
      else {
        txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                     stub, sep = "")
        eval(parse(text = txt))
      }
    }
  }
  else {
    if (color == "bw" || color == "gray") {
      subs <- paste("grey are +/-", se, "s.e.")
      lo.col <- "gray"
      hi.col <- "gray"
    }
    else {
      subs <- paste("red/green are +/-", se, "s.e.")
      lo.col <- "green"
      hi.col <- "red"
    }
    if (!is.null(zlim)) {
      if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
        stop("Something wrong with zlim")
      min.z <- zlim[1]
      max.z <- zlim[2]
    }
    else {
      z.max <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
      z.min <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
    }
    zlim <- c(z.min, z.max)
    z <- fv$fit - fv$se.fit * se
    z <- matrix(z, n.grid, n.grid)
    if (plot.type == "contour") 
      warning("sorry no option for contouring with errors: try plot.gam")
    stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                  ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                                                                         dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                                                                                                        "", ",sub=subs"), ",...)", sep = "")
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=lo.col"), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=\"black\""), stub, sep = "")
    eval(parse(text = txt))
    par(new = TRUE)
    z <- fv$fit + se * fv$se.fit
    z <- matrix(z, n.grid, n.grid)
    txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
                                                             dnm, "", ",border=hi.col"), stub, sep = "")
    eval(parse(text = txt))
  }
}


# For usage just use color = 'jet':

myvis.gam(g.selec,
          view = c("X_3035", "Y_3035"),
          plot.type = "contour",
          contour.col = "black",labcex=0.8,
          color = "jet",
          too.far = 0.04,
          type = "response",
          main="TP reference")
points(x=df.ref$X_3035, y=df.ref$Y_3035, type="p", pch=1, col="grey", bg=NA, cex=0.5)

# save as tiff 400 x 300


##########################################
#####    PREDICTIONS TP REFFERENCE   #####
##########################################


# TP reference conditions for the 1081 lakes
newdata = data.frame(Alk = df$Alk,
                     Area_km2 = df$Area_km2,
                     Colour = df$Colour,
                     X_3035 = df$X_3035,
                     Y_3035 = df$Y_3035)
head(newdata)
colSums(is.na(newdata)) # no NAs

pgam <- predict(g.selec, newdata,
                se.fit = T, type = "response")

# Add the fit and se to the new data frame
newdata$TP.fit <- exp(pgam$fit)
# 95% Wald confidence limits calculated on the log-odds scale and 
# backtransformed to response scale with plogis 
newdata$TP.se_lwr <- exp(pgam$fit + qnorm(0.025)* pgam$se.fit) # qnorm(0.025) = -1.96
newdata$TP.se_upr <- exp(pgam$fit + qnorm(0.975)* pgam$se.fit) # qnorm(0.975) = +1.96
newdata$TP.se_range <- newdata$TP.se_upr - newdata$TP.se_lwr

head(newdata)
newdata <- newdata %>% select(-(Alk:Y_3035))
df.ref.pred <- cbind(df, newdata)

head(df.ref.pred)
write.table(df.ref.pred,"~/TP reference/Fennoscandia_1081lakes_ref.txt",sep="\t",row.names=FALSE)

