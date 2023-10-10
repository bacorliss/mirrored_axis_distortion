



library(stringr)

# Convert fold change units to mirrored fold change
fc_to_mfc <- function(x) {x[x < 1 & !is.na(x)]<- -1/x[x < 1 & !is.na(x)]; return(x)}
# Mirrored fold change units to fold change
mfc_to_fc <- function(x) {x[x < 0 & !is.na(x)]<- -1/x[x < 0 & !is.na(x)]; return(x)}


#' Contraction transform
#' 
#' @description moves all points 1 unit closer to origin (points with (-1,1) 
#' become NA)
#' 
#' @param x numeric vector x
#' @returns numeric vector with transform applied.
#' 
contract1 <- function(x) {
  x1<-x
  x1[x <= -1 & !is.na(x)] <- x1[x <= -1 & !is.na(x)] + 1
  x1[x >=  1 & !is.na(x)] <- x1[x >=  1 & !is.na(x)] - 1
  x1[(x >  -1 & !is.na(x)) & (x <  1  & !is.na(x))] <- NaN
  return(x1)
}

#' Reverse contraction transform
#' 
#' @description moves all points 1 unit further from origin 
#' 
#' @param x numeric vector x
#' @returns numeric vector with transform applied.
#' 
rev_contract1 <- function(x) {
  x1<-x
  x1[x < 0 & !is.na(x)] <- x1[x < 0 & !is.na(x)] - 1
  x1[x >= 0 & !is.na(x)] <- x1[x >= 0 & !is.na(x)] + 1
  return(x1)
}



#' Mirrored relative change transform
#' 
#' @description converts elements of numeric vector x from units of relative 
#' change to units of mirrored relative change.
#' 
#' The equation for rc is:
#' rc = (y-x)/x
#' 
#' @param x numeric vector x
#' @param forward boolean for direction of conversion. TRUE: rc to mrc. 
#' FALSE: mrc to rc.
#' 
#' @return x modified vector with converted elements. 
#'
mirror_rc <- function(x, forward = TRUE) {

  # Define conversion functions
  rc_pos_to_neg <- function(x)  { - 1/(x + 1) + 1}
  rc_neg_to_pos <- function(x)  { - 1/(x - 1) - 1}
  
  if (forward) {
    x[x<0] = rc_pos_to_neg(x[x<0]); 
  } else if (!forward) {
    x[x<0] = rc_neg_to_pos(x[x<0]);
  } else {
    stop(sprintf("Argument for direction: must be boolean", direction))
  }
  
  return(x);
}


#' Mirrored fold change
#'
#' @description converts elements of numeric vector x from units of fold 
#' change (fc) to units of mirrored fold change (mfc).
#' 
#' The equation for fc is:
#' fc = y/x
#' 
#' fc is a measure of amount rather than change. An FC of 1 is no change.
#' 
#' @param x numeric vector x
#' @param forward boolean for direction of conversion. TRUE: fc to mfc. 
#' FALSE: mfc to fc.
#' 
#' @return x modified vector with converted elements.
#' 
mirror_fc <- function(x, forward = TRUE) {

  # Transforms to convert between fold change and mirrored fold change.
  fc_to_mfc <- function(x) {x[x < 1]<- -1/x[x < 1]; return(x)}
  mfc_to_fc <- function(x) {x[x < 0]<- -1/x[x < 0]; return(x)}
  
  # Cpnvert depending on direction of transform
  if (forward) {
    x1 <- fc_to_mfc(x)
  } else if (!forward) {
    x1 <- mfc_to_fc(x)
  }
  
  return(x1);
  
}


#' Reverses mirrored fold change transform on axis tick labels of ggplot2 
#' object to complete the mad-fc visualization.
#'
#' @description given a ggplot object where the specified axis has undergone a
#' mad-fc transform, reverses the axis transform on the axis labels (this
#' completes the visualization)
#' 
#' Essentially only the negative axis labels change, and are converted to
#' either decimal, power, or fraction format.
#' 
#' @param gg ggplot object that has one axis with mad-fc transform
#' @param ax axis that has the mad-fc trasnformed data (values either x or y)
#' @param num_format string specifying the axis tick label format, (values:
#' either "decimal", "fraction", or "power")
#' 
#' @return gg2 returns modified ggplot2 object with axis tick marks relabeled.
#' 
gg_revaxis_mfc<- function(gg, ax = "y", num_format = "decimal") {

  # browser()
  
  xlabs <- ggplot_build(gg)$layout$panel_params[[1]][[ax]]$get_labels()
  # Remove NAs (sometimes there are hidden empty ticks)
  xlabs <- xlabs[!is.na(xlabs)]
  
  x_breaks <- ggplot_build(gg)$layout$panel_params[[1]][[ax]]$get_breaks()
  x_breaks <- x_breaks[!is.na(x_breaks)]
  x_sigdigs = unname(sapply(xlabs, function(x) length(str_replace(x, "^[-0.]*",""))))
  
  new_xlabs_num <- rev_contract1(as.numeric(xlabs))
  new_xlabs_num <- as.numeric(xlabs)
  new_xlabs_num[new_xlabs_num==0] <- 1
  
  new_xbreaks <- contract1(x_breaks)
  new_xbreaks[is.nan(new_xbreaks)]<-0
  
  # Convert axis labels to specified format output
  if (num_format == "decimal") {
    new_xlabs_num[new_xlabs_num<1] <-  -1/new_xlabs_num[new_xlabs_num<1] 
    new_x_sigdigs = unname(sapply(as.character(new_xlabs_num), function(x) length(str_replace(x, "^[-0.]*",""))))
    new_xlabs = rep("",length(new_xlabs_num))
    for (n in seq_along(new_xlabs_num)) {
      new_xlabs[n] <- sprintf(paste0("%.",x_sigdigs[n]+ as.numeric(x_breaks<0)[n],"g"),
                              new_xlabs_num[n])
    }
    
  } else if (num_format == "power") {
    new_xlabs <- as.character(new_xlabs_num)
    new_xlabs[new_xlabs_num<1] <- paste(as.character(abs(new_xlabs_num[new_xlabs_num<1])),"^-1",sep="")
    

  } else if (num_format == "fraction") {
    new_xlabs <- as.character(new_xlabs_num)
    new_xlabs[new_xlabs_num<1] <- paste("1/", as.character(abs(new_xlabs_num[new_xlabs_num<1])),sep="")
    
  } else { stop("num_format only supports: decimal, power, or fraction") }
  
 
  
  if (ax=='x') {
    gg2<- gg + scale_x_continuous(labels = parse(text = new_xlabs), breaks = new_xbreaks)
  } else if (ax=='y'){
    gg2<- gg + scale_y_continuous(labels = parse(text = new_xlabs), breaks = new_xbreaks)
  } else (stop("gg_revaxis_mfc:only x and y axis supported for input argument"))
  
  
  return(gg2)
}




#' ggplot2 relabeling function for mad-fc transform
#'
#' @description 
#' 
#' @param 
#' @param 
#' 
#' @return returns new axis labels with mirrored fc transform reversed to complete 
#' mad-fc visualization.
#' 
mad_fc_labeller <- function(breaks, num_format = "fraction") {
  
  labels = rep("", length(breaks))
  labels[breaks == 0] <- as.character(1)
  labels[breaks > 0] <- as.character(breaks[breaks>0]+1)
  
  if (num_format == "decimal") {
    labels[breaks<0] <- as.character(1 /(breaks[breaks<0]-1))
  }  else if (num_format == "fraction") {
    labels[breaks<0] <- paste0("-1/",as.character(abs(breaks[breaks<0]-1)))
  }  else if (num_format == "power") {
    labels[breaks<0] <- paste0(as.character(abs(breaks[breaks<0]-1)),"^-1")
  } else {stop("num_format only supports: decimal, power, or fraction") }
  
  return(labels)
}
