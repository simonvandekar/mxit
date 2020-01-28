#' Calculates density of image intensity.
#'
#' @param dat list of imageDensity values
#' @keywords plot density of image intensity
#' @return Returns the density plot of image intensity for supplied `data`
#' @importFrom ggplot2 ggplot
#' @importFrom viridis viridis
#' @export
plotDensity = function(dat){
  if(!is.null(attr(dat,"logged"))){
    axis.txt = "Log Image Intensity"
  } else{
    axis.txt = "Image Intensity"
  }
  if(length(dat) > 1 & is.null(attr(dat,"class"))){
    #combine data into one nice ggplot set
    plot_dat = data.frame()
    for(i in 1:length(dat)){
      p = data.frame(cbind(dat[[i]]$x,
                dat[[i]]$y,
                attr(dat[[i]],"bias_correct")))
      colnames(p) = c("x",
                      "y",
                      "bias_correct")
      p$name = attr(dat[[i]],"name")
      p$bias_correct = factor(ifelse(p$bias_correct==1,
                                     "Yes",
                                     "No"),
                              levels = c("Yes","No"))
      plot_dat = rbind(plot_dat, p)
    }
  } else{
    plot_dat = data.frame(cbind(dat$x,
                     dat$y,
                     attr(dat,"bias_correct")))
    
    colnames(plot_dat) = c("x",
                           "y",
                           "bias_correct")
    plot_dat$name = attr(dat, "name")
    plot_dat$bias_correct = factor(ifelse(plot_dat$bias_correct==1,
                                          "Yes",
                                          "No"),
                                   levels = c("Yes","No"))
  }
  
  ggplot(plot_dat) +
    geom_line(aes(x=x,y=y,color=bias_correct),size=1.2) +
    theme_minimal()+
    scale_x_continuous(axis.txt)+
    scale_y_continuous("Density") +
    facet_wrap(~name) +
    scale_color_manual("Bias Corrected",
                       values = viridis(length(dat)))
}