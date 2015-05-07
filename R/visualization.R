# visualization for TFBS analysis

# COLOR SCHEMES - from Paul Tol
# Qualitative color schemes by Paul Tol
tol1qualitative=c("#4477AA")
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

col_list = list(tol1qualitative,tol2qualitative, tol3qualitative, tol4qualitative, tol5qualitative, 
                tol6qualitative, tol7qualitative, tol8qualitative, tol9qualitative, tol10qualitative, tol11qualitative, tol12qualitative)

pal <- function(col, border = "light gray", ...){
  n <- length(col)
  plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, xlab = "", ylab = "", ...)
  rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
}

#pal(tol7qualitative)

basic_hist <- function(counts, xlab = "Counts", main = "Histogram of Counts", col = "cyan"){
  
  # makes a histogram of simulation outcomes and plots observed counts as dotted black line
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  breaks = seq(min(counts) - 0.5, max(counts) + 0.5,1)
  h = hist(counts, xlab=xlab, main=main, breaks = breaks, col=col, xaxt="n")
  axis(side=1, at = breaks+0.5, labels=breaks + 0.5)
}

coverage_hist <- function(coverages, xlab = "% of sequence bp covered by TFBS", main = "Histogram of Coverage", col = "cyan"){
  
  # makes a histogram of percentages binned by every 5%
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  breaks = seq(0,1,by=0.05)
  h = hist(coverages, xlab=xlab, main=main, breaks = breaks, col=col, xaxt="n")
  axis(side=1, at = seq(0,1,by=0.25), labels=seq(0,1,by=0.25))
}

two_panel_hist <- function(counts1, counts2, xlab = "Counts", mains = c("Category 1", "Category 2"), cols = tol2qualitative){
  
  # makes  a two color histogram comparing counts1, counts2, two different categories
  
  layout(cbind(1,2))
  breaks = seq(min(c(counts1, counts2)) - 0.5, max(c(counts1, counts2)) + 0.5, 1)

  
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h1 = hist(counts1, xlab=xlab, main=mains[1], breaks = breaks, col=cols[1], xaxt="n", yaxt="n", freq=FALSE )
  axis(side=1, at = breaks+0.5, labels=breaks + 0.5)
  axis(side=2, at = seq(0,0.6,0.1), labels=seq(0,0.6,0.1))
  
  
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h2 = hist(counts2, xlab=xlab, main=mains[2], breaks = breaks, col=cols[2], xaxt="n", yaxt="n", freq=FALSE)
  axis(side=1, at = breaks+0.5, labels=breaks + 0.5)
  axis(side=2, at = seq(0,0.6,0.1), labels=seq(0,0.6,0.1))
  
}

stacked_bar <- function(counts, labels, group_names, ylab = "Counts", beside = FALSE) {
  
  # plot a stacked bar graph
  par(mar=c(5,4,5,2) + 0.5)   # extra large bottom margin
  if (is.null(nrow(counts))){
    col = tol1qualitative
  } else {
    col = col_list[[nrow(counts)]]
  }
  
  bar <- barplot(counts, las=1, col = col, names.arg = labels, legend.text = group_names, ylab = ylab, beside = beside)
  return(bar)  
}

sim_hist <- function(counts, observed, xlab = "Simulation Outcomes", main = "Comparing Observation to Simulation", col = "cyan"){
  
  # makes a histogram of simulation outcomes and plots observed counts as dotted black line
  
  breaks = seq(min(min(counts), observed) - 0.5, max(max(counts), observed) + 0.5,1)
  h = hist(counts, xlab=xlab, main=main, breaks = breaks, col=col, xaxt="n")
  axis(side=1, at = breaks + 0.5, labels=breaks + 0.5)
  abline(v=observed, col="black", lty=3, lw=5)
}
