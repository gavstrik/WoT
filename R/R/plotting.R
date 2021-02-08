add_alpha <- function(col, alpha=1){
  ## Add an alpha value to a color to make it transparent

  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

wbeta_colors <- function(wB, neutral = 0.4, cols = c("red","chartreuse3","blue")){
  # Make colors of weighted betas wB, with neutral at 0.4

  cols = colorRampPalette(cols, space="rgb", interpolate="linear")(200)
  midx = seq(0,2*neutral, length=length(cols))

  out = numeric(length(wB))
  for(i in 1:length(wB)){
    out[i] = cols[(midx-wB[i])^2 == min((midx-wB[i])^2)]
  }
  return(out)
}

