
#' Measure the time and display original message.
#'
#'@param StopTime How long time you want to measure?
#'@param MESSAGE What message U want to display when stop timer.
#'@param disp logical. If FALSE, don't display the message.
#'@param units time units. auto", "secs", "mins", "hours", "days", "weeks". The default is "secs".
#'@param cex graphics parameter. expansion rate.
#'@param font graphics parameter. font number.
#'@param col graphics parameter. font colour.
#'@param family graohics parameter.
#'
#'@export
timer <- function(StopTime, MESSAGE = "STOP!!", disp = T, units = "secs",cex = 5,font = 1,col = 1,family=""){
  time <- 0
  if(disp == T){
    repeat {
      if(time == 0) START <- Sys.time()
      t0 <- Sys.time()
      DIFF <- 0
      while(DIFF == 0){
        t1 <- Sys.time()
        DIFF <- t1 - t0
        if(DIFF < 0.01) DIFF <- 0
      }
      time <- t1 - START
      if(time >= StopTime) break
      time <- round(as.numeric(time, units = units), digits = 3)
      cat(time, units,".\r")
    }

    cat("\n",MESSAGE)

    grDevices::dev.new()

    graphics::plot(NA, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    graphics::par(cex =cex, col = col, font = font, family = family)
    graphics::text(0.5,0.5,MESSAGE)
  }
}
