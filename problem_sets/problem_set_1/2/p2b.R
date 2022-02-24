## script to generate plot for PS1 2b
library(ggplot2)

dat <- read.csv("L2_h.csv", header = FALSE)

colnames(dat) = c("h", "x_norm", "y_norm")
dat$x_norm <- log(dat$x_norm)
dat$y_norm <- log(dat$y_norm)
dat$h <- log(dat$h)

out <- ggplot(dat, aes(x = h)) +
    geom_point(aes(y = x_norm), color = "red") +
    geom_line(aes(y = x_norm), color = "red") +
    geom_point(aes(y = y_norm), color = "black") + 
    geom_line(aes(y = y_norm), color = "black") +
    xlab("h") +
    ylab("L2 norm (log scale)") +
    scale_y_continuous(expand = c(0,0), limits = c(NA, NA)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, NA)) +
    theme_bw()

is.vector(out)

filepath <- "p2b.png"
ggsave(filepath, out)


 



