## script to generate plot for PS1 2b
library(ggplot2)

dat <- read.csv("L2_h.csv", header = FALSE)

colnames(dat) <- c("h", "x_norm", "y_norm")
dat$x_norm <- log(dat$x_norm)
dat$y_norm <- log(dat$y_norm)
dat$h <- log(dat$h)

ndat <- rbind(data.frame(h = dat$h, norm = dat$x_norm, norm_var = rep("x", times = nrow(dat))),
              data.frame(h = dat$h, norm = dat$y_norm, norm_var = rep("y", times = nrow(dat))))


out <- ggplot(ndat, aes(x = h, y = norm, linetype = norm_var)) +
    geom_point(size = 2) +
    geom_line(size = 1.5) +
    xlab("h (log scale)") +
    ylab("L2 norm (log scale)") +
    scale_y_continuous(expand = c(0,0), limits = c(-20, 5)) +
    scale_x_continuous(expand = c(0,0), limits = c(NA, -2)) +
    theme_bw()

is.vector(out)

filepath <- "p2b.png"
ggsave(filepath, out)


 



