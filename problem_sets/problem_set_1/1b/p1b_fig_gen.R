## file to generate plot for p1b

library(ggplot2)

data <- read.csv("time.csv")
colnames(data) <- c("num_cores", "runtime")

out <- ggplot(data, aes(x = num_cores, y = runtime) +
    geom_point(size = 2) +
    geom_line() +
    theme_bw +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))

filepath <- "omp_thread.png"
ggsave(out, filepath)