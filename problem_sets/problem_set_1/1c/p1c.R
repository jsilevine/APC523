## file to generate plot for p1b

library(ggplot2)

data <- read.csv("time.csv", header = FALSE)
colnames(data) <- c("num_cores", "runtime")

out <- ggplot(data, aes(x = num_cores, y = runtime)) +
    geom_point(size = 2) +
    geom_line() +
    theme_bw() +
    scale_y_continuous(expand = c(0,0), limits = c(0, 6)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 18)) +
  ylab("runtime (s)")

filepath <- "mpi_thread.png"
ggsave(filepath, out)
