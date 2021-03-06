## file to make figure for problem set 1.1
require("ggplot2")

## read data
dat <- read.csv("final_phi.csv", header = FALSE)

n <- ncol(dat)-1

## convert to long format and add x and y vals for plotting
dat <- matrix(as.matrix(dat[,1:n]), ncol = 1)
dat <- cbind(dat, matrix(c(rep(seq(1, n), each = n), rep(seq(1, n), times = n)), ncol = 2))
colnames(dat) <- c("phi", "x", "y")
dat <- as.data.frame(dat) ## change to data.frame

## generate plot
out <- ggplot(dat, aes(x = x, y = y, fill = phi)) +
  geom_raster() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  ##coord_quickmap() +
  theme_bw()

filepath <- "phi_128.png"
ggsave(filepath, out)

print(paste0("output saved to ", filepath))
