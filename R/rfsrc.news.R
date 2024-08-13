irfsrc.news <- function(...) {
  newsfile <- file.path(system.file(package="incrementalrandomForestSRC"), "NEWS")
  file.show(newsfile)
}
