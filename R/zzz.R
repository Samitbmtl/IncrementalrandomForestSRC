.onAttach <- function(libname, pkgname) {
  irfsrc.version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                            fields="Version")
  packageStartupMessage(paste("\n",
                              pkgname,
                              irfsrc.version,
                              "\n",
                              "\n",
                              "Type irfsrc.news() to see new features, changes, and bug fixes.",
                              "\n",
                              "\n"))
}
