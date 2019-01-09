# Auxiliary binary package rhlra2

## Installation instructions for linux/macOS

Requirements:
1. Build environment (working `GCC`/`clang` compiler together with `dev`-packages for `R` programming language (it depends on package management in your OS)).
2. `fftw` library (install it using package manager (`apt-get`, `brew`, etc.)).

Installation (two alternatives):
1. Using RStudio (New project -> Existing directory -> choose rhlra2 directory. Then click `Clean and rebuild` button (or `Build and reload`) for building. It is placed in top right position inside `Build` tab).
2. Using console: run `R CMD INSTALL .` inside `rhlra2` directory. 

## Installation instructions for Windows

Requirements:
1. [Rtools](https://cran.r-project.org/bin/windows/Rtools/) build environment. Please install it into `C:\RBuildTools` if possible. 
2. Download `fftw` build from Windows [from official site](http://www.fftw.org/install/windows.html). Choose a needed version of library (32bit/64bit) from "Precompiled FFTW 3.3.5 Windows DLLs" section. **It is required** to unpack `fftw` into `C:\RBuildTools\fftw` and **put this directory into PATH environment variable** (see instructions for Windows 10 [here](https://superuser.com/questions/949560/how-do-i-set-system-environment-variables-in-windows-10); **It is required** to do logout/login after editing PATH variable.
3. `devtools` package for `R` (`install.packages("devtools")`).

Installation (by steps):
1. Open RStudio, load `devtools` package (`library(devtools)`).
2. Run `find_rtools()` and make sure it returns `TRUE`.
3. Then install package in RStudio as usual: (New project -> Existing directory -> choose rhlra2 directory. Then click `Clean and rebuild` button (or `Build and reload`) for building. It is placed in top right position inside `Build` tab).