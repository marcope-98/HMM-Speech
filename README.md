# Hidden Markov Models applied to spoken digits recognition

C++ project on the use case of Hidden Markov Models to recognize spoken digits between different speakers.

## Dependencies

This library depends on R and Eigen. For now the dependencies are required especially for the Eigen3 library even though it is not currently in use.

As for the R/Rcpp/RInside dependencies they are required as of now. In the future a compile definition/option might be instantiated to turn on and off such dependency.

### R
The project uses R (version 4.1.2) to display plots and images, hence it needs a base installation of R and the Rcpp and RInside packages to be installed in the system.

To install R on Ubuntu, open a terminal and execute the command
```console
$ sudo apt update && sudo apt upgrade
$ sudo apt install r-base r-base-dev
```

Once R is installed in your system open a R session with sudo privileges and install the Rcpp and RInside packages as follows:
```console
$ sudo R
...
> install.packages(c("Rcpp", "RInside"))
```

## TODO:
- [x] WAV file format parser
- [x] Include Eigen for matrix operations
- [x] Fast Fourier Transform with a window
- [ ] Feature extraction pipeline
    - [ ] Preprocessing: Preemphesis, Voice Activation Detection
    - [ ] Frame Blocking and Windowing
    - [ ] Feature Extraction: Melcepstra coefficients
    - [ ] Post processing: Weight Function and Normalization
- [ ] Hidden Markov Models
- [ ] CMake compile option to turn on and off R dependency
- [ ] Find a way to deal with SIMD intrinsics: either with a compile option or a compile definition or by simply checking compatibility