#ifndef HMM_GRAPHICS_HPP_
#define HMM_GRAPHICS_HPP_

#include <RInside.h>
#include <Rcpp.h>

namespace hmm
{
  struct graphics
  {
    static RInside R;

    static void open_window()
    {
      Rcpp::Function x11("x11");
      x11();
    }

    static void wait() { R.parseEval("while(names(dev.cur()) != 'null device') Sys.sleep(1)"); }

    static void plot(const float *const x, const float *const y, const std::size_t &size,
                     const char *x_label, const char *y_label)
    {
      Rcpp::Function      plot("plot");
      Rcpp::NumericVector x_rcpp;
      Rcpp::NumericVector y_rcpp;
      for (std::size_t i = 0; i < size; ++i)
      {
        x_rcpp.push_back(x[i]);
        y_rcpp.push_back(y[i]);
      }
      plot(
          Rcpp::Named("x")    = x_rcpp,
          Rcpp::Named("y")    = y_rcpp,
          Rcpp::Named("type") = "l",
          Rcpp::Named("xlab") = x_label,
          Rcpp::Named("ylab") = y_label);
    }
  };

  // WARN: This should be in a source file to avoid new instantiation every time this file is included
  RInside graphics::R;
} // namespace hmm

#endif // HMM_GRAPHICS_HPP_