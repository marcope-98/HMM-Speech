#ifndef HMM_IO_GRAPHICS_HPP_
#define HMM_IO_GRAPHICS_HPP_
#include <RInside.h>
#include <Rcpp.h>
#include <cmath> // TODO: remove this later

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
                     const char *x_label = "", const char *y_label = "")
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

    static void matrix(const float *const x, const float *const y, const float *const z,
                       const std::size_t &rows, const std::size_t &cols,
                       const char *xlab = "", const char *ylab = "")
    {
      Rcpp::Function      image("image");
      Rcpp::Function      hcl("hcl.colors");
      Rcpp::NumericVector x_out;
      Rcpp::NumericVector y_out;
      for (std::size_t i = 0; i < cols; ++i)
        x_out.push_back(x[i]);
      for (std::size_t i = 0; i < rows; ++i)
        y_out.push_back(y[i]);
      Rcpp::NumericMatrix z_out(cols, rows, z);
      image(Rcpp::Named("x")    = x_out,
            Rcpp::Named("y")    = y_out,
            Rcpp::Named("z")    = z_out,
            Rcpp::Named("xlab") = xlab,
            Rcpp::Named("ylab") = ylab,
            Rcpp::Named("col")  = hcl(
                Rcpp::Named("n")       = 100,
                Rcpp::Named("palette") = "inferno",
                Rcpp::Named("rev")     = false));
    }
  };

} // namespace hmm

#endif // HMM_GRAPHICS_HPP_