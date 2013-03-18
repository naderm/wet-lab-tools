#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import sys
from scipy.odr.odrpack import Model, RealData, ODR
from numpy import *
from pylab import *
# from scipy import optimize

def read_cd_data(cd_file):
  data = []

  for line in cd_file:
    if "$" in line:
      continue

    cells = line.split()
    if len(cells) != 5:
      continue

    try:
      data += [float(i) for i in cells],
    except ValueError:
      continue

  data = array(data)
  return data[:,0], data[:,1], data[:,2]

def gibbs_free_energy(h, cp, tm, t):
  return h * (1 - t / tm) - cp * ((tm - t) + t * log(t / tm))

def fit_cd_melt(T, sig, error, show_graph = True, name = "melt"):
  R = 8.314
  T += 273

  def k(h, cp, tm, t):
    return exp(gibbs_free_energy(h, cp, tm, t) / (R * t))

  def alpha(h, cp, tm, t):
    k_calc = k(h, cp, tm, t)
    return k_calc / (1 + k_calc)

  def fit_func(h, cp, tm, sig_f, sig_u, t):
    return alpha(h, cp, tm, t) * (sig_f - sig_u) + sig_u

  # def residuals(p, t, signal):
  #   return fit_func(p[0], p[1], p[2], p[3], p[4], t) - signal

  def fit_func_2(B, t):
    return fit_func(B[0], B[1], B[2], B[3], B[4], t)

  sig_f_guess, sig_u_guess = min(sig), max(sig)
  sig_mid = (sig_f_guess + sig_u_guess) / 2
  tm_guess = min(enumerate(T), key = lambda x: abs(sig[x[0]] - sig_mid))[1]

  guesses = [0, 0, tm_guess, sig_f_guess, sig_u_guess]
  # p = optimize.leastsq(residuals, guesses[:], args = (T, sig))[0]
  # res_var = sum(i ** 2 for i in residuals(p, T, sig)) / (len(sig) - len(guesses))
  # ss_tot = sum((sig - sig.mean()) ** 2)
  # res_var = 1 - ss_err / ss_tot

  # Instead of using least squares, use orthogonal distance regression. This
  # lets us account for errors in the measurements of the data.
  # See: http://docs.scipy.org/doc/scipy/reference/odr.html
  linear = Model(fit_func_2)
  data = RealData(T, sig, sy = error)
  odr = ODR(data, linear, beta0 = guesses)
  output = odr.run()

  if show_graph:
    temp = linspace(T.min(), T.max(), 100)
    plot(T, sig, "ro", temp, fit_func_2(output.beta, temp), "r-")
    title("Temperature melt of {}".format(name))
    xlabel("Temperature (K)")
    ylabel("CD Signal (millidegrees)")
    show()
    savefig("{}.png".format(os.path.splitext(name)[0]))

  return output.beta, output.sd_beta, output.res_var
  # return p, [0.] * len(p), res_var

def main(args):
  for arg in args:
    print("{}:".format(arg))

    with open(arg) as cd_input:
      T, sig, error = read_cd_data(cd_input)
      p, p_sd, res_var = fit_cd_melt(T, sig, error, name = arg)
      h, cp, tm = p[:3]
      h_sd, cp_sd, tm_sd = p_sd[:3]

      # Covert the enthalpy from J/mol to kJ/mol
      h, h_sd = h / 1000, h_sd / 1000

      delta = u"\N{GREEK CAPITAL LETTER DELTA}".encode("utf-8")

      print("  {}H: {:.6} +/- {:.4} kJ/mol".format(delta, h, h_sd))
      print("  Cp: {:.6} +/- {:.4} J/mol/K".format(cp, cp_sd))
      print("  Tm: {:.6} +/- {:.4} K".format(tm, tm_sd))

      dg_t = 25 + 273
      dg = gibbs_free_energy(h, cp, tm, dg_t) / 1000
      print("  {}G @ {} K: {:.5} kJ/mol".format(delta, dg_t, dg))
      print("  Residual variance: {:.3}".format(res_var))

if __name__ == "__main__":
  main(sys.argv[1:])
