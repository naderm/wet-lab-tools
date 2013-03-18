#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import sys
from scipy.odr.odrpack import Model, RealData, ODR
from numpy import *
from pylab import *

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
  T = data[:,0] + 273
  signal = data[:,1]
  error = data[:,2]
  return T, signal, error

def _gibbs_free_energy(dH, C_p, T_m, t):
  return dH * (1 - t / T_m) - C_p * ((T_m - t) + t * log(t / T_m))

def _k(dH, C_p, T_m, t):
  R = 8.314
  return exp(_gibbs_free_energy(dH, C_p, T_m, t) / (R * t))

def _alpha(dH, C_p, T_m, t):
  k_calc = _k(dH, C_p, T_m, t)
  return k_calc / (1 + k_calc)

def _expected_signal(B, t):
  dH, C_p, T_m, sig_f, sig_u = B
  return _alpha(dH, C_p, T_m, t) * (sig_f - sig_u) + sig_u

def fit_cd_melt(T, sig, error):
  # Set up the guesses for the sig_f, sig_u, and T_m, all easy to find
  sig_f_guess, sig_u_guess = min(sig), max(sig)
  sig_mid = (sig_f_guess + sig_u_guess) / 2
  T_m_guess = min(enumerate(T), key = lambda x: abs(sig[x[0]] - sig_mid))[1]
  dH_guess = 0
  C_p_guess = 0
  guesses = [dH_guess, C_p_guess, T_m_guess, sig_f_guess, sig_u_guess]

  # Instead of using least squares, use orthogonal distance regression. This
  # lets us account for errors in the measurements of the data.
  # See: http://docs.scipy.org/doc/scipy/reference/odr.html
  linear = Model(_expected_signal)
  data = RealData(T, sig, sy = error)
  odr = ODR(data, linear, beta0 = guesses)
  output = odr.run()

  return output.beta, output.sd_beta, output.res_var

def main(args, show_graph = True):
  for arg in args:
    print("{}:".format(arg))

    with open(arg) as cd_input:
      T, sig, error = read_cd_data(cd_input)
      p, p_sd, res_var = fit_cd_melt(T, sig, error)
      dH, C_p, T_m = p[:3]
      dH_sd, C_p_sd, T_m_sd = p_sd[:3]

      # Covert the enthalpy from J/mol to kJ/mol
      dH, dH_sd = dH / 1000, dH_sd / 1000

      delta = u"\N{GREEK CAPITAL LETTER DELTA}".encode("utf-8")

      print("  {}H: {:.6} +/- {:.4} kJ/mol".format(delta, dH, dH_sd))
      print("  C_p: {:.6} +/- {:.4} J/mol/K".format(C_p, C_p_sd))
      print("  T_m: {:.6} +/- {:.4} K".format(T_m, T_m_sd))

      dg_t = 25 + 273
      dg = _gibbs_free_energy(dH, C_p, T_m, dg_t) / 1000
      print("  {}G @ {} K: {:.5} kJ/mol".format(delta, dg_t, dg))
      print("  Residual variance: {:.3}".format(res_var))

      if show_graph:
        temp = linspace(T.min(), T.max(), 100)
        plot(T, sig, "ro", temp, _expected_signal(p, temp), "k-")
        title("Temperature Melt of {}".format(arg))
        xlabel("Temperature (K)")
        ylabel("CD Signal (millidegrees)")
        show()
        savefig("{}.png".format(os.path.splitext(arg)[0]))

if __name__ == "__main__":
  main(sys.argv[1:])
