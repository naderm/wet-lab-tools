#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import sys
from numpy import *
from pylab import *
from scipy import optimize

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

def fit_cd_melt(T, sig, error):
  R = 8.314

  def k(h, cp, tm, t):
    return exp(gibbs_free_energy(h, cp, tm, t) / (R * t))

  def alpha(h, cp, tm, t):
    k_calc = k(h, cp, tm, t)
    return k_calc / (1 + k_calc)

  def fit_func(h, cp, tm, sig_f, sig_u, t):
    return alpha(h, cp, tm, t) * (sig_f - sig_u) + sig_u

  def err_func(p, t, signal):
    return fit_func(p[0], p[1], p[2], p[3], p[4], t) - signal

  sig_f_guess, sig_u_guess = min(sig), max(sig)
  sig_mid = (sig_f_guess + sig_u_guess) / 2
  tm_guess = min(enumerate(T), key = lambda x: abs(sig[x[0]] - sig_mid))[1]

  p0 = [0, 0, tm_guess, sig_f_guess, sig_u_guess]
  p2 = optimize.leastsq(err_func, p0[:], args = (T, sig))[0]
  return p2

def main(args):
  for arg in args:
    print("{}:".format(arg))
    with open(arg) as f:
      T, sig, error = read_cd_data(f)
      h, cp, tm = fit_cd_melt(T, sig, error)[:3]
      print("  dH: {:.5}".format(h))
      print("  Cp: {:.5}".format(cp))
      print("  Tm: {:.5}".format(tm))

      dg_t = 25
      dg = gibbs_free_energy(h, cp, tm, dg_t)
      degree_sign= u'\N{DEGREE SIGN}'.encode("utf-8")
      print("  dG @ {}{}C: {:.5}".format(dg_t, degree_sign, dg))



if __name__ == "__main__":
  main(sys.argv[1:])
