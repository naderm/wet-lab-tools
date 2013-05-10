#!bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import argparse
import os
import sys
from scipy.odr.odrpack import Model, RealData, ODR
from numpy import array, exp, log, gradient
from pylab import errorbar, linspace, title, xlabel, ylabel, show, savefig
from scipy import optimize

def read_cd_data(cd_file):
    """
    Reads a file containing CD data as output by Aviv's software.

    Returns an array of the temperatures, signals, and errors for the melt.
    """
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
    T = data[:,0] + 273.15
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

def fit_cd_melt(T, sig, error, odr = True):
    """
    Estimate the dH, C_p, T_m, sig_f, and sig_u associated with a temperature
    melt experiment. Uses scipy.odr to estimate error on the values.

    Returns an array of the estimatations, an array of the standard deviations
    associated with the estimates, and the residual variance of the fit.
    """
    # Set up the guesses for the sig_f, sig_u, and T_m, all easy to find
    sig_f_guess, sig_u_guess = min(sig), max(sig)
    sig_mid = (sig_f_guess + sig_u_guess) / 2
    T_m_guess = max(zip(T, gradient(T)), key = lambda x: x[1])[0]
    dH_guess = 0
    C_p_guess = 0
    guesses = [dH_guess, C_p_guess, T_m_guess, sig_f_guess, sig_u_guess]

    if odr:
        # Instead of using least squares, use orthogonal distance regression.
        # This allows us to account for errors in the measurements of the data.
        # See: http://docs.scipy.org/doc/scipy/reference/odr.html
        linear = Model(_expected_signal)
        data = RealData(T, sig, sy = error)
        odr = ODR(data, linear, beta0 = guesses)
        output = odr.run()
        output.pprint()
        return output.beta, output.sd_beta, output.res_var
    else:
        def err_func(p, t, signal):
            return _expected_signal(p, t) - signal

        p2 = optimize.leastsq(err_func, guesses[:], args = (T, sig))
        return p2[0], [None] * len(p2[0]), p2[1]

def _create_parser():
    parser = argparse.ArgumentParser(
    description = "Process data for CD temperature melts")
    parser.add_argument("file", nargs = "+",
                        help = "File containing CD data")
    parser.add_argument("-o", "--odr", action = "store_true",
                        help = "Use ODRs to calculate the fit, giving error "
                        "estimates")
    parser.add_argument("--room_temp", type = float, default = 273.15 + 25,
                        help = "Room temperature to use when calculating dG")

    return parser

def main(args, show_graph = True):
    parser = _create_parser()
    args = parser.parse_args(args)
    for f in args.file:
        print("{}:".format(f))

        with open(f) as cd_input:
            T, sig, error = read_cd_data(cd_input)
            p, p_sd, res_var = fit_cd_melt(T, sig, error, odr = args.odr)
            dH, C_p, T_m = p[:3]
            dH_sd, C_p_sd, T_m_sd = p_sd[:3]

            # Covert the enthalpy from J/mol to kJ/mol
            dH /= 1000

            if dH_sd is not None:
                dH_sd /= 1000

            delta = u"\N{GREEK CAPITAL LETTER DELTA}".encode("utf-8")

            print(("  {}H: {:.6}" + (" +/- {:.4}" if dH_sd is not None else "")
                   + " kJ/mol").format(delta, dH, dH_sd))
            print(("  C_p: {:.6}" + (" +/- {:.4}" if C_p_sd is not None else "")
                   + " J/mol/K").format(C_p, C_p_sd))
            print(("  T_m: {:.6}" + (" +/- {:.4}" if T_m_sd is not None else "")
                   + " K").format(T_m, T_m_sd))

            dg_t = args.room_temp
            dg = _gibbs_free_energy(dH, C_p, T_m, dg_t) / 1000
            print("  {}G @ {} K: {:.5} kJ/mol".format(delta, dg_t, dg))
            print("  Residual variance: {}".format(res_var))

            if show_graph:
                temp = linspace(T.min(), T.max(), 100)
                errorbar(T, sig, yerr = error, fmt = "ro")
                errorbar(temp, _expected_signal(p, temp), fmt = "k-")
                title("Temperature Melt of {}".format(f))
                xlabel("Temperature (K)")
                ylabel("CD Signal (millidegrees)")
                show()
                savefig("{}.png".format(os.path.splitext(f)[0]))

if __name__ == "__main__":
    main(sys.argv[1:])
