import math

ep = 1e-5

def get_s_of_pm(x: float, m0: float) -> float:
    x_sq = x * x
    x4 = x_sq * x_sq
    return (4.0 / (x4 * x)) * math.exp(-1.0 / x4) * m0

def get_s_of_js(x: float, m0: float, l: float) -> float:
    xm = math.pow(4.0 / 5.0, 0.25)
    xPM = xm + l * (x - xm);
    return l * get_s_of_pm(xPM, m0)


# ------------------------------------------------------------------
#  Main routine – translated from C++ to Python
# ------------------------------------------------------------------
def calc_freq_spectrum(
    jsGamma: float,
    spectrum_A: float,
    spectrum_B: float,
    spectrum_b: float,      # note the name in the original was `spectrum_BB`
) -> tuple[list[float], list[float]]:

    # ------------------------------------------------------------------
    #  Prepare the vectors (Python lists are dynamic; no reserve needed)
    # ------------------------------------------------------------------
    s_vec = []
    w_vec = []

    n_bins = 80

    gamma = -1.0
    use_pm = True

    # ------------------------------------------------------------------
    #  Gamma selection – copy of the C++ logic
    # ------------------------------------------------------------------
    if jsGamma > 1 + ep:
        gamma = jsGamma
        use_pm = False
        lgG = math.log10(gamma);
        _lambda = 1 + 0.261 * lgG - 0.0525 * lgG * lgG

    # ------------------------------------------------------------------
    #  Core calculation
    # ------------------------------------------------------------------
    upper_bound = 4
    x_delta = upper_bound / n_bins

    PM_m0 = 0.25 * spectrum_A / spectrum_B
    minS = 2e-4;

    # Fill the vectors
    for i in range(1, n_bins + 1):
        x = i * x_delta
        if use_pm:
            s_val = get_s_of_pm(x, PM_m0)
        else:
            s_val = get_s_of_js(x, PM_m0, _lambda)

        if s_val > minS:
            s_vec.append(s_val)
            w_vec.append(x / spectrum_b)   # Eq 16/18

    return s_vec, w_vec

def area(xVec, yVec):
    sum = 0
    for i in range(0, len(xVec) - 1):
        dx = xVec[i + 1] - xVec[i]
        y = (yVec[i + 1] + yVec[i]) / 2
        sum += (y * dx)
    return sum

pp1d = 7
swh = 2

jsGamma = 3.3
spec_B = math.pow(2 * math.pi / pp1d, 4) / 0.8

spec_A = 0.25 * spec_B * swh * swh
spec_b = math.pow(spec_B, -0.25)

s_valsPM, w_valsPM = calc_freq_spectrum(-1, spec_A, spec_B, spec_b)
s_valsJS, w_valsJS = calc_freq_spectrum(jsGamma, spec_A, spec_B, spec_b)

#sPM_max = max(s_valsPM)
#graphArea = area(w_valsPM, s_valsPM)

import matplotlib.pyplot as plt
import numpy as np

plt.scatter(w_valsPM, s_valsPM, label="PM")
plt.scatter(w_valsJS, s_valsJS, label=f"JS gamma {jsGamma}")

plt.xlabel("w")
plt.ylabel("S")
plt.legend()

plt.show()
