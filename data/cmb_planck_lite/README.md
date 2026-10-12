# Planck 2018 CMB power spectra (plik-lite, binned TT/TE/EE)

Files from the Planck 2018 data release (Planck Legacy Archive), as packaged in
the planck-lite-py repository (Prince & Dunkley 2019, arXiv:1909.05869;
github.com/heatherprince/planck-lite-py), unchanged:

* `planck2018_plik_lite/`: high-ℓ "plik-lite" binned spectra, foregrounds
  marginalised: 215 TT bins (ℓ 30–2508), 199 TE and 199 EE bins
  (ℓ 30–1996); `cl_cmb_plik_v22.dat` (bin centre, C_ℓ, σ),
  `c_matrix_plik_v22.dat` (613×613 covariance, Fortran binary), bin edges
  and weights.
* `planck2018_low_ell/`: two Gaussian low-ℓ TT bins (ℓ 2–29).

Cite: Planck Collaboration 2020, A&A 641, A5 ("Planck 2018 results. V. CMB
power spectra and likelihoods", arXiv:1907.12875) and Prince & Dunkley 2019.
The benchmark's own reader is `ubench/cmb_lite.py`, an independent
implementation of the same plik-lite logic. Planck calibration is fixed to 1.
Reionisation optical depth: Planck 2018 low-ℓ EE gives τ = 0.0544 ± 0.0073
(Planck 2018 VI, arXiv:1807.06209), used as one extra Gaussian data point.
