import subprocess
import sys
import unittest

import isitgr


def _sigma8_and_tt(**mg_params):
    pars = isitgr.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122, mnu=0.06, omk=0, tau=0.06, **mg_params)
    pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
    pars.set_for_lmax(1000, lens_potential_accuracy=0)
    pars.set_matter_power(redshifts=[0.0], kmax=2.0)
    results = isitgr.get_results(pars)
    tt = results.get_cmb_power_spectra(pars, CMB_unit="muK")["total"][100, 0]
    return results.get_sigma8_0(), tt


class ISiTGRTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sigma8_GR, cls.tt_GR = _sigma8_and_tt()

    def test_no_fast_math(self):
        # importing the library must not switch on flush-to-zero (e.g. from -ffast-math)
        code = "before = str(1e-323); import isitgr; assert str(1e-323) == before, str(1e-323)"
        subprocess.run([sys.executable, "-c", code], check=True)

    def test_gr(self):
        self.assertTrue(0.80 < self.sigma8_GR < 0.805, self.sigma8_GR)
        self.assertGreater(self.tt_GR, 0)

    def test_mg_gr_limits(self):
        for mg in (
            dict(MG_parameterization="muSigma", mu0=0.0, Sigma0=0.0),
            dict(MG_parameterization="mueta", E11=0.0, E22=0.0),
        ):
            with self.subTest(**mg):
                sigma8, tt = _sigma8_and_tt(**mg)
                self.assertAlmostEqual(sigma8, self.sigma8_GR, places=5)
                self.assertAlmostEqual(tt / self.tt_GR, 1, places=5)

    def test_mg_enhanced_growth(self):
        for mg in (
            dict(MG_parameterization="mueta", E11=1.0, E22=1.0),
            dict(MG_parameterization="muSigma", mu0=0.2, Sigma0=0.0),
        ):
            with self.subTest(**mg):
                sigma8, _ = _sigma8_and_tt(**mg)
                self.assertGreater(sigma8, self.sigma8_GR * 1.01)


if __name__ == "__main__":
    unittest.main()
