//! Spin-restricted exchange-correlation functionals.
//!
//! For each point: e = energy density per volume, vrho = de/drho, vsigma = de/dsigma with
//! sigma = |grad rho|^2. Matches libxc's lda_x + lda_c_vwn (VWN5) and gga_x_pbe + gga_c_pbe.

use std::f64::consts::{LN_2, PI};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Functional {
    /// Slater exchange + VWN5 correlation.
    Lda,
    /// PBE exchange + PBE correlation.
    Pbe,
}

impl Functional {
    pub fn is_gga(self) -> bool {
        matches!(self, Functional::Pbe)
    }

    /// Fraction of exact (Hartree-Fock) exchange.
    pub fn exact_exchange(self) -> f64 {
        0.0
    }
}

const RHO_MIN: f64 = 1e-14;

fn cx() -> f64 {
    -0.75 * (3.0 / PI).cbrt()
}

/// Slater exchange: (e, vrho).
pub fn slater(rho: f64) -> (f64, f64) {
    let r13 = rho.cbrt();
    (cx() * rho * r13, 4.0 / 3.0 * cx() * r13)
}

/// VWN5 paramagnetic correlation: (e, vrho).
pub fn vwn5(rho: f64) -> (f64, f64) {
    let (a, b, c, x0) = (0.031_090_7, 3.727_44, 12.935_2, -0.104_98);
    let rs = (3.0 / (4.0 * PI * rho)).cbrt();
    let x = rs.sqrt();
    let xx = x * x + b * x + c;
    let xx0 = x0 * x0 + b * x0 + c;
    let q = (4.0 * c - b * b).sqrt();
    let at = (q / (2.0 * x + b)).atan();
    let eps = a
        * ((x * x / xx).ln() + 2.0 * b / q * at
            - b * x0 / xx0 * (((x - x0) * (x - x0) / xx).ln() + 2.0 * (b + 2.0 * x0) / q * at));
    let deps = a
        * (2.0 / x
            - (2.0 * x + b) / xx
            - b / xx
            - b * x0 / xx0 * (2.0 / (x - x0) - (2.0 * x + b) / xx - (b + 2.0 * x0) / xx));
    (rho * eps, eps - x / 6.0 * deps)
}

const KAPPA: f64 = 0.804;
const MU: f64 = 0.219_514_972_764_517_1;

/// PBE exchange: (e, vrho, vsigma).
pub fn pbe_x(rho: f64, sigma: f64) -> (f64, f64, f64) {
    let s2_pref = 1.0 / (4.0 * (3.0 * PI * PI).powf(2.0 / 3.0));
    let (el, vl) = slater(rho);
    let r83 = rho.powf(8.0 / 3.0);
    let s2 = s2_pref * sigma / r83;
    let den = 1.0 + MU * s2 / KAPPA;
    let f = 1.0 + KAPPA - KAPPA / den;
    let df = MU / (den * den);
    (
        el * f,
        vl * f + el * df * (-8.0 / 3.0 * s2) / rho,
        el * df * s2_pref / r83,
    )
}

const BETA: f64 = 0.066_724_550_603_149_22;

/// PW92 (libxc PW_MOD) unpolarized correlation energy per particle and d/drs.
pub fn pw92(rs: f64) -> (f64, f64) {
    let (a, a1, b1, b2, b3, b4) = (0.031_090_7, 0.2137, 7.5957, 3.5876, 1.6382, 0.49294);
    let srs = rs.sqrt();
    let q0 = -2.0 * a * (1.0 + a1 * rs);
    let q1 = 2.0 * a * (b1 * srs + b2 * rs + b3 * rs * srs + b4 * rs * rs);
    let dq1 = a * (b1 / srs + 2.0 * b2 + 3.0 * b3 * srs + 4.0 * b4 * rs);
    let lg = (1.0 + 1.0 / q1).ln();
    (q0 * lg, -2.0 * a * a1 * lg - q0 * dq1 / (q1 * q1 + q1))
}

/// PBE correlation: (e, vrho, vsigma).
pub fn pbe_c(rho: f64, sigma: f64) -> (f64, f64, f64) {
    let gamma = (1.0 - LN_2) / (PI * PI);
    let t2_pref = PI / (16.0 * (3.0 * PI * PI).cbrt());
    let rs = (3.0 / (4.0 * PI * rho)).cbrt();
    let (eps, deps) = pw92(rs);
    let drs = -rs / (3.0 * rho);
    let ex = (-eps / gamma).exp();
    let a = BETA / gamma / (ex - 1.0);
    let da_deps = a * a * ex / BETA;
    let r73 = rho.powf(7.0 / 3.0);
    let y = t2_pref * sigma / r73;
    let nn = 1.0 + a * y;
    let dd = 1.0 + a * y + a * a * y * y;
    let g = BETA / gamma * y * nn / dd;
    let h = gamma * (1.0 + g).ln();
    let dh_dg = gamma / (1.0 + g);
    let dg_dy = BETA / gamma * (nn * dd + y * (a * dd - nn * (a + 2.0 * a * a * y))) / (dd * dd);
    let dg_da = -(BETA / gamma) * a * y * y * y * (2.0 + a * y) / (dd * dd);
    let dh_drho = dh_dg * (dg_da * da_deps * deps * drs + dg_dy * (-7.0 / 3.0 * y) / rho);
    (
        rho * (eps + h),
        eps + h + rho * (deps * drs + dh_drho),
        rho * dh_dg * dg_dy * t2_pref / r73,
    )
}

/// Evaluate a functional at one point: (e, vrho, vsigma). vsigma is 0 for LDA.
#[inline]
pub fn eval_point(xc: Functional, rho: f64, sigma: f64) -> (f64, f64, f64) {
    if rho < RHO_MIN {
        return (0.0, 0.0, 0.0);
    }
    match xc {
        Functional::Lda => {
            let (ex, vx) = slater(rho);
            let (ec, vc) = vwn5(rho);
            (ex + ec, vx + vc, 0.0)
        }
        Functional::Pbe => {
            let s = sigma.max(0.0);
            let (ex, vx, sx) = pbe_x(rho, s);
            let (ec, vc, sc) = pbe_c(rho, s);
            (ex + ec, vx + vc, sx + sc)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Potentials must be the derivatives of the energy density (central differences).
    #[test]
    fn potentials_are_derivatives_of_energy() {
        for xc in [Functional::Lda, Functional::Pbe] {
            for &(rho, sigma) in &[
                (1e-3, 1e-7),
                (0.05, 0.002),
                (0.7, 0.3),
                (12.0, 400.0),
                (300.0, 1e4),
            ] {
                let (_, vr, vs) = eval_point(xc, rho, sigma);
                let h = rho * 1e-5;
                let num_r = (eval_point(xc, rho + h, sigma).0 - eval_point(xc, rho - h, sigma).0)
                    / (2.0 * h);
                assert!(
                    (vr - num_r).abs() < 1e-6 * num_r.abs().max(1.0),
                    "{xc:?} vrho {vr} vs {num_r}"
                );
                if xc.is_gga() {
                    // e is nearly linear in sigma, so a large step avoids cancellation error.
                    let hs = sigma * 1e-3;
                    let num_s = (eval_point(xc, rho, sigma + hs).0
                        - eval_point(xc, rho, sigma - hs).0)
                        / (2.0 * hs);
                    assert!(
                        (vs - num_s).abs() < 2e-5 * num_s.abs(),
                        "{xc:?} vsigma {vs} vs {num_s}"
                    );
                }
            }
        }
    }

    /// (rho, sigma, e, vrho, vsigma) from libxc 'lda,vwn' and 'pbe' via PySCF 2.14 eval_xc.
    const LIBXC_LDA: [(f64, f64, f64, f64, f64); 4] = [
        (
            0.001,
            1e-07,
            -9.872067156718418e-05,
            -0.12819269645829556,
            0.0,
        ),
        (0.05, 0.002, -0.01602293898963633, -0.41823756163336406, 0.0),
        (0.7, 0.3, -0.5070853659079907, -0.9512160386446162, 0.0),
        (12.0, 400.0, -21.410008837154535, -2.3568666127091586, 0.0),
    ];
    const LIBXC_PBE: [(f64, f64, f64, f64, f64); 4] = [
        (
            0.001,
            1e-07,
            -9.881540390601881e-05,
            -0.12820208617569598,
            -0.7507018370001077,
        ),
        (
            0.05,
            0.002,
            -0.016059842393440122,
            -0.41551073125614657,
            -0.03707249296661849,
        ),
        (
            0.7,
            0.3,
            -0.5069132076014515,
            -0.9503654214471243,
            -0.0005126242297795866,
        ),
        (
            12.0,
            400.0,
            -21.40807480457341,
            -2.3550827559211243,
            -2.0732383123747602e-05,
        ),
    ];

    #[test]
    fn matches_libxc_reference_values() {
        for (xc, table) in [(Functional::Lda, LIBXC_LDA), (Functional::Pbe, LIBXC_PBE)] {
            for (rho, sigma, e, vr, vs) in table {
                let (e2, vr2, vs2) = eval_point(xc, rho, sigma);
                assert!(
                    (e2 - e).abs() < 1e-10 * e.abs().max(1.0),
                    "{xc:?} e {e2} vs {e}"
                );
                assert!(
                    (vr2 - vr).abs() < 1e-10 * vr.abs().max(1.0),
                    "{xc:?} vrho {vr2} vs {vr}"
                );
                assert!(
                    (vs2 - vs).abs() < 1e-10 * vs.abs().max(1.0),
                    "{xc:?} vsigma {vs2} vs {vs}"
                );
            }
        }
    }

    #[test]
    fn negligible_density_gives_zero() {
        assert_eq!(eval_point(Functional::Pbe, 1e-20, 1e-30), (0.0, 0.0, 0.0));
    }

    #[test]
    fn energy_densities_match_libxc_limits() {
        // Uniform electron gas: PBE reduces to LDA exchange and PW92 correlation at sigma = 0.
        let rho = 0.3;
        let (ex, _, _) = pbe_x(rho, 0.0);
        assert!((ex - slater(rho).0).abs() < 1e-15);
        let (ec, _, _) = pbe_c(rho, 0.0);
        let rs = (3.0 / (4.0 * PI * rho)).cbrt();
        assert!((ec - rho * pw92(rs).0).abs() < 1e-15);
        // PBE exchange enhancement is bounded by 1 + kappa.
        let (exl, _, _) = pbe_x(rho, 1e12);
        assert!(exl / slater(rho).0 <= 1.0 + KAPPA + 1e-12);
    }
}
