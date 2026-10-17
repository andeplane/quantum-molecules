//! Boys function F_n(T) = ∫_0^1 t^(2n) exp(-T t^2) dt.
//!
//! Below `T_MAX` a table on a 0.05 grid plus a 6-term Taylor expansion gives F_nmax to about
//! 1e-14 and downward recursion fills the lower orders. Above `T_MAX` the asymptotic F_0 with
//! upward recursion is exact to double precision.

use std::sync::OnceLock;

pub const BOYS_NMAX: usize = 16;
const TAYLOR: usize = 6;
const STEP: f64 = 0.05;
const T_MAX: f64 = 40.0;
const NPTS: usize = (T_MAX / STEP) as usize + 1;
const NTAB: usize = BOYS_NMAX + TAYLOR + 1;
const INV_FACT: [f64; 7] = [
    1.0,
    1.0,
    0.5,
    1.0 / 6.0,
    1.0 / 24.0,
    1.0 / 120.0,
    1.0 / 720.0,
];

/// Reference evaluation by the convergent series
/// F_n(T) = exp(-T) Σ_k (2T)^k / ((2n+1)(2n+3)...(2n+2k+1)).
pub fn boys_series(n: usize, t: f64) -> f64 {
    let mut term = 1.0 / (2 * n + 1) as f64;
    let mut sum = term;
    for k in 1..2000 {
        term *= 2.0 * t / (2 * n + 2 * k + 1) as f64;
        sum += term;
        if term < sum * 1e-17 {
            break;
        }
    }
    (-t).exp() * sum
}

fn table() -> &'static [f64] {
    static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut tab = vec![0.0; NPTS * NTAB];
        for i in 0..NPTS {
            for n in 0..NTAB {
                tab[i * NTAB + n] = boys_series(n, i as f64 * STEP);
            }
        }
        tab
    })
}

/// Fill `out[0..=nmax]` with F_0(T)..F_nmax(T).
pub fn boys(nmax: usize, t: f64, out: &mut [f64]) {
    debug_assert!(nmax <= BOYS_NMAX);
    if t < T_MAX {
        let tab = table();
        let i = (t / STEP).round() as usize;
        let d = i as f64 * STEP - t;
        let base = i * NTAB + nmax;
        let mut f = 0.0;
        let mut dp = 1.0;
        for k in 0..TAYLOR {
            f += tab[base + k] * dp * INV_FACT[k];
            dp *= d;
        }
        out[nmax] = f;
        if nmax > 0 {
            let e = (-t).exp();
            for n in (0..nmax).rev() {
                out[n] = (2.0 * t * out[n + 1] + e) / (2 * n + 1) as f64;
            }
        }
    } else {
        out[0] = 0.5 * (std::f64::consts::PI / t).sqrt();
        if nmax > 0 {
            let e = (-t).exp();
            let inv2t = 1.0 / (2.0 * t);
            for n in 0..nmax {
                out[n + 1] = ((2 * n + 1) as f64 * out[n] - e) * inv2t;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_series_everywhere() {
        let mut out = [0.0; 20];
        for &t in &[
            0.0, 1e-7, 0.024, 0.3, 2.71, 11.2, 25.0, 39.97, 40.0, 40.1, 75.0, 110.0,
        ] {
            boys(BOYS_NMAX, t, &mut out);
            for n in 0..=BOYS_NMAX {
                let r = boys_series(n, t);
                assert!(
                    (out[n] - r).abs() <= 1e-12 * r,
                    "n={n} t={t} {} vs {r}",
                    out[n]
                );
            }
        }
    }

    #[test]
    fn known_values() {
        let mut out = [0.0; 4];
        boys(0, 0.0, &mut out);
        assert!((out[0] - 1.0).abs() < 1e-15);
        boys(3, 0.0, &mut out);
        assert!((out[3] - 1.0 / 7.0).abs() < 1e-15);
        boys(0, 300.0, &mut out);
        assert!((out[0] - 0.5 * (std::f64::consts::PI / 300.0).sqrt()).abs() < 1e-15);
        // F_0(T) = sqrt(pi/T) erf(sqrt(T)) / 2 with erf(1) = 0.8427007929497149
        boys(0, 1.0, &mut out);
        assert!(
            (out[0] - 0.5 * std::f64::consts::PI.sqrt() * 0.842_700_792_949_714_9).abs() < 1e-14
        );
    }
}
