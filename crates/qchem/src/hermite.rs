//! McMurchie-Davidson Hermite expansion coefficients and Hermite Coulomb integrals.

use crate::boys::boys;

/// Hermite expansion coefficients E^{ij}_t for one cartesian direction, i <= imax, j <= jmax,
/// t <= i + j. E^{00}_0 = exp(-mu X_AB^2).
#[derive(Clone, Debug)]
pub struct HermiteE {
    pub imax: usize,
    pub jmax: usize,
    pub data: Vec<f64>,
}

impl HermiteE {
    pub fn new(imax: usize, jmax: usize, a: f64, b: f64, xab: f64) -> Self {
        let p = a + b;
        let mu = a * b / p;
        let xpa = -b / p * xab;
        let xpb = a / p * xab;
        let oo2p = 1.0 / (2.0 * p);
        let tdim = imax + jmax + 1;
        let jdim = jmax + 1;
        let mut e = vec![0.0; (imax + 1) * jdim * tdim];
        let idx = |i: usize, j: usize, t: usize| (i * jdim + j) * tdim + t;
        e[0] = (-mu * xab * xab).exp();
        for i in 0..imax {
            for t in 0..=i + 1 {
                let mut v = 0.0;
                if t > 0 {
                    v += oo2p * e[idx(i, 0, t - 1)];
                }
                if t <= i {
                    v += xpa * e[idx(i, 0, t)];
                }
                if t < i {
                    v += (t + 1) as f64 * e[idx(i, 0, t + 1)];
                }
                e[idx(i + 1, 0, t)] = v;
            }
        }
        for i in 0..=imax {
            for j in 0..jmax {
                for t in 0..=i + j + 1 {
                    let mut v = 0.0;
                    if t > 0 {
                        v += oo2p * e[idx(i, j, t - 1)];
                    }
                    if t <= i + j {
                        v += xpb * e[idx(i, j, t)];
                    }
                    if t < i + j {
                        v += (t + 1) as f64 * e[idx(i, j, t + 1)];
                    }
                    e[idx(i, j + 1, t)] = v;
                }
            }
        }
        HermiteE {
            imax,
            jmax,
            data: e,
        }
    }

    #[inline(always)]
    pub fn tdim(&self) -> usize {
        self.imax + self.jmax + 1
    }

    /// Offset of E^{ij}_0; E^{ij}_t is at offset + t.
    #[inline(always)]
    pub fn offset(&self, i: usize, j: usize) -> usize {
        (i * (self.jmax + 1) + j) * self.tdim()
    }

    /// E^{ij}_t with zero outside the valid range (negative indices allowed).
    #[inline(always)]
    pub fn get(&self, i: isize, j: isize, t: usize) -> f64 {
        if i < 0
            || j < 0
            || i as usize > self.imax
            || j as usize > self.jmax
            || t > (i + j) as usize
        {
            0.0
        } else {
            self.data[self.offset(i as usize, j as usize) + t]
        }
    }

    /// dE^{ij}_t / dA = 2a E^{i+1,j}_t - i E^{i-1,j}_t
    #[inline(always)]
    pub fn d_a(&self, i: usize, j: usize, t: usize, a: f64) -> f64 {
        let (i, j) = (i as isize, j as isize);
        2.0 * a * self.get(i + 1, j, t) - i as f64 * self.get(i - 1, j, t)
    }

    /// dE^{ij}_t / dB = 2b E^{i,j+1}_t - j E^{i,j-1}_t
    #[inline(always)]
    pub fn d_b(&self, i: usize, j: usize, t: usize, b: f64) -> f64 {
        let (i, j) = (i as isize, j as isize);
        2.0 * b * self.get(i, j + 1, t) - j as f64 * self.get(i, j - 1, t)
    }
}

/// Hermite Coulomb integrals R^0_{tuv}(alpha, X, Y, Z) for t + u + v <= L, stored in a cube
/// `out[(t * (L+1) + u) * (L+1) + v]`.
pub struct HermiteR {
    work: Vec<f64>,
    fbuf: [f64; 32],
    pub l: usize,
    pub out: Vec<f64>,
}

impl Default for HermiteR {
    fn default() -> Self {
        Self::new()
    }
}

impl HermiteR {
    pub fn new() -> Self {
        HermiteR {
            work: Vec::new(),
            fbuf: [0.0; 32],
            l: 0,
            out: Vec::new(),
        }
    }

    #[inline(always)]
    pub fn at(&self, t: usize, u: usize, v: usize) -> f64 {
        let l1 = self.l + 1;
        self.out[(t * l1 + u) * l1 + v]
    }

    pub fn compute(&mut self, l: usize, alpha: f64, x: f64, y: f64, z: f64) {
        let l1 = l + 1;
        let cube = l1 * l1 * l1;
        if self.work.len() < 2 * cube {
            self.work.resize(2 * cube, 0.0);
        }
        if self.out.len() < cube {
            self.out.resize(cube, 0.0);
        }
        self.l = l;
        boys(l, alpha * (x * x + y * y + z * z), &mut self.fbuf);
        let m2a = -2.0 * alpha;
        let mut pw = m2a.powi(l as i32);
        let (mut prev, mut cur) = (0usize, cube);
        let s1 = l1 * l1;
        for n in (0..=l).rev() {
            let lim = l - n;
            for t in 0..=lim {
                for u in 0..=lim - t {
                    for v in 0..=lim - t - u {
                        let k = (t * l1 + u) * l1 + v;
                        let r = if t > 0 {
                            let mut r = x * self.work[prev + k - s1];
                            if t > 1 {
                                r += (t - 1) as f64 * self.work[prev + k - 2 * s1];
                            }
                            r
                        } else if u > 0 {
                            let mut r = y * self.work[prev + k - l1];
                            if u > 1 {
                                r += (u - 1) as f64 * self.work[prev + k - 2 * l1];
                            }
                            r
                        } else if v > 0 {
                            let mut r = z * self.work[prev + k - 1];
                            if v > 1 {
                                r += (v - 1) as f64 * self.work[prev + k - 2];
                            }
                            r
                        } else {
                            pw * self.fbuf[n]
                        };
                        self.work[cur + k] = r;
                    }
                }
            }
            pw /= m2a;
            std::mem::swap(&mut prev, &mut cur);
        }
        self.out[..cube].copy_from_slice(&self.work[prev..prev + cube]);
    }
}
