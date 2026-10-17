//! Small dense linear algebra on row-major `Vec<f64>` matrices.

/// (n x m) * (m x k)
pub fn matmul(a: &[f64], b: &[f64], n: usize, m: usize, k: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * k];
    for i in 0..n {
        for p in 0..m {
            let av = a[i * m + p];
            if av == 0.0 {
                continue;
            }
            let bo = &b[p * k..p * k + k];
            let co = &mut c[i * k..i * k + k];
            for j in 0..k {
                co[j] += av * bo[j];
            }
        }
    }
    c
}

pub fn transpose(a: &[f64], n: usize, m: usize) -> Vec<f64> {
    let mut t = vec![0.0; n * m];
    for i in 0..n {
        for j in 0..m {
            t[j * n + i] = a[i * m + j];
        }
    }
    t
}

/// X^T A X for A (n x n) and X (n x m).
pub fn congruence(x: &[f64], a: &[f64], n: usize, m: usize) -> Vec<f64> {
    matmul(&transpose(x, n, m), &matmul(a, x, n, n, m), m, n, m)
}

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

pub struct Eigen {
    /// Ascending eigenvalues.
    pub values: Vec<f64>,
    /// Eigenvectors as columns (row-major n x n).
    pub vectors: Vec<f64>,
}

/// Symmetric eigen-decomposition: Householder tridiagonalization + implicit QL (after JAMA).
pub fn eigh(a: &[f64], n: usize) -> Eigen {
    let mut v: Vec<Vec<f64>> = (0..n).map(|i| a[i * n..i * n + n].to_vec()).collect();
    let mut d = vec![0.0; n];
    let mut e = vec![0.0; n];
    if n > 0 {
        tred2(n, &mut v, &mut d, &mut e);
        tql2(n, &mut v, &mut d, &mut e);
    }
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        out[i * n..i * n + n].copy_from_slice(&v[i]);
    }
    Eigen {
        values: d,
        vectors: out,
    }
}

fn tred2(n: usize, v: &mut [Vec<f64>], d: &mut [f64], e: &mut [f64]) {
    for j in 0..n {
        d[j] = v[n - 1][j];
    }
    for i in (1..n).rev() {
        let mut scale = 0.0;
        let mut h = 0.0;
        for k in 0..i {
            scale += d[k].abs();
        }
        if scale == 0.0 {
            e[i] = d[i - 1];
            for j in 0..i {
                d[j] = v[i - 1][j];
                v[i][j] = 0.0;
                v[j][i] = 0.0;
            }
        } else {
            for k in 0..i {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            let mut f = d[i - 1];
            let mut g = h.sqrt();
            if f > 0.0 {
                g = -g;
            }
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for j in 0..i {
                e[j] = 0.0;
            }
            for j in 0..i {
                f = d[j];
                v[j][i] = f;
                g = e[j] + v[j][j] * f;
                for k in j + 1..i {
                    g += v[k][j] * d[k];
                    e[k] += v[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for j in 0..i {
                e[j] /= h;
                f += e[j] * d[j];
            }
            let hh = f / (h + h);
            for j in 0..i {
                e[j] -= hh * d[j];
            }
            for j in 0..i {
                f = d[j];
                g = e[j];
                for k in j..i {
                    v[k][j] -= f * e[k] + g * d[k];
                }
                d[j] = v[i - 1][j];
                v[i][j] = 0.0;
            }
        }
        d[i] = h;
    }
    for i in 0..n - 1 {
        v[n - 1][i] = v[i][i];
        v[i][i] = 1.0;
        let h = d[i + 1];
        if h != 0.0 {
            for k in 0..=i {
                d[k] = v[k][i + 1] / h;
            }
            for j in 0..=i {
                let mut g = 0.0;
                for k in 0..=i {
                    g += v[k][i + 1] * v[k][j];
                }
                for k in 0..=i {
                    v[k][j] -= g * d[k];
                }
            }
        }
        for k in 0..=i {
            v[k][i + 1] = 0.0;
        }
    }
    for j in 0..n {
        d[j] = v[n - 1][j];
        v[n - 1][j] = 0.0;
    }
    v[n - 1][n - 1] = 1.0;
    e[0] = 0.0;
}

fn tql2(n: usize, v: &mut [Vec<f64>], d: &mut [f64], e: &mut [f64]) {
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    let mut f = 0.0;
    let mut tst1: f64 = 0.0;
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            loop {
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for i in l + 2..n {
                    d[i] -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for row in v.iter_mut() {
                        h = row[i + 1];
                        row[i + 1] = s * row[i] + c * h;
                        row[i] = c * row[i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for j in i + 1..n {
            if d[j] < p {
                k = j;
                p = d[j];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for row in v.iter_mut() {
                row.swap(i, k);
            }
        }
    }
}

/// Solve A x = b by Gaussian elimination with partial pivoting. Returns None if singular.
pub fn solve(a: &[f64], b: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    for c in 0..n {
        let mut piv = c;
        for r in c + 1..n {
            if a[r * n + c].abs() > a[piv * n + c].abs() {
                piv = r;
            }
        }
        if piv != c {
            for k in 0..n {
                a.swap(c * n + k, piv * n + k);
            }
            b.swap(c, piv);
        }
        let d = a[c * n + c];
        if d.abs() < 1e-300 {
            return None;
        }
        for r in c + 1..n {
            let f = a[r * n + c] / d;
            if f == 0.0 {
                continue;
            }
            for k in c..n {
                a[r * n + k] -= f * a[c * n + k];
            }
            b[r] -= f * b[c];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_symmetric(n: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        let mut rnd = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = rnd();
                a[i * n + j] = v;
                a[j * n + i] = v;
            }
        }
        a
    }

    #[test]
    fn eigh_reconstructs_matrix() {
        for &n in &[1usize, 2, 5, 17, 40] {
            let a = random_symmetric(n, n as u64 + 3);
            let Eigen { values, vectors } = eigh(&a, n);
            for i in 1..n {
                assert!(values[i] >= values[i - 1]);
            }
            // V^T V = I and V diag(w) V^T = A
            let vtv = congruence(&vectors, &vec_identity(n), n, n);
            for i in 0..n {
                for j in 0..n {
                    let expect = if i == j { 1.0 } else { 0.0 };
                    assert!((vtv[i * n + j] - expect).abs() < 1e-12);
                    let mut s = 0.0;
                    for k in 0..n {
                        s += vectors[i * n + k] * values[k] * vectors[j * n + k];
                    }
                    assert!((s - a[i * n + j]).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn eigh_handles_degenerate_values() {
        let n = 4;
        let mut a = vec![0.0; n * n];
        for i in 0..n {
            a[i * n + i] = 2.0;
        }
        let e = eigh(&a, n);
        assert!(e.values.iter().all(|v| (v - 2.0).abs() < 1e-14));
    }

    #[test]
    fn solve_matches_known_system() {
        let a = [2.0, 1.0, -1.0, -3.0, -1.0, 2.0, -2.0, 1.0, 2.0];
        let b = [8.0, -11.0, -3.0];
        let x = solve(&a, &b, 3).unwrap();
        assert!(
            (x[0] - 2.0).abs() < 1e-12 && (x[1] - 3.0).abs() < 1e-12 && (x[2] + 1.0).abs() < 1e-12
        );
        assert!(solve(&[1.0, 2.0, 2.0, 4.0], &[1.0, 2.0], 2).is_none());
    }

    fn vec_identity(n: usize) -> Vec<f64> {
        let mut m = vec![0.0; n * n];
        for i in 0..n {
            m[i * n + i] = 1.0;
        }
        m
    }
}
