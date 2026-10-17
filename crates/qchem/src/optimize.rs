//! Geometry optimization by BFGS in Cartesian coordinates with a step-length cap and backtracking.

#[derive(Clone, Debug)]
pub struct OptimizerOptions {
    /// Convergence threshold on the largest gradient component (Hartree/bohr).
    pub max_force: f64,
    /// Convergence threshold on the RMS gradient (Hartree/bohr).
    pub rms_force: f64,
    /// Largest displacement of any atom in one step (bohr).
    pub max_step: f64,
    /// Initial Hessian guess (Hartree/bohr^2) times the identity.
    pub hessian_guess: f64,
}

impl Default for OptimizerOptions {
    fn default() -> Self {
        OptimizerOptions {
            max_force: 4.5e-4,
            rms_force: 3.0e-4,
            max_step: 0.3,
            hessian_guess: 0.4,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum OptStatus {
    /// Evaluate the energy and gradient at these coordinates next.
    Continue(Vec<f64>),
    Converged,
}

/// Step-wise BFGS driver: feed it energy and gradient at the coordinates it asks for.
pub struct Bfgs {
    pub options: OptimizerOptions,
    /// Inverse Hessian approximation.
    hinv: Vec<f64>,
    /// Last accepted point.
    best: Option<(Vec<f64>, f64, Vec<f64>)>,
    /// Step taken from `best` to the point being evaluated.
    last_step: Option<Vec<f64>>,
    pub steps: usize,
}

pub fn force_stats(g: &[f64]) -> (f64, f64) {
    let max = g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let rms = (g.iter().map(|v| v * v).sum::<f64>() / g.len().max(1) as f64).sqrt();
    (max, rms)
}

impl Bfgs {
    pub fn new(ncoord: usize, options: OptimizerOptions) -> Self {
        let mut hinv = vec![0.0; ncoord * ncoord];
        for i in 0..ncoord {
            hinv[i * ncoord + i] = 1.0 / options.hessian_guess;
        }
        Bfgs {
            options,
            hinv,
            best: None,
            last_step: None,
            steps: 0,
        }
    }

    pub fn is_converged(&self, g: &[f64]) -> bool {
        let (max, rms) = force_stats(g);
        max < self.options.max_force && rms < self.options.rms_force
    }

    /// Report energy and gradient at `x` (the point returned by the previous call, or the start).
    pub fn step(&mut self, x: &[f64], energy: f64, g: &[f64]) -> OptStatus {
        let n = x.len();
        self.steps += 1;
        if let (Some((bx, be, bg)), Some(step)) = (&self.best, &self.last_step) {
            if energy > be + 1e-7 {
                // Uphill: halve the step from the last good point and keep the Hessian.
                let half: Vec<f64> = step.iter().map(|v| 0.5 * v).collect();
                if half.iter().map(|v| v.abs()).fold(0.0, f64::max) < 1e-5 {
                    return OptStatus::Converged;
                }
                let next: Vec<f64> = bx.iter().zip(&half).map(|(a, b)| a + b).collect();
                let _ = bg;
                self.last_step = Some(half);
                return OptStatus::Continue(next);
            }
            // BFGS update of the inverse Hessian with s = x - bx, y = g - bg.
            let s: Vec<f64> = x.iter().zip(bx).map(|(a, b)| a - b).collect();
            let y: Vec<f64> = g.iter().zip(bg).map(|(a, b)| a - b).collect();
            let sy: f64 = s.iter().zip(&y).map(|(a, b)| a * b).sum();
            if sy > 1e-10 {
                let rho = 1.0 / sy;
                let hy: Vec<f64> = (0..n)
                    .map(|i| (0..n).map(|j| self.hinv[i * n + j] * y[j]).sum())
                    .collect();
                let yhy: f64 = y.iter().zip(&hy).map(|(a, b)| a * b).sum();
                for i in 0..n {
                    for j in 0..n {
                        self.hinv[i * n + j] += (1.0 + rho * yhy) * rho * s[i] * s[j]
                            - rho * (hy[i] * s[j] + s[i] * hy[j]);
                    }
                }
            }
        }
        self.best = Some((x.to_vec(), energy, g.to_vec()));
        if self.is_converged(g) {
            return OptStatus::Converged;
        }
        let mut p: Vec<f64> = (0..n)
            .map(|i| -(0..n).map(|j| self.hinv[i * n + j] * g[j]).sum::<f64>())
            .collect();
        // Fall back to steepest descent if the quasi-Newton step goes uphill.
        if p.iter().zip(g).map(|(a, b)| a * b).sum::<f64>() >= 0.0 {
            p = g.iter().map(|v| -v / self.options.hessian_guess).collect();
        }
        let longest = (0..n / 3)
            .map(|a| (p[3 * a].powi(2) + p[3 * a + 1].powi(2) + p[3 * a + 2].powi(2)).sqrt())
            .fold(0.0, f64::max);
        if longest > self.options.max_step {
            let f = self.options.max_step / longest;
            p.iter_mut().for_each(|v| *v *= f);
        }
        let next = x.iter().zip(&p).map(|(a, b)| a + b).collect();
        self.last_step = Some(p);
        OptStatus::Continue(next)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 2D Rosenbrock-like quadratic valley padded to 3 coordinates.
    fn f(x: &[f64]) -> (f64, Vec<f64>) {
        let (a, b, c) = (x[0] - 1.0, x[1] + 0.5, x[2]);
        let e = 2.0 * a * a + 0.5 * b * b + 0.3 * a * b + c * c;
        (e, vec![4.0 * a + 0.3 * b, b + 0.3 * a, 2.0 * c])
    }

    #[test]
    fn converges_on_a_quadratic() {
        let mut opt = Bfgs::new(
            3,
            OptimizerOptions {
                max_force: 1e-8,
                rms_force: 1e-8,
                ..Default::default()
            },
        );
        let mut x = vec![-1.0, 2.0, 0.7];
        for _ in 0..200 {
            let (e, g) = f(&x);
            match opt.step(&x, e, &g) {
                OptStatus::Continue(next) => x = next,
                OptStatus::Converged => break,
            }
        }
        assert!(
            (x[0] - 1.0).abs() < 1e-7 && (x[1] + 0.5).abs() < 1e-7 && x[2].abs() < 1e-7,
            "{x:?}"
        );
        assert!(opt.steps < 60);
    }
}
