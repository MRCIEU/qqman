//! `qqman::qq()`.

use anyhow::{Result, bail};
use rs_graphics::prelude::*;

use crate::draw;
use crate::plot::Common;

pub struct Prepared {
    pub expected: Vec<f64>,
    pub observed: Vec<f64>,
    pub xlim: (f64, f64),
    pub ylim: (f64, f64),
}

/// `ppoints(n)`.
fn ppoints(n: usize) -> impl Iterator<Item = f64> {
    let a = if n <= 10 { 3.0 / 8.0 } else { 0.5 };
    (1..=n).map(move |i| (i as f64 - a) / (n as f64 + 1.0 - 2.0 * a))
}

pub fn prepare(mut p: Vec<f64>, common: &Common) -> Result<Prepared> {
    p.retain(|&v| v.is_finite() && v > 0.0 && v < 1.0);
    if p.is_empty() {
        bail!("no p-values in (0, 1)");
    }
    p.sort_unstable_by(f64::total_cmp);
    let observed: Vec<f64> = p.iter().map(|v| -v.log10()).collect();
    let expected: Vec<f64> = ppoints(p.len()).map(|v| -v.log10()).collect();
    // Both are decreasing, so the maxima come first.
    Ok(Prepared {
        xlim: common.xlim.unwrap_or((0.0, expected[0])),
        ylim: common.ylim.unwrap_or((0.0, observed[0])),
        expected,
        observed,
    })
}

pub fn draw<B: Backend>(dev: &mut Device<B>, d: &Prepared, common: &Common) {
    // plot(e, o, pch = 20, ...): points, then axes, box and titles.
    dev.plot_new();
    dev.plot_window()
        .xlim(d.xlim.0, d.xlim.1)
        .ylim(d.ylim.0, d.ylim.1)
        .draw();
    draw::points_dedup(
        dev,
        &d.expected,
        &d.observed,
        Color::BLACK,
        common.cex,
        common.dedup_grid,
    );
    dev.axis().side(1).draw();
    dev.axis().side(2).draw();
    dev.box_().draw();
    dev.title()
        .with(&TitleOptions {
            main: common.main.as_deref().map(Label::from),
            xlab: Some(Label::math("Expected ~ ~-log[10](italic(p))")),
            ylab: Some(Label::math("Observed ~ ~-log[10](italic(p))")),
            ..TitleOptions::default()
        })
        .draw();
    // abline(0, 1, col = "red")
    dev.abline().a(0.0).b(1.0).col("red").draw();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ppoints_matches_r() {
        // ppoints(3); ppoints(12)[1]
        let small: Vec<f64> = ppoints(3).collect();
        let expect = [0.1923076923076923, 0.5, 0.8076923076923077];
        for (a, b) in small.iter().zip(expect) {
            assert!((a - b).abs() < 1e-15);
        }
        assert!((ppoints(12).next().unwrap() - 0.5 / 12.0).abs() < 1e-15);
    }

    #[test]
    fn filters_and_sorts() {
        let d = prepare(vec![0.5, f64::NAN, 1.0, 0.0, 0.01, 0.1], &Common::default()).unwrap();
        assert_eq!(d.observed, vec![2.0, 1.0, -(0.5f64.log10())]);
        assert_eq!(d.ylim, (0.0, 2.0));
    }
}
