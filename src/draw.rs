//! `points_dedup`, a drawing helper that qqman needs but rs-graphics doesn't provide.

use std::collections::HashSet;

use rs_graphics::prelude::*;

/// `points(x, y, pch = 20, col = col, cex = cex)`, skipping any point that
/// lands in the same cell of a `grid`-per-inch grid as one already drawn by
/// this call. Such a point adds (almost) no ink, and dropping it keeps plots
/// of millions of SNPs fast. See `plot::dedup_grid` for the cell size.
pub fn points_dedup<B: Backend>(
    dev: &mut Device<B>,
    xs: &[f64],
    ys: &[f64],
    col: Color,
    cex: f64,
    grid: f64,
) {
    let mut seen = HashSet::new();
    let (xs, ys): (Vec<f64>, Vec<f64>) = xs
        .iter()
        .zip(ys)
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .filter(|&(&x, &y)| {
            let p = dev.user_to_inches(x, y);
            seen.insert(((p.x * grid).round() as i64, (p.y * grid).round() as i64))
        })
        .unzip();
    dev.points(&xs, &ys).pch(20).col(col).cex(cex).draw();
}
