//! `qqman::manhattan()`.

use std::collections::HashSet;
use std::ops::Range;

use anyhow::{Result, bail};
use rs_graphics::prelude::*;

use crate::draw;
use crate::plot::{Common, resolve_color};

pub struct ManhattanInput {
    pub chr: Vec<f64>,
    pub bp: Vec<f64>,
    pub p: Vec<f64>,
    /// Only read when highlighting or annotating.
    pub snp: Option<Vec<String>>,
}

pub struct ManhattanOptions {
    pub col: Vec<ColorSpec>,
    pub chrlabs: Option<Vec<String>>,
    pub suggestiveline: Option<f64>,
    pub genomewideline: Option<f64>,
    pub highlight: Option<HashSet<String>>,
    pub logp: bool,
    pub annotate_pval: Option<f64>,
    pub annotate_top: bool,
    pub common: Common,
}

/// Everything `manhattan()` computes before drawing, in qqman's
/// `order(CHR, BP)` order.
pub struct Prepared {
    pub pos: Vec<f64>,
    pub logp: Vec<f64>,
    pub p: Vec<f64>,
    pub chr: Vec<f64>,
    pub snp: Option<Vec<String>>,
    /// Row ranges of each chromosome (qqman's `index`).
    pub groups: Vec<Range<usize>>,
    pub ticks: Vec<f64>,
    pub labs: Vec<String>,
    pub xlabel: String,
    pub xlim: (f64, f64),
    pub ylim: (f64, f64),
}

/// R's `as.character()` for a chromosome code.
fn format_chr(c: f64) -> String {
    if c.fract() == 0.0 && c.abs() < 1e15 {
        format!("{}", c as i64)
    } else {
        format!("{c}")
    }
}

/// `quantile(x, 0.5)` (type 7) of an already-sorted slice.
fn sorted_median(x: &[f64]) -> f64 {
    let h = (x.len() - 1) as f64 * 0.5;
    let lo = h.floor() as usize;
    x[lo] + (h - lo as f64) * (x[h.ceil() as usize] - x[lo])
}

pub fn prepare(input: ManhattanInput, opts: &ManhattanOptions) -> Result<Prepared> {
    let ManhattanInput { chr, bp, p, snp } = input;
    let n = chr.len();

    // qqman would fail on NA/Inf; drop those rows with a warning instead.
    let keep: Vec<usize> = (0..n)
        .filter(|&i| {
            chr[i].is_finite()
                && bp[i].is_finite()
                && p[i].is_finite()
                && (!opts.logp || p[i] > 0.0)
        })
        .collect();
    if keep.len() < n {
        eprintln!(
            "warning: dropped {} rows with missing, non-finite or zero values",
            n - keep.len()
        );
    }
    if keep.is_empty() {
        bail!("no rows with valid CHR, BP and P values");
    }

    // d[order(d$CHR, d$BP), ] - a stable sort, skipped if already sorted.
    let mut order = keep;
    let key = |i: usize| (chr[i], bp[i]);
    if !order.windows(2).all(|w| key(w[0]) <= key(w[1])) {
        order.sort_by(|&a, &b| key(a).partial_cmp(&key(b)).unwrap());
    }
    let chr: Vec<f64> = order.iter().map(|&i| chr[i]).collect();
    let mut bp: Vec<f64> = order.iter().map(|&i| bp[i]).collect();
    let p: Vec<f64> = order.iter().map(|&i| p[i]).collect();
    let snp = snp.map(|s| order.iter().map(|&i| s[i].clone()).collect::<Vec<_>>());
    let logp: Vec<f64> = if opts.logp {
        p.iter().map(|v| -v.log10()).collect()
    } else {
        p.clone()
    };

    let mut groups = Vec::new();
    let mut start = 0;
    for i in 1..=chr.len() {
        if i == chr.len() || chr[i] != chr[start] {
            groups.push(start..i);
            start = i;
        }
    }

    let (pos, ticks, labs, xlabel);
    if groups.len() == 1 {
        pos = bp;
        ticks = Vec::new();
        labs = Vec::new();
        xlabel = format!("Chromosome {} position", format_chr(chr[0]));
    } else {
        // qqman rebases every chromosome after the first to start at 1,
        // and offsets it by the summed (rebased) maxima of the earlier ones.
        // Rows are sorted by BP within each chromosome, so a group's max BP
        // is its last row.
        let mut lastbase = 0.0;
        let mut prev_max = bp[groups[0].end - 1];
        for g in &groups[1..] {
            lastbase += prev_max;
            let min = bp[g.start];
            for v in &mut bp[g.clone()] {
                *v = *v - min + 1.0 + lastbase;
            }
            prev_max = bp[g.end - 1] - lastbase;
        }
        pos = bp;
        ticks = groups
            .iter()
            .map(|g| sorted_median(&pos[g.clone()]))
            .collect();
        labs = groups.iter().map(|g| format_chr(chr[g.start])).collect();
        xlabel = "Chromosome".to_string();
    }

    let labs = match &opts.chrlabs {
        Some(custom) if !labs.is_empty() && custom.len() == labs.len() => custom.clone(),
        Some(_) if !labs.is_empty() => {
            eprintln!(
                "warning: You're trying to specify chromosome labels but the number of labels != number of chromosomes."
            );
            labs
        }
        _ => labs,
    };

    let max_pos = pos.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let max_logp = logp.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok(Prepared {
        xlim: opts
            .common
            .xlim
            .unwrap_or(((max_pos * -0.03).floor(), (max_pos * 1.03).ceil())),
        ylim: opts.common.ylim.unwrap_or((0.0, max_logp.ceil())),
        pos,
        logp,
        p,
        chr,
        snp,
        groups,
        ticks,
        labs,
        xlabel,
    })
}

/// Fails where R would stop: a page too small for the margins.
pub fn draw<B: Backend>(
    dev: &mut Device<B>,
    d: &Prepared,
    opts: &ManhattanOptions,
) -> Result<(), rs_graphics::Error> {
    let c = &opts.common;

    // plot(NA, xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", las = 1, ...)
    dev.par().xaxs(AxisStyle::Exact).yaxs(AxisStyle::Exact).draw();
    dev.try_plot_new()?;
    dev.plot_window()
        .xlim(d.xlim.0, d.xlim.1)
        .ylim(d.ylim.0, d.ylim.1)
        .draw();
    dev.axis().side(2).las(1).draw();
    dev.title()
        .with(&TitleOptions {
            main: c.main.as_deref().map(Label::from),
            xlab: Some(d.xlabel.as_str().into()),
            ylab: Some(Label::math("-log[10](italic(p))")),
            ..TitleOptions::default()
        })
        .draw();

    // qqman forwards `cex` to axis(1, ...) too, but C_axis() sizes labels
    // from cex.axis alone, so it has no effect there.
    if d.groups.len() == 1 {
        dev.axis().side(1).draw();
    } else {
        dev.axis().side(1).at(&d.ticks).labels(&d.labs).draw();
    }

    let palette = dev.current_par().palette.clone();
    let cols: Vec<Color> = opts
        .col
        .iter()
        .map(|&s| resolve_color(s, &palette))
        .collect();
    for (gi, g) in d.groups.iter().enumerate() {
        let col = cols[gi % cols.len()];
        draw::points_dedup(
            dev,
            &d.pos[g.clone()],
            &d.logp[g.clone()],
            col,
            c.cex,
            c.dedup_grid,
        );
    }

    if let Some(h) = opts.suggestiveline.filter(|&h| h != 0.0) {
        dev.abline().h(h).col("blue").draw();
    }
    if let Some(h) = opts.genomewideline.filter(|&h| h != 0.0) {
        dev.abline().h(h).col("red").draw();
    }

    if let (Some(highlight), Some(snp)) = (&opts.highlight, &d.snp) {
        if highlight.iter().any(|h| !snp.contains(h)) {
            eprintln!("warning: You're trying to highlight SNPs that don't exist in your results.");
        }
        let (xs, ys): (Vec<f64>, Vec<f64>) = (0..snp.len())
            .filter(|&i| highlight.contains(&snp[i]))
            .map(|i| (d.pos[i], d.logp[i]))
            .unzip();
        draw::points_dedup(
            dev,
            &xs,
            &ys,
            Color::parse("green3").unwrap(),
            c.cex,
            c.dedup_grid,
        );
    }

    if let Some(threshold) = opts.annotate_pval {
        annotate(dev, d, opts, threshold);
    }
    Ok(())
}

/// The `annotatePval` block, including `calibrate::textxy()`.
fn annotate<B: Backend>(
    dev: &mut Device<B>,
    d: &Prepared,
    opts: &ManhattanOptions,
    threshold: f64,
) {
    let Some(snp) = &d.snp else { return };
    let hit = |i: usize| {
        if opts.logp {
            d.p[i] <= threshold
        } else {
            d.p[i] >= threshold
        }
    };
    let mut hits: Vec<usize> = (0..d.p.len()).filter(|&i| hit(i)).collect();
    let cex = if opts.annotate_top {
        // Best hit per chromosome, chromosomes in order of their best P.
        hits.sort_by(|&a, &b| d.p[a].partial_cmp(&d.p[b]).unwrap());
        let mut seen = Vec::new();
        hits.retain(|&i| {
            let fresh = !seen.contains(&d.chr[i]);
            if fresh {
                seen.push(d.chr[i]);
            }
            fresh
        });
        0.5
    } else {
        0.45
    };
    // textxy(X, Y, labs, offset = 0.625): labels sit up and to the right of
    // points in the positive quadrant (all of a Manhattan plot).
    let offset = 0.625;
    for i in hits {
        let (x, y) = (d.pos[i], d.logp[i]);
        let xc = if x >= 0.0 { 0.5 - offset } else { 0.5 + offset };
        let yc = if y >= 0.0 { 0.5 - offset } else { 0.5 + offset };
        let label = TextOptions::default().adj(xc).vadj(yc).cex(cex).xpd(true);
        dev.text()
            .x(&[x])
            .y(&[y])
            .labels(&[snp[i].as_str()])
            .with(&label)
            .draw();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts() -> ManhattanOptions {
        ManhattanOptions {
            col: vec![],
            chrlabs: None,
            suggestiveline: None,
            genomewideline: None,
            highlight: None,
            logp: true,
            annotate_pval: None,
            annotate_top: true,
            common: Common::default(),
        }
    }

    #[test]
    fn positions_match_qqman() {
        // Expected values computed with qqman's own loop in R.
        let input = ManhattanInput {
            chr: vec![2.0, 1.0, 1.0, 3.0, 2.0, 3.0],
            bp: vec![150.0, 10.0, 100.0, 1000.0, 50.0, 2000.0],
            p: vec![0.1, 0.01, 0.5, 1e-8, 0.2, 0.3],
            snp: None,
        };
        let d = prepare(input, &opts()).unwrap();
        // chr1: 10, 100; chr2: 50,150 -> 1,101 +100; chr3: 1000,2000 -> 1,1001 + 100 + 101
        assert_eq!(d.pos, vec![10.0, 100.0, 101.0, 201.0, 202.0, 1202.0]);
        assert_eq!(d.ticks, vec![55.0, 151.0, 702.0]);
        assert_eq!(d.labs, vec!["1", "2", "3"]);
        assert_eq!(
            d.xlim,
            ((1202.0f64 * -0.03).floor(), (1202.0f64 * 1.03).ceil())
        );
        assert_eq!(d.ylim, (0.0, 8.0));
    }

    #[test]
    fn single_chromosome_uses_bp() {
        let input = ManhattanInput {
            chr: vec![7.0, 7.0],
            bp: vec![20.0, 10.0],
            p: vec![0.5, 0.01],
            snp: None,
        };
        let d = prepare(input, &opts()).unwrap();
        assert_eq!(d.pos, vec![10.0, 20.0]);
        assert_eq!(d.xlabel, "Chromosome 7 position");
    }
}
