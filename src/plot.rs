//! Options shared by both plots, and writing a figure to SVG/PDF/PNG.

use std::path::Path;

use anyhow::{Context, Result, anyhow, bail};
use rs_graphics::prelude::*;

/// The subset of R's `...` graphical arguments the CLI exposes.
#[derive(Clone, Debug)]
pub struct Common {
    pub main: Option<String>,
    pub xlim: Option<(f64, f64)>,
    pub ylim: Option<(f64, f64)>,
    pub cex: f64,
    pub cex_axis: f64,
    pub cex_lab: f64,
    pub pointsize: f64,
    pub width: f64,
    pub height: f64,
    pub dpi: f64,
    /// Grid cells per inch for `draw::points_dedup`; see [`dedup_grid`].
    pub dedup_grid: f64,
}

/// How finely `draw::points_dedup` tells points apart. Vector output can be
/// zoomed, so it uses 1/1000 inch. PNG uses half a pixel: about 4x faster
/// for millions of points, and on a 5M-SNP plot only a handful of pixels
/// change noticeably (at the soft edges of the densest clusters).
pub fn dedup_grid(out: &Path, dpi: f64) -> f64 {
    let is_png = out
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("png"));
    if is_png { 2.0 * dpi } else { 1000.0 }
}

impl Default for Common {
    fn default() -> Self {
        Common {
            main: None,
            xlim: None,
            ylim: None,
            cex: 1.0,
            cex_axis: 1.0,
            cex_lab: 1.0,
            pointsize: 12.0,
            width: 7.0,
            height: 7.0,
            dpi: 150.0,
            dedup_grid: 1000.0,
        }
    }
}

pub fn resolve_color(spec: ColorSpec, palette: &[Color]) -> Color {
    spec.resolve(palette).unwrap_or(Color::TRANSPARENT)
}

/// Parses an R colour: a name, `#RRGGBB[AA]`, or a palette index.
pub fn parse_color(s: &str) -> Result<ColorSpec> {
    if let Ok(i) = s.parse::<i32>() {
        return Ok(ColorSpec::Index(i));
    }
    match Color::parse(s) {
        Some(c) => Ok(ColorSpec::Color(c)),
        None if s.eq_ignore_ascii_case("NA") || s.eq_ignore_ascii_case("transparent") => {
            Ok(ColorSpec::None)
        }
        None => bail!("invalid color name '{s}'"),
    }
}

/// A plot that can be drawn onto any backend.
pub trait Figure {
    fn draw<B: Backend>(&self, dev: &mut Device<B>);
}

/// Writes `fig` to `out`, choosing the format from the file extension.
pub fn save<F: Figure>(fig: &F, out: &Path, c: &Common) -> Result<()> {
    let backend = AnyBackend::for_path(out, c.width, c.height, c.dpi).map_err(|_| {
        anyhow!(
            "can't tell output format from {}: use .png, .pdf or .svg",
            out.display()
        )
    })?;
    let par = Par {
        ps: c.pointsize,
        cex_axis: c.cex_axis,
        cex_lab: c.cex_lab,
        ..Par::default()
    };
    let mut dev = Device::new(backend, par);
    fig.draw(&mut dev);
    dev.save(out)
        .with_context(|| format!("cannot write {}", out.display()))
}
