//! Options shared by both plots, and writing a figure to SVG/PDF/PNG.

use std::path::Path;

use anyhow::{Context, Result, anyhow};
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

/// Parses an R colour: a name, `#RRGGBB[AA]`, or a palette index (>= 0,
/// as in R). Used as a clap `value_parser`, so a bad `--col` fails before
/// the input is read.
pub fn parse_color(s: &str) -> std::result::Result<ColorSpec, String> {
    if let Ok(i) = s.parse::<i32>() {
        if i < 0 {
            return Err(format!("numerical color values must be >= 0, found {i}"));
        }
        return Ok(ColorSpec::Index(i));
    }
    match Color::parse(s) {
        Some(c) => Ok(ColorSpec::Color(c)),
        None if s.eq_ignore_ascii_case("NA") || s.eq_ignore_ascii_case("transparent") => {
            Ok(ColorSpec::None)
        }
        None => Err(format!("invalid color name '{s}'")),
    }
}

/// A plot that can be drawn onto any backend. Fails with rs-graphics'
/// error where R would stop (a page too small for the margins).
pub trait Figure {
    fn draw<B: Backend>(&self, dev: &mut Device<B>) -> std::result::Result<(), rs_graphics::Error>;
}

/// Writes `fig` to `out`, choosing the format from the file extension.
/// When the page is too small for the plot's margins, fails with R's
/// message (e.g. "figure margins too large") and removes `out`.
pub fn save<F: Figure>(fig: &F, out: &Path, c: &Common) -> Result<()> {
    let format = out
        .extension()
        .and_then(|e| e.to_str())
        .and_then(Format::parse)
        .ok_or_else(|| {
            anyhow!(
                "can't tell output format from {}: use .png, .pdf or .svg",
                out.display()
            )
        })?;
    let mut dev = Device::format(format, out)
        .width(c.width)
        .height(c.height)
        .res(c.dpi)
        .pointsize(c.pointsize)
        .open()?;
    dev.par().cex_axis(c.cex_axis).cex_lab(c.cex_lab).draw();
    if let Err(e) = fig.draw(&mut dev) {
        // `open()` created the file; nothing usable will be written to it.
        drop(dev);
        let _ = std::fs::remove_file(out);
        if e.kind() == ErrorKind::DoesNotFit {
            return Err(anyhow!(
                "{}: use a larger --width/--height or a smaller --pointsize",
                e.message()
            ));
        }
        return Err(e.into());
    }
    dev.close()
        .with_context(|| format!("cannot write {}", out.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Points;

    impl Figure for Points {
        fn draw<B: Backend>(
            &self,
            dev: &mut Device<B>,
        ) -> std::result::Result<(), rs_graphics::Error> {
            dev.plot().x(&[1, 2]).y(&[1, 2]).try_draw()
        }
    }

    #[test]
    fn too_small_a_page_fails_with_r_message() {
        let dir = std::env::temp_dir().join(format!("qqman-small-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let out = dir.join("small.pdf");
        let small = Common {
            width: 1.0,
            height: 1.0,
            ..Common::default()
        };
        let err = save(&Points, &out, &small).unwrap_err().to_string();
        assert!(err.starts_with("figure margins too large: "), "{err}");
        assert!(!out.exists());
        // The default size works.
        save(&Points, &out, &Common::default()).unwrap();
        assert!(out.exists());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn parse_color_follows_r() {
        assert_eq!(parse_color("2"), Ok(ColorSpec::Index(2)));
        assert_eq!(parse_color("0"), Ok(ColorSpec::Index(0)));
        // R ignores spaces in colour names.
        assert_eq!(
            parse_color("light blue"),
            Ok(ColorSpec::Color(Color::rgb(173, 216, 230)))
        );
        assert_eq!(parse_color("NA"), Ok(ColorSpec::None));
        assert_eq!(
            parse_color("-1"),
            Err("numerical color values must be >= 0, found -1".to_string())
        );
        assert_eq!(
            parse_color("notacolour"),
            Err("invalid color name 'notacolour'".to_string())
        );
    }
}
