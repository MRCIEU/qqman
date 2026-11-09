//! `qqman`: Manhattan and Q-Q plots for GWAS results, after the R package
//! qqman (https://github.com/stephenturner/qqman). Licensed under GPL-3.0-only;
//! see LICENSE.

mod draw;
mod input;
mod manhattan;
mod plot;
mod qq;

use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Args, CommandFactory, Parser, Subcommand};
use rs_graphics::backend::Backend;
use rs_graphics::color::ColorSpec;
use rs_graphics::device::Device;

use input::{ColKind, ColSpec, read_columns};
use plot::{Common, Figure, parse_color, save};

#[derive(Parser)]
#[command(version, about, arg_required_else_help = true)]
struct Cli {
    /// Print a summary of licence information
    #[arg(short = 'l', exclusive = true)]
    licences_summary: bool,
    /// Print full licence information, including all third-party notices
    #[arg(long = "licences", alias = "licenses", exclusive = true)]
    licences: bool,
    #[command(subcommand)]
    command: Option<Command>,
}

/// The NOTICE file: this program's licence, and the notices of its third-party dependencies
const NOTICE: &str = include_str!("../NOTICE");

/// The opening statement of NOTICE, i.e. everything before the first third-party notice
fn get_notice_summary() -> &'static str {
    NOTICE
        .split_once("\n---")
        .map_or(NOTICE, |(summary, _)| summary)
}

// The variants differ in size, but a `Command` is built once per run, so
// boxing the larger one would buy nothing.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
enum Command {
    /// Manhattan plot of GWAS results (qqman::manhattan).
    Manhattan(ManhattanArgs),
    /// Q-Q plot of p-values (qqman::qq).
    Qq(QqArgs),
}

/// Graphical parameters shared by both plots.
#[derive(Args)]
struct CommonArgs {
    /// Output file; the format comes from the extension (.png, .pdf, .svg).
    #[arg(short, long)]
    output: PathBuf,
    /// Plot title.
    #[arg(long)]
    main: Option<String>,
    /// X-axis limits, as `min,max`.
    #[arg(long, value_parser = parse_range, allow_hyphen_values = true)]
    xlim: Option<(f64, f64)>,
    /// Y-axis limits, as `min,max`.
    #[arg(long, value_parser = parse_range, allow_hyphen_values = true)]
    ylim: Option<(f64, f64)>,
    /// Point size multiplier (R's `cex`).
    #[arg(long, default_value_t = 1.0)]
    cex: f64,
    /// Axis label size multiplier (R's `cex.axis`).
    #[arg(long, default_value_t = 1.0)]
    cex_axis: f64,
    /// Axis title size multiplier (R's `cex.lab`).
    #[arg(long, default_value_t = 1.0)]
    cex_lab: f64,
    /// Base font size in points.
    #[arg(long, default_value_t = 12.0)]
    pointsize: f64,
    /// Width in inches [default: 14 for manhattan, 7 for qq].
    #[arg(long)]
    width: Option<f64>,
    /// Height in inches.
    #[arg(long, default_value_t = 7.0)]
    height: f64,
    /// Resolution of PNG output, in pixels per inch.
    #[arg(long, default_value_t = 150.0)]
    dpi: f64,
}

impl CommonArgs {
    /// `default_width` is the subcommand's default for `--width`.
    fn common(&self, default_width: f64) -> Common {
        Common {
            main: self.main.clone(),
            xlim: self.xlim,
            ylim: self.ylim,
            cex: self.cex,
            cex_axis: self.cex_axis,
            cex_lab: self.cex_lab,
            pointsize: self.pointsize,
            width: self.width.unwrap_or(default_width),
            height: self.height,
            dpi: self.dpi,
            dedup_grid: plot::dedup_grid(&self.output, self.dpi),
        }
    }
}

#[derive(Args)]
struct ManhattanArgs {
    /// Results table (tab/comma/whitespace-delimited, optionally gzipped; `-` for stdin).
    input: PathBuf,
    /// Chromosome column (must be numeric).
    #[arg(long, default_value = "CHR")]
    chr: String,
    /// Base-pair position column.
    #[arg(long, default_value = "BP")]
    bp: String,
    /// P-value column.
    #[arg(long, default_value = "P")]
    p: String,
    /// SNP name column (needed for --highlight and --annotate-pval).
    #[arg(long, default_value = "SNP")]
    snp: String,
    /// Colours to cycle through by chromosome: names, `#RRGGBB[AA]` or
    /// palette numbers (>= 0).
    #[arg(
        long,
        value_delimiter = ',',
        value_parser = parse_color,
        allow_hyphen_values = true,
        default_value = "gray10,gray60"
    )]
    col: Vec<ColorSpec>,
    /// Chromosome labels, one per chromosome, comma-separated.
    #[arg(long, value_delimiter = ',')]
    chrlabs: Option<Vec<String>>,
    /// Where to draw the blue suggestive line, or `false`.
    #[arg(long, value_parser = parse_line, default_value = "5")]
    suggestiveline: Line,
    /// Where to draw the red genome-wide line, or `false`.
    #[arg(long, value_parser = parse_line, default_value = "7.301029995663981")]
    genomewideline: Line,
    /// SNPs to highlight in green, comma-separated.
    #[arg(long, value_delimiter = ',')]
    highlight: Option<Vec<String>>,
    /// File of SNPs to highlight, one per line.
    #[arg(long)]
    highlight_file: Option<PathBuf>,
    /// Plot the P column as-is instead of -log10(P).
    #[arg(long)]
    no_logp: bool,
    /// Label SNPs with P at or below this threshold.
    #[arg(long)]
    annotate_pval: Option<f64>,
    /// With --annotate-pval, label every SNP past the threshold rather than
    /// only the top SNP on each chromosome.
    #[arg(long)]
    annotate_all: bool,
    #[command(flatten)]
    common: CommonArgs,
}

#[derive(Args)]
struct QqArgs {
    /// Results table (tab/comma/whitespace-delimited, optionally gzipped; `-` for stdin).
    input: PathBuf,
    /// P-value column.
    #[arg(long, default_value = "P")]
    p: String,
    #[command(flatten)]
    common: CommonArgs,
}

/// Default widths in inches. Manhattan plots are usually wide; a Q-Q plot is
/// read against its y = x line, so it stays square (R's 7 x 7 `pdf()` size).
const MANHATTAN_WIDTH: f64 = 14.0;
const QQ_WIDTH: f64 = 7.0;

fn parse_range(s: &str) -> Result<(f64, f64), String> {
    let (a, b) = s.split_once(',').ok_or("expected `min,max`")?;
    let num = |v: &str| v.trim().parse::<f64>().map_err(|e| format!("{v:?}: {e}"));
    Ok((num(a)?, num(b)?))
}

/// A threshold line's position, or `None` for R's `FALSE`.
#[derive(Clone, Copy)]
struct Line(Option<f64>);

fn parse_line(s: &str) -> Result<Line, String> {
    match s.to_ascii_lowercase().as_str() {
        "false" | "none" | "f" => Ok(Line(None)),
        v => v
            .parse()
            .map(|h| Line(Some(h)))
            .map_err(|e| format!("{s:?}: {e}")),
    }
}

struct ManhattanFigure(manhattan::Prepared, manhattan::ManhattanOptions);

impl Figure for ManhattanFigure {
    fn draw<B: Backend>(&self, dev: &mut Device<B>) -> Result<(), rs_graphics::Error> {
        manhattan::draw(dev, &self.0, &self.1)
    }
}

struct QqFigure(qq::Prepared, Common);

impl Figure for QqFigure {
    fn draw<B: Backend>(&self, dev: &mut Device<B>) -> Result<(), rs_graphics::Error> {
        qq::draw(dev, &self.0, &self.1)
    }
}

fn run_manhattan(args: ManhattanArgs) -> Result<()> {
    let mut highlight: Option<HashSet<String>> = args.highlight.map(|v| v.into_iter().collect());
    if let Some(path) = &args.highlight_file {
        let text =
            fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
        highlight
            .get_or_insert_default()
            .extend(text.split_whitespace().map(String::from));
    }
    let need_snp = highlight.is_some() || args.annotate_pval.is_some();

    let t = Instant::now();
    let specs = [
        ColSpec {
            name: &args.chr,
            kind: ColKind::Num,
            required: true,
        },
        ColSpec {
            name: &args.bp,
            kind: ColKind::Num,
            required: true,
        },
        ColSpec {
            name: &args.p,
            kind: ColKind::Num,
            required: true,
        },
        ColSpec {
            name: &args.snp,
            kind: ColKind::Str,
            required: false,
        },
    ];
    let mut cols =
        read_columns(&args.input, if need_snp { &specs } else { &specs[..3] })?.into_iter();
    let chr = cols.next().flatten().unwrap().into_num();
    let bp = cols.next().flatten().unwrap().into_num();
    let p = cols.next().flatten().unwrap().into_num();
    let snp = cols.next().flatten().map(|c| c.into_str());
    if need_snp && snp.is_none() {
        eprintln!("warning: No SNP column found. OK unless you're trying to highlight.");
    }
    eprintln!("read {} rows in {:.2?}", chr.len(), t.elapsed());

    let opts = manhattan::ManhattanOptions {
        col: args.col,
        chrlabs: args.chrlabs,
        suggestiveline: args.suggestiveline.0,
        genomewideline: args.genomewideline.0,
        highlight,
        logp: !args.no_logp,
        annotate_pval: args.annotate_pval,
        annotate_top: !args.annotate_all,
        common: args.common.common(MANHATTAN_WIDTH),
    };
    let prepared = manhattan::prepare(manhattan::ManhattanInput { chr, bp, p, snp }, &opts)?;
    let t = Instant::now();
    let common = opts.common.clone();
    save(
        &ManhattanFigure(prepared, opts),
        &args.common.output,
        &common,
    )?;
    eprintln!(
        "wrote {} in {:.2?}",
        args.common.output.display(),
        t.elapsed()
    );
    Ok(())
}

fn run_qq(args: QqArgs) -> Result<()> {
    let t = Instant::now();
    let spec = [ColSpec {
        name: &args.p,
        kind: ColKind::Num,
        required: true,
    }];
    let p = read_columns(&args.input, &spec)?
        .remove(0)
        .unwrap()
        .into_num();
    eprintln!("read {} rows in {:.2?}", p.len(), t.elapsed());
    let common = args.common.common(QQ_WIDTH);
    let prepared = qq::prepare(p, &common)?;
    let t = Instant::now();
    save(
        &QqFigure(prepared, common.clone()),
        &args.common.output,
        &common,
    )?;
    eprintln!(
        "wrote {} in {:.2?}",
        args.common.output.display(),
        t.elapsed()
    );
    Ok(())
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    if cli.licences {
        print!("{NOTICE}");
        return Ok(());
    }
    if cli.licences_summary {
        println!("{}", get_notice_summary());
        println!(
            "Use \"qqman --licences\" to print the full text, including all third-party notices."
        );
        return Ok(());
    }
    match cli.command {
        Some(Command::Manhattan(args)) => run_manhattan(args),
        Some(Command::Qq(args)) => run_qq(args),
        None => Ok(Cli::command().print_help()?),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_notice_summary_stops_before_first_third_party_notice() {
        let summary = get_notice_summary();
        assert!(summary.starts_with("qqman is licensed under the GNU General Public License"));
        assert!(!summary.contains("---"));
        assert!(summary.len() < NOTICE.len());
    }

    #[test]
    fn negative_col_is_rejected_while_parsing_arguments() {
        let parse = |col: &str| {
            Cli::try_parse_from(["qqman", "manhattan", "in.tsv", "-o", "out.png", "--col", col])
        };
        let err = parse("1,-2").err().expect("a negative colour is an error");
        assert!(
            err.to_string()
                .contains("numerical color values must be >= 0, found -2"),
            "{err}"
        );
        assert!(parse("1,2").is_ok());
    }
}
