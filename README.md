# qqman-rs

`qqman` is a command-line tool for drawing Manhattan and Q-Q plots of GWAS
results. It is a Rust port of the R package
[qqman](https://github.com/stephenturner/qqman): the plots look the same as
qqman's, but there is no R to install, and it is much faster. A Manhattan plot
of 5 million SNPs takes about a second as a PNG, where R/qqman takes around
35 seconds.

> Status: early development. Options may change.

## Installation

Install the [Rust toolchain](https://www.rust-lang.org/tools/install), then:

```
cargo install --git ssh://git@github.com/MRCIEU/qqman-rs.git
```

This puts the `qqman` binary in Cargo's bin directory (usually
`~/.cargo/bin`). The repositories are fetched over SSH, so you need a GitHub
SSH key with access to them. If Cargo can't authenticate, set
`CARGO_NET_GIT_FETCH_WITH_CLI=true` so it uses your system `git`.

To build from a clone instead:

```
git clone git@github.com:MRCIEU/qqman-rs.git
cd qqman-rs
cargo build --release   # binary at target/release/qqman
```

## Usage

```
qqman manhattan results.tsv -o manhattan.png
qqman qq results.tsv -o qq.pdf
```

- **Input** is a table with a header row: tab-, comma- or
  whitespace-separated, and optionally gzipped. Both are detected
  automatically, so `results.tsv.gz` works as is. Use `-` to read from stdin:

  ```
  # keep the header row, drop chromosomes above 22 (CHR is column 2 here)
  awk 'NR == 1 || $2 <= 22' results.tsv | qqman manhattan - -o manhattan.png
  ```

- **Columns** default to qqman's names: `CHR`, `BP`, `P` and `SNP`. Use
  `--chr`, `--bp`, `--p` and `--snp` for other names, e.g. PLINK 2's
  `--chr '#CHROM' --bp POS --snp ID`. Chromosomes must be numeric, as in qqman.
- **Output format** comes from the `-o` file extension: `.png`, `.pdf` or
  `.svg`. Manhattan plots default to 14 × 7 inches and Q-Q plots to 7 × 7;
  change this with `--width`/`--height`, and PNG resolution with `--dpi`
  (default 150).
- **Options follow qqman's arguments**, in kebab-case: `col`, `chrlabs`,
  `suggestiveline`, `genomewideline`, `highlight`, `annotatePval`, plus
  graphical parameters such as `main`, `ylim`, `cex` and `cex.axis`:

  ```
  qqman manhattan results.tsv -o manhattan.pdf \
      --main "Manhattan Plot" --col blue4,orange3 --cex 0.6 \
      --highlight-file snps.txt --annotate-pval 1e-6
  ```

  `--suggestiveline false` or `--genomewideline false` turns a line off.

- **Missing values:** Manhattan plots skip rows with a missing, non-finite
  or zero value, with a warning (R/qqman would fail on them). Q-Q plots use
  only P values strictly between 0 and 1, as qqman's `qq()` does.

## Help

<details>
<summary><code>qqman --help</code></summary>

```
Usage: qqman [OPTIONS] [COMMAND]

Commands:
  manhattan  Manhattan plot of GWAS results (qqman::manhattan)
  qq         Q-Q plot of p-values (qqman::qq)
  help       Print this message or the help of the given subcommand(s)

Options:
  -l              Print a summary of licence information
      --licences  Print full licence information, including all third-party notices
  -h, --help      Print help
  -V, --version   Print version
```
</details>

<details>
<summary><code>qqman manhattan --help</code></summary>

```
Usage: qqman manhattan [OPTIONS] --output <OUTPUT> <INPUT>

Arguments:
  <INPUT>  Results table (tab/comma/whitespace-delimited, optionally gzipped; `-` for stdin)

Options:
      --chr <CHR>
          Chromosome column (must be numeric) [default: CHR]
      --bp <BP>
          Base-pair position column [default: BP]
      --p <P>
          P-value column [default: P]
      --snp <SNP>
          SNP name column (needed for --highlight and --annotate-pval) [default: SNP]
      --col <COL>
          Colours to cycle through by chromosome [default: gray10,gray60]
      --chrlabs <CHRLABS>
          Chromosome labels, one per chromosome, comma-separated
      --suggestiveline <SUGGESTIVELINE>
          Where to draw the blue suggestive line, or `false` [default: 5]
      --genomewideline <GENOMEWIDELINE>
          Where to draw the red genome-wide line, or `false` [default: 7.301029995663981]
      --highlight <HIGHLIGHT>
          SNPs to highlight in green, comma-separated
      --highlight-file <HIGHLIGHT_FILE>
          File of SNPs to highlight, one per line
      --no-logp
          Plot the P column as-is instead of -log10(P)
      --annotate-pval <ANNOTATE_PVAL>
          Label SNPs with P at or below this threshold
      --annotate-all
          With --annotate-pval, label every SNP past the threshold rather than only the top SNP on each chromosome
  -o, --output <OUTPUT>
          Output file; the format comes from the extension (.png, .pdf, .svg)
      --main <MAIN>
          Plot title
      --xlim <XLIM>
          X-axis limits, as `min,max`
      --ylim <YLIM>
          Y-axis limits, as `min,max`
      --cex <CEX>
          Point size multiplier (R's `cex`) [default: 1]
      --cex-axis <CEX_AXIS>
          Axis label size multiplier (R's `cex.axis`) [default: 1]
      --cex-lab <CEX_LAB>
          Axis title size multiplier (R's `cex.lab`) [default: 1]
      --pointsize <POINTSIZE>
          Base font size in points [default: 12]
      --width <WIDTH>
          Width in inches [default: 14 for manhattan, 7 for qq]
      --height <HEIGHT>
          Height in inches [default: 7]
      --dpi <DPI>
          Resolution of PNG output, in pixels per inch [default: 150]
  -h, --help
          Print help
```
</details>

<details>
<summary><code>qqman qq --help</code></summary>

```
Usage: qqman qq [OPTIONS] --output <OUTPUT> <INPUT>

Arguments:
  <INPUT>  Results table (tab/comma/whitespace-delimited, optionally gzipped; `-` for stdin)

Options:
      --p <P>                  P-value column [default: P]
  -o, --output <OUTPUT>        Output file; the format comes from the extension (.png, .pdf, .svg)
      --main <MAIN>            Plot title
      --xlim <XLIM>            X-axis limits, as `min,max`
      --ylim <YLIM>            Y-axis limits, as `min,max`
      --cex <CEX>              Point size multiplier (R's `cex`) [default: 1]
      --cex-axis <CEX_AXIS>    Axis label size multiplier (R's `cex.axis`) [default: 1]
      --cex-lab <CEX_LAB>      Axis title size multiplier (R's `cex.lab`) [default: 1]
      --pointsize <POINTSIZE>  Base font size in points [default: 12]
      --width <WIDTH>          Width in inches [default: 14 for manhattan, 7 for qq]
      --height <HEIGHT>        Height in inches [default: 7]
      --dpi <DPI>              Resolution of PNG output, in pixels per inch [default: 150]
  -h, --help                   Print help
```
</details>

## Licence

GPL-3.0-only, like R/qqman; see [LICENSE](LICENSE) and [NOTICE](NOTICE).
`qqman -l` prints a summary.
