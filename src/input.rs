//! Reading GWAS results tables.
//!
//! Accepts tab-, comma- or whitespace-delimited text with a header row
//! (the delimiter is detected from the header), optionally gzipped
//! (detected from the magic bytes, not the file name). `-` reads stdin.
//! Only the requested columns are parsed, so wide PLINK/REGENIE output
//! costs little more than a four-column file.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ColKind {
    Num,
    Str,
}

pub struct ColSpec<'a> {
    pub name: &'a str,
    pub kind: ColKind,
    pub required: bool,
}

#[derive(Debug)]
pub enum Column {
    Num(Vec<f64>),
    Str(Vec<String>),
}

impl Column {
    pub fn into_num(self) -> Vec<f64> {
        match self {
            Column::Num(v) => v,
            Column::Str(_) => unreachable!("column was requested as numeric"),
        }
    }

    pub fn into_str(self) -> Vec<String> {
        match self {
            Column::Str(v) => v,
            Column::Num(_) => unreachable!("column was requested as string"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Delim {
    Tab,
    Comma,
    Whitespace,
}

impl Delim {
    fn detect(header: &str) -> Delim {
        if header.contains('\t') {
            Delim::Tab
        } else if header.contains(',') {
            Delim::Comma
        } else {
            Delim::Whitespace
        }
    }

    /// Calls `f(index, field)` for the first `n` fields of `line`, returning
    /// how many there were (at most `n`).
    fn for_each_field<E>(
        self,
        line: &str,
        n: usize,
        mut f: impl FnMut(usize, &str) -> Result<(), E>,
    ) -> Result<usize, E> {
        let mut count = 0;
        let mut visit = |field: &str| {
            let r = f(count, field);
            count += 1;
            r
        };
        match self {
            Delim::Tab => line.split('\t').take(n).try_for_each(|x| visit(x.trim()))?,
            Delim::Comma => line
                .split(',')
                .take(n)
                .try_for_each(|x| visit(x.trim().trim_matches('"')))?,
            Delim::Whitespace => line.split_whitespace().take(n).try_for_each(visit)?,
        }
        Ok(count)
    }
}

fn open(path: &Path) -> Result<Box<dyn BufRead>> {
    let raw: Box<dyn Read> = if path == Path::new("-") {
        Box::new(io::stdin().lock())
    } else {
        Box::new(File::open(path).with_context(|| format!("cannot open {}", path.display()))?)
    };
    let mut reader = BufReader::with_capacity(1 << 20, raw);
    let is_gzip = reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    Ok(if is_gzip {
        Box::new(BufReader::with_capacity(
            1 << 20,
            MultiGzDecoder::new(reader),
        ))
    } else {
        Box::new(reader)
    })
}

/// R's `read.table()` treats these as missing; they become `NaN`.
fn parse_num(field: &str) -> Option<f64> {
    match field {
        "" | "NA" | "NaN" | "nan" | "." => Some(f64::NAN),
        _ => field.parse().ok(),
    }
}

/// Reads the requested columns from `path`. The result has one entry per
/// spec, `None` for an optional column that isn't in the header.
pub fn read_columns(path: &Path, specs: &[ColSpec]) -> Result<Vec<Option<Column>>> {
    let mut reader = open(path)?;
    let mut line = String::new();
    if reader.read_line(&mut line)? == 0 {
        bail!("{} is empty", path.display());
    }
    let header = line.trim_end_matches(['\n', '\r']).to_string();
    let delim = Delim::detect(&header);
    let mut fields: Vec<String> = Vec::new();
    delim
        .for_each_field::<()>(&header, usize::MAX, |_, f| {
            fields.push(f.to_string());
            Ok(())
        })
        .unwrap();

    let mut wanted = Vec::new(); // (field index, spec index)
    let mut columns: Vec<Option<Column>> = Vec::with_capacity(specs.len());
    for (si, spec) in specs.iter().enumerate() {
        match fields.iter().position(|f| *f == spec.name) {
            Some(fi) => {
                wanted.push((fi, si));
                columns.push(Some(match spec.kind {
                    ColKind::Num => Column::Num(Vec::new()),
                    ColKind::Str => Column::Str(Vec::new()),
                }));
            }
            None if spec.required => bail!("Column {} not found!", spec.name),
            None => columns.push(None),
        }
    }
    let n_needed = wanted.iter().map(|&(fi, _)| fi + 1).max().unwrap_or(0);
    let mut slot = vec![None; n_needed]; // field index -> spec index
    for &(fi, si) in &wanted {
        slot[fi] = Some(si);
    }

    let mut lineno = 1;
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        lineno += 1;
        let row = line.trim_end_matches(['\n', '\r']);
        if row.trim().is_empty() {
            continue;
        }
        let found = delim.for_each_field(row, n_needed, |fi, field| {
            let Some(si) = slot[fi] else { return Ok(()) };
            match columns[si].as_mut().unwrap() {
                Column::Num(v) => match parse_num(field) {
                    Some(x) => v.push(x),
                    None => bail!(
                        "line {lineno}: {} column should be numeric, found {field:?}",
                        specs[si].name
                    ),
                },
                Column::Str(v) => v.push(field.to_string()),
            }
            Ok(())
        })?;
        if found < n_needed {
            bail!("line {lineno}: expected at least {n_needed} fields, found {found}");
        }
    }
    Ok(columns)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_delimiters() {
        assert_eq!(Delim::detect("SNP\tCHR\tBP"), Delim::Tab);
        assert_eq!(Delim::detect("SNP,CHR,BP"), Delim::Comma);
        assert_eq!(Delim::detect(" CHR   SNP  BP"), Delim::Whitespace);
    }

    #[test]
    fn missing_values_are_nan() {
        assert!(parse_num("NA").unwrap().is_nan());
        assert_eq!(parse_num("1e-8"), Some(1e-8));
        assert_eq!(parse_num("X"), None);
    }
}
