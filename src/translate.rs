//! DNA -> protein translation (standard genetic code), mirroring gclib's
//! `translateDNA()`. Stop codons are emitted as `.` (or `*` with `star_stop`).

/// Standard codon table: index = (nt1, nt2, nt3) packed as base-4 with
/// A=0,C=1,G=2,T/U=3. Unknown codons translate to `X`.
fn codon_index(c: &[u8]) -> usize {
    let mut idx = 0usize;
    for &b in c {
        let v = match b.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' | b'U' => 3,
            _ => return usize::MAX,
        };
        idx = idx * 4 + v;
    }
    idx
}

/// Lookup the amino acid for a 3-letter codon (standard table).
/// Handles IUPAC ambiguous bases: if a codon contains ambiguous bases
/// (N, R, Y, S, W, K, M, B, D, H, V), expands all possible unambiguous
/// codons and returns the amino acid if all expansions agree; otherwise 'X'.
/// This mirrors gclib's codon table which has explicit entries for
/// ambiguous codons like GGN→G (all Gly), GTN→V (all Val), etc.
fn aa_for_codon(cod: &[u8]) -> u8 {
    // Direct codon string match is simplest and avoids table-sizing issues.
    let c = [cod[0].to_ascii_uppercase(), cod[1].to_ascii_uppercase(), cod[2].to_ascii_uppercase()];
    match c {
        [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => b'F',
        [b'T', b'T', b'A'] | [b'T', b'T', b'G'] | [b'C', b'T', b'A'] | [b'C', b'T', b'G']
        | [b'C', b'T', b'T'] | [b'C', b'T', b'C'] => b'L',
        [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => b'I',
        [b'A', b'T', b'G'] => b'M',
        [b'G', b'T', b'T'] | [b'G', b'T', b'C'] | [b'G', b'T', b'A'] | [b'G', b'T', b'G'] => b'V',
        [b'T', b'C', b'T'] | [b'T', b'C', b'C'] | [b'T', b'C', b'A'] | [b'T', b'C', b'G'] => b'S',
        [b'C', b'C', b'T'] | [b'C', b'C', b'C'] | [b'C', b'C', b'A'] | [b'C', b'C', b'G'] => b'P',
        [b'A', b'C', b'T'] | [b'A', b'C', b'C'] | [b'A', b'C', b'A'] | [b'A', b'C', b'G'] => b'T',
        [b'G', b'C', b'T'] | [b'G', b'C', b'C'] | [b'G', b'C', b'A'] | [b'G', b'C', b'G'] => b'A',
        [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => b'Y',
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] => b'.', // stop
        [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => b'H',
        [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => b'Q',
        [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => b'N',
        [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => b'K',
        [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => b'D',
        [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => b'E',
        [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => b'C',
        [b'T', b'G', b'A'] => b'.', // stop
        [b'T', b'G', b'G'] => b'W',
        [b'C', b'G', b'T'] | [b'C', b'G', b'C'] | [b'C', b'G', b'A'] | [b'C', b'G', b'G']
        | [b'A', b'G', b'A'] | [b'A', b'G', b'G'] => b'R',
        [b'A', b'G', b'T'] | [b'A', b'G', b'C'] => b'S',
        [b'G', b'G', b'T'] | [b'G', b'G', b'C'] | [b'G', b'G', b'A'] | [b'G', b'G', b'G'] => b'G',
        _ => {
            // Codon contains an IUPAC ambiguous base. Expand all possible
            // unambiguous codons; if they all translate to the same amino
            // acid, return it, otherwise return 'X'. This mirrors gclib's
            // codon table behavior (e.g. GGN→G because all GG? are Gly).
            let expansions = expand_ambiguous_codon(&c);
            if expansions.is_empty() {
                return b'X';
            }
            let first_aa = aa_for_codon_simple(&expansions[0]);
            for e in &expansions[1..] {
                if aa_for_codon_simple(e) != first_aa {
                    return b'X';
                }
            }
            first_aa
        }
    }
}

/// Expand an IUPAC ambiguous codon into all possible unambiguous codons.
fn expand_ambiguous_codon(c: &[u8; 3]) -> Vec<[u8; 3]> {
    fn expand_base(b: u8) -> Vec<u8> {
        match b {
            b'A' => vec![b'A'],
            b'C' => vec![b'C'],
            b'G' => vec![b'G'],
            b'T' | b'U' => vec![b'T'],
            b'N' => vec![b'A', b'C', b'G', b'T'],
            b'R' => vec![b'A', b'G'],
            b'Y' => vec![b'C', b'T'],
            b'S' => vec![b'G', b'C'],
            b'W' => vec![b'A', b'T'],
            b'K' => vec![b'G', b'T'],
            b'M' => vec![b'A', b'C'],
            b'B' => vec![b'C', b'G', b'T'],
            b'D' => vec![b'A', b'G', b'T'],
            b'H' => vec![b'A', b'C', b'T'],
            b'V' => vec![b'A', b'C', b'G'],
            _ => vec![],
        }
    }
    let b1 = expand_base(c[0]);
    let b2 = expand_base(c[1]);
    let b3 = expand_base(c[2]);
    let mut out = Vec::new();
    for &x in &b1 {
        for &y in &b2 {
            for &z in &b3 {
                out.push([x, y, z]);
            }
        }
    }
    out
}

/// Strict unambiguous codon lookup (no IUPAC expansion). Callers must ensure
/// the codon contains only A/C/G/T.
fn aa_for_codon_simple(c: &[u8; 3]) -> u8 {
    match c {
        [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => b'F',
        [b'T', b'T', b'A'] | [b'T', b'T', b'G'] | [b'C', b'T', b'A'] | [b'C', b'T', b'G']
        | [b'C', b'T', b'T'] | [b'C', b'T', b'C'] => b'L',
        [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => b'I',
        [b'A', b'T', b'G'] => b'M',
        [b'G', b'T', b'T'] | [b'G', b'T', b'C'] | [b'G', b'T', b'A'] | [b'G', b'T', b'G'] => b'V',
        [b'T', b'C', b'T'] | [b'T', b'C', b'C'] | [b'T', b'C', b'A'] | [b'T', b'C', b'G'] => b'S',
        [b'C', b'C', b'T'] | [b'C', b'C', b'C'] | [b'C', b'C', b'A'] | [b'C', b'C', b'G'] => b'P',
        [b'A', b'C', b'T'] | [b'A', b'C', b'C'] | [b'A', b'C', b'A'] | [b'A', b'C', b'G'] => b'T',
        [b'G', b'C', b'T'] | [b'G', b'C', b'C'] | [b'G', b'C', b'A'] | [b'G', b'C', b'G'] => b'A',
        [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => b'Y',
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] => b'.',
        [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => b'H',
        [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => b'Q',
        [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => b'N',
        [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => b'K',
        [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => b'D',
        [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => b'E',
        [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => b'C',
        [b'T', b'G', b'A'] => b'.',
        [b'T', b'G', b'G'] => b'W',
        [b'C', b'G', b'T'] | [b'C', b'G', b'C'] | [b'C', b'G', b'A'] | [b'C', b'G', b'G']
        | [b'A', b'G', b'A'] | [b'A', b'G', b'G'] => b'R',
        [b'A', b'G', b'T'] | [b'A', b'G', b'C'] => b'S',
        [b'G', b'G', b'T'] | [b'G', b'G', b'C'] | [b'G', b'G', b'A'] | [b'G', b'G', b'G'] => b'G',
        _ => b'X',
    }
}

/// Translate a DNA byte string. Returns the amino-acid sequence.
/// `len` is the number of nucleotides to translate (must be a multiple of 3 for
/// full codons; trailing bases are ignored). Stop codons are `.`.
pub fn translate_dna(nt: &[u8], len: usize) -> Vec<u8> {
    let n = len.min(nt.len());
    let n = n - (n % 3);
    let mut out = Vec::with_capacity(n / 3);
    let mut i = 0;
    while i + 3 <= n {
        let cod = &nt[i..i + 3];
        out.push(aa_for_codon(cod));
        i += 3;
    }
    let _ = codon_index; // referenced to silence dead_code if table helper unused
    out
}

/// Write a FASTA sequence (70 chars/line, matching gffread) to `out`. `use_star`
/// replaces `.` (stop) with `*`.
pub fn print_fasta(
    out: &mut dyn std::io::Write,
    defline: Option<&str>,
    seq: &[u8],
    use_star: bool,
) -> std::io::Result<()> {
    if let Some(d) = defline {
        writeln!(out, ">{}", d)?;
    }
    let mut ilen = 0usize;
    for &b in seq {
        if ilen == 70 {
            writeln!(out)?;
            ilen = 0;
        }
        let c = if use_star && b == b'.' { b'*' } else { b };
        out.write_all(&[c])?;
        ilen += 1;
    }
    writeln!(out)?;
    Ok(())
}
