//! Output formatting: GFF3 / GTF / BED / TLF / table printing.

use std::io::Write;

use crate::config::{AppConfig, OutFormat};
use crate::types::{GffObj, GffPrintMode, GSeg, TableField};

/// Decode `%XX` hex-escaped characters in an attribute value (-D).
pub fn decode_hex_chars(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' && i + 2 < bytes.len() {
            let hi = hex_val(bytes[i + 1]);
            let lo = hex_val(bytes[i + 2]);
            if let (Some(h), Some(l)) = (hi, lo) {
                out.push((h << 4) | l);
                i += 3;
                continue;
            }
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn hex_val(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        b'A'..=b'F' => Some(b - b'A' + 10),
        _ => None,
    }
}

/// Public accessor for FASTA defline attribute printing (used by filter.rs).
///
/// In the C++ gclib, `extractGFFAttr` removes `ID=`, `Parent=` (GFF3 input)
/// and `transcript_id`, `gene_id` (GTF input) from the attribute string
/// *before* `parseAttrs` runs, so they never appear in the attrs list and are
/// never printed in FASTA deflines. The Rust parser stores all attributes,
/// so we replicate the consumption here by filtering them out.
pub fn __attrs_for_fasta<'a>(cfg: &AppConfig, t: &'a GffObj) -> Vec<(&'a str, &'a str)> {
    let is_gff3 = t.get_attr("ID").is_some();
    let consumed: &[&str] = if is_gff3 {
        &["ID", "Parent"]
    } else {
        &["transcript_id", "gene_id"]
    };

    if !cfg.attr_list.is_empty() {
        let mut out = Vec::new();
        for (k, v) in &t.attrs {
            if consumed.contains(&k.as_str()) {
                continue;
            }
            if cfg.attr_list.contains(k) {
                out.push((k.as_str(), v.as_str()));
            }
        }
        return out;
    }

    if cfg.loader.full_attributes {
        return t.attrs.iter()
            .filter(|(k, _)| !consumed.contains(&k.as_str()))
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
    }

    let mut prog: Vec<&str> = Vec::new();
    if cfg.add_cds_attrs || cfg.adjust_stop {
        prog.extend_from_slice(&["partialness", "InFrameStop", "CDstopAdjusted"]);
    }
    if cfg.add_has_cds {
        prog.push("hasCDS");
    }
    let mut out = Vec::new();
    for (k, v) in &t.attrs {
        if prog.contains(&k.as_str()) {
            out.push((k.as_str(), v.as_str()));
        }
    }
    out
}

/// Choose the attributes to print for a record, honoring -F / --attrs / -G.
fn attrs_to_print<'a>(cfg: &AppConfig, t: &'a GffObj) -> Vec<(&'a str, &'a str)> {
    if !cfg.attr_list.is_empty() {
        // --attrs filter: only listed attributes
        let mut out = Vec::new();
        for (k, v) in &t.attrs {
            if cfg.attr_list.contains(k) {
                out.push((k.as_str(), v.as_str()));
            }
        }
        return out;
    }
    if cfg.loader.full_attributes {
        return t.attrs.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
    }
    // minimal attribute set
    let mut keep: Vec<&str> =
        vec!["Name", "gene_name", "geneID", "gene_id", "note", "product", "Parent"];
    // Attributes added by -P / --add-hasCDS / --adj-stop should always survive
    // the minimal-attribute filter so the user can see them.
    if cfg.add_cds_attrs || cfg.adjust_stop {
        keep.push("partialness");
        keep.push("InFrameStop");
        keep.push("CDstopAdjusted");
    }
    if cfg.add_has_cds {
        keep.push("hasCDS");
    }
    let mut out = Vec::new();
    for (k, v) in &t.attrs {
        if keep.contains(&k.as_str()) {
            out.push((k.as_str(), v.as_str()));
        }
    }
    out
}

fn fmt_attr_value(v: &str, decode: bool) -> String {
    if decode {
        decode_hex_chars(v)
    } else {
        v.to_string()
    }
}

/// Build the GFF3 attribute column (ID=...;Parent=...;...).
fn gff3_attrs(cfg: &AppConfig, t: &GffObj, include_id: bool, include_parent: bool) -> String {
    let decode = cfg.decode_chars;
    let mut parts: Vec<String> = Vec::new();
    if include_id && !t.id.is_empty() {
        parts.push(format!("ID={}", t.id));
    }
    // Native gclib: for transcript records, when the parent gene is NOT
    // printed (the default, unless --keep-genes), the gene association is
    // emitted as `geneID=` instead of `Parent=`. This matches the behavior
    // where genes are discarded during finalize() so `parent->isDiscarded()`
    // is true and `Parent=` is suppressed.
    let use_geneid_for_parent = t.is_mrna && !cfg.loader.keep_genes;
    if include_parent {
        if let Some(p) = &t.parent_id {
            if use_geneid_for_parent {
                parts.push(format!("geneID={}", p));
            } else {
                parts.push(format!("Parent={}", p));
            }
        }
    }
    let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
    if include_id {
        seen.insert("ID");
    }
    if include_parent {
        // The Parent link is always consumed (either emitted as Parent= or
        // geneID=), so never re-print a stored Parent attribute from attrs.
        seen.insert("Parent");
        if use_geneid_for_parent {
            seen.insert("geneID");
        }
    }
    for (k, v) in attrs_to_print(cfg, t) {
        if seen.contains(k) {
            continue;
        }
        // Skip a gene_id attr that duplicates the geneID we just printed.
        if use_geneid_for_parent && k == "gene_id" {
            if let Some(p) = &t.parent_id {
                if v == p {
                    continue;
                }
            }
        }
        seen.insert(k);
        parts.push(format!("{}={}", k, fmt_attr_value(v, decode)));
    }
    parts.join(";")
}

/// Native gclib resolves the GTF `gene_id` as: the `geneID` attribute, else
/// the `gene_id` attribute, else the GFF3 `Parent` link, else the transcript's
/// own ID.
fn gtf_gene_id(t: &GffObj) -> String {
    t.get_attr("geneID")
        .or_else(|| t.get_attr("gene_id"))
        .or_else(|| t.parent_id.as_deref())
        .unwrap_or(t.id.as_str())
        .to_string()
}

/// gclib `GffObj::geneID` member: only set from a `geneID`/`gene_id`
/// attribute or the GFF3 `Parent` link. Used for the `gene_id` value on
/// exon/CDS lines (no `getAttr` fallback there).
fn gtf_gene_id_member(t: &GffObj) -> Option<&str> {
    t.get_attr("geneID")
        .or_else(|| t.get_attr("gene_id"))
        .or_else(|| t.parent_id.as_deref())
}

/// gclib `GffObj::gene_name` member: first of the `gene_name`/`geneName`/
/// `gene_sym`/`gene` attributes found on the transcript.
fn gtf_gene_name_member(t: &GffObj) -> Option<&str> {
    ["gene_name", "geneName", "gene_sym", "gene"]
        .iter()
        .find_map(|k| t.get_attr(k))
}

/// Attributes that surface in native GTF output. gclib only stores parsed
/// attributes with `-F`/`--attrs`; without them, only the attributes that
/// gffread adds programmatically (-P / --adj-stop / --add-hasCDS) survive.
fn gtf_extra_attrs<'a>(cfg: &AppConfig, t: &'a GffObj) -> Vec<(&'a str, &'a str)> {
    if !cfg.attr_list.is_empty() {
        return t
            .attrs
            .iter()
            .filter(|(k, _)| cfg.attr_list.contains(k.as_str()))
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
    }
    if cfg.loader.full_attributes {
        return t.attrs.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
    }
    let mut out = Vec::new();
    for (k, v) in &t.attrs {
        let prog = (cfg.add_cds_attrs || cfg.adjust_stop)
            && matches!(k.as_str(), "partialness" | "InFrameStop" | "CDstopAdjusted");
        if prog || (cfg.add_has_cds && k == "hasCDS") {
            out.push((k.as_str(), v.as_str()));
        }
    }
    out
}

/// Attribute column for the GTF transcript line, mirroring gclib's
/// `GffObj::printGxf` GTF branch + `printAttrs`:
/// `transcript_id "X"; gene_id "Y"; gene_name "Z"; a1 "v1"; ...`
/// (separator `; `, NO trailing semicolon).
fn gtf_transcript_attrs(cfg: &AppConfig, t: &GffObj) -> String {
    let decode = cfg.decode_chars;
    let visible = gtf_extra_attrs(cfg, t);
    let mut parts: Vec<String> = Vec::new();
    parts.push(format!("transcript_id \"{}\"", t.id));
    let gid = gtf_gene_id(t);
    // gclib prints the gene_id/gene_name member values verbatim (no -D decode)
    parts.push(format!("gene_id \"{}\"", gid));
    // gene_name is printed from the member only when no gene_name/GENE_NAME
    // attribute is visible (it will then be emitted by the attrs loop below)
    if let Some(gn) = gtf_gene_name_member(t) {
        let shadowed = visible
            .iter()
            .any(|(k, _)| *k == "gene_name" || *k == "GENE_NAME");
        if !shadowed {
            parts.push(format!("gene_name \"{}\"", gn));
        }
    }
    // gclib consumes ID/Parent on GFF3 lines and transcript_id/gene_id on GTF
    // lines, so they never reach the printed attributes
    let gff3_input = t.get_attr("ID").is_some();
    let mut tr_id_seen = false;
    for (k, v) in visible {
        let mut name = k;
        if gff3_input {
            if name == "ID" || name == "Parent" {
                continue;
            }
        } else if name == "transcript_id" || name == "gene_id" {
            continue;
        }
        if name == "transcriptID" {
            if tr_id_seen {
                continue;
            }
            tr_id_seen = true;
        } else if name == "transcript_id" {
            if tr_id_seen {
                continue;
            }
            name = "transcriptID";
            tr_id_seen = true;
        }
        if name == "geneID" && v == gid {
            continue;
        }
        if name == "gene_id" {
            continue;
        }
        parts.push(format!("{} \"{}\"", name, fmt_attr_value(v, decode)));
    }
    parts.join("; ")
}

/// Attribute column for GTF exon/CDS lines, mirroring gclib `printGxfExon`:
/// `transcript_id "X"; gene_id "Y"; gene_name "Z";` (each attr ends in `;`).
fn gtf_exon_attrs(_cfg: &AppConfig, t: &GffObj) -> String {
    let mut s = format!("transcript_id \"{}\";", t.id);
    if let Some(gid) = gtf_gene_id_member(t) {
        s.push_str(&format!(" gene_id \"{}\";", gid));
    }
    if let Some(gn) = gtf_gene_name_member(t) {
        s.push_str(&format!(" gene_name \"{}\";", gn));
    }
    s
}

fn track_for(cfg: &AppConfig, t: &GffObj) -> String {
    cfg.track_label
        .clone()
        .unwrap_or_else(|| if t.track.is_empty() { ".".to_string() } else { t.track.clone() })
}

/// Print the GFF3 header.
pub fn print_gff3_header<W: Write>(out: &mut W, cfg: &AppConfig) -> std::io::Result<()> {
    if cfg.loader.keep_gff3_comments && !cfg.header_lines.is_empty() {
        for h in &cfg.header_lines {
            writeln!(out, "{}", h)?;
        }
    } else {
        writeln!(out, "##gff-version 3")?;
        writeln!(out, "# gffread v{}", env!("CARGO_PKG_VERSION"))?;
        if !cfg.cmd_line.is_empty() {
            writeln!(out, "# {}", cfg.cmd_line)?;
        }
    }
    Ok(())
}

pub fn print_seq_region<W: Write>(
    out: &mut W,
    cfg: &AppConfig,
    name: &str,
    start: u32,
    end: u32,
) -> std::io::Result<()> {
    if cfg.loader.keep_gff3_comments && start > 0 && end > 0 {
        writeln!(out, "##sequence-region {} {} {}", name, start, end)?;
    }
    Ok(())
}

/// Print a record in the configured format.
pub fn print_gxf<W: Write>(
    out: &mut W,
    t: &GffObj,
    cfg: &AppConfig,
    mode: GffPrintMode,
) -> std::io::Result<()> {
    match cfg.out_format {
        OutFormat::Gff3 => print_gff3(out, t, cfg, mode),
        OutFormat::Gtf => print_gtf(out, t, cfg),
        OutFormat::Bed => print_bed(out, t, cfg),
        OutFormat::Tlf => print_tlf(out, t, cfg),
        OutFormat::Table => print_table_data(out, t, cfg, false),
    }
}

fn print_gff3<W: Write>(
    out: &mut W,
    t: &GffObj,
    cfg: &AppConfig,
    mode: GffPrintMode,
) -> std::io::Result<()> {
    let track = track_for(cfg, t);
    let feat = if t.feature.is_empty() { "mRNA".to_string() } else { t.feature.clone() };
    // transcript/feature line
    let attrs = gff3_attrs(cfg, t, true, t.parent_id.is_some());
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        t.gseq_name, track, feat, t.start, t.end, t.score, t.strand, t.phase, attrs
    )?;

    // sub-features: print exons/CDS for transcript features in "Any" (default)
    // or "Both" mode. "Both" additionally prints them for non-transcript features.
    let print_subs = (mode == GffPrintMode::Any || mode == GffPrintMode::Both)
        && t.is_mrna;
    if print_subs {
        let print_exons = !t.exons.is_empty();
        if print_exons {
            for ex in t.exons.iter() {
                // Native gclib printGxfExon: GFF3 exon lines only carry
                // `Parent=<transcript_id>`; the exon ID is consumed and not
                // re-emitted in default (non -F) mode.
                let eattrs = format!("Parent={}", t.id);
                writeln!(
                    out,
                    "{}\t{}\texon\t{}\t{}\t.\t{}\t.\t{}",
                    t.gseq_name, track, ex.start, ex.end, t.strand, eattrs
                )?;
            }
        }
        if t.has_cds() {
            let cds_segs = cds_segments(t);
            let phases = compute_cds_phases(t, &cds_segs);
            for (i, c) in cds_segs.iter().enumerate() {
                let cattrs = format!("Parent={}", t.id);
                let phase = phases[i];
                writeln!(
                    out,
                    "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\t{}",
                    t.gseq_name, track, c.start, c.end, t.strand, phase, cattrs
                )?;
            }
        }
    }
    Ok(())
}

fn print_gtf<W: Write>(out: &mut W, t: &GffObj, cfg: &AppConfig) -> std::io::Result<()> {
    let source = if t.track.is_empty() { ".".to_string() } else { t.track.clone() };
    let score = if t.score.is_empty() { ".".to_string() } else { t.score.clone() };

    // transcript line: gclib hardcodes the feature type "transcript"
    writeln!(
        out,
        "{}\t{}\ttranscript\t{}\t{}\t{}\t{}\t.\t{}",
        t.gseq_name, source, t.start, t.end, score, t.strand,
        gtf_transcript_attrs(cfg, t)
    )?;

    // exons
    for ex in &t.exons {
        writeln!(
            out,
            "{}\t{}\texon\t{}\t{}\t.\t{}\t.\t{}",
            t.gseq_name, source, ex.start, ex.end, t.strand,
            gtf_exon_attrs(cfg, t)
        )?;
    }
    // CDS
    if t.has_cds() {
        let cds_segs = cds_segments(t);
        let phases = compute_gtf_cds_phases(t, &cds_segs);
        for (i, c) in cds_segs.iter().enumerate() {
            writeln!(
                out,
                "{}\t{}\tCDS\t{}\t{}\t.\t{}\t{}\t{}",
                t.gseq_name, source, c.start, c.end, t.strand, phases[i],
                gtf_exon_attrs(cfg, t)
            )?;
        }
    }
    Ok(())
}

fn print_bed<W: Write>(out: &mut W, t: &GffObj, _cfg: &AppConfig) -> std::io::Result<()> {
    let chrom = &t.gseq_name;
    let cstart = if t.start > 0 { t.start - 1 } else { 0 };
    let cend = t.end;
    let name = &t.id;
    let score = if t.score == "." { "0".to_string() } else { t.score.clone() };
    let strand = t.strand;
    let (tstart, tend) = if t.has_cds() {
        (if t.cds_start > 0 { t.cds_start - 1 } else { cstart }, t.cds_end)
    } else {
        (cstart, cstart)
    };
    let block_count = t.exons.len() as u32;
    let mut sizes = Vec::new();
    let mut starts = Vec::new();
    // BED exons must be sorted by start; exons are already sorted.
    // blockStarts are 0-based offsets from chromStart (also 0-based).
    for ex in &t.exons {
        sizes.push(format!("{}", ex.len()));
        // ex.start is 1-based; convert to 0-based then subtract cstart (0-based).
        starts.push(format!("{}", (ex.start as i64 - 1) - cstart as i64));
    }
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t{}\t{}",
        chrom,
        cstart,
        cend,
        name,
        score,
        strand,
        tstart,
        tend,
        block_count,
        sizes.join(","),
        starts.join(",")
    )
}

fn print_tlf<W: Write>(out: &mut W, t: &GffObj, cfg: &AppConfig) -> std::io::Result<()> {
    let track = track_for(cfg, t);
    let feat = if t.feature.is_empty() { "mRNA".to_string() } else { t.feature.clone() };
    let mut attrs = gff3_attrs(cfg, t, true, t.parent_id.is_some());
    // structural attributes
    let exon_str = t
        .exons
        .iter()
        .map(|e| format!("{}-{}", e.start, e.end))
        .collect::<Vec<_>>()
        .join(",");
    let attr_extra = format!("exoncount={};exons={}", t.exons.len(), exon_str);
    if t.has_cds() {
        let cds_str = format!(";CDSphase={};CDS={}:{}", phase_char(t.cds_phase), t.cds_start, t.cds_end);
        attrs.push_str(&format!(";{}{}", attr_extra, cds_str));
    } else {
        attrs.push_str(&format!(";{}", attr_extra));
    }
    writeln!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t{}",
        t.gseq_name, track, feat, t.start, t.end, t.score, t.strand, attrs
    )
}

fn phase_char(p: char) -> char {
    if p == '.' { '0' } else { p }
}

/// CDS segments in genomic coordinates, sorted by start.
fn cds_segments(t: &GffObj) -> Vec<GSeg> {
    if !t.has_cds() {
        return Vec::new();
    }
    let mut segs = Vec::new();
    for ex in &t.exons {
        if ex.end < t.cds_start || ex.start > t.cds_end {
            continue;
        }
        let s = ex.start.max(t.cds_start);
        let e = ex.end.min(t.cds_end);
        if e >= s {
            segs.push(GSeg::new(s, e));
        }
    }
    segs.sort_by_key(|s| s.start);
    segs
}

/// Compute the GFF3 phase column (one per CDS segment).
///
/// For `+` strand features, translation begins at the leftmost CDS segment;
/// for `-` strand features, translation begins at the rightmost segment.
/// The phase of the first translated segment is the record's stored phase
/// (defaulting to `0`). Subsequent segments have phase `(3 - L%3) % 3` where
/// `L` is the cumulative CDS length up to (and excluding) the current segment.
fn compute_cds_phases(t: &GffObj, segs: &[GSeg]) -> Vec<char> {
    let n = segs.len();
    if n == 0 {
        return Vec::new();
    }
    let mut out = vec!['0'; n];
    let initial = if t.cds_phase == '.' { '0' } else { t.cds_phase };
    if t.strand == '-' {
        // translation proceeds from the highest-coordinate segment down
        out[n - 1] = initial;
        let mut cum: u32 = 0;
        for i in (0..n.saturating_sub(1)).rev() {
            cum += segs[i + 1].len();
            let p = (3 - (cum % 3)) % 3;
            out[i] = char::from_digit(p, 10).unwrap_or('0');
        }
    } else {
        out[0] = initial;
        let mut cum: u32 = 0;
        for i in 1..n {
            cum += segs[i - 1].len();
            let p = (3 - (cum % 3)) % 3;
            out[i] = char::from_digit(p, 10).unwrap_or('0');
        }
    }
    out
}

/// CDS phase per segment for GTF output, mirroring gclib's
/// `GffObj::updateCDSPhase()`.
///
/// Unlike the GFF3 printer, gclib ALWAYS recomputes the per-segment phases for
/// the GTF output: translation begins at the leftmost CDS segment for `+`
/// strand features and at the rightmost for `-` strand features, starting with
/// the record's initial phase (`t.cds_phase`, the first/last CDS segment's
/// input phase). That initial phase consumes `3 - phase` bases of the first
/// codon, so the running counter is seeded with it; every segment then gets
/// phase `(3 - acc % 3) % 3` before its length is added to `acc`. (Verified
/// against native gffread for every strand/phase combination.)
fn compute_gtf_cds_phases(t: &GffObj, segs: &[GSeg]) -> Vec<char> {
    let n = segs.len();
    if n == 0 {
        return Vec::new();
    }
    let mut out = vec!['0'; n];
    let mut acc: u32 = match t.cds_phase {
        '1' => 2,
        '2' => 1,
        _ => 0,
    };
    if t.strand == '-' {
        for i in (0..n).rev() {
            out[i] = char::from_digit((3 - acc % 3) % 3, 10).unwrap_or('0');
            acc += segs[i].len();
        }
    } else {
        for i in 0..n {
            out[i] = char::from_digit((3 - acc % 3) % 3, 10).unwrap_or('0');
            acc += segs[i].len();
        }
    }
    out
}

/// Table output (one line per record).
pub fn print_table_data<W: Write>(
    out: &mut W,
    t: &GffObj,
    cfg: &AppConfig,
    in_fasta: bool,
) -> std::io::Result<()> {
    for (i, col) in cfg.table_cols.iter().enumerate() {
        if i > 0 || in_fasta {
            if !in_fasta || !matches!(col, TableField::Id) {
                out.write_all(b"\t")?;
            }
        }
        match col {
            TableField::Attr(name) => {
                let v = t.get_attr(name);
                if let Some(v) = v {
                    let s = if cfg.decode_chars { decode_hex_chars(v) } else { v.to_string() };
                    out.write_all(s.as_bytes())?;
                } else {
                    out.write_all(b".")?;
                }
            }
            TableField::Chr => out.write_all(t.gseq_name.as_bytes())?,
            TableField::Track => out.write_all(track_for(cfg, t).as_bytes())?,
            TableField::Id => {
                if !in_fasta {
                    out.write_all(t.id.as_bytes())?;
                }
            }
            TableField::GeneId => {
                let g = t.get_gene_id().unwrap_or(".");
                out.write_all(g.as_bytes())?;
            }
            TableField::GeneName => {
                let g = t.get_gene_name().unwrap_or(".");
                out.write_all(g.as_bytes())?;
            }
            TableField::Parent => {
                let p = t.parent_id.as_deref().unwrap_or(".");
                out.write_all(p.as_bytes())?;
            }
            TableField::Feature => out.write_all(t.feature.as_bytes())?,
            TableField::Start => write!(out, "{}", t.start)?,
            TableField::End => write!(out, "{}", t.end)?,
            TableField::Strand => out.write_all(&[t.strand as u8])?,
            TableField::NumExons => write!(out, "{}", t.exons.len())?,
            TableField::Exons => {
                if t.exons.is_empty() {
                    out.write_all(b".")?;
                } else {
                    for (k, e) in t.exons.iter().enumerate() {
                        if k > 0 {
                            out.write_all(b",")?;
                        }
                        write!(out, "{}-{}", e.start, e.end)?;
                    }
                }
            }
            TableField::Introns => {
                if t.exons.len() < 2 {
                    out.write_all(b".")?;
                } else {
                    for k in 0..t.exons.len() - 1 {
                        if k > 0 {
                            out.write_all(b",")?;
                        }
                        write!(out, "{}-{}", t.exons[k].end + 1, t.exons[k + 1].start - 1)?;
                    }
                }
            }
            TableField::Cds => {
                if t.has_cds() {
                    let segs = cds_segments(t);
                    for (k, c) in segs.iter().enumerate() {
                        if k > 0 {
                            out.write_all(b",")?;
                        }
                        write!(out, "{}-{}", c.start, c.end)?;
                    }
                } else {
                    out.write_all(b".")?;
                }
            }
            TableField::CovLen => write!(out, "{}", t.covlen)?,
            TableField::CdsLen => {
                if t.has_cds() {
                    let segs = cds_segments(t);
                    let l: u32 = segs.iter().map(|s| s.len()).sum();
                    write!(out, "{}", l)?;
                } else {
                    out.write_all(b"0")?;
                }
            }
            TableField::AllAttrs => {
                let s: Vec<String> = t.attrs.iter().map(|(k, v)| format!("{}={}", k, v)).collect();
                out.write_all(s.join(";").as_bytes())?;
            }
        }
    }
    if !in_fasta {
        writeln!(out)?;
    }
    Ok(())
}

/// Print a locus feature line (for -M output).
pub fn print_locus<W: Write>(
    out: &mut W,
    gseq_name: &str,
    loctrack: &str,
    locname: &str,
    start: u32,
    end: u32,
    strand: char,
    gene_names: &[String],
    gene_ids: &[String],
    rna_ids: &[String],
) -> std::io::Result<()> {
    write!(
        out,
        "{}\t{}\tlocus\t{}\t{}\t.\t{}\t.\tID={}",
        gseq_name, loctrack, start, end, strand, locname
    )?;
    if !gene_names.is_empty() {
        write!(out, ";genes={}", gene_names.join(","))?;
    }
    if !gene_ids.is_empty() {
        write!(out, ";geneIDs={}", gene_ids.join(","))?;
    }
    if !rna_ids.is_empty() {
        write!(out, ";transcripts={}", rna_ids.join(","))?;
    }
    writeln!(out)?;
    Ok(())
}
