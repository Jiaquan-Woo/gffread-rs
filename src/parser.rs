//! GFF3 / GTF / BED / TLF parser.
//!
//! Builds `GffObj` transcripts (with exons/CDS) from feature lines, mirroring
//! the assembly logic of gclib's `GffReader::readAll()`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use anyhow::{anyhow, Result};

use crate::config::AppConfig;
use crate::types::{
    feature_id, is_gene_type, is_subfeature_type, is_transcript_type, GenomicSeqData, GffObj, GSeg,
    GFF_FID_EXON,
};

/// Open a (possibly gzipped) text reader.
pub fn open_reader(path: &str) -> Result<Box<dyn BufRead>> {
    if path == "-" || path == "stdin" {
        let r = Box::new(BufReader::new(std::io::stdin()));
        return Ok(r);
    }
    let f = File::open(path).map_err(|e| anyhow!("Error: cannot open input file {}: {}", path, e))?;
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".gz") {
        let dec = flate2::read::GzDecoder::new(f);
        Ok(Box::new(BufReader::new(dec)))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn file_ext(path: &str) -> &str {
    let base = Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    // strip .gz
    let base = base.strip_suffix(".gz").unwrap_or(base);
    match base.rsplit_once('.') {
        Some((_, ext)) => ext,
        None => "",
    }
}

/// Parsed comment metadata.
#[derive(Debug, Default)]
pub struct Comments {
    pub header_lines: Vec<String>,
    /// gseq_name -> (start, end)
    pub seq_regions: HashMap<String, (u32, u32)>,
}

/// Parse a single 9-column feature line into raw fields.
struct RawLine {
    seqname: String,
    source: String,
    feature: String,
    start: u32,
    end: u32,
    score: String,
    strand: char,
    phase: char,
    attrs: Vec<(String, String)>,
}

fn split_line(line: &str) -> Option<RawLine> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 8 {
        return None;
    }
    let start: u32 = parts[3].parse().ok()?;
    let end: u32 = parts[4].parse().ok()?;
    let strand = parts[6].chars().next().unwrap_or('.');
    let phase = parts[7].chars().next().unwrap_or('.');
    let attrs = if parts.len() >= 9 {
        parse_attrs(parts[8])
    } else {
        Vec::new()
    };
    Some(RawLine {
        seqname: parts[0].to_string(),
        source: parts[1].to_string(),
        feature: parts[2].to_string(),
        start,
        end,
        score: parts[5].to_string(),
        strand,
        phase,
        attrs,
    })
}

/// Attribute parser handling both GFF3 (`k=v;k=v`) and GTF (`k "v"; k "v"`).
fn parse_attrs(s: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for piece in s.split(';') {
        let p = piece.trim();
        if p.is_empty() {
            continue;
        }
        if let Some(eq) = p.find('=') {
            let k = p[..eq].trim().to_string();
            let v = strip_quotes(p[eq + 1..].trim());
            if !k.is_empty() {
                out.push((k, v));
            }
        } else {
            // split on first whitespace
            let mut iter = p.splitn(2, char::is_whitespace);
            let k = iter.next().unwrap_or("").trim().to_string();
            let v = iter.next().unwrap_or("").trim();
            if !k.is_empty() {
                out.push((k, strip_quotes(v)));
            }
        }
    }
    out
}

fn strip_quotes(s: &str) -> String {
    let s = s.trim();
    if s.len() >= 2 && s.starts_with('"') && s.ends_with('"') {
        s[1..s.len() - 1].to_string()
    } else {
        s.to_string()
    }
}

fn attr_get<'a>(attrs: &'a [(String, String)], key: &str) -> Option<&'a str> {
    attrs.iter().find(|(k, _)| k == key).map(|(_, v)| v.as_str())
}

/// Merge exon list, optionally merging exons separated by tiny introns (-Z).
fn merge_close_exons(exons: &mut Vec<GSeg>, threshold: u32) {
    if exons.len() < 2 {
        return;
    }
    exons.sort_by_key(|e| e.start);
    let mut merged: Vec<GSeg> = Vec::with_capacity(exons.len());
    merged.push(exons[0]);
    for &e in &exons[1..] {
        let last = merged.last_mut().unwrap();
        if e.start <= last.end + 1 + threshold {
            if e.end > last.end {
                last.end = e.end;
            }
        } else {
            merged.push(e);
        }
    }
    *exons = merged;
}

/// Collapse overlapping / adjacent exon segments into a single interval per
/// cluster, mirroring gclib `GffObj::addExon()` which absorbs any exon/CDS/UTR
/// segment that overlaps or is adjacent to an existing exon (fully contained
/// segments simply vanish). Adjacent exons (gap == 0) are merged too.
fn merge_overlapping_exons(exons: &mut Vec<GSeg>) {
    if exons.len() < 2 {
        return;
    }
    exons.sort_by_key(|e| e.start);
    let mut merged: Vec<GSeg> = Vec::with_capacity(exons.len());
    merged.push(exons[0]);
    for &e in &exons[1..] {
        let last = merged.last_mut().unwrap();
        if e.start <= last.end + 1 {
            if e.end > last.end {
                last.end = e.end;
            }
        } else {
            merged.push(e);
        }
    }
    *exons = merged;
}

/// Parse BED12 lines into transcripts.
fn parse_bed_lines(cfg: &AppConfig, lines: impl BufRead, comments: &mut Comments) -> Result<Vec<GffObj>> {
    let _ = cfg;
    let mut rnas = Vec::new();
    for line in lines.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 6 {
            continue;
        }
        let chrom = cols[0].to_string();
        let cstart: u32 = cols[1].parse().unwrap_or(0); // 0-based
        let cend: u32 = cols[2].parse().unwrap_or(0);
        let name = cols.get(3).map(|s| s.to_string()).unwrap_or_default();
        let score = cols.get(4).map(|s| s.to_string()).unwrap_or_else(|| ".".to_string());
        let strand = cols.get(5).and_then(|s| s.chars().next()).unwrap_or('.');

        let mut t = GffObj::new();
        t.gseq_name = chrom.clone();
        t.track = "bed".to_string();
        t.feature = "mRNA".to_string();
        t.strand = strand;
        t.score = score;
        t.id = name;
        t.is_mrna = true;
        t.start = if cend > cstart { cstart + 1 } else { cstart + 1 };
        t.end = cend;
        t.phase = '.';

        // CDS from thickStart/thickEnd
        if cols.len() >= 8 {
            let tstart: u32 = cols[6].parse().unwrap_or(0);
            let tend: u32 = cols[7].parse().unwrap_or(0);
            if tend > tstart {
                t.cds_start = tstart + 1;
                t.cds_end = tend;
                t.cds_phase = '0';
            }
        }

        // blocks (exons)
        if cols.len() >= 12 {
            let block_count: u32 = cols[9].parse().unwrap_or(0);
            let sizes: Vec<&str> = cols[10].split(',').filter(|s| !s.is_empty()).collect();
            let starts: Vec<&str> = cols[11].split(',').filter(|s| !s.is_empty()).collect();
            for i in 0..block_count as usize {
                let sz: u32 = sizes.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
                let st: u32 = starts.get(i).and_then(|s| s.parse().ok()).unwrap_or(0);
                if sz == 0 {
                    continue;
                }
                let estart = cstart + st + 1;
                let eend = estart + sz - 1;
                t.exons.push(GSeg::new(estart, eend));
            }
        } else {
            // single block
            t.exons.push(GSeg::new(t.start, t.end));
        }
        t.finalize_exons();
        rnas.push(t);
    }
    let _ = comments;
    Ok(rnas)
}

/// Parse TLF (transcript-line-format): one GFF line per transcript with
/// `exoncount=N;exons=...;CDSphase=N;CDS=...` attributes.
fn parse_tlf_lines(cfg: &AppConfig, lines: impl BufRead, comments: &mut Comments) -> Result<Vec<GffObj>> {
    let _ = cfg;
    let mut rnas = Vec::new();
    for line in lines.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(raw) = split_line(line) else {
            continue;
        };
        let mut t = GffObj::new();
        t.gseq_name = raw.seqname;
        t.track = raw.source;
        t.feature = raw.feature;
        t.start = raw.start;
        t.end = raw.end;
        t.score = raw.score;
        t.strand = raw.strand;
        t.phase = raw.phase;
        t.is_mrna = true;
        if let Some(id) = attr_get(&raw.attrs, "ID") {
            t.id = id.to_string();
        }
        if let Some(pid) = attr_get(&raw.attrs, "Parent") {
            t.parent_id = Some(pid.to_string());
        }
        // exons
        if let Some(exstr) = attr_get(&raw.attrs, "exons") {
            for ep in exstr.split(',') {
                let ep = ep.trim();
                if let Some((a, b)) = ep.split_once('-') {
                    if let (Ok(s), Ok(e)) = (a.trim().parse::<u32>(), b.trim().parse::<u32>()) {
                        t.exons.push(GSeg::new(s, e));
                    }
                }
            }
        }
        if t.exons.is_empty() {
            t.exons.push(GSeg::new(t.start, t.end));
        }
        // CDS
        if let Some(cdsstr) = attr_get(&raw.attrs, "CDS") {
            // could be start:end or list like exons
            let parts: Vec<&str> = cdsstr.split(',').collect();
            let mut cstart = u32::MAX;
            let mut cend = 0u32;
            for p in parts {
                let p = p.trim();
                let (a, b) = if let Some((a, b)) = p.split_once('-') {
                    (a, b)
                } else if let Some((a, b)) = p.split_once(':') {
                    (a, b)
                } else {
                    continue;
                };
                if let (Ok(s), Ok(e)) = (a.trim().parse::<u32>(), b.trim().parse::<u32>()) {
                    if s < cstart {
                        cstart = s;
                    }
                    if e > cend {
                        cend = e;
                    }
                }
            }
            if cstart <= cend {
                t.cds_start = cstart;
                t.cds_end = cend;
            }
        }
        if let Some(ph) = attr_get(&raw.attrs, "CDSphase") {
            t.cds_phase = ph.chars().next().unwrap_or('.');
        }
        // copy remaining attrs except the structural ones
        for (k, v) in &raw.attrs {
            if matches!(k.as_str(), "exoncount" | "exons" | "CDS" | "CDSphase") {
                continue;
            }
            t.attrs.push((k.clone(), v.clone()));
        }
        t.finalize_exons();
        rnas.push(t);
    }
    let _ = comments;
    Ok(rnas)
}

/// Parse GFF3/GTF lines into transcripts + gene records.
fn parse_gxf_lines(
    cfg: &AppConfig,
    lines: impl BufRead,
    comments: &mut Comments,
) -> Result<(Vec<GffObj>, Vec<GffObj>)> {
    let mut rnas_by_id: HashMap<String, GffObj> = HashMap::new();
    let mut rna_order: Vec<String> = Vec::new();
    let mut genes: Vec<GffObj> = Vec::new();
    // sub-features awaiting parent resolution: (parent_id, exon_seg, is_cds, phase)
    let mut subs: Vec<(String, GSeg, bool, char)> = Vec::new();
    // standalone sub-features whose parent is a gene (for gene2exon / implicit transcripts)
    // raw feature lines that are neither transcript nor gene nor subfeature, kept for -O
    let mut others: Vec<GffObj> = Vec::new();

    let keep_all_attrs = cfg.loader.full_attributes;
    let gather_exon_attrs = cfg.loader.gather_exon_attrs;

    for line in lines.lines() {
        let line = line?;
        let line = line.trim_end_matches(['\r', '\n']);
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            // comment
            if line.starts_with("##") {
                let body = &line[2..];
                let mut tk = body.split_whitespace();
                if tk.next() == Some("sequence-region") {
                    let chr = tk.next().unwrap_or("");
                    let s: u32 = tk.next().and_then(|x| x.parse().ok()).unwrap_or(0);
                    let e: u32 = tk.next().and_then(|x| x.parse().ok()).unwrap_or(0);
                    if !chr.is_empty() {
                        comments.seq_regions.insert(chr.to_string(), (s, e));
                    }
                    continue;
                }
            }
            if cfg.loader.keep_gff3_comments {
                comments.header_lines.push(line.to_string());
            }
            continue;
        }

        let Some(raw) = split_line(line) else {
            continue;
        };

        let feat_l = raw.feature.to_ascii_lowercase();
        let is_sub = is_subfeature_type(&feat_l) || feat_l == "exon" || feat_l == "cds";
        let is_t = is_transcript_type(&feat_l);
        let is_g = is_gene_type(&feat_l);

        if is_sub {
            // exon / CDS / UTR subfeature
            let is_cds = feat_l == "cds";
            let seg = GSeg::new(raw.start, raw.end);
            // resolve parent(s): GTF lines link sub-features via transcript_id
            // (native gclib ignores a stray Parent attr in GTF mode); GFF3
            // lines use Parent.
            let parents: Vec<String> = if let Some(p) = attr_get(&raw.attrs, "transcript_id") {
                vec![p.to_string()]
            } else if let Some(p) = attr_get(&raw.attrs, "Parent") {
                p.split(',').map(|s| s.trim().to_string()).collect()
            } else {
                vec![]
            };
            for pid in parents {
                subs.push((pid, seg, is_cds, raw.phase));
            }
            continue;
        }

        if is_t || is_g {
            let mut o = GffObj::new();
            o.gseq_name = raw.seqname;
            o.track = raw.source;
            o.ftype_id = feature_id(&raw.feature);
            o.feature = raw.feature;
            o.start = raw.start;
            o.end = raw.end;
            o.score = raw.score;
            o.strand = raw.strand;
            o.phase = raw.phase;
            if is_t {
                o.is_mrna = true;
            }
            if is_g {
                o.is_gene = true;
            }
            // ID
            if let Some(id) = attr_get(&raw.attrs, "ID") {
                o.id = id.to_string();
            } else if let Some(id) = attr_get(&raw.attrs, "transcript_id") {
                o.id = id.to_string();
            }
            if let Some(p) = attr_get(&raw.attrs, "Parent") {
                o.parent_id = Some(p.to_string());
            }
            if let Some(gid) = attr_get(&raw.attrs, "gene_id").or_else(|| attr_get(&raw.attrs, "geneID")) {
                o.gene_id = Some(gid.to_string());
            }
            if let Some(gn) = attr_get(&raw.attrs, "gene_name").or_else(|| attr_get(&raw.attrs, "Name")) {
                if is_t {
                    o.gene_name = Some(gn.to_string());
                }
            }
            // copy attributes (full vs reduced handled later at print time)
            if keep_all_attrs || gather_exon_attrs || is_t {
                for (k, v) in &raw.attrs {
                    o.attrs.push((k.clone(), v.clone()));
                }
            } else {
                // minimal attributes
                for key in ["ID", "Name", "Parent", "gene_id", "gene_name", "geneID"] {
                    if let Some(v) = attr_get(&raw.attrs, key) {
                        o.attrs.push((key.to_string(), v.to_string()));
                    }
                }
            }

            if is_t {
                if o.id.is_empty() {
                    // no ID, synthesize
                    o.id = format!("{}_{}_{}", o.gseq_name, o.start, o.end);
                }
                if !rnas_by_id.contains_key(&o.id) {
                    rna_order.push(o.id.clone());
                    rnas_by_id.insert(o.id.clone(), o);
                } else {
                    // merge attrs into existing
                    let existing = rnas_by_id.get_mut(&o.id).unwrap();
                    for (k, v) in &o.attrs {
                        if existing.get_attr(k).is_none() {
                            existing.attrs.push((k.clone(), v.clone()));
                        }
                    }
                    if existing.gene_id.is_none() && o.gene_id.is_some() {
                        existing.gene_id = o.gene_id.clone();
                    }
                    if existing.gene_name.is_none() && o.gene_name.is_some() {
                        existing.gene_name = o.gene_name.clone();
                    }
                }
            } else {
                // gene
                if o.id.is_empty() {
                    o.id = format!("gene_{}_{}_{}", o.gseq_name, o.start, o.end);
                }
                genes.push(o);
            }
            continue;
        }

        // other feature types (only kept with -O)
        if !cfg.loader.transcripts_only {
            let mut o = GffObj::new();
            o.gseq_name = raw.seqname;
            o.track = raw.source;
            o.ftype_id = feature_id(&raw.feature);
            o.feature = raw.feature;
            o.start = raw.start;
            o.end = raw.end;
            o.score = raw.score;
            o.strand = raw.strand;
            o.phase = raw.phase;
            if let Some(id) = attr_get(&raw.attrs, "ID") {
                o.id = id.to_string();
            }
            if let Some(p) = attr_get(&raw.attrs, "Parent") {
                o.parent_id = Some(p.to_string());
            }
            for (k, v) in &raw.attrs {
                o.attrs.push((k.clone(), v.clone()));
            }
            others.push(o);
        }
    }

    // resolve sub-features into transcripts
    let mut implicit: HashMap<String, GffObj> = HashMap::new();
    for (pid, seg, is_cds, phase) in subs {
        if let Some(t) = rnas_by_id.get_mut(&pid) {
            // gclib stores every exon/CDS/UTR segment in the exon list and
            // merges overlapping or adjacent ones, so a UTR or CDS that is
            // fully contained in an exon never appears as a separate exon.
            t.exons.push(seg);
            if is_cds {
                if t.cds_start == 0 || seg.start < t.cds_start {
                    t.cds_start = seg.start;
                }
                if seg.end > t.cds_end {
                    t.cds_end = seg.end;
                }
                // gclib: CDphase = phase of the first CDS segment for `+`
                // strands (min start) and of the last one for `-` strands
                // (max start in sorted order, i.e. cdss->Last()->phase).
                // For `-` strand, only update phase when this segment has
                // the highest start coordinate seen so far.
                if t.strand == '-' {
                    if seg.start >= t.cds_max_start {
                        t.cds_max_start = seg.start;
                        t.cds_phase = phase;
                    }
                } else {
                    if seg.start == t.cds_start {
                        t.cds_phase = phase;
                    }
                }
            }
        } else {
            // parent is a gene or unknown -> create implicit transcript (unless parent is a gene and !gene2exon)
            let entry = implicit.entry(pid.clone()).or_insert_with(|| {
                let mut t = GffObj::new();
                t.id = pid.clone();
                t.is_mrna = true;
                t.feature = "mRNA".to_string();
                t
            });
            entry.exons.push(seg);
            if is_cds {
                if entry.cds_start == 0 || seg.start < entry.cds_start {
                    entry.cds_start = seg.start;
                }
                if seg.end > entry.cds_end {
                    entry.cds_end = seg.end;
                }
                // Strand is unknown here; track phases for both strands.
                // cds_phase = phase of min-start CDS (for `+` strand).
                // cds_phase_neg = phase of max-start CDS (for `-` strand).
                if seg.start == entry.cds_start {
                    entry.cds_phase = phase;
                }
                if seg.start >= entry.cds_max_start {
                    entry.cds_max_start = seg.start;
                    entry.cds_phase_neg = phase;
                }
            }
        }
    }

    // promote implicit transcripts
    for (id, mut t) in implicit {
        if rnas_by_id.contains_key(&id) {
            continue;
        }
        // try to inherit gene info from a gene record matching this id as Parent
        // (i.e. the id was actually a gene id and exons parented to gene)
        if let Some(g) = genes.iter().find(|g| g.id == id) {
            t.gseq_name = g.gseq_name.clone();
            t.track = g.track.clone();
            t.strand = g.strand;
            t.start = g.start;
            t.end = g.end;
            t.gene_id = Some(g.id.clone());
            t.parent_id = Some(g.id.clone());
            if let Some(n) = g.get_attr("Name") {
                t.gene_name = Some(n.to_string());
            }
            // Strand now known: for '-' strand, use the max-start CDS phase.
            if t.strand == '-' && t.cds_phase_neg != '.' {
                t.cds_phase = t.cds_phase_neg;
            }
        }
        if t.gseq_name.is_empty() {
            // no info; skip
            continue;
        }
        t.finalize_exons();
        rna_order.push(id.clone());
        rnas_by_id.insert(id, t);
    }

    // finalize transcripts
    let mut rnas: Vec<GffObj> = Vec::with_capacity(rna_order.len());
    for id in rna_order {
        if let Some(mut t) = rnas_by_id.remove(&id) {
            // gclib merges overlapping/adjacent exon/CDS/UTR segments while
            // parsing; replicate that before any downstream use.
            merge_overlapping_exons(&mut t.exons);
            if t.exons.is_empty() {
                // transcript line with no exons; treat whole span as one exon if gene2exon or force
                if cfg.loader.gene2exon || cfg.loader.force_exons {
                    t.exons.push(GSeg::new(t.start, t.end));
                }
            }
            t.finalize_exons();
            if cfg.loader.merge_close_exons {
                merge_close_exons(&mut t.exons, 3);
                t.finalize_exons();
            }
            if cfg.loader.force_exons {
                t.subftype_id = GFF_FID_EXON;
            }
            rnas.push(t);
        }
    }

    // attach others as gfs too
    let mut gfs = genes;
    gfs.extend(others);

    Ok((rnas, gfs))
}

/// Load a single input file into per-reference `GenomicSeqData`.
pub fn load_file(cfg: &mut AppConfig, path: &str) -> Result<Vec<GenomicSeqData>> {
    let mut comments = Comments::default();
    let bed = cfg.loader.bed_input || file_ext(path).eq_ignore_ascii_case("bed");
    let tlf = cfg.loader.tlf_input || file_ext(path).eq_ignore_ascii_case("tlf");

    let (rnas, gfs) = if bed {
        let reader = open_reader(path)?;
        (parse_bed_lines(cfg, reader, &mut comments)?, Vec::new())
    } else if tlf {
        let reader = open_reader(path)?;
        (parse_tlf_lines(cfg, reader, &mut comments)?, Vec::new())
    } else {
        let reader = open_reader(path)?;
        parse_gxf_lines(cfg, reader, &mut comments)?
    };

    // store header lines into config (first file only)
    if cfg.header_lines.is_empty() {
        cfg.header_lines = comments.header_lines.clone();
    }

    // group by genomic sequence
    let mut by_gseq: HashMap<String, GenomicSeqData> = HashMap::new();
    let mut gseq_order: Vec<String> = Vec::new();

    for mut t in rnas {
        // apply ref name mapping
        if !cfg.ref_map.is_empty() {
            if let Some(new) = cfg.ref_map.get(&t.gseq_name) {
                t.gseq_name = new.clone();
            } else {
                // try stripping version suffix (.N)
                if t.gseq_name.len() > 2 {
                    let bytes = t.gseq_name.as_bytes();
                    if bytes[bytes.len() - 2] == b'.' && bytes[bytes.len() - 1].is_ascii_digit() {
                        let base = &t.gseq_name[..t.gseq_name.len() - 2];
                        if let Some(new) = cfg.ref_map.get(base) {
                            t.gseq_name = new.clone();
                        }
                    }
                }
            }
        }
        let gname = t.gseq_name.clone();
        if !by_gseq.contains_key(&gname) {
            gseq_order.push(gname.clone());
            let gid = cfg.names.add_gseq(&gname);
            let mut gd = GenomicSeqData::new(gid, gname.clone());
            if let Some((s, e)) = comments.seq_regions.get(&gname) {
                gd.seqreg_start = *s;
                gd.seqreg_end = *e;
            }
            by_gseq.insert(gname.clone(), gd);
        }
        let gd = by_gseq.get_mut(&gname).unwrap();
        t.gseq_id = gd.gseq_id;
        gd.rnas.push(t);
    }

    for mut g in gfs {
        if !cfg.ref_map.is_empty() {
            if let Some(new) = cfg.ref_map.get(&g.gseq_name) {
                g.gseq_name = new.clone();
            }
        }
        let gname = g.gseq_name.clone();
        if !by_gseq.contains_key(&gname) {
            gseq_order.push(gname.clone());
            let gid = cfg.names.add_gseq(&gname);
            let gd = GenomicSeqData::new(gid, gname.clone());
            by_gseq.insert(gname.clone(), gd);
        }
        let gd = by_gseq.get_mut(&gname).unwrap();
        g.gseq_id = gd.gseq_id;
        gd.gfs.push(g);
    }

    // sort rnas/gfs by start coordinate
    for gd in by_gseq.values_mut() {
        gd.rnas.sort_by_key(|t| t.start);
        gd.gfs.sort_by_key(|t| t.start);
    }

    let mut result: Vec<GenomicSeqData> = gseq_order
        .into_iter()
        .map(|n| by_gseq.remove(&n).unwrap())
        .collect();

    if cfg.loader.sort_refs_alpha {
        result.sort_by(|a, b| a.gseq_name.cmp(&b.gseq_name));
    } else if let Some(lst) = &cfg.sort_by.clone() {
        // order by the ref list
        let order: Vec<String> = load_ref_list(lst).unwrap_or_default();
        let idx = |name: &str| -> usize {
            order.iter().position(|x| x == name).unwrap_or(usize::MAX)
        };
        result.sort_by_key(|gd| idx(&gd.gseq_name));
    }

    Ok(result)
}

/// Load a whitespace/comma delimited reference name list.
fn load_ref_list(path: &str) -> Result<Vec<String>> {
    let mut out = Vec::new();
    let reader = open_reader(path)?;
    for line in reader.lines() {
        let line = line?;
        for tok in line.split(|c: char| c.is_whitespace() || c == ',' || c == ';') {
            let t = tok.trim();
            if !t.is_empty() {
                out.push(t.to_string());
            }
        }
    }
    Ok(out)
}
