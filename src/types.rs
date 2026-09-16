//! Core data types for the gffread Rust rewrite.
//!
//! These mirror the gclib structures used by the original C++ implementation
//! (`GffObj`, `GSeg`, `GffLocus`, `GenomicSeqData`, etc.) but are owned Rust
//! values without raw pointers.

use std::collections::HashMap;

/// Genomic segment / interval (1-based, inclusive coordinates).
#[derive(Debug, Clone, Copy)]
pub struct GSeg {
    pub start: u32,
    pub end: u32,
}

impl GSeg {
    pub fn new(start: u32, end: u32) -> Self {
        GSeg { start, end }
    }
    pub fn len(&self) -> u32 {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }
    pub fn is_empty(&self) -> bool {
        self.end < self.start
    }
    pub fn overlap(&self, other: &GSeg) -> bool {
        self.start <= other.end && self.end >= other.start
    }
    pub fn overlap_len(&self, other: &GSeg) -> u32 {
        let s = self.start.max(other.start);
        let e = self.end.min(other.end);
        if s <= e {
            e - s + 1
        } else {
            0
        }
    }
    pub fn overlap_len_coords(&self, ostart: u32, oend: u32) -> u32 {
        let s = self.start.max(ostart);
        let e = self.end.min(oend);
        if s <= e {
            e - s + 1
        } else {
            0
        }
    }
}

impl Default for GSeg {
    fn default() -> Self {
        GSeg { start: 0, end: 0 }
    }
}

/// Exon coordinate alias.
pub type GffExon = GSeg;

/// udata bit flags (match the C++ T_PRINTABLE / T_NO_PRINT / T_DUPSHOW macros).
pub const UDATA_PRINTABLE: u32 = 0x100;
pub const UDATA_DUPSHOW: u32 = 0x200;
pub const UDATA_OSTRAND_MASK: u32 = 0xFF;

/// A single GFF/GTF/BED record. Transcript features carry exon/CDS info.
#[derive(Debug, Clone)]
pub struct GffObj {
    /// Reference sequence name (column 1).
    pub gseq_name: String,
    /// Source / method / track (column 2).
    pub track: String,
    /// Feature type (column 3), e.g. "gene", "mRNA", "exon", "CDS".
    pub feature: String,
    /// Start coordinate (1-based, inclusive).
    pub start: u32,
    /// End coordinate (1-based, inclusive).
    pub end: u32,
    /// Score (column 6, kept as string).
    pub score: String,
    /// Strand: '+', '-' or '.'.
    pub strand: char,
    /// Phase (column 8) for CDS: '0','1','2' or '.'.
    pub phase: char,

    /// Feature ID (ID attribute / transcript_id).
    pub id: String,
    /// Parent ID, if any.
    pub parent_id: Option<String>,
    /// Gene ID.
    pub gene_id: Option<String>,
    /// Gene name.
    pub gene_name: Option<String>,
    /// All attributes, preserved in input order.
    pub attrs: Vec<(String, String)>,

    /// Exons (sorted, for transcript features).
    pub exons: Vec<GSeg>,
    /// CDS start (CDstart in C++).
    pub cds_start: u32,
    /// CDS end (CDend in C++).
    pub cds_end: u32,
    /// CDS phase (CDphase in C++).
    pub cds_phase: char,
    /// Max start coordinate among CDS segments (for '-' strand phase resolution).
    pub cds_max_start: u32,
    /// Phase of max-start CDS segment (used when strand is unknown during parsing).
    pub cds_phase_neg: char,
    /// Covered length (sum of exon lengths).
    pub covlen: u32,

    /// Whether this is a transcript feature.
    pub is_mrna: bool,
    /// Whether this is a gene feature.
    pub is_gene: bool,
    /// Numeric genomic-sequence id (index into name table).
    pub gseq_id: i32,
    /// Numeric track id.
    pub track_id: i32,
    /// Feature-type id.
    pub ftype_id: i32,
    /// Sub-feature type id.
    pub subftype_id: i32,
    /// User flags (printability, original-strand storage for -Y, etc.).
    pub udata: u32,
}

impl GffObj {
    pub fn new() -> Self {
        GffObj {
            gseq_name: String::new(),
            track: String::new(),
            feature: String::new(),
            start: 0,
            end: 0,
            score: ".".to_string(),
            strand: '.',
            phase: '.',
            id: String::new(),
            parent_id: None,
            gene_id: None,
            gene_name: None,
            attrs: Vec::new(),
            exons: Vec::new(),
            cds_start: 0,
            cds_end: 0,
            cds_phase: '.',
            cds_max_start: 0,
            cds_phase_neg: '.',
            covlen: 0,
            is_mrna: false,
            is_gene: false,
            gseq_id: 0,
            track_id: 0,
            ftype_id: 0,
            subftype_id: 0,
            udata: 0,
        }
    }

    /// Linear attribute lookup.
    pub fn get_attr(&self, key: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    pub fn add_attr(&mut self, key: &str, value: &str) {
        if let Some(slot) = self.attrs.iter_mut().find(|(k, _)| k == key) {
            slot.1 = value.to_string();
        } else {
            self.attrs.push((key.to_string(), value.to_string()));
        }
    }

    pub fn remove_attr(&mut self, key: &str) {
        self.attrs.retain(|(k, _)| k != key);
    }

    pub fn has_cds(&self) -> bool {
        self.cds_start > 0 && self.cds_end > 0
    }

    pub fn is_transcript(&self) -> bool {
        self.is_mrna
    }

    pub fn is_printable(&self) -> bool {
        (self.udata & UDATA_PRINTABLE) == 0
    }

    pub fn no_print(&mut self) {
        self.udata |= UDATA_PRINTABLE;
    }

    pub fn is_dupshowable(&self) -> bool {
        (self.udata & UDATA_DUPSHOW) == 0
    }

    pub fn no_dupshow(&mut self) {
        self.udata |= UDATA_DUPSHOW;
    }

    /// Stored original strand (low byte of udata), used by -Y.
    pub fn orig_strand(&self) -> char {
        let b = (self.udata & UDATA_OSTRAND_MASK) as u8;
        if b == 0 {
            0 as char
        } else {
            b as char
        }
    }

    pub fn set_orig_strand(&mut self, s: char) {
        self.udata = (self.udata & !UDATA_OSTRAND_MASK) | (s as u32 & UDATA_OSTRAND_MASK);
    }

    pub fn get_gene_id(&self) -> Option<&str> {
        self.gene_id.as_deref().or_else(|| self.get_attr("geneID"))
    }

    pub fn get_gene_name(&self) -> Option<&str> {
        self.gene_name
            .as_deref()
            .or_else(|| self.get_attr("gene_name"))
            .or_else(|| self.get_attr("Name"))
    }

    /// Recompute covlen/start/end from the exon list.
    pub fn finalize_exons(&mut self) {
        if self.exons.is_empty() {
            return;
        }
        self.exons.sort_by_key(|e| e.start);
        // merge identical/adjacent is left to callers when needed
        let mut cov = 0u32;
        let mut s = u32::MAX;
        let mut e = 0u32;
        for ex in &self.exons {
            cov += ex.len();
            if ex.start < s {
                s = ex.start;
            }
            if ex.end > e {
                e = ex.end;
            }
        }
        self.covlen = cov;
        if self.start == 0 {
            self.start = s;
        } else {
            self.start = self.start.min(s);
        }
        self.end = self.end.max(e);
    }

    /// CDS segments projected onto the spliced transcript coordinates.
    pub fn get_cds_segs(&self) -> Vec<GSeg> {
        let mut cds = Vec::new();
        if !self.has_cds() {
            return cds;
        }
        let mut offset = 0u32;
        for exon in &self.exons {
            let exon_len = exon.len();
            let exon_start = offset + 1;
            if exon.end < self.cds_start || exon.start > self.cds_end {
                offset += exon_len;
                continue;
            }
            let cds_in_exon_start = if exon.start < self.cds_start {
                exon_start + (self.cds_start - exon.start)
            } else {
                exon_start
            };
            let cds_in_exon_end = if exon.end > self.cds_end {
                exon_start + (self.cds_end - exon.start)
            } else {
                exon_start + exon_len - 1
            };
            if cds_in_exon_start <= cds_in_exon_end {
                cds.push(GSeg::new(cds_in_exon_start, cds_in_exon_end));
            }
            offset += exon_len;
        }
        cds
    }
}

impl Default for GffObj {
    fn default() -> Self {
        GffObj::new()
    }
}

/// Nucleotide complement.
pub fn nt_complement(c: char) -> char {
    match c {
        'A' | 'a' => 'T',
        'T' | 't' => 'A',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        'U' | 'u' => 'A',
        'N' | 'n' => 'N',
        // IUPAC ambiguity codes
        'R' | 'r' => 'Y',
        'Y' | 'y' => 'R',
        'S' | 's' => 'S',
        'W' | 'w' => 'W',
        'K' | 'k' => 'M',
        'M' | 'm' => 'K',
        'B' | 'b' => 'V',
        'V' | 'v' => 'B',
        'D' | 'd' => 'H',
        'H' | 'h' => 'D',
        _ => 'N',
    }
}

/// Reverse-complement a byte sequence in place semantics (returns new Vec).
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| nt_complement(b as char) as u8)
        .collect()
}

/// Splice site (donor/acceptor dinucleotide).
#[derive(Debug, Clone, Copy)]
pub struct GSpliceSite {
    pub nt: [char; 2],
}

impl GSpliceSite {
    pub fn new(c1: char, c2: char) -> Self {
        GSpliceSite {
            nt: [c1.to_ascii_uppercase(), c2.to_ascii_uppercase()],
        }
    }

    /// Build from an intron sequence. `get_acceptor` selects the acceptor
    /// (3') site instead of the donor (5') site. `revc` reverse-complements
    /// (for transcripts on the minus strand).
    pub fn from_intron(intron: &[u8], get_acceptor: bool, revc: bool) -> Self {
        let mut nt = ['N', 'N'];
        if intron.len() < 2 {
            return GSpliceSite { nt };
        }
        if revc {
            if get_acceptor {
                nt[0] = nt_complement(intron[0] as char);
                nt[1] = nt_complement(intron[1] as char);
            } else {
                let idx = intron.len() - 2;
                nt[0] = nt_complement(intron[idx] as char);
                nt[1] = nt_complement(intron[idx + 1] as char);
            }
        } else if get_acceptor {
            let idx = intron.len() - 2;
            nt[0] = (intron[idx] as char).to_ascii_uppercase();
            nt[1] = (intron[idx + 1] as char).to_ascii_uppercase();
        } else {
            nt[0] = (intron[0] as char).to_ascii_uppercase();
            nt[1] = (intron[1] as char).to_ascii_uppercase();
        }
        GSpliceSite { nt }
    }

    pub fn canonical_donor(&self) -> bool {
        self.nt[0] == 'G' && (self.nt[1] == 'T' || self.nt[1] == 'C')
    }

    pub fn as_str(&self) -> [u8; 2] {
        [self.nt[0] as u8, self.nt[1] as u8]
    }
}

impl PartialEq for GSpliceSite {
    fn eq(&self, other: &Self) -> bool {
        self.nt == other.nt
    }
}

/// Intron record for `-j` junction output.
#[derive(Debug, Clone)]
pub struct CIntronData {
    pub start: u32,
    pub end: u32,
    pub strand: char,
    pub ts: Vec<String>,
}

impl CIntronData {
    pub fn new(start: u32, end: u32, strand: char, t_id: &str) -> Self {
        CIntronData {
            start,
            end,
            strand,
            ts: vec![t_id.to_string()],
        }
    }
    pub fn add(&mut self, t_id: &str) {
        self.ts.push(t_id.to_string());
    }
}

/// Intron list for `-j` output (per genomic sequence).
#[derive(Debug, Clone, Default)]
pub struct CIntronList {
    pub gseq_id: i32,
    pub last_t_start: u32,
    pub introns: Vec<CIntronData>,
}

impl CIntronList {
    pub fn new() -> Self {
        CIntronList::default()
    }

    /// Add all introns of transcript `t` (id, strand, exons).
    pub fn add(&mut self, t_id: &str, strand: char, exons: &[GSeg]) {
        if exons.len() < 2 {
            return;
        }
        for i in 1..exons.len() {
            let istart = exons[i - 1].end + 1;
            let iend = if exons[i].start > 0 {
                exons[i].start - 1
            } else {
                0
            };
            if let Some(slot) = self.introns.iter_mut().find(|d| {
                d.start == istart && d.end == iend && d.strand == strand
            }) {
                slot.add(t_id);
            } else {
                self.introns.push(CIntronData::new(istart, iend, strand, t_id));
            }
            self.last_t_start = exons[0].start;
        }
    }

    pub fn print_to(&self, out: &mut dyn std::io::Write, gseqname: &str) -> std::io::Result<()> {
        let mut introns = self.introns.clone();
        introns.sort_by(|a, b| a.start.cmp(&b.start).then(a.end.cmp(&b.end)));
        for idata in &introns {
            let mut ts = idata.ts.clone();
            ts.sort();
            let joined = ts.join(",");
            writeln!(out, "{}\t{}\t{}\t{}\t{}", gseqname, idata.start, idata.end, idata.strand, joined)?;
        }
        Ok(())
    }
}

/// Sequence info from the `-s` option.
#[derive(Debug, Clone)]
pub struct SeqInfo {
    pub length: u32,
    pub description: Option<String>,
}

/// Reference name mapping entry for `-m`.
#[derive(Debug, Clone)]
pub struct RefTran {
    pub new_name: String,
}

/// Range filter (`-r` / `--jmatch`).
#[derive(Debug, Clone)]
pub struct RangeFilter {
    pub ref_name: Option<String>,
    pub strand: char,
    pub start: u32,
    pub end: u32,
}

impl RangeFilter {
    pub fn new() -> Self {
        RangeFilter {
            ref_name: None,
            strand: 0 as char,
            start: 0,
            end: u32::MAX,
        }
    }
}

/// ID filter mode.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum IDFltType {
    None,
    Only,
    Exclude,
}

/// GFF/GTF/BED print mode (exon printing style).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GffPrintMode {
    Any,
    Both,
    Gtf,
    Bed,
    Tlf,
}

/// Table output column descriptor.
#[derive(Debug, Clone)]
pub enum TableField {
    Attr(String),
    Chr,
    Track,
    Id,
    GeneId,
    GeneName,
    Parent,
    Feature,
    Start,
    End,
    Strand,
    NumExons,
    Exons,
    Introns,
    Cds,
    CovLen,
    CdsLen,
    AllAttrs,
}

/// Genomic-sequence aggregate (all records on one reference).
#[derive(Debug, Clone)]
pub struct GenomicSeqData {
    pub gseq_id: i32,
    pub gseq_name: String,
    pub gfs: Vec<GffObj>, // gene / non-transcript features
    pub rnas: Vec<GffObj>, // transcripts
    pub loci: Vec<GffLocus>,
    pub seqreg_start: u32,
    pub seqreg_end: u32,
    pub f_bases: u64,
    pub r_bases: u64,
    pub u_bases: u64,
}

impl GenomicSeqData {
    pub fn new(gseq_id: i32, gseq_name: String) -> Self {
        GenomicSeqData {
            gseq_id,
            gseq_name,
            gfs: Vec::new(),
            rnas: Vec::new(),
            loci: Vec::new(),
            seqreg_start: 0,
            seqreg_end: 0,
            f_bases: 0,
            r_bases: 0,
            u_bases: 0,
        }
    }
}

/// A genomic locus (cluster of transcripts / genes), for `-M`.
#[derive(Debug, Clone)]
pub struct GffLocus {
    pub gseq_id: i32,
    pub locus_num: u32,
    pub strand: char,
    pub is_mrna: bool,
    /// index into the locus rnas of the max-coverage transcript.
    pub t_maxcov: Option<usize>,
    pub gfs: Vec<usize>, // indices into the parent gdata.gfs (not used heavily here)
    pub rna_ids: Vec<String>,
    pub rna_covlens: Vec<u32>,
    pub mexons: Vec<GSeg>,
    pub gene_names: Vec<String>,
    pub gene_ids: Vec<String>,
    pub start: u32,
    pub end: u32,
}

impl GffLocus {
    pub fn new() -> Self {
        GffLocus {
            gseq_id: 0,
            locus_num: 0,
            strand: '.',
            is_mrna: false,
            t_maxcov: None,
            gfs: Vec::new(),
            rna_ids: Vec::new(),
            rna_covlens: Vec::new(),
            mexons: Vec::new(),
            gene_names: Vec::new(),
            gene_ids: Vec::new(),
            start: 0,
            end: 0,
        }
    }
}

/// Global name table (mirrors GffNames in gclib).
#[derive(Debug, Default, Clone)]
pub struct GffNames {
    pub gseq_names: Vec<String>,
    pub gseq_index: HashMap<String, i32>,
    pub feat_names: Vec<String>,
    pub attr_names: Vec<String>,
}

impl GffNames {
    pub fn add_gseq(&mut self, name: &str) -> i32 {
        if let Some(&i) = self.gseq_index.get(name) {
            return i;
        }
        let i = self.gseq_names.len() as i32;
        self.gseq_names.push(name.to_string());
        self.gseq_index.insert(name.to_string(), i);
        i
    }
    pub fn gseq_name(&self, id: i32) -> Option<&str> {
        self.gseq_names.get(id as usize).map(|s| s.as_str())
    }
}

/// Maximum locus span when grouping transcripts (`GFF_MAX_LOCUS`).
pub const GFF_MAX_LOCUS: u32 = 100_000;

/// Feature-type id constants (mirrors gclib).
pub const GFF_FID_MRNA: i32 = 1;
pub const GFF_FID_EXON: i32 = 2;
pub const GFF_FID_CDS: i32 = 3;
pub const GFF_FID_GENE: i32 = 4;
pub const GFF_FID_RNA: i32 = 5;
pub const GFF_FID_MISC: i32 = 6;

/// Classify a feature name into a coarse id.
pub fn feature_id(name: &str) -> i32 {
    let l = name.to_ascii_lowercase();
    match l.as_str() {
        "mrna" | "transcript" => GFF_FID_MRNA,
        "exon" => GFF_FID_EXON,
        "cds" => GFF_FID_CDS,
        "gene" | "pseudogene" => GFF_FID_GENE,
        _ => {
            if l.ends_with("gene") {
                GFF_FID_GENE
            } else if l.contains("rna") || l == "transcript" {
                GFF_FID_RNA
            } else {
                GFF_FID_MISC
            }
        }
    }
}

/// Whether a feature type is a transcript/RNA type.
pub fn is_transcript_type(name: &str) -> bool {
    let l = name.to_ascii_lowercase();
    matches!(
        l.as_str(),
        "mrna" | "transcript" | "trna" | "rrna" | "mirna"
        | "snorna" | "snrna" | "ncrna" | "lncrna" | "sirna"
        | "pirna" | "tmrna" | "rnase_mrp_rna" | "srp_rna"
        | "misc_rna" | "primary_transcript"
        | "processed_transcript" | "pseudogene_transcript"
        | "pseudogenic_transcript" | "pseudogenic_trna"
    ) || l.ends_with("rna")
        || (l.ends_with("transcript") && l != "transcript_region")
}

/// Whether a feature type is a gene/locus container.
pub fn is_gene_type(name: &str) -> bool {
    let l = name.to_ascii_lowercase();
    l == "gene" || l.ends_with("gene") || l == "locus"
}

/// Whether a feature type is a sub-feature (exon/CDS/UTR/...).
pub fn is_subfeature_type(name: &str) -> bool {
    let l = name.to_ascii_lowercase();
    matches!(
        l.as_str(),
        "exon" | "cds" | "utr" | "three_prime_utr" | "five_prime_utr"
        | "3utr" | "5utr" | "intron" | "sig_peptide" | "tss" | "tts"
        | "start_codon" | "stop_codon"
    ) || l.ends_with("_utr")
        || l == "utr"
        || l.ends_with("exon")
}

/// Merge two sorted exon lists into a non-overlapping merged list (used for
/// locus mexons).
pub fn merge_segments(segs: &mut Vec<GSeg>) {
    if segs.len() < 2 {
        return;
    }
    segs.sort_by_key(|s| s.start);
    let mut out: Vec<GSeg> = Vec::with_capacity(segs.len());
    out.push(segs[0]);
    for &s in &segs[1..] {
        let last = out.last_mut().unwrap();
        if s.start <= last.end + 1 {
            if s.end > last.end {
                last.end = s.end;
            }
        } else {
            out.push(s);
        }
    }
    *segs = out;
}
