//! gffread - filter, convert or cluster GFF/GTF/BED records and extract
//! transcript sequences. Rust rewrite of gffread v0.12.9.

mod cluster;
mod config;
mod fasta;
mod filter;
mod output;
mod parser;
mod translate;
mod types;

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use anyhow::{anyhow, bail, Result};
use clap::Parser;

use cluster::{cluster_gdata, restore_strands};
use config::{AppConfig, LoaderOpts, OutFormat, SeqInfoEntry};
use fasta::FastaDb;
use filter::{check_filters_mut, process_transcript, Outputs};
use parser::load_file;
use types::{GffObj, GffPrintMode, IDFltType, RangeFilter, TableField};

const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Parser, Debug)]
#[command(
    name = "gffread",
    version = VERSION,
    disable_version_flag = true,
    about = "Filter, convert or cluster GFF/GTF/BED records, extract transcript sequences."
)]
struct Args {
    /// Input GFF/GTF/BED file(s). stdin if "-" or omitted.
    #[arg(name = "input")]
    inputs: Vec<String>,

    // ---- output format ----
    #[arg(short = 'T', long = "gtf", help = "Output GTF instead of GFF3")]
    gtf: bool,
    #[arg(long = "bed", help = "Output BED format")]
    bed: bool,
    #[arg(long = "tlf", help = "Output TLF (transcript line format)")]
    tlf: bool,
    #[arg(long = "table", value_name = "ATTRS", help = "Tab-delimited table output")]
    table: Option<String>,

    // ---- output files ----
    #[arg(short = 'o', value_name = "FILE", help = "Write output records to FILE")]
    out: Option<String>,
    #[arg(short = 'w', value_name = "FILE", help = "Write spliced exons (transcripts) FASTA")]
    w: Option<String>,
    #[arg(short = 'x', value_name = "FILE", help = "Write spliced CDS FASTA")]
    x: Option<String>,
    #[arg(short = 'y', value_name = "FILE", help = "Write translated protein FASTA")]
    y: Option<String>,
    #[arg(short = 'u', value_name = "FILE", help = "Write unspliced transcript FASTA")]
    u: Option<String>,
    #[arg(short = 'j', value_name = "FILE", help = "Write junctions (introns)")]
    j: Option<String>,
    #[arg(short = 'd', value_name = "FILE", help = "Write duplicate info (for -M)")]
    d: Option<String>,

    // ---- FASTA / sequence ----
    #[arg(short = 'g', value_name = "FASTA|DIR", help = "Genomic sequences FASTA or directory")]
    g: Option<String>,
    #[arg(long = "w-add", value_name = "N", help = "Additional bases around -w/-u boundaries")]
    w_add: Option<i64>,
    #[arg(long = "w-nocds", help = "For -w, do not output CDS info")]
    w_nocds: bool,
    #[arg(short = 'W', help = "Write exon coordinates projected onto spliced sequence")]
    w_proj: bool,
    #[arg(short = 'S', help = "Use '*' instead of '.' for stop codons")]
    star_stop: bool,
    #[arg(short = 'P', help = "Add coding-status GFF attributes (requires -g)")]
    add_cds_attrs: bool,
    #[arg(long = "add-hasCDS", help = "Add a hasCDS=true attribute")]
    add_hascds: bool,
    #[arg(long = "adj-stop", help = "Automatic CDS stop codon adjustment (enables -P)")]
    adj_stop: bool,

    // ---- filters ----
    #[arg(short = 'i', value_name = "N", help = "Discard transcripts with intron > N")]
    maxintron: Option<i64>,
    #[arg(short = 'l', value_name = "N", help = "Discard transcripts shorter than N bases")]
    minlen: Option<u32>,
    #[arg(short = 'r', value_name = "RANGE", help = "Only transcripts overlapping RANGE (chr:start-end)")]
    range: Option<String>,
    #[arg(short = 'R', help = "For -r, require full containment within range")]
    within: bool,
    #[arg(long = "jmatch", value_name = "RANGE", help = "Only transcripts matching given junction")]
    jmatch: Option<String>,
    #[arg(short = 'U', help = "Discard single-exon transcripts")]
    multi_exon: bool,
    #[arg(short = 'C', help = "Discard mRNAs without CDS")]
    coding_only: bool,
    #[arg(long = "nc", help = "Discard mRNAs with CDS")]
    noncoding_only: bool,
    #[arg(short = 'V', help = "Discard mRNAs with in-frame stop codons (requires -g)")]
    valid_cds: bool,
    #[arg(short = 'H', help = "For -V, try alternate CDS phases")]
    alt_phases: bool,
    #[arg(short = 'B', help = "For -V, also check single-exon on opposite strand")]
    both_strands: bool,
    #[arg(short = 'N', help = "Discard mRNAs with non-canonical splice sites")]
    splice_check: bool,
    #[arg(short = 'J', help = "Discard mRNAs without complete CDS (START+STOP)")]
    full_cds: bool,
    #[arg(long = "no-pseudo", help = "Filter out pseudo records")]
    no_pseudo: bool,
    #[arg(long = "ids", value_name = "FILE", help = "Only keep IDs listed in FILE")]
    ids: Option<String>,
    #[arg(long = "nids", value_name = "FILE", help = "Discard IDs listed in FILE")]
    nids: Option<String>,

    // ---- attributes / formatting ----
    #[arg(short = 'F', help = "Keep all GFF attributes")]
    full_attrs: bool,
    #[arg(long = "keep-exon-attrs", help = "For -F, keep exon/CDS attributes")]
    keep_exon_attrs: bool,
    #[arg(short = 'G', help = "Move exon attributes to the transcript feature")]
    move_exon_attrs: bool,
    #[arg(long = "attrs", value_name = "LIST", help = "Only output listed attributes")]
    attrs: Option<String>,
    #[arg(long = "keep-genes", help = "Also preserve gene records in transcript-only mode")]
    keep_genes: bool,
    #[arg(long = "keep-comments", help = "Preserve GFF3 comment lines")]
    keep_comments: bool,
    #[arg(short = 'O', help = "Process non-transcript GFF records too")]
    process_others: bool,
    #[arg(short = 'D', help = "Decode url-encoded characters in attributes")]
    decode: bool,
    #[arg(short = 't', value_name = "NAME", help = "Use NAME as the track (column 2)")]
    track: Option<String>,
    #[arg(short = 'L', help = "Ensembl GTF to GFF3 conversion")]
    ensembl: bool,
    #[arg(short = 'A', help = "Add description from seq-info as a 'descr' attribute")]
    add_descr: bool,
    #[arg(short = 'Z', help = "Merge very close exons (intron < 4)")]
    merge_close: bool,

    // ---- input format ----
    #[arg(long = "in-bed", help = "Input is BED format")]
    in_bed: bool,
    #[arg(long = "in-tlf", help = "Input is TLF format")]
    in_tlf: bool,

    // ---- clustering ----
    #[arg(short = 'M', long = "merge", help = "Cluster transcripts into loci, discard redundant")]
    merge: bool,
    #[arg(long = "cluster-only", help = "Like -M but do not discard redundant transcripts")]
    cluster_only: bool,
    #[arg(short = 'K', help = "For -M, also discard contained transcripts")]
    discard_contained: bool,
    #[arg(long = "cset", help = "For -K, collapse single-exon transcripts into multi-exon exons")]
    cset: bool,
    #[arg(short = 'Q', help = "For -M, do not require boundary containment")]
    no_contain: bool,
    #[arg(short = 'Y', help = "For -M, discard overlapping single-exon transcripts (any strand)")]
    discard_set: bool,

    // ---- misc ----
    #[arg(long = "sort-alpha", help = "Sort chromosomes alphabetically")]
    sort_alpha: bool,
    #[arg(long = "sort-by", value_name = "FILE", help = "Sort chromosomes by order in FILE")]
    sort_by: Option<String>,
    #[arg(short = 'm', value_name = "FILE", help = "Reference name mapping table")]
    ref_map: Option<String>,
    #[arg(short = 's', value_name = "FILE", help = "Sequence info file (name len description)")]
    seq_info: Option<String>,
    #[arg(long = "stream", help = "Stream processing (no sorting/clustering)")]
    stream: bool,
    #[arg(long = "cov-info", help = "Report genome coverage by transcripts")]
    cov_info: bool,
    #[arg(long = "force-exons", help = "Treat lowest features as exon")]
    force_exons: bool,
    #[arg(long = "gene2exon", help = "Add exon feature spanning single-line genes")]
    gene2exon: bool,
    #[arg(long = "t-adopt", help = "Adopt orphan transcripts into overlapping genes")]
    t_adopt: bool,
    #[arg(long = "ignore-locus", help = "Discard locus features and attributes")]
    ignore_locus: bool,
    #[arg(short = 'v', help = "Verbose")]
    verbose: bool,
    #[arg(short = 'E', help = "Verbose (same as -v)")]
    expose: bool,
    /// Print version and exit
    #[arg(long = "version", action = clap::ArgAction::Version)]
    version: (),
}

fn open_writer(path: &str) -> Result<Box<dyn Write>> {
    if path == "-" || path == "stdout" {
        Ok(Box::new(BufWriter::new(io::stdout())))
    } else {
        let f = File::create(path).map_err(|e| anyhow!("Error creating file {}: {}", path, e))?;
        Ok(Box::new(BufWriter::new(f)))
    }
}

fn parse_range(s: &str) -> Result<RangeFilter> {
    // formats: [strand]chr:start-end  or  chr:start-end  or  chr:start  or  start-end
    let mut rf = RangeFilter::new();
    let mut rest = s;
    // optional leading strand
    let first = rest.chars().next();
    if first == Some('+') || first == Some('-') {
        rf.strand = first.unwrap();
        rest = &rest[1..];
    }
    if let Some(colon) = rest.rfind(':') {
        // could be chr:start-end
        let (chr, coords) = rest.split_at(colon);
        rf.ref_name = Some(chr.to_string());
        rest = &coords[1..];
    }
    if let Some(dash) = rest.rfind('-') {
        let a = &rest[..dash];
        let b = &rest[dash + 1..];
        if !a.is_empty() {
            rf.start = a.parse().unwrap_or(0);
        }
        if !b.is_empty() {
            rf.end = b.parse().unwrap_or(u32::MAX);
        }
    } else if !rest.is_empty() {
        // single coordinate
        let v: u32 = rest.parse().unwrap_or(0);
        rf.start = v;
        rf.end = v;
    }
    if rf.end == 0 {
        rf.end = u32::MAX;
    }
    Ok(rf)
}

fn parse_table_format(s: &str) -> Vec<TableField> {
    let mut cols = Vec::new();
    for raw in s.split(|c: char| c.is_whitespace() || c == ',' || c == ';' || c == ':') {
        let w = raw.trim();
        if w.is_empty() {
            continue;
        }
        if let Some(at) = w.strip_prefix('@') {
            let lw = at.to_ascii_lowercase();
            let f = match lw.as_str() {
                "chr" => TableField::Chr,
                "track" => TableField::Track,
                "id" => TableField::Id,
                "geneid" => TableField::GeneId,
                "genename" => TableField::GeneName,
                "parent" => TableField::Parent,
                "feature" => TableField::Feature,
                "start" => TableField::Start,
                "end" => TableField::End,
                "strand" => TableField::Strand,
                "numexons" => TableField::NumExons,
                "exons" => TableField::Exons,
                "introns" => TableField::Introns,
                "cds" => TableField::Cds,
                "covlen" => TableField::CovLen,
                "cdslen" => TableField::CdsLen,
                "attrs" => TableField::AllAttrs,
                _ => continue,
            };
            cols.push(f);
        } else {
            match w {
                "ID" | "transcript_id" => cols.push(TableField::Id),
                "geneID" | "gene_id" => cols.push(TableField::GeneId),
                "Parent" => cols.push(TableField::Parent),
                _ => cols.push(TableField::Attr(w.to_string())),
            }
        }
    }
    cols
}

fn load_id_list(path: &str) -> Result<HashSet<String>> {
    let mut set = HashSet::new();
    let reader = parser::open_reader(path)?;
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        for tok in line.split_whitespace() {
            if !tok.is_empty() {
                set.insert(tok.to_string());
            }
        }
    }
    Ok(set)
}

fn load_ref_map(path: &str) -> Result<std::collections::HashMap<String, String>> {
    let mut map = std::collections::HashMap::new();
    let reader = parser::open_reader(path)?;
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        let mut iter = line.split_whitespace();
        let Some(orig) = iter.next() else { continue };
        let Some(new) = iter.next() else { continue };
        map.insert(orig.to_string(), new.to_string());
    }
    Ok(map)
}

fn load_seq_info(path: &str) -> Result<std::collections::HashMap<String, SeqInfoEntry>> {
    let mut map = std::collections::HashMap::new();
    let reader = parser::open_reader(path)?;
    for line in std::io::BufRead::lines(reader) {
        let line = line?;
        let mut parts = line.splitn(3, |c: char| c == '\t' || c == ' ');
        let Some(name) = parts.next() else { continue };
        let name = name.trim();
        if name.is_empty() {
            continue;
        }
        let len: u32 = parts.next().and_then(|s| s.trim().parse().ok()).unwrap_or(0);
        let descr = parts.next().map(|s| s.trim().to_string()).filter(|s| !s.is_empty());
        map.insert(name.to_string(), SeqInfoEntry { length: len, description: descr });
    }
    Ok(map)
}

fn build_config(args: &Args) -> Result<AppConfig> {
    let mut cfg = AppConfig::default();

    // output format
    cfg.out_format = if args.table.is_some() {
        OutFormat::Table
    } else if args.gtf {
        OutFormat::Gtf
    } else if args.bed {
        OutFormat::Bed
    } else if args.tlf {
        OutFormat::Tlf
    } else {
        OutFormat::Gff3
    };

    // loader options
    let mut lo = LoaderOpts::default();
    lo.transcripts_only = !args.process_others;
    lo.gene2exon = args.gene2exon;
    lo.full_attributes = args.full_attrs;
    lo.keep_all_exon_attrs = args.keep_exon_attrs;
    lo.gather_exon_attrs = !args.full_attrs; // will refine below
    lo.merge_close_exons = args.merge_close;
    lo.ignore_locus = args.ignore_locus;
    lo.no_pseudo = args.no_pseudo;
    lo.bed_input = args.in_bed;
    lo.tlf_input = args.in_tlf;
    lo.keep_genes = args.keep_genes;
    lo.tr_adoption = args.t_adopt;
    lo.keep_gff3_comments = args.keep_comments;
    lo.sort_refs_alpha = args.sort_alpha;
    lo.force_exons = args.force_exons;
    lo.stream_in = args.stream;
    lo.ensembl_proc = args.ensembl;
    lo.match_all_introns = !args.discard_contained;
    lo.nc_span = args.no_contain;
    lo.d_ovl_set = args.discard_set;
    lo.cset_merge = args.cset;

    // -G: move exon attrs to transcript
    if args.move_exon_attrs {
        lo.gather_exon_attrs = true;
        lo.full_attributes = true;
    }
    if args.no_pseudo && !lo.full_attributes {
        lo.gather_exon_attrs = true;
        lo.full_attributes = true;
    }
    if lo.ensembl_proc {
        lo.full_attributes = true;
        lo.gather_exon_attrs = false;
    }
    // clustering
    if args.merge {
        lo.do_cluster = true;
        lo.collapse_redundant = true;
    }
    if args.cluster_only {
        lo.do_cluster = true;
        lo.collapse_redundant = false;
    }
    if args.cov_info {
        lo.do_cluster = true;
    }
    if args.discard_set {
        lo.nc_span = true; // -Y enforces -Q
    }
    cfg.loader = lo;

    // validation: -K/-Q/-Y require -M or --cluster-only
    let needs_merge = args.discard_contained || args.no_contain || args.discard_set;
    if needs_merge && !(cfg.loader.do_cluster) {
        bail!("Error: options -K,-Q,-Y require -M/--merge option!");
    }
    if args.cset && cfg.loader.match_all_introns {
        bail!("Error: option --cset requires option -K");
    }
    if args.keep_exon_attrs && !args.full_attrs {
        bail!("Error: option --keep-exon-attrs requires option -F!");
    }
    if args.sort_by.is_some() && args.sort_alpha {
        bail!("Error: options --sort-by and --sort-alpha are mutually exclusive!");
    }
    if args.within && args.range.is_none() {
        bail!("Error: option -R requires -r!");
    }

    // filters
    cfg.maxintron = args.maxintron.unwrap_or(999_000_000);
    cfg.min_len = args.minlen.unwrap_or(0);
    cfg.multi_exon = args.multi_exon;
    cfg.wcds_only = args.coding_only;
    cfg.wnc_only = args.noncoding_only;
    cfg.valid_cds_only = args.valid_cds;
    cfg.full_cds_only = args.full_cds;
    cfg.alt_phases = args.alt_phases;
    cfg.both_strands = args.both_strands;
    cfg.splice_check = args.splice_check;
    if cfg.full_cds_only {
        cfg.valid_cds_only = true;
    }
    cfg.rflt_within = args.within;
    cfg.flt_range = args.range.as_deref().map(|s| parse_range(s)).transpose()?;
    if let Some(jm) = &args.jmatch {
        let mut rf = parse_range(jm)?;
        if rf.strand == '.' {
            rf.strand = 0 as char;
        }
        cfg.flt_junction = Some(rf);
    }

    // ID filter
    if let Some(p) = &args.ids {
        cfg.id_flt = IDFltType::Only;
        cfg.flt_ids = load_id_list(p)?;
    } else if let Some(p) = &args.nids {
        cfg.id_flt = IDFltType::Exclude;
        cfg.flt_ids = load_id_list(p)?;
    }
    if cfg.flt_ids.is_empty() {
        cfg.id_flt = IDFltType::None;
    }

    // attrs / table
    if let Some(a) = &args.attrs {
        for tok in a.split(|c: char| c == ',' || c == ';' || c == ':') {
            let t = tok.trim();
            if !t.is_empty() {
                cfg.attr_list.insert(t.to_string());
            }
        }
        cfg.loader.attrs_filter = cfg.attr_list.len() > 1;
        cfg.loader.full_attributes = true;
    }
    if let Some(t) = &args.table {
        cfg.table_cols = parse_table_format(t);
    }
    if cfg.out_format == OutFormat::Table {
        cfg.loader.full_attributes = true;
    }

    // FASTA
    cfg.fasta_path = args.g.clone();
    cfg.w_padding = args.w_add.unwrap_or(0);
    cfg.wfa_no_cds = args.w_nocds;
    cfg.write_exon_segs = args.w_proj;
    cfg.star_stop = args.star_stop;
    cfg.add_cds_attrs = args.add_cds_attrs;
    cfg.add_has_cds = args.add_hascds;
    cfg.adjust_stop = args.adj_stop;
    if cfg.adjust_stop {
        cfg.add_cds_attrs = true;
    }

    // output files
    cfg.out_file = args.out.clone();
    cfg.w_file = args.w.clone();
    cfg.u_file = args.u.clone();
    cfg.x_file = args.x.clone();
    cfg.y_file = args.y.clone();
    cfg.j_file = args.j.clone();
    cfg.dup_info_file = args.d.clone();

    cfg.decode_chars = args.decode;
    cfg.cov_info = args.cov_info;
    cfg.add_descr = args.add_descr;
    cfg.track_label = args.track.clone();
    cfg.verbose = args.verbose || args.expose;

    // ref map / seq info / sort-by
    if let Some(p) = &args.ref_map {
        cfg.ref_map = load_ref_map(p)?;
    }
    if let Some(p) = &args.seq_info {
        cfg.seq_info = load_seq_info(p)?;
    }
    cfg.sort_by = args.sort_by.clone();

    // -g required for these
    if cfg.fasta_path.is_none() && cfg.needs_fasta() {
        bail!("Error: -g option is required for options -w/x/y/u/V/N/M !");
    }

    Ok(cfg)
}

fn print_mode_for(cfg: &AppConfig) -> GffPrintMode {
    match cfg.out_format {
        OutFormat::Gtf => {
            if cfg.loader.force_exons {
                GffPrintMode::Both
            } else {
                GffPrintMode::Any
            }
        }
        OutFormat::Bed => GffPrintMode::Bed,
        OutFormat::Tlf => GffPrintMode::Tlf,
        _ => {
            if cfg.loader.force_exons {
                GffPrintMode::Both
            } else {
                GffPrintMode::Any
            }
        }
    }
}

fn run() -> Result<()> {
    let args = Args::parse();
    let mut cfg = build_config(&args)?;
    // Capture the original command line for the GFF3 header comment.
    cfg.cmd_line = std::env::args().collect::<Vec<String>>().join(" ");

    // open FASTA db if needed
    let mut fasta_db = if let Some(p) = &cfg.fasta_path {
        Some(FastaDb::new(p)?)
    } else {
        None
    };

    // open output sinks
    let mut outs = Outputs {
        out: None,
        w: None,
        u: None,
        x: None,
        y: None,
        j: None,
        dup: None,
    };
    // default: if no output file at all, write records to stdout
    let any_out = cfg.out_file.is_some()
        || cfg.w_file.is_some()
        || cfg.u_file.is_some()
        || cfg.x_file.is_some()
        || cfg.y_file.is_some()
        || cfg.j_file.is_some()
        || cfg.cov_info;
    if !any_out {
        cfg.out_file = Some("-".to_string());
    }
    if let Some(p) = &cfg.out_file {
        outs.out = Some(open_writer(p)?);
    }
    if let Some(p) = &cfg.w_file {
        outs.w = Some(open_writer(p)?);
    }
    if let Some(p) = &cfg.u_file {
        outs.u = Some(open_writer(p)?);
    }
    if let Some(p) = &cfg.x_file {
        outs.x = Some(open_writer(p)?);
    }
    if let Some(p) = &cfg.y_file {
        outs.y = Some(open_writer(p)?);
    }
    if let Some(p) = &cfg.j_file {
        outs.j = Some(open_writer(p)?);
    }
    if let Some(p) = &cfg.dup_info_file {
        outs.dup = Some(open_writer(p)?);
    }

    let inputs: Vec<String> = if args.inputs.is_empty() {
        vec!["-".to_string()]
    } else {
        args.inputs.clone()
    };

    // ---- streaming mode ----
    if cfg.loader.stream_in {
        return run_stream(&mut cfg, &inputs, &mut fasta_db, &mut outs);
    }

    // ---- bulk load ----
    let mut g_data: Vec<types::GenomicSeqData> = Vec::new();
    for infile in &inputs {
        // extension-based BED/TLF detection
        let lower = infile.to_ascii_lowercase();
        let base = lower.strip_suffix(".gz").unwrap_or(&lower);
        if base.ends_with(".bed") {
            cfg.loader.bed_input = true;
        }
        if base.ends_with(".tlf") {
            cfg.loader.tlf_input = true;
        }
        let mut file_data = load_file(&mut cfg, infile)?;
        g_data.append(&mut file_data);
    }

    // process pass: filters + FASTA writing + CDS validation
    let mut intron_lists: std::collections::HashMap<String, types::CIntronList> =
        std::collections::HashMap::new();
    for gd in g_data.iter_mut() {
        for t in gd.rnas.iter_mut() {
            if !check_filters_mut(&cfg, t) {
                t.no_print();
                continue;
            }
            if cfg.fasta_path.is_some() && t.is_mrna {
                if let Some(db) = fasta_db.as_mut() {
                    if !process_transcript(&cfg, t, db, &mut outs) {
                        t.no_print();
                        continue;
                    }
                }
            }
            // collect introns for -j
            if cfg.j_file.is_some() && t.is_mrna && t.exons.len() > 1 {
                let list = intron_lists.entry(t.gseq_name.clone()).or_default();
                list.gseq_id = t.gseq_id;
                list.add(&t.id, t.strand, &t.exons);
            }
        }
        // gfs filter: drop genes with no surviving children when TFilters active
        if cfg.has_tfilters() {
            let live_ids: HashSet<String> = gd
                .rnas
                .iter()
                .filter(|t| t.is_printable())
                .filter_map(|t| t.parent_id.clone())
                .collect();
            gd.gfs.retain(|g| g.is_gene && (cfg.loader.keep_genes || live_ids.contains(&g.id)));
        }
    }

    // flush -j intron lists
    if let Some(fj) = outs.j.as_mut() {
        for (gname, list) in &intron_lists {
            list.print_to(fj, gname)?;
        }
    }

    // clustering
    let mut all_replaced: std::collections::HashMap<String, String> =
        std::collections::HashMap::new();
    if cfg.loader.do_cluster {
        for gd in g_data.iter_mut() {
            let res = cluster_gdata(&cfg, gd);
            for (k, v) in res.replaced_by {
                all_replaced.insert(k, v);
            }
            restore_strands(&cfg, gd);
        }
    }

    // coverage info
    if cfg.cov_info {
        let mut f_bases = 0u64;
        let mut r_bases = 0u64;
        let mut u_bases = 0u64;
        for gd in &g_data {
            f_bases += gd.f_bases;
            r_bases += gd.r_bases;
            u_bases += gd.u_bases;
        }
        let stdout = io::stdout();
        let mut o = stdout.lock();
        writeln!(o, "Total bases covered by transcripts:")?;
        if f_bases > 0 {
            writeln!(o, "\t{} on + strand", f_bases)?;
        }
        if r_bases > 0 {
            writeln!(o, "\t{} on - strand", r_bases)?;
        }
        if u_bases > 0 {
            writeln!(o, "\t{} on . strand", u_bases)?;
        }
    }

    // ---- print records ----
    let mode = print_mode_for(&cfg);
    let loctrack = cfg.track_label.clone().unwrap_or_else(|| "gffcl".to_string());

    // Only print GFF/GTF/BED/TLF/table records when an -o sink (or stdout
    // default) exists. If only FASTA outputs (-w/-x/-y/-u) were requested,
    // skip record printing entirely (matches original gffread behavior).
    if let Some(mut out) = outs.out.take() {
        let mut first_gff3 = cfg.out_format == OutFormat::Gff3;

        // build a gene id -> gf lookup per gseq for parent printing
        for gd in &g_data {
            let mut first_gseq_header = cfg.out_format == OutFormat::Gff3;
            let mut printed_genes: HashSet<String> = HashSet::new();
            // O(1) id lookups during printing (avoids O(n²) scans on large inputs)
            let rnas_by_id: std::collections::HashMap<&str, &GffObj> =
                gd.rnas.iter().map(|t| (t.id.as_str(), t)).collect();
            let gfs_by_id: std::collections::HashMap<&str, &GffObj> =
                gd.gfs.iter().map(|g| (g.id.as_str(), g)).collect();

            if cfg.out_format == OutFormat::Gff3 {
                output::print_seq_region(&mut out, &cfg, &gd.gseq_name, gd.seqreg_start, gd.seqreg_end)?;
            }

            if cfg.loader.do_cluster {
                for loc in &gd.loci {
                    let locname = format!("RLOC_{:08}", loc.locus_num);
                    // print locus line in GFF3
                    if cfg.out_format == OutFormat::Gff3 {
                        if first_gff3 {
                            output::print_gff3_header(&mut out, &cfg)?;
                            first_gff3 = false;
                        }
                        if first_gseq_header {
                            first_gseq_header = false;
                        }
                        output::print_locus(
                            &mut out,
                            &gd.gseq_name,
                            &loctrack,
                            &locname,
                            loc.start,
                            loc.end,
                            loc.strand,
                            &loc.gene_names,
                            &loc.gene_ids,
                            &loc.rna_ids,
                        )?;
                    }
                    // print rnas in the locus
                    for rid in &loc.rna_ids {
                        if all_replaced.contains_key(rid) {
                            // write dup info
                            if let Some(fd) = outs.dup.as_mut() {
                                let repl = &all_replaced[rid];
                                writeln!(fd, "{} => {}", rid, repl)?;
                            }
                            if cfg.verbose {
                                eprintln!("Info: {} discarded: superseded by {}", rid, all_replaced[rid]);
                            }
                            continue;
                        }
                        if let Some(&t) = rnas_by_id.get(rid.as_str()) {
                            if !t.is_printable() {
                                continue;
                            }
                            if cfg.out_format == OutFormat::Gff3 {
                                // print parent gene first
                                if let Some(pid) = &t.parent_id {
                                    if !printed_genes.contains(pid) {
                                        if let Some(&g) = gfs_by_id.get(pid.as_str()) {
                                            output::print_gxf(&mut out, g, &cfg, mode)?;
                                            printed_genes.insert(pid.clone());
                                        }
                                    }
                                }
                            }
                            if cfg.out_format == OutFormat::Table {
                                output::print_table_data(&mut out, t, &cfg, false)?;
                            } else {
                                output::print_gxf(&mut out, t, &cfg, mode)?;
                            }
                        }
                    }
                }
            } else {
                // non-clustered: print rnas interleaved with gfs
                let mut gfs_i = 0usize;
                for t in &gd.rnas {
                    if !t.is_printable() || all_replaced.contains_key(&t.id) {
                        continue;
                    }
                    // print genes with start <= t.start first
                    while gfs_i < gd.gfs.len() && gd.gfs[gfs_i].start <= t.start {
                        let g = &gd.gfs[gfs_i];
                        // Native gclib discards non-transcript records during
                        // finalize() unless --keep-genes (genes) or -O (others).
                        let printable = !cfg.loader.transcripts_only
                            || (g.is_gene && cfg.loader.keep_genes);
                        if printable && !printed_genes.contains(&g.id) && g.is_printable() {
                            if cfg.out_format == OutFormat::Gff3 {
                                if first_gff3 {
                                    output::print_gff3_header(&mut out, &cfg)?;
                                    first_gff3 = false;
                                }
                                output::print_gxf(&mut out, g, &cfg, mode)?;
                                printed_genes.insert(g.id.clone());
                            } else if cfg.out_format == OutFormat::Table {
                                output::print_table_data(&mut out, g, &cfg, false)?;
                                printed_genes.insert(g.id.clone());
                            }
                        }
                        gfs_i += 1;
                    }
                    if cfg.out_format == OutFormat::Gff3 {
                        if first_gff3 {
                            output::print_gff3_header(&mut out, &cfg)?;
                            first_gff3 = false;
                        }
                        if let Some(pid) = &t.parent_id {
                            if !printed_genes.contains(pid) {
                                if let Some(&g) = gfs_by_id.get(pid.as_str()) {
                                    let printable = !cfg.loader.transcripts_only
                                        || (g.is_gene && cfg.loader.keep_genes);
                                    if printable {
                                        output::print_gxf(&mut out, g, &cfg, mode)?;
                                        printed_genes.insert(pid.clone());
                                    }
                                }
                            }
                        }
                    }
                    if cfg.out_format == OutFormat::Table {
                        output::print_table_data(&mut out, t, &cfg, false)?;
                    } else {
                        output::print_gxf(&mut out, t, &cfg, mode)?;
                    }
                }
                // remaining gfs
                while gfs_i < gd.gfs.len() {
                    let g = &gd.gfs[gfs_i];
                    let printable = !cfg.loader.transcripts_only
                        || (g.is_gene && cfg.loader.keep_genes);
                    if printable && !printed_genes.contains(&g.id) && g.is_printable() {
                        if cfg.out_format == OutFormat::Gff3 {
                            if first_gff3 {
                                output::print_gff3_header(&mut out, &cfg)?;
                                first_gff3 = false;
                            }
                            output::print_gxf(&mut out, g, &cfg, mode)?;
                        } else if cfg.out_format == OutFormat::Table {
                            output::print_table_data(&mut out, g, &cfg, false)?;
                        }
                    }
                    gfs_i += 1;
                }
            }
        }

        out.flush()?;
    }

    // flush outputs
    for s in [&mut outs.w, &mut outs.u, &mut outs.x, &mut outs.y, &mut outs.j, &mut outs.dup] {
        if let Some(w) = s {
            w.flush()?;
        }
    }
    Ok(())
}

/// Streaming mode: process each transcript as it's parsed and print immediately.
fn run_stream(
    cfg: &mut AppConfig,
    inputs: &[String],
    fasta_db: &mut Option<FastaDb>,
    outs: &mut Outputs,
) -> Result<()> {
    let mode = print_mode_for(cfg);
    // Only open a record sink if -o was given or no FASTA outputs at all.
    let mut out = outs.out.take();
    let mut first_gff3 = cfg.out_format == OutFormat::Gff3;

    for infile in inputs {
        let lower = infile.to_ascii_lowercase();
        let base = lower.strip_suffix(".gz").unwrap_or(&lower);
        if base.ends_with(".bed") {
            cfg.loader.bed_input = true;
        }
        if base.ends_with(".tlf") {
            cfg.loader.tlf_input = true;
        }
        let file_data = load_file(cfg, infile)?;
        for gd in file_data {
            for mut t in gd.rnas {
                if !check_filters_mut(cfg, &mut t) {
                    continue;
                }
                if cfg.fasta_path.is_some() && t.is_mrna {
                    if let Some(db) = fasta_db.as_mut() {
                        if !process_transcript(cfg, &mut t, db, outs) {
                            continue;
                        }
                    }
                }
                if let Some(o) = out.as_mut() {
                    if first_gff3 && cfg.out_format == OutFormat::Gff3 {
                        output::print_gff3_header(o, cfg)?;
                        first_gff3 = false;
                    }
                    if cfg.out_format == OutFormat::Table {
                        output::print_table_data(o, &t, cfg, false)?;
                    } else {
                        output::print_gxf(o, &t, cfg, mode)?;
                    }
                }
            }
        }
    }
    if let Some(o) = out.as_mut() {
        o.flush()?;
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{}", e);
        std::process::exit(1);
    }
}

#[allow(dead_code)]
fn unused_path_check(_p: &Path) {}
