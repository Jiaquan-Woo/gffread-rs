//! Runtime configuration mirroring the C++ globals in `gff_utils.cpp` and
//! the `GffLoader` option bitfield.

use std::collections::{HashMap, HashSet};

use crate::types::{GffNames, IDFltType, RangeFilter, TableField};

/// Loader / processing options (subset of the C++ `GffLoader` bitfield that
/// affects parsing and clustering).
#[derive(Debug, Clone)]
pub struct LoaderOpts {
    pub transcripts_only: bool,
    pub gene2exon: bool,
    pub full_attributes: bool,
    pub keep_all_exon_attrs: bool,
    pub gather_exon_attrs: bool,
    pub merge_close_exons: bool,
    pub ignore_locus: bool,
    pub no_pseudo: bool,
    pub bed_input: bool,
    pub tlf_input: bool,
    pub keep_genes: bool,
    pub tr_adoption: bool,
    pub keep_gff3_comments: bool,
    pub sort_refs_alpha: bool,
    pub do_cluster: bool,
    pub collapse_redundant: bool,
    pub match_all_introns: bool,
    pub nc_span: bool,
    pub d_ovl_set: bool,
    pub force_exons: bool,
    pub stream_in: bool,
    pub ensembl_proc: bool,
    pub attrs_filter: bool,
    pub cset_merge: bool,
}

impl Default for LoaderOpts {
    fn default() -> Self {
        LoaderOpts {
            transcripts_only: true,
            gene2exon: false,
            full_attributes: false,
            keep_all_exon_attrs: false,
            gather_exon_attrs: true,
            merge_close_exons: false,
            ignore_locus: false,
            no_pseudo: false,
            bed_input: false,
            tlf_input: false,
            keep_genes: false,
            tr_adoption: false,
            keep_gff3_comments: false,
            sort_refs_alpha: false,
            do_cluster: false,
            collapse_redundant: false,
            match_all_introns: true,
            nc_span: false,
            d_ovl_set: false,
            force_exons: false,
            stream_in: false,
            ensembl_proc: false,
            attrs_filter: false,
            cset_merge: false,
        }
    }
}

/// Output-format selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OutFormat {
    Gff3,
    Gtf,
    Bed,
    Tlf,
    Table,
}

/// All runtime configuration assembled from CLI args.
#[derive(Debug, Clone)]
pub struct AppConfig {
    pub loader: LoaderOpts,
    pub out_format: OutFormat,

    // filters / validation
    pub maxintron: i64,
    pub min_len: u32,
    pub multi_exon: bool, // -U
    pub wcds_only: bool,  // -C
    pub wnc_only: bool,   // --nc
    pub valid_cds_only: bool, // -V
    pub full_cds_only: bool,  // -J
    pub alt_phases: bool,     // -H
    pub both_strands: bool,   // -B
    pub splice_check: bool,   // -N

    pub flt_range: Option<RangeFilter>,
    pub rflt_within: bool,
    pub flt_junction: Option<RangeFilter>,
    pub id_flt: IDFltType,
    pub flt_ids: HashSet<String>,
    pub attr_list: HashSet<String>,
    pub table_cols: Vec<TableField>,

    // sequence / fasta
    pub fasta_path: Option<String>,
    pub w_padding: i64,
    pub wfa_no_cds: bool,
    pub write_exon_segs: bool, // -W
    pub star_stop: bool,       // -S
    pub add_cds_attrs: bool,   // -P
    pub add_has_cds: bool,
    pub adjust_stop: bool,

    // output files
    pub out_file: Option<String>,
    pub w_file: Option<String>,
    pub u_file: Option<String>,
    pub x_file: Option<String>,
    pub y_file: Option<String>,
    pub j_file: Option<String>,
    pub dup_info_file: Option<String>,

    pub decode_chars: bool,
    pub cov_info: bool,
    pub add_descr: bool,
    pub track_label: Option<String>,
    pub verbose: bool,

    // ref name mapping (-m) and seq info (-s)
    pub ref_map: HashMap<String, String>,
    pub seq_info: HashMap<String, SeqInfoEntry>,
    pub sort_by: Option<String>,

    pub names: GffNames,
    pub header_lines: Vec<String>,
    pub cmd_line: String,
}

#[derive(Debug, Clone)]
pub struct SeqInfoEntry {
    pub length: u32,
    pub description: Option<String>,
}

impl Default for AppConfig {
    fn default() -> Self {
        AppConfig {
            loader: LoaderOpts::default(),
            out_format: OutFormat::Gff3,
            maxintron: 999_000_000,
            min_len: 0,
            multi_exon: false,
            wcds_only: false,
            wnc_only: false,
            valid_cds_only: false,
            full_cds_only: false,
            alt_phases: false,
            both_strands: false,
            splice_check: false,
            flt_range: None,
            rflt_within: false,
            flt_junction: None,
            id_flt: IDFltType::None,
            flt_ids: HashSet::new(),
            attr_list: HashSet::new(),
            table_cols: Vec::new(),
            fasta_path: None,
            w_padding: 0,
            wfa_no_cds: false,
            write_exon_segs: false,
            star_stop: false,
            add_cds_attrs: false,
            add_has_cds: false,
            adjust_stop: false,
            out_file: None,
            w_file: None,
            u_file: None,
            x_file: None,
            y_file: None,
            j_file: None,
            dup_info_file: None,
            decode_chars: false,
            cov_info: false,
            add_descr: false,
            track_label: None,
            verbose: false,
            ref_map: HashMap::new(),
            seq_info: HashMap::new(),
            sort_by: None,
            names: GffNames::default(),
            header_lines: Vec::new(),
            cmd_line: String::new(),
        }
    }
}

impl AppConfig {
    /// True when any transcript-level filter is active.
    pub fn has_tfilters(&self) -> bool {
        self.multi_exon || self.wcds_only || self.wnc_only
    }

    pub fn needs_fasta(&self) -> bool {
        self.valid_cds_only
            || self.splice_check
            || self.w_file.is_some()
            || self.x_file.is_some()
            || self.y_file.is_some()
            || self.u_file.is_some()
            || self.full_cds_only
    }
}
