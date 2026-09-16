//! Transcript clustering into loci (`-M` / `--merge` / `--cluster-only`).
//!
//! Implements the core of `GffLoader::placeGf()`, `redundantTranscripts()`
//! and `collectLocusData()` from `gff_utils.cpp`.

use std::collections::{HashMap, HashSet};

use crate::config::AppConfig;
use crate::types::{merge_segments, GenomicSeqData, GffLocus, GffObj, GSeg, GFF_MAX_LOCUS};

/// Result of clustering one genomic sequence.
#[derive(Debug, Default)]
pub struct ClusterResult {
    /// replaced_transcript_id -> surviving_transcript_id
    pub replaced_by: HashMap<String, String>,
}

/// Check whether any of `t`'s exons overlaps any segment in `mexons`.
fn exon_overlaps_mexons(t: &GffObj, mexons: &[GSeg]) -> bool {
    if t.exons.is_empty() {
        // gene-like: span overlap
        return mexons.iter().any(|m| m.overlap(&GSeg::new(t.start, t.end)));
    }
    for ex in &t.exons {
        for m in mexons {
            if ex.overlap(m) {
                return true;
            }
        }
    }
    false
}

/// Redundancy test between two transcripts.
///
/// Returns `Some(true)` if `a` is the "bigger" (surviving) transcript and `b`
/// should be replaced by `a`, `Some(false)` if the reverse, or `None` if they
/// are not redundant.
fn redundant(a: &GffObj, b: &GffObj, cfg: &AppConfig) -> Option<bool> {
    let adj = if cfg.loader.d_ovl_set { 1u32 } else { 0u32 };
    // span overlap check
    if a.start > b.end + adj || b.start > a.end + adj {
        return None;
    }
    if b.strand != '.' && a.strand != '.' && a.strand != b.strand {
        return None;
    }

    let ai = a.exons.len();
    let bj = b.exons.len();
    let a_bigger = a.covlen >= b.covlen;

    if cfg.loader.match_all_introns {
        // must have the same number of exons
        if ai != bj {
            return None;
        }
        if !cfg.loader.nc_span {
            // require containment of the smaller span
            if a_bigger {
                if a.start > b.start || a.end < b.end {
                    return None;
                }
            } else if b.start > a.start || b.end < a.end {
                return None;
            }
        }
        // all introns must match
        for i in 0..ai.saturating_sub(1) {
            if a.exons[i].end != b.exons[i].end || a.exons[i + 1].start != b.exons[i + 1].start {
                return None;
            }
        }
        // single-exon: require 80% overlap if nc_span, else containment
        if ai == 1 {
            if cfg.loader.nc_span {
                let minlen = a.covlen.min(b.covlen);
                let ovl = a.exons[0].overlap_len(&b.exons[0]);
                if (ovl as f64) < 0.8 * (minlen as f64) {
                    return None;
                }
            } else {
                let (big, small) = if a_bigger { (a, b) } else { (b, a) };
                if small.start < big.start || small.end > big.end {
                    return None;
                }
            }
        }
        return Some(a_bigger);
    }

    // match_all_introns == false: intron-chain containment
    let (bigger, smaller) = if a_bigger { (a, b) } else { (b, a) };
    if ai == 1 && bj == 1 {
        if cfg.loader.nc_span {
            let minlen = a.covlen.min(b.covlen);
            let ovl = a.exons[0].overlap_len(&b.exons[0]);
            return if (ovl as f64) >= 0.8 * (minlen as f64) {
                Some(a_bigger)
            } else {
                None
            };
        }
        return if smaller.start >= bigger.start && smaller.end <= bigger.end {
            Some(a_bigger)
        } else {
            None
        };
    }

    if cfg.loader.cset_merge && smaller.exons.len() == 1 {
        return if unspl_contained(smaller, bigger, cfg) {
            Some(a_bigger)
        } else {
            None
        };
    }

    // both multi-exon: check intron-chain containment
    if ai == 0 || bj == 0 {
        return None;
    }
    if a.exons[ai - 1].start < b.exons[0].end || b.exons[bj - 1].start < a.exons[0].end {
        return None; // intron chains do not overlap
    }
    if !cfg.loader.nc_span && (bigger.start > smaller.start || bigger.end < smaller.end) {
        return None;
    }
    // find first intron overlap
    let mut i = 1usize;
    let mut j = 1usize;
    loop {
        if i >= ai || j >= bj {
            break;
        }
        let eis = a.exons[i - 1].end;
        let eie = a.exons[i].start;
        let ejs = b.exons[j - 1].end;
        let eje = b.exons[j].start;
        if eje < eis {
            j += 1;
            continue;
        }
        if eie < ejs {
            i += 1;
            continue;
        }
        break;
    }
    if i > 1 && j > 1 {
        return None;
    }
    if i >= ai || j >= bj {
        return None;
    }
    if a.exons[i - 1].end != b.exons[j - 1].end || a.exons[i].start != b.exons[j].start {
        return None;
    }
    // remaining introns must match in sequence
    let mut i = i + 1;
    let mut j = j + 1;
    while i < ai && j < bj {
        if a.exons[i - 1].end != b.exons[j - 1].end || a.exons[i].start != b.exons[j].start {
            return None;
        }
        i += 1;
        j += 1;
    }
    Some(a_bigger)
}

/// Check if a single-exon transcript is contained in an exon of a multi-exon
/// one without crossing intron-exon boundaries (`unsplContained`).
fn unspl_contained(ti: &GffObj, tj: &GffObj, cfg: &AppConfig) -> bool {
    let max_intron_ovl = if cfg.loader.d_ovl_set { 25u32 } else { 0u32 };
    if ti.exons.len() != 1 {
        return false;
    }
    for (j, exj) in tj.exons.iter().enumerate() {
        let exon_overlap = if cfg.loader.nc_span {
            if cfg.loader.d_ovl_set {
                ti.exons[0].overlap_len_coords(exj.start.saturating_sub(1), exj.end + 1) > 0
            } else {
                ti.exons[0].overlap_len(exj) as f64 >= 0.8 * ti.covlen as f64
            }
        } else {
            ti.exons[0].end <= exj.end && ti.exons[0].start >= exj.start
        };
        if exon_overlap {
            if (j > 0 && ti.start + max_intron_ovl < tj.exons[j].start)
                || (j < tj.exons.len() - 1 && ti.end > tj.exons[j].end + max_intron_ovl)
            {
                return false;
            }
            return true;
        }
    }
    false
}

/// Cluster all transcripts of one genomic sequence into loci.
pub fn cluster_gdata(cfg: &AppConfig, gd: &mut GenomicSeqData) -> ClusterResult {
    let mut result = ClusterResult::default();
    if gd.rnas.is_empty() && gd.gfs.is_empty() {
        return result;
    }

    let mut loci: Vec<GffLocus> = Vec::new();
    // Only cluster transcripts that survived filtering (printable).
    let rnas_snapshot: Vec<GffObj> = gd.rnas.iter().filter(|t| t.is_printable()).cloned().collect();

    // map id -> index into rnas_snapshot for redundancy lookups
    let id_index: HashMap<String, usize> = rnas_snapshot
        .iter()
        .enumerate()
        .map(|(i, r)| (r.id.clone(), i))
        .collect();

    for t in &rnas_snapshot {
        // dOvlSET: temporarily neutralize strand for single-exon
        let t_strand = if cfg.loader.d_ovl_set && t.exons.len() == 1 {
            '.'
        } else {
            t.strand
        };

        // find overlapping loci (immutable scan)
        let mut found: Vec<usize> = Vec::new();
        for (li, loc) in loci.iter().enumerate() {
            if (loc.strand == '+' || loc.strand == '-')
                && t_strand != '.'
                && loc.strand != t_strand
            {
                continue;
            }
            let t_start = if cfg.loader.d_ovl_set { t.start.saturating_sub(1) } else { t.start };
            let t_end = if cfg.loader.d_ovl_set { t.end + 1 } else { t.end };
            if t_start > loc.end {
                // loci sorted by start; allow far gaps up to GFF_MAX_LOCUS
                if t.start.saturating_sub(loc.start) > GFF_MAX_LOCUS {
                    // can't break globally since not strictly sorted by end; just skip
                }
                continue;
            }
            if t_end < loc.start {
                continue;
            }
            if exon_overlaps_mexons(t, &loc.mexons) {
                found.push(li);
            }
        }

        if found.is_empty() {
            // create a new locus
            let mut loc = GffLocus::new();
            loc.gseq_id = t.gseq_id;
            loc.strand = t_strand;
            loc.start = t.start;
            loc.end = t.end;
            loc.mexons = if t.exons.is_empty() {
                vec![GSeg::new(t.start, t.end)]
            } else {
                t.exons.clone()
            };
            merge_segments(&mut loc.mexons);
            loc.rna_ids.push(t.id.clone());
            loc.rna_covlens.push(t.covlen);
            loc.is_mrna = t.is_mrna;
            loc.t_maxcov = Some(0);
            loci.push(loc);
            continue;
        }

        // collapse-redundant check against all rnas in found loci
        if cfg.loader.collapse_redundant && t.is_mrna {
            let mut candidates: Vec<String> = Vec::new();
            for &li in &found {
                for rid in &loci[li].rna_ids {
                    if rid != &t.id && !result.replaced_by.contains_key(rid) {
                        candidates.push(rid.clone());
                    }
                }
            }
            for rid in candidates {
                let Some(&oi) = id_index.get(&rid) else {
                    continue;
                };
                let other = &rnas_snapshot[oi];
                if let Some(a_bigger) = redundant(t, other, cfg) {
                    if a_bigger {
                        result.replaced_by.insert(other.id.clone(), t.id.clone());
                    } else {
                        result.replaced_by.insert(t.id.clone(), other.id.clone());
                    }
                }
            }
        }

        // add t to the first found locus, merge the rest into it
        let primary = found[0];
        {
            let loc = &mut loci[primary];
            loc.rna_ids.push(t.id.clone());
            loc.rna_covlens.push(t.covlen);
            if t.covlen > loc.rna_covlens.iter().copied().max().unwrap_or(0) {
                loc.t_maxcov = Some(loc.rna_ids.len() - 1);
            }
            if t.is_mrna {
                loc.is_mrna = true;
            }
            // merge exons
            for ex in &t.exons {
                loc.mexons.push(*ex);
            }
            if t.exons.is_empty() {
                loc.mexons.push(GSeg::new(t.start, t.end));
            }
            merge_segments(&mut loc.mexons);
            if t.start < loc.start {
                loc.start = t.start;
            }
            if t.end > loc.end {
                loc.end = t.end;
            }
            if (loc.strand == '.' || loc.strand == 0 as char) && t_strand != '.' {
                loc.strand = t_strand;
            }
        }
        // merge additional loci into primary (highest indices first)
        for &li in found.iter().skip(1).rev() {
            let other = loci[li].clone();
            let loc = &mut loci[primary];
            for rid in &other.rna_ids {
                loc.rna_ids.push(rid.clone());
            }
            loc.rna_covlens.extend(other.rna_covlens.iter().copied());
            for ex in &other.mexons {
                loc.mexons.push(*ex);
            }
            merge_segments(&mut loc.mexons);
            if other.start < loc.start {
                loc.start = other.start;
            }
            if other.end > loc.end {
                loc.end = other.end;
            }
            if other.is_mrna {
                loc.is_mrna = true;
            }
            loci.remove(li);
        }
        // keep loci sorted by start (bubble the primary down if needed)
        let mut k = primary;
        while k > 0 && loci[k].start < loci[k - 1].start {
            loci.swap(k, k - 1);
            k -= 1;
        }
    }

    // also place gene features (gfs) into overlapping loci (for gene_names/ids)
    for g in &gd.gfs {
        for loc in loci.iter_mut() {
            if g.start <= loc.end && g.end >= loc.start {
                if let Some(n) = g.get_attr("Name") {
                    if !loc.gene_names.contains(&n.to_string()) {
                        loc.gene_names.push(n.to_string());
                    }
                }
                if !loc.gene_ids.contains(&g.id) {
                    loc.gene_ids.push(g.id.clone());
                }
            }
        }
    }

    // collect gene names/ids from transcripts
    for loc in loci.iter_mut() {
        let mut gnames: HashSet<String> = loc.gene_names.iter().cloned().collect();
        let mut gids: HashSet<String> = loc.gene_ids.iter().cloned().collect();
        for rid in &loc.rna_ids {
            if let Some(&idx) = id_index.get(rid) {
                let t = &rnas_snapshot[idx];
                if let Some(n) = t.get_gene_name() {
                    gnames.insert(n.to_string());
                }
                if let Some(gid) = t.get_gene_id() {
                    gids.insert(gid.to_string());
                }
            }
        }
        loc.gene_names = gnames.into_iter().collect();
        loc.gene_ids = gids.into_iter().collect();
    }

    // assign locus numbers and decide strand by majority
    let mut n = 0u32;
    for loc in loci.iter_mut() {
        n += 1;
        loc.locus_num = n;
        let mut fstrand = 0i32;
        let mut rstrand = 0i32;
        let mut ustrand = 0i32;
        for rid in &loc.rna_ids {
            if let Some(&idx) = id_index.get(rid) {
                let t = &rnas_snapshot[idx];
                let s = if cfg.loader.d_ovl_set && t.exons.len() == 1 {
                    t.orig_strand()
                } else {
                    t.strand
                };
                let s = if s == 0 as char { '.' } else { s };
                match s {
                    '+' => fstrand += 1,
                    '-' => rstrand += 1,
                    _ => ustrand += 1,
                }
            }
        }
        if (fstrand > 0 && rstrand > 0) || (fstrand == 0 && rstrand == 0) {
            loc.strand = '.';
        } else if rstrand > 0 {
            loc.strand = '-';
        } else {
            loc.strand = '+';
        }
    }

    gd.loci = loci;
    result
}

/// Restore original strands on transcripts after clustering (-Y).
pub fn restore_strands(cfg: &AppConfig, gd: &mut GenomicSeqData) {
    if !cfg.loader.d_ovl_set {
        return;
    }
    for t in gd.rnas.iter_mut() {
        if t.exons.len() == 1 {
            let o = t.orig_strand();
            if o != 0 as char {
                t.strand = o;
            }
        }
    }
}
