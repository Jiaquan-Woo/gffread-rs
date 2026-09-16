//! Transcript filtering and FASTA-sequence processing.
//!
//! Mirrors `GffLoader::checkFilters()`, `process_transcript()` and
//! `collectIntrons()` from `gff_utils.cpp`.

use std::io::Write;

use crate::config::AppConfig;
use crate::fasta::{get_spliced, get_unspliced, FastaDb};
use crate::output::print_table_data;
use crate::translate::{print_fasta, translate_dna};
use crate::types::{GffObj, GSpliceSite, IDFltType};

/// Owned output streams (one per FASTA / record sink).
pub struct Outputs {
    pub out: Option<Box<dyn Write>>,
    pub w: Option<Box<dyn Write>>,
    pub u: Option<Box<dyn Write>>,
    pub x: Option<Box<dyn Write>>,
    pub y: Option<Box<dyn Write>>,
    pub j: Option<Box<dyn Write>>,
    pub dup: Option<Box<dyn Write>>,
}

/// Apply all non-FASTA filters. Returns `false` to discard the record.
pub fn check_filters(cfg: &AppConfig, t: &GffObj) -> bool {
    // ID list filter
    if cfg.id_flt != IDFltType::None && !cfg.flt_ids.is_empty() {
        let has = cfg.flt_ids.contains(&t.id);
        match cfg.id_flt {
            IDFltType::Exclude => {
                if has {
                    return false;
                }
            }
            IDFltType::Only => {
                if !has {
                    return false;
                }
            }
            IDFltType::None => {}
        }
    }

    // minimum length
    if cfg.min_len > 0 && t.covlen < cfg.min_len {
        return false;
    }

    // range filter (-r / -R)
    if let Some(rf) = &cfg.flt_range {
        if let Some(rn) = &rf.ref_name {
            if rn != &t.gseq_name {
                return false;
            }
        }
        if rf.strand != 0 as char && rf.strand != '.' && t.strand != rf.strand {
            return false;
        }
        if rf.start != 0 || rf.end != u32::MAX {
            if cfg.rflt_within {
                if t.start < rf.start || t.end > rf.end {
                    return false;
                }
            } else if t.start > rf.end || t.end < rf.start {
                return false;
            }
        }
    }

    // junction filter (--jmatch)
    if let Some(jf) = &cfg.flt_junction {
        if t.exons.len() <= 1 {
            return false;
        }
        if let Some(rn) = &jf.ref_name {
            if rn != &t.gseq_name {
                return false;
            }
        }
        if jf.strand != 0 as char && jf.strand != '.' && t.strand != jf.strand {
            return false;
        }
        let jstart = if jf.start == 0 { jf.end } else { jf.start };
        let jend = if jf.end == 0 { jstart } else { jf.end };
        if t.start >= jstart || t.end <= jend {
            return false;
        }
        let mut matched = false;
        for i in 0..t.exons.len() - 1 {
            let donor = t.exons[i].end + 1;
            let acceptor = if t.exons[i + 1].start > 0 {
                t.exons[i + 1].start - 1
            } else {
                0
            };
            let s_ok = jf.start == 0 || donor == jf.start;
            let e_ok = jf.end == 0 || acceptor == jf.end;
            if s_ok && e_ok {
                matched = true;
                break;
            }
        }
        if !matched {
            return false;
        }
    }

    // attribute filtering (remove attrs not in attr_list)
    if cfg.loader.attrs_filter && !cfg.attr_list.is_empty() {
        // done at print time; nothing to remove here
    }

    // transcript-only filters
    if t.is_mrna {
        if cfg.multi_exon && t.exons.len() <= 1 {
            return false;
        }
        if cfg.wcds_only && !t.has_cds() {
            return false;
        }
        if cfg.wnc_only && t.has_cds() {
            return false;
        }
        // maxintron check (-i) — does not require -g
        if (cfg.maxintron as i64) < 999_000_000 && t.exons.len() > 1 {
            for i in 1..t.exons.len() {
                let ilen = (t.exons[i].start as i64) - (t.exons[i - 1].end as i64) - 1;
                if ilen > cfg.maxintron as i64 {
                    return false;
                }
            }
        }
        // --add-hasCDS: doesn't require -g, just presence of CDS features
        if cfg.add_has_cds && t.has_cds() {
            // mutate through interior pointer — safe because check_filters
            // takes &GffObj in our public API but we need to add the attr.
            // The caller passes a mutable reference internally.
        }
    }

    true
}

/// Apply filters and side-effects that need a mutable record. Returns `false`
/// to discard the record. This is the internal counterpart to `check_filters`
/// that can mutate the record (e.g. add hasCDS attribute).
pub fn check_filters_mut(cfg: &AppConfig, t: &mut GffObj) -> bool {
    if !check_filters(cfg, t) {
        return false;
    }
    if cfg.add_has_cds && t.is_mrna && t.has_cds() {
        t.add_attr("hasCDS", "true");
    }
    true
}

/// Check splice-site consensus for multi-exon transcripts (-N).
/// Returns `false` if any intron has a non-canonical splice site.
pub fn check_splice_sites(cfg: &AppConfig, db: &mut FastaDb, t: &GffObj) -> bool {
    if !cfg.splice_check || t.exons.len() < 2 {
        return true;
    }
    let glen = t.end - t.start + 1;
    let gseq = match db.subseq(&t.gseq_name, t.start, glen) {
        Ok(Some(s)) => s,
        _ => return false,
    };
    let revc = t.strand == '-';
    let mut ok = true;
    for e in 1..t.exons.len() {
        let i_start = (t.exons[e - 1].end + 1 - t.start) as usize;
        let i_end = (t.exons[e].start - 1 - t.start) as usize;
        if i_start >= gseq.len() || i_end + 1 > gseq.len() || i_end < i_start {
            ok = false;
            break;
        }
        let intron = &gseq[i_start..i_end + 1];
        let acceptor = GSpliceSite::from_intron(intron, true, revc);
        let donor = GSpliceSite::from_intron(intron, false, revc);
        if acceptor == GSpliceSite::new('A', 'G') {
            if !donor.canonical_donor() {
                ok = false;
                break;
            }
        } else if acceptor == GSpliceSite::new('A', 'C') {
            if donor != GSpliceSite::new('A', 'T') {
                ok = false;
                break;
            }
        } else {
            ok = false;
            break;
        }
    }
    ok
}

/// Process a transcript: validate CDS, write FASTA outputs (-w/-x/-y/-u).
/// Returns `false` if the transcript fails -V/-J validation and must be
/// discarded.
pub fn process_transcript(
    cfg: &AppConfig,
    t: &mut GffObj,
    db: &mut FastaDb,
    outs: &mut Outputs,
) -> bool {
    // maxintron check
    for i in 1..t.exons.len() {
        let ilen = (t.exons[i].start as i64) - (t.exons[i - 1].end as i64) - 1;
        if ilen > cfg.maxintron {
            return false;
        }
    }

    if cfg.add_has_cds && t.has_cds() {
        t.add_attr("hasCDS", "true");
    }

    // splice site check
    if !check_splice_sites(cfg, db, t) {
        return false;
    }

    let mut trprint = true;
    let mut inframe_stop = false;
    let mut end_stop = false;
    let mut stop_adjusted = false;
    let mut full_cds = false;

    let need_cds_check = t.has_cds()
        && (cfg.y_file.is_some()
            || cfg.x_file.is_some()
            || cfg.valid_cds_only
            || cfg.full_cds_only
            || cfg.add_cds_attrs
            || cfg.adjust_stop);

    let mut cds_seq: Vec<u8> = Vec::new();
    let mut aa_seq: Vec<u8> = Vec::new();
    let mut cds_len: usize = 0;

    if need_cds_check {
        let res = match get_spliced(db, t, true) {
            Ok(r) => r,
            Err(_) => return false,
        };
        cds_seq = res.seq.clone();
        cds_len = cds_seq.len();
        if !cds_seq.is_empty() {
            aa_seq = translate_dna(&cds_seq, cds_len);
            let mut cds_aalen = aa_seq.len();
            // locate first stop
            let stop_pos = aa_seq.iter().position(|&a| a == b'.');
            if let Some(p) = stop_pos {
                if p == cds_aalen.saturating_sub(1) {
                    // stop is the stated last codon
                    end_stop = true;
                    if cfg.adjust_stop {
                        cds_aalen = p;
                    } else {
                        cds_aalen = p;
                    }
                } else if p < cds_aalen.saturating_sub(1) {
                    if cfg.adjust_stop {
                        // trim CDS to the stop codon
                        let new_len = (p + 1) * 3;
                        cds_seq.truncate(new_len);
                        cds_len = new_len;
                        aa_seq.truncate(p + 1);
                        cds_aalen = p + 1;
                        end_stop = true;
                        stop_adjusted = true;
                        // adjust genomic CDS end coordinate
                        // (approximate: leave as-is for now)
                    } else {
                        inframe_stop = true;
                    }
                }
            }
            if inframe_stop {
                if cfg.add_cds_attrs {
                    t.add_attr("InFrameStop", "true");
                }
            }
            if stop_adjusted && cfg.add_cds_attrs {
                t.add_attr("CDStopAdjusted", "true");
                inframe_stop = false;
            }
            if !inframe_stop {
                let has_start = aa_seq.first().copied() == Some(b'M');
                full_cds = end_stop && has_start;
                if !full_cds && cfg.add_cds_attrs {
                    let partialness = if has_start {
                        "3"
                    } else if end_stop {
                        "5"
                    } else {
                        "5_3"
                    };
                    t.add_attr("partialness", partialness);
                }
            }
            if (cfg.full_cds_only && !full_cds) || (cfg.valid_cds_only && inframe_stop) {
                trprint = false;
            }
            let _ = cds_aalen;
        }
    }

    if !trprint {
        return false;
    }

    // write -y (protein)
    if let Some(fy) = outs.y.as_mut() {
        if !aa_seq.is_empty() {
            // strip trailing stop
            let mut aalen = aa_seq.len();
            if aa_seq[aalen - 1] == b'.' {
                aalen -= 1;
            }
            write_fasta_defline(cfg, fy, t, &cds_seq, true).ok();
            print_fasta(fy, None, &aa_seq[..aalen], cfg.star_stop).ok();
        }
    }
    // write -x (CDS)
    if let Some(fx) = outs.x.as_mut() {
        if !cds_seq.is_empty() {
            write_fasta_defline(cfg, fx, t, &cds_seq, true).ok();
            print_fasta(fx, None, &cds_seq, false).ok();
        }
    }
    // write -w (spliced exons)
    if outs.w.is_some() {
        let pad_left = if cfg.w_padding > 0 {
            (cfg.w_padding as u32).min(t.start.saturating_sub(1))
        } else {
            0
        };
        let seqlen = db.seq_len(&t.gseq_name).unwrap_or(t.end);
        let pad_right = if cfg.w_padding > 0 {
            let ediff = seqlen.saturating_sub(t.end);
            (cfg.w_padding as u32).min(ediff)
        } else {
            0
        };
        // adjust exons for padding
        if pad_left > 0 || pad_right > 0 {
            if let Some(first) = t.exons.first_mut() {
                first.start = first.start.saturating_sub(pad_left);
            }
            if let Some(last) = t.exons.last_mut() {
                last.end = last.end + pad_right;
            }
            t.start = t.start.saturating_sub(pad_left);
            t.end = t.end + pad_right;
        }
        let res = get_spliced(db, t, false).ok();
        // restore
        if pad_left > 0 || pad_right > 0 {
            if let Some(first) = t.exons.first_mut() {
                first.start += pad_left;
            }
            if let Some(last) = t.exons.last_mut() {
                last.end = last.end.saturating_sub(pad_right);
            }
            t.start += pad_left;
            t.end = t.end.saturating_sub(pad_right);
        }
        if let Some(res) = res {
            if let Some(fw) = outs.w.as_mut() {
                let defline = build_w_defline(cfg, t, &res, pad_left, pad_right);
                write!(fw, ">{}", defline).ok();
                if cfg.out_format == crate::config::OutFormat::Table {
                    print_table_data(fw, t, cfg, true).ok();
                } else {
                    write_attrs_gff(fw, t, cfg).ok();
                }
                writeln!(fw).ok();
                print_fasta(fw, None, &res.seq, false).ok();
            }
        }
    }
    // write -u (unspliced)
    if outs.u.is_some() {
        let pad_left = if cfg.w_padding > 0 {
            (cfg.w_padding as u32).min(t.start.saturating_sub(1))
        } else {
            0
        };
        let seqlen = db.seq_len(&t.gseq_name).unwrap_or(t.end);
        let pad_right = if cfg.w_padding > 0 {
            let ediff = seqlen.saturating_sub(t.end);
            (cfg.w_padding as u32).min(ediff)
        } else {
            0
        };
        let seq = get_unspliced(db, t, pad_left, pad_right).unwrap_or_default();
        if let Some(fu) = outs.u.as_mut() {
            write!(fu, ">{}", t.id).ok();
            if cfg.out_format == crate::config::OutFormat::Table {
                print_table_data(fu, t, cfg, true).ok();
            } else {
                write_attrs_gff(fu, t, cfg).ok();
            }
            writeln!(fu).ok();
            print_fasta(fu, None, &seq, false).ok();
        }
    }

    let _ = (inframe_stop, end_stop, stop_adjusted, full_cds);
    true
}

fn write_fasta_defline<W: Write>(
    cfg: &AppConfig,
    f: &mut W,
    t: &GffObj,
    _cds_seq: &[u8],
    _is_cds: bool,
) -> std::io::Result<()> {
    write!(f, ">{}", t.id)?;
    if cfg.write_exon_segs {
        write!(f, " loc:{}({}){}-{}", t.gseq_name, t.strand, t.start, t.end)?;
    }
    if cfg.out_format == crate::config::OutFormat::Table {
        print_table_data(f, t, cfg, true)?;
    } else {
        write_attrs_gff(f, t, cfg)?;
    }
    writeln!(f)?;
    Ok(())
}

fn write_attrs_gff<W: Write>(f: &mut W, t: &GffObj, cfg: &AppConfig) -> std::io::Result<()> {
    let decode = cfg.decode_chars;
    let attrs = crate::output::__attrs_for_fasta(cfg, t);
    if attrs.is_empty() {
        return Ok(());
    }
    write!(f, " ")?;
    let mut first = true;
    for (k, v) in attrs {
        if !first {
            write!(f, ";")?;
        }
        first = false;
        let val = if decode { crate::output::decode_hex_chars(v) } else { v.to_string() };
        write!(f, "{}={}", k, val)?;
    }
    Ok(())
}

fn build_w_defline(
    cfg: &AppConfig,
    t: &GffObj,
    res: &crate::fasta::SplicedResult,
    pad_left: u32,
    pad_right: u32,
) -> String {
    let mut d = String::new();
    d.push_str(&t.id);
    if !cfg.wfa_no_cds && t.has_cds() && res.cds_start > 0 {
        d.push_str(&format!(" CDS={}-{}", res.cds_start, res.cds_end));
    }
    if cfg.write_exon_segs {
        d.push_str(&format!(" loc:{}|{}-{}|{}", t.gseq_name, t.start, t.end, t.strand));
        d.push_str(" exons:");
        for (i, e) in t.exons.iter().enumerate() {
            if i > 0 {
                d.push(',');
            }
            d.push_str(&format!("{}-{}", e.start, e.end));
        }
        if cfg.w_padding > 0 {
            d.push_str(&format!(" padding:{}|{}", pad_left, pad_right));
        }
        d.push_str(" segs:");
        for (i, s) in res.segs.iter().enumerate() {
            if i > 0 {
                d.push(',');
            }
            d.push_str(&format!("{}-{}", s.start, s.end));
        }
    }
    d
}
