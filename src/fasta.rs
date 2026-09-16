//! FASTA database with offset indexing and spliced-sequence extraction
//! (replaces gclib's `GFastaDb` / `GFaSeqGet` / `GffObj::getSpliced`).

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;

use anyhow::{anyhow, Result};

use crate::types::{reverse_complement, GSeg, GffObj};

/// Indexed location of one sequence inside a FASTA file.
#[derive(Debug, Clone, Copy)]
struct SeqEntry {
    /// byte offset in the file where the sequence data begins (right after the
    /// header line's terminating newline)
    offset: u64,
    /// number of bytes of sequence (including newlines) until the next header.
    raw_len: u64,
}

/// A FASTA database backed by a single file or a directory of single-fasta files.
pub struct FastaDb {
    path: PathBuf,
    is_dir: bool,
    /// file index (only for single-file mode).
    index: HashMap<String, SeqEntry>,
    indexed: bool,
    /// whole raw file contents (single-file mode), used for lazy per-sequence
    /// extraction. Dropped once it is no longer needed to save memory.
    data: Option<Vec<u8>>,
    /// cache of already-loaded, newline-stripped sequences.
    cache: HashMap<String, Vec<u8>>,
    seq_lens: HashMap<String, u32>,
}

impl FastaDb {
    pub fn new(path: &str) -> Result<Self> {
        let p = PathBuf::from(path);
        let is_dir = p.is_dir();
        if !is_dir && !p.exists() && path != "-" {
            return Err(anyhow!("FASTA path does not exist: {}", path));
        }
        Ok(FastaDb {
            path: p,
            is_dir,
            index: HashMap::new(),
            indexed: false,
            data: None,
            cache: HashMap::new(),
            seq_lens: HashMap::new(),
        })
    }

    /// Build the index for single-file mode (scans headers only).
    ///
    /// Reads the whole file into memory once — gffread's `-g` FASTA inputs are
    /// typically modest in size, and this avoids subtle seek/buffer bugs.
    fn ensure_indexed(&mut self) -> Result<()> {
        if self.indexed || self.is_dir {
            self.indexed = true;
            return Ok(());
        }
        let mut f = File::open(&self.path)?;
        let mut data = Vec::new();
        f.read_to_end(&mut data)?;

        let n = data.len();
        let mut pos: u64 = 0;
        let mut cur_name: Option<String> = None;
        let mut cur_start: u64 = 0;
        let mut i = 0usize;
        let mut at_line_start = true;
        while i < n {
            let b = data[i];
            if at_line_start && b == b'>' {
                // close previous record
                if let Some(name) = cur_name.take() {
                    let raw_len = pos - cur_start;
                    self.index.insert(name, SeqEntry { offset: cur_start, raw_len });
                }
                // read header line until newline
                let header_start = i + 1;
                let mut j = header_start;
                while j < n && data[j] != b'\n' {
                    j += 1;
                }
                let header = &data[header_start..j];
                let header_str = String::from_utf8_lossy(header);
                let name = header_str.split_whitespace().next().unwrap_or("").to_string();
                // skip past the newline (if any)
                if j < n {
                    j += 1;
                }
                i = j;
                pos = i as u64;
                cur_name = Some(name);
                cur_start = pos;
                at_line_start = true;
                continue;
            }
            at_line_start = b == b'\n';
            pos += 1;
            i += 1;
        }
        if let Some(name) = cur_name.take() {
            let raw_len = pos - cur_start;
            self.index.insert(name, SeqEntry { offset: cur_start, raw_len });
        }
        self.data = Some(data);
        self.indexed = true;
        Ok(())
    }

    /// Make sure the (newline-stripped) sequence `name` is present in the
    /// cache. Returns `false` if the sequence does not exist.
    ///
    /// This is the only place that pays the full-sequence cost (one read +
    /// one newline-strip per distinct sequence name). All `subseq` calls then
    /// borrow from the cache instead of re-cloning the whole sequence, which
    /// is critical for large genomes with many transcripts.
    fn load_seq_to_cache(&mut self, name: &str) -> Result<bool> {
        if self.cache.contains_key(name) {
            return Ok(true);
        }
        if self.is_dir {
            return self.get_seq_from_dir(name).map(|o| o.is_some());
        }
        self.ensure_indexed()?;
        let Some(entry) = self.index.get(name).copied() else {
            return Ok(false);
        };
        let data = self.data.as_ref().unwrap();
        let raw = &data[entry.offset as usize..(entry.offset + entry.raw_len) as usize];
        let seq: Vec<u8> = raw.iter().copied().filter(|&b| b != b'\n' && b != b'\r').collect();
        self.seq_lens.insert(name.to_string(), seq.len() as u32);
        self.cache.insert(name.to_string(), seq);
        Ok(true)
    }

    /// Load a sequence by name (1-based genomic coordinates externally).
    pub fn get_seq(&mut self, name: &str) -> Result<Option<Vec<u8>>> {
        if self.load_seq_to_cache(name)? {
            Ok(self.cache.get(name).cloned())
        } else {
            Ok(None)
        }
    }

    fn get_seq_from_dir(&mut self, name: &str) -> Result<Option<Vec<u8>>> {
        let candidates = [
            format!("{}.fa", name),
            format!("{}.fasta", name),
            format!("{}.fa.gz", name),
            format!("{}.fasta.gz", name),
            format!("{}.fna", name),
            name.to_string(),
        ];
        for c in &candidates {
            let p = self.path.join(c);
            if p.exists() {
                let f = File::open(&p)?;
                let lower = c.to_ascii_lowercase();
                let seq = if lower.ends_with(".gz") {
                    let mut dec = flate2::read::GzDecoder::new(f);
                    let mut buf = Vec::new();
                    dec.read_to_end(&mut buf)?;
                    parse_first_seq(&buf)
                } else {
                    let mut buf = Vec::new();
                    let mut f = f;
                    f.read_to_end(&mut buf)?;
                    parse_first_seq(&buf)
                };
                self.seq_lens.insert(name.to_string(), seq.len() as u32);
                self.cache.insert(name.to_string(), seq.clone());
                return Ok(Some(seq));
            }
        }
        Ok(None)
    }

    pub fn seq_len(&mut self, name: &str) -> Result<u32> {
        if let Some(&l) = self.seq_lens.get(name) {
            return Ok(l);
        }
        if self.load_seq_to_cache(name)? {
            return Ok(*self.seq_lens.get(name).unwrap_or(&0));
        }
        Ok(0)
    }

    /// Genomic subsequence, 1-based `start`, length `len`. Returns raw bytes
    /// (forward strand, no orientation change).
    ///
    /// Only the requested window is copied; the full sequence lives in the
    /// cache and is loaded (read + newline-strip) at most once per name.
    pub fn subseq(&mut self, name: &str, start: u32, len: u32) -> Result<Option<Vec<u8>>> {
        if !self.load_seq_to_cache(name)? {
            return Ok(None);
        }
        let seq = self.cache.get(name).unwrap();
        if start == 0 || start > seq.len() as u32 {
            return Ok(Some(Vec::new()));
        }
        let s0 = (start - 1) as usize;
        let end = (s0 + len as usize).min(seq.len());
        Ok(Some(seq[s0..end].to_vec()))
    }
}

/// Parse the first sequence from an in-memory FASTA buffer.
fn parse_first_seq(buf: &[u8]) -> Vec<u8> {
    let mut seq = Vec::new();
    let mut in_seq = false;
    for &b in buf {
        if b == b'\n' || b == b'\r' {
            continue;
        }
        if !in_seq {
            if b == b'>' {
                in_seq = true;
                // skip rest of header line
                continue;
            }
            // bytes before any header — ignore
            continue;
        }
        if b == b'>' {
            break;
        }
        seq.push(b);
    }
    seq
}

/// Result of spliced sequence extraction.
#[derive(Debug, Clone, Default)]
pub struct SplicedResult {
    pub seq: Vec<u8>,
    /// CDS start/end offsets in the spliced sequence (1-based, 0 if none).
    pub cds_start: u32,
    pub cds_end: u32,
    /// Segments in spliced coordinates (exons for -w, CDS for -x).
    pub segs: Vec<GSeg>,
}

/// Extract spliced sequence. When `cds_only`, extract only the CDS (mRNA
/// orientation). Otherwise extract all exons and project CDS coords.
pub fn get_spliced(db: &mut FastaDb, t: &GffObj, cds_only: bool) -> Result<SplicedResult> {
    if t.exons.is_empty() {
        return Ok(SplicedResult::default());
    }
    let revc = t.strand == '-';
    let mut exons = t.exons.clone();
    exons.sort_by_key(|e| e.start);

    // gather forward-strand exon sequences
    let mut fwd: Vec<u8> = Vec::new();
    // map: (exon_idx) -> (offset_start in fwd, length)
    let mut exon_offsets: Vec<(usize, u32)> = Vec::new();
    for ex in &exons {
        let off = fwd.len();
        if let Some(s) = db.subseq(&t.gseq_name, ex.start, ex.len())? {
            fwd.extend_from_slice(&s);
        }
        exon_offsets.push((off, ex.len()));
    }

    if cds_only {
        if !t.has_cds() {
            return Ok(SplicedResult::default());
        }
        // extract CDS portions in genomic order, then RC if needed
        let mut cds_fwd: Vec<u8> = Vec::new();
        let mut segs: Vec<GSeg> = Vec::new();
        let mut offset = 0u32;
        for ex in &exons {
            if ex.end < t.cds_start || ex.start > t.cds_end {
                continue;
            }
            let cstart = ex.start.max(t.cds_start);
            let cend = ex.end.min(t.cds_end);
            let clen = if cend >= cstart { cend - cstart + 1 } else { 0 };
            if clen == 0 {
                continue;
            }
            if let Some(s) = db.subseq(&t.gseq_name, cstart, clen)? {
                cds_fwd.extend_from_slice(&s);
            }
            segs.push(GSeg::new(offset + 1, offset + clen));
            offset += clen;
        }
        let seq = if revc { reverse_complement(&cds_fwd) } else { cds_fwd };
        let total = seq.len() as u32;
        let segs = if revc {
            segs.into_iter()
                .map(|s| GSeg::new(total - s.end + 1, total - s.start + 1))
                .collect()
        } else {
            segs
        };
        let len = seq.len() as u32;
        Ok(SplicedResult {
            seq,
            cds_start: 1,
            cds_end: len,
            segs,
        })
    } else {
        // full spliced exons
        let segs: Vec<GSeg> = exon_offsets
            .iter()
            .map(|(off, len)| GSeg::new(*off as u32 + 1, *off as u32 + *len))
            .collect();
        // project CDS genomic coords onto forward spliced coords
        let mut cds_fwd_start = 0u32;
        let mut cds_fwd_end = 0u32;
        if t.has_cds() {
            for (i, ex) in exons.iter().enumerate() {
                let (off, _len) = exon_offsets[i];
                let off = off as u32;
                if t.cds_end >= ex.start && t.cds_start <= ex.end {
                    let cs = t.cds_start.max(ex.start);
                    let ce = t.cds_end.min(ex.end);
                    let local_s = off + (cs - ex.start) + 1;
                    let local_e = off + (ce - ex.start) + 1;
                    if cds_fwd_start == 0 {
                        cds_fwd_start = local_s;
                    }
                    cds_fwd_end = local_e;
                }
            }
        }
        let total = fwd.len() as u32;
        // Map forward-strand spliced coordinates to mRNA-orientation coordinates
        // (reverse strand flips the axis: fwd_pos -> total - fwd_pos + 1).
        let map_pos = |p: u32| -> u32 {
            if p == 0 {
                0
            } else if revc {
                total - p + 1
            } else {
                p
            }
        };
        let segs_final: Vec<GSeg> = segs
            .into_iter()
            .map(|s| GSeg::new(map_pos(s.start), map_pos(s.end)))
            .collect();
        // Project CDS bounds (forward -> mRNA orientation), then order them.
        let cs = map_pos(cds_fwd_start);
        let ce = map_pos(cds_fwd_end);
        let (cds_start, cds_end) = match (cs, ce) {
            (0, 0) => (0u32, 0u32),
            (a, 0) | (0, a) => (a, a),
            (a, b) => (a.min(b), a.max(b)),
        };
        let seq = if revc { reverse_complement(&fwd) } else { fwd };
        Ok(SplicedResult {
            seq,
            cds_start,
            cds_end,
            segs: segs_final,
        })
    }
}

/// Extract the unspliced transcript span (genomic, including introns), with
/// optional padding. Returns the forward-strand sequence.
pub fn get_unspliced(db: &mut FastaDb, t: &GffObj, pad_left: u32, pad_right: u32) -> Result<Vec<u8>> {
    let start = t.start.saturating_sub(pad_left).max(1);
    let len = t.end - start + 1 + pad_right;
    Ok(db.subseq(&t.gseq_name, start, len)?.unwrap_or_default())
}
