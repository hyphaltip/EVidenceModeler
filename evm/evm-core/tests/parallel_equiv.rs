//! The parallel convert and recombine steps must give the same bytes as a
//! serial run, in partitions-listing order, including across batch boundaries.

use std::fs;
use std::path::{Path, PathBuf};

use evm_core::gff3_convert::evm_to_gff3::{
    convert_all_to_gff3_opts, evm_output_to_gff3, ConvertOptions,
};
use evm_core::io::partitions::PartitionEntry;
use evm_core::recombine::recombine::recombine_outputs;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../tests/fixtures")
        .join(name)
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("evm_par_{}_{}", name, std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn entry(acc: &str, dir: &Path) -> PartitionEntry {
    PartitionEntry {
        accession: acc.to_string(),
        base_dir: dir.to_string_lossy().into_owned(),
        is_partitioned: false,
        partition_dir: None,
    }
}

/// More entries than one parallel batch (4096), so batch order is exercised.
const N: usize = 4300;

#[test]
fn parallel_convert_matches_serial_in_listing_order() {
    let root = scratch("convert");
    let evm_text = fs::read_to_string(fixture("Contig1.perl.evm.out")).unwrap();
    let mut entries = Vec::new();
    for i in 0..N {
        let d = root.join(format!("c{:05}", i));
        fs::create_dir_all(&d).unwrap();
        // Every 7th contig is empty; the rest hold the fixture model.
        let body = if i % 7 == 0 { "" } else { evm_text.as_str() };
        fs::write(d.join("evm.out"), body).unwrap();
        entries.push(entry(&format!("Contig{}", i), &d));
    }

    // Serial reference: convert each entry alone, concatenate in order,
    // dropping the empty ones (what --skip_empty is meant to do).
    let mut want = String::new();
    for (i, e) in entries.iter().enumerate() {
        if i % 7 == 0 {
            continue;
        }
        let out = root.join(format!("ref{}.gff3", i));
        evm_output_to_gff3(
            &format!("{}/evm.out", e.base_dir),
            &e.accession,
            out.to_str().unwrap(),
        )
        .unwrap();
        want.push_str(&fs::read_to_string(&out).unwrap());
        fs::remove_file(&out).unwrap();
    }

    let concat = root.join("all.gff3");
    let opts = ConvertOptions {
        skip_empty: true,
        concat_out: Some(concat.to_string_lossy().into_owned()),
    };
    convert_all_to_gff3_opts(&entries, "evm.out", &opts).unwrap();

    assert_eq!(
        fs::read_to_string(&concat).unwrap(),
        want,
        "concatenated GFF3 differs from the serial reference"
    );
    // Per-directory files: written for non-empty, absent for empty.
    for (i, e) in entries.iter().enumerate() {
        let gff = Path::new(&e.base_dir).join("evm.out.gff3");
        assert_eq!(gff.exists(), i % 7 != 0, "entry {}", i);
    }
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn duplicate_accession_converts_once_with_last_base_dir() {
    let root = scratch("dup");
    let evm_text = fs::read_to_string(fixture("Contig1.perl.evm.out")).unwrap();
    let (a, b) = (root.join("a"), root.join("b"));
    for d in [&a, &b] {
        fs::create_dir_all(d).unwrap();
        fs::write(d.join("evm.out"), &evm_text).unwrap();
    }
    let entries = vec![entry("Contig1", &a), entry("Contig1", &b)];
    convert_all_to_gff3_opts(&entries, "evm.out", &ConvertOptions::default()).unwrap();
    assert!(!a.join("evm.out.gff3").exists());
    assert!(b.join("evm.out.gff3").exists());
    let _ = fs::remove_dir_all(&root);
}

#[test]
fn parallel_recombine_matches_one_contig_at_a_time() {
    let root = scratch("recombine");
    let evm_text = fs::read_to_string(fixture("Contig1.perl.evm.out")).unwrap();
    let mut all = Vec::new();
    let mut per_contig: Vec<Vec<PartitionEntry>> = Vec::new();
    for c in 0..64 {
        let base = root.join(format!("ctg{:03}", c));
        let mut mine = Vec::new();
        for (lend, rend) in [(1u32, 50_000u32), (50_001, 100_000)] {
            let pdir = base.join(format!("ctg{:03}_{}-{}", c, lend, rend));
            fs::create_dir_all(&pdir).unwrap();
            fs::write(pdir.join("evm.out"), &evm_text).unwrap();
            mine.push(PartitionEntry {
                accession: format!("ctg{:03}", c),
                base_dir: base.to_string_lossy().into_owned(),
                is_partitioned: true,
                partition_dir: Some(pdir.to_string_lossy().into_owned()),
            });
        }
        all.extend(mine.iter().cloned());
        per_contig.push(mine);
    }

    // Reference: recombine each contig on its own, keep the bytes, delete.
    let mut want = Vec::new();
    for mine in &per_contig {
        recombine_outputs(mine, "evm.out").unwrap();
        let p = Path::new(&mine[0].base_dir).join("evm.out");
        want.push(fs::read(&p).unwrap());
        fs::remove_file(&p).unwrap();
    }
    recombine_outputs(&all, "evm.out").unwrap();
    for (mine, want_bytes) in per_contig.iter().zip(&want) {
        let got = fs::read(Path::new(&mine[0].base_dir).join("evm.out")).unwrap();
        assert_eq!(&got, want_bytes, "contig {}", mine[0].accession);
    }
    assert!(!want[0].is_empty());
    let _ = fs::remove_dir_all(&root);
}
