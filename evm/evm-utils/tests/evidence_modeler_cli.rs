//! The evidence_modeler binary writes the same evm.out to stdout (how
//! funannotate calls it) and to `-o`, and both match the Perl golden file.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../tests/fixtures")
        .join(name)
}

fn args() -> Vec<String> {
    vec![
        "-G".into(),
        fixture("elm.genome.fasta").to_string_lossy().into_owned(),
        "-w".into(),
        fixture("elm.weights.txt").to_string_lossy().into_owned(),
        "-g".into(),
        fixture("elm.gene_predictions.gff3")
            .to_string_lossy()
            .into_owned(),
        "--report_ELM".into(),
    ]
}

#[test]
fn stdout_and_output_file_match_perl_golden() {
    let bin = env!("CARGO_BIN_EXE_evidence_modeler");
    let want = fs::read(fixture("Contig_elm.perl.evm.out")).unwrap();

    let out = Command::new(bin).args(args()).output().unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(out.stdout, want, "stdout differs from the Perl golden");

    let path = std::env::temp_dir().join(format!("evm_cli_{}.out", std::process::id()));
    let status = Command::new(bin)
        .args(args())
        .arg("-o")
        .arg(&path)
        .status()
        .unwrap();
    assert!(status.success());
    assert_eq!(
        fs::read(&path).unwrap(),
        want,
        "-o file differs from golden"
    );
    let _ = fs::remove_file(&path);
}
