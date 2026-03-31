process buildCombinedAnnotIndex {
    tag "${sample_id}"
    cpus 8
    memory '8 GB'
    publishDir { "${params.outdir}/annotation_indices/combined" }, mode: 'symlink', overwrite: false

    input:
    // tuple( sample_id, [ fasta1, fasta2, ... ] )
    tuple val(sample_id), path(fastas)

    output:
    // Emit the sample_id, labeled combined FASTA (needed for label parsing),
    // and all bowtie index files
     tuple val(sample_id),
        path("${sample_id}_combined.fa"),
        path("index/"),
        emit: index

    script:
    """
    # For each feature FASTA, prefix sequence headers with the feature type.
    # Feature type is derived from the filename: <sample_id>_<feature>.fa
    # e.g. ath_miRNA.fa -> feature type = "miRNA"
    # The labeled header format is: >featureType|originalId
    # This allows the quantification script to parse feature type directly from
    # the refence name in the BAM.

    > ${sample_id}_combined.fa

for fasta in ${fastas}; do

    [ -s "\$fasta" ] || continue

    feature=\$(basename \$fasta .fa | sed 's/^${sample_id}_//' | sed 's/^genome_//')

    python3 - "\$fasta" "\$feature" "${sample_id}_combined.fa" <<PYEOF
import sys

fasta_in  = sys.argv[1]
feature   = sys.argv[2]
fasta_out = sys.argv[3]

header    = None
seq_lines = []
written   = 0
skipped   = 0

def is_valid(seq_lines):
    seq = "".join(seq_lines).upper().strip()
    clean = seq.replace("N", "")
    return len(clean) > 0

with open(fasta_in) as fin, open(fasta_out, "a") as fout:
    for line in fin:
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                if is_valid(seq_lines):
                    fout.write(f">{feature}|{header}\\n")
                    fout.write("\\n".join(seq_lines) + "\\n")
                    written += 1
                else:
                    print(f"[SKIP] All-N sequence removed: {header}",
                          file=sys.stderr)
                    skipped += 1
            header    = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)

    # Final record handler is now INSIDE the with block -- fout is still open
    if header is not None:
        if is_valid(seq_lines):
            fout.write(f">{feature}|{header}\\n")
            fout.write("\\n".join(seq_lines) + "\\n")
            written += 1
        else:
            print(f"[SKIP] All-N sequence removed: {header}",
                  file=sys.stderr)
            skipped += 1

print(f"[{feature}] Written: {written}  Skipped (all-N): {skipped}",
      file=sys.stderr)
PYEOF

done

if [ ! -s "${sample_id}_combined.fa" ]; then
    echo "ERROR: Combined FASTA for ${sample_id} is empty after filtering."
    exit 1
fi

echo "Combined FASTA sequence count:"
grep -c "^>" ${sample_id}_combined.fa

# Attempt standard index build first, fall back to large index if it fails.

echo "Attempting standard bowtie-build..."

mkdir -p index

if bowtie2-build \
    --threads ${task.cpus} \
    ${sample_id}_combined.fa \
    index/${sample_id}_combined \
    2>${sample_id}_bowtie_build.log; then
    echo "Standard index built successfully."
else
    echo "bowtie2-build failed for ${sample_id}"
    cat ${sample_id}_bowtie_build.log
    exit 1
fi
    """
}