process mergeSizeDistributions {
    cpus 1
    memory '2 GB'
    publishDir { "${params.outdir}/04_size_distribution" }, mode: 'symlink', overwrite: false

    input:
    path size_dist_counts_files
    path size_dist_rpm_files

    output:
    path "all_samples_size_distribution_counts.csv"
    path "all_samples_size_distribution_rpm.csv"

    script:
    """
    #!/usr/bin/env python3
    import csv
    import glob
    from collections import defaultdict

    def merge_files(pattern, output_file):
        all_rows  = []
        all_sizes = set()

        for csv_file in sorted(glob.glob(pattern)):
            with open(csv_file, newline="") as f:
                reader = csv.reader(f)
                header = next(reader)

                # CHANGED: header format is now:
                # library, genome, tissue, feature_type, 10, 11, ..., 50
                # size columns start at index 4 (previously index 2)
                if len(header) < 5:
                    continue

                sizes = [int(x) for x in header[4:]]    # CHANGED: offset 4
                all_sizes.update(sizes)

                for data_row in reader:
                    if len(data_row) < 5:
                        continue
                    library_name = data_row[0]
                    genome       = data_row[1]           # CHANGED: new column
                    tissue       = data_row[2]           # CHANGED: new column
                    feature_type = data_row[3]           # CHANGED: shifted index
                    values       = {}
                    for i, size in enumerate(sizes):
                        try:
                            values[size] = float(data_row[i + 4])  # CHANGED: offset 4
                        except (ValueError, IndexError):
                            values[size] = 0
                    all_rows.append((library_name, genome, tissue,
                                     feature_type, values))

        if not all_rows:
            with open(output_file, "w", newline="") as f:
                csv.writer(f).writerow(["library", "genome", "tissue",
                                        "feature_type"])
            return

        all_sizes_sorted = sorted(all_sizes)

        # header now includes genome and tissue columns
        header_out = ["library", "genome", "tissue", "feature_type"] + \
                     [str(s) for s in all_sizes_sorted]

        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header_out)
            # CHANGED: unpack five-element tuples including genome and tissue
            for library_name, genome, tissue, feature_type, values in all_rows:
                row = [library_name, genome, tissue, feature_type]
                row += [values.get(s, 0) for s in all_sizes_sorted]
                writer.writerow(row)

    merge_files("*_size_dist_counts.csv",
                "all_samples_size_distribution_counts.csv")
    merge_files("*_size_dist_rpm.csv",
                "all_samples_size_distribution_rpm.csv")
    """
}