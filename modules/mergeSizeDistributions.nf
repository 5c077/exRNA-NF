process mergeSizeDistributions {
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/04_size_distribution", mode: 'symlink', overwrite: false

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
        all_rows  = []   # list of (library, category, {size: value}) dicts
        all_sizes = set()

        for csv_file in sorted(glob.glob(pattern)):
            with open(csv_file, newline="") as f:
                reader = csv.reader(f)
                header = next(reader)

                # Header format: library, Category, 10, 11, 12, ...
                if len(header) < 3:
                    continue

                sizes = [int(x) for x in header[2:]]
                all_sizes.update(sizes)

                for data_row in reader:
                    if len(data_row) < 3:
                        continue
                    library_name = data_row[0]
                    category     = data_row[1]
                    values       = {}
                    for i, size in enumerate(sizes):
                        try:
                            values[size] = float(data_row[i + 2])
                        except (ValueError, IndexError):
                            values[size] = 0
                    all_rows.append((library_name, category, values))

        if not all_rows:
            with open(output_file, "w", newline="") as f:
                csv.writer(f).writerow(["library", "Category"])
            return

        all_sizes_sorted = sorted(all_sizes)
        header_out       = ["library", "Category"] + \
                           [str(s) for s in all_sizes_sorted]

        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header_out)
            for library_name, category, values in all_rows:
                row = [library_name, category]
                row += [values.get(s, 0) for s in all_sizes_sorted]
                writer.writerow(row)

    merge_files("*_size_dist_counts.csv", "all_samples_size_distribution_counts.csv")
    merge_files("*_size_dist_rpm.csv",    "all_samples_size_distribution_rpm.csv")
    """
}