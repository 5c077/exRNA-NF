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
    from collections import defaultdict
    import glob

    def merge_files(pattern, output_file):
        # Read all individual CSV files
        all_data = []
        all_sizes = set()
        
        for csv_file in glob.glob(pattern):
            with open(csv_file, 'r') as f:
                reader = csv.reader(f)
                header = next(reader)
                
                if len(header) > 1:  # Not empty
                    sizes = [int(x) for x in header[1:]]
                    all_sizes.update(sizes)
                    
                    data_row = next(reader)
                    library_name = data_row[0]
                    
                    # Handle both integer counts and float RPM values
                    values = {}
                    for i in range(len(sizes)):
                        try:
                            values[sizes[i]] = float(data_row[i+1])
                        except (ValueError, IndexError):
                            values[sizes[i]] = 0
                    
                    all_data.append((library_name, values))
        
        # Write merged CSV
        if all_data:
            all_sizes_sorted = sorted(list(all_sizes))
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['library'] + all_sizes_sorted)
                
                for library_name, values in all_data:
                    row = [library_name]
                    for size in all_sizes_sorted:
                        row.append(values.get(size, 0))
                    writer.writerow(row)
        else:
            # Create empty file
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['library'])
    
    # Merge counts files
    merge_files("*_size_dist_counts.csv", "all_samples_size_distribution_counts.csv")
    
    # Merge RPM files
    merge_files("*_size_dist_rpm.csv", "all_samples_size_distribution_rpm.csv")
    """
}