process extractSizeDistribution {
    tag "${bam.simpleName}"
    cpus 2
    memory '4 GB'
    publishDir "${params.outdir}/04_size_distribution", mode: 'symlink', overwrite: false

    input:
    path bam

    output:
    path "${bam.simpleName}_size_dist_counts.csv", emit: size_dist_counts
    path "${bam.simpleName}_size_dist_rpm.csv", emit: size_dist_rpm

    script:
    """
    #!/usr/bin/env python3
    import pysam
    import csv
    from collections import defaultdict

    # Extract read lengths from BAM
    size_counts = defaultdict(int)
    
    try:
        bam_file = pysam.AlignmentFile("${bam}", "rb", check_sq=False)
        
        for read in bam_file:
            if not read.is_unmapped:
                read_length = read.query_length
                size_counts[read_length] += 1
        
        bam_file.close()
    except Exception as e:
        print(f"Error reading BAM file: {e}")
        # Create empty outputs so pipeline doesn't fail
        with open("${bam.simpleName}_size_dist_counts.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['library'])
            writer.writerow(["${bam.simpleName}"])
        with open("${bam.simpleName}_size_dist_rpm.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['library'])
            writer.writerow(["${bam.simpleName}"])
        raise
    
    # Calculate total reads for RPM normalization
    total_reads = sum(size_counts.values())
    library_name = "${bam.simpleName}"
    
    # Write counts CSV
    with open(f"{library_name}_size_dist_counts.csv", 'w', newline='') as csvfile:
        if size_counts:
            min_size = min(size_counts.keys())
            max_size = max(size_counts.keys())
            
            # Write header
            writer = csv.writer(csvfile)
            writer.writerow(['library'] + list(range(min_size, max_size + 1)))
            
            # Write data row
            row = [library_name]
            for size in range(min_size, max_size + 1):
                row.append(size_counts.get(size, 0))
            writer.writerow(row)
        else:
            # Empty file case
            writer = csv.writer(csvfile)
            writer.writerow(['library'])
            writer.writerow([library_name])
    
    # Write RPM CSV
    with open(f"{library_name}_size_dist_rpm.csv", 'w', newline='') as csvfile:
        if size_counts and total_reads > 0:
            min_size = min(size_counts.keys())
            max_size = max(size_counts.keys())
            
            # Write header
            writer = csv.writer(csvfile)
            writer.writerow(['library'] + list(range(min_size, max_size + 1)))
            
            # Write data row with RPM values
            row = [library_name]
            for size in range(min_size, max_size + 1):
                count = size_counts.get(size, 0)
                rpm = (count / total_reads) * 1_000_000 if total_reads > 0 else 0
                row.append(round(rpm, 2))
            writer.writerow(row)
        else:
            # Empty file case
            writer = csv.writer(csvfile)
            writer.writerow(['library'])
            writer.writerow([library_name])
    """
}