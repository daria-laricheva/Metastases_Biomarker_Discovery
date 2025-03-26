# Define functions
def perform_fastqc(input_files_path, fastqc_out):
    """
    Perform quality control on all files in subdirectories of input_files_path using FastQC.
    
    Parameters:
    input_files_path (str): Path to the directory containing subdirectories of files.
    fastqc_out (str): Path to the directory where FastQC outputs will be saved.
    """
    # Iterate over each subdirectory in input_files_path
    for subdir in os.listdir(input_files_path):
        subdir_path = os.path.join(input_files_path, subdir)
        if not os.path.isdir(subdir_path):
            continue 
        
        # Create an output directory for the FastQC reports of this subdirectory
        output_subdir = os.path.join(fastqc_out, subdir)
        if not os.path.exists(output_subdir):
            os.makedirs(output_subdir)
        
        # Iterate over all files in the current subdirectory
        for filename in os.listdir(subdir_path):
            file_path = os.path.join(subdir_path, filename)
            # Run FastQC for the current file
            fastqc_cmd = ["fastqc", file_path, "-o", output_subdir]
            try:
                subprocess.run(fastqc_cmd, check=True)
                print(f"FastQC analysis completed for {filename}. Report saved in {output_subdir}")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while running FastQC on {filename}: {e}")

def perform_multiqc(fastqc_out, multiqc_out):
    """
    Run MultiQC on the FastQC output directory using Docker.
    
    Parameters:
    fastqc_out (str): Path to the directory containing FastQC outputs.
    multiqc_out (str): Path to the directory where MultiQC report will be saved.
    """
    # Make sure the output directory exists
    if not os.path.exists(multiqc_out):
        os.makedirs(multiqc_out)
    
    # Construct the Docker command
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{os.path.abspath(fastqc_out)}:/data",
        "-v", f"{os.path.abspath(multiqc_out)}:/output",
        "ewels/multiqc", "multiqc", "/data", "-o", "/output"
    ]
    
    # Execute the Docker command
    try:
        subprocess.run(docker_cmd, check=True)
        print(f"MultiQC analysis completed. Report saved in {multiqc_out}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running MultiQC: {e}")

def build_index(reference_genome_path, index_prefix, data_type, aligner):
    """
    Build index files for a reference genome using the specified aligner and data type.
    
    Parameters:
    reference_genome_path (str): Path to the reference genome fasta file.
    index_prefix (str): The prefix for the generated index files.
    data_type (str): Type of data ('dna' or 'rna').
    aligner (str): The aligner tool to use for indexing ('bowtie2', 'hisat2').
    """
    index_dir = os.path.join("indexes", data_type)  # Organize indexes by data type
    full_index_path = os.path.join(index_dir, index_prefix)

    # Ensure the index directory exists
    if not os.path.exists(index_dir):
        os.makedirs(index_dir)
        print(f"Created directory for {data_type} index files at: {index_dir}")

    if aligner == "bowtie2" and data_type == "dna":
        cmd = f"bowtie2-build {reference_genome_path} {full_index_path}"
    elif aligner == "hisat2" and data_type == "rna":
        cmd = f"hisat2-build {reference_genome_path} {full_index_path}"
    else:
        raise ValueError(f"Unsupported aligner/data type combination specified: {aligner}/{data_type}")

    # Execute the command
    try:
        print(f"Building {aligner} index for {data_type} data...")
        subprocess.run(cmd, shell=True, check=True)
        print(f"Indexing with {aligner} for {data_type} data completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during {aligner} indexing for {data_type} data: {e}")

def align_rna_hisat2(input_files_path, aligned_files_base_path, index_prefix):
    """
    Align paired-end RNA-seq reads using HISAT2.

    Parameters:
    input_files_path (str): Path to the directory containing paired-end RNA-seq FASTQ files.
    aligned_files_base_path (str): Base path for where aligned RNA-seq files will be saved.
    index_prefix (str): Full path with the prefix to the HISAT2 index for RNA.
    """
    # Create the rna directory in the aligned_files_base_path if it doesn't exist
    rna_aligned_files_path = os.path.join(aligned_files_base_path, "rna")
    if not os.path.exists(rna_aligned_files_path):
        os.makedirs(rna_aligned_files_path)
        
    # List all files in the input_files_path
    for filename in os.listdir(input_files_path):
        if filename.endswith("_1.fastq.gz"):
            # Construct the paths for paired-end FASTQ files
            file_path_1 = os.path.join(input_files_path, filename)
            file_path_2 = file_path_1.replace("_1.fastq.gz", "_2.fastq.gz")
            # Define the output SAM file path
            output_file = os.path.join(rna_aligned_files_path, filename.replace("_1.fastq.gz", ".sam"))
            
            # Prepare and execute the HISAT2 alignment command
            align_cmd = f"hisat2 -x {index_prefix} -1 {file_path_1} -2 {file_path_2} -S {output_file}"
            
            try:
                subprocess.run(align_cmd, shell=True, check=True)
                print(f"Alignment completed for {filename}. Output saved to {output_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while aligning {filename}: {e}")

def sam_to_sorted_bam(sam_file_path, output_dir):
    """
    Convert a SAM file to a sorted BAM file and index it.
    
    Parameters:
    sam_file_path (str): The path to the SAM file.
    output_dir (str): The directory where the BAM file will be saved.
    """
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Extract the base name for the file (without extension)
    base_name = os.path.basename(sam_file_path).replace(".sam", "")
    
    # Construct file paths
    bam_file_path = os.path.join(output_dir, base_name + ".bam")
    sorted_bam_file_path = os.path.join(output_dir, base_name + ".sorted.bam")
    
    # Convert SAM to BAM
    cmd_convert = f"samtools view -bS {sam_file_path} > {bam_file_path}"
    # Sort BAM file
    cmd_sort = f"samtools sort {bam_file_path} -o {sorted_bam_file_path}"
    # Index sorted BAM file
    cmd_index = f"samtools index {sorted_bam_file_path}"
    
    try:
        print(f"Converting {sam_file_path} to BAM...")
        subprocess.run(cmd_convert, shell=True, check=True)
        print(f"Sorting BAM file {bam_file_path}...")
        subprocess.run(cmd_sort, shell=True, check=True)
        print(f"Indexing sorted BAM file {sorted_bam_file_path}...")
        subprocess.run(cmd_index, shell=True, check=True)
        print(f"Process completed for {sam_file_path}. Sorted and indexed BAM file is {sorted_bam_file_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")

def process_all_sam_files_in_subdirectories(parent_dir, output_base_dir):
    """
    Process all SAM files within each subdirectory of a given parent directory,
    converting them to sorted and indexed BAM files.
    
    Parameters:
    parent_dir (str): Parent directory containing subdirectories with SAM files.
    output_base_dir (str): Base directory to save sorted and indexed BAM files.
    """
    for subdir in os.listdir(parent_dir):
        input_dir = os.path.join(parent_dir, subdir)
        if os.path.isdir(input_dir):
            output_dir = os.path.join(output_base_dir, subdir)
            print(f"Processing SAM files in {input_dir}...")
            
            for filename in os.listdir(input_dir):
                if filename.endswith(".sam"):
                    sam_file_path = os.path.join(input_dir, filename)
                    sam_to_sorted_bam(sam_file_path, output_dir)


import os
import subprocess

# Define necessary paths
working_dir = "."
input_files_path = os.path.join(working_dir, "Input_Files", "Bulk_RNA-Seq_input_files")

fastqc_out = os.path.join(working_dir, "fastqc_results")
multiqc_out = os.path.join(working_dir, "multiqc_results")

trimmed_files_path = os.path.join(working_dir, "Trimmed_Files")

trimmed_files_path_rna = os.path.join(trimmed_files_path, "Bulk_RNA-Seq_input_files")


fastqc_trimmed = os.path.join(working_dir, "fastqc_trimmed")
multiqc_trimmed = os.path.join(working_dir, "multiqc_trimmed")

reference_genome_17 = os.path.join(working_dir, "reference_genome_files", "chr17.fa.gz")
reference_genome_17_decomp = os.path.join(working_dir, "reference_genome_files", "chr17.fa")

aligned_files = os.path.join(working_dir, "aligned_files")

indexes = os.path.join(working_dir, "indexes")

sorted_indexed_files = os.path.join(working_dir, "sorted_indexed_files")

sorted_indexed_files_rna = os.path.join(sorted_indexed_files, "rna", "rna")

gtf_chr_17 = os.path.join(working_dir, "reference_genome_files", "chr17.gtf")

counts_data = os.path.join(working_dir, "counts_data")

# Command to subset reference genome file to chr17 (because of memory constraints)
command = """
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | \
awk '/^>17/{flag=1} flag; /^>18/{flag=0; exit}' | \
gzip > chr17.fa.gz
"""

# Execute the command
try:
    subprocess.run(command, shell=True, check=True)
    print("Successfully extracted chr17.")
except subprocess.CalledProcessError as e:
    print(f"Error occurred: {e}")

# Perform QC with fastqc
perform_fastqc(input_files_path, fastqc_out)

# Perform reports' aggregation with multiqc
perform_multiqc(fastqc_out, multiqc_out)

# Build index using hisat2 (for further alignmnet of rna sequences)
build_index(reference_genome_17_decomp, "hisat2_index", "rna", "hisat2")

# Perform alignment of RNA sequences
align_rna_hisat2(input_files_path, aligned_files, os.path.join(indexes, "rna", "hisat2_index"))

# Convert .sam files to sorted and indexed .bam files
process_all_sam_files_in_subdirectories(aligned_files, sorted_indexed_files_rna)

# Create the featureCounts command
feature_counts_command = f"""
featureCounts -a {gtf_chr_17} \
               -o {counts_data}/combined_counts.txt \
               -g gene_name \
               -p \
               -t exon \
               -s 2 \
               -T 4 \
               {sorted_indexed_files_rna}/*.sorted.bam
"""

# Execute the featureCounts command
try:
    subprocess.run(feature_counts_command, shell=True, check=True)
    print(f"featureCounts has successfully generated the counts file at {counts_data}.")
except subprocess.CalledProcessError as e:
    print(f"An error occurred while running featureCounts: {e}")
