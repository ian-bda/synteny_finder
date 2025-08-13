import os
import re
import subprocess
import shutil
import pandas as pd
import json
import gzip

# Set up working directories
WORK_DIR = os.path.join(os.path.dirname(__file__), 'synteny_work')
BEDS_DIR = os.path.join(WORK_DIR, 'beds')
CDS_DIR = os.path.join(WORK_DIR, 'cds')
ANCHORS_DIR = os.path.join(WORK_DIR, 'anchors')
CSVS_DIR = os.path.join(WORK_DIR, 'csvs')
FINAL_OUTPUT_DIR = os.path.join(CSVS_DIR, 'final_synteny_csvs')
for d in [WORK_DIR, BEDS_DIR, CDS_DIR, ANCHORS_DIR, CSVS_DIR, FINAL_OUTPUT_DIR]:
    os.makedirs(d, exist_ok=True)

# Input directories
GFF_DIR = os.path.join(os.path.dirname(__file__), 'gff_files')
ORIG_CDS_DIR = os.path.join(os.path.dirname(__file__), 'cds_files')

# Helper: Extract species name from GFF filename (before first dot)
def extract_species_name(gff_filename):
    return gff_filename.split('.')[0]

def gff_to_bed(gff_path, bed_path):
    # Smart key selection based on GFF content
    # First check if 'gene' field exists and has meaningful values
    if has_gene_field(gff_path):
        try:
            subprocess.run([
                'python', '-m', 'jcvi.formats.gff', 'bed', '--type=mRNA', '--key=gene', gff_path, '-o', bed_path
            ], check=True)
            return
        except subprocess.CalledProcessError:
            pass
    
    # Try locus_tag if gene field doesn't work
    if has_locus_tag_field(gff_path):
        try:
            subprocess.run([
                'python', '-m', 'jcvi.formats.gff', 'bed', '--type=mRNA', '--key=locus_tag', gff_path, '-o', bed_path
            ], check=True)
            return
        except subprocess.CalledProcessError:
            pass
    
    # Final fallback to ID
    subprocess.run([
        'python', '-m', 'jcvi.formats.gff', 'bed', '--type=mRNA', '--key=ID', gff_path, '-o', bed_path
    ], check=True)

def has_gene_field(gff_path):
    """Check if GFF has meaningful gene field values"""
    open_func = gzip.open if gff_path.endswith('.gz') else open
    with open_func(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            # Check if gene field has actual gene names, not just empty or generic values
            gene_match = re.search(r'gene=([^;\s]+)', line)
            if gene_match and gene_match.group(1) and gene_match.group(1) != '.':
                # Skip mitochondrial genes and other non-standard identifiers
                gene_name = gene_match.group(1)
                if not gene_name.startswith(('ND', 'COX', 'ATP', 'CYTB', 'MT-', 'MT')):
                    return True
    return False

def has_locus_tag_field(gff_path):
    """Check if GFF has locus_tag field"""
    open_func = gzip.open if gff_path.endswith('.gz') else open
    with open_func(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if 'locus_tag=' in line:
                return True
    return False

def format_cds_fasta(cds_path, formatted_cds_path):
    subprocess.run([
        'python', '-m', 'jcvi.formats.fasta', 'format', cds_path, formatted_cds_path
    ], check=True)

def run_ortholog_search(species1, species2, bed1, bed2, cds1, cds2, anchors_dir):
    # Copy files to cwd with expected names for jcvi
    bed1_target = f"{species1}_parsed.bed"
    bed2_target = f"{species2}_parsed.bed"
    cds1_target = f"{species1}_parsed.cds"
    cds2_target = f"{species2}_parsed.cds"
    shutil.copyfile(bed1, bed1_target)
    shutil.copyfile(bed2, bed2_target)
    shutil.copyfile(cds1, cds1_target)
    shutil.copyfile(cds2, cds2_target)
    anchors_path = None
    try:
        subprocess.run([
            'python', '-m', 'jcvi.compara.catalog', 'ortholog', f'{species1}_parsed', f'{species2}_parsed', '--no_strip_names'
        ], check=True)
        # Check for anchors file in current directory first
        anchors_file = f"{species1}_parsed.{species2}_parsed.anchors"
        if not os.path.exists(anchors_file):
            anchors_file = f"{species2}_parsed.{species1}_parsed.anchors"
        if os.path.exists(anchors_file):
            anchors_path = anchors_file  # Return the path to the file in current directory
        else:
            # If not found in current directory, check if it was moved to anchors_dir
            anchors_file_in_dir = os.path.join(anchors_dir, f"{species1}_parsed.{species2}_parsed.anchors")
            if not os.path.exists(anchors_file_in_dir):
                anchors_file_in_dir = os.path.join(anchors_dir, f"{species2}_parsed.{species1}_parsed.anchors")
            if os.path.exists(anchors_file_in_dir):
                anchors_path = anchors_file_in_dir
            else:
                anchors_path = None
    except subprocess.CalledProcessError:
        print(f"Warning: No anchors found for {species1} vs {species2}, skipping.")
        anchors_path = None
    finally:
        # Clean up temp files
        for f in [bed1_target, bed2_target, cds1_target, cds2_target]:
            if os.path.exists(f):
                os.remove(f)
    return anchors_path

def build_id_to_name_map(gff_file):
    id_to_name = {}
    open_func = gzip.open if gff_file.endswith('.gz') else open
    with open_func(gff_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            feature_type = fields[2]
            if feature_type not in ('mRNA', 'transcript'):
                continue
            attr = fields[8]
            # Try to get ID and gene_name or Name
            id_match = re.search(r'ID=([^;]+)', attr)
            name_match = re.search(r'gene_name=([^;]+)', attr)
            if not name_match:
                name_match = re.search(r'Name=([^;]+)', attr)
            if id_match and name_match:
                ensembl_id = id_match.group(1).replace('transcript:', '').split('.')[0]
                gene_name = name_match.group(1)
                id_to_name[ensembl_id] = gene_name
    return id_to_name

def parse_anchors_to_csv(anchors_file, bed1_file, bed2_file, output_dir, id_to_name1=None, id_to_name2=None):
    bed1_df = pd.read_csv(bed1_file, sep="\t", header=None, names=["chrom", "start", "end", "gene", "score", "strand"])
    bed2_df = pd.read_csv(bed2_file, sep="\t", header=None, names=["chrom", "start", "end", "gene", "score", "strand"])
    bed1_df = bed1_df[["chrom", "start", "end", "gene"]]
    bed2_df = bed2_df[["chrom", "start", "end", "gene"]]
    anchors_df = pd.read_csv(anchors_file, sep="\t", header=None, names=["gene1", "gene2", "score"])
    anchors_df['gene1'] = anchors_df['gene1'].str.strip()
    anchors_df['gene2'] = anchors_df['gene2'].str.strip()
    bed1_df['gene'] = bed1_df['gene'].str.strip()
    bed2_df['gene'] = bed2_df['gene'].str.strip()
    merged_df1 = pd.merge(anchors_df, bed1_df, left_on="gene1", right_on="gene", how="left")
    merged_df1 = merged_df1.rename(columns={"chrom": "chrom1", "start": "start1", "end": "end1"}).drop(columns=["gene"])
    merged_df2 = pd.merge(merged_df1, bed2_df, left_on="gene2", right_on="gene", how="left")
    merged_df2 = merged_df2.rename(columns={"chrom": "chrom2", "start": "start2", "end": "end2"}).drop(columns=["gene"])
    merged_df2 = merged_df2.dropna(subset=["chrom1", "chrom2"])
    # Replace Ensembl IDs with gene names if mapping provided
    if id_to_name1 is not None:
        merged_df2['gene1'] = merged_df2['gene1'].map(lambda x: id_to_name1.get(x, x))
    if id_to_name2 is not None:
        merged_df2['gene2'] = merged_df2['gene2'].map(lambda x: id_to_name2.get(x, x))
    output_file = os.path.join(output_dir, f"synteny_{os.path.basename(anchors_file).replace('.anchors', '.csv')}")
    merged_df2.to_csv(output_file, sep="\t", index=False)
    print(f"Saved synteny file: {output_file}")
    return output_file

def filter_csv_by_chromosomes(parsed_csv, sp1, sp2, chrom_map, output_dir):
    df = pd.read_csv(parsed_csv, sep='\t')
    chrom1 = chrom_map.get(sp1)
    chrom2 = chrom_map.get(sp2)
    if chrom1 is None or chrom2 is None:
        print(f"Warning: No chromosome mapping for {sp1} or {sp2}")
        return None
    filtered = df[(df['chrom1'].astype(str) == str(chrom1)) & (df['chrom2'].astype(str) == str(chrom2))]
    outname = os.path.join(output_dir, os.path.basename(parsed_csv))
    filtered.to_csv(outname, sep='\t', index=False)
    print(f"Filtered synteny file: {outname}")
    return outname

def normalize_bed_ids(bed_path):
    tmp_path = bed_path + ".tmp"
    with open(bed_path) as infile, open(tmp_path, "w") as outfile:
        for line in infile:
            fields = line.rstrip().split('\t')
            if len(fields) >= 4:
                # Extract the locus tag from the gene ID (e.g., rna-KC01_LOCUS1 -> KC01_LOCUS1)
                gene_id = fields[3]
                if 'rna-' in gene_id:
                    gene_id = gene_id.replace('rna-', '')
                elif 'transcript:' in gene_id:
                    gene_id = gene_id.replace('transcript:', '').split('.')[0]
                fields[3] = gene_id
            outfile.write('\t'.join(fields) + '\n')
    os.replace(tmp_path, bed_path)

def normalize_cds_ids(cds_path):
    tmp_path = cds_path + ".tmp"
    with open(cds_path) as infile, open(tmp_path, "w") as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract the gene name from the header
                header = line[1:]
                gene_match = re.search(r'\[gene=([^\]]+)\]', header)
                if gene_match:
                    gene_id = gene_match.group(1)
                else:
                    # Fallback: try locus_tag
                    locus_match = re.search(r'\[locus_tag=([^\]]+)\]', header)
                    if locus_match:
                        gene_id = locus_match.group(1)
                    else:
                        # Final fallback: extract from the first part
                        gene_id = header.split()[0].split('|')[-1].split('.')[0]
                
                # Keep the original header but replace the ID
                rest = ' '.join(header.split()[1:])
                outfile.write(f'>{gene_id} {rest}\n' if rest else f'>{gene_id}\n')
            else:
                outfile.write(line)
    os.replace(tmp_path, cds_path)

def check_bed_cds_overlap(bed_path, cds_path, species):
    bed_ids = set()
    with open(bed_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                bed_ids.add(fields[3])
    cds_ids = set()
    with open(cds_path) as f:
        for line in f:
            if line.startswith('>'):
                cds_ids.add(line[1:].split()[0])
    overlap = bed_ids & cds_ids
    print(f"[ID CHECK] {species}: BED IDs: {len(bed_ids)}, CDS IDs: {len(cds_ids)}, Overlap: {len(overlap)}")

if __name__ == '__main__':
    # Load chromosome map
    with open('chromosome_map.json') as f:
        chrom_map = json.load(f)

    # 1. Find all GFF and CDS files, match by species
    gff_files = sorted([f for f in os.listdir(GFF_DIR) if f.endswith('.gff')])
    cds_files = sorted([f for f in os.listdir(ORIG_CDS_DIR) if f.endswith('_cds.fna')])
    species_info = {}
    id_to_name_map = {}
    for gff in gff_files:
        species = extract_species_name(gff)
        cds = next((c for c in cds_files if species in c), None)
        if not cds:
            print(f"Warning: No CDS file found for {species}")
            continue
        species_info[species] = {
            'gff': os.path.join(GFF_DIR, gff),
            'cds': os.path.join(ORIG_CDS_DIR, cds)
        }
        # Build ID to gene name map for this species
        id_to_name_map[species] = build_id_to_name_map(os.path.join(GFF_DIR, gff))

    # 2. For each species, process GFF to BED, format CDS
    for species, info in species_info.items():
        print(f"\nProcessing {species} ...")
        bed = os.path.join(BEDS_DIR, f"{species}_parsed.bed")
        formatted_cds = os.path.join(CDS_DIR, f"{species}.cds")
        if not os.path.exists(bed):
            gff_to_bed(info['gff'], bed)
        normalize_bed_ids(bed)
        if not os.path.exists(formatted_cds):
            format_cds_fasta(info['cds'], formatted_cds)
        normalize_cds_ids(formatted_cds)
        check_bed_cds_overlap(bed, formatted_cds, species)
        info['bed'] = bed
        info['cds_filtered'] = formatted_cds

    # 3. For all species pairs, run ortholog search and parse anchors
    species_list = list(species_info.keys())
    for i in range(len(species_list)):
        for j in range(i+1, len(species_list)):
            sp1 = species_list[i]
            sp2 = species_list[j]
            print(f"\nRunning ortholog search: {sp1} vs {sp2}")
            anchors_path = run_ortholog_search(
                sp1, sp2,
                species_info[sp1]['bed'], species_info[sp2]['bed'],
                species_info[sp1]['cds_filtered'], species_info[sp2]['cds_filtered'],
                ANCHORS_DIR
            )
            # Move all relevant output files to anchors/<sp1>__<sp2>/
            pair_dir = os.path.join(ANCHORS_DIR, f"{sp1}__{sp2}")
            os.makedirs(pair_dir, exist_ok=True)
            base_patterns = [
                f"{sp1}_parsed.{sp2}_parsed",
                f"{sp2}_parsed.{sp1}_parsed"
            ]
            for pat in base_patterns:
                for ext in [".last", ".last.filtered", ".anchors", ".pdf", ".lifted.anchors"]:
                    fname = f"{pat}{ext}"
                    if os.path.exists(fname):
                        shutil.move(fname, os.path.join(pair_dir, fname))
            # Check if anchors file exists in the pair-specific subdirectory
            pair_dir = os.path.join(ANCHORS_DIR, f"{sp1}__{sp2}")
            anchors_file_in_pair = os.path.join(pair_dir, f"{sp1}_parsed.{sp2}_parsed.anchors")
            if not os.path.exists(anchors_file_in_pair):
                anchors_file_in_pair = os.path.join(pair_dir, f"{sp2}_parsed.{sp1}_parsed.anchors")
            
            if os.path.exists(anchors_file_in_pair):
                # Parse the anchors file from the pair directory
                parsed_csv = parse_anchors_to_csv(
                    anchors_file_in_pair,
                    species_info[sp1]['bed'],
                    species_info[sp2]['bed'],
                    CSVS_DIR,
                    id_to_name_map[sp1],
                    id_to_name_map[sp2]
                )
                filter_csv_by_chromosomes(parsed_csv, sp1, sp2, chrom_map, FINAL_OUTPUT_DIR)
            else:
                print(f"Warning: No anchors file found for {sp1} vs {sp2}") 
