import argparse
import pysam
import pandas as pd
from collections import defaultdict
import os
import concurrent.futures
import textwrap

def load_gene_annotations_with_transcript_ends(gtf_file):
    genes = {}
    transcript_ends = {}
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)
            if feature == "gene":
                gene_id = get_attribute_value(attributes, "gene_id")
                key = (chrom, strand)
                genes.setdefault(key, []).append((start, end, gene_id))
            elif feature == "mRNA":
                transcript_id = get_attribute_value(attributes, "transcript_id")
                key = (chrom, strand)
                transcript_ends.setdefault(key, []).append((end, transcript_id))
    return genes, transcript_ends

def get_attribute_value(attributes, key):
    for attr in attributes.split(";"):
        if key in attr:
            return attr.split('"')[1]
    return None

def merge_polyA_sites(polyA_sites, merge_window, min_reads):
    merged_sites = []
    for (chrom, strand), positions in polyA_sites.items():
        positions.sort()
        start = positions[0]
        count = 1
        for i in range(1, len(positions)):
            if positions[i] - positions[i - 1] <= merge_window:
                count += 1
            else:
                if count >= min_reads:
                    merged_sites.append((chrom, start - 1, positions[i - 1], count, strand))
                start = positions[i]
                count = 1
        if count >= min_reads:
            merged_sites.append((chrom, start - 1, positions[-1], count, strand))
    return merged_sites

def process_bam_file(bam_file):
    polyA_sites = defaultdict(list)
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam.fetch(until_eof=True):
        if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
            if read.is_reverse:
                site = read.reference_start + 1
                strand = "-"
            else:
                site = read.reference_end
                strand = "+"
            polyA_sites[(read.reference_name, strand)].append(site)
    return polyA_sites

def identify_polyA_sites(bam_files, merge_window, min_reads, num_threads):
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        results = list(executor.map(process_bam_file, bam_files))
    merged_polyA_sites = defaultdict(list)
    for result in results:
        for key, sites in result.items():
            merged_polyA_sites[key].extend(sites)
    return merge_polyA_sites(merged_polyA_sites, merge_window, min_reads)

def assign_polyA_to_genes(polyA_regions, genes, apa_extension):
    assigned_regions = []
    for chrom, start, end, count, strand in polyA_regions:
        if (chrom, strand) not in genes:
            assigned_regions.append((chrom, start, end, count, strand, "no_gene", "no_relation"))
            continue
        closest_gene = None
        closest_distance = float("inf")
        relative_position = "no_relation"
        for gene_start, gene_end, gene_id in genes[(chrom, strand)]:
            polyA_site_start = gene_end - apa_extension if strand == "+" else gene_end - apa_extension
            polyA_site_end = gene_end + apa_extension if strand == "+" else gene_start + apa_extension

            if end < polyA_site_start:
                distance = polyA_site_start - end
                relation = "upstream" if strand == "+" else "downstream"
            elif start > polyA_site_end:
                distance = start - polyA_site_end
                relation = "downstream" if strand == "+" else "upstream"
            else:
                distance = 0
                relation = "within"
            if distance < closest_distance:
                closest_distance = distance
                closest_gene = gene_id
                relative_position = relation
        assigned_regions.append((chrom, start, end, count, strand, closest_gene or "no_gene", relative_position))
    return assigned_regions

def reassign_polyA_classes(assigned_regions):
    """
    Reassign poly(A) sites into three classes: most abundant, upstream, downstream.

    Parameters:
        assigned_regions (list): List of poly(A) regions with gene assignments.

    Returns:
        list: Updated regions with new classifications.
    """
    gene_to_sites = defaultdict(list)
    for region in assigned_regions:
        chrom, start, end, count, strand, gene_id, relation = region
        if gene_id != "no_gene":
            gene_to_sites[gene_id].append(region)

    updated_regions = []
    for gene_id, sites in gene_to_sites.items():
        # Sort by abundance (count) and position
        sites.sort(key=lambda x: (-x[3], x[1]))
        most_abundant = sites[0]
        updated_regions.append((*most_abundant[:6], "most_abundant"))

        for site in sites[1:]:
            _, start, end, _, strand, _, _ = site
            _, m_start, m_end, _, _, m_strand, _ = most_abundant

            
            relation = "upstream" if (end < m_start and strand == "+") or (end > m_end and strand == "-") else "downstream"



            updated_regions.append((*site[:6], relation))

    return updated_regions

def main():
    parser = argparse.ArgumentParser(description="Identify and analyze poly(A) sites from BAM files from transcriptomic Nanopore experiments (either DRS or cDNA).", epilog=textwrap.dedent("epilog"))
    parser.add_argument("-gtf", "--gtf_file", required=True, help="Path to the gene annotation GTF file.")
    parser.add_argument("-bam", "--bam_files", nargs='+', required=True, help="List of BAM files separated by space.")
    parser.add_argument("-o", "--output_bed", required=True, help="Path to the output BED file.")
    parser.add_argument("-w", "--merge_window", type=int, default=5, help="Maximum distance to merge sites.")
    parser.add_argument("-m", "--min_reads", type=int, default=5, help="Minimum reads to retain a poly(A) site.")
    parser.add_argument("-ae", "--annotation_extension", type=int, default=10, help="Number of bases for enlarging the annotation APA site in assign_polyA_to_genes.")
    parser.add_argument("-t", "--num_threads", type=int, default=1, help="Number of threads, for tasks that are multi-threaded.")
    parser.add_argument("-head", "--bed_header", action="store_true", help="Include header in the output BED file.")
    parser.add_argument("-nohead", "--no_bed_header", action="store_false", dest="bed_header", help="Exclude header from the output BED file.")
    parser.set_defaults(bed_header=True)


    args = parser.parse_args()

    merged_sites = identify_polyA_sites(args.bam_files, args.merge_window, args.min_reads, args.num_threads)
    print("Poly-A sites identified.")

    genes, transcript_ends = load_gene_annotations_with_transcript_ends(args.gtf_file)
    print("Annotation file with transcript ends loaded.")

    assigned_regions = assign_polyA_to_genes(merged_sites, genes, args.annotation_extension)
    reassigned_regions = reassign_polyA_classes(assigned_regions)
    print("Poly-A sites assigned and reassigned relative to the most abundant site.")

    bam_specific_counts = {bam: process_bam_file(bam) for bam in args.bam_files}

    extended_regions = []
    for region in reassigned_regions:
        chrom, start, end, count, strand, gene_id, relation = region

        overlapping_transcripts = [
            transcript_id for (trans_end, transcript_id)
            in transcript_ends.get((chrom, strand), [])
            if trans_end >= start and trans_end <= end
        ]
        transcript = overlapping_transcripts[0] if overlapping_transcripts else "none"

        bam_counts = [
            sum(1 for site in bam_specific_counts[bam][(chrom, strand)] if start <= site <= end)
            for bam in args.bam_files
        ]

        extended_regions.append((chrom, start, end, gene_id, count, strand, relation, transcript, *bam_counts))

    bam_file_columns = [os.path.basename(bam).replace(".bam", "") for bam in args.bam_files]
    columns = ["chrom", "start", "end", "gene_id", "total_count", "strand", "relation", "transcript", *bam_file_columns]
    counts_df = pd.DataFrame(extended_regions, columns=columns)

    counts_df.to_csv(args.output_bed, sep="\t", index=False, header=args.bed_header)
    print(f"Results saved to {args.output_bed}")

if __name__ == "__main__":
    main()
