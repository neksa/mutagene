import os
from mutagene.io.context_window import read_MAF_with_context_window
from mutagene.motifs.motifs import identify_motifs
import collections

data_path = os.path.dirname(os.path.realpath(__file__))


def test_motif():
    MAF_file = open(data_path + "/test.maf")
    genome = open(data_path + "/test_genome.2bit")
    window_size = 40
    custom_motif = "A[A>G]A"  # Lys to Arg
    mutations, mutations_with_context, processing_stats = read_MAF_with_context_window(MAF_file, genome, window_size)
    matching_motifs = identify_motifs(mutations_with_context, custom_motif) if mutations_with_context is not None else []
    print(matching_motifs)


def generate_maf():
    fields = "Hugo_Symbol Entrez_Gene_Id  Center  NCBI_Build  Chromosome  Start_position  End_position    Strand  Variant_Classification  Variant_Type    Reference_Allele    Tumor_Seq_Allele1   Tumor_Seq_Allele2   dbSNP_RS    dbSNP_Val_Status    Tumor_Sample_Barcode    Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1  Match_Norm_Seq_Allele2  Tumor_Validation_Allele1    Tumor_Validation_Allele2    Match_Norm_Validation_Allele1   Match_Norm_Validation_Allele2   Verification_Status Validation_Status   Mutation_Status Sequencing_Phase    Sequence_Source Validation_Method   Score   BAM_file    Sequencer   Tumor_Sample_UUID   Matched_Norm_Sample_UUID    Genome_Change   Annotation_Transcript   Transcript_Strand   Transcript_Exon Transcript_Position cDNA_Change Codon_Change    Protein_Change  Other_Transcripts   Refseq_mRNA_Id  Refseq_prot_Id  SwissProt_acc_Id    SwissProt_entry_Id  Description UniProt_AApos   UniProt_Region  UniProt_Site    UniProt_Natural_Variations  UniProt_Experimental_Info   GO_Biological_Process   GO_Cellular_Component   GO_Molecular_Function   COSMIC_overlapping_mutations    COSMIC_fusion_genes COSMIC_tissue_types_affected    COSMIC_total_alterations_in_gene    Tumorscape_Amplification_Peaks  Tumorscape_Deletion_Peaks   TCGAscape_Amplification_Peaks   TCGAscape_Deletion_Peaks    DrugBank    ref_context gc_content  CCLE_ONCOMAP_overlapping_mutations  CCLE_ONCOMAP_total_mutations_in_gene    CGC_Mutation_Type   CGC_Translocation_Partner   CGC_Tumor_Types_Somatic CGC_Tumor_Types_Germline    CGC_Other_Diseases  DNARepairGenes_Role FamilialCancerDatabase_Syndromes    MUTSIG_Published_Results    OREGANNO_ID OREGANNO_Values t_alt_count t_ref_count validation_alt_allele   validation_method   validation_status   validation_tumor_sample pox qox pox_cutoff  isArtifactMode  oxoGCut
    MAF = collections.namedtuple('MAF', fields.split())




if __name__ == '__main__':
    generate_maf()
