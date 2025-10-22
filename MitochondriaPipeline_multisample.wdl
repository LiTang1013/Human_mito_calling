version 1.0

import "AlignAndCall.wdl" as AlignAndCall

workflow MitochondriaMultiSamplePipeline {

  meta {
    description: "Takes in an hg38 bam or cram and outputs VCF of SNP/Indel calls on the mitochondria."
    allowNestedInputs: true
  }

  input {
    #File wgs_aligned_input_bam_or_cram
    #File wgs_aligned_input_bam_or_cram_index

    File inputSamplesFile

    String contig_name = "chrM"
    Float autosomal_coverage = 30

    String gatk
    String picard
    String haplocheckCLI
    Boolean compress_output_vcf

    # Full reference is only requred if starting with a CRAM (BAM doesn't need these files)
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa

    #Shifted reference is used for calling the control region (edge of mitochondria reference).
    #This solves the problem that BWA doesn't support alignment to circular contigs.
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa

    File shift_back_chain
    File blacklisted_sites
    File blacklisted_sites_index

    File non_control_region_interval_list    
    File control_region_shifted_reference_interval_list

  }

  parameter_meta {
    wgs_aligned_input_bam_or_cram: "Full WGS hg38 bam or cram"
    autosomal_coverage: "Median coverage of full input bam"
    out_vcf: "Final VCF of mitochondrial SNPs and INDELs"
    vaf_filter_threshold: "Hard threshold for filtering low VAF sites"
    f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
    contig_name: "Name of mitochondria contig in reference that wgs_aligned_input_bam_or_cram is aligned to"
  }


  Array[Array[String]] inputSamples = read_tsv(inputSamplesFile)

  scatter (sample in inputSamples) {
    call SubsetBamToChrM {
      input:
        gatk = gatk,
        input_bam = sample[0],
        input_bai = sample[1],
        contig_name = contig_name,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict
    }

    call RevertSam {
      input:
        picard = picard,
        input_bam = SubsetBamToChrM.output_bam,
    }

    String base_name = basename(SubsetBamToChrM.output_bam, ".bam")

    call AlignAndCall.AlignAndCall as AlignAndCall {
      input:
        unmapped_bam = RevertSam.unmapped_bam,
        autosomal_coverage = autosomal_coverage,
        base_name = base_name,
        picard = picard,
        gatk = gatk,
        haplocheckCLI = haplocheckCLI,
        mt_dict = mt_dict,
        mt_fasta = mt_fasta,
        mt_fasta_index = mt_fasta_index,
        mt_amb = mt_amb,
        mt_ann = mt_ann,
        mt_bwt = mt_bwt,
        mt_pac = mt_pac,
        mt_sa = mt_sa,
        mt_shifted_dict = mt_shifted_dict,
        mt_shifted_fasta = mt_shifted_fasta,
        mt_shifted_fasta_index = mt_shifted_fasta_index,
        mt_shifted_amb = mt_shifted_amb,
        mt_shifted_ann = mt_shifted_ann,
        mt_shifted_bwt = mt_shifted_bwt,
        mt_shifted_pac = mt_shifted_pac,
        mt_shifted_sa = mt_shifted_sa,
        compress_output_vcf = compress_output_vcf,
        shift_back_chain = shift_back_chain,
        blacklisted_sites = blacklisted_sites,
        blacklisted_sites_index = blacklisted_sites_index
    }

    # This is a temporary task to handle "joint calling" until Mutect2 can produce a GVCF.
    # This proivdes coverage at each base so low coverage sites can be considered ./. rather than 0/0.
    call CoverageAtEveryBase {
      input:
        picard = picard,
        input_bam_regular_ref = AlignAndCall.mt_aligned_bam,
        input_bam_regular_ref_index = AlignAndCall.mt_aligned_bai,
        input_bam_shifted_ref = AlignAndCall.mt_aligned_shifted_bam,
        input_bam_shifted_ref_index = AlignAndCall.mt_aligned_shifted_bai,
        shift_back_chain = shift_back_chain,
        control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list,
        non_control_region_interval_list = non_control_region_interval_list,
        ref_fasta = mt_fasta,
        ref_fasta_index = mt_fasta_index,
        ref_dict = mt_dict,
        shifted_ref_fasta = mt_shifted_fasta,
        shifted_ref_fasta_index = mt_shifted_fasta_index,
        shifted_ref_dict = mt_shifted_dict
    }
    
    call SplitMultiAllelicSites {
      input:
        gatk = gatk,
        input_vcf = AlignAndCall.out_vcf,
        input_vcf_index = AlignAndCall.out_vcf_index,
        base_name = base_name,
        ref_fasta = mt_fasta,
        ref_fasta_index = mt_fasta_index,
        ref_dict = mt_dict
    }

  }

  output {
    #Array[File] subset_bam = SubsetBamToChrM.output_bam
    #Array[File] subset_bam = SubsetBamToChrM.output_bam
    #Array[File] subset_bai = SubsetBamToChrM.output_bai
    #Array[File] mt_aligned_bam = AlignAndCall.mt_aligned_bam
    #Array[File] mt_aligned_bai = AlignAndCall.mt_aligned_bai
    #Array[File] duplicate_metrics = AlignAndCall.duplicate_metrics    
    #Array[File] coverage_metrics = AlignAndCall.coverage_metrics
    #Array[File] theoretical_sensitivity_metrics = AlignAndCall.theoretical_sensitivity_metrics
    #Array[Int] mean_coverage = AlignAndCall.mean_coverage

    Array[File] contamination_metrics = AlignAndCall.contamination_metrics
    #Array[Float] contamination = AlignAndCall.contamination
    #Array[String] major_haplogroup = AlignAndCall.major_haplogroup

    #Array[File] input_vcf_for_haplochecker = AlignAndCall.input_vcf_for_haplochecker
    #Array[File] out_vcf = AlignAndCall.out_vcf
    #Array[File] out_vcf_index = AlignAndCall.out_vcf_index

    Array[File] base_level_coverage_metrics = CoverageAtEveryBase.table
    Array[File] split_vcf = SplitMultiAllelicSites.split_vcf
    Array[File] split_vcf_index = SplitMultiAllelicSites.split_vcf_index
  }

}

task SubsetBamToChrM {
  input {
    String gatk
    File input_bam
    File input_bai
    String contig_name
    String basename = basename(basename(input_bam, ".cram"), ".bam")
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict

  }

  meta {
    description: "Subsets a whole genome bam to just Mitochondria reads"
  }
  parameter_meta {
    ref_fasta: "Reference is only required for cram input. If it is provided ref_fasta_index and ref_dict are also required."
    input_bam: {
      localization_optional: true
    }
    input_bai: {
      localization_optional: true
    }
  }
  command <<<
    java -Xmx4G -jar ~{gatk} PrintReads \
      ~{"-R " + ref_fasta} \
      -L ~{contig_name} \
      --read-filter MateOnSameContigOrNoMappedMateReadFilter \
      --read-filter MateUnmappedAndUnmappedReadFilter \
      -I ~{input_bam} \
      --read-index ~{input_bai} \
      -O ~{basename}.bam
  >>>

  output {
    File output_bam = "~{basename}.bam"
    File output_bai = "~{basename}.bai"
  }
}

task RevertSam {
  input {
    String picard
    File input_bam
    String basename = basename(input_bam, ".bam")
  }

  meta {
    description: "Removes alignment information while retaining recalibrated base qualities and original alignment tags"
  }
  parameter_meta {
    input_bam: "aligned bam"
  }
  command {
    java -Xmx4G -jar ~{picard} \
    RevertSam \
    INPUT=~{input_bam} \
    OUTPUT_BY_READGROUP=false \
    OUTPUT=~{basename}.bam \
    VALIDATION_STRINGENCY=LENIENT \
    ATTRIBUTE_TO_CLEAR=FT \
    ATTRIBUTE_TO_CLEAR=CO \
    SORT_ORDER=queryname \
    RESTORE_ORIGINAL_QUALITIES=false
  }

  output {
    File unmapped_bam = "~{basename}.bam"
  }
}


task CoverageAtEveryBase {
  input {
    String picard
    File input_bam_regular_ref
    File input_bam_regular_ref_index
    File input_bam_shifted_ref
    File input_bam_shifted_ref_index
    File shift_back_chain
    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File shifted_ref_fasta
    File shifted_ref_fasta_index
    File shifted_ref_dict

  }

  meta {
    description: "Remove this hack once there's a GVCF solution."
  }

  String basename = basename(input_bam_regular_ref, ".realigned.bam")
  String coverage_file = basename + ".per_base_coverage.tsv"

  command <<<
    set -e

    java -Xmx4G -jar ~{picard} CollectHsMetrics \
      I=~{input_bam_regular_ref} \
      R=~{ref_fasta} \
      PER_BASE_COVERAGE=non_control_region.tsv \
      O=non_control_region.metrics \
      TI=~{non_control_region_interval_list} \
      BI=~{non_control_region_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    java -Xmx4G -jar ~{picard} CollectHsMetrics \
      I=~{input_bam_shifted_ref} \
      R=~{shifted_ref_fasta} \
      PER_BASE_COVERAGE=control_region_shifted.tsv \
      O=control_region_shifted.metrics \
      TI=~{control_region_shifted_reference_interval_list} \
      BI=~{control_region_shifted_reference_interval_list} \
      COVMAX=20000 \
      SAMPLE_SIZE=1

    R --vanilla <<CODE
      shift_back = function(x) {
        if (x < 8570) {
          return(x + 8000)
        } else {
          return (x - 8569)
        }
      }

      control_region_shifted = read.table("control_region_shifted.tsv", header=T)
      shifted_back = sapply(control_region_shifted[,"pos"], shift_back)
      control_region_shifted[,"pos"] = shifted_back

      beginning = subset(control_region_shifted, control_region_shifted[,'pos']<8000)
      end = subset(control_region_shifted, control_region_shifted[,'pos']>8000)

      non_control_region = read.table("non_control_region.tsv", header=T)
      combined_table = rbind(beginning, non_control_region, end)
      #write.table(combined_table, "per_base_coverage.tsv", row.names=F, col.names=T, quote=F, sep="\t")
      write.table(combined_table, "~{coverage_file}", row.names=F, col.names=T, quote=F, sep="\t")

    CODE
  >>>

  output {
    #File table = "per_base_coverage.tsv"
    File table = "~{coverage_file}"

  }
}

task SplitMultiAllelicSites {
  input {
    String gatk
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File input_vcf
    File input_vcf_index
    String base_name
  }

  String output_vcf = base_name + ".final.split.vcf"
  String output_vcf_index = output_vcf + ".idx"

  command {
    set -e

    java -Xmx4G -jar ~{gatk} LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{input_vcf} \
      -O ~{output_vcf} \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac
  }
  output {
    File split_vcf = "~{output_vcf}"
    File split_vcf_index = "~{output_vcf}"
  }

}


