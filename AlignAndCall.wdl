version 1.0

import "AlignmentPipeline.wdl" as AlignAndMarkDuplicates

workflow AlignAndCall {
  meta {
    description: "Takes in unmapped bam and outputs VCF of SNP/Indel calls on the mitochondria."
  }

  input {
    File unmapped_bam
    Float? autosomal_coverage
    String base_name

    String picard
    String gatk
    String haplocheckCLI

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

    Boolean compress_output_vcf

    Float? verifyBamID
    Int? max_low_het_sites
    String? m2_extra_args
    String? m2_filter_extra_args
    Float? vaf_filter_threshold
    Float? f_score_beta

    # Read length used for optimization only. If this is too small CollectWgsMetrics might fail, but the results are not
    # affected by this number. Default is 151.
    Int? max_read_length


  }

  parameter_meta {
    unmapped_bam: "Unmapped and subset bam, optionally with original alignment (OA) tag"
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToMt {
    input:
      input_bam = unmapped_bam,
      picard = picard,
      mt_dict = mt_dict,
      mt_fasta = mt_fasta,
      mt_fasta_index = mt_fasta_index,
      mt_amb = mt_amb,
      mt_ann = mt_ann,
      mt_bwt = mt_bwt,
      mt_pac = mt_pac,
      mt_sa = mt_sa
  }

  call AlignAndMarkDuplicates.AlignmentPipeline as AlignToShiftedMt {
    input:
      input_bam = unmapped_bam,
      picard = picard,
      mt_dict = mt_shifted_dict,
      mt_fasta = mt_shifted_fasta,
      mt_fasta_index = mt_shifted_fasta_index,
      mt_amb = mt_shifted_amb,
      mt_ann = mt_shifted_ann,
      mt_bwt = mt_shifted_bwt,
      mt_pac = mt_shifted_pac,
      mt_sa = mt_shifted_sa
  }

  call CollectWgsMetrics {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bam_index = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      picard = picard,
      read_length = max_read_length,
      coverage_cap = 100000
  }

  #Int? M2_mem = if CollectWgsMetrics.mean_coverage > 25000 then 14 else 7

  call M2 as CallMt {
    input:
      input_bam = AlignToMt.mt_aligned_bam,
      input_bai = AlignToMt.mt_aligned_bai,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      gatk = gatk,
      compress = compress_output_vcf,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:576-16024 "
  }

  call M2 as CallShiftedMt {
    input:
      input_bam = AlignToShiftedMt.mt_aligned_bam,
      input_bai = AlignToShiftedMt.mt_aligned_bai,
      ref_fasta = mt_shifted_fasta,
      ref_fai = mt_shifted_fasta_index,
      ref_dict = mt_shifted_dict,
      gatk = gatk,
      compress = compress_output_vcf,
      # Everything is called except the control region.
      m2_extra_args = select_first([m2_extra_args, ""]) + " -L chrM:8025-9144 "
  }

  call LiftoverAndCombineVcfs {
    input:
      picard = picard,
      shifted_vcf = CallShiftedMt.raw_vcf,
      vcf = CallMt.raw_vcf,
      ref_fasta = mt_fasta,
      ref_fasta_index = mt_fasta_index,
      ref_dict = mt_dict,
      shift_back_chain = shift_back_chain
  }

  call MergeStats {
    input:
      gatk = gatk,
      shifted_stats = CallShiftedMt.stats,
      non_shifted_stats = CallMt.stats
  }

  call Filter as InitialFilter {
    input:
      gatk = gatk,
      raw_vcf = LiftoverAndCombineVcfs.merged_vcf,
      raw_vcf_index = LiftoverAndCombineVcfs.merged_vcf_index,
      raw_vcf_stats = MergeStats.stats,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = 0,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta,
      run_contamination = false
  }

 
  call SplitMultiAllelicsAndRemoveNonPassSites {
    input:
      gatk = gatk,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      filtered_vcf = InitialFilter.filtered_vcf,
      filtered_vcf_idx = InitialFilter.filtered_vcf_idx
  }

  call GetContamination {
    input:
      haplocheckCLI = haplocheckCLI,
      input_vcf = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
  }

  call Filter as FilterContamination {
    input:
      gatk = gatk,
      raw_vcf = InitialFilter.filtered_vcf,
      raw_vcf_index = InitialFilter.filtered_vcf_idx,
      raw_vcf_stats = MergeStats.stats,
      run_contamination = true,
      hasContamination = GetContamination.hasContamination,
      contamination_major = GetContamination.major_level,
      contamination_minor = GetContamination.minor_level,
      verifyBamID = verifyBamID,
      base_name = base_name,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      compress = compress_output_vcf,
      m2_extra_filtering_args = m2_filter_extra_args,
      max_alt_allele_count = 4,
      vaf_filter_threshold = vaf_filter_threshold,
      blacklisted_sites = blacklisted_sites,
      blacklisted_sites_index = blacklisted_sites_index,
      f_score_beta = f_score_beta
 }


  if ( defined(autosomal_coverage) ) {
    call FilterNuMTs {
      input:
        gatk = gatk,
        filtered_vcf = FilterContamination.filtered_vcf,
        filtered_vcf_index = FilterContamination.filtered_vcf_idx,
        ref_fasta = mt_fasta,
        ref_fai = mt_fasta_index,
        ref_dict = mt_dict,
        autosomal_coverage = autosomal_coverage,
        compress = compress_output_vcf
    }
  }

  File low_het_vcf = select_first([FilterNuMTs.numt_filtered_vcf, FilterContamination.filtered_vcf])
  File low_het_vcf_index = select_first([FilterNuMTs.numt_filtered_vcf_idx, FilterContamination.filtered_vcf_idx])

  call FilterLowHetSites {
    input:
      gatk = gatk,
      filtered_vcf = low_het_vcf,
      filtered_vcf_index = low_het_vcf_index,
      ref_fasta = mt_fasta,
      ref_fai = mt_fasta_index,
      ref_dict = mt_dict,
      max_low_het_sites = max_low_het_sites,
      compress = compress_output_vcf,
      base_name = base_name
  }


  output {
    File mt_aligned_bam = AlignToMt.mt_aligned_bam
    File mt_aligned_bai = AlignToMt.mt_aligned_bai
    File mt_aligned_shifted_bam = AlignToShiftedMt.mt_aligned_bam
    File mt_aligned_shifted_bai = AlignToShiftedMt.mt_aligned_bai
    File duplicate_metrics = AlignToMt.duplicate_metrics
    File coverage_metrics = CollectWgsMetrics.metrics
    File theoretical_sensitivity_metrics = CollectWgsMetrics.theoretical_sensitivity
    Int mean_coverage = CollectWgsMetrics.mean_coverage

    File contamination_metrics = GetContamination.contamination_file
    String major_haplogroup = GetContamination.major_hg
    Float contamination = FilterContamination.contamination

    File input_vcf_for_haplochecker = SplitMultiAllelicsAndRemoveNonPassSites.vcf_for_haplochecker
    File out_vcf = FilterLowHetSites.final_filtered_vcf
    File out_vcf_index = FilterLowHetSites.final_filtered_vcf_idx

  }
}




task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String picard
    Int? read_length
    Int? coverage_cap
  }

  Int read_length_for_optimization = select_first([read_length, 151])

  meta {
    description: "Collect coverage metrics"
  }
  parameter_meta {
    read_length: "Read length used for optimization only. If this is too small CollectWgsMetrics might fail. Default is 151."
  }

  command <<<
    set -e

    java -Xms2000m -jar ~{picard} \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=metrics.txt \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length_for_optimization} \
      ~{"COVERAGE_CAP=" + coverage_cap} \
      INCLUDE_BQ_HISTOGRAM=true \
      THEORETICAL_SENSITIVITY_OUTPUT=theoretical_sensitivity.txt

    R --vanilla <<CODE
      df = read.table("metrics.txt",skip=6,header=TRUE,stringsAsFactors=FALSE,sep='\t',nrows=1)
      write.table(floor(df[,"MEAN_COVERAGE"]), "mean_coverage.txt", quote=F, col.names=F, row.names=F)
    CODE
  >>>

  output {
    File metrics = "metrics.txt"
    File theoretical_sensitivity = "theoretical_sensitivity.txt"
    Int mean_coverage = read_int("mean_coverage.txt")
  }
}


task M2 {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    File input_bam
    File input_bai
    String gatk
    Int? max_reads_per_alignment_start
    String? m2_extra_args
    Boolean? make_bamout
    Boolean compress
  }

  Int max_reads_per_alignment_start_arg = select_first([max_reads_per_alignment_start, 75])
  String output_vcf = "raw" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  meta {
    description: "Mutect2 for calling Snps and Indels"
  }
  parameter_meta {
    input_bam: "Aligned Bam"
  }
  command <<<
      set -e

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      java -Xmx4G -jar ~{gatk} Mutect2 \
        -R ~{ref_fasta} \
        -I ~{input_bam} \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O ~{output_vcf} \
        ~{true='--bam-output bamout.bam' false='' make_bamout} \
        ~{m2_extra_args} \
        --annotation StrandBiasBySample \
        --mitochondria-mode \
        --max-reads-per-alignment-start ~{max_reads_per_alignment_start_arg} \
        --max-mnp-distance 0
  >>>

  output {
      File raw_vcf = "~{output_vcf}"
      File raw_vcf_idx = "~{output_vcf_index}"
      File stats = "~{output_vcf}.stats"
      File output_bamOut = "bamout.bam"
  }
}

task LiftoverAndCombineVcfs {
  input {
    String picard
    File shifted_vcf
    File vcf
    String basename = basename(shifted_vcf, ".vcf")

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File shift_back_chain

  }


  meta {
    description: "Lifts over shifted vcf of control region and combines it with the rest of the chrM calls."
  }
  parameter_meta {
    shifted_vcf: "VCF of control region on shifted reference"
    vcf: "VCF of the rest of chrM on original reference"
    ref_fasta: "Original (not shifted) chrM reference"
    shift_back_chain: "Chain file to lift over from shifted reference to original chrM"
  }
  command<<<
    set -e

    java -Xms2000m -jar ~{picard} LiftoverVcf \
      I=~{shifted_vcf} \
      O=~{basename}.shifted_back.vcf \
      R=~{ref_fasta} \
      CHAIN=~{shift_back_chain} \
      REJECT=~{basename}.rejected.vcf

    java -Xms2000m -jar ~{picard} MergeVcfs \
      I=~{basename}.shifted_back.vcf \
      I=~{vcf} \
      O=~{basename}.merged.vcf
    >>>
    output{
        # rejected_vcf should always be empty
        File rejected_vcf = "~{basename}.rejected.vcf"
        File merged_vcf = "~{basename}.merged.vcf"
        File merged_vcf_index = "~{basename}.merged.vcf.idx"
    }
}

task MergeStats {
  input {
    String gatk
    File shifted_stats
    File non_shifted_stats
  }

  command{
    set -e

    java -Xmx4G -jar ~{gatk} MergeMutectStats --stats ~{shifted_stats} --stats ~{non_shifted_stats} -O raw.combined.stats
  }
  output {
    File stats = "raw.combined.stats"
  }

}


task Filter {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File raw_vcf
    File raw_vcf_index
    File raw_vcf_stats
    Boolean compress
    Float? vaf_cutoff
    String base_name

    String? m2_extra_filtering_args
    Int max_alt_allele_count
    Float? autosomal_coverage
    Float? vaf_filter_threshold
    Float? f_score_beta

    Boolean run_contamination
    String? hasContamination
    Float? contamination_major
    Float? contamination_minor
    Float? verifyBamID
     
    File blacklisted_sites
    File blacklisted_sites_index

  }

  String output_vcf = base_name + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"

  Float hc_contamination = if run_contamination && hasContamination == "YES" then (if contamination_major == 0.0 then contamination_minor else 1.0 - contamination_major) else 0.0
  Float max_contamination = if defined(verifyBamID) && verifyBamID > hc_contamination then verifyBamID else hc_contamination

  meta {
    description: "Mutect2 Filtering for calling Snps and Indels"
  }
  parameter_meta {
      vaf_filter_threshold: "Hard cutoff for minimum allele fraction. All sites with VAF less than this cutoff will be filtered."
      f_score_beta: "F-Score beta balances the filtering strategy between recall and precision. The relative weight of recall to precision."
  }
  command <<<
      set -e

      # We need to create these files regardless, even if they stay empty
      touch bamout.bam

      java -Xmx4G -jar ~{gatk} FilterMutectCalls -V ~{raw_vcf} \
        -R ~{ref_fasta} \
        -O filtered.vcf \
        --stats ~{raw_vcf_stats} \
        ~{m2_extra_filtering_args} \
        --max-alt-allele-count ~{max_alt_allele_count} \
        --mitochondria-mode \
        ~{"--min-allele-fraction " + vaf_filter_threshold} \
        ~{"--f-score-beta " + f_score_beta} \
        ~{"--contamination-estimate " + max_contamination}

      java -Xmx4G -jar ~{gatk} VariantFiltration -V filtered.vcf \
        -O ~{output_vcf} \
        --apply-allele-specific-filters \
        --mask ~{blacklisted_sites} \
        --mask-name "blacklisted_site"

  >>>

  output {
      File filtered_vcf = "~{output_vcf}"
      File filtered_vcf_idx = "~{output_vcf_index}"
      Float contamination = "~{hc_contamination}"
  }
}

task SplitMultiAllelicsAndRemoveNonPassSites {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File filtered_vcf_idx
  }

  String basename = basename(filtered_vcf, ".vcf.gz")
  String output_vcf = basename + ".splitAndPassOnly.vcf"

  command {
    set -e

    java -Xmx4G -jar ~{gatk} LeftAlignAndTrimVariants \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O split.vcf \
      --split-multi-allelics \
      --dont-trim-alleles \
      --keep-original-ac

      java -Xmx4G -jar ~{gatk} SelectVariants \
        -V split.vcf \
        -O ~{output_vcf} \
        --exclude-filtered
  
  }
  output {
    #File vcf_for_haplochecker = "splitAndPassOnly.vcf"
    File vcf_for_haplochecker = "~{output_vcf}"
  }
}


task GetContamination {
  input {
    String haplocheckCLI
    File input_vcf
  }

  meta {
    description: "Uses new Haplochecker to estimate levels of contamination in mitochondria"
  }
  parameter_meta {
    input_vcf: "Filtered and split multi-allelic sites VCF for mitochondria"
  }

  String basename = basename(input_vcf, ".splitAndPassOnly.vcf")
  String output_file = basename + ".haplocheck_contamination.txt"

  command <<<
  set -e
  PARENT_DIR="$(dirname "~{input_vcf}")"
  java -Xmx4G -jar ~{haplocheckCLI} "${PARENT_DIR}"

  #sed 's/\"//g' output > output-noquotes
  #grep "SampleID" output-noquotes > headers
  
  sed 's/\"//g' output > ~{output_file}
  grep "SampleID" ~{output_file} > headers


  FORMAT_ERROR="Bad contamination file format"
  if [ `awk '{print $2}' headers` != "Contamination" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $6}' headers` != "HgMajor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $8}' headers` != "HgMinor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $14}' headers` != "MeanHetLevelMajor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi
  if [ `awk '{print $15}' headers` != "MeanHetLevelMinor" ]; then
    echo $FORMAT_ERROR; exit 1
  fi

  #grep -v "SampleID" output-noquotes > output-data
  
  grep -v "SampleID" ~{output_file} > output-data
  awk -F "\t" '{print $2}' output-data > contamination.txt
  awk -F "\t" '{print $6}' output-data > major_hg.txt
  awk -F "\t" '{print $8}' output-data > minor_hg.txt
  awk -F "\t" '{print $14}' output-data > mean_het_major.txt
  awk -F "\t" '{print $15}' output-data > mean_het_minor.txt
  >>>

  output {
    #File contamination_file = "output-noquotes"
    File contamination_file = "~{output_file}"
    String hasContamination = read_string("contamination.txt") 
    String major_hg = read_string("major_hg.txt")
    String minor_hg = read_string("minor_hg.txt")
    Float major_level = read_float("mean_het_major.txt")
    Float minor_level = read_float("mean_het_minor.txt")
  }
}

task FilterNuMTs {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File filtered_vcf_index
    Float? autosomal_coverage
    Boolean compress
  }
  
  String basename = basename(filtered_vcf, ".vcf")
  String output_vcf = basename + ".numt" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  
  parameter_meta {
    autosomal_coverage: "Median coverage of the autosomes for filtering potential polymorphic NuMT variants"
  }

  command {
    set -e

    java -Xmx4G -jar ~{gatk} NuMTFilterTool \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O ~{output_vcf} \
      --autosomal-coverage ~{autosomal_coverage}
  
  }
  output {
    File numt_filtered_vcf = "~{output_vcf}"
    File numt_filtered_vcf_idx = "~{output_vcf_index}"
  }
}

task FilterLowHetSites {
  input {
    String gatk
    File ref_fasta
    File ref_fai
    File ref_dict
    File filtered_vcf
    File filtered_vcf_index
    String base_name
    Int? max_low_het_sites
    Boolean compress
  }

  String output_vcf = base_name + ".final" + if compress then ".vcf.gz" else ".vcf"
  String output_vcf_index = output_vcf + if compress then ".tbi" else ".idx"
  Int max_sites = select_first([max_low_het_sites, 3])
  
  command {
    set -e

    java -Xmx4G -jar ~{gatk} MTLowHeteroplasmyFilterTool \
      -R ~{ref_fasta} \
      -V ~{filtered_vcf} \
      -O ~{output_vcf} \
      --max-allowed-low-hets ~{max_sites}
  
  }
  output {
    File final_filtered_vcf = "~{output_vcf}"
    File final_filtered_vcf_idx = "~{output_vcf_index}"
  }
}





