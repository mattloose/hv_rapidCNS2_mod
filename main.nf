#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// set up and create an output directory
outdir = file(params.outdir)
outdir.mkdir()

log.info """\

        ================================================================
        Human-variation -> Rapid_CNS2 - Nextflow P I P E L I N E
        ================================================================

        INPUTS
        ================================================================
        sample                  : ${params.sample}
        patient                 : ${params.patient}
        input_bam               : ${params.bam}
        outdir                  : ${params.outdir}
        reference               : ${params.reference}
        threads                 : ${params.threads}
        bam_min_coverage        : ${params.bam_min_coverage}
        min_mgmt_coverage	: ${params.minimum_mgmt_cov}
        ================================================================
        To run with SLURM, uncomment the "process.executor" line in nextflow.config
        ================================================================

        """
        .stripIndent()

process check_bam_has_meth_data {
    input:
        path(input_bam)
    
    output:
        stdout emit: meth_check

    script:
        """
        samtools view ${input_bam} | grep -m 1 MM:Z
        """
}

process index_input_bam {
    input:
        path(input_bam)

    output:
        path "*.bai", emit: indexed_bam

    script:
        """
        samtools index ${input_bam}
        """
}

process check_mgmt_coverage {
    input:
        path(input_bam)
	path(mgmt_bed)
	val(minimum_mgmt_coverage)

    publishDir("${params.outdir}")

    output:
	val true
	path "*.txt", emit: mgmt_avg_cov_file
        stdout emit: mgmt_avg_cov

    script:
        """
        cov="\$(bedtools coverage -a ${mgmt_bed} -b ${input_bam} | awk '{print \$4}')"
        echo \${cov}
	if [ \${cov} -ge ${minimum_mgmt_coverage} ]; then
		echo \${cov} > mgmt_avg_cov.txt
        else 
                echo \${cov} > mgmt_cov_below_thresh.txt
	fi
	"""
}

process draw_mgmt_methylartist {
    input:
        path(indexed_bam)
        path(bam)
        path(reference)
        val(params.outdir)
	val ready // mgmt_coverage has run

    publishDir("${params.outdir}")
    
    output:
        val true
        path "*.png", emit: mgmt_plot optional true
    
    script:
        cov_file = file("${params.outdir}/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )    
            """
            methylartist locus -i chr10:129466536-129467536 -b ${bam} --ref ${reference} --motif CG
            """
        else
            """
            """
}

process mosdepth {
    input:
        val(threads)
        path(targets)
        path(input_bam)
        val(sample)
        path(indexed_bam)
    
    publishDir("${params.outdir}")
	
    output:
        path "*.mosdepth.summary.txt", emit: mosdepth_out

    script:
        """
        /mosdepth -t ${threads} -n --by ${targets} --fast-mode ${sample} ${input_bam}
        """
}

process index_subsetted_bam {
    input:
        path(input_bam)
    
    output:
        val true

    script:
        """
        samtools index ${input_bam}
        """
}

process human_variation_sv_methyl_cnv {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        path(indexed_bam)
    
    publishDir("${params.outdir}")

    output:
        val true
    
    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
            -with-report ${PWD}/${params.outdir}/human_variation_sv_methyl_cnv_nextflow_report.html \
            -profile standard \
            -w ${PWD}/${params.outdir}/workspace \
            --ref ${reference} \
            --sv \
            --methyl \
            --cnv \
            --bam ${input_bam} \
            --bed ${targets_bed} \
            --sample_name ${sample} \
            --bam_min_coverage ${bam_min_coverage} \
            --out_dir ${PWD}/${params.outdir}
        """ 
}
 process subset_bam_by_bed {
    input:
        path(input_bam)
        path(input_bed)
        val(sample)

    output:
        path "*subset.bam", emit: subsetted_bam

    script:
        """
        samtools view \
        -b \
        -h \
        -L ${input_bed} ${input_bam} > ${sample}_subset.bam
        """
 }

 process human_variation_snp {
    input:
        path(input_bam)
        path(targets_bed)
        path(reference)
        val(sample)
        val(outdir)
        val(bam_min_coverage)
        val ready

    output:
        val true

    script:
        """
        nextflow run epi2me-labs/wf-human-variation \
            -with-report ${params.outdir}/human_variation_snp_nextflow_report.html \
            -profile standard \
            -w ${PWD}/${params.outdir}/workspace \
            --ref ${reference} \
            --snp \
            --bam ${input_bam} \
            --bed ${targets_bed} \
            --sample_name ${sample} \
            --bam_min_coverage ${bam_min_coverage} \
            --out_dir ${PWD}/${params.outdir}
        """
 }

 process cnvpytor{
    input:
        val(sample)
        path(input_bam)
        val(threads)

    output:
        val true
    
    publishDir("${params.outdir}")

    script:
        """
        cnvpytor -root ${sample}_CNV.pytor -rd ${input_bam.toRealPath()} -j ${threads}
        cnvpytor -root ${sample}_CNV.pytor -his 1000 10000 100000 -j ${threads} # 4 mins
        cnvpytor -root ${sample}_CNV.pytor -partition 1000 10000 100000 -j ${threads} # SLOW
        cnvpytor -root ${sample}_CNV.pytor -call 1000 -j ${threads} > ${sample}.cnvpytor.calls.1000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 10000 -j ${threads} > ${sample}.cnvpytor.calls.10000.tsv
        cnvpytor -root ${sample}_CNV.pytor -call 100000 -j ${threads} > ${sample}.cnvpytor.calls.100000.tsv
        cnvpytor -root ${sample}_CNV.pytor -plot manhattan 100000 -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY -o ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.png
        """
}

process vcftools {
    input:
        val(ready) // human-variation sv, cnv, methyl
        val(sample)

    output:
        path "*.vcf", emit: variants_out

    script:
        """
        vcftools --gzvcf ${PWD}/${params.outdir}/${sample}.wf_snp.vcf.gz \
        --remove-filtered-all \
        --recode \
        --out ${sample}
        """
}

process gzip {
    input:
        path(input_file)

    output:
        path "*.gz", emit: compressed_out

    script:
        """
        gzip -c ${input_file} > ${input_file}.gz
        """
}

process bedtools_intersect {
    input:
        path(input1)
        path(input2)
        val(output_file)
        val(ext)

    output:   
        path "*.vcf", emit: intersect_vcf
        
    script:
        """
        bedtools intersect \
        -a ${input1} \
        -b ${input2} > ${output_file}${ext}
        """
}

process convert2annovar{
    input:
        path(input)
        val(output_file)
        val(ext)

    output:
        path "*.avinput", emit: annovar_input

    script:
        """
        /annovar/convert2annovar.pl \
        -format vcf4 ${input} \
        -withfreq \
        -includeinfo > ${output_file}${ext}
        """
}

process table_annovar {
    input:
        path(input)
        val(annovar_ver)
        val(output_file)
        val(ext)
    
    output:
        path "*_multianno.csv", emit: clair3_output
      
    script:
        """
        /annovar/table_annovar.pl ${input} \
        /annovar/humandb/ \
        -buildver ${annovar_ver} \
        -out ${output_file}${ext} \
        -protocol refGene,cytoBand,avsnp147,dbnsfp30a,1000g2015aug_eur,cosmic70 \
        -operation gx,r,f,f,f,f \
        -nastring . \
        -csvout \
        -polish \
        -otherinfo
        """
}

process bedtools_intersect2 {
    input:
        val(ready) // bedmethyl file from human-variation_sv_methyl_cnv
        path(input2)
        val(output_file)
        val(ext)
        val(sample)

    publishDir("${params.outdir}")

    output:   
        path "*.bed", emit: intersect_bed
        
    script:
        """
        bedtools intersect \
        -a ${PWD}/${params.outdir}/${sample}.methyl.cpg.bed.gz \
        -b ${input2} > ${output_file}.${ext}
        """
}

process mgmt_pred {
    input:
        path(mgmt_pred)
        file(intersect_bed)
        path(probes)
        path(model)
        val(sample)
        val(params.outdir)

    output:
        val true

    script:
        cov_file = file("${params.outdir}/mgmt_avg_cov.txt")
        if ( cov_file.exists() == true )
        """
            Rscript ${mgmt_pred} \
            --input ${intersect_bed} \
            --probes ${probes} \
            --model ${model} \
            --out_dir ${PWD}/${params.outdir} \
            --sample ${sample}
        """
        else
            """
            """
}

process meth_classification {
    input:
        path(meth_class)
        val(sample)
        val(params.outdir)
        path(topprobes)
        path(trainingdata)
        path(arrayfile)
        val(threads)
        val(ready) // human_var_sv_cnv_meth

    publishDir("${params.outdir}")

    output:
        val true

    script:
        """
        Rscript ${meth_class} \
        --sample ${sample} \
        --out_dir ${PWD}/${params.outdir} \
        --in_file ${PWD}/${params.outdir}/${sample}.methyl.cpg.bed.gz \
        --probes ${topprobes} \
        --training_data ${trainingdata} \
        --array_file ${arrayfile} \
        --threads ${threads}
        """
}
process filter_report {
    input:
        path(filterreport)
        path(clair3_multianno)
        val(sample)
        val(params.outdir)
    
    publishDir("${params.outdir}")

    output:
        val true
    
    script:
        """
        Rscript ${filterreport} \
        --input ${clair3_multianno} \
        --output ${PWD}/${params.outdir}/${sample}_clair3_report.csv \
        --sample ${sample}
        """        
}

process make_report {
    input:
        path(makereport)
        val(ready) // cnvpytor
        val(ready) // mgmt_pred
        val(ready) // meth_classification
        val(ready) // filter_report
        val(sample)
        path(report_UKHD)
        path(mosdepth_plot_data) // mosdepth
        val(methylartist_plot)
        path(report_UKHD_no_mgmt)
        val(mgmt_cov)
    
    output:
	val true

    script:
        // if mgmt coverage was above thresh (and we have the data)
        mgmt_status_file = file("${params.outdir}/${sample}_mgmt_status.csv")
        if ( mgmt_status_file.exists() == true )
        
        """
        Rscript ${makereport} \
        --prefix ${sample} \
        --mutations ${PWD}/${params.outdir}/${sample}_clair3_report.csv \
        --cnv_plot ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.global.0000.png \
        --rf_details ${PWD}/${params.outdir}/${sample}_rf_details.tsv \
        --votes ${PWD}/${params.outdir}/${sample}_votes.tsv \
        --output_dir ${PWD}/${params.outdir} \
        --patient ${sample} \
        --coverage ${PWD}/${params.outdir}/${sample}.mosdepth.summary.txt \
        --sample ${sample} \
        --report_UKHD ${report_UKHD} \
        --methylartist ${PWD}/${params.outdir}/*.locus.meth.png \
        --mgmt ${PWD}/${params.outdir}/${sample}_mgmt_status.csv \
        --report_UKHD_no_mgmt ${report_UKHD_no_mgmt} \
        --promoter_mgmt_coverage ${mgmt_cov}
        """

        // else, mgmt was below thresh and there is no data
        else
        """
        Rscript ${makereport} \
        --prefix ${sample} \
        --mutations ${PWD}/${params.outdir}/${sample}_clair3_report.csv \
        --cnv_plot ${PWD}/${params.outdir}/${sample}_cnvpytor_100k.global.0000.png \
        --rf_details ${PWD}/${params.outdir}/${sample}_rf_details.tsv \
        --votes ${PWD}/${params.outdir}/${sample}_votes.tsv \
        --output_dir ${PWD}/${params.outdir} \
        --patient ${sample} \
        --coverage ${PWD}/${params.outdir}/${sample}.mosdepth.summary.txt \
        --sample ${sample} \
        --report_UKHD ${report_UKHD} \
        --report_UKHD_no_mgmt ${report_UKHD_no_mgmt} \
        --promoter_mgmt_coverage ${mgmt_cov} \
        """
}

///////////////////////////
// MAIN WORKFLOW SECTION //
///////////////////////////

workflow {

///////////////////// - read all input parameters

///////////////////// - read user defined parameters
    Channel.fromPath(params.bam, checkIfExists: true)
    .set {input_bam}

    Channel.fromPath(params.reference, checkIfExists: true)
    .set {reference}

    Channel.from(params.sample)
    .set {sample}

    Channel.from(params.outdir)
    .set {outdir}

    Channel.from(params.patient)
    .set {patient}

    Channel.from(params.threads)
    .set {threads}

    Channel.from(params.minimum_mgmt_cov)
    .set {minimum_mgmt_cov}

///////////////////// - read required scripts and files from bin/

    Channel.fromPath("${projectDir}/bin/NPHD_panel_hg38_clean.bed", checkIfExists: true)
    .set {targets}

    Channel.fromPath("${projectDir}/bin/mgmt_hg38.bed", checkIfExists: true)
    .set {mgmt_bed}

    Channel.fromPath("${projectDir}/bin/mgmt_probes.Rdata", checkIfExists: true)
    .set {probes}

    Channel.fromPath("${projectDir}/bin/mgmt_137sites_mean_model.Rdata", checkIfExists: true)
    .set {model}

    Channel.fromPath("${projectDir}/bin/mgmt_pred_v0.2.R", checkIfExists: true)
    .set {mgmt_pred}

    Channel.fromPath("${projectDir}/bin/methylation_classification_nanodx_v0.1.R", checkIfExists: true)
    .set {meth_class}

    Channel.fromPath("${projectDir}/bin/top_probes_hm450.Rdata", checkIfExists: true)
    .set {topprobes}
    
    Channel.fromPath("${projectDir}/bin/capper_top_100k_betas_binarised.Rdata", checkIfExists: true)
    .set {trainingdata}

    Channel.fromPath("${projectDir}/bin/HM450.hg38.manifest.gencode.v22.Rdata", checkIfExists: true)
    .set {arrayfile}

    Channel.fromPath("${projectDir}/bin/filter_report_v0.1.R", checkIfExists: true)
    .set {filterreport}

    Channel.fromPath("${projectDir}/bin/make_report_v0.2.R", checkIfExists: true)
    .set {makereport}

    Channel.fromPath("${projectDir}/bin/Rapid_CNS2_report_UKHD_v0.2.Rmd", checkIfExists: true)
    .set {report_UKHD}

    Channel.fromPath("${projectDir}/bin/Rapid_CNS2_report_UKHD_no_mgmt_v0.2.Rmd", checkIfExists: true)
    .set {report_UKHD_no_mgmt}

///////////////////// - run the workflow
    // check input bam file for methylation tags
    check_ch = check_bam_has_meth_data(input_bam)

    // index the input bam file 
    index_ch = index_input_bam(input_bam)

    // check the coverage of the mgmt region in the sequence data
    mgmt_coverage_ch = check_mgmt_coverage(input_bam, mgmt_bed, minimum_mgmt_cov)

    // methylartist mgmt plot
    // methylartist only runs if the coverage is detected to be high enough (> minimum_mgmt_cov)
    methyl_artist_ch = draw_mgmt_methylartist(index_ch.indexed_bam, input_bam, reference, outdir, check_mgmt_coverage.out[0])

    // mosdepth coverage plots
    mosdepth_ch = mosdepth(threads, targets, input_bam, sample, index_ch.indexed_bam)

    // call and run the epi2me-labs/wf-human-variation
    human_variation_sv_methyl_cnv(input_bam, targets, reference, sample, outdir, 1, index_ch.indexed_bam)

    // subset the input bam file by the bed file - to speed up the SNP portion of human-variation
    subsetted_bam_ch = subset_bam_by_bed(input_bam, targets, sample)

    // index the subsetted bam generated above
    index_subsetted_bam(subsetted_bam_ch.subsetted_bam)

    // run the SNP human variation workflow on the subsetted bam
    human_variation_snp(subsetted_bam_ch.subsetted_bam, targets, reference, sample, outdir, 1, index_subsetted_bam.out)

    // run cnvPYTOR
    cnvpytor(sample, input_bam, threads)

    // run vcftools
    variants_ch = vcftools(human_variation_snp.out, sample)

    // compress output from vcftools
    compressed_variants_ch = gzip(variants_ch.variants_out)

    // run bedtools intersect on output of vcftools
    vcf_intersect_ch = bedtools_intersect(compressed_variants_ch.compressed_out, targets, sample, '_clair3_panel.vcf')

    // convert vcf to annovar input
    converted_annovar_ch = convert2annovar(vcf_intersect_ch.intersect_vcf, sample, '_clair3_panel.avinput')

    // run table_annovar on the converted annovar input
    clair3_annovar_ch = table_annovar(converted_annovar_ch.annovar_input, 'hg38', sample, '_clair3_panel')

    // run bedtools intersect
    intersect_bed_ch = bedtools_intersect2(human_variation_sv_methyl_cnv.out, mgmt_bed, 'mgmt_5mC.hg38', 'bed', sample)

    // run the mgmt_pred script
    mgmt_pred_ch = mgmt_pred(mgmt_pred, intersect_bed_ch.intersect_bed, probes, model, sample, params.outdir)
   
    // run the meth classification script
    meth_classification(meth_class, sample, params.outdir, topprobes, trainingdata, arrayfile, threads, human_variation_sv_methyl_cnv.out)
    
    // collect report data and generate report
    filter_report(filterreport, clair3_annovar_ch.clair3_output, sample, params.outdir)  

    // collect data and generate final report    
    make_report(makereport, cnvpytor.out, mgmt_pred.out, meth_classification.out, filter_report.out, sample, report_UKHD, mosdepth_ch.mosdepth_out, draw_mgmt_methylartist.out[0], report_UKHD_no_mgmt, mgmt_coverage_ch.mgmt_avg_cov)
}
