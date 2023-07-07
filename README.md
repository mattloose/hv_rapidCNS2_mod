### Nextflow wf-human-variation and rapidCNS2 combined workflow

It takes an aligned bam of ONT data, runs epi2me-labs wf-human-variation to generate SNP, SV, and CNV variants and also aggregate methylation data (the --snp, --sv, --cnv, and --methyl options).

Outputs from wf-human-variation are then analysed using Rapid-CNS2 to generate the final report.

Input BAM file **MUST** contain methylation data (MM:Z tags) otherwise methylation related outputs of Rapid-CNS2 will be incorrect/gibberish. The workflow should ERROR if MM:Z tags are not present.

The provided reference should be the same reference used to align reads in the input BAM.

At present, it relies upon the *outdir* being within the current working directory. Setting a full path to some other location will cause errors.

### Requirements:
```
docker
nextflow
```

### Pull latest docker image:
```
docker pull graefox/rapid_cns2:latest
```

### Pull latest workflow
```
nextflow pull graemefox/hv_rapidCNS2
```


### Example command:
```
## define sample name, ID and output directory, input BAM and reference genome:

SAMPLE=sample_01
PATIENT=JohnDoe
OUTPUT_DIR=${SAMPLE}_output
BAM=my_data.bam
REFERENCE=my_reference.fa.gz

## run the pipeline

nextflow run graemefox/hv_rapidCNS2 \
-with-docker graefox/rapid_cns2:latest \
-with-report ${OUTPUT_DIR}/${SAMPLE}_nextflow_report.html \
--sample ${SAMPLE} \
--patient ${PATIENT} \
--bam ${BAM} \
--outdir ${OUTPUT_DIR} \
--reference ${REFERENCE} \
--threads 64

```