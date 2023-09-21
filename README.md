### Nextflow wf-human-variation and rapidCNS2 combined workflow

It takes an aligned bam of ONT data along with the .bai index and runs epi2me-labs wf-human-variation to generate SNP, SV, and CNV variants and also aggregate methylation data (the --snp, --sv, --cnv, and --methyl options).

Outputs from wf-human-variation are then analysed using Rapid-CNS2 to generate the final report.

Input BAM file **MUST** contain methylation data (MM:Z tags). The workflow should ERROR if MM:Z tags are not present.
The .bai index must be in the same location as the input bam.

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

### Clone workflow repository
```
git clone https://github.com/graemefox/hv_rapidCNS2.git
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

nextflow run hv_rapidCNS2/main.nf \
-with-docker graefox/rapid_cns2:latest \
-with-report ${OUTPUT_DIR}/${SAMPLE}_nextflow_report.html \
--sample ${SAMPLE} \
--patient ${PATIENT} \
--bam ${BAM} \
--outdir ${OUTPUT_DIR} \
--reference ${REFERENCE}

```


### optional extra parameters (with their default values)
```
--threads 16 (CPUs to use [default: 64]) 
--bam_min_coverage (minimum coverage required to run the epi2melabs/wf-human-variation stages [ default: 5]) 
--minimum_mgmt_cov (minimum avg coverage at the mgmt promoter. Coverage must be greater than this to run the analysis of mgmt methylation)

```

### Troubleshooting tips
If the run seems to hang forever at the cnvpytor step, it may be that you have not indexed your input bam.

Docker Error: "docker: permission denied while trying to connect to the docker daemon socket"
On *buntu systems, you need to add your user to the docker group. 
Follow the instruction here: (https://www.digitalocean.com/community/questions/how-to-fix-docker-got-permission-denied-while-trying-to-connect-to-the-docker-daemon-socket)



