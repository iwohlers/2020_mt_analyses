# Getting sequencing reads mapped to the mitochondrium from the WGS HGDP data
# set of (Bergström et al., 2020)
# Before executing this workflow, run 
# 'snakemake BERGSTROEM2020/populations_indiviuals.txt' to get the IDs and 
# populations of the samples to be downloaded

rule mt_sam_to_bam:
    input: "{dataset}/{population}/{individual}_mt.sam"
    output: "{dataset}/{population}/{individual}_mt.bam"
    conda: "envs/samtools.yaml"
    shell: "samtools view -S -b {input} > {output}"



################################################################################
########################### Sudanese MT sequencing files #######################
################################################################################

rule get_indv_pop_sudan:
    input: "data/raw_sudanese/SUDAN_populations_individuals.txt"
    output: "SUDAN2020/populations_individuals.txt"
    shell: "cp {input} {output}"

INDIVIDUALS_SUDAN = []
POPULATION_SUDAN = {}
with open("SUDAN2020/populations_individuals.txt","r") as f_in:
    for line in f_in:
        pop,indv = line.strip().split("\t")
        # Three individuals have so little reads that haplogroups cannot be 
        # computed (2, 14, 2 reads)
        if not indv in ["Psor-092_S15","Psor-028_S13","Psor-102_S25"]:
            INDIVIDUALS_SUDAN.append(indv)
            POPULATION_SUDAN[indv] = pop

rule fastqc:
    input: "data/raw_sudanese/{sample}_L001_{mate}_001.fastq.gz"
    output: "fastqc/{sample}_L001_{mate}_001_fastqc.html"
    conda: "envs/fastqc.yaml"
    shell: "fastqc --extract --outdir=fastqc/ {input}"
    
rule fastqc_all:
    input: expand("fastqc/{sample}_L001_{mate}_001_fastqc.html",sample=INDIVIDUALS_SUDAN,mate=["R1","R2"])

rule fastqc_summary:
    input: expand("fastqc/{sample}_L001_{mate}_001_fastqc/summary.txt",sample=INDIVIDUALS_SUDAN,mate=["R1","R2"])
    output: "fastqc/fastqc_summary.txt"
    run:
        with open(output[0],"w") as f_out:
            header = ["","Basic Statistics","Per base sequence quality",\
            "Per tile sequence quality","Per sequence quality scores", \
            "Per base sequence content","Per sequence GC content", \
            "Per base N content","Sequence Length Distribution", \
            "Sequence Duplication Levels","Overrepresented sequences", \
            "Adapter Content"]
            f_out.write("\t".join(header)+"\n")
            for filename in input:
                f_out.write(filename.split("/")[1]+"\t")
                with open(filename,"r") as f_in:
                    i = 0
                    for line in f_in:
                        if i<10:
                            f_out.write(line.split("\t")[0]+"\t")
                            i += 1
                        else:
                            f_out.write(line.split("\t")[0]+"\n")
                           

rule get_mt_reference:
    output: "ref/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"
    shell: "wget -P ref ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"

# Mapping to reference/assembly using bwa
rule bwa_index:
    input: "ref/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz"
    output: "bwa_index/Homo_sapiens.GRCh38.dna.chromosome.MT.sa"
    conda: "envs/bwa.yaml"
    shell: "bwa index -p bwa_index/Homo_sapiens.GRCh38.dna.chromosome.MT {input}"

rule bwa_mem:
    input: index = "bwa_index/Homo_sapiens.GRCh38.dna.chromosome.MT.sa",
           fastq_r1 = "data/raw_sudanese/{sample}_L001_R1_001.fastq.gz",
           fastq_r2 = "data/raw_sudanese/{sample}_L001_R2_001.fastq.gz",
    output: "SUDAN2020/mapped/{sample}.sam"
    conda: "envs/bwa.yaml"
    shell: "bwa mem -t 8 " + \
           "bwa_index/Homo_sapiens.GRCh38.dna.chromosome.MT " + \
           "{input.fastq_r1} {input.fastq_r2} > {output}"

rule samtools_sort_to_bam:
    input: "SUDAN2020/mapped/{sample}.sam"
    output: "SUDAN2020/mapped/{sample}.bam"
    conda: "envs/samtools.yaml"
    shell: "samtools sort -O BAM {input} > {output}"

rule samtools_stats:
    input: "SUDAN2020/mapped/{sample}.bam"
    output: "SUDAN2020/mapped/{sample}.stats"
    conda: "envs/samtools.yaml"
    shell: "samtools stats {input} > {output}"

rule mapping_stats_all:
    input: expand("SUDAN2020/mapped/{sample}.stats",sample=INDIVIDUALS_SUDAN)

rule prepare_sudan:
    input: "SUDAN2020/mapped/{sample}.bam"
    output: "SUDAN2020/Sudanese/{sample}_mt.sam"
    shell: "samtools view -h {input} > {output}"

rule summarize_mapping_stats:
    input: expand("SUDAN2020/mapped/{sample}.stats", sample=INDIVIDUALS_SUDAN)
    output: "SUDAN2020/mapped/summary.stats"
    run:
        with open(output[0],"w") as f_out:
            f_out.write("sample\treads\treads_mapped\tbases_mapped\taverage_length\tapprox_mt_coverage\n")
            for filename in input:
                with open(filename,"r") as f_in:
                    f_out.write(filename.split("/")[-1].split(".")[0]+"\t")
                    for line in f_in:
                        if line[:13] == "SN\tsequences:":
                            f_out.write(line.strip("\n").split("\t")[2]+"\t")
                        if line[:16] == "SN\treads mapped:":
                            f_out.write(line.strip("\n").split("\t")[2]+"\t")     
                        if line[:16] == "SN\tbases mapped:":
                            bases_mapped = line.strip("\n").split("\t")[2]
                            f_out.write(bases_mapped+"\t")
                        if line[:18] == "SN\taverage length:":
                            f_out.write(line.strip("\n").split("\t")[2]+"\t") 
                    f_out.write(str(round(int(bases_mapped)/16569,3))+"\n")

rule get_coverage:
    input: "SUDAN2020/mapped/{sample}.stats"
    output: "SUDAN2020/mapped/{sample}.cov"
    shell: "cat {input} | grep ^COV | cut -f 2- > {output}"

rule get_coverage_all:
    input: expand("SUDAN2020/mapped/{sample}.cov",sample=INDIVIDUALS_SUDAN)

rule prepare_sudan_mt_all:
     input: expand("SUDAN2020/Sudanese/{individual}_mt.sam", \
            individual=INDIVIDUALS_SUDAN)

################################################################################
####################### Bergstroem MT reads from bam files #####################
################################################################################

rule get populations_and_individual_ids:
    input: "BERGSTROEM2020/hgdp.high_coverage.GRCh38DH.alignment.index"
    output: "BERGSTROEM2020/populations_indiviuals.txt"
    run:
        with open(input[0],"r") as f_in, open(output[0],"w") as f_out:
            for line in f_in:
                if line[0] == "#":
                    continue
                s = line.split("/")
                population = s[8]
                individual = s[9]
                f_out.write(population+"\t"+individual+"\n")

INDIVIDUALS_HGDP = []
POPULATION_HGDP = {}
with open("BERGSTROEM2020/populations_indiviuals.txt","r") as f_in:
    for line in f_in:
        pop,indv = line.strip().split("\t")
        INDIVIDUALS_HGDP.append(indv)
        POPULATION_HGDP[indv] = pop
        
rule download_readme_bergstroem:
     output: "BERGSTROEM2020/README_HGDP.md"
     shell: "wget -P BERGSTROEM2020 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/README_HGDP.md"

rule download_reusestatement_bergstroem:
     output: "BERGSTROEM2020/README_HGDP_datareuse_statement.md"
     shell: "wget -P BERGSTROEM2020 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/README_HGDP_datareuse_statement.md"

rule download_alignmentindex_bergstroem:
     output: "BERGSTROEM2020/hgdp.high_coverage.GRCh38DH.alignment.index"
     shell: "wget -P BERGSTROEM2020 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/hgdp.high_coverage.GRCh38DH.alignment.index"

rule download_sequenceindex_bergstroem:
     output: "BERGSTROEM2020/hgdp_wgs.sequence.index"
     shell: "wget -P BERGSTROEM2020 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/hgdp_wgs.sequence.index"

rule download_bergstroem_mt:
     output: "BERGSTROEM2020/{population}/{individual}_mt.sam"
     shell: "samtools view -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGDP/data/{wildcards.population}/{wildcards.individual}/alignment/{wildcards.individual}.alt_bwamem_GRCh38DH.20181023.{wildcards.population}.cram chrM > BERGSTROEM2020/{wildcards.population}/{wildcards.individual}_mt.sam"

rule download_bergstroem_mt_all:
     input: expand("BERGSTROEM2020/{population}/{individual}_mt.sam", zip, \
                    individual=INDIVIDUALS_HGDP, \
                    population=[POPULATION_HGDP[indv] for indv in INDIVIDUALS_HGDP])


################################################################################
####################### 1000G MT reads from bam files ##########################
################################################################################

POP2POPCODE = {
"African Ancestry in Southwest US": "ASW",
"African Caribbean in Barbados": "ACB",
"Bengali in Bangladesh": "BEB",
"British in England and Scotland": "GBR",
"Chinese Dai in Xishuangbanna, China": "CFX",
"Colombian in Medellin, Colombia": "CLM",
"Esan in Nigeria": "ESN",
"Finnish in Finland": "FIN",
"Gambian in Western Division, The Gambia - Mandinka": "GWD",
"Gujarati Indians in Houston, TX": "GIH",
"Han Chinese in Beijing, China": "CHB",
"Han Chinese South": "CHS",
"Iberian populations in Spain": "IBS",
"Indian Telugu in the UK": "ITU",
"Japanese in Tokyo, Japan": "JPT",
"Kinh in Ho Chi Minh City, Vietnam": "KHV",
"Luhya in Webuye, Kenya": "LWK",
"Mende in Sierra Leone": "MSL",
"Mexican Ancestry in Los Angeles, California": "MXL",
"Peruvian in Lima, Peru": "PEL",
"Puerto Rican in Puerto Rico": "PUR",
"Punjabi in Lahore, Pakistan": "PJL",
"Sri Lankan Tamil in the UK": "STU",
"Toscani in Italy": "TSI",
"Utah residents (CEPH) with Northern and Western European ancestry": "CEU",
"Yoruba in Ibadan, Nigeria": "YRI"
}

# The file with the details on where to download the cram files was downloaded
# here: https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
# https://www.internationalgenome.org/api/beta/file/_search/igsr_30x%20GRCh38.tsv.tsv
INDIVIDUALS_1000G = []
POPULATION_1000G = {}
FTP_1000G = {}
with open("data/igsr_30x GRCh38.tsv.tsv","r") as f_in:
    for line in f_in:
        if line[:3] == "url":
            continue
        s = line.strip().split("\t")
        indv = s[5]
        pop = s[6]
        INDIVIDUALS_1000G.append(indv)
        POPULATION_1000G[indv] = POP2POPCODE[pop]
        FTP_1000G[indv] = s[0]


# These files have been remapped against GRCh38, and the list of files has been
# obtained from 
# https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
# These are the cram files of 2504 individuals from 26 populations

rule get_datareuse_1000g:
    output: "1000G/20190405_1000G_2504_high_cov_data_reuse_README.md"
    shell: "wget -P 1000G http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20190405_1000G_2504_high_cov_data_reuse_README.md"

rule download_1000g_mt:
    output: "1000G/{population}/{individual}_mt.sam"
    run: 
        shell("samtools view -h "+FTP_1000G[wildcards.individual]+" chrM > "
        "1000G/"+wildcards.population+"/"+wildcards.individual+"_mt.sam")

rule download_1000g_mt_all:
     input: expand("1000G/{population}/{individual}_mt.sam", zip, \
                    individual=INDIVIDUALS_1000G, \
                    population=[POPULATION_1000G[indv] for indv in INDIVIDUALS_1000G])


################################################################################
############################### Variant calling ################################
################################################################################

# Puttick et al., 2019: Using default FreeBayes settings, the reproducibility from 13 WGS replicates of the NA12878 cell line was poor (Supp. Fig 3a). We thus optimised the mapping quality (MQ≥30), and base quality (BQ≥24) filters, resulting in an average of 17.5±0.5 variants per replicate, with variant allele frequency (VAF) >0.01 (Supp. Fig 3a). 


################################################################################
################# Haplogroups from BAM files usinng Haplocheck #################
################################################################################

# Download and install
#mkdir haplocheck 
#cd haplocheck
#wget https://github.com/genepi/cloudgene/releases/download/v2.2.0/cloudgene-installer.sh
#bash cloudgene-installer.sh
#./cloudgene install https://github.com/genepi/haplocheck/releases/download/v1.1.2/haplocheck.zip 

rule run_haplocheck:
    input: "{dataset}/{population}/{individual}_mt.bam"
    output: "haplocheck/results/{dataset}_{population}_{individual}/haplogroups/haplogroups.txt"
    wildcard_constraints: 
        dataset="[A-Z,0-9]+",
        population="[A-Z,a-z]+"
    params: prefix=lambda wildcards, output: "/".join(output[0].split("/")[1:-2])
    shell: "./haplocheck/cloudgene run haplocheck@1.1.2 " + \
                                            "--files ../{input} " +\
                                            "--output {params.prefix} "

rule run_haplocheck_1000g:
    input: expand("haplocheck/results/1000G_{population}_{individual}/haplogroups/haplogroups.txt", \
                    zip, \
                    individual=INDIVIDUALS_1000G, \
                    population=[POPULATION_1000G[indv] for indv in INDIVIDUALS_1000G])

rule run_haplocheck_bergstroem:
    input: expand("haplocheck/results/BERGSTROEM2020_{population}_{individual}/haplogroups/haplogroups.txt", \
                    zip, \
                    individual=INDIVIDUALS_HGDP, \
                    population=[POPULATION_HGDP[indv] for indv in INDIVIDUALS_HGDP])

rule run_haplocheck_sudan:
    input: expand("haplocheck/results/SUDAN2020_Sudanese_{individual}/haplogroups/haplogroups.txt", \
                    individual=INDIVIDUALS_SUDAN)

rule combined_haplogroup_file:
    input: expand("haplocheck/results/1000G_{population}_{individual}/haplogroups/haplogroups.txt", \
                    zip, \
                    individual=INDIVIDUALS_1000G, \
                    population=[POPULATION_1000G[indv] for indv in INDIVIDUALS_1000G]), \
           expand("haplocheck/results/BERGSTROEM2020_{population}_{individual}/haplogroups/haplogroups.txt", \
                    zip, \
                    individual=INDIVIDUALS_HGDP, \
                    population=[POPULATION_HGDP[indv] for indv in INDIVIDUALS_HGDP]),
           expand("haplocheck/results/SUDAN2020_Sudanese_{individual}/haplogroups/haplogroups.txt", \
                    individual=INDIVIDUALS_SUDAN) 
    output: "haplocheck/results/all_haplogroups.txt"
    shell: "cat {input[0]} | head -n 1 > {output}; " + \
            "cat haplocheck/results/*/haplogroups/haplogroups.txt | " + \
            "grep -v SampleID >> {output} " 

rule combined_sudan_haplogroup_file:
    input: expand("haplocheck/results/SUDAN2020_Sudanese_{individual}/haplogroups/haplogroups.txt", \
                    individual=INDIVIDUALS_SUDAN) 
    output: "haplocheck/results/sudan_haplogroups.txt"
    shell: "cat {input[0]} | head -n 1 > {output}; " + \
            "cat haplocheck/results/SUDAN2020*/haplogroups/haplogroups.txt | " + \
            "grep -v SampleID >> {output} " 

rule combined_sudan_contamination_file:
    input: expand("haplocheck/results/SUDAN2020_Sudanese_{individual}/contamination/contamination.txt", \
                    individual=INDIVIDUALS_SUDAN) 
    output: "haplocheck/results/sudan_contaminations.txt"
    shell: "cat {input[0]} | head -n 1 > {output}; " + \
            "cat haplocheck/results/SUDAN2020*/contamination/contamination.txt | " + \
            "grep -v Sample >> {output} " 


