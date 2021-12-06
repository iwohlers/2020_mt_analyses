in_config_file = snakemake.input[0]

in_1000g = snakemake.input[1]
in_bergstroem = snakemake.input[2]
in_egypt = snakemake.input[3]
in_german = snakemake.input[4]
in_schuenemann = snakemake.input[5]
in_sudan = snakemake.input[6]
in_wohlers = snakemake.input[7]

cont_1000g = snakemake.input[8]
cont_bergstroem = snakemake.input[9]
cont_egypt = snakemake.input[10]
cont_german = snakemake.input[11]
cont_schuenemann = snakemake.input[12]
cont_sudan = snakemake.input[13]
cont_wohlers = snakemake.input[14]

out_file = snakemake.output[0]

in_files = [in_1000g,in_bergstroem,in_egypt,in_german,in_schuenemann,in_sudan,in_wohlers]
cont_files = [cont_1000g,cont_bergstroem,cont_egypt,cont_german,cont_schuenemann,cont_sudan,cont_wohlers]
dataset_handle = ["1000G","BERGSTROEM2020","EGYPT2020","GERMAN2021","SCHUENEMANN2017","SUDAN2020","WOHLERS2020"]


# Read in the input files with population and individuals
samples = {}
for filename,handle in zip(in_files,dataset_handle):
    with open(filename,"r") as f_in:
        for line in f_in:
            s = line.strip().split("\t")
            sample = {}
            sample["dataset"] = handle
            sample["population"] = s[0]
            sample["individual"] = s[1]
            samples[s[1]] = sample

# Read in the haplocheck input files
for filename,handle in zip(cont_files,dataset_handle):
    with open(filename,"r") as f_in:
        for line in f_in:
            if line[:7] == "\"Sample":
                header = line
                continue
            s = line.strip().split("\t")
            info = [x.replace('"','') for x in s]
            sample_id = info[0].replace("_mt.bam","")
            sample = samples[sample_id]
            sample["dataset"] = handle
            sample["contamination"] = info[1]
            sample["coverage"] = info[4]
            sample["major_haplogroup"] = info[9]
            sample["minor_haplogroup"] = info[11]
            sample["haplocheck_line"] = line
            samples[s[0]] = sample

#print(samples)

# Go through the config file line by line and remember all samples that
# fulfill at least one config line
samples_included = {}
with open(in_config_file,"r") as f_in, open(out_file,"w") as f_out:
    for line in f_in:
        if line[0] == '#':
            f_out.write("#dataset"+"\t"+"individual"+"\t"+"population"+"\t"+header.replace("\"",""))
            continue
        s = line.strip("\n").split("\t")
        dataset = s[0]
        population = s[1]
        individual = s[2]
        haplogroup = s[3]
        min_coverage = s[4]
        contamination = s[5]
        major_minor_same = s[6]
        for sample_index in samples:
            sample = samples[sample_index]
            #print(dataset)
            # From the Sudanese data set, some samples have no haplocheck output
            # Skip those
            if not "haplocheck_line" in sample:
                print("Sample "+sample["individual"]+" has no haplocheck output")
                continue
            if not dataset == "" and not dataset == sample["dataset"]:
                continue
            if not population == "" and not population == sample["population"]:
                continue
            if not individual == "" and not individual == sample["individual"]:
                continue
            if not haplogroup == "" and not haplogroup == sample["major_haplogroup"][:len(haplogroup)]:
                continue
            if not min_coverage == "" and not int(sample["coverage"])>=int(min_coverage):
                continue
            if not contamination == "" and not contamination == sample["contamination"]:
                continue
            if not major_minor_same == "" and not sample["major_haplogroup"]==sample["minor_haplogroup"] and not major_minor_same == "NO":
                continue    
            if not major_minor_same == "" and sample["major_haplogroup"]==sample["minor_haplogroup"] and not major_minor_same == "YES":
                continue
            if not sample["individual"] in samples_included:
                f_out.write(sample["dataset"]+ "\t"+sample["individual"]+ "\t"+sample["population"]+"\t"+sample["haplocheck_line"].replace("\"",""))
                samples_included[sample["individual"]] = True
                
                
                    

