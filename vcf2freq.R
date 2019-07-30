vcf2freq <- function(input,high_parent,low_parent,output,cycle)
{
    if(missing(input))
    {
        stop("ERROR: An input vcf file should be provided")
    }
    if(missing(output))
    {
        out_file = paste0(input,".freq")
    }else
    {
        out_file = output
    }
    con1 = file(input,"read")
    while(TRUE)
    {
        line = readLines(con1,1)
        if(grepl("#CHROM",line)){break}
    }
    header = strsplit(line,"\t")[[1]]
    close(con1)
    ## by default the high perent is the first genotype in the file and the low parent is the second one,
    ## otherwise, their were specified by the user.
    if(missing(high_parent))
    {
        high_p = 10
    }else
    {
        high_p =  which(header == high_parent)
    }
    
    if(missing(low_parent))
    {
        low_p = 11
    }else
    {
        low_p =  which(header == low_parent)
    }
    ##Reading the genotypic data from the vcf file
    vcf = read.table(input)
    OUTPUT = data.frame()
    high_geno = vcf[,high_p]
    low_geno = vcf[,low_p]
    format_geno = vcf[,9]
    count = 1
    for(i in 1:length(format_geno))
    {
        if(grepl("AD",as.character(format_geno[i])))
        {
            AD_index = which(strsplit(as.character(format_geno[i]),':')[[1]]=="AD")
            ##getting high freq
            high = strsplit(as.character(high_geno[i]),':')[[1]][AD_index]
            high_ref = as.numeric(strsplit(high,",")[[1]][1])
            high_alt = as.numeric(strsplit(high,",")[[1]][2])
            high_freq = high_ref/(high_ref+high_alt)
            ##getting low freq
            low = strsplit(as.character(low_geno[i]),':')[[1]][AD_index]
            low_ref = as.numeric(strsplit(low,",")[[1]][1])
            low_alt = as.numeric(strsplit(low,",")[[1]][2])
            low_freq = low_ref/(low_ref+low_alt)
            ### getting delta SNP
            deltaSNP = high_freq - low_freq
            ### filling the data frame
            OUTPUT = rbind.data.frame(OUTPUT,cbind.data.frame(vcf[i,1:2],high_freq,low_freq,deltaSNP,cycle))
            
            if(i > 1 && (vcf[i,1] != vcf[i-1,1] || i == length(format_geno)))
            {
                print(paste0("SNPs and delta SNPs of ",vcf[i,1]," were calculated >> ",count," SNPs were extracted"))
                count=0
            }
            count = count+1
        }
    }
	colnames(OUTPUT) = c("chr","position","highBulkFreq","lowBulkFreq","deltaSNP","cycle")
    write.table(OUTPUT,out_file,col.names=T,row.names=FALSE,sep='\t',quote=FALSE)
}
args <- commandArgs(TRUE)
vcf2freq(args[1],args[2],args[3],args[4],args[5])