####Downlaod and super fast retrive of ALL TCGA RNAseq data in h5 format
####Author: Ni Shuai and Bernd Fischer

######load packages
library(rhdf5)
required_lib=(c("HGNChelper", "RCurl", "httr", 
    "stringr", "digest", "bitops"))
lapply(required_lib, require, character.only = TRUE)
library(org.Hs.eg.db)
source('Module_A.r')
source('Module_B.r')

#database_date='Feb-26-2016'
######traverseALLdirectories TCGA
TraverseAllDirectories(
  entryPoint= 
  "htts://tcga-data.nci.nih.gov/tcgafiles/
  ftp_auth/distro_ftpusers/anonymous/tumor/", 
  fileLabel ="DirectoryTraverseResult");
##it is also possible to traverse cancer type specific URLs only. 
###TCGA assembler automatically name the traverse file in the same manner as below
###date is the local date
database_date=paste(unlist(strsplit(format(Sys.Date(), '%b, %d, %Y'), ', ')), collapse='-')
TraverseResult=paste("DirectoryTraverseResult_", database_date, '.rda', sep='')


######Download all RNAseqV2 data for all cancer types
DirRnaseqV2='RNAseq'
dir.create(file.path(DirRnaseqV2), showWarnings = FALSE)
for (cancer_type in c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
            'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 
            'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 
            'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM', 
            'STAD', 'THCA', 'UCEC', 'UCS')) {
  DownloadRNASeqData(
  traverseResultFile = file.path('.', TraverseResult), 
  saveFolderName = file.path(DirRnaseqV2), 
  cancerType =cancer_type, assayPlatform = "RNASeqV2", 
  dataType = c("rsem.genes.results", "rsem.genes.normalized_results"));
}


######Read and render datasets into a data.frames and .rda files
OutputDir='RNAseqV2Counts_Rda'
dir.create(file.path('.', OutputDir), showWarnings = FALSE)
File_list=list.files(file.path('.', DirRnaseqV2))
File_list_result=File_list[grep('genes.results',File_list)]
File_list_norm_result=File_list[grep('genes.normalized',File_list)]

for (file in File_list_result){
  cancer_type=unlist(strsplit(file, '__'))[1]
  pltfm=unlist(strsplit(file, '__'))[3]
  rnaseqv2=read.table(file.path(DirRnaseqV2, file), 
  sep='\t', header=TRUE, row.names=1,stringsAsFactors = FALSE)
  rnaseqv2=rnaseqv2[,-1]
  raw_count=rnaseqv2[,grep('raw', rnaseqv2[1,])]
  raw_count=raw_count[-1,]
  scaled_estimate=rnaseqv2[,grep('scaled',rnaseqv2[1,])]
  scaled_estimate=scaled_estimate[-1,]
  names(scaled_estimate)=sapply(
  names(scaled_estimate), function(x) 
    paste(head(unlist(strsplit(x, "[.]")),-1), collapse ='.'))
  
  saveRDS(raw_count, file=file.path(
     OutputDir, paste('raw_count',
                        cancer_type, 
                        pltfm,
                        'rda',sep='.')))
  
  print(paste('saved raw counts for', cancer_type ))
    saveRDS(scaled_estimate, file=file.path(
     OutputDir, paste('scaled_estimate',
                                   cancer_type,
                                   pltfm,
                        'rda',sep='.')))
    
  print(paste('saved scaled estimate for', cancer_type ))
}

for (file in File_list_norm_result){
   cancer_type=unlist(strsplit(file, '__'))[1]
   pltfm=unlist(strsplit(file, '__'))[3]
   norm_rnaseqv2=read.table(file.path(DirRnaseqV2, file),
  sep='\t', header=TRUE, row.names=1,stringsAsFactors = FALSE)
  print(file)
  norm_rnaseqv2=norm_rnaseqv2[-1,]
  saveRDS(norm_rnaseqv2, file=file.path(
    OutputDir, paste('normalized_count',
                                  cancer_type,
                                  pltfm,
                       'rda',sep='.')))
  print(paste('saved normalized count for', cancer_type ))
} 



# #Convert gene names into ENSEMBL names
    # x <- org.Hs.egENSEMBL
    # mapped_genes <- mappedkeys(x)
    # xx <- as.list(x[mapped_genes])
    # unlist(xx[c('145376','145389')]) 
    # ####Convert the Entrez ID to ENSEMBL ID
    # ENSGid=xx[sapply(row.names(raw_count), 
    #   function(x) unlist(strsplit(x, "[|]"))[2])]
    # ENSGid[listLen(ENSGid)!=1]=NA
    # ENSGid=unname(unlist(ENSGid))
    # row.names(raw_count)=ENSGid
    # 
    # row.names(scaled_estimate)=xx[
    #   sapply(row.names(scaled_estimate), 
    #          function(x) unlist(strsplit(x, "[|]"))[2])]
    # scaled_estimate[1:8,1:8]


#####Gene identifiers are the same for all 27 datasets
#####Only sequencing data from Hiseq is used to aviod RNA-seq batch effects of sequencing platforms bwtween Hiseq series and Genome Analyzer II 

DirRnaseqV2='RNAseq'
OutputDir='RNAseqV2Counts_Rda'

CountType=c('raw_count', 'normalized_count','scaled_estimate')
for (counttype in CountType){ 
  files=list.files(file.path(OutputDir))
  files=files[grep('hiseq', files)]
  files=files[grep(counttype, files)]
  ####read an example file to initialize the matrices
  
  example=readRDS(file.path(OutputDir,files[1]))
  example=as.matrix(example)
  class(example)='numeric'
  Genenames=vapply(
    strsplit(row.names(example),split='[|]'), 
    function(x) x[1], 'ni')
  example=example[Genenames!='?',]
  Genenames=vapply(
    strsplit(row.names(example),split='[|]'), 
    function(x) x[1], 'ni')
  
  RNAseq_Datamaatrix=matrix(,nrow=nrow(example), ncol=0)
  RNAseq_Gene_info=Genenames
  RNAseq_Patient_info=matrix(, nrow=0, ncol=2)
  colnames(RNAseq_Patient_info)=c('PatientID', 'TumorType')
  
  ###fill in the initiated data matrices
  for (i in files){
    cancer=unlist(strsplit(i, '[.]'))[2]
    example=readRDS(file.path(OutputDir,i))
    Genenames=vapply(
      strsplit(row.names(example),split='[|]'), 
        function(x) x[1], 'ni')
    example=example[Genenames!='?',]
    Genenames=vapply(
      strsplit(row.names(example),split='[|]'), 
        function(x) x[1], 'ni')
    RNAseq_Patient_info=rbind(
      RNAseq_Patient_info,
      matrix(c(colnames(example), rep(cancer, ncol(example))), ncol=2))
    
    example=as.matrix(example)
    row.names(example)=NULL; colnames(example)=NULL;
    class(example)='numeric' 
    RNAseq_Datamaatrix=cbind(RNAseq_Datamaatrix,example)
  }
 
  ####Add explainary information to the patient IDs
  center_code=read.table('code_sampleID/centerCode_TCGA.txt',sep='\t', header=TRUE,quote = "", na.string='z1x2cv4b81nM')
  sample_code=read.table('code_sampleID/sampleTypeCode_TCGA.txt',sep='\t', header=TRUE,quote = "", na.string='NiS2zx235tPx')
  tissue_code=read.table('code_sampleID/tissueSourceSiteCode_TCGA.txt',sep='\t', header=TRUE,quote = "", na.string='nasS8nAg55Asa')
  
  PatientIDs= strsplit(RNAseq_Patient_info[,1], split='[.]')
  TissueName=as.character(sapply(PatientIDs, 
         function(x) {
            tissue_code$Study.Name[tissue_code$TSS.Code %in% x[2]]
         }))
  CenterName=as.character(sapply(PatientIDs, 
         function(x) {
            center_code$Display.Name[center_code$Code %in%  as.numeric(x[7])]
         }))
  SampleName=as.character(sapply(PatientIDs, 
         function(x) {
           sample_code$Definition[sample_code$Code %in%  as.numeric(substr(x[4], 1,2))]
         }))
  ShortSampleName=as.character(sapply(PatientIDs, 
                               function(x) {
                                 sample_code$Short.Letter.Code[sample_code$Code %in%  as.numeric(substr(x[4], 1,2))]
                               }))
  RNAseq_Patient_info=cbind(RNAseq_Patient_info,TissueName)
  RNAseq_Patient_info=cbind(RNAseq_Patient_info,CenterName)
  RNAseq_Patient_info=cbind(RNAseq_Patient_info,SampleName)
  RNAseq_Patient_info=cbind(RNAseq_Patient_info,ShortSampleName)
  
  ####create hdr5 dataset
  h5RNAseqname=paste(paste('TCGA_RNAseq',counttype, database_date,sep='_'),'h5',sep='.')
  h5createFile(h5RNAseqname)
  h5createGroup(h5RNAseqname,'patient_info')
  h5createGroup(h5RNAseqname,'Gene_info')
  h5createGroup(h5RNAseqname,'Data_matrix')
  h5createGroup(h5RNAseqname,'Version')
  h5write(RNAseq_Datamaatrix,h5RNAseqname,'Data_matrix/RNAseq_Datamatrix', level=0)
  h5write(RNAseq_Gene_info,h5RNAseqname,'Gene_info/RNAseq_Gene_info')
  h5write(RNAseq_Patient_info,h5RNAseqname,'patient_info/RNAseq_Patient_info')
  h5write(database_date, h5RNAseqname,'Version/Datebase_Version')
  h5write(R.Version()[[1]],h5RNAseqname,'Version/R_Version')
  h5ls(h5RNAseqname)
  H5close()
}  


#######################
######Download CNV data
#######################
###Cancer type 'LAML' could not be downloaded
DirCNV='CNV'
dir.create(file.path(DirCNV), showWarnings = FALSE)
for (cancer_type in c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
                      'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                      'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 
                      'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM', 
                      'STAD', 'THCA', 'UCEC', 'UCS')) {
  DownloadCNAData(
    traverseResultFile =file.path(TraverseResult), 
    saveFolderName = file.path(DirCNV), 
    assayPlatform = 'genome_wide_snp_6',
    cancerType =cancer_type)
}


######process CNV files using TCGA assembler build-in function

####It takes a little while!
DIRCNVgene='GeneLevelCNV'
DirCNV='CNV'
dir.create(file.path(DIRCNVgene), showWarnings = FALSE)
###rename function
rename=function(filename){
  aa=unlist(strsplit(filename, '__'))
  paste('GeneLevelCNA',
        paste(as.character(aa[c(-2, -length(aa))]), 
              collapse ='__'), database_date,sep='__')
}

files=list.files(file.path(DirCNV))
files=files[grep('nocnv_hg19', files)]
for (filename in files){
  ProcessCNAData(inputFilePath = file.path(DirCNV,filename),
                 outputFileName =rename(filename),
                 outputFileFolder = file.path(DIRCNVgene),
                 refGenomeFile =file.path('SupportingFiles/Hg19GenePosition.txt')
               )
  }


######read the rda files into a hdf5 object

#####The gene names are identical for all cancer types, therefore just save one copy of that
#####First try to read one rda file to initialize the matrices

files=list.files(file.path(DIRCNVgene))
files=files[grep('.rda', files)]
load(file.path(DIRCNVgene, files[1]))
CNV_Datamaatrix=matrix( , nrow=nrow(Data), ncol=0)
CNV_Gene_info=Des
CNV_Patient_info=matrix( , nrow=0, ncol=2)
colnames(CNV_Patient_info)=c('PatientID', 'TumorType')

for (i in files){
  load(file.path(DIRCNVgene, i))
  cancer=unlist(strsplit(i, '__'))[2]
  CNV_Patient_info=rbind(
    CNV_Patient_info,
    matrix(c(colnames(Data), 
             rep(cancer, ncol(Data))), ncol=2))
  Data=as.matrix(Data)
  colnames(Data)=NULL; row.names(Data)=NULL
  CNV_Datamaatrix=cbind(CNV_Datamaatrix,Data)
}

####Add explainary information to the patient IDs
center_code=read.table('code_sampleID/centerCode_TCGA.txt',sep='\t', header=TRUE,quote = "", na.string='z1x2cv4b81nM')
sample_code=read.table('code_sampleID/sampleTypeCode_TCGA.txt',sep='\t', header=TRUE,quote = "", na.string='NiS2zx235tPx')
tissue_code=read.table('code_sampleID/tissueSourceSiteCode_TCGA.txt',sep='\t', header=TRUE,quote = "", na.string='nasS8nAg55Asa')

PatientIDs= strsplit(CNV_Patient_info[,1], split='-')
TissueName=as.character(sapply(PatientIDs, 
                               function(x) {
                                 tissue_code$Study.Name[tissue_code$TSS.Code %in% x[2]]
                               }))
CenterName=as.character(sapply(PatientIDs, 
                               function(x) {
                                 center_code$Display.Name[center_code$Code %in%  as.numeric(x[7])]
                               }))
SampleName=as.character(sapply(PatientIDs, 
                               function(x) {
                                 sample_code$Definition[sample_code$Code %in%  as.numeric(substr(x[4], 1,2))]
                               }))
ShortSampleName=as.character(sapply(PatientIDs, 
                               function(x) {
                                 sample_code$Short.Letter.Code[sample_code$Code %in%  as.numeric(substr(x[4], 1,2))]
                               }))
                                 
CNV_Patient_info=cbind(CNV_Patient_info,TissueName)
CNV_Patient_info=cbind(CNV_Patient_info,CenterName)
CNV_Patient_info=cbind(CNV_Patient_info,SampleName)
CNV_Patient_info=cbind(CNV_Patient_info,ShortSampleName)

####Put data matrices into a single h5 file
h5CNVname=paste(paste('TCGA_CNV',database_date,sep='_'),'h5',sep='.')
h5createFile(h5CNVname)
h5createGroup(h5CNVname,'patient_info')
h5createGroup(h5CNVname,'Gene_info')
h5createGroup(h5CNVname,'Data_matrix')
h5createGroup(h5CNVname,'Version')

h5write(CNV_Datamaatrix,h5CNVname,'Data_matrix/CNV_Datamatrix')
h5write(CNV_Gene_info,h5CNVname,'Gene_info/CNV_Gene_info')
h5write(CNV_Patient_info,h5CNVname,'patient_info/CNV_Patient_info')
h5write(database_date, h5CNVname,'Version/Datebase_Version')
h5write(R.Version()[[1]],h5CNVname,'Version/R_Version')
H5close()

######Download clinical data

DirClinic='Clinic'
dir.create(file.path(DirClinic), showWarnings = FALSE)
for (cancer_type in c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
                      'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 
                      'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 
                      'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM', 
                      'STAD', 'THCA', 'UCEC', 'UCS')) {
  DownloadClinicalData(
    traverseResultFile = file.path(TraverseResult), 
    saveFolderName = file.path(DirClinic), 
    cancerType =cancer_type, 
    clinicalDataType= c('patient','drug', 'follow_up','radiation'))
}

####Render clinical data into big matrices
# 
# 
# files=list.files(file.path(DirClinic))
# files=files[grep('patient', files)]
# file=files[2]
# for (file in files){
#   ss=read.table(file.path(DirClinic, file), quote = "\"'",skip = 0,check.names = TRUE,
#   dec='.',   sep='\t',stringsAsFactors =TRUE)
#   aa=file(description=file.path(DirClinic, file), open='r')
#   readLines(con = file, n = -1L, ok = TRUE, warn = TRUE,
#             encoding = "unknown", skipNul = FALSE)
#   
#   names(ss)=ss[1,]
#   ss=ss[-1:-2, ]
#   names(ss)
#   print(dim(ss))}
#   ss$bcr_patient_barcode
#   sss=H5readCNV(Dir='.',h5file = 'TCGA_CNV_Feb-26-2016.h5',TumorType = 'ACC')
#   str(sss)
#   order(ss$bcr_patient_barcode)
#   length(ss$bcr_patient_barcode)
#   sss$PatientID
#   substr(sss$PatientID, 1, 12) %in% ss$bcr_patient_barcode
# ###180 records in cnv data, 92 records in patient data. 
# ###there are 88 duplicated records in cnv data, which are tumatched data, these can be ideitified by the 5th position of the conjuction. 
# ###Task1, identify tumor matched data
# ###Task2, match cnv data to patient information
