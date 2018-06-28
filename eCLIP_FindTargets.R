library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(org.Hs.eg.db)
library(GenomicAlignments)
library(gdata)
library(idr)
####author: Ni Shuai and Bernd Fischer
###read the metadata for all downloaded data
load('~/wrk/bernd/eCLIP/eCLIPannotation.rda')
###remove columns that are the same for every row, they are not informative
rm_redunt_col=function(Df){
    I=apply(Df, 2, function(x) length(unique(x)))
    Df[,I!=1]
}
###clean up the meta file for peaks
peaks=eCLIPannotation
peaks=peaks[peaks$Biosample.term.name !='adrenal gland',]
peaks=rm_redunt_col(peaks)
peaks$sample_name=unlist(strsplit(peaks$Experiment.target, split='-'))[c(TRUE, FALSE)]
peaks$sample_name=paste(peaks$sample_name, peaks$Biosample.term.name, sep='-')

####get peaks into one single target file
peak_files=list.files('~/wrk/bernd/eCLIP/bed/')
peak_files_accession=unlist(strsplit(peak_files, split='[.]'))[c(TRUE,FALSE,FALSE)]

bam_files=list.files('~/wrk/bernd/eCLIP/bam/')
bam_files_accession=unlist(strsplit(bam_files, split='[.]'))[c(TRUE,FALSE)]
##we have 168 experiments in total
length(bam_files_accession)/3
length(peak_files_accession)/2

####for each protein-cellline combination find the corresponding peak and bam file

extraCols_narrowPeak <- c(sigvalue= "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer") 
import.narrowPeak <- function(... ) {
    import(..., format="bed", extraCols=extraCols_narrowPeak)
}

###function that finds corresponding file accession for each protein-cellline combination

getBam=function(sample_name){
    peaks$File.accession[peaks$File.format=='bam' & peaks$sample_name==sample_name]
}
getBed=function(sample_name){
    peaks$File.accession[peaks$File.format=='bed narrowPeak' & peaks$sample_name==sample_name]
}
getMockbam=function(sample_name){
    protein=unlist(strsplit(sample_name, split='-'))[1]
    cellline=unlist(strsplit(sample_name, split='-'))[2]
    peaks$File.accession[peaks$File.format=='bam' & 
                             grepl('eCLIP mock', peaks$sample_name) & grepl(protein, peaks$sample_name) & 
                             grepl(cellline, peaks$sample_name)]
}

##get all experiment names
Exps=unique(peaks$sample_name[-grep('input', peaks$sample_name)])

for (i in 1:length(Exps)){
    ##define sample to work with
    sample_name=Exps[i]
    ##import 2 bed files for that sample
    x1 <-import.narrowPeak(file.path('~/wrk/bernd/eCLIP/bed', paste0(getBed(sample_name)[1], '.bed.gz')))
    x2 <-import.narrowPeak(file.path('~/wrk/bernd/eCLIP/bed', paste0(getBed(sample_name)[2], '.bed.gz')))
    ###find overlapped regions
    x=findOverlapsOfPeaks(x1, x2, minoverlap = 0, maxgap = 0, ignore.strand = FALSE)
    overlaps=x$peaklist[["x1///x2"]]
    ###get the coverage from corresponding bam files
    bamFileA <-file.path('~/wrk/bernd/eCLIP/bam', paste0(getBam(sample_name)[1], '.bam'))
    bamFileB <-file.path('~/wrk/bernd/eCLIP/bam', paste0(getBam(sample_name)[2], '.bam'))
    bamFileMock <-file.path('~/wrk/bernd/eCLIP/bam', paste0(getMockbam(sample_name), '.bam'))
    ###get coverage from 2 biological replicates
    coverage <- summarizeOverlaps(features = overlaps, 
                  reads = c(bamFileA, bamFileB),
                  mode=Union, 
                  ignore.strand = FALSE, 
                  singleEnd=TRUE)
    
    ###get coverage from Mock input
    Mockcoverage <- summarizeOverlaps(features = overlaps, 
                  reads = bamFileMock,
                  mode=Union, 
                  ignore.strand = FALSE, 
                  singleEnd=TRUE)
 
    ###remove entries that is irreproducible
    idr=est.IDR(assay(coverage)/ width(rowRanges(coverage)),
             mu=2.07, sigma=1.34, rho=0.89, p=0.84)
    overlaps=overlaps[idr$IDR<=0.3]
    coverage=coverage[idr$IDR<=0.3]
    Mockcoverage=Mockcoverage[idr$IDR<=0.3]
    #annotate those overlaps
    ucsc.hg38.knownGene <-genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
    anno <- annotatePeakInBatch(overlaps, 
                                     AnnotationData=ucsc.hg38.knownGene)
    anno=anno[!duplicated(anno),]
    coverage=coverage[anno$insideFeature=='inside' & !is.na(anno$insideFeature),]
    Mockcoverage=Mockcoverage[anno$insideFeature=='inside' & !is.na(anno$insideFeature),]
    anno=anno[anno$insideFeature=='inside' & !is.na(anno$insideFeature),]
    anno=anno[,2:3]
    anno <- addGeneIDs(annotatedPeak=anno, 
                            orgAnn="org.Hs.eg.db", 
                            feature_id_type="entrez_id",
                            IDs2Add="symbol")
    ##remove rare cases that genes overlap with each other, in this case only one gene is taken
    anno <- addGeneIDs(annotatedPeak=anno, 
                            orgAnn="org.Hs.eg.db", 
                            feature_id_type="entrez_id",
                            IDs2Add="ensembl")
    

    ###do Deseq analysis
    counttable=data.frame(coverage1=assay(coverage)[,1],coverage2=assay(coverage)[,2], coverageMock=assay(Mockcoverage)[,1])
    
    rownames(counttable)=names(anno)
    coldata=data.frame(condition=c('treated','treated','untreated'), type=c('paired-end','paired-end','paired-end'))
    rownames(coldata)=names(counttable)
    
    dds <- DESeqDataSetFromMatrix(countData =counttable,
                                  colData = coldata,
                                  design = ~ condition)
    
    
    dds=DESeq(dds)
    res=results(dds)
    pack=list(dds=dds, res=res, peaks=anno)
    saveRDS(pack,file.path('results', paste0(sample_name, '.rds')))
}

####after having significant peaks, get all interacting RNA for each experiment

##get unique peaks inside the the feature and with p value smaller than or equal 0.01
targets=list.files('~/wrk/bernd/eCLIP/results/')
Target=lapply(1:length(targets), function(i){
    target=readRDS(file.path('~/wrk/bernd/eCLIP/results/', targets[i]))
    peaknames=rownames(target$res)[target$res$padj<=0.01 & !is.na(target$res$padj)]
    unique(target$peaks$ensembl[names(target$peaks) %in% peaknames])
})

names(Target)=unlist(strsplit(targets, split='[.]'))[c(TRUE, FALSE)]
saveRDS(Target, '~/wrk/bernd/eCLIP/AllEclipTargets.rds')

####testing
load(file = file.path('/Users/sni/wrk/bernd/MouseNeuralDevelopment/eCLIP/AUH-HepG2.rda'))
bb=peaks$geneID[res$padj<=0.1 & !is.na(res$padj)]
bb %in% anno$ensembl
