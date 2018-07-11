######Function to read h5 file given specification like cancer type, Gene names or patient ID
######Author: Ni Shuai and Bernd Fischer (n.shuai[at]dkfz.de)
library(rhdf5)

##ExampleDB=H5readRNAseq(h5file = 'TCGA_RNAseq_raw_count_Feb-26-2016.h5',TumorType = 'BLCA')
H5readRNAseq=function(h5file, TumorType, PatientID, Genename){
  if (nargs()<1) {
    
    print ("Need to specify at leat one parameter from TumprType, PatientID and Genename.");
    if (!exists('TCGA_Tumor_Type')) {assign('TCGA_Tumor_Type',c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
                                                                'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                                                                'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV',
                                                                'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM',
                                                                'STAD', 'THCA', 'UCEC', 'UCS'),envir = .GlobalEnv)
      print ("The possible input tumor types has been saved in Variable 'TCGA_Tumor_Type'" )}
    else { 
      print ("The possible input tumor types are: " )
      print (c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
               'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 
               'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 
               'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM', 
               'STAD', 'THCA', 'UCEC', 'UCS'))}
    
    return ()
  }
  
  if (missing(h5file)) print ('No h5 file specified')
  if (!file.exists(file.path(h5file))) stop('H5 file not exist')
  R_version=h5read(h5file,'Version/R_Version')
  Datebase_date=h5read(h5file,'Version/Datebase_Version')
  Patient_Tumor=h5read(h5file,'patient_info/RNAseq_Patient_info')
  Patients=Patient_Tumor[,1]
  Tumors=Patient_Tumor[,2]
  Genes=h5read(h5file,'Gene_info/RNAseq_Gene_info')
  dim1=1:length(Genes) 
  dim2=1:length(Patients)
  Selected_gene=rep(TRUE,length(Genes))
  Selected_patient=rep(TRUE,length(Patients))
  
  
  if (!missing(TumorType)) Selected_patient=Tumors %in% TumorType
  if (!missing(PatientID)) Selected_patient=Selected_patient & (Patients %in% PatientID)
  if (!missing(Genename))  Selected_gene=Genes %in% Genename
  if (sum(Selected_gene)==0 | sum(Selected_patient)==0) stop('There is no record matches your specification')
  sss=h5read(h5file,'Data_matrix/RNAseq_Datamatrix', index=list(dim1[Selected_gene],dim2[Selected_patient]))
  output=list()
  output$Patient_info=Patient_Tumor[Selected_patient,]
  colnames(output$Patient_info)=c('PatientID','TumorType','Tissue','Center','Site','AbbrSite')
  output$TumorType=Tumors[Selected_patient]
  output$GeneName=Genes[Selected_gene]
  row.names(sss)=output$GeneName
  colnames(sss)=Patients[Selected_patient]
  output$Data=sss
  output$R_version=R_version
  output$Datebase_date=Datebase_date
  
  H5close()
  output 
  #sss=h5read(h5file,'Data_matrix/RNAseq_Datamatrix', index=list(dim1[Selected_gene], dim2[Selected_patient]))
}



######Function to read back the h5 file given specification like cancer type, Gene names or patient ID
H5readCNV=function(h5file, TumorType, PatientID, Genename){
  if (nargs()<1) {
    print ("Need to specify at leat one parameter from TumprType, PatientID and Genename.");
    if (!exists('TCGA_Tumor_Type')) {assign('TCGA_Tumor_Type',c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
                                                                'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP',
                                                                'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV',
                                                                'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM',
                                                                'STAD', 'THCA', 'UCEC', 'UCS'),envir = .GlobalEnv)
      print ("The possible input tumor types has been saved in Variable 'TCGA_Tumor_Type'" )}
    else { 
      print ("The possible input tumor types are: " )
      print (c('ACC','BLCA', 'BRCA', 'CESC', 'COAD', 'DLBC', 
               'ESCA', 'GBM', 'HNSC', 'KICH', 'KIRC', 'KIRP', 
               'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'OV', 
               'PAAD', 'PRAD', 'READ', 'SARC', 'SKCM', 
               'STAD', 'THCA', 'UCEC', 'UCS'))}
    
    return ()
  }
  
  if (missing(h5file)) stop ('No h5 file specified')
  
  R_version=h5read(h5file,'Version/R_Version')
  Datebase_date=h5read(h5file,'Version/Datebase_Version')
  Patient_Tumor=h5read(h5file,'patient_info/CNV_Patient_info')
  Patients=Patient_Tumor[,1]
  Tumors=Patient_Tumor[,2]
  Genes=h5read(h5file,'Gene_info/CNV_Gene_info')[,1]
  dim1=1:length(Genes) 
  dim2=1:length(Patients)
  Selected_gene=rep(TRUE,length(Genes))
  Selected_patient=rep(TRUE,length(Patients))
  
  
  if (!missing(TumorType)) Selected_patient=Tumors %in% TumorType
  if (!missing(PatientID)) Selected_patient=Selected_patient & (Patients %in% PatientID)
  if (!missing(Genename))  Selected_gene= Genes %in% Genename
  if (sum(Selected_gene)==0 | sum(Selected_patient)==0) stop('There is no record matches your specification')
  sss=h5read(h5file,'Data_matrix/CNV_Datamatrix', index=list(dim1[Selected_gene],dim2[Selected_patient]))
  output=list()
  output$PatientID=Patients[Selected_patient]
  output$TumorType=Tumors[Selected_patient]
  output$GeneName=Genes[Selected_gene]
  output$Data=sss
  output$R_version=R_version
  output$Datebase_date=Datebase_date
  H5close()
  output 
}

#ss=H5readCNV('TCGA_CNV_Feb-26-2016.h5', TumorType='BLCA')
#ss$PatientID
#str(ss)
######Download clinical data
