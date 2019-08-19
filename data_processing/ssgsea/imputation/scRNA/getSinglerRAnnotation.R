#' Wrapper function to create a SingleR object
#'
#' @param counts single cell expression matrix with gene row and sample column. 
#' @param project.name the project name
#' @param min.genes Include cells where at least this many genes are detected (number non-zero genes).
#' @param technology The technology used for creating the single-cell data.
#' @param species The species of the sample ('Human' or 'Mouse').
#' @param citation a citation for the project.
#' @param ref.list a list of reference objects. If NULL uses the predefined reference objects - Mouse: ImmGen and Mouse.RNAseq, Human: HPCA and Blueprint+Encode. 
#' @param normalize.gene.length if a full-length method set to TRUE, if a 3' method set to FALSE.
#' @param variable.genes variable gene method to use - 'sd' or 'de'. Default is 'de'.
#' @param fine.tune perform fine tuning. Default is TRUE. Fine-tuning may take long to run.
#' @param do.signatures create signatures data
#' @param clusters input cluster id for each of the cells with at least min.genes, if NULL uses SingleR clusterings.
#' @param do.main.types run the SingleR pipeline for main cell types (cell types grouped together) as well.
#' @param reduce.file.size remove less used SingleR fields that increase the object size.
#' @param temp.dir used by the SingleR webtool.
#' @param numCores Number of cores to use.
#'
#' @return a SingleR object
mySCObject = function(counts, project.name,
                      min.genes=0,technology='10X',
                      species='Human',citation='',
                      ref.list=list(),normalize.gene.length=F,
                      variable.genes='de',fine.tune=T,
                      do.signatures=F,clusters=NULL,
                      do.main.types=T,reduce.file.size=T,
                      temp.dir=NULL,numCores = SingleR.numCores) {
    
    require(SingleR)
    
    
    print(paste0('Dimensions of counts data: ',
                 nrow(counts),'x',ncol(counts)))
    
    singler = list()
    
    
    N = colSums(counts>0)
    counts = counts[,N>=min.genes]
    orig.ident = colnames(counts)
    
    sc.data.gl = counts 
    
    if (length(ref.list)==0) {
        if (species == 'Mouse') {
            #if (!exists('immgen'))
            #  data('Immgen')
            #if (!exists('mouse.rnaseq'))
            #  data('Mouse-RNAseq')
            res = list(SingleR.CreateObject(sc.data.gl,immgen,clusters,species,
                                            citation,technology,
                                            do.main.types=do.main.types,
                                            variable.genes=variable.genes,
                                            fine.tune=fine.tune,numCores = numCores),
                       SingleR.CreateObject(sc.data.gl,mouse.rnaseq,clusters,
                                            species,citation,technology,
                                            do.main.types=do.main.types,
                                            variable.genes=variable.genes,
                                            fine.tune=fine.tune,numCores = numCores)
            )
        } else if (species == 'Human') {
            #if(!exists('hpca'))
            #  data ('HPCA')
            #if (!exists('blueprint_encode'))
            #  data('Blueprint_Encode')
            res = list(SingleR.CreateObject(sc.data.gl,hpca,clusters,species,
                                            citation,technology,
                                            do.main.types = do.main.types,
                                            variable.genes=variable.genes,
                                            fine.tune=fine.tune,numCores = numCores),
                       SingleR.CreateObject(sc.data.gl,blueprint_encode,
                                            clusters,species,citation,technology,
                                            do.main.types = do.main.types,
                                            variable.genes=variable.genes,
                                            fine.tune=fine.tune,numCores = numCores))
        }
    } else {
        res = lapply(ref.list, FUN=function(x) {
            SingleR.CreateObject(sc.data.gl,x,clusters,species,citation,technology,
                                 do.main.types=do.main.types,
                                 variable.genes=variable.genes,fine.tune=fine.tune,
                                 numCores = numCores)
        })
    }
    
    singler$singler = res
    
    if (do.signatures==TRUE) {
        signatures = calculateSingScores(sc.data.gl,species=species)
        singler$signatures = signatures
        
    }
    
    if (species == 'Human') {
        kang = SingleR.CreateKangAnnotations(sc.data.gl)
        singler$other = kang$kang_annotation
    }
    
    singler$meta.data = list(project.name=project.name,orig.ident=orig.ident)
    
    if (reduce.file.size==T) {
        singler = remove.Unnecessary.Data.single(singler)
    }
    
    singler
    
}

library(SingleR)
load("/liulab/asahu/data/ssgsea/xiaoman/getz/icb.expression.matched.RData")
singler = mySCObject(counts=t(icb.expression.matched), 
                     project.name = 'GSE120575', min.genes = 100,
                     technology= 'Smart-seq', species = 'Human', citation = 'Moshe et al. 2018', reduce.file.size = T, variable.genes = 'de', normalize.gene.length = F)
save(file="/liulab/asahu/data/ssgsea/xiaoman/getz/singler.RData", singler)

