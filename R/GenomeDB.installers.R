# ASSEMBLIES[hs37d5]="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
# ASSEMBLIES[hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz"
# ASSEMBLIES[GRCh37]="ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
# ASSEMBLIES[hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromFa.tar.gz"
# ASSEMBLIES[GRCh38]="ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
# ASSEMBLIES[mm10]="http://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/chromFa.tar.gz"
# ASSEMBLIES[GRCm38]="ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
#
# declare -A ANNOTATIONS
# ANNOTATIONS[GENCODE19]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
# ANNOTATIONS[RefSeq_hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz"
# ANNOTATIONS[ENSEMBL87]="ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
# ANNOTATIONS[GENCODE28]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
# ANNOTATIONS[RefSeq_hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz"
# ANNOTATIONS[ENSEMBL93]="ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.chr.gtf.gz"
# ANNOTATIONS[ENSEMBL104]="ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz"
# ANNOTATIONS[GENCODEM25]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
# ANNOTATIONS[RefSeq_mm10]="http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/refGene.txt.gz"

# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
###dependencies
library(R.utils)

#' SetDBdirectory
#' Set GenomeDB directory for the RSAP family.
#' It will create a directory .../main/GenomesDB, where main is "my/favorite/directory"
#' @param main directory path name
#' @export
#' @usage SetDBdirectory("my/preferred/path")
#'
#
SetDBdirectory  <- function(main){
  software <- .OpenConfigFile()
  if(length(software)==0){##no hay directorio previo
    software$main <- main
  }
  software$GenomesDB$main <- file.path(software$main,"GenomesDB")
  if(dir.create(software$GenomesDB$main)==FALSE){
    stop("FAIL")
  }else{
    message(paste0("\n",software$GenomesDB$main, " created"))
  }
  .OpenConfigFile(software)
}

.removeConfigFile <- function(){
  file.remove(file.path( .libPaths()[1],"GenomeDB.ElmerLab.RData"))
}
#' .OpenConfigFile
#' Internal. save and restore config info
#' @param config config list structure with software$softName...
#'
.OpenConfigFile <- function(config){
  config.file <- file.path( .libPaths()[1],"GenomeDB.ElmerLab.RData")
  if(missing(config)){##si lo llama vacio, devuelve lo que hay en el file
    if(file.exists(config.file)){
      config <- readRDS(config.file)
      return(invisible(config))
    }else{
      return(list())
    }
    }else{#si config no esta vacio, entonces hay que guardar
      saveRDS(config,config.file)
      return(invisible(config))
    }



}


#' AddGenome
#' This function will install a Genome into your database
#'@param species string
#'@param urlfasta string, the url of the genome assembly
#'@param urlGTF string, the url of the genome GTF annotation file
#'@param version string denoting the version name to identify the assembly+annotation pair
#'@export
AddGenome <- function(species, urlFasta, urlGTF, version){


  software <- .OpenConfigFile()

    # assembly.destination <- "/media/respaldo4t/Assemblies"
  # annotation.destination <- "/media/respaldo4t/Annotation"

  ##aca verificamos si existe la version de la especie
  if(any(stringr::str_detect(software$GenomesDB,species))==FALSE){
    message(paste0("\nThe ",species, " specie does not exist in your Genomes Database"))
  }
  ##aqui agrego la especie a la lista, luego ya la puedo agregar con $
  software$GenomesDB[[species]]$main <- file.path(software$GenomesDB,species)
  #-------------------------
  if(dir.exists(software$GenomesDB[[species]]$main)==FALSE){
    if(dir.create(software$GenomesDB[[species]]$main)==FALSE){
      stop("FAIL create")
    }
  }else{
    message(paste0("\n",software$GenomesDB[[species]]$main," already exists"))
  }


  if(dir.exists(file.path(software$GenomesDB[[species]],version))==FALSE){
    if(is.null(software$GenomesDB[[species]]$version)){
      software$GenomesDB[[species]]$version <- file.path(software$GenomesDB[[species]],version)
      names(software$GenomesDB[[species]]$version) <- version
    }else{
      onames <- names(software$GenomesDB[[species]]$version)
      software$GenomesDB[[species]]$version <- c(software$GenomesDB[[species]]$version, file.path(software$GenomesDB[[species]],version) )
      names(software$GenomesDB[[species]]$version) <- c(onames, version)
    }

    dir.create(software$GenomesDB[[species]]$version[version])
  }else{
    #la version ya existe
    stop(paste0("\nThe " , species,":",version, " already exists, pls use a different version"))
  }

## we will identify the species by specie and version
  download.file(url = urlFasta,
                method = "wget",
                destfile = file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlFasta))),
                extra = "--quiet"
  )
  if(file.exists(file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlFasta))))){
    if(stringr::str_detect(file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlFasta))), ".gz")){
      R.utils::gunzip(filename = file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlFasta))),
                      remove = T)
    }
    cat(paste0("\n",paste0(version,"_",basename(urlFasta))," Installed"))
  }else{
    stop("FAIL genome instalation")
  }
  # R.utils::gunzip(filename = file.path(software$assembly,basename(urlFasta)),
  #       remove = T)
  download.file(url = urlGTF,
                method = "wget",
                destfile = file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlGTF))),
                extra = "--quiet"
  )
  if(file.exists(file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlGTF))))){
    if(stringr::str_detect(file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlGTF))), ".gz")){
      R.utils::gunzip(filename = file.path(software$GenomesDB[[species]]$version[version],paste0(version,"_",basename(urlGTF))),
                      remove = T)
    }
    cat(paste0("\n",paste0(version,"_",basename(urlGTF))," Installed"))
  }
  .OpenConfigFile(software)

}

#' AddHumanGenomes
#' @param urlFasta url for the genome fasta file
#' @param urlGTF the url of the annotation file
#' @param version the name of the version identifying assembly+annotation sources
#' @usage if any is missing, the GRCh38 assembly and the GENCODE release 38 annotation GTF is installed with version name "GRCh38+GENCODE". Else
#' all the parameters are required
#' @export
AddHumanGenomes <- function(urlFasta, urlGTF,version){
  if(any(c(missing(urlFasta),missing(urlGTF),missing(version)))==TRUE){
    message("\n adding Human Genome GRCh38 release 104")
    AddGenome(species = "Human",
              urlFasta = "ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
              urlGTF ="ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz",
              version = "GRCh38+GENECODE")
    gen.files <- GetGenome("GRCh38+GENECODE")

  }else{
    AddGenome(species = "Human",
              urlFasta = urlFasta,
              urlGTF = urlGTF,
              version = version)
  }

}

#' GetGenome
#' return a list with the slots fasta and gtf with the full path to the genomes and annotation files
#' @param specie the stored specie (Human" etc..)
#' @param version "the assembly+annotation version/code
#' @export
GetGenome <- function(species, version){
  databases <- GenomeDB:::.OpenConfigFile()
  if(all(stringr::str_detect(names(databases$GenomesDB[[species]]$version),version))==FALSE){

    stop(paste0("\nThe ", specie, " and ", version, " is not in the Database\nAvailable versions are: ", names(database$GenomesDB[[species]]$version)))
  }
  if(dir.exists(file.path(databases$GenomesDB[[species]]$main,version))==FALSE){
    stop(paste0("\nGenome for ",species, " version ", version, " NOT FOUND"))
  }
  fasta.gtf.files <- list.files(file.path(databases$GenomesDB[[species]]$main,version), full.names = TRUE)
  fasta <- fasta.gtf.files[stringr::str_detect(fasta.gtf.files,".fasta|.fa")]
  gtf <- fasta.gtf.files[stringr::str_detect(fasta.gtf.files,".gtf")]
  return(list(fasta=fasta,gtf=gtf))
}

