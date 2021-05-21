##Aligners
UpdateSTARindex <- function(){
  databases <- GenomeDB:::.OpenConfigFile()
  species <- names(databases$GenomesDB)[-1]

  software <- ArribaR:::.OpenConfigFile()##debera cambiar, solo para usar la version instalada como prueba
  if(is.null(software$star$path)){
      stop("STAR not installed")
    }
  cat("\nBuilding STAR 2.7.6a - ARRIBA Genome Index...This may take some time")
  thr <- ifelse(parallel::detectCores() > 3, parallel::detectCores()-1, parallel::detectCores())
  cat(paste0("\nRunning STAR genomeGenerate with ",thr," CPU cores"))
  for(sp in species){

    versiones <- basename(databases$GenomesDB[[sp]]$version)
    for(vv in versiones){
      if(dir.exists(file.path(software$star$path,versiones))==FALSE){
        if(dir.create(file.path(software$star$path,versiones))==FALSE){
          stop("error creating index directory")
        }
      }
        ##se creo el directorio
        software$path[[versiones]]$main <- file.path(software$star$path,versiones)
        fasta.gtf.files <- GenomeDB::GetGenome(sp, versiones)
        genome.fasta <- file.path(databases$GenomesDB[[sp]],paste0(versiones,"/"))
        system2(command = software$star$command,
                args = c("--runMode genomeGenerate",
                         paste0("--genomeDir ",software$path[[versiones]]$main),
                         paste0("--genomeFastaFiles ",fasta.gtf.files$fasta),
                         paste0("--sjdbGTFfile ",fasta.gtf.files$gtf),
                         paste0("--runThreadN ",thr),
                         "--sjdbOverhang 250") )
      }
    }
    # .OpenConfigFile(software)

}


UpdateRsubreadindex <- function(){
  databases <- GenomeDB:::.OpenConfigFile()
  species <- names(databases$GenomesDB)[-1]

  software <- ArribaR:::.OpenConfigFile()##debera cambiar, solo para usar la version instalada como prueba

  require(Rsubread)

  thr <- ifelse(parallel::detectCores() > 3, parallel::detectCores()-1, parallel::detectCores())
  cat(paste0("\nRunning Rsubread genomeGenerate with ",thr," CPU cores"))
  if(dir.exists(file.path(software$path,"Rsubread"))==FALSE){
    dir.create(file.path(software$path,"Rsubread"))
    software$Rsubread$main <- file.path(software$path,"Rsubread")
  }

  for(sp in species){

    versiones <- basename(databases$GenomesDB[[sp]]$version)
    for(vv in versiones){
      if(dir.exists(file.path(file.path(software$Rsubread$main,versiones))==FALSE){
        if(dir.create(file.path(file.path(software$Rsubread$main,versiones))==FALSE){
          stop("error creating Rsubread index directory")
        }
      }
      ##se creo el directorio
      software$Rsubread$main[[versiones]] <- file.path(file.path(software$Rsubread$main,versiones)
      fasta.gtf.files <- GenomeDB::GetGenome(sp, versiones)
      genome.fasta <- file.path(databases$GenomesDB[[sp]],paste0(versiones,"/"))
      Rsubread::buildindex(basename = software$Rsubread$main[[versiones]],
                           reference = fasta.gtf.files$fasta)

    }
  }
  # .OpenConfigFile(software)

}
