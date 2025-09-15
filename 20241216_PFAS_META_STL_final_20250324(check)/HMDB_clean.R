
library(XML)
library("methods")
# # 
setwd("D:/OneDrive/Postdoctoral_SCI/Sperm_bank/hmdb_metabolites")
file <- dir()
hmdbLib <- vector(mode = "list", length = length(file))
for(i in 1:length(file)){
  cat(i, " ")
  temp.file <- file[i]
  result <- xmlParse(file = temp.file)
  rootnode <- xmlRoot(result)
  info <- xmlSApply(rootnode, function(x) xmlSApply(x, xmlValue))
  names(info) <- unname(names(rootnode))
  info <- lapply(info, unname)
  ms2.mz <- sapply(getNodeSet(doc = rootnode[[26]], path = "//mass-charge"), xmlValue)
  ms2.int <- sapply(getNodeSet(doc = rootnode[[26]], path = "//intensity"), xmlValue)
  ms2 <- data.frame("mz" = ms2.mz, "intensity" = ms2.int,
                    stringsAsFactors = FALSE)
  info[[26]] <- ms2
  
  ms2.info <- info[c(10, 17, 18, 24, 25)]
  ms2.info <- lapply(ms2.info, function(x){
    if(length(x) == 0){
      x <- NA
    }
    x
  })
  ms2.spec <- info[[26]]
  
  ms2.info <- do.call(rbind, ms2.info)
  ms2.info <- data.frame("name" = rownames(ms2.info), "value" = ms2.info[,1],
                         stringsAsFactors = FALSE)
  ms2.info <- tibble::as_tibble(ms2.info)
  
  temp.spec <- list(ms2.info = ms2.info, ms2.spec = ms2.spec)
  hmdbLib[[i]] <- temp.spec
  
}

hmdb <- hmdbLib




convert_hmdb2metid <-
  function(data,
           path = ".",
           threads = 5,
           ms1_or_ms2 = c("ms1", "ms2")) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
    ms1_or_ms2 <- match.arg(ms1_or_ms2)
    
    if (ms1_or_ms2 == "ms1") {
      data <-
        data %>%
        dplyr::filter(!is.na(monisotopic_molecular_weight))
      
      data[which(data == "", arr.ind = TRUE)] <- NA
      
      data <-
        data %>%
        dplyr::rename(
          mz = average_molecular_weight,
          Create_date = creation_date,
          Updated_date = update_date,
          Lab.ID = accession,
          Compound.name = name,
          Description = description,
          Synonyms = synonyms,
          Formula = chemical_formula,
          IUPAC_name = iupac_name,
          Traditional_IUPAC_name = traditional_iupac,
          CAS.ID = cas_registry_number,
          SMILES.ID = smiles,
          INCHI.ID = inchi,
          INCHIKEY.ID = inchikey,
          Kingdom = kingdom,
          Super_class = super_class,
          Class = class,
          Sub_class = sub_class,
          State = state,
          Biospecimen_locations = biospecimen_locations,
          Cellular_locations = cellular_locations,
          Tissue_locations = tissue_locations,
          CHEMSPIDER.ID = chemspider_id,
          DRUGBANK.ID = drugbank_id,
          FOODB.ID = foodb_id,
          PUBCHEM.ID = pubchem_compound_id,
          CHEBI.ID = chebi_id,
          KEGG.ID = kegg_id,
          BIOCYC.ID = biocyc_id,
          BIGG.ID = bigg_id,
          WIKIPEDIA.ID = wikipedia_id,
          METLIN.ID = metlin_id
        ) %>%
        dplyr::mutate(
          HMDB.ID = Lab.ID,
          RT = NA,
          mz.pos = NA,
          mz.neg = NA,
          Submitter = "HMDB",
          From_human = "Yes"
        ) %>%
        dplyr::select(
          Lab.ID,
          Compound.name,
          mz,
          RT,
          CAS.ID,
          HMDB.ID,
          KEGG.ID,
          Formula,
          mz.pos,
          mz.neg,
          Submitter,
          everything()
        )
      
      temp_file <- tempfile()
      dir.create(temp_file, showWarnings = FALSE)
      
      readr::write_csv(x = data,
                       file = file.path(temp_file, "data.csv"))
      
      hmdb_ms1 <-
        metid::construct_database(
          path = temp_file,
          version =  as.character(Sys.Date()),
          metabolite.info.name = "data.csv",
          source = "HMDB",
          link = "https://hmdb.ca/",
          creater = "Xiaotao Shen",
          email = "shenxt@stanford.edu",
          rt = FALSE,
          threads = threads
        )
      
      save(hmdb_ms1, file = file.path(path, "hmdb_ms1"))
      invisible(hmdb_ms1)
    } else{
      remove_idx <-
        data %>%
        lapply(function(x) {
          nrow(x$ms2)
        }) %>%
        unlist() %>%
        `==`(0) %>%
        which()
      
      if (length(remove_idx) > 0) {
        data <-
          data[-remove_idx]
      }
      
      spectra_info <-
        data %>%
        purrr::map(function(x) {
          x$ms1_info
        }) %>%
        dplyr::bind_rows() %>%
        as.data.frame()
      
      spectra_data <-
        data %>%
        purrr::map(function(x) {
          x$ms2
        })
      
      spectra_info[which(spectra_info == "NA", arr.ind = TRUE)] <-
        NA
      spectra_info[which(spectra_info == "n/a", arr.ind = TRUE)] <-
        NA
      spectra_info[which(spectra_info == "N/A", arr.ind = TRUE)] <-
        NA
      
      spectra_info <-
        spectra_info %>%
        dplyr::select(HMDB.ID,
                      Instrument_type,
                      Polarity,
                      collision_energy_voltage,
                      adduct)
      
      remove_idx <-
        which(is.na(spectra_info$Polarity))
      
      if (length(remove_idx) > 0) {
        spectra_info <-
          spectra_info[-remove_idx, ]
        
        spectra_data <-
          spectra_data[-remove_idx]
      }
      
      spectra_info <-
        spectra_info %>%
        dplyr::mutate(
          Polarity = case_when(
            Polarity == "positive" ~ "Positive",
            Polarity == "negative" ~ "Negative",
            TRUE ~ Polarity
          )
        )
      
      spectra_info$Lab.ID <-
        masstools::name_duplicated(spectra_info$HMDB.ID) %>%
        paste("shen", sep = "_")
      
      spectra_info2 <-
        spectra_info %>%
        plyr::dlply(.variables = .(HMDB.ID)) %>%
        purrr::map(function(y) {
          if (sum(is.na(y$collision_energy_voltage)) > 0) {
            y$collision_energy_voltage[is.na(y$collision_energy_voltage)] <-
              paste("Unknown",
                    1:length(y$collision_energy_voltage[is.na(y$collision_energy_voltage)]),
                    sep = "_")
          }
          y
        }) %>%
        dplyr::bind_rows() %>%
        as.data.frame()
      
      spectra_info2 <-
        spectra_info2[match(spectra_info$Lab.ID, spectra_info2$Lab.ID), ]
      
      spectra_data2 <-
        1:length(spectra_data) %>%
        purrr::map(function(i) {
          x <- spectra_data[[i]]
          x <- list(x)
          names(x) <-
            spectra_info2$collision_energy_voltage[i]
          x
        })
      
      names(spectra_data2) <- spectra_info2$Lab.ID
      
      ######positive mode
      spectra_info2$Lab.ID == names(spectra_data2)
      
      index_pos <- which(spectra_info2$Polarity == "Positive")
      index_neg <- which(spectra_info2$Polarity == "Negative")
      
      spectra_info_pos <- spectra_info2[index_pos, ]
      spectra_data_pos <- spectra_data2[index_pos]
      
      spectra_info_neg <- spectra_info2[index_neg, ]
      spectra_data_neg <- spectra_data2[index_neg]
      
      colnames(spectra_info2)
      colnames(hmdb_ms1@spectra.info)
      
      spectra_info2 <-
        spectra_info2 %>%
        dplyr::rename(CE = "collision_energy_voltage")
      
      spectra_info2 <-
        spectra_info2 %>%
        dplyr::left_join(hmdb_ms1@spectra.info %>% dplyr::select(-Lab.ID),
                         by = c("HMDB.ID"))
      
      temp_file <- tempfile()
      dir.create(temp_file, showWarnings = FALSE)
      
      readr::write_csv(x = spectra_info2,
                       file = file.path(temp_file, "spectra_info2.csv"))
      
      hmdb_ms2 <-
        metid::construct_database(
          path = temp_file,
          version =  as.character(Sys.Date()),
          metabolite.info.name = "spectra_info2.csv",
          source = "HMDB",
          link = "https://hmdb.ca/",
          creater = "Xiaotao Shen",
          email = "shenxt@stanford.edu",
          rt = FALSE,
          threads = threads
        )
      
      hmdb_ms2@spectra.data$Spectra.positive <- spectra_data_pos
      hmdb_ms2@spectra.data$Spectra.negative <- spectra_data_neg
      
      save(hmdb_ms2, file = file.path(path, "hmdb_ms2"))
      message("Done.")
      invisible(hmdb_ms2)
    }
    
  }







temp_file <- tempfile()
dir.create(temp_file, showWarnings = FALSE)

readr::write_csv(x = data,
                 file = file.path(temp_file, "data.csv"))

hmdb_ms1 <-
  metid::construct_database(
    path = temp_file,
    version =  as.character(Sys.Date()),
    metabolite.info.name = "data.csv",
    source = "HMDB",
    link = "https://hmdb.ca/",
    # creater = "Xiaotao Shen",
    # email = "shenxt@stanford.edu",
    rt = FALSE,
    threads = threads
  )

save(hmdb_ms1, file = file.path(path, "hmdb_ms1"))
invisible(hmdb_ms1) 





setClass(
  Class = "pathway_database",
  representation(
    database_info = "list",
    pathway_id = "vector",
    pathway_name = "vector",
    describtion = "list",
    pathway_class = "list",
    gene_list = "list",
    compound_list = "list",
    protein_list = "list",
    reference_list = "list",
    related_disease = "list",
    related_module = "list"
  ),
  prototype = list(
    database_info = list(),
    pathway_id = list(),
    pathway_name = list(),
    describtion = list(),
    pathway_class = list(),
    gene_list = list(),
    compound_list = list(),
    protein_list = list(),
    reference_list = list(),
    related_disease = list(),
    related_module = list()
  )
)


setMethod(
  f = "show",
  signature = "pathway_database",
  definition = function(object) {
    version <- try(object@database_info$version, silent = TRUE)
    source <- try(object@database_info$source, silent = TRUE)
    if (!is(version, "try-error")) {
      message(crayon::green("---------Pathway source&version---------"))
      message(crayon::green(source, " & ", version))
    }
    message(crayon::green("-----------Pathway information------------"))
    message(crayon::green(length(object@pathway_id), " pathways"))
    message(
      crayon::green(
        object@gene_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        " pathways have genes"
      )
    )
    
    message(
      crayon::green(
        object@protein_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        " pathways have proteins"
      )
    )
    
    message(
      crayon::green(
        object@compound_list %>%
          lapply(nrow) %>%
          unlist() %>%
          `!=`(0) %>%
          sum(),
        " pathways have compounds"
      )
    )
    
    message(crayon::green("Pathway class (top 10):",
                          paste(unique(head(
                            unlist(object@pathway_class), 10
                          )), collapse = ";")))
  }
)


library(massdatabase)
library(metpath)

get_hmdb_pathway() -> hmdb_pathway



get_pathway_class(hmdb_pathway)
#get the class of pathways
pathway_class =
  metpath::pathway_class(hmdb_pathway)

remain_idx = which(unlist(pathway_class) == "Metabolic;primary_pathway")

hmdb_pathway =
  hmdb_pathway[remain_idx]

hmdb_pathway@pathway_name


download_smpdb_pathway(path = ".")

# convert_smpdb2metpath()

data <- read_smpdb_pathway(path = ".", only_primarity_pathway = TRUE)

smpdb_pathway_database <- convert_smpdb2metpath(data = data, path = ".")





if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('biodbHmdb')

# 加载包
library(biodbHmdb)

# 获取所有代谢通路
pathways <- getPathways()

# 查看结果
print(pathways)

browseVignettes("biodbHmdb")   ###查看HMDB文档

mybiodb <- biodb::newInst()


conn <- mybiodb$getFactory()$createConn('hmdb.metabolites')

dbExtract <- system.file("extdata", 'generated', "hmdb_extract.zip",
                         package="biodbHmdb")
conn$setPropValSlot('urls', 'db.zip.url', dbExtract)

conn$getNbEntries()

mybiodb$terminate()








library(httr)

HMDBDownloader <- list(
  BASE_URL = "https://hmdb.ca/unearth/q",
  
  download_metabolite = function(metabolite_id) {
    # 构建请求的 URL
    url <- paste0(HMDBDownloader$BASE_URL, "/", metabolite_id)
    response <- GET(url)
    
    # 检查响应状态
    if (status_code(response) == 200) {
      content(response, "text", encoding = "UTF-8")
    } else {
      stop(paste("Failed to fetch data for metabolite ID", metabolite_id))
    }
  },
  
  save_metabolite_data = function(metabolite_id, file_path) {
    # 下载数据
    data <- HMDBDownloader$download_metabolite(metabolite_id)
    
    # 将数据写入文件
    writeLines(data, file_path)
  }
)



# HMDBDownloader$download_metabolite

HMDBDownloader$save_metabolite_data("HMDB00001", "metabolite_HMDB00001.txt")




# 
# # 加载必要的库
# library(httr)
# library(xml2)
# library(dplyr)
# library(stringr)
# 
# # 1. 批量爬取 HMDB 每页的网页
# base_url <- "https://hmdb.ca/hml/metabolites?page="
# for (i in 1:46) {
#   url <- paste0(base_url, i)
#   response <- GET(url)
#   writeBin(content(response, "raw"), paste0("page", i, ".html"))
#   Sys.sleep(15)  # 防止请求过快
#   cat("Downloaded page", i, "\n")
# }
# 
# # 2. 解析本地的 HTML 文件，提取代谢物的 ID 和名称
# metabolite_dict <- data.frame(id = character(), name = character(), stringsAsFactors = FALSE)
# 
# for (i in 1:46) {
#   html_content <- readLines(paste0("page", i, ".html"), warn = FALSE)
#   metabolite_lines <- grep("class=\"metabolite-name\">", html_content, value = TRUE)
#   
#   for (line in metabolite_lines) {
#     id <- str_extract(line, "HMDB[0-9]{7}")
#     name <- str_extract(line, "(?<=<strong>)[^<]+")
#     if (!is.na(id) && !is.na(name)) {
#       metabolite_dict <- rbind(metabolite_dict, data.frame(id = id, name = name, stringsAsFactors = FALSE))
#     }
#   }
#   cat("Processed page", i, "\n")
# }
# 
# # 保存代谢物 ID 和名称的字典
# write.table(metabolite_dict, "metabolite_dict.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# 
# getwd()
# 
# # 3. 批量爬取每个代谢物的详细信息（XML 文件）
# base_metabolite_url <- "https://hmdb.ca/metabolites/"
# 
# for (i in 1:nrow(metabolite_dict)) {
#   id <- metabolite_dict$id[i]
#   url <- paste0(base_metabolite_url, id, ".xml")
#   response <- GET(url)
#   writeBin(content(response, "raw"), paste0(id, ".xml"))
#   Sys.sleep(10)  # 防止请求过快
#   cat("Downloaded XML for", id, "\n")
# }
# 
# # 4. 解析 XML 文件，提取代谢物-通路信息
# meta2pathway <- data.frame(
#   HMDB_ID = character(),
#   HMDB_Name = character(),
#   Pathway_Name = character(),
#   SMPDB_ID = character(),
#   KEGG_Map_ID = character(),
#   stringsAsFactors = FALSE
# )
# 
# for (i in 1:nrow(metabolite_dict)) {
#   id <- metabolite_dict$id[i]
#   name <- metabolite_dict$name[i]
#   xml_file <- paste0(id, ".xml")
#   
#   if (file.exists(xml_file)) {
#     xml_content <- read_xml(xml_file)
#     pathway_names <- xml_text(xml_find_all(xml_content, "//pathway/name"))
#     smpdb_ids <- xml_text(xml_find_all(xml_content, "//pathway/smpdb_id"))
#     kegg_map_ids <- xml_text(xml_find_all(xml_content, "//pathway/kegg_map_id"))
#     
#     if (length(pathway_names) > 0) {
#       for (j in 1:length(pathway_names)) {
#         meta2pathway <- rbind(meta2pathway, data.frame(
#           HMDB_ID = id,
#           HMDB_Name = name,
#           Pathway_Name = pathway_names[j],
#           SMPDB_ID = smpdb_ids[j],
#           KEGG_Map_ID = kegg_map_ids[j],
#           stringsAsFactors = FALSE
#         ))
#       }
#       cat("Processed XML for", id, "\n")
#     } else {
#       cat("Skipped XML for", id, "(no pathway data)\n")
#     }
#   } else {
#     cat("Skipped XML for", id, "(file not found)\n")
#   }
# }
# 
# # 保存代谢物-通路信息
# write.table(meta2pathway, "Meta2Pathway.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# 

