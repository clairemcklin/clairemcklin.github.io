library(neonUtilities)
library(dplyr)
library(readr)
library(tidyverse)
library(vegan)
library(stringr)
library(ggplot2) 
library(data.table)

# load in microbe data from 19 sites
microbes <- loadByProduct(dpID = "DP1.10081.001", 
                          site = c("GUAN","HEAL","HARV","CPER","TALL","SCBI",
                                   "OSBS","UNDE","UKFS","GRSM","WOOD","OAES",
                                   "YELL","NIWO","SRER","ONAQ","ABBY","SJER",
                                   "BARR"),
                          release = "RELEASE-2024", 
                          check.size = F)

# accesses the 16S (bacteria and archaea) data from microbe NEON dataset and 
##loads it as a data table
table1 <- as.data.frame(microbes$mcc_soilSeqVariantMetadata_16S)

# accesses ITS (fungi) data from microbe NEON dataset as a data table
table2 <- as.data.frame(microbes$mcc_soilSeqVariantMetadata_ITS)

# excludes rows containing "below threshold" from 16S data table
filtered_table1 <- table1[table1$sequenceCountQF == "OK", ]

# merging data files keeping plotIDs and siteIDs
microbe_links <- filtered_table1$downloadFileUrl # creates object with all links to CSV files

read_csv_from_url <- function(url, plotID, siteID) {
  temp_file <- tempfile()  # create a temporary file
  download.file(url, temp_file, quiet = TRUE)  # downloads the file
  data <- read.csv(temp_file)  # reads the CSV
  data$plotID <- plotID  # adds plotID column to the dataset
  data$siteID <- siteID # adds siteID column to dataset
  return(data)  # returns the modified dataset
}

# processes each link, downloads the contents of CSV file, and merges all files 
# into one data table
merged_data <- bind_rows(mapply(read_csv_from_url, 
                                microbe_links, 
                                filtered_table1$plotID, 
                                filtered_table1$siteID, 
                                SIMPLIFY = FALSE))

head(merged_data)

#excludes rows containing "below threshold" from ITS data table
filtered_table2 <- table2[table1$sequenceCountQF == "OK", ]

#merge data files keeping plotIDs and siteIDs
microbe_links2 <- filtered_table2$downloadFileUrl # creates object with all links to CSV files

read_csv_from_url2 <- function(url, plotID, siteID) {
  temp_file <- tempfile()  # create a temporary file
  download.file(url, temp_file, quiet = TRUE)  # downloads the file
  data <- read.csv(temp_file)  # reads the CSV
  data$plotID <- plotID  # adds plotID column to the dataset
  data$siteID <- siteID # adds siteID column to dataset
  return(data)  # returns the modified dataset
}

# processes each link, downloads the contents of CSV file, and merges all files 
# into one data table
merged_data2 <- bind_rows(mapply(read_csv_from_url2, 
                                 microbe_links2, 
                                 filtered_table2$plotID, 
                                 filtered_table2$siteID, 
                                 SIMPLIFY = FALSE))

head(merged_data2)

# combines 16S data table and ITS data table
microbes_total <- bind_rows(merged_data,merged_data2)

# make separate data tables for each site
microbe_GUAN <- filter(microbes_total, siteID=="GUAN")
microbe_HEAL <- filter(microbes_total, siteID=="HEAL")
microbe_HARV <- filter(microbes_total, siteID=="HARV")
microbe_CPER <- filter(microbes_total, siteID=="CPER")
microbe_TALL <- filter(microbes_total, siteID=="TALL")
microbe_SCBI <- filter(microbes_total, siteID=="SCBI")
microbe_OSBS <- filter(microbes_total, siteID=="OSBS")
microbe_UNDE <- filter(microbes_total, siteID=="UNDE")
microbe_UKFS <- filter(microbes_total, siteID=="UKFS")
microbe_GRSM <- filter(microbes_total, siteID=="GRSM")
microbe_WOOD <- filter(microbes_total, siteID=="WOOD")
microbe_OAES <- filter(microbes_total, siteID=="OAES")
microbe_YELL <- filter(microbes_total, siteID=="YELL")
microbe_NIWO <- filter(microbes_total, siteID=="NIWO")
microbe_SRER <- filter(microbes_total, siteID=="SRER")
microbe_ONAQ <- filter(microbes_total, siteID=="ONAQ")
microbe_ABBY <- filter(microbes_total, siteID=="ABBY")
microbe_SJER <- filter(microbes_total, siteID=="SJER")
microbe_BARR <- filter(microbes_total, siteID=="BARR")


## ## make data tables for each site indicating the number of species per plot, make NA 
### values equal 0 , then remove the plot names in order to calculate SHD and richness
microbe_matrix_GUAN <- microbe_GUAN %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_GUAN[is.na(microbe_matrix_GUAN)] <- 0
microbe_matrix_GUAN <- microbe_matrix_GUAN[,-1]

microbe_matrix_HEAL <- microbe_HEAL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_HEAL[is.na(microbe_matrix_HEAL)] <- 0
microbe_matrix_HEAL <- microbe_matrix_HEAL[,-1]

microbe_matrix_HARV <- microbe_HARV %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_HARV[is.na(microbe_matrix_HARV)] <- 0
microbe_matrix_HARV <- microbe_matrix_HARV[,-1]

microbe_matrix_CPER <- microbe_CPER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_CPER[is.na(microbe_matrix_CPER)] <- 0
microbe_matrix_CPER <- microbe_matrix_CPER[,-1]

microbe_matrix_TALL <- microbe_TALL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_TALL[is.na(microbe_matrix_TALL)] <- 0
microbe_matrix_TALL <- microbe_matrix_TALL[,-1]

microbe_matrix_SCBI <- microbe_SCBI %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_SCBI[is.na(microbe_matrix_SCBI)] <- 0
microbe_matrix_SCBI <- microbe_matrix_SCBI[,-1]

microbe_matrix_OSBS <- microbe_OSBS %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_OSBS[is.na(microbe_matrix_OSBS)] <- 0
microbe_matrix_OSBS <- microbe_matrix_OSBS[,-1]

microbe_matrix_UNDE <- microbe_UNDE %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_UNDE[is.na(microbe_matrix_UNDE)] <- 0
microbe_matrix_UNDE <- microbe_matrix_UNDE[,-1]

microbe_matrix_UKFS <- microbe_UKFS %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_UKFS[is.na(microbe_matrix_UKFS)] <- 0
microbe_matrix_UKFS <- microbe_matrix_UKFS[,-1]

microbe_matrix_GRSM <- microbe_GRSM %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_GRSM[is.na(microbe_matrix_GRSM)] <- 0
microbe_matrix_GRSM <- microbe_matrix_GRSM[,-1]

microbe_matrix_WOOD <- microbe_WOOD %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_WOOD[is.na(microbe_matrix_WOOD)] <- 0
microbe_matrix_WOOD <- microbe_matrix_WOOD[,-1]

microbe_matrix_OAES <- microbe_OAES %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_OAES[is.na(microbe_matrix_OAES)] <- 0
microbe_matrix_OAES <- microbe_matrix_OAES[,-1]

microbe_matrix_YELL <- microbe_YELL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_YELL[is.na(microbe_matrix_YELL)] <- 0
microbe_matrix_YELL <- microbe_matrix_YELL[,-1]

microbe_matrix_NIWO <- microbe_NIWO %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_NIWO[is.na(microbe_matrix_NIWO)] <- 0
microbe_matrix_NIWO <- microbe_matrix_NIWO[,-1]

microbe_matrix_SRER <- microbe_SRER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_SRER[is.na(microbe_matrix_SRER)] <- 0
microbe_matrix_SRER <- microbe_matrix_SRER[,-1]

microbe_matrix_ONAQ <- microbe_ONAQ %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_ONAQ[is.na(microbe_matrix_ONAQ)] <- 0
microbe_matrix_ONAQ <- microbe_matrix_ONAQ[,-1]

microbe_matrix_ABBY <- microbe_ABBY %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_ABBY[is.na(microbe_matrix_ABBY)] <- 0
microbe_matrix_ABBY <- microbe_matrix_ABBY[,-1]

microbe_matrix_SJER <- microbe_SJER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_SJER[is.na(microbe_matrix_SJER)] <- 0
microbe_matrix_SJER <- microbe_matrix_SJER[,-1]

microbe_matrix_BARR <- microbe_BARR %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
microbe_matrix_BARR[is.na(microbe_matrix_BARR)] <- 0
microbe_matrix_BARR <- microbe_matrix_BARR[,-1]

##calculating shannon diversity and richness for each site
shd_GUANmicrobes <- diversity(microbe_matrix_GUAN, index = "shannon", MARGIN = 1, base = exp(1))
shd_HEALmicrobes <- diversity(microbe_matrix_HEAL, index = "shannon", MARGIN = 1, base = exp(1))
shd_HARVmicrobes <- diversity(microbe_matrix_HARV, index = "shannon", MARGIN = 1, base = exp(1))
shd_CPERmicrobes <- diversity(microbe_matrix_CPER, index = "shannon", MARGIN = 1, base = exp(1))
shd_TALLmicrobes <- diversity(microbe_matrix_TALL, index = "shannon", MARGIN = 1, base = exp(1))
shd_SCBImicrobes <- diversity(microbe_matrix_SCBI, index = "shannon", MARGIN = 1, base = exp(1))
shd_OSBSmicrobes <- diversity(microbe_matrix_OSBS, index = "shannon", MARGIN = 1, base = exp(1))
shd_UNDEmicrobes <- diversity(microbe_matrix_UNDE, index = "shannon", MARGIN = 1, base = exp(1))
shd_UKFSmicrobes <- diversity(microbe_matrix_UKFS, index = "shannon", MARGIN = 1, base = exp(1))
shd_GRSMmicrobes <- diversity(microbe_matrix_GRSM, index = "shannon", MARGIN = 1, base = exp(1))
shd_WOODmicrobes <- diversity(microbe_matrix_WOOD, index = "shannon", MARGIN = 1, base = exp(1))
shd_OAESmicrobes <- diversity(microbe_matrix_OAES, index = "shannon", MARGIN = 1, base = exp(1))
shd_YELLmicrobes <- diversity(microbe_matrix_YELL, index = "shannon", MARGIN = 1, base = exp(1))
shd_NIWOmicrobes <- diversity(microbe_matrix_NIWO, index = "shannon", MARGIN = 1, base = exp(1))
shd_SRERmicrobes <- diversity(microbe_matrix_SRER, index = "shannon", MARGIN = 1, base = exp(1))
shd_ONAQmicrobes <- diversity(microbe_matrix_ONAQ, index = "shannon", MARGIN = 1, base = exp(1))
shd_ABBYmicrobes <- diversity(microbe_matrix_ABBY, index = "shannon", MARGIN = 1, base = exp(1))
shd_SJERmicrobes <- diversity(microbe_matrix_SJER, index = "shannon", MARGIN = 1, base = exp(1))
shd_BARRmicrobes <- diversity(microbe_matrix_BARR, index = "shannon", MARGIN = 1, base = exp(1))

rich_GUANmicrobes <- specnumber(microbe_matrix_GUAN, MARGIN = 1)
rich_HEALmicrobes <- specnumber(microbe_matrix_HEAL, MARGIN = 1)
rich_HARVmicrobes <- specnumber(microbe_matrix_HARV, MARGIN = 1)
rich_CPERmicrobes <- specnumber(microbe_matrix_CPER, MARGIN = 1)
rich_TALLmicrobes <- specnumber(microbe_matrix_TALL, MARGIN = 1)
rich_SCBImicrobes <- specnumber(microbe_matrix_SCBI, MARGIN = 1)
rich_OSBSmicrobes <- specnumber(microbe_matrix_OSBS, MARGIN = 1)
rich_UNDEmicrobes <- specnumber(microbe_matrix_UNDE, MARGIN = 1)
rich_UKFSmicrobes <- specnumber(microbe_matrix_UKFS, MARGIN = 1)
rich_GRSMmicrobes <- specnumber(microbe_matrix_GRSM, MARGIN = 1)
rich_WOODmicrobes <- specnumber(microbe_matrix_WOOD, MARGIN = 1)
rich_OAESmicrobes <- specnumber(microbe_matrix_OAES, MARGIN = 1)
rich_YELLmicrobes <- specnumber(microbe_matrix_YELL, MARGIN = 1)
rich_NIWOmicrobes <- specnumber(microbe_matrix_NIWO, MARGIN = 1)
rich_SRERmicrobes <- specnumber(microbe_matrix_SRER, MARGIN = 1)
rich_ONAQmicrobes <- specnumber(microbe_matrix_ONAQ, MARGIN = 1)
rich_ABBYmicrobes <- specnumber(microbe_matrix_ABBY, MARGIN = 1)
rich_SJERmicrobes <- specnumber(microbe_matrix_SJER, MARGIN = 1)
rich_BARRmicrobes <- specnumber(microbe_matrix_BARR, MARGIN = 1)


### make data tables for each site indicating the number of species per plot, make NA 
### values equal 0 , then remove the plot names in order to calculate SHD and richness
plotID <- c("GUAN_001","GUAN_002","GUAN_003","GUAN_004","GUAN_006","GUAN_007",
            "GUAN_042","GUAN_043","GUAN_048","GUAN_049")
siteID <- rep("GUAN", length(plotID))
species <- rep("microbe", length(plotID))
x_GUAN_microbes <- data.frame(plotID, siteID, shd_GUANmicrobes, rich_GUANmicrobes,
                              species)
colnames(x_GUAN_microbes)<- c("plotID","siteID", "shd", "rich","species")

### HEAL data table with shannon and diversity
plotID <- c("HEAL_001","HEAL_004","HEAL_005","HEAL_010","HEAL_011","HEAL_025",
            "HEAL_045","HEAL_046","HEAL_047","HEAL_048")
siteID <- rep("HEAL", length(plotID))
species <- rep("microbe", length(plotID))
x_HEAL_microbes <- data.frame(plotID, siteID, shd_HEALmicrobes, rich_HEALmicrobes,
                              species)
colnames(x_HEAL_microbes)<- c("plotID","siteID", "shd", "rich","species")

### HARV data table with shannon and diversity
plotID <- c("HARV_001","HARV_002","HARV_004","HARV_005","HARV_010","HARV_013",
            "HARV_016","HARV_020","HARV_021","HARV_033","HARV_034","HARV_035",
            "HARV_037")
siteID <- rep("HARV", length(plotID))
species <- rep("microbe", length(plotID))
x_HARV_microbes <- data.frame(plotID, siteID, shd_HARVmicrobes, rich_HARVmicrobes,
                              species)
colnames(x_HARV_microbes)<- c("plotID","siteID", "shd", "rich","species")

### CPER data table with shannon and diversity
plotID <- c("CPER_001","CPER_002","CPER_003","CPER_004","CPER_005","CPER_006",
            "CPER_045","CPER_046","CPER_047","CPER_048")
siteID <- rep("CPER", length(plotID))
species <- rep("microbe", length(plotID))
x_CPER_microbes <- data.frame(plotID, siteID, shd_CPERmicrobes, rich_CPERmicrobes,
                              species)
colnames(x_CPER_microbes)<- c("plotID","siteID", "shd", "rich","species")

### TALL data table with shannon and diversity
plotID <- c("TALL_001","TALL_002","TALL_003","TALL_004","TALL_006","TALL_007",
            "TALL_044","TALL_047","TALL_051","TALL_054")
siteID <- rep("TALL", length(plotID))
species <- rep("microbe", length(plotID))
x_TALL_microbes <- data.frame(plotID, siteID, shd_TALLmicrobes, rich_TALLmicrobes,
                              species)
colnames(x_TALL_microbes)<- c("plotID","siteID", "shd", "rich","species")

### SCBI data table with shannon and diversity
plotID <- c("SCBI_002","SCBI_003","SCBI_004","SCBI_005","SCBI_006","SCBI_008",
            "SCBI_012","SCBI_045","SCBI_046","SCBI_047","SCBI_049","SCBI_067")
siteID <- rep("SCBI", length(plotID))
species <- rep("microbe", length(plotID))
x_SCBI_microbes <- data.frame(plotID, siteID, shd_SCBImicrobes, rich_SCBImicrobes,
                              species)
colnames(x_SCBI_microbes)<- c("plotID","siteID", "shd", "rich","species")

### OSBS data table with shannon and diversity
plotID <- c("OSBS_001","OSBS_002","OSBS_003","OSBS_004","OSBS_005","OSBS_022",
            "OSBS_023","OSBS_026","OSBS_027","OSBS_029","OSBS_031")
siteID <- rep("OSBS", length(plotID))
species <- rep("microbe", length(plotID))
x_OSBS_microbes <- data.frame(plotID, siteID, shd_OSBSmicrobes, rich_OSBSmicrobes,
                              species)
colnames(x_OSBS_microbes)<- c("plotID","siteID", "shd", "rich","species")

### UNDE data table with shannon and diversity
plotID <- c("UNDE_001","UNDE_002","UNDE_003","UNDE_006","UNDE_007","UNDE_008",
            "UNDE_010","UNDE_013","UNDE_014","UNDE_017","UNDE_019","UNDE_027",
            "UNDE_034","UNDE_035","UNDE_037","UNDE_038","UNDE_043","UNDE_044")
siteID <- rep("UNDE", length(plotID))
species <- rep("microbe", length(plotID))
x_UNDE_microbes <- data.frame(plotID, siteID, shd_UNDEmicrobes, rich_UNDEmicrobes,
                              species)
colnames(x_UNDE_microbes)<- c("plotID","siteID", "shd", "rich","species")

### UKFS data table with shannon and diversity
plotID <- c("UKFS_001","UKFS_002","UKFS_003","UKFS_004","UKFS_005","UKFS_009",
            "UKFS_031","UKFS_032","UKFS_043","UKFS_044")
siteID <- rep("UKFS", length(plotID))
species <- rep("microbe", length(plotID))
x_UKFS_microbes <- data.frame(plotID, siteID, shd_UKFSmicrobes, rich_UKFSmicrobes,
                              species)
colnames(x_UKFS_microbes)<- c("plotID","siteID", "shd", "rich","species")

### GRSM data table with shannon and diversity
plotID <- c("GRSM_001","GRSM_002","GRSM_003","GRSM_006","GRSM_007","GRSM_016",
            "GRSM_055","GRSM_058","GRSM_059","GRSM_060")
siteID <- rep("GRSM", length(plotID))
species <- rep("microbe", length(plotID))
x_GRSM_microbes <- data.frame(plotID, siteID, shd_GRSMmicrobes, rich_GRSMmicrobes,
                              species)
colnames(x_GRSM_microbes)<- c("plotID","siteID", "shd", "rich","species")

### WOOD data table with shannon and diversity
plotID <- c("WOOD_001","WOOD_002","WOOD_003","WOOD_004","WOOD_005","WOOD_022",
            "WOOD_024","WOOD_042","WOOD_043","WOOD_044","WOOD_045")
siteID <- rep("WOOD", length(plotID))
species <- rep("microbe", length(plotID))
x_WOOD_microbes <- data.frame(plotID, siteID, shd_WOODmicrobes, rich_WOODmicrobes,
                              species)
colnames(x_WOOD_microbes)<- c("plotID","siteID", "shd", "rich","species")

### OAES data table with shannon and diversity
plotID <- c("OAES_001","OAES_002","OAES_003","OAES_004","OAES_007","OAES_009",
            "OAES_042","OAES_043","OAES_044","OAES_045")
siteID <- rep("OAES", length(plotID))
species <- rep("microbe", length(plotID))
x_OAES_microbes <- data.frame(plotID, siteID, shd_OAESmicrobes, rich_OAESmicrobes,
                              species)
colnames(x_OAES_microbes)<- c("plotID","siteID", "shd", "rich","species")

### YELL data table with shannon and diversity
plotID <- c("YELL_001","YELL_002","YELL_003","YELL_009","YELL_012","YELL_016",
            "YELL_046","YELL_048","YELL_051","YELL_052")
siteID <- rep("YELL", length(plotID))
species <- rep("microbe", length(plotID))
x_YELL_microbes <- data.frame(plotID, siteID, shd_YELLmicrobes, rich_YELLmicrobes,
                              species)
colnames(x_YELL_microbes)<- c("plotID","siteID", "shd", "rich","species")

### NIWO data table with shannon and diversity
plotID <- c("NIWO_001","NIWO_002","NIWO_003","NIWO_004","NIWO_005","NIWO_006",
            "NIWO_008","NIWO_040","NIWO_041","NIWO_042","NIWO_043")
siteID <- rep("NIWO", length(plotID))
species <- rep("microbe", length(plotID))
x_NIWO_microbes <- data.frame(plotID, siteID, shd_NIWOmicrobes, rich_NIWOmicrobes,
                              species)
colnames(x_NIWO_microbes)<- c("plotID","siteID", "shd", "rich","species")

### SRER data table with shannon and diversity
plotID <- c("SRER_001","SRER_002","SRER_003","SRER_004","SRER_005","SRER_006",
            "SRER_043","SRER_047","SRER_052","SRER_053")
siteID <- rep("SRER", length(plotID))
species <- rep("microbe", length(plotID))
x_SRER_microbes <- data.frame(plotID, siteID, shd_SRERmicrobes, rich_SRERmicrobes,
                              species)
colnames(x_SRER_microbes)<- c("plotID","siteID", "shd", "rich","species")

### ONAQ data table with shannon and diversity
plotID <- c("ONAQ_002","ONAQ_003","ONAQ_004","ONAQ_005","ONAQ_006","ONAQ_007",
            "ONAQ_008","ONAQ_009","ONAQ_010","ONAQ_012","ONAQ_017","ONAQ_041",
            "ONAQ_042","ONAQ_043","ONAQ_044")
siteID <- rep("ONAQ", length(plotID))
species <- rep("microbe", length(plotID))
x_ONAQ_microbes <- data.frame(plotID, siteID, shd_ONAQmicrobes, rich_ONAQmicrobes,
                              species)
colnames(x_ONAQ_microbes)<- c("plotID","siteID", "shd", "rich","species")

### ABBY data table with shannon and diversity
plotID <- c("ABBY_001","ABBY_002","ABBY_003","ABBY_004","ABBY_006","ABBY_023",
            "ABBY_061","ABBY_062","ABBY_063","ABBY_070")
siteID <- rep("ABBY", length(plotID))
species <- rep("microbe", length(plotID))
x_ABBY_microbes <- data.frame(plotID, siteID, shd_ABBYmicrobes, rich_ABBYmicrobes,
                              species)
colnames(x_ABBY_microbes)<- c("plotID","siteID", "shd", "rich","species")

### SJER data table with shannon and diversity
plotID <- c("SJER_001","SJER_002","SJER_003","SJER_004","SJER_005","SJER_025",
            "SJER_045","SJER_046","SJER_047","SJER_048")
siteID <- rep("SJER", length(plotID))
species <- rep("microbe", length(plotID))
x_SJER_microbes <- data.frame(plotID, siteID, shd_SJERmicrobes, rich_SJERmicrobes,
                              species)
colnames(x_SJER_microbes)<- c("plotID","siteID", "shd", "rich","species")

### BARR data table with shannon and diversity
plotID <- c("BARR_001","BARR_002","BARR_003","BARR_004","BARR_005","BARR_006",
            "BARR_051","BARR_052","BARR_053","BARR_054")
siteID <- rep("BARR", length(plotID))
species <- rep("microbe", length(plotID))
x_BARR_microbes <- data.frame(plotID, siteID, shd_BARRmicrobes, rich_BARRmicrobes,
                              species)
colnames(x_BARR_microbes)<- c("plotID","siteID", "shd", "rich","species")

# combine data tables from each site
x_microbes <- bind_rows(x_GUAN_microbes,x_HEAL_microbes,x_HARV_microbes,
                        x_CPER_microbes,x_TALL_microbes,x_SCBI_microbes,x_OSBS_microbes,
                        x_UNDE_microbes,x_UKFS_microbes,x_GRSM_microbes,x_WOOD_microbes,
                        x_OAES_microbes,x_YELL_microbes,x_NIWO_microbes,x_SRER_microbes,
                        x_ONAQ_microbes,x_ABBY_microbes,x_SJER_microbes,x_BARR_microbes)







# bring in plant data from 19 sites
plants <- loadByProduct(dpID = "DP1.10058.001", 
                        site = c("GUAN","HEAL","HARV","CPER","TALL","SCBI",
                                 "OSBS","UNDE","UKFS","GRSM","WOOD","OAES",
                                 "YELL","NIWO","SRER","ONAQ","ABBY","SJER",
                                 "BARR"),
                        release = "RELEASE-2024", 
                        check.size = F)

# make a data table of the plant data from 10m2 plots
table_plants <- as.data.frame(plants$div_10m2Data100m2Data)


#make separate data tables for each site
plant_GUAN <- filter(table_plants, siteID=="GUAN")
plant_HEAL <- filter(table_plants, siteID=="HEAL")
plant_HARV <- filter(table_plants, siteID=="HARV")
plant_CPER <- filter(table_plants, siteID=="CPER")
plant_TALL <- filter(table_plants, siteID=="TALL")
plant_SCBI <- filter(table_plants, siteID=="SCBI")
plant_OSBS <- filter(table_plants, siteID=="OSBS")
plant_UNDE <- filter(table_plants, siteID=="UNDE")
plant_UKFS <- filter(table_plants, siteID=="UKFS")
plant_GRSM <- filter(table_plants, siteID=="GRSM")
plant_WOOD <- filter(table_plants, siteID=="WOOD")
plant_OAES <- filter(table_plants, siteID=="OAES")
plant_YELL <- filter(table_plants, siteID=="YELL")
plant_NIWO <- filter(table_plants, siteID=="NIWO")
plant_SRER <- filter(table_plants, siteID=="SRER")
plant_ONAQ <- filter(table_plants, siteID=="ONAQ")
plant_ABBY <- filter(table_plants, siteID=="ABBY")
plant_SJER <- filter(table_plants, siteID=="SJER")
plant_BARR <- filter(table_plants, siteID=="BARR")


## make data tables for each site indicating the number of species per plot, make NA 
### values equal 0 , then remove the plot names in order to calculate SHD and richness
plant_matrix_GUAN <- plant_GUAN %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_GUAN[is.na(plant_matrix_GUAN)] <- 0
plant_matrix_GUAN <- plant_matrix_GUAN[,-1]

plant_matrix_HEAL <- plant_HEAL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_HEAL[is.na(plant_matrix_HEAL)] <- 0
plant_matrix_HEAL <- plant_matrix_HEAL[,-1]

plant_matrix_HARV <- plant_HARV %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_HARV[is.na(plant_matrix_HARV)] <- 0
plant_matrix_HARV <- plant_matrix_HARV[,-1]

plant_matrix_CPER <- plant_CPER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_CPER[is.na(plant_matrix_CPER)] <- 0
plant_matrix_CPER <- plant_matrix_CPER[,-1]

plant_matrix_TALL <- plant_TALL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_TALL[is.na(plant_matrix_TALL)] <- 0
plant_matrix_TALL <- plant_matrix_TALL[,-1]

plant_matrix_SCBI <- plant_SCBI %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_SCBI[is.na(plant_matrix_SCBI)] <- 0
plant_matrix_SCBI <- plant_matrix_SCBI[,-1]

plant_matrix_OSBS <- plant_OSBS %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_OSBS[is.na(plant_matrix_OSBS)] <- 0
plant_matrix_OSBS <- plant_matrix_OSBS[,-1]

plant_matrix_UNDE <- plant_UNDE %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_UNDE[is.na(plant_matrix_UNDE)] <- 0
plant_matrix_UNDE <- plant_matrix_UNDE[,-1]

plant_matrix_UKFS <- plant_UKFS %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_UKFS[is.na(plant_matrix_UKFS)] <- 0
plant_matrix_UKFS <- plant_matrix_UKFS[,-1]

plant_matrix_GRSM <- plant_GRSM %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_GRSM[is.na(plant_matrix_GRSM)] <- 0
plant_matrix_GRSM <- plant_matrix_GRSM[,-1]

plant_matrix_WOOD <- plant_WOOD %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_WOOD[is.na(plant_matrix_WOOD)] <- 0
plant_matrix_WOOD <- plant_matrix_WOOD[,-1]

plant_matrix_OAES <- plant_OAES %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_OAES[is.na(plant_matrix_OAES)] <- 0
plant_matrix_OAES <- plant_matrix_OAES[,-1]

plant_matrix_YELL <- plant_YELL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_YELL[is.na(plant_matrix_YELL)] <- 0
plant_matrix_YELL <- plant_matrix_YELL[,-1]

plant_matrix_NIWO <- plant_NIWO %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_NIWO[is.na(plant_matrix_NIWO)] <- 0
plant_matrix_NIWO <- plant_matrix_NIWO[,-1]

plant_matrix_SRER <- plant_SRER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_SRER[is.na(plant_matrix_SRER)] <- 0
plant_matrix_SRER <- plant_matrix_SRER[,-1]

plant_matrix_ONAQ <- plant_ONAQ %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_ONAQ[is.na(plant_matrix_ONAQ)] <- 0
plant_matrix_ONAQ <- plant_matrix_ONAQ[,-1]

plant_matrix_ABBY <- plant_ABBY %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_ABBY[is.na(plant_matrix_ABBY)] <- 0
plant_matrix_ABBY <- plant_matrix_ABBY[,-1]

plant_matrix_SJER <- plant_SJER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_SJER[is.na(plant_matrix_SJER)] <- 0
plant_matrix_SJER <- plant_matrix_SJER[,-1]

plant_matrix_BARR <- plant_BARR %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
plant_matrix_BARR[is.na(plant_matrix_BARR)] <- 0
plant_matrix_BARR <- plant_matrix_BARR[,-1]

##calculating shannon diversity and richness based on each site
shd_GUANplants <- diversity(plant_matrix_GUAN, index = "shannon", MARGIN = 1, base = exp(1))
shd_HEALplants <- diversity(plant_matrix_HEAL, index = "shannon", MARGIN = 1, base = exp(1))
shd_HARVplants <- diversity(plant_matrix_HARV, index = "shannon", MARGIN = 1, base = exp(1))
shd_CPERplants <- diversity(plant_matrix_CPER, index = "shannon", MARGIN = 1, base = exp(1))
shd_TALLplants <- diversity(plant_matrix_TALL, index = "shannon", MARGIN = 1, base = exp(1))
shd_SCBIplants <- diversity(plant_matrix_SCBI, index = "shannon", MARGIN = 1, base = exp(1))
shd_OSBSplants <- diversity(plant_matrix_OSBS, index = "shannon", MARGIN = 1, base = exp(1))
shd_UNDEplants <- diversity(plant_matrix_UNDE, index = "shannon", MARGIN = 1, base = exp(1))
shd_UKFSplants <- diversity(plant_matrix_UKFS, index = "shannon", MARGIN = 1, base = exp(1))
shd_GRSMplants <- diversity(plant_matrix_GRSM, index = "shannon", MARGIN = 1, base = exp(1))
shd_WOODplants <- diversity(plant_matrix_WOOD, index = "shannon", MARGIN = 1, base = exp(1))
shd_OAESplants <- diversity(plant_matrix_OAES, index = "shannon", MARGIN = 1, base = exp(1))
shd_YELLplants <- diversity(plant_matrix_YELL, index = "shannon", MARGIN = 1, base = exp(1))
shd_NIWOplants <- diversity(plant_matrix_NIWO, index = "shannon", MARGIN = 1, base = exp(1))
shd_SRERplants <- diversity(plant_matrix_SRER, index = "shannon", MARGIN = 1, base = exp(1))
shd_ONAQplants <- diversity(plant_matrix_ONAQ, index = "shannon", MARGIN = 1, base = exp(1))
shd_ABBYplants <- diversity(plant_matrix_ABBY, index = "shannon", MARGIN = 1, base = exp(1))
shd_SJERplants <- diversity(plant_matrix_SJER, index = "shannon", MARGIN = 1, base = exp(1))
shd_BARRplants <- diversity(plant_matrix_BARR, index = "shannon", MARGIN = 1, base = exp(1))

rich_GUANplants <- specnumber(plant_matrix_GUAN, MARGIN = 1)
rich_HEALplants <- specnumber(plant_matrix_HEAL, MARGIN = 1)
rich_HARVplants <- specnumber(plant_matrix_HARV, MARGIN = 1)
rich_CPERplants <- specnumber(plant_matrix_CPER, MARGIN = 1)
rich_TALLplants <- specnumber(plant_matrix_TALL, MARGIN = 1)
rich_SCBIplants <- specnumber(plant_matrix_SCBI, MARGIN = 1)
rich_OSBSplants <- specnumber(plant_matrix_OSBS, MARGIN = 1)
rich_UNDEplants <- specnumber(plant_matrix_UNDE, MARGIN = 1)
rich_UKFSplants <- specnumber(plant_matrix_UKFS, MARGIN = 1)
rich_GRSMplants <- specnumber(plant_matrix_GRSM, MARGIN = 1)
rich_WOODplants <- specnumber(plant_matrix_WOOD, MARGIN = 1)
rich_OAESplants <- specnumber(plant_matrix_OAES, MARGIN = 1)
rich_YELLplants <- specnumber(plant_matrix_YELL, MARGIN = 1)
rich_NIWOplants <- specnumber(plant_matrix_NIWO, MARGIN = 1)
rich_SRERplants <- specnumber(plant_matrix_SRER, MARGIN = 1)
rich_ONAQplants <- specnumber(plant_matrix_ONAQ, MARGIN = 1)
rich_ABBYplants <- specnumber(plant_matrix_ABBY, MARGIN = 1)
rich_SJERplants <- specnumber(plant_matrix_SJER, MARGIN = 1)
rich_BARRplants <- specnumber(plant_matrix_BARR, MARGIN = 1)

# GUAN data table with shannon and diversity
plotID <- c("GUAN_001","GUAN_002","GUAN_003","GUAN_004","GUAN_005","GUAN_006",
            "GUAN_007","GUAN_008","GUAN_009","GUAN_010","GUAN_011","GUAN_012",
            "GUAN_013","GUAN_014","GUAN_015","GUAN_016","GUAN_017","GUAN_018",
            "GUAN_019","GUAN_020","GUAN_021","GUAN_022","GUAN_023","GUAN_024",
            "GUAN_042","GUAN_043","GUAN_048")
siteID <- rep("GUAN", length(plotID))
species <- rep("plant", length(plotID))
x_GUAN_plants <- data.frame(plotID, siteID, shd_GUANplants, rich_GUANplants,species)
colnames(x_GUAN_plants)<- c("plotID","siteID", "shd", "rich","species")

# HEAL data table with shannon and diversity
plotID <- c("HEAL_001","HEAL_002","HEAL_003","HEAL_004","HEAL_005","HEAL_006",
            "HEAL_007","HEAL_008","HEAL_009","HEAL_010","HEAL_011","HEAL_012",
            "HEAL_013","HEAL_014","HEAL_015","HEAL_016","HEAL_017","HEAL_018",
            "HEAL_019","HEAL_020","HEAL_021","HEAL_022","HEAL_023","HEAL_024",
            "HEAL_025","HEAL_026","HEAL_027","HEAL_028","HEAL_029","HEAL_045",
            "HEAL_046","HEAL_047","HEAL_077")
siteID <- rep("HEAL", length(plotID))
species <- rep("plant", length(plotID))
x_HEAL_plants <- data.frame(plotID, siteID, shd_HEALplants, rich_HEALplants,species)
colnames(x_HEAL_plants)<- c("plotID","siteID", "shd", "rich","species")

# HARV data table with shannon and diversity
plotID <- c("HARV_001","HARV_002","HARV_004","HARV_005","HARV_006","HARV_008",
            "HARV_010","HARV_011","HARV_012","HARV_013","HARV_014","HARV_015",
            "HARV_016","HARV_017","HARV_018","HARV_020","HARV_021","HARV_022",
            "HARV_023","HARV_024","HARV_025","HARV_026","HARV_027","HARV_028",
            "HARV_029","HARV_030","HARV_031","HARV_033","HARV_034","HARV_035",
            "HARV_058","HARV_059","HARV_063")
siteID <- rep("HARV", length(plotID))
species <- rep("plant", length(plotID))
x_HARV_plants <- data.frame(plotID, siteID, shd_HARVplants, rich_HARVplants,species)
colnames(x_HARV_plants)<- c("plotID","siteID", "shd", "rich","species")

# CPER data table with shannon and diversity
plotID <- c("CPER_001","CPER_002","CPER_003","CPER_004","CPER_005","CPER_006",
            "CPER_007","CPER_008","CPER_009","CPER_010","CPER_011","CPER_012",
            "CPER_013","CPER_014","CPER_015","CPER_016","CPER_017","CPER_018",
            "CPER_019","CPER_020","CPER_021","CPER_022","CPER_023","CPER_024",
            "CPER_025","CPER_026","CPER_027","CPER_028","CPER_029","CPER_030",
            "CPER_031","CPER_033","CPER_034","CPER_035","CPER_036","CPER_037",
            "CPER_038","CPER_039","CPER_040","CPER_041","CPER_045","CPER_046",
            "CPER_047","CPER_048","CPER_049","CPER_050","CPER_051","CPER_052",
            "CPER_053","CPER_054","CPER_055","CPER_056","CPER_057","CPER_058",
            "CPER_059","CPER_060","CPER_061","CPER_062","CPER_063","CPER_064",
            "CPER_065","CPER_066","CPER_067","CPER_068","CPER_069","CPER_070",
            "CPER_071","CPER_072","CPER_073","CPER_074")
siteID <- rep("CPER", length(plotID))
species <- rep("plant", length(plotID))
x_CPER_plants <- data.frame(plotID, siteID, shd_CPERplants, rich_CPERplants,species)
colnames(x_CPER_plants)<- c("plotID","siteID", "shd", "rich","species")

# TALL data with shannon and diversity
plotID <- c("TALL_001","TALL_002","TALL_003","TALL_004","TALL_005","TALL_006",
            "TALL_007","TALL_008","TALL_009","TALL_010","TALL_011","TALL_012",
            "TALL_013","TALL_015","TALL_016","TALL_017","TALL_018","TALL_019",
            "TALL_020","TALL_021","TALL_022","TALL_023","TALL_024","TALL_025",
            "TALL_026","TALL_027","TALL_029","TALL_030","TALL_031","TALL_032",
            "TALL_044","TALL_051","TALL_054","TALL_064")
siteID <- rep("TALL", length(plotID))
species <- rep("plant", length(plotID))
x_TALL_plants <- data.frame(plotID, siteID, shd_TALLplants, rich_TALLplants,species)
colnames(x_TALL_plants)<- c("plotID","siteID", "shd", "rich","species")

# SCBI data with shannon and diversity
plotID <- c("SCBI_002","SCBI_003","SCBI_004","SCBI_005","SCBI_006","SCBI_007",
            "SCBI_008","SCBI_010","SCBI_011","SCBI_012","SCBI_013","SCBI_014",
            "SCBI_015","SCBI_016","SCBI_017","SCBI_018","SCBI_019","SCBI_021",
            "SCBI_022","SCBI_023","SCBI_033","SCBI_034","SCBI_035","SCBI_037",
            "SCBI_038","SCBI_039","SCBI_040","SCBI_041","SCBI_042","SCBI_043",
            "SCBI_044","SCBI_045","SCBI_047","SCBI_067")
siteID <- rep("SCBI", length(plotID))
species <- rep("plant", length(plotID))
x_SCBI_plants <- data.frame(plotID, siteID, shd_SCBIplants, rich_SCBIplants,species)
colnames(x_SCBI_plants)<- c("plotID","siteID", "shd", "rich","species")

# OSBS data with shannon and diversity
plotID <- c("OSBS_001","OSBS_002","OSBS_003","OSBS_004","OSBS_005","OSBS_006",
            "OSBS_007","OSBS_008","OSBS_009","OSBS_010","OSBS_011","OSBS_012",
            "OSBS_013","OSBS_014","OSBS_015","OSBS_016","OSBS_017","OSBS_018",
            "OSBS_019","OSBS_020","OSBS_021","OSBS_022","OSBS_023","OSBS_024",
            "OSBS_025","OSBS_026","OSBS_027","OSBS_029","OSBS_048","OSBS_049",
            "OSBS_050","OSBS_051","OSBS_063","OSBS_065","OSBS_066","OSBS_067",
            "OSBS_068","OSBS_069","OSBS_070","OSBS_071","OSBS_072","OSBS_073",
            "OSBS_074","OSBS_075","OSBS_076")
siteID <- rep("OSBS", length(plotID))
species <- rep("plant", length(plotID))
x_OSBS_plants <- data.frame(plotID, siteID, shd_OSBSplants, rich_OSBSplants,species)
colnames(x_OSBS_plants)<- c("plotID","siteID", "shd", "rich","species")

# UNDE data with shannon and diversity
plotID <- c("UNDE_001","UNDE_002","UNDE_003","UNDE_006","UNDE_007","UNDE_008",
            "UNDE_010","UNDE_011","UNDE_012","UNDE_013","UNDE_014","UNDE_015",
            "UNDE_016","UNDE_017","UNDE_018","UNDE_019","UNDE_020","UNDE_021",
            "UNDE_022","UNDE_023","UNDE_024","UNDE_025","UNDE_027","UNDE_028",
            "UNDE_029","UNDE_030","UNDE_032","UNDE_034","UNDE_035","UNDE_036",
            "UNDE_037","UNDE_038","UNDE_043","UNDE_044","UNDE_077")
siteID <- rep("UNDE", length(plotID))
species <- rep("plant", length(plotID))
x_UNDE_plants <- data.frame(plotID, siteID, shd_UNDEplants, rich_UNDEplants,species)
colnames(x_UNDE_plants)<- c("plotID","siteID", "shd", "rich","species")

# UKFS data with shannon and diversity
plotID <- c("UKFS_001","UKFS_002","UKFS_003","UKFS_004","UKFS_005","UKFS_006",
            "UKFS_007","UKFS_008","UKFS_009","UKFS_010","UKFS_011","UKFS_012",
            "UKFS_013","UKFS_014","UKFS_015","UKFS_016","UKFS_017","UKFS_018",
            "UKFS_019","UKFS_020","UKFS_021","UKFS_022","UKFS_023","UKFS_024",
            "UKFS_025","UKFS_026","UKFS_027","UKFS_028","UKFS_029","UKFS_030",
            "UKFS_031","UKFS_032","UKFS_043")
siteID <- rep("UKFS", length(plotID))
species <- rep("plant", length(plotID))
x_UKFS_plants <- data.frame(plotID, siteID, shd_UKFSplants, rich_UKFSplants,species)
colnames(x_UKFS_plants)<- c("plotID","siteID", "shd", "rich","species")

# GRSM data with shannon and diversity
plotID <- c("GRSM_001","GRSM_002","GRSM_003","GRSM_004","GRSM_005","GRSM_006",
            "GRSM_007","GRSM_008","GRSM_009","GRSM_010","GRSM_011","GRSM_012",
            "GRSM_013","GRSM_014","GRSM_015","GRSM_016","GRSM_017","GRSM_018",
            "GRSM_019","GRSM_020","GRSM_021","GRSM_022","GRSM_023","GRSM_024",
            "GRSM_025","GRSM_026","GRSM_027","GRSM_028","GRSM_029","GRSM_030",
            "GRSM_047","GRSM_048","GRSM_050","GRSM_055","GRSM_058","GRSM_059",
            "GRSM_072")
siteID <- rep("GRSM", length(plotID))
species <- rep("plant", length(plotID))
x_GRSM_plants <- data.frame(plotID, siteID, shd_GRSMplants, rich_GRSMplants,species)
colnames(x_GRSM_plants)<- c("plotID","siteID", "shd", "rich","species")

# WOOD data with shannon and diversity
plotID <- c("WOOD_001","WOOD_002","WOOD_003","WOOD_004","WOOD_005","WOOD_006",
            "WOOD_007","WOOD_008","WOOD_009","WOOD_010","WOOD_011","WOOD_012",
            "WOOD_013","WOOD_014","WOOD_015","WOOD_016","WOOD_017","WOOD_018",
            "WOOD_019","WOOD_020","WOOD_021","WOOD_022","WOOD_023","WOOD_024",
            "WOOD_025","WOOD_026","WOOD_027","WOOD_028","WOOD_029","WOOD_030",
            "WOOD_042","WOOD_043","WOOD_044")
siteID <- rep("WOOD", length(plotID))
species <- rep("plant", length(plotID))
x_WOOD_plants <- data.frame(plotID, siteID, shd_WOODplants, rich_WOODplants,species)
colnames(x_WOOD_plants)<- c("plotID","siteID", "shd", "rich","species")

# OAES data with shannon and diversity
plotID <- c("OAES_001","OAES_002","OAES_003","OAES_004","OAES_005","OAES_006",
            "OAES_007","OAES_008","OAES_009","OAES_010","OAES_011","OAES_012",
            "OAES_013","OAES_014","OAES_015","OAES_016","OAES_017","OAES_018",
            "OAES_019","OAES_020","OAES_021","OAES_022","OAES_023","OAES_024",
            "OAES_025","OAES_026","OAES_027","OAES_028","OAES_029","OAES_030",
            "OAES_042","OAES_043","OAES_044")
siteID <- rep("OAES", length(plotID))
species <- rep("plant", length(plotID))
x_OAES_plants <- data.frame(plotID, siteID, shd_OAESplants, rich_OAESplants,species)
colnames(x_OAES_plants)<- c("plotID","siteID", "shd", "rich","species")

# YELL data with shannon and diversity
plotID <- c("YELL_001","YELL_002","YELL_003","YELL_004","YELL_005","YELL_006",
            "YELL_007","YELL_008","YELL_009","YELL_010","YELL_011","YELL_012",
            "YELL_013","YELL_014","YELL_015","YELL_016","YELL_017","YELL_018",
            "YELL_019","YELL_020","YELL_021","YELL_022","YELL_023","YELL_024",
            "YELL_025","YELL_026","YELL_027","YELL_028","YELL_029","YELL_030",
            "YELL_046","YELL_048","YELL_051","YELL_071")
siteID <- rep("YELL", length(plotID))
species <- rep("plant", length(plotID))
x_YELL_plants <- data.frame(plotID, siteID, shd_YELLplants, rich_YELLplants,species)
colnames(x_YELL_plants)<- c("plotID","siteID", "shd", "rich","species")

# NIWO data with shannon and diversity
plotID <- c("NIWO_001","NIWO_002","NIWO_003","NIWO_004","NIWO_005","NIWO_006",
            "NIWO_007","NIWO_008","NIWO_009","NIWO_010","NIWO_011","NIWO_012",
            "NIWO_013","NIWO_014","NIWO_015","NIWO_016","NIWO_017","NIWO_018",
            "NIWO_019","NIWO_020","NIWO_021","NIWO_022","NIWO_023","NIWO_024",
            "NIWO_025","NIWO_026","NIWO_027","NIWO_028","NIWO_029","NIWO_030",
            "NIWO_040","NIWO_041","NIWO_042")
siteID <- rep("NIWO", length(plotID))
species <- rep("plant", length(plotID))
x_NIWO_plants <- data.frame(plotID, siteID, shd_NIWOplants, rich_NIWOplants,species)
colnames(x_NIWO_plants)<- c("plotID","siteID", "shd", "rich","species")

# SRER data with shannon and diversity
plotID <- c("SRER_001","SRER_002","SRER_003","SRER_004","SRER_005","SRER_006",
            "SRER_007","SRER_008","SRER_009","SRER_010","SRER_011","SRER_012",
            "SRER_013","SRER_014","SRER_015","SRER_016","SRER_017","SRER_018",
            "SRER_019","SRER_020","SRER_021","SRER_022","SRER_023","SRER_024",
            "SRER_025","SRER_026","SRER_027","SRER_028","SRER_029","SRER_030",
            "SRER_043","SRER_047","SRER_052")
siteID <- rep("SRER", length(plotID))
species <- rep("plant", length(plotID))
x_SRER_plants <- data.frame(plotID, siteID, shd_SRERplants, rich_SRERplants,species)
colnames(x_SRER_plants)<- c("plotID","siteID", "shd", "rich","species")

# ONAQ data with shannon and diversity
plotID <- c("ONAQ_002","ONAQ_003","ONAQ_004","ONAQ_005","ONAQ_006","ONAQ_007",
            "ONAQ_008","ONAQ_009","ONAQ_010","ONAQ_011","ONAQ_012","ONAQ_014",
            "ONAQ_015","ONAQ_016","ONAQ_017","ONAQ_018","ONAQ_019","ONAQ_020",
            "ONAQ_021","ONAQ_022","ONAQ_023","ONAQ_024","ONAQ_025","ONAQ_026",
            "ONAQ_027","ONAQ_028","ONAQ_029","ONAQ_030","ONAQ_031","ONAQ_032",
            "ONAQ_041","ONAQ_042","ONAQ_043","ONAQ_073")
siteID <- rep("ONAQ", length(plotID))
species <- rep("plant", length(plotID))
x_ONAQ_plants <- data.frame(plotID, siteID, shd_ONAQplants, rich_ONAQplants,species)
colnames(x_ONAQ_plants)<- c("plotID","siteID", "shd", "rich","species")

# ABBY data with shannon and diversity
plotID <- c("ABBY_001","ABBY_002","ABBY_003","ABBY_004","ABBY_005","ABBY_006",
            "ABBY_007","ABBY_008","ABBY_009","ABBY_010","ABBY_011","ABBY_012",
            "ABBY_013","ABBY_014","ABBY_015","ABBY_016","ABBY_017","ABBY_018",
            "ABBY_019","ABBY_020","ABBY_021","ABBY_022","ABBY_023","ABBY_024",
            "ABBY_025","ABBY_026","ABBY_027","ABBY_028","ABBY_029","ABBY_030",
            "ABBY_061","ABBY_063","ABBY_070","ABBY_077")
siteID <- rep("ABBY", length(plotID))
species <- rep("plant", length(plotID))
x_ABBY_plants <- data.frame(plotID, siteID, shd_ABBYplants, rich_ABBYplants,species)
colnames(x_ABBY_plants)<- c("plotID","siteID", "shd", "rich","species")

# SJER data with shannon and diversity
plotID <- c("SJER_001","SJER_002","SJER_003","SJER_004","SJER_005","SJER_006",
            "SJER_007","SJER_008","SJER_009","SJER_010","SJER_011","SJER_012",
            "SJER_013","SJER_014","SJER_015","SJER_016","SJER_017","SJER_018",
            "SJER_019","SJER_020","SJER_021","SJER_022","SJER_023","SJER_024",
            "SJER_025","SJER_026","SJER_027","SJER_028","SJER_029","SJER_030",
            "SJER_045","SJER_046","SJER_047")
siteID <- rep("SJER", length(plotID))
species <- rep("plant", length(plotID))
x_SJER_plants <- data.frame(plotID, siteID, shd_SJERplants, rich_SJERplants,species)
colnames(x_SJER_plants)<- c("plotID","siteID", "shd", "rich","species")

# BARR data with shannon and diversity
plotID <- c("BARR_001","BARR_002","BARR_003","BARR_004","BARR_005","BARR_006",
            "BARR_007","BARR_008","BARR_009","BARR_011","BARR_012","BARR_013",
            "BARR_014","BARR_015","BARR_016","BARR_017","BARR_018","BARR_019",
            "BARR_020","BARR_021","BARR_022","BARR_023","BARR_024","BARR_025",
            "BARR_026","BARR_027","BARR_028","BARR_029","BARR_030","BARR_032",
            "BARR_051","BARR_052","BARR_053")
siteID <- rep("BARR", length(plotID))
species <- rep("plant", length(plotID))
x_BARR_plants <- data.frame(plotID, siteID, shd_BARRplants, rich_BARRplants,species)
colnames(x_BARR_plants)<- c("plotID","siteID", "shd", "rich","species")



# combine plant data tables from each site
x_plants <- bind_rows(x_GUAN_plants, x_HEAL_plants, x_HARV_plants, x_CPER_plants,
                      x_TALL_plants, x_SCBI_plants, x_OSBS_plants, x_UNDE_plants,
                      x_UKFS_plants, x_GRSM_plants, x_WOOD_plants, x_OAES_plants,
                      x_YELL_plants, x_NIWO_plants, x_SRER_plants, x_ONAQ_plants,
                      x_ABBY_plants, x_SJER_plants, x_BARR_plants)


# combine plant and microbe data
x <- bind_rows(x_microbes,x_plants)


# make a new data table with both plant and microbe data in each row by plotID
merged_x <- merge(x_microbes, x_plants, by = "plotID")

# clean up names of columns
setnames(merged_x, old = "rich.x", new = "rich.mic")
setnames(merged_x, old = "rich.y", new = "rich.pla")
setnames(merged_x, old = "shd.x", new = "shd.mic")
setnames(merged_x, old = "shd.y", new = "shd.pla")
setnames(merged_x, old = "siteID.x", new = "siteID")

# make a graph that compares plant and microbial species richness by site ID
ggplot(merged_x, aes(x = rich.mic, y = rich.pla, color = siteID)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    title = "Relationship Between Microbial and Plant Richness",
    x = "Microbe Richness",
    y = "Plant Richness"
  ) +
  theme_minimal()

# run statistics on richness graph
model <- lm(rich.pla ~ rich.mic, data = merged_x)
summary(model)

# make a graph that compares plant and microbial SHD by site ID
ggplot(merged_x, aes(x = shd.mic, y = shd.pla, color = siteID)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(
    title = "Relationship Between Microbial and Plant Shannon",
    x = "Microbe Shannon",
    y = "Plant Shannon"
  ) +
  theme_minimal()

# run statistics on shannon graph
model <- lm(shd.pla ~ shd.mic, data = merged_x)
summary(model)











# bring in soil temperature and moisture data
soil <- loadByProduct(dpID = "DP1.10086.001", 
                      site = c("GUAN","HEAL","HARV","CPER","TALL","SCBI",
                               "OSBS","UNDE","UKFS","GRSM","WOOD","OAES",
                               "YELL","NIWO","SRER","ONAQ","ABBY","SJER",
                               "BARR"),
                      release = "RELEASE-2024", 
                      check.size = F)

temp_table <- as.data.frame(soil$sls_soilCoreCollection)

#make a separate temperature data table with only site, plot, and soil temp data
temp_data <- temp_table %>% select(siteID, plotID, soilTemp)
setDT(temp_data)
temp_data <- temp_data[!is.na(soilTemp)]
temp_data <- temp_data[order(plotID)]

#find the average soil temperature for each site, round it to nearest whole number
avg_temp <- temp_data[, .(mean_temp = round(mean(soilTemp, na.rm = TRUE), 1)), by = siteID]

# merge soil temperature data with plant and microbe data
temp <- merge(avg_temp, x, by = "siteID")

# make the soil temperature column into a factor rather than a numeric value
## to graph shannon and richness by temperature
temp$mean_temp <- factor(temp$mean_temp)


# graph shannon and diversity as boxplots by soil temperature rather than site/plot
ggplot(temp, aes(x = mean_temp, y = shd, fill = species)) +
  geom_boxplot(position = position_dodge()) +
  labs(
    title = "Shannon Diversity Index by Soil Temperature",
    x = "Soil Temperature (˚C)",
    y = "Shannon Diversity"
  ) +
  theme_classic()

model <- lm(shd ~ mean_temp, data = temp)
summary(model)

ggplot(temp, aes(x = mean_temp, y = rich, fill = species)) +
  geom_boxplot(position = position_dodge()) +
  labs(
    title = "Species Richness by Soil Temperature",
    x = "Soil Temperature (˚C)",
    y = "Species Richness"
  ) +
  theme_classic()

model <- lm(rich ~ mean_temp, data = temp)
summary(model)



x_temp <- merge(merged_x, avg_temp, by = "siteID")







# make a graph that compares plant and microbial species richness by temp
ggplot(x_temp, aes(x = rich.mic, y = rich.pla, color = mean_temp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "gold", high = "red2") +
  labs(
    title = "Relationship Between Microbial and Plant Richness",
    x = "Microbe Richness",
    y = "Plant Richness",
    color = "Soil Temp (˚C)"
  ) +
  theme_minimal()


# split this table into three categories by temperature: low, medium and high
x_temp <- x_temp %>%
  mutate(temp_group = case_when(
    mean_temp <= quantile(mean_temp, 0.33) ~ "Low Temp",
    mean_temp <= quantile(mean_temp, 0.66) ~ "Medium Temp",
    TRUE ~ "High Temp"
  ))

# graph the table to show correlations by each temperature group
ggplot(x_temp, aes(x = rich.mic, y = rich.pla, color = temp_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Low Temp" = "purple", "Medium Temp" = "orange", "High Temp" = "red")) +
  labs(
    title = "Relationship Between Microbe and Plant Richness",
    x = "Microbe Richness",
    y = "Plant Richness",
    color = "Temperature"
  ) +
  theme_minimal()

# run statistics on richness graph
model <- lm(rich.pla ~ rich.mic, data = x_temp)
summary(model)
model <- lm(rich.mic ~ rich.pla + mean_temp + rich.pla*mean_temp, data = x_temp)
summary(model)

# statistics on richness graph by soil temperature group
model <- lm(rich.mic ~ rich.pla + temp_group + rich.pla*temp_group, data = x_temp)
summary(model)
model <- lm(rich.mic ~ rich.pla, data = subset(x_temp, temp_group == "Low Temp"))
summary(model)
model <- lm(rich.mic ~ rich.pla, data = subset(x_temp, temp_group == "Medium Temp"))
summary(model)
model <- lm(rich.mic ~ rich.pla, data = subset(x_temp, temp_group == "High Temp"))
summary(model)

# graph the table to show correlations by each temperature group
ggplot(x_temp, aes(x = shd.mic, y = shd.pla, color = temp_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Low Temp" = "purple", "Medium Temp" = "orange", "High Temp" = "red")) +
  labs(
    title = "Relationship Between Microbe and Plant Shannon",
    x = "Microbe Shannon",
    y = "Plant Shannon",
    color = "Temperature"
  ) +
  theme_minimal()

# statistics on richness graph by soil temperature group
model <- lm(shd.mic ~ shd.pla + temp_group + shd.pla*temp_group, data = x_temp)
summary(model)
model <- lm(shd.mic ~ shd.pla, data = subset(x_temp, temp_group == "Low Temp"))
summary(model)
model <- lm(shd.mic ~ shd.pla, data = subset(x_temp, temp_group == "Medium Temp"))
summary(model)
model <- lm(shd.mic ~ shd.pla, data = subset(x_temp, temp_group == "High Temp"))
summary(model)



# split this table into three categories by moisture: low, medium and high
x_soil <- x_soil %>%
  mutate(moist_group = case_when(
    mean_moist <= quantile(mean_moist, 0.33) ~ "Dry",
    mean_moist <= quantile(mean_moist, 0.66) ~ "Medium",
    TRUE ~ "Wet"
  ))

# graph the table to show correlations by each temperature group
ggplot(x_soil, aes(x = rich.mic, y = rich.pla, color = moist_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Dry" = "skyblue1", "Medium" = "royalblue", "Wet" = "navy")) +
  labs(
    title = "Relationship Between Microbe and Plant Richness",
    x = "Microbe Richness",
    y = "Plant Richness",
    color = "Moisture"
  ) +
  theme_minimal()


# statistics on richness graph by soil temperature group
model <- lm(rich.mic ~ rich.pla + moist_group + rich.pla*moist_group, data = x_soil)
summary(model)
model <- lm(rich.mic ~ rich.pla, data = subset(x_soil, moist_group == "Dry"))
summary(model)
model <- lm(rich.mic ~ rich.pla, data = subset(x_soil, moist_group == "Medium"))
summary(model)
model <- lm(rich.mic ~ rich.pla, data = subset(x_soil, moist_group == "Wet"))
summary(model)

# graph the table to show correlations by each temperature group
ggplot(x_soil, aes(x = shd.mic, y = shd.pla, color = moist_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Dry" = "skyblue1", "Medium" = "royalblue", "Wet" = "navy")) +
  labs(
    title = "Relationship Between Microbe and Plant Shannon",
    x = "Microbe Shannon",
    y = "Plant Shannon",
    color = "Moisture"
  ) +
  theme_minimal()


# statistics on richness graph by soil temperature group
model <- lm(shd.mic ~ shd.pla + moist_group + shd.pla*moist_group, data = x_soil)
summary(model)
model <- lm(shd.mic ~ shd.pla, data = subset(x_soil, moist_group == "Dry"))
summary(model)
model <- lm(shd.mic ~ shd.pla, data = subset(x_soil, moist_group == "Medium"))
summary(model)
model <- lm(shd.mic ~ shd.pla, data = subset(x_soil, moist_group == "Wet"))
summary(model)












# make a graph that compares plant and microbial SHD by site ID
ggplot(x_temp, aes(x = shd.mic, y = shd.pla, color = mean_temp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "gold", high = "red2") +
  labs(
    title = "Relationship Between Microbial and Plant Shannon",
    x = "Microbe Shannon",
    y = "Plant Shannon",
    color = "Soil Temp"
  ) +
  theme_minimal()

# run statistics on richness graph
model <- lm(shd.mic ~ shd.pla, data = x_temp)
summary(model)
model <- lm(shd.mic ~ shd.pla + mean_temp + shd.pla*mean_temp, data = x_temp)
model <- lm(shd.mic ~ shd.pla + mean_temp, data = x_temp)
summary(model)



# moisture data
moisture_table <- as.data.frame(soil$sls_soilMoisture)

#make a separate temperature data table with only site, plot, and soil moisture data
moist_data <- moisture_table %>% select(siteID, plotID, soilMoisture)
setDT(moist_data)
moist_data <- moist_data[!is.na(soilMoisture)]
moist_data <- moist_data[order(plotID)]

#find the average soil moisture for each site, round it to three decimal places
avg_moist <- moist_data[, .(mean_moist = round(mean(soilMoisture, na.rm = TRUE), 3)), by = siteID]

# merge soil temperature data with plant and microbe data
moist <- merge(avg_moist, x, by = "siteID")

# make the soil temperature column into a factor rather than a numeric value
## to graph shannon and richness by temperature
moist$mean_moist <- factor(moist$mean_moist)




# graph shannon and diversity as boxplots by soil temperature rather than site/plot
ggplot(moist, aes(x = mean_moist, y = shd, fill = species)) +
  geom_boxplot(position = position_dodge()) +
  labs(
    title = "Shannon Diversity Index by Soil Moisture",
    x = "Soil Moisture (g water/g dry soil)",
    y = "Shannon Diversity"
  ) +
  theme_classic()

model <- lm(shd ~ mean_moist, data = moist)
summary(model)

ggplot(moist, aes(x = mean_moist, y = rich, fill = species)) +
  geom_boxplot(position = position_dodge()) +
  labs(
    title = "Species Richness by Soil Moisture",
    x = "Soil Moisture (g water/g dry soil)",
    y = "Species Richness"
  ) +
  theme_classic()

model <- lm(rich ~ mean_moist, data = moist)
summary(model)



x_soil <- merge(x_temp, avg_moist, by = "siteID")


# make a graph that compares plant and microbial species richness by soil moisture
ggplot(x_soil, aes(x = rich.mic, y = rich.pla, color = mean_moist)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Relationship Between Microbial and Plant Richness",
    x = "Microbe Richness",
    y = "Plant Richness",
    color = "Soil Moisture
(g water/g dry soil)"
  ) +
  theme_minimal()

model <- lm(rich.mic ~ rich.pla + mean_moist + rich.pla*mean_moist, data = x_soil)
summary(model)

# run statistics on richness graph
model <- lm(rich.pla ~ rich.mic, data = x_soil)
summary(model)
model <- lm(rich.pla + mean_moist ~ rich.mic, data = x_soil)
summary(model)
model <- lm(rich.mic ~ rich.pla + mean_moist + mean_temp, data = x_soil)
summary(model)

# make a graph that compares plant and microbial SHD by soil moisture
ggplot(x_soil, aes(x = shd.mic, y = shd.pla, color = mean_moist)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Relationship Between Microbial and Plant Shannon",
    x = "Microbe Shannon",
    y = "Plant Shannon",
    color = "Soil Moisture
(g water/g dry soil)"
  ) +
  theme_minimal()


model <- lm(shd.mic ~ shd.pla + mean_moist + shd.pla*mean_moist, data = x_soil)
summary(model)

# run statistics on shannon graph
model <- lm(shd.pla ~ shd.mic, data = x_soil)
summary(model)
model <- lm(shd.pla + mean_moist ~ shd.mic, data = x_soil)
summary(model)
model <- lm(shd.pla + mean_moist + mean_temp ~ shd.mic, data = x_soil)
summary(model)


model <- lm(rich.pla + mean_moist ~ rich.mic, data = x_soil)
summary(model)
model <- lm(rich.pla + mean_moist + mean_temp ~ rich.mic, data = x_soil)
summary(model)



# just fungi data (merged_data2)

# make separate data tables for each site
fungi_GUAN <- filter(merged_data2, siteID=="GUAN")
fungi_HEAL <- filter(merged_data2, siteID=="HEAL")
fungi_HARV <- filter(merged_data2, siteID=="HARV")
fungi_CPER <- filter(merged_data2, siteID=="CPER")
fungi_TALL <- filter(merged_data2, siteID=="TALL")
fungi_SCBI <- filter(merged_data2, siteID=="SCBI")
fungi_OSBS <- filter(merged_data2, siteID=="OSBS")
fungi_UNDE <- filter(merged_data2, siteID=="UNDE")
fungi_UKFS <- filter(merged_data2, siteID=="UKFS")
fungi_GRSM <- filter(merged_data2, siteID=="GRSM")
fungi_WOOD <- filter(merged_data2, siteID=="WOOD")
fungi_OAES <- filter(merged_data2, siteID=="OAES")
fungi_YELL <- filter(merged_data2, siteID=="YELL")
fungi_NIWO <- filter(merged_data2, siteID=="NIWO")
fungi_SRER <- filter(merged_data2, siteID=="SRER")
fungi_ONAQ <- filter(merged_data2, siteID=="ONAQ")
fungi_ABBY <- filter(merged_data2, siteID=="ABBY")
fungi_SJER <- filter(merged_data2, siteID=="SJER")
fungi_BARR <- filter(merged_data2, siteID=="BARR")


## make data tables for each site indicating the number of species per plot, make NA 
### values equal 0 , then remove the plot names in order to calculate SHD and richness
fungi_matrix_GUAN <- fungi_GUAN %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_GUAN[is.na(fungi_matrix_GUAN)] <- 0
fungi_matrix_GUAN <- fungi_matrix_GUAN[,-1]

fungi_matrix_HEAL <- fungi_HEAL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_HEAL[is.na(fungi_matrix_HEAL)] <- 0
fungi_matrix_HEAL <- fungi_matrix_HEAL[,-1]

fungi_matrix_HARV <- fungi_HARV %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_HARV[is.na(fungi_matrix_HARV)] <- 0
fungi_matrix_HARV <- fungi_matrix_HARV[,-1]

fungi_matrix_CPER <- fungi_CPER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_CPER[is.na(fungi_matrix_CPER)] <- 0
fungi_matrix_CPER <- fungi_matrix_CPER[,-1]

fungi_matrix_TALL <- fungi_TALL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_TALL[is.na(fungi_matrix_TALL)] <- 0
fungi_matrix_TALL <- fungi_matrix_TALL[,-1]

fungi_matrix_SCBI <- fungi_SCBI %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_SCBI[is.na(fungi_matrix_SCBI)] <- 0
fungi_matrix_SCBI <- fungi_matrix_SCBI[,-1]

fungi_matrix_OSBS <- fungi_OSBS %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_OSBS[is.na(fungi_matrix_OSBS)] <- 0
fungi_matrix_OSBS <- fungi_matrix_OSBS[,-1]

fungi_matrix_UNDE <- fungi_UNDE %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_UNDE[is.na(fungi_matrix_UNDE)] <- 0
fungi_matrix_UNDE <- fungi_matrix_UNDE[,-1]

fungi_matrix_UKFS <- fungi_UKFS %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_UKFS[is.na(fungi_matrix_UKFS)] <- 0
fungi_matrix_UKFS <- fungi_matrix_UKFS[,-1]

fungi_matrix_GRSM <- fungi_GRSM %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_GRSM[is.na(fungi_matrix_GRSM)] <- 0
fungi_matrix_GRSM <- fungi_matrix_GRSM[,-1]

fungi_matrix_WOOD <- fungi_WOOD %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_WOOD[is.na(fungi_matrix_WOOD)] <- 0
fungi_matrix_WOOD <- fungi_matrix_WOOD[,-1]

fungi_matrix_OAES <- fungi_OAES %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_OAES[is.na(fungi_matrix_OAES)] <- 0
fungi_matrix_OAES <- fungi_matrix_OAES[,-1]

fungi_matrix_YELL <- fungi_YELL %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_YELL[is.na(fungi_matrix_YELL)] <- 0
fungi_matrix_YELL <- fungi_matrix_YELL[,-1]

fungi_matrix_NIWO <- fungi_NIWO %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_NIWO[is.na(fungi_matrix_NIWO)] <- 0
fungi_matrix_NIWO <- fungi_matrix_NIWO[,-1]

fungi_matrix_SRER <- fungi_SRER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_SRER[is.na(fungi_matrix_SRER)] <- 0
fungi_matrix_SRER <- fungi_matrix_SRER[,-1]

fungi_matrix_ONAQ <- fungi_ONAQ %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_ONAQ[is.na(fungi_matrix_ONAQ)] <- 0
fungi_matrix_ONAQ <- fungi_matrix_ONAQ[,-1]

fungi_matrix_ABBY <- fungi_ABBY %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_ABBY[is.na(fungi_matrix_ABBY)] <- 0
fungi_matrix_ABBY <- fungi_matrix_ABBY[,-1]

fungi_matrix_SJER <- fungi_SJER %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_SJER[is.na(fungi_matrix_SJER)] <- 0
fungi_matrix_SJER <- fungi_matrix_SJER[,-1]

fungi_matrix_BARR <- fungi_BARR %>% group_by(plotID, scientificName) %>% 
  summarise(count=n()) %>%
  spread(scientificName,count)
fungi_matrix_BARR[is.na(fungi_matrix_BARR)] <- 0
fungi_matrix_BARR <- fungi_matrix_BARR[,-1]

##calculating shannon diversity and richness for each site
shd_GUANfungi <- diversity(fungi_matrix_GUAN, index = "shannon", MARGIN = 1, base = exp(1))
shd_HEALfungi <- diversity(fungi_matrix_HEAL, index = "shannon", MARGIN = 1, base = exp(1))
shd_HARVfungi <- diversity(fungi_matrix_HARV, index = "shannon", MARGIN = 1, base = exp(1))
shd_CPERfungi <- diversity(fungi_matrix_CPER, index = "shannon", MARGIN = 1, base = exp(1))
shd_TALLfungi <- diversity(fungi_matrix_TALL, index = "shannon", MARGIN = 1, base = exp(1))
shd_SCBIfungi <- diversity(fungi_matrix_SCBI, index = "shannon", MARGIN = 1, base = exp(1))
shd_OSBSfungi <- diversity(fungi_matrix_OSBS, index = "shannon", MARGIN = 1, base = exp(1))
shd_UNDEfungi <- diversity(fungi_matrix_UNDE, index = "shannon", MARGIN = 1, base = exp(1))
shd_UKFSfungi <- diversity(fungi_matrix_UKFS, index = "shannon", MARGIN = 1, base = exp(1))
shd_GRSMfungi <- diversity(fungi_matrix_GRSM, index = "shannon", MARGIN = 1, base = exp(1))
shd_WOODfungi <- diversity(fungi_matrix_WOOD, index = "shannon", MARGIN = 1, base = exp(1))
shd_OAESfungi <- diversity(fungi_matrix_OAES, index = "shannon", MARGIN = 1, base = exp(1))
shd_YELLfungi <- diversity(fungi_matrix_YELL, index = "shannon", MARGIN = 1, base = exp(1))
shd_NIWOfungi <- diversity(fungi_matrix_NIWO, index = "shannon", MARGIN = 1, base = exp(1))
shd_SRERfungi <- diversity(fungi_matrix_SRER, index = "shannon", MARGIN = 1, base = exp(1))
shd_ONAQfungi <- diversity(fungi_matrix_ONAQ, index = "shannon", MARGIN = 1, base = exp(1))
shd_ABBYfungi <- diversity(fungi_matrix_ABBY, index = "shannon", MARGIN = 1, base = exp(1))
shd_SJERfungi <- diversity(fungi_matrix_SJER, index = "shannon", MARGIN = 1, base = exp(1))
shd_BARRfungi <- diversity(fungi_matrix_BARR, index = "shannon", MARGIN = 1, base = exp(1))

rich_GUANfungi <- specnumber(fungi_matrix_GUAN, MARGIN = 1)
rich_HEALfungi <- specnumber(fungi_matrix_HEAL, MARGIN = 1)
rich_HARVfungi <- specnumber(fungi_matrix_HARV, MARGIN = 1)
rich_CPERfungi <- specnumber(fungi_matrix_CPER, MARGIN = 1)
rich_TALLfungi <- specnumber(fungi_matrix_TALL, MARGIN = 1)
rich_SCBIfungi <- specnumber(fungi_matrix_SCBI, MARGIN = 1)
rich_OSBSfungi <- specnumber(fungi_matrix_OSBS, MARGIN = 1)
rich_UNDEfungi <- specnumber(fungi_matrix_UNDE, MARGIN = 1)
rich_UKFSfungi <- specnumber(fungi_matrix_UKFS, MARGIN = 1)
rich_GRSMfungi <- specnumber(fungi_matrix_GRSM, MARGIN = 1)
rich_WOODfungi <- specnumber(fungi_matrix_WOOD, MARGIN = 1)
rich_OAESfungi <- specnumber(fungi_matrix_OAES, MARGIN = 1)
rich_YELLfungi <- specnumber(fungi_matrix_YELL, MARGIN = 1)
rich_NIWOfungi <- specnumber(fungi_matrix_NIWO, MARGIN = 1)
rich_SRERfungi <- specnumber(fungi_matrix_SRER, MARGIN = 1)
rich_ONAQfungi <- specnumber(fungi_matrix_ONAQ, MARGIN = 1)
rich_ABBYfungi <- specnumber(fungi_matrix_ABBY, MARGIN = 1)
rich_SJERfungi <- specnumber(fungi_matrix_SJER, MARGIN = 1)
rich_BARRfungi <- specnumber(fungi_matrix_BARR, MARGIN = 1)


### GUAN data table with shannon and diversity
plotID <- c("GUAN_001","GUAN_002","GUAN_003","GUAN_004","GUAN_006","GUAN_007",
            "GUAN_042","GUAN_043","GUAN_048","GUAN_049")
siteID <- rep("GUAN", length(plotID))
species <- rep("fungi", length(plotID))
x_GUAN_fungi <- data.frame(plotID, siteID, shd_GUANfungi, rich_GUANfungi,
                              species)
colnames(x_GUAN_fungi)<- c("plotID","siteID", "shd", "rich","species")

### HEAL data table with shannon and diversity
plotID <- c("HEAL_001","HEAL_004","HEAL_005","HEAL_010","HEAL_011","HEAL_025",
            "HEAL_045","HEAL_046","HEAL_047","HEAL_048")
siteID <- rep("HEAL", length(plotID))
species <- rep("fungi", length(plotID))
x_HEAL_fungi <- data.frame(plotID, siteID, shd_HEALfungi, rich_HEALfungi,
                              species)
colnames(x_HEAL_fungi)<- c("plotID","siteID", "shd", "rich","species")

### HARV data table with shannon and diversity
plotID <- c("HARV_001","HARV_002","HARV_004","HARV_005","HARV_010","HARV_013",
            "HARV_016","HARV_020","HARV_021","HARV_033","HARV_034","HARV_035",
            "HARV_037")
siteID <- rep("HARV", length(plotID))
species <- rep("fungi", length(plotID))
x_HARV_fungi <- data.frame(plotID, siteID, shd_HARVfungi, rich_HARVfungi,
                              species)
colnames(x_HARV_fungi)<- c("plotID","siteID", "shd", "rich","species")

### CPER data table with shannon and diversity
plotID <- c("CPER_001","CPER_002","CPER_003","CPER_004","CPER_005","CPER_006",
            "CPER_045","CPER_046","CPER_047","CPER_048")
siteID <- rep("CPER", length(plotID))
species <- rep("fungi", length(plotID))
x_CPER_fungi <- data.frame(plotID, siteID, shd_CPERfungi, rich_CPERfungi,
                              species)
colnames(x_CPER_fungi)<- c("plotID","siteID", "shd", "rich","species")

### TALL data table with shannon and diversity
plotID <- c("TALL_001","TALL_002","TALL_003","TALL_004","TALL_006","TALL_007",
            "TALL_044","TALL_047","TALL_051","TALL_054")
siteID <- rep("TALL", length(plotID))
species <- rep("fungi", length(plotID))
x_TALL_fungi <- data.frame(plotID, siteID, shd_TALLfungi, rich_TALLfungi,
                              species)
colnames(x_TALL_fungi)<- c("plotID","siteID", "shd", "rich","species")

### SCBI data table with shannon and diversity
plotID <- c("SCBI_002","SCBI_003","SCBI_004","SCBI_005","SCBI_006","SCBI_008",
            "SCBI_012","SCBI_045","SCBI_046","SCBI_047","SCBI_049","SCBI_067")
siteID <- rep("SCBI", length(plotID))
species <- rep("fungi", length(plotID))
x_SCBI_fungi <- data.frame(plotID, siteID, shd_SCBIfungi, rich_SCBIfungi,
                              species)
colnames(x_SCBI_fungi)<- c("plotID","siteID", "shd", "rich","species")

### OSBS data table with shannon and diversity
plotID <- c("OSBS_001","OSBS_002","OSBS_003","OSBS_004","OSBS_005","OSBS_022",
            "OSBS_023","OSBS_026","OSBS_027","OSBS_029","OSBS_031")
siteID <- rep("OSBS", length(plotID))
species <- rep("fungi", length(plotID))
x_OSBS_fungi <- data.frame(plotID, siteID, shd_OSBSfungi, rich_OSBSfungi,
                              species)
colnames(x_OSBS_fungi)<- c("plotID","siteID", "shd", "rich","species")

### UNDE data table with shannon and diversity
plotID <- c("UNDE_001","UNDE_002","UNDE_003","UNDE_006","UNDE_007","UNDE_008",
            "UNDE_010","UNDE_013","UNDE_014","UNDE_017","UNDE_019","UNDE_027",
            "UNDE_034","UNDE_035","UNDE_037","UNDE_038","UNDE_043","UNDE_044")
siteID <- rep("UNDE", length(plotID))
species <- rep("fungi", length(plotID))
x_UNDE_fungi <- data.frame(plotID, siteID, shd_UNDEfungi, rich_UNDEfungi,
                              species)
colnames(x_UNDE_fungi)<- c("plotID","siteID", "shd", "rich","species")

### UKFS data table with shannon and diversity
plotID <- c("UKFS_001","UKFS_002","UKFS_003","UKFS_004","UKFS_005","UKFS_009",
            "UKFS_031","UKFS_032","UKFS_043","UKFS_044")
siteID <- rep("UKFS", length(plotID))
species <- rep("fungi", length(plotID))
x_UKFS_fungi <- data.frame(plotID, siteID, shd_UKFSfungi, rich_UKFSfungi,
                              species)
colnames(x_UKFS_fungi)<- c("plotID","siteID", "shd", "rich","species")

### GRSM data table with shannon and diversity
plotID <- c("GRSM_001","GRSM_002","GRSM_003","GRSM_006","GRSM_007","GRSM_016",
            "GRSM_055","GRSM_058","GRSM_059","GRSM_060")
siteID <- rep("GRSM", length(plotID))
species <- rep("fungi", length(plotID))
x_GRSM_fungi <- data.frame(plotID, siteID, shd_GRSMfungi, rich_GRSMfungi,
                              species)
colnames(x_GRSM_fungi)<- c("plotID","siteID", "shd", "rich","species")

### WOOD data table with shannon and diversity
plotID <- c("WOOD_001","WOOD_002","WOOD_003","WOOD_004","WOOD_005","WOOD_022",
            "WOOD_024","WOOD_042","WOOD_043","WOOD_044","WOOD_045")
siteID <- rep("WOOD", length(plotID))
species <- rep("fungi", length(plotID))
x_WOOD_fungi <- data.frame(plotID, siteID, shd_WOODfungi, rich_WOODfungi,
                              species)
colnames(x_WOOD_fungi)<- c("plotID","siteID", "shd", "rich","species")

### OAES data table with shannon and diversity
plotID <- c("OAES_001","OAES_002","OAES_003","OAES_004","OAES_007","OAES_009",
            "OAES_042","OAES_043","OAES_044","OAES_045")
siteID <- rep("OAES", length(plotID))
species <- rep("fungi", length(plotID))
x_OAES_fungi <- data.frame(plotID, siteID, shd_OAESfungi, rich_OAESfungi,
                              species)
colnames(x_OAES_fungi)<- c("plotID","siteID", "shd", "rich","species")

### YELL data table with shannon and diversity
plotID <- c("YELL_001","YELL_002","YELL_003","YELL_009","YELL_012","YELL_016",
            "YELL_046","YELL_048","YELL_051","YELL_052")
siteID <- rep("YELL", length(plotID))
species <- rep("fungi", length(plotID))
x_YELL_fungi <- data.frame(plotID, siteID, shd_YELLfungi, rich_YELLfungi,
                              species)
colnames(x_YELL_fungi)<- c("plotID","siteID", "shd", "rich","species")

### NIWO data table with shannon and diversity
plotID <- c("NIWO_001","NIWO_002","NIWO_003","NIWO_004","NIWO_005","NIWO_006",
            "NIWO_008","NIWO_040","NIWO_041","NIWO_042","NIWO_043")
siteID <- rep("NIWO", length(plotID))
species <- rep("fungi", length(plotID))
x_NIWO_fungi <- data.frame(plotID, siteID, shd_NIWOfungi, rich_NIWOfungi,
                              species)
colnames(x_NIWO_fungi)<- c("plotID","siteID", "shd", "rich","species")

### SRER data table with shannon and diversity
plotID <- c("SRER_001","SRER_002","SRER_003","SRER_004","SRER_005","SRER_006",
            "SRER_043","SRER_047","SRER_052","SRER_053")
siteID <- rep("SRER", length(plotID))
species <- rep("fungi", length(plotID))
x_SRER_fungi <- data.frame(plotID, siteID, shd_SRERfungi, rich_SRERfungi,
                              species)
colnames(x_SRER_fungi)<- c("plotID","siteID", "shd", "rich","species")

### ONAQ data table with shannon and diversity
plotID <- c("ONAQ_002","ONAQ_003","ONAQ_004","ONAQ_005","ONAQ_006","ONAQ_007",
            "ONAQ_008","ONAQ_009","ONAQ_010","ONAQ_012","ONAQ_017","ONAQ_041",
            "ONAQ_042","ONAQ_043","ONAQ_044")
siteID <- rep("ONAQ", length(plotID))
species <- rep("fungi", length(plotID))
x_ONAQ_fungi <- data.frame(plotID, siteID, shd_ONAQfungi, rich_ONAQfungi,
                              species)
colnames(x_ONAQ_fungi)<- c("plotID","siteID", "shd", "rich","species")

### ABBY data table with shannon and diversity
plotID <- c("ABBY_001","ABBY_002","ABBY_003","ABBY_004","ABBY_006","ABBY_023",
            "ABBY_061","ABBY_062","ABBY_063","ABBY_070")
siteID <- rep("ABBY", length(plotID))
species <- rep("fungi", length(plotID))
x_ABBY_fungi <- data.frame(plotID, siteID, shd_ABBYfungi, rich_ABBYfungi,
                              species)
colnames(x_ABBY_fungi)<- c("plotID","siteID", "shd", "rich","species")

### SJER data table with shannon and diversity
plotID <- c("SJER_001","SJER_002","SJER_003","SJER_004","SJER_005","SJER_025",
            "SJER_045","SJER_046","SJER_047","SJER_048")
siteID <- rep("SJER", length(plotID))
species <- rep("fungi", length(plotID))
x_SJER_fungi <- data.frame(plotID, siteID, shd_SJERfungi, rich_SJERfungi,
                              species)
colnames(x_SJER_fungi)<- c("plotID","siteID", "shd", "rich","species")

### BARR data table with shannon and diversity
plotID <- c("BARR_001","BARR_002","BARR_003","BARR_004","BARR_005","BARR_006",
            "BARR_051","BARR_052","BARR_053","BARR_054")
siteID <- rep("BARR", length(plotID))
species <- rep("fungi", length(plotID))
x_BARR_fungi <- data.frame(plotID, siteID, shd_BARRfungi, rich_BARRfungi,
                              species)
colnames(x_BARR_fungi)<- c("plotID","siteID", "shd", "rich","species")

# combine data tables from each site
x_fungi <- bind_rows(x_GUAN_fungi,x_HEAL_fungi,x_HARV_fungi,
                     x_CPER_fungi,x_TALL_fungi,x_SCBI_fungi,x_OSBS_fungi,
                     x_UNDE_fungi,x_UKFS_fungi,x_GRSM_fungi,x_WOOD_fungi,
                     x_OAES_fungi,x_YELL_fungi,x_NIWO_fungi,x_SRER_fungi,
                     x_ONAQ_fungi,x_ABBY_fungi,x_SJER_fungi,x_BARR_fungi)




# combine plant and microbe data
x2 <- bind_rows(x_fungi,x_plants)


# make a new data table with both plant and microbe data in each row by plotID
merged_x2 <- merge(x_fungi, x_plants, by = "plotID")

# clean up names of columns
setnames(merged_x2, old = "rich.x", new = "rich.fun")
setnames(merged_x2, old = "rich.y", new = "rich.pla")
setnames(merged_x2, old = "shd.x", new = "shd.fun")
setnames(merged_x2, old = "shd.y", new = "shd.pla")
setnames(merged_x2, old = "siteID.x", new = "siteID")
setDT(merged_x2)
merged_x2 <- merged_x2[-14]

# add soil data
x2_temp <- merge(merged_x2, avg_temp, by = "siteID")
x2_soil <- merge(x2_temp, avg_moist, by = "siteID")
x_fun_temp <- merge(x_fungi, avg_temp, by = "siteID")

# make a graph that compares plant and fungi species richness by mean temp
ggplot(x2_temp, aes(x = rich.fun, y = rich.pla, color = mean_temp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "gold", high = "red2") +
  labs(
    title = "Relationship Between Fungi and Plant Richness",
    x = "Fungi Richness",
    y = "Plant Richness",
    color = "Soil Temp (˚C)"
  ) +
  theme_minimal()

# run statistics on richness graph
model <- lm(rich.fun ~ rich.pla + mean_temp + rich.pla*mean_temp, data = x2_temp)
summary(model)

# split this table into three categories by temperature: low, medium and high
x2_temp <- x2_temp %>%
  mutate(temp_group = case_when(
    mean_temp <= quantile(mean_temp, 0.33) ~ "Low Temp",
    mean_temp <= quantile(mean_temp, 0.66) ~ "Medium Temp",
    TRUE ~ "High Temp"
  ))

# graph the table to show correlations by each temperature group
ggplot(x2_temp, aes(x = rich.fun, y = rich.pla, color = temp_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Low Temp" = "purple", "Medium Temp" = "orange", "High Temp" = "red")) +
  labs(
    title = "Relationship Between Fungi and Plant Richness",
    x = "Fungi Richness",
    y = "Plant Richness",
    color = "Temperature"
  ) +
  theme_minimal()

# statistics on richness graph by soil temperature group
model <- lm(rich.fun ~ rich.pla + temp_group + rich.pla*temp_group, data = x2_temp)
summary(model)
model <- lm(rich.fun ~ rich.pla, data = subset(x2_temp, temp_group == "Low Temp"))
summary(model)
model <- lm(rich.fun ~ rich.pla, data = subset(x2_temp, temp_group == "Medium Temp"))
summary(model)
model <- lm(rich.fun ~ rich.pla, data = subset(x2_temp, temp_group == "High Temp"))
summary(model)


# make a graph that compares plant and microbial SHD by site ID
ggplot(x2_temp, aes(x = shd.fun, y = shd.pla, color = mean_temp)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "gold", high = "red2") +
  labs(
    title = "Relationship Between Fungi and Plant Shannon",
    x = "Fungi Shannon",
    y = "Plant Shannon",
    color = "Soil Temp (˚C)"
  ) +
  theme_minimal()


# run statistics on shannon graph
model <- lm(shd.fun ~ shd.pla + mean_temp + shd.pla*mean_temp, data = x2_temp)
summary(model)


# split this table into three categories by moisture: low, medium and high
x2_soil <- x2_soil %>%
  mutate(moist_group = case_when(
    mean_moist <= quantile(mean_moist, 0.33) ~ "Dry",
    mean_moist <= quantile(mean_moist, 0.66) ~ "Medium",
    TRUE ~ "Wet"
  ))

# graph the table to show correlations by each temperature group
ggplot(x2_soil, aes(x = rich.fun, y = rich.pla, color = moist_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Dry" = "skyblue1", "Medium" = "royalblue", "Wet" = "navy")) +
  labs(
    title = "Relationship Between Fungi and Plant Richness",
    x = "Fungi Richness",
    y = "Plant Richness",
    color = "Moisture"
  ) +
  theme_minimal()


# statistics on richness graph by soil temperature group
model <- lm(rich.fun ~ rich.pla + moist_group + rich.pla*moist_group, data = x2_soil)
summary(model)
model <- lm(rich.fun ~ rich.pla, data = subset(x2_soil, moist_group == "Dry"))
summary(model)
model <- lm(rich.fun ~ rich.pla, data = subset(x2_soil, moist_group == "Medium"))
summary(model)
model <- lm(rich.fun ~ rich.pla, data = subset(x2_soil, moist_group == "Wet"))
summary(model)




#merge fungi, plant, and soil data (long)
fun_temp <- merge(x2, avg_temp, by = "siteID")
fun_soil <- merge(fun_temp, avg_moist, by = "siteID")
setDT(fun_soil)
fun_soil <- fun_soil[order(plotID)]
fun_soil$mean_temp <- factor(fun_soil$mean_temp)
fun_soil$mean_moist <- factor(fun_soil$mean_moist)


# graph shannon and diversity as boxplots by soil temperature rather than site/plot
ggplot(fun_soil, aes(x = mean_temp, y = shd, fill = species)) +
  geom_boxplot(position = position_dodge()) +
  labs(
    title = "Shannon Diversity Index by Soil Temperature",
    x = "Soil Temperature (˚C)",
    y = "Shannon Diversity"
  ) +
  theme_classic()

model <- lm(shd ~ mean_temp, data = fun_soil)
summary(model)

ggplot(fun_soil, aes(x = mean_temp, y = rich, fill = species)) +
  geom_boxplot(position = position_dodge()) +
  labs(
    title = "Species Richness by Soil Temperature",
    x = "Soil Temperature (˚C)",
    y = "Species Richness"
  ) +
  theme_classic()

model <- lm(rich ~ mean_temp, data = fun_soil)
summary(model)


model <- lm(shd.pla + mean_moist + mean_temp ~ shd.fun, data = x2_soil)
summary(model)


# make a graph that compares plant and fungi species richness by soil moisture
ggplot(x2_soil, aes(x = rich.fun, y = rich.pla, color = mean_moist)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Relationship Between Fungi and Plant Richness",
    x = "Fungi Richness",
    y = "Plant Richness",
    color = "Soil Moisture
(g water/g dry soil)"
  ) +
  theme_minimal()

# run statistics on richness graph
model <- lm(rich.fun ~ rich.pla + mean_moist + rich.pla*mean_moist, data = x2_soil)
summary(model)


# make a graph that compares plant and microbial SHD by soil moisture
ggplot(x2_soil, aes(x = shd.fun, y = shd.pla, color = mean_moist)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(
    title = "Relationship Between Fungi and Plant Shannon",
    x = "Fungi Shannon",
    y = "Plant Shannon",
    color = "Soil Moisture
(g water/g dry soil)"
  ) +
  theme_minimal()


# run statistics on shannon graph
model <- lm(shd.fun ~ shd.pla + mean_moist + shd.pla*mean_moist, data = x2_soil)
summary(model)


# graph the table to show correlations by each temperature group
ggplot(x2_soil, aes(x = shd.fun, y = shd.pla, color = moist_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Dry" = "skyblue1", "Medium" = "royalblue", "Wet" = "navy")) +
  labs(
    title = "Relationship Between Fungi and Plant Shannon",
    x = "Fungi Shannon",
    y = "Plant Shannon",
    color = "Moisture"
  ) +
  theme_minimal()

# statistics on shannon graph by soil moisture group
model <- lm(shd.fun ~ shd.pla + moist_group + shd.pla*moist_group, data = x2_soil)
summary(model)
model <- lm(shd.fun ~ shd.pla, data = subset(x2_soil, moist_group == "Dry"))
summary(model)
model <- lm(shd.fun ~ shd.pla, data = subset(x2_soil, moist_group == "Medium"))
summary(model)
model <- lm(shd.fun ~ shd.pla, data = subset(x2_soil, moist_group == "Wet"))
summary(model)


# graph the table to show correlations by each temperature group
ggplot(x2_temp, aes(x = shd.fun, y = shd.pla, color = temp_group)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("Low Temp" = "purple", "Medium Temp" = "orange", "High Temp" = "red")) +
  labs(
    title = "Relationship Between Fungi and Plant Shannon",
    x = "Fungi Shannon",
    y = "Plant Shannon",
    color = "Temperature"
  ) +
  theme_minimal()

# statistics on shannon graph by soil moisture group
model <- lm(shd.fun ~ shd.pla + temp_group + shd.pla*temp_group, data = x2_temp)
summary(model)
model <- lm(shd.fun ~ shd.pla, data = subset(x2_temp, temp_group == "Low Temp"))
summary(model)
model <- lm(shd.fun ~ shd.pla, data = subset(x2_temp, temp_group == "Medium Temp"))
summary(model)
model <- lm(shd.fun ~ shd.pla, data = subset(x2_temp, temp_group == "High Temp"))
summary(model)

