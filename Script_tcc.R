
## ------------------------------------------- ##
##               Script para o TCC             ##
##                                             ##
##              Construção de árvore           ##
##                 de genes                    ##
##                                             ##
##                                             ##
## ------------------------------------------- ## 


# Alinhamento das sequências de D.melanogaster

R

arquivo_entrada <- "/mnt/datos/Pasantia_Lais/nova_abordagem/Ir-Dm-prot_filtrado.fasta"
arquivo_saida <- "alinhamento_irs_dm"

comando_mafft <- "mafft"
args <- c("--auto", "--maxiterate", "1000", arquivo_entrada, ">", arquivo_saida)
system2(comando_mafft, args = args)

#71seqs

# Alinhamento das sequências de C.hominivorax

arquivo_entrada1 <- "/mnt/datos/Pasantia_Lais/nova_abordagem/Ir-Ch88-prot.fasta"
arquivo_saida1 <- "alinhamento_irs_ch"

comando_mafft <- "mafft"
args <- c("--auto", "--maxiterate", "1000", arquivo_entrada1, ">", arquivo_saida1)
system2(comando_mafft, args = args)

#88 seqs


# Alinhamento das sequências de M.domestica

arquivo_entrada2 <- "/mnt/datos/Pasantia_Lais/nova_abordagem/irs_mdom.fasta"
arquivo_saida2 <- "alinhamento_irs_md"

comando_mafft <- "mafft"
args <- c("--auto", "--maxiterate", "1000", arquivo_entrada2, ">", arquivo_saida2)
system2(comando_mafft, args = args)


# Alinhamento das sequências de Anopheles

arquivo_entrada3 <- "/mnt/datos/Pasantia_Lais/nova_abordagem/irs_Ano.fasta"
arquivo_saida3 <- "alinhamento_irs_ano"

comando_mafft <- "mafft"
args <- c("--auto", "--maxiterate", "1000", arquivo_entrada3, ">", arquivo_saida3)
system2(comando_mafft, args = args)


# Alinhamento das sequências de Scal

arquivo_entrada4 <- "/mnt/datos/Pasantia_Lais/nova_abordagem/irs_Scal.fasta"
arquivo_saida4 <- "alinhamento_irs_scal"

comando_mafft <- "mafft"
args <- c("--auto", "--maxiterate", "1000", arquivo_entrada4, ">", arquivo_saida4)
system2(comando_mafft, args = args)

# --------------------------------------------------------------------------------------------------

# Mandar os arquivos dos irs de glossina (antenais e divergentes) para o outro computador

#rsync -avzu arquivo.txt usuario@ip_do_destino:/caminho/de/destino/

rsync -avzu Antennal_irs.aln lcardoso@ubi12:/mnt/datos/Pasantia_Lais/nova_abordagem
rsync -avzu Divergent_IRs.aln lcardoso@ubi12:/mnt/datos/Pasantia_Lais/nova_abordagem

#Transformar o arquivo clustal em fasta 
R
setwd("/mnt/datos/Pasantia_Lais/nova_abordagem")
getwd()

library(Biostrings)

#Antenais
aln_file <- "/mnt/datos/Pasantia_Lais/nova_abordagem/Antennal_irs.aln" # Substitua pelo caminho do arquivo .aln
alignment <- readAAMultipleAlignment(filepath = aln_file, format = "clustal")

fasta_file <- "Antennal_irs.fasta" # Nome do arquivo de saída em FASTA
writeXStringSet(as(alignment, "AAStringSet"), filepath = fasta_file, format = "fasta")

#Divergentes

aln_file <- "/mnt/datos/Pasantia_Lais/nova_abordagem/Divergent_IRs.aln" # Substitua pelo caminho do arquivo .aln
alignment <- readAAMultipleAlignment(filepath = aln_file, format = "clustal")

fasta_file <- "Divergent_irs.fasta" # Nome do arquivo de saída em FASTA
writeXStringSet(as(alignment, "AAStringSet"), filepath = fasta_file, format = "fasta")

# Vou pegar esses arquivos .fasta e pedir que deixe apenas aquelas de glossina
# Depois eu vou juntar esses dois arquivos para deixar tudo de glossina

fasta_file <- "Antennal_irs.fasta"
sequences <- readAAStringSet(fasta_file, format = "fasta") # Lê o arquivo FASTA
filtered_sequences <- sequences[grepl("^G", names(sequences))] # Filtra os cabeçalhos que começam com "G"
writeXStringSet(filtered_sequences, "Antenais_glossina.fasta", format = "fasta") # Salva o arquivo com as sequências filtradas


fasta_file <- "Divergent_irs.fasta"
sequences <- readAAStringSet(fasta_file, format = "fasta") 
filtered_sequences <- sequences[grepl("^G", names(sequences))] 
writeXStringSet(filtered_sequences, "Divergentes_glossina.fasta", format = "fasta") 

antenais <- "Antenais_glossina.fasta"
divergentes <- "Divergentes_glossina.fasta"

seqs1 <- readAAStringSet(antenais, format = "fasta")
seqs2 <- readAAStringSet(divergentes, format = "fasta")

merged_seqs <- c(seqs1, seqs2)
writeXStringSet(merged_seqs, "glossina_irs.fasta", format = "fasta")

# Juntar todos os arquivos de todas as espécies com o intuito de alinhar 
# movi todos os arquivos das seqs para uma nova pasta em bash 
cp glossina_irs.fasta /mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs
cp Ir-Ch88-prot.fasta /mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs
cp Ir-Dm-prot_filtrado.fasta /mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs
cp irs_Ano.fasta /mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs
cp irs_mdom.fasta /mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs
cp irs_Scal.fasta /mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs

pasta_fasta <- "/mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs" 
arquivos <- list.files(path = pasta_fasta, pattern = "\\.fasta$", full.names = TRUE)

todas_sequencias <- do.call(c, lapply(arquivos, readAAStringSet))
writeXStringSet(todas_sequencias, "all_irs.fasta", format = "fasta")
# ---------------------------------------------------------------------------------------


# Alinhamento 

# Arquivo de entrada com sequências de aminoácidos
arquivo_entrada <- "all_irs.fasta"

# Arquivo de saída com alinhamento
arquivo_saida <- "alinhamento.fasta"

comando_mafft <- "mafft"
args <- c("--auto", "--anysymbol", "--maxiterate", "1000", arquivo_entrada, ">", arquivo_saida)
system2(comando_mafft, args = args)


# Construção da árvore - IQ-TREE


# caminho do executável do IQ-TREE # /mnt/datos/programas/iqtree-1.6.12-Linux/bin/iqtree
iqtree_path <- "/mnt/datos/programas/iqtree-1.6.12-Linux/bin/iqtree"  
input_file <- "/mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs/alinhamento.fasta" # Caminho para o arquivo de entrada

# Comando para rodar o IQ-TREE
system(paste(iqtree_path, "-s", input_file, "-bb", 1000, "-alrt", 1000)) 


# Construção da Tabela metadados
##    ---------------------------

setwd("/mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs")
library(ggtree)

carrega a árvore 
#tree <- read.tree("sua_arvore.treefile")
arvore <- read.tree("/mnt/datos/Pasantia_Lais/nova_abordagem/sequencias_irs/alinhamento.fasta.treefile")

tip_labels <- arvore$tip.label
tip_data <- data.frame(TipLabel = tip_labels)
write.csv(tip_data, "metadados.csv", row.names = FALSE)


# Edição da árvore de genes
##    ---------------------------

# R studio


setwd("/home/lais-cardoso-valadao/Area_de_trabalho/Pasantia/irs")
getwd()
list.files()

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggplot2)
library(ggtree) # Pacote ggtree carregado
#library(dplyr)
arvore <- read.tree("alinhamento.fasta.treefile") # arquivo da árvore carregada
arvoreR <- phangorn::midpoint(arvore)
ggtree(arvoreR)
#class(arvore)
metadados <- read.csv("metadados.csv")
metadados
# Tamanho da fonte 
ggtree(arvoreR) + 
  geom_tiplab(size = 1)

# aumentar ou diminuir o espaçamento entre os rótulos e os ramos da árvore
ggtree(arvoreR) +
  geom_tiplab(size = 2, offset = 0.1) +  # Aumenta o espaço entre o rótulo e o ramo
  theme_tree2()

# open.angle: Define o ângulo de abertura da árvore circular (valores como 45, 90, ou outros podem ser experimentados)

ggtree(arvoreR, layout = "circular", open.angle = 180) + 
  geom_tiplab(size = 1, offset = 0.5, align = TRUE) +
  theme_tree2()

# Adiciona pontos e cores nos ramos da árvore 
gg <- ggtree(arvoreR, layout = "circular")
gg <- gg %<+% metadados
gg <- gg + geom_tippoint(aes(color = factor(Gene)))
#gg <- gg + geom_tippoint(aes(color = (Gene)))
gg

gg + theme(legend.position = "none")  # Remove a legenda se não for essencial

# Criando o gráfico com cores manuais
gg <- ggtree(arvoreR, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Gene))) +
  scale_color_manual(values = c( 
    "Ir8a" = "#32CD32",
    "Ir25a" = "#303F9F",
    "Ir21a" = "#E54100",
    "Ir93a" = "#FF69B4",
    "Ir76b" = "#D32F2F")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Árvore de Nematocera e Brachycera

library(ape)
library(ggtree)
#install.packages("phytools")
library(phytools)

gg <- ggtree(arvoreR, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Especie))) +
  geom_nodelab(aes(label = label), hjust = 1.5, vjust = -0.5, size = 1) +  # Ajuste do tamanho da fonte
  scale_color_manual(values = c( 
    "Glossina" = "#ffad08",
    "Cochliomyia hominivorax" = "#edd75a",
    "Stomoxys calcitrans" = "#f01880",
    "Drosophila melanogaster" = "#73b06f",
    "Musca domestica" ="#0c8f8f",
    "Aedes aegypti" ="#3a4c38"
  )) +
  theme(legend.position = "right")  

print(gg)

ggr <- ggtree(arvoreR, layout = "daylight") %<+% metadados +
  geom_tippoint(aes(color = factor(Gene))) +
  scale_color_manual(values = c( 
    "Ir8a" = "#32CD32",
    "Ir25a" = "#303F9F",
    "Ir21a" = "#E54100",
    "Ir93a" = "#FF69B4",
    "Ir76b" = "#D32F2F")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(ggr)


# Antenais e Divergentes
gg <- ggtree(arvoreR, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Tipcolors))) +
  scale_color_manual(values = c( 
    "green" = "#b50a04",
    "blue" = "#f58207",
    "lightblue" = "#16a0c9", 
    "black" = "#000000",
    "grey" = "#8a9294")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Criação da árvore sem nematocera
"AaegIr7a","AaegIr103","AaegIr41h", "AaegIr75b", "AaegIr75c",
"AaegIr7k","AaegIr127","AaegIr112","AaegIr133","AaegIr139",
"AaegIr7e","AaegIr123","AaegIr41c","AaegIr161","AaegIr41m",
"AaegIr132","AaegIr145","AaegIr146","AaegIr41p","AaegIr148",
"AaegIr149", "AaegIr150", "AaegIr75a", "AaegIr75k", "AaegIr122",
"AaegIr31a1","AaegIr128", "AaegIr7l", "AaegIr93a", "AaegIr160",
"AaegIr104","AaegIr7m", "AaegIr64a", "AaegIr102", "AaegIr60a",
"AaegIr147","AaegIr21a", "AaegIr41g", "AaegIr118", "AaegIr7o", 
"AaegIr8a","AaegIr25a", "AaegIr75d", "AaegIr105", "AaegIr41n",
"AaegIr68a", "AaegIr106", "AaegIr7n", "AaegIr7p", "AaegIr7q",
"AaegIr7r", "AaegIr31a2", "AaegIr41l", "AaegIr41b", "AaegIr110",
"AaegIr136", "AaegIr121", "AaegIr131", "AaegIr164", "AaegIr143",
"AaegIr152", "AaegIr153", "AaegIr154", "AaegIr41e", "AaegIr75i",
"AaegIr75j","AaegIr107", "AaegIr41f", "AaegIr157","AaegIr41d", 
"AaegIr41i", "AaegIr75e","AaegIr75f","AaegIr75g", "AaegIr138",
"AaegIr142","AaegIr144", "AaegIr141", "AaegIr156", "AaegIr41o",
"AaegIr124","AaegIr125", "AaegIr108", "AaegIr163","AaegIr111",
"AaegIr130", "AaegIr155","AaegIr7b", "AaegIr169", "AaegIr170", 
"AaegIr140", "AaegIr7j","AaegIr117", "AaegIr100c", "AaegIr41a",
"AaegIr100a", "AaegIr41j", "AaegIr113", "AaegIr100b", "AaegIr7c",
"AaegIr7d", "AaegIr167","AaegIr168","AaegIr7h","AaegIr171","AaegIr115",
"AaegIr75l","AaegIr126","AaegIr120", "AaegIr137", "AaegIr7g",
"AaegIr151","AaegIr7i","AaegIr40a","AaegIr135","AaegIr172",
"AaegIr75h","AaegIr100d","AaegIr116","AaegIr165","AaegIr76b",
"AaegIr41k", "AaegIr166", "AaegIr87a2", "AaegIr134","AaegIr119",
"AaegIr7f","AaegIr114","AaegIr158","AaegIr159","AaegIr162","AaegIr101",
"AaegIr87a1","AaegIr129", "AaegIr109" 
library(ape)
mrca_node <- getMRCA(arvoreR, c("AaegIr7a","AaegIr103","AaegIr41h", "AaegIr75b", "AaegIr75c",
                                "AaegIr7k","AaegIr127","AaegIr112","AaegIr133","AaegIr139",
                                "AaegIr7e","AaegIr123","AaegIr41c","AaegIr161","AaegIr41m",
                                "AaegIr132","AaegIr145","AaegIr146","AaegIr41p","AaegIr148",
                                "AaegIr149", "AaegIr150", "AaegIr75a", "AaegIr75k", "AaegIr122",
                                "AaegIr31a1","AaegIr128", "AaegIr7l", "AaegIr93a", "AaegIr160",
                                "AaegIr104","AaegIr7m", "AaegIr64a", "AaegIr102", "AaegIr60a",
                                "AaegIr147","AaegIr21a", "AaegIr41g", "AaegIr118", "AaegIr7o", 
                                "AaegIr8a","AaegIr25a", "AaegIr75d", "AaegIr105", "AaegIr41n",
                                "AaegIr68a", "AaegIr106", "AaegIr7n", "AaegIr7p", "AaegIr7q",
                                "AaegIr7r", "AaegIr31a2", "AaegIr41l", "AaegIr41b", "AaegIr110",
                                "AaegIr136", "AaegIr121", "AaegIr131", "AaegIr164", "AaegIr143",
                                "AaegIr152", "AaegIr153", "AaegIr154", "AaegIr41e", "AaegIr75i",
                                "AaegIr75j","AaegIr107", "AaegIr41f", "AaegIr157","AaegIr41d", 
                                "AaegIr41i", "AaegIr75e","AaegIr75f","AaegIr75g", "AaegIr138",
                                "AaegIr142","AaegIr144", "AaegIr141", "AaegIr156", "AaegIr41o",
                                "AaegIr124","AaegIr125", "AaegIr108", "AaegIr163","AaegIr111",
                                "AaegIr130", "AaegIr155","AaegIr7b", "AaegIr169", "AaegIr170", 
                                "AaegIr140", "AaegIr7j","AaegIr117", "AaegIr100c", "AaegIr41a",
                                "AaegIr100a", "AaegIr41j", "AaegIr113", "AaegIr100b", "AaegIr7c",
                                "AaegIr7d", "AaegIr167","AaegIr168","AaegIr7h","AaegIr171","AaegIr115",
                                "AaegIr75l","AaegIr126","AaegIr120", "AaegIr137", "AaegIr7g",
                                "AaegIr151","AaegIr7i","AaegIr40a","AaegIr135","AaegIr172",
                                "AaegIr75h","AaegIr100d","AaegIr116","AaegIr165","AaegIr76b",
                                "AaegIr41k", "AaegIr166", "AaegIr87a2", "AaegIr134","AaegIr119",
                                "AaegIr7f","AaegIr114","AaegIr158","AaegIr159","AaegIr162","AaegIr101",
                                "AaegIr87a1","AaegIr129", "AaegIr109")) 
print(mrca_node)

# Define the node of the clade to prune
node_to_prune <- 671

# Get descendants (tips) of that node
desc_tips <- extract.clade(arvoreR, node_to_prune)$tip.label

# Prune the clade
pruned_tree <- drop.tip(arvoreR, desc_tips)

#metadados <- read.csv("metadados.csv")

gg <- ggtree(pruned_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Classe))) +
  scale_color_manual(values = c( 
  "Antenais D.melanogaster" = "#b50a04",
    "Divergente D.melanogaster" = "#f58207",
    "C.hominivorax" = "#16a0c9", 
    "Glossina" = "#b009ed",
    "M.domestica" = "#e6ed09",
  "S.calcitrans" = "#ed09c7",
  "Co-receptor D.melanogaster"= "#0911ed")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Raiz

"ScalIr25a2","ScalIr25a3N","ScalIr25a4N", "ScalIr25a1"
"MdIR25a","ChomIr25a-2","ChomIr25a-1","Dm_Ir25a_isoC",
"Dm_Ir25a_isoD","Dm_Ir25a_isoE","Dm_Ir25a_isoB", "GbrIr25a_1",
"GffIr25a_1","GaIr25a_16","GpdIr25a_1","GmmIr25a_9","GbrIr8a_1-",
"GffIr8a_12","GaIr8a_1-3","GpdIr8a_14","GmmIr8a_1-",
"ChomIr8a","MdIR8a","ScalIr8a","Dm_Ir8a"




mrca_node <- getMRCA(pruned_tree, c("ScalIr25a2","ScalIr25a3N","ScalIr25a4N", "ScalIr25a1",
                                "MdIR25a","ChomIr25a-2","ChomIr25a-1","Dm_Ir25a_isoC",
                                "Dm_Ir25a_isoD","Dm_Ir25a_isoE","Dm_Ir25a_isoB", "GbrIr25a_1",
                                "GffIr25a_1","GaIr25a_16","GpdIr25a_1","GmmIr25a_9","GbrIr8a_1-",
                                "GffIr8a_12","GaIr8a_1-3","GpdIr8a_14","GmmIr8a_1-",
                                "ChomIr8a","MdIR8a","ScalIr8a","Dm_Ir8a"))
  
print(mrca_node)

new_root_tree <- ape::root(pruned_tree, node = 493)

gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Classe))) +
  scale_color_manual(values = c( 
    "Antenais D.melanogaster" = "#b50a04",
    "Divergente D.melanogaster" = "#f58207",
    "C.hominivorax" = "#16a0c9", 
    "Glossina" = "#b009ed",
    "M.domestica" = "#e6ed09",
    "S.calcitrans" = "#ed09c7",
    "Co-receptor D.melanogaster"= "#0911ed")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Antenais
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Gene))) +
  scale_color_manual(values = c( 
    "Ir25a" = "#b50a04",
    "Ir8a" = "#b50a04",
    "Ir21a" ="#b50a04",
    "Ir31a" = "#b50a04",
    "Ir40a" = "#b50a04",
    "Ir64a" = "#b50a04",
    "Ir75a"="#b50a04",
    "Ir75b" = "#b50a04",
    "Ir75c" ="#b50a04",
    "Ir75d"="#b50a04",
    "Ir76a"="#b50a04",
    "Ir76b"="#b50a04",
    "Ir84a"="#b50a04",
    "Ir92a"="#b50a04",
    "Ir93a"="#b50a04",
    "Ir68a" ="#4a2d57",
    "Ir60a"="#2ab55f",
    "Ir41a"="#f5de14",
    "Ir67a"= "#1835f0"))+
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)



# Derivados
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Gene))) +
  scale_color_manual(values = c( 
    "Ir182" = "#b50a04",
    "Ir181" = "#731428",
    "Ir180" = "#cf2147", 
    "Ir179" = "#e65574",
    "Ir172" = "#70061d",
    "Ir171" = "#e0778e",
    "Ir173" = "#4110e3",
    "Ir168" = "#1f0966",
    "Ir178" = "#734cf5",
    "Ir169" = "#150059",
    "Ir170" = "#50eb26",
    "Ir166" = "#207508",
    "Ir174" = "#0a2e00",
    "Ir169" = "#81e366",
    "Ir164" = "#c9eb09",
    "Ir163" = "#e3ff47",
    "Ir162" = "#586318",
    "Ir160"= "#d9ff00",
    "Ir161"= "#a624c7",
    "Ir168" = "#592666",
    "Ir167" = "#dd88f2",
    "Ir172" = "#17011c",
    "Ir184" = "#471800",
    "Ir183" = "#f58953",
    "Ir170" = "#8a2e00",
    "Ir171" = "#ff5500"
    )) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Genes por espécie
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Especie))) +
  scale_color_manual(values = c( 
    "Glossina" = "#b50a04",
    "Cochliomyia hominivorax" = "#4110e3",
    "Musca domestica" = "#50eb26",
    "Drosophila melanogaster" = "#c9eb09",
    "Stomoxys calcitrans" = "#ff5500")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Musca domestica
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Especie))) +
  scale_color_manual(values = c( 
    "Musca domestica" = "#50eb26")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Drosophila melanogaster
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Especie))) +
  scale_color_manual(values = c( 
    "Drosophila melanogaster" = "#c9eb09" )) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)


# Dados só de D.melanogaster
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Classe))) +
  scale_color_manual(values = c( 
    "Antenais D.melanogaster" = "#b50a04",
    "Divergente D.melanogaster" = "#f58207",
    "Co-receptor D.melanogaster"= "#b50a04")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Stomoxys calcitrans
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Especie))) +
  scale_color_manual(values = c( 
    "Stomoxys calcitrans" = "#ff5500")) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)

# Genes 7 divergentes
gg <- ggtree(new_root_tree, layout = "circular") %<+% metadados +
  geom_tippoint(aes(color = factor(Gene))) +
  scale_color_manual(values = c( 
    "Ir7d" = "#b50a04",
    "Ir7e"= "#b50a04",
    "Ir7f" ="#b50a04",
    "Ir7g"= "#b50a04",
    "Ir7b" ="#b50a04",
    "Ir7c"="#b50a04",
    "Ir7a" ="#b50a04",
    "Ir7b" ="#b50a04",
    "Ir71h" ="#b50a04"
    )) +
  theme(legend.position = "right") # Ajuste a posição da legenda, se necessário
# Exibindo o gráfico
print(gg)










































