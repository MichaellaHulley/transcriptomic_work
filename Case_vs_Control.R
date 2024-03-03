## LOAD LIBRARIES PACKAGES
source('expression_libraries.R')

## OPTIONS
options(repr.matrix.max.rows = 100)
columns(org.Hs.eg.db)

## SET WORKING DIRECTORY (ASSUMING YOU HAVE EVERYTHING SETUP IN THE FOLDER YOU WORKING IN)
setwd(getwd())

## FOR HTSEQ ANALYSIS - COMMENT BELOW
counts.dir <- paste0(getwd(),"/expression_output")
if (file.exists(counts.dir)){
    in.file <- "output/4_Read_Counts/htseqCounts/gene_counts_final.txt"
    out.dir <- paste0(counts.dir,'/')
}else {
    dir.create(file.path(counts.dir), showWarnings = FALSE)
    in.file <- "output/4_Read_Counts/htseqCounts/gene_counts_final.txt"
    out.dir <- paste0(counts.dir,'/')
}

## PATH FOR FIGURES
figures.out.dir <- paste0(out.dir, "figures/")
if (file.exists(figures.out.dir)){
    system(paste0('rm -r ','"',figures.out.dir,'"',"*"))
}else {
    dir.create(file.path(figures.out.dir), showWarnings = TRUE)
}

## PATH FOR TABLES
tables.out.dir <- paste0(out.dir, "tables/")
if (file.exists(tables.out.dir)){
    system(paste0('rm -r ','"',tables.out.dir,'"',"*"))
}else {
    dir.create(file.path(tables.out.dir), showWarnings = TRUE)
}

## CREATE ENRICHR DIR
enrichr.out.dir <- paste0(out.dir, "enrichr/")
if (file.exists(enrichr.out.dir)){
    system(paste0('rm -r ','"',enrichr.out.dir,'"',"*"))
}else {
    dir.create(file.path(enrichr.out.dir), showWarnings = TRUE)
}

## CREATE OUTPUT DIRECTORY CALLED "PATHWAYS" INSIDE THE FOLDER OF OUTPUT
path.out.dir <- paste0(out.dir, "pathways/")
if (file.exists(path.out.dir)){
    system(paste0('rm -r ','"',path.out.dir,'"',"*"))
}else {
    dir.create(file.path(path.out.dir), showWarnings = TRUE)
}

## KEGG DOWNLOADED PATHWAYS
kegg.out.dir <- paste0(out.dir, "pathways/kegg_downloads/")
if (file.exists(kegg.out.dir)){
    system(paste0('rm -r ','"',kegg.out.dir,'"',"*"))
}else {
    dir.create(file.path(kegg.out.dir), showWarnings = TRUE)
}

## LOAD FILE WITH THE FUNCTIONS USED IN THIS SCRIPT
source('functions.R')

## ## FOR FEATURESEQ ANALYSIS - COMMENT ABOVE
## feature.counts.dir <- paste0(getwd(),"/results_feature_counts")
## if (file.exists(feature.counts.dir)){
##     in.file <- "nf-rnaSeqCount-scleroderma/featureCounts/gene_counts_final.txt"
##     out.dir <- paste0(feature.counts.dir,'/')
## }else {
##     dir.create(file.path(feature.counts.dir), showWarnings = FALSE)
##     in.file <- "nf-rnaSeqCount-scleroderma/featureCounts/gene_counts_final.txt"
##     out.dir <- paste0(feature.counts.dir,'/')
## }

## ------------------------------ LOAD COUNTS MATRIX DATA & ASSIGN NAMES, CONDITION, SITE AND GENDER ETC.
the.data <- read.csv(file=in.file, header=TRUE, row.names=1, sep = '\t')
the.data <- the.data[,order(names(the.data))]
the.data <- the.data[!row.names(the.data) %in% c("__no_feature", "__ambiguous", "__too_low_aQual","__not_aligned","__alignment_not_unique"),]

## ASSIGN SAMPLE NAMES
new_names <- c('CASE06_LE-NOXFS','CASE09_LE-NOXFS','CTRL03','CTRL04','CTRL05','CASE06','CASE07','CASE08','CASE09',
               'CTRL10','CASE11','CASE12','CTRL13','CTRL14','CASE15','CTRL16','CASE17','CTRL18',
               'CASE19','CASE20','CASE21','CASE22_RE_XFS','CASE22','CTRL24','CASE25','CASE07_RE_NOXFS','CASE27',
               'CASE28','CTRL29','CASE15_RE_XFS','CTRL31','CASE32','CASE33','CTRL34','CTRL35','CASE36',
               'CTRL37','CASE38','CASE39','CTRL40') 


colnames(the.data) <- new_names
the.data <- the.data[,order(names(the.data))]
the.data <- as.matrix(the.data)
compare.COMP <- list(c('CTRL03','CTRL04','CTRL05','CTRL10','CTRL13','CTRL14','CTRL16','CTRL18','CTRL24','CTRL29','CTRL31','CTRL34','CTRL35','CTRL37','CTRL40'),
                     c('CASE06_LE-NOXFS','CASE09_LE-NOXFS','CASE06','CASE07_RE_NOXFS','CASE07','CASE08','CASE09','CASE11','CASE12','CASE15','CASE17','CASE19','CASE20','CASE21','CASE22_RE_XFS','CASE22','CASE25','CASE27','CASE28','CASE15_RE_XFS','CASE32','CASE33','CASE36','CASE38','CASE39'))

## SAMPLE CONRITION
condition <- factor(c(rep('case',25),rep('ctrl',15))) ## MUST AUTOMATE - USER INPUT MAYBE
sex <- factor(c('M','M','F','M','F','F','M','F','M','M','M','F','F','M','F','M','M','M','F','M',
                'M','F','M','M','M','F','F','F','F','F','F','F','M','M','F','M','M','F','F','M'))
age <- c('53','63','45','53','76','67','80','80','69','62','58','84','84','71','87','85','80','73','74','66',
         '74','77','53','53','51','78','58','75','67','58','81','83','78','63','69','71','78','55','67','64')
the.coldata <- DataFrame(row.names=colnames(the.data), condition, sex, age)

## SUBSET
the.data <- the.data[, !colnames(the.data) %in% c('CASE06_LE-NOXFS','CASE09_LE-NOXFS','CASE15_RE_XFS','CASE22_RE_XFS','CASE07_RE_NOXFS')]
the.coldata <- the.coldata[!row.names(the.coldata) %in% c('CASE06_LE-NOXFS','CASE09_LE-NOXFS','CASE15_RE_XFS','CASE22_RE_XFS','CASE07_RE_NOXFS'),]

the.data <- the.data[which(rowSums(the.data) >= 1),]
the.data <- the.data[which(rowSums(the.data >= 5) >= 3),]

## LOAD ENSEMBL DATASETS FOR HUMAN
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## GET ENSEMBLE ANNOTATIONS
ensembl.annotations <- get_ensembl_annotations(row.names(the.data))
#ensembl.annotations <- ensembl.annotations[which(!ensembl.annotations$HGNC %in% c("C12orf74","LINC00595","CCL3L3")),]

## PATH FOR TABLES
tables.out.dir <- paste0(out.dir, "tables/")
if (file.exists(tables.out.dir)){
    system(paste0('rm -r ','"',tables.out.dir,'"',"*"))
}else {
    dir.create(file.path(tables.out.dir), showWarnings = TRUE)
}

## DO DIFF!
COMP <- DESeqDataSetFromMatrix(countData=the.data, colData=the.coldata, design = ~ sex + condition)
COMP$condition <- relevel(COMP$condition, ref="ctrl")
COMP <- estimateSizeFactors(COMP)
COMP <- COMP[ rowSums(counts(COMP, normalized = TRUE) >= 5) > 3,]
dds.COMP <- DESeq(COMP)

## EVALUATE THE BEST LFC THRESHOLD TO USE BY CREATING A TABLE FOR LFC AT 0.05 FDR FOR ALL COMPARISONS
check_lfc <- function(dds) {
    the_list <- vector()
    for (lfc in c(seq(0,2,0.1))) {
        res <- results(dds, alpha=0.05, lfcThreshold=lfc, altHypothesis="greaterAbs")
        summary(res)
        sig <- sum(res$padj < 0.05, na.rm=TRUE)
        the_list <- append(the_list, sig)
    }
    return(the_list)
}

lfc.table <- data.frame(row.names=(seq(0,2,0.1)))
lfc.table$COMP <- check_lfc(dds.COMP)

lfc.table

res.COMP <- results(dds.COMP, alpha=0.05, lfcThreshold=0.5, altHypothesis="greaterAbs")
summary(res.COMP)
sum(res.COMP$padj < 0.05, na.rm=TRUE)

sig.COMP <- get_sig(res.COMP, 0.05, 0.5)
sig_genes <- rownames(sig.COMP)
summary(sig.COMP)

rld.COMP <- rlog(COMP, blind=FALSE)

selected.COMP <- row.names(sig.COMP)
genes.COMP <- get_gene_names(selected.COMP)

write.csv(genes.COMP, paste0(tables.out.dir,'genes.csv'))

## PLOTS
pdf(paste0(figures.out.dir,'plot.pca.pdf'), width = 7, height = 6, onefile=FALSE)
plot_pca(rld.COMP)
dev.off()

pdf(paste0(figures.out.dir,'plot.volcano.pdf'), width = 7, height = 6, onefile=FALSE)
my_volcano(res.COMP,0.5,0.05)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.pdf'),onefile=FALSE)
plot_heat(rld.COMP,selected.COMP,genes.COMP,title.COMP)
dev.off()

## FUNCTION FOR CREATING THE TABLE AND ALL THE STATS
diff_expr_results <- function(sig.res,genes,dds){
    ## GET DATA FROM SIG
    sig <- as.data.frame(sig.res)[order(sig.res$padj),]

    ## GET DATA FROM DDS (NORMALIZED COUNTS)
    dds <- as.data.frame(counts(dds, normalized=TRUE))

    ## MERGE THE TWO
    expr_table <- merge(sig,dds,by="row.names",sort=FALSE)
    colnames(expr_table)[1] <- "GeneID"
    colnames(expr_table)[3] <- "log2FC"

    ## GET THE DESCRIPTION, CROMOSOME AND HGNC ID FROM THE ANNOTATION TABLE OF THE DATASET
    expr_table$Description <- genes$description[match(expr_table$GeneID, genes$ensembl_id)]
    expr_table$Chr <- genes$chromosome[match(expr_table$GeneID, genes$ensembl_id)]
    expr_table$GeneID <- ifelse(genes$HGNC[match(expr_table$GeneID, genes$ensembl_id)] == "",
                                expr_table$GeneID,
                                genes$HGNC[match(expr_table$GeneID, genes$ensembl_id)])
    ## ORDER THE TABLE
    len <- length(colnames(expr_table))
    expr_table <- expr_table[c(1, (len-1):len, 2:(len-2))]

    ## FORMAT THE NUMBERS
    expr_table$pvalue <- formatC(expr_table$pvalue, format = "e", digits = 2)
    expr_table$padj <- formatC(expr_table$padj, format = "e", digits = 2)
    expr_table[,c(4:7,10:len)] <- sapply(expr_table[,c(4:7,10:len)], function(number) format(round(number, 2), nsmall = 2))

    ## WRITE OUT TABLE FOR EXPRESSION AND ALL
    write.table(expr_table,file=paste0(tables.out.dir,"stats.csv"), row.names=FALSE,sep="|")
    return(expr_table)
}

table.COMP <- (diff_expr_results(sig.COMP,genes.COMP,dds.COMP))

##### ========== END DIFFERENTIAL EXPRESSION ANALYSIS ========== #####

analysis.genes <- unique(selected.COMP)
the_table <- as.data.frame(analysis.genes)
rownames(the_table) <- the_table$analysis.genes
the_table <- the_table[,-1]

the_names <- get_gene_names(analysis.genes)
the_table$HGNC <- the_names$HGNC[match(row.names(the_table), the_names$ensembl_id)]
the_table$NAME <- the_names$description[match(row.names(the_table), the_names$ensembl_id)]
write.csv(the_table, file=paste0(tables.out.dir,"galbladder.csv"))

final_genes <- the_table
final_genes.ensembl <- unique(row.names(final_genes))

final_genes.hgnc <- unique(the_names[the_names$ensembl_id %in% final_genes.ensembl,]$HGNC)
final_genes.hgnc <- final_genes.hgnc[final_genes.hgnc != ""]
write.csv(final_genes[,colnames(final_genes) %in% c("HGNC","NAME")], file=paste0(tables.out.dir,"final.gene.list.csv"), row.names=TRUE,col.names=TRUE)

## ENRICHR
dbs <- listEnrichrDbs()$libraryName
dbs <- dbs[dbs != "NCI-Nature_2015"]
my.db.list=c("GO_Molecular_Function_2013",
             "GO_Biological_Process_2013",
             "GO_Cellular_Component_2013",
             "GO_Molecular_Function_2015",
             "GO_Biological_Process_2015",
             "GO_Cellular_Component_2015",
             "GO_Molecular_Function_2017",
             "GO_Biological_Process_2017",
             "GO_Cellular_Component_2017",
             "GO_Molecular_Function_2017b",
             "GO_Biological_Process_2017b",
             "GO_Cellular_Component_2017b",
             "GO_Molecular_Function_2018",
             "GO_Biological_Process_2018",
             "GO_Cellular_Component_2018",
             "WikiPathways_2013",
             "WikiPathways_2015",
             "WikiPathways_2016",
             "KEGG_2013",
             "KEGG_2015",
             "KEGG_2016",
             "BioCarta_2013",
             "BioCarta_2015",
             "BioCarta_2016",
             "Panther_2015",
             "Panther_2016",
             "Reactome_2013",
             "Reactome_2015",
             "Reactome_2016",
             "OMIM_Disease",
             "OMIM_Expanded",
             "Jensen_DISEASES",
             "Jensen_COMPARTMENTS",
             "Jensen_TISSUES")

dbs <- dbs[which(dbs %in% my.db.list)]

## ALL SIGNIFICANT GENES COMBINED!
enriched.COMP <- enrichr(final_genes.hgnc, dbs)

put_enrichment_into_files <- function(set) {
    set.name <- deparse(substitute(set))
    set.name <- substr(set.name,10,nchar(set.name))
    lapply(names(set), function(db) 
        if (nrow(set[[db]]) >= 1 & nrow(set[[db]][which(set[[db]]$Adjusted.P.value <= 0.05),]) >= 1) {
            print(c("There are ",nrow(set[[db]])," in this set"))
            write.csv(as.data.frame(set[[db]][,c("Term","Overlap","P.value","Adjusted.P.value","Combined.Score","Genes")][which(set[[db]]$Adjusted.P.value <= 0.05),]),
                      file=paste0(enrichr.out.dir,"enrichr_annotation_",db,"_",set.name,".csv"),
                      row.names=FALSE,col.names=TRUE)
        } else {}
        )
}

put_enrichment_into_files(enriched.COMP)

## ATTEMPT AT VISUALIZING GO TERMS
plot.enrichr <- function(set,db.in){
    if (nrow(set[[db.in]]) >= 1 & nrow(set[[db.in]][which(set[[db.in]]$Adjusted.P.value <= 0.05),]) >= 1) {
        the.set <- set[[db.in]][order(set[[db.in]]$Adjusted.P.value, decreasing=FALSE),]
        if (nrow(the.set[which(the.set$Adjusted.P.value < 0.05 & str_detect(the.set$Term,"_Mus musculus") == FALSE),]) >= 1) {
            the.set <- head(the.set[which(the.set$Adjusted.P.value < 0.05 & str_detect(the.set$Term,"_Mus musculus") == FALSE),],20)
            the.set$gene.counts <- count.fields(textConnection(the.set$Genes),sep=";")
            the.set <- the.set[,which(colnames(the.set) %in% c("Term","P.value","Combined.Score","gene.counts"))]
            if (str_detect(db.in,"GO_") == TRUE) {
                ## the.set$Term <- gsub(".*\\(|\\)", "", the.set$Term)
                the.set$Term <- gsub(" \\(GO:.*","",the.set$Term)
            } else if (str_detect(db.in,"WikiPathways_") == TRUE) {
                the.set$Term <- gsub("_Homo.*","",the.set$Term)
            } else if (str_detect(db.in,"Reactome_") == TRUE) {
                the.set$Term <- gsub("_Homo.*","",the.set$Term)
            } else if (str_detect(db.in,"KEGG_")) {
                the.set$Term <- substr(the.set$Term,0,nchar(the.set$Term)-22)
            } else {}
        } else {
            the.set <- data.frame(Term="No Significant Results",P.value=0.1,Combined.Score=1,gene.counts=1)
        }
    } else {
        the.set <- data.frame(Term="NO SIGNIFICANT RESULTS",P.value=0.1,Combined.Score=1,gene.counts=1)
    }
    ## PLOT
    ggplot(data=the.set, aes(x=-log10(P.value), y=gene.counts, color=-log10(P.value), size=Combined.Score)) +
        geom_point(alpha=0.5,shape=1) +
        guides(size = "none", colour = "legend") +
        geom_text_repel(data=head(the.set,5), aes(label=Term),
                        size=1.5, min.segment.length = 0.4,
                        segment.size = 0.1, box.padding = 0.5,
                        ## label.padding = 0.1, 
                        ## nudge_y = 0.1,
                        ## nudge_x = 0.1,
                        direction = "both",
                        hjust = 1
                        ) +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Gene Counts") +
        scale_x_continuous(breaks = seq(0,15,by=0.5)) +
        scale_y_continuous(breaks = seq(0,15,1)) +
        scale_fill_manual(values="white") +
        scale_color_gradient(low='black',high='red') +
        theme(panel.grid.major = element_blank(),
              legend.position="none",
              ## legend.position = "bottom",
              ## legend.justification = "top",
              ## legend.key = element_rect(fill = "white", size=0.1),
              ## legend.text = element_text(size = 3),
              ## legend.title = element_text(size = 4,face = 'bold'),
              ## legend.key.size = unit(0.25, "cm"),
              ## legend.margin = margin(t = 0, unit='cm'),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=5),
              axis.text = element_text(size=4),
              plot.title = element_text(size=7, hjust = 0.5),
              plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"))
}


## enriched.COMP
go.biol.COMP <- plot.enrichr(enriched.COMP,"GO_Biological_Process_2018")
go.mol.COMP <- plot.enrichr(enriched.COMP,"GO_Molecular_Function_2018")
go.cell.COMP <- plot.enrichr(enriched.COMP,"GO_Cellular_Component_2018")
WIKI.COMP <- plot.enrichr(enriched.COMP,"WikiPathways_2016")
KEGG.COMP <- plot.enrichr(enriched.COMP,"KEGG_2016")
REACTOME.COMP <- plot.enrichr(enriched.COMP,"Reactome_2016")
pdf(paste0(figures.out.dir,'enrichment.plot.COMP_Terms.pdf'), paper='a4', width = 0, height = 0, onefile=FALSE)
grid.arrange(go.biol.COMP + ggtitle("GO Biological Process 2018"),
             go.mol.COMP + ggtitle("GO Molecular Function 2018"),
             go.cell.COMP + ggtitle("GO Cellular Component 2018"),
             WIKI.COMP + ggtitle("WikiPathways 2016"),
             KEGG.COMP + ggtitle("KEGG 2016"),
             REACTOME.COMP + ggtitle("Reactome 2016"),
             ## + ggtitle(title.MILD.BACK),
             ## + ggtitle(title.SEVERE.MILD),
             layout_matrix = rbind(c(1,2), c(3,4), c(5,6)),
             widths = c(1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for WIKI Pathways of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## ------------------------------ Pathway Analysis ------------------------------#####
## Map ENSEMBL IDs to Entrez for pathway analysis
## Order by P-value

## Detach the PROPRER PACKAGE - INTERFERES WITH PATHWAY ANALYSIS
## detach("package:PROPER", unload=TRUE)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
columns(org.Hs.eg.db)

mapping <- data.frame(row.names=final_genes.hgnc)
mapping$ensembl <- mapIds(org.Hs.eg.db, keys=row.names(mapping), column="SYMBOL", keytype="SYMBOL", multiVals="first")
mapping$entrez <- mapIds(org.Hs.eg.db, keys=row.names(mapping), column="ENTREZID", keytype="SYMBOL", multiVals="first")
mapping$name   <- mapIds(org.Hs.eg.db, keys=row.names(mapping), column="GENENAME", keytype="SYMBOL", multiVals="first")

do_pathway <- function(res) {
    ## GET THE NAME OF THE COMPARISON SET WE ARE WORKING WITH FROM RES NAME
    path.name <- deparse(substitute(res))
    path.name <- substring(path.name,5,nchar(path.name))
    
    ## ORDER THE DIFFERENTIAL EXPRESSION RESULTS BY ADJUSTED PVALUE
    resPath <- res[order(res$padj),]
    
    ## GET THE ENTREZ ID'S FOR MAPPING TO KEGG PATHWAYS
    resPath$symbol = mapIds(org.Hs.eg.db, keys=row.names(resPath), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    resPath$entrez = mapIds(org.Hs.eg.db, keys=row.names(resPath), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
    resPath$name =   mapIds(org.Hs.eg.db, keys=row.names(resPath), column="GENENAME", keytype="ENSEMBL", multiVals="first")
    
    ## EXTRACT FOLD-CHANGE INFORMATION FROM THE DE SET
    foldchanges = resPath$log2FoldChange
    names(foldchanges) = resPath$entrez
    
    ## RUN PATHWAY ANALYSIS WITH GAGE
    keggres <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE) #, set.size=c(10,50))
    
    ## GET THE PATHWAYS UP-REGULARED
    select.greater <- keggres$greater[, "q.val"] < 0.05 & !is.na(keggres$greater[,"q.val"])
    up.reg <- rownames(keggres$greater)[select.greater]
    lapply(up.reg, function(x) print(paste0('upreg_',x)))
    
    ## GET THE PATHWAYS DOWN-REGULARED
    select.less <- keggres$less[, "q.val"] < 0.05 & !is.na(keggres$less[,"q.val"])
    down.reg <- rownames(keggres$less)[select.less]
    lapply(down.reg, function(x) print(paste0('downreg_',x)))
    
    ## GET THE TABLE FOR THE VALUES
    if (nrow(keggres$greater[select.greater,]) >=1) {
        ## print(head(keggres$greater[select.greater,]))
        path.table <- as.data.frame(keggres$greater[select.greater,])
        path.table$p.geomean <- formatC(path.table$p.geomean, format = "e", digits = 2)
        path.table$stat.mean <- formatC(path.table$stat.mean, format = "e", digits = 2)
        path.table$p.val <- formatC(path.table$p.val, format = "e", digits = 2)
        path.table$q.val <- formatC(path.table$q.val, format = "e", digits = 2)
        path.table$exp1 <- formatC(path.table$exp1, format = "e", digits = 2)
        print(path.table)
        write.csv(path.table,
                  file=paste0(path.out.dir,"pathway_table_",path.name,"_greater.csv"),
                  row.names=TRUE,col.names=TRUE)

        out.table = paste0(path.out.dir,"pathway_table_",path.name,"_greater.tex")
        table.head <- "\\begin{longtable}{ L{28em} L{5em} R{5em} R{5em} R{5em} R{5em} R{5em} }"
        table.caption <- paste0("\\caption[Up-regulated pathways in the ",path.name," comparison identified by \\gage{}]{\\textbf{Up-regulated pathways in the ",path.name," comparison identified by \\gage{}.}}")
        table.label <- paste0("\\label{tab:gage.up.",tolower(path.name),"}\\\\")
        ##
        write.table(table.head, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
        write.table(table.caption, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        write.table(table.label, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        ##
        write.table("\\toprule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(paste0("\\rowcolor{gray!25} Pathway name & ", paste(names(path.table),collapse=" & ")),
                    file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\midrule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(path.table, file=out.table, row.names=TRUE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\bottomrule\n\\end{longtable}", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")

    }else {}
    
    if (nrow(keggres$less[select.less,]) >=1) {
        ## print(head(keggres$less[select.less,]))
        path.table <- as.data.frame(keggres$less[select.less,])
        path.table$p.geomean <- formatC(path.table$p.geomean, format = "e", digits = 2)
        path.table$stat.mean <- formatC(path.table$stat.mean, format = "e", digits = 2)
        path.table$p.val <- formatC(path.table$p.val, format = "e", digits = 2)
        path.table$q.val <- formatC(path.table$q.val, format = "e", digits = 2)
        path.table$exp1 <- formatC(path.table$exp1, format = "e", digits = 2)
        print(path.table)
        write.csv(path.table,
                  file=paste0(path.out.dir,"pathway_table_",path.name,"_less.csv"),
                  row.names=TRUE,col.names=TRUE)
        
        out.table = paste0(path.out.dir,"pathway_table_",path.name,"_less.tex")
        table.head <- "\\begin{longtable}{ L{28em} L{5em} R{5em} R{5em} R{5em} R{5em} R{5em} }"
        table.caption <- paste0("\\caption[Down-regulated pathways in the ",path.name," comparison identified by \\gage{}]{\\textbf{Down-regulated pathways in the ",path.name," comparison identified by \\gage{}.}}")
        table.label <- paste0("\\label{tab:gage.down.",tolower(path.name),"}\\\\")
        ##
        write.table(table.head, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
        write.table(table.caption, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        write.table(table.label, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        ##
        write.table("\\toprule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(paste0("\\rowcolor{gray!25} Pathway name & ", paste(names(path.table),collapse=" & ")),
                    file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\midrule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(path.table, file=out.table, row.names=TRUE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\bottomrule\n\\end{longtable}", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    }else {}
    
    ## GET THE LIST
    keggresnames <- c(down.reg, up.reg)
    print(keggresnames)
    
    keggresids <- substr(c(up.reg, down.reg), 1, 8)
    keggresids <- substr(keggresnames, 1, 8)
    
    keggresnames.table <- data.frame(row.names=keggresids)
    keggresnames.table$names <- substr(keggresnames,10,nchar(keggresnames))
    
    print(keggresids)
    print(keggresnames.table)
    
    ## PLOT PATHWAYS USING PATHVIEW
    sapply(keggresids, function(pid) pathview(gene.data=foldchanges,
                                              pathway.id=pid,
                                              species="hsa" ,
                                              kegg.dir = kegg.out.dir,
                                              limit = list(gene = 10, cpd = 10),
                                              out.suffix = paste0(gsub(' ','.',keggresnames.table$names[match(pid, row.names(keggresnames.table))]),"_",path.name)
                                              )
           )

    ## MOVE PATHWAYS TO PATHWAY FOLDER
    system(paste0('mv *',path.name,".png ",'"',path.out.dir,'"'))
}

try(do_pathway(res.COMP), TRUE)
## try(do_pathway(res.CASE.ARM), TRUE)
## try(do_pathway(res.CASE.BACK), TRUE)
## try(do_pathway(res.CASE), TRUE)
## try(do_pathway(res.SEVERE.ARM), TRUE)
## try(do_pathway(res.SEVERE.BACK), TRUE)
## try(do_pathway(res.MILD.ARM), TRUE)
## try(do_pathway(res.MILD.BACK), TRUE)
## try(do_pathway(res.SEVERE.MILD), TRUE)

library(PROPER)
