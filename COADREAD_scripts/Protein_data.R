RPPA


RPPA <- read.delim("~/Downloads/gdac.broadinstitute.org_COADREAD.RPPA_AnnotateWithGene.Level_3.2016012800.0.0/COADREAD.rppa.txt")
converter <- read.delim("~/Downloads/gdac.broadinstitute.org_COADREAD.RPPA_AnnotateWithGene.Level_3.2016012800.0.0/COADREAD.antibody_annotation.txt")
RPPA1 <- separate(RPPA, "Composite.Element.REF", c("Gene.Name", "Composite.Element.REF"), sep = "\\|")

head(RPPA1)[, 1:10]

unique(RPPA1$Gene.Name)
