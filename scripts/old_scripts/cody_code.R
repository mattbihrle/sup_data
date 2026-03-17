####Running an NMDS on the data
GLNPO_4_NMDS <- fread("~/Google Drive/My Drive/Projects/Collaborations/USEPA_GLNPO/Combined_Sylph_secondround.tsv") %>%
  separate_wider_delim(clade_name, names = c("Domain", "Phylum", "Class", "Order", "Family","Genus", "Species","GenomeHit"), delim ="|" , too_few = "align_start") %>%
  replace_na(list(Phylum = "unknown", Class = "unknown", Order= "unknown", Family= "unknown",Genus= "unknown", Species= "unknown", GenomeHit= "unknown")) %>%
  filter(GenomeHit != "unknown") %>%
  select(-Domain,-Phylum, -Class, -Order, -Family, -Genus, -Species, -GenomeHit)

NMDS <- metaMDS(decostand(t(GLNPO_4_NMDS), "hellinger"), distance="bray", autotransform=FALSE, binary=FALSE, noshare=TRUE, zerodist="ignore")

NMDS_scores_MAG<-scores(NMDS, display = "sites") %>%
  data.frame() %>%
  rownames_to_column("Sample") %>%
  left_join(TotalSizes, by="Sample") %>%
  select(-V3) %>%
  na.omit()

adonis2(decostand(t(GLNPO_4_NMDS), "hellinger") ~NMDS_scores_MAG$Lake, strata=NMDS_scores_MAG$Filter_size)

#### All samples
NMDS_scores_MAG %>%
  #filter( Lake == c("Erie", "Ontario")) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, fill=Lake, shape = Date)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_shape_manual(values=c(21,22,23,24,25), name = "Season") +
  #scale_fill_manual(values= c("white", "#3182BD"), name="", labels= c("Deep Chlorophyll Max", "Integrate"))+
  guides(fill = guide_legend("Legend fill", override.aes = list(shape = 21)))