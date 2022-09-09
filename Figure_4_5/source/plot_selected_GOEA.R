for (i in names(clust_prof)) {
  
  for (j in names(clust_prof[[i]])[grep("_table", names(clust_prof[[i]]))]) {
    
    clust_prof[[i]][[j]]$module <- i
    
  } 
  
}

collector <- list()

for (i in names(clust_prof)) {
  
  for (j in names(clust_prof[[i]])[grep("GO_table", names(clust_prof[[i]]))]) {
    
    collector[[i]] <- clust_prof[[i]][[j]]
    
  }
  
}

GO <-  do.call("rbind", collector)

steelblue <- c("protein deubiquitination","regulation of cytokine−mediated signaling pathway","interferon−gamma production")
sendybrown <- c("negative regulation of growth","retinoid metabolic process","viral budding")
indianred <- c("ribonucleoprotein complex biogenesis","ncRNA metabolic process","RNA splicing")
darkgreen <- c("granulocyte activation","regulation of inflammatory response","tumor necrosis factor production")
pink <- c("autophagy","Fc−gamma receptor signaling pathway","chromatin silencing")
lightblue <- c("T cell activation","lymphocyte differentiation","lymphocyte costimulation")
darkgrey <- c("B cell activation","leukocyte proliferation","B cell receptor signaling pathway")
darkorange <- c("negative regulation of immune response","respiratory burst","macrophage migration")
orchid <- c("immune response−regulating signaling pathway","phagocytosis","Notch signaling pathway")
maroon <- c("neutrophil mediated immunity","phospholipid metabolic process","regulation of viral life cycle")
khaki <- c("defense response to other organism","response to interferon−gamma","type I interferon production")
turquoise <- c("chromosome segregation","adaptive immune response","lymphocyte apoptotic process")
gold <- c("cell−cell recognition","leukocyte differentiation","regulation of regulatory T cell differentiation")
lightgreen <- c("cell killing","natural killer cell mediated immunity","positive regulation of interferon−gamma production")

GO_plot <- c(steelblue,sendybrown,indianred,darkgreen,pink,lightblue,darkgrey,darkorange,orchid,maroon,khaki,turquoise,gold,lightgreen)

plot <- GO[GO$Description %in% GO_plot,]

plot <- plot[plot$p.adjust<0.05,]

plot$module <- factor(plot$module, levels = rev(c("steelblue", "darkorange", "lightgreen", "orchid", "pink", "darkgreen", "maroon", "khaki", "sandybrown", "darkgrey", "indianred", "lightblue", "gold", "turquoise")))

plot$Description <- factor(plot$Description, levels = GO_plot)

gg <- ggplot(plot, aes(x=Description, y=module, size=Count, fill=-log10(p.adjust))) +
        geom_point(shape=21, color="black") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        scale_fill_gradient(limits = c(-log10(0.05),13), low = "yellow", high = "red", space = "Lab") +
        scale_size_continuous(limits = c(0, 40), range = c(2, 8)) +
        ggtitle("Fig. 4g - GOEA - selected terms")

print(gg)
