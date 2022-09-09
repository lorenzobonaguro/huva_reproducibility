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

steelblue <- c("protein localization to membrane", "post-translational protein modification", "carboxylic acid biosynthetic process")
darkgreen <- c("neutrophil activation", "immune response-activating signal transduction", "blood coagulation")
darkgrey <- c("antigen processing and presentation of peptide antigen", "regulation of cell cycle G2/M phase transition", "T cell receptor signaling pathway")
darkorange <- c("sensory system development", "DNA geometric change", "V(D)J recombination")
orchid <- c("lymphocyte differentiation", "T cell differentiation", "B cell proliferation")
maroon <- c("activation of immune response", "positive regulation of defense response", "activation of innate immune response")
gold <- c("regulation of cellular protein catabolic process", "regulation of cellular protein catabolic process", "tRNA metabolic process")

GO_plot <- c(steelblue,darkgreen,darkgrey,darkorange,orchid,maroon,gold)

plot <- GO[GO$Description %in% GO_plot,]

plot <- plot[plot$p.adjust<0.05,]

plot$module <- factor(plot$module, levels = rev(c("maroon", "darkgreen", "darkgrey", "steelblue", "gold", "darkorange", "orchid")))

plot$Description <- factor(plot$Description)

gg <- ggplot(plot, aes(x=Description, y=module, size=Count, fill=-log10(p.adjust))) +
        geom_point(shape=21, color="black") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
        scale_fill_gradient(limits = c(-log10(0.05),16), low = "yellow", high = "red", space = "Lab") +
        scale_size_continuous(limits = c(0, 60), range = c(2, 8)) +
        ggtitle("Fig. S10g - GOEA - selected terms")

print(gg)
