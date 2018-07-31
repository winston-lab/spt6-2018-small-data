library(tidyverse)
library(forcats)
library(ggthemes)

df = read_tsv('tfiib-qpcr.tsv') %>%
    mutate(condition=fct_inorder(condition, ordered=TRUE))

plot_qpcr = function(df, gene_id, norm="input"){
    subdf = df %>% filter(gene==gene_id)
    if (norm != "input"){
        subdf = df %>% filter(gene==norm) %>%
            select(condition, replicate, spikein=value) %>%
            left_join(subdf, ., by=c("condition", "replicate")) %>%
            mutate(value=value/spikein)
    }

    strand = subdf %>% distinct(strand) %>% pull(strand)

    if (strand=="+"){
        subdf = subdf %>%
            mutate_at(vars(amplicon_start, amplicon_end, transcript_end, orf_start, orf_end, transcript_start),
                      funs(.-transcript_start))
    } else if (strand=="-"){
        subdf = subdf %>% mutate_at(vars(amplicon_start, amplicon_end, transcript_start, orf_start, orf_end, transcript_end),
                                  funs(transcript_end-.))
    }

    max_val = subdf %>% pull(value) %>% max()
    txn_end = max(c(subdf[["transcript_start"]], subdf[["transcript_end"]]))
    orf_start = min(c(subdf[["orf_start"]], subdf[["orf_end"]]))
    orf_end = max(c(subdf[["orf_start"]], subdf[["orf_end"]]))

    plot = ggplot(data = subdf, aes(x=(amplicon_start+amplicon_end)/2, y=value,
                            group=fct_inorder(paste(amplicon_start, condition), ordered=TRUE), fill=condition)) +
        geom_vline(aes(xintercept = amplicon_start), linetype="dashed", alpha=0.25) +
        geom_vline(aes(xintercept = amplicon_end), linetype="dashed", alpha=0.25) +
        annotate(geom="segment", color="black", size=2,
                 x=0, xend=txn_end, y=-.04*max_val, yend=-.04*max_val) +
        annotate(geom="rect", color="black", fill="grey95", size=0.5,
                 xmin=orf_start, xmax=orf_end, ymin=-.08*max_val, ymax=0) +
        annotate(geom="text", label=gene_id, x=(orf_start+orf_end)/2,
                 y = -.035*max_val, size=4, fontface="bold") +
        annotate(geom="segment", y=-0.06*max_val, yend=-0.06*max_val,
                 x=0.2*(orf_end-orf_start) + orf_start,
                 xend=0.8*(orf_end-orf_start) + orf_start,
                 size=0.5, arrow=arrow(length=unit(0.03, "npc"))) +
        geom_boxplot(position=position_dodge(450), width=500) +
        geom_point(shape=21, size=1, stroke=1, color="black",
                   position=position_jitterdodge(jitter.width=100, dodge.width = 450)) +
        scale_fill_ptol() +
        scale_color_ptol() +
        scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                           labels= function(x){if_else(x==0, "TSS",
                                                       if(abs(txn_end)>500){as.character(x/1000)}
                                                       else {as.character(x)})},
                           name=paste("distance from TSS", if(abs(txn_end)>500){"(kb)"}
                                      else {"(nt)"})) +
        ylab("enrichment (AU)") +
        ggtitle(paste("TFIIB ChIP-qPCR at", gene_id),
                subtitle = if(norm != "input"){paste("normalized to S. pombe", norm, "promoter")} else {"input normalized"}) +
        theme_light() +
        theme(text = element_text(size=12, face="bold"),
              legend.title = element_blank(),
              legend.text = element_text(size=12),
              axis.text = element_text(size=12, color="black"),
              axis.title.x = element_text(size=10, face="plain"),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=10, face="plain"))
}

pma = plot_qpcr(df, gene_id="pma1+", norm="input")
act = plot_qpcr(df, gene_id="act1+", norm="input")
vam_input = plot_qpcr(df, gene_id="VAM6", norm="input")
vam_pma = plot_qpcr(df, gene_id="VAM6", norm="pma1+")
vam_act = plot_qpcr(df, gene_id="VAM6", norm="act1+")
hsp_input = plot_qpcr(df, gene_id="HSP82", norm="input")
hsp_pma = plot_qpcr(df, gene_id="HSP82", norm="pma1+")
hsp_act = plot_qpcr(df, gene_id="HSP82", norm="act1+")
pma_input = plot_qpcr(df, gene_id="PMA1", norm="input")
pma_pma = plot_qpcr(df, gene_id="PMA1", norm="pma1+")
pma_act = plot_qpcr(df, gene_id="PMA1", norm="act1+")
flo_input = plot_qpcr(df, gene_id="FLO8", norm="input")
flo_pma = plot_qpcr(df, gene_id="FLO8", norm="pma1+")
flo_act = plot_qpcr(df, gene_id="FLO8", norm="act1+")
avt_input = plot_qpcr(df, gene_id="AVT2", norm="input")
avt_pma = plot_qpcr(df, gene_id="AVT2", norm="pma1+")
avt_act = plot_qpcr(df, gene_id="AVT2", norm="act1+")
ypt_input = plot_qpcr(df, gene_id="YPT52", norm="input")
ypt_pma = plot_qpcr(df, gene_id="YPT52", norm="pma1+")
ypt_act = plot_qpcr(df, gene_id="YPT52", norm="act1+")

ggsave('tfiib-chip-qpcr_spom-pma1-input-norm.svg', plot=pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_spom-act1-input-norm.svg', plot=act, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-vam6-input-norm.svg', plot=vam_input, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-vam6-spom-pma1-norm.svg', plot=vam_pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-vam6-spom-act1-norm.svg', plot=vam_act, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-hsp82-input-norm.svg', plot=hsp_input, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-hsp82-spom-pma1-norm.svg', plot=hsp_pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-hsp82-spom-act1-norm.svg', plot=hsp_act, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-pma1-input-norm.svg', plot=pma_input, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-pma1-spom-pma1-norm.svg', plot=pma_pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-pma1-spom-act1-norm.svg', plot=pma_act, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-flo8-input-norm.svg', plot=flo_input, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-flo8-spom-pma1-norm.svg', plot=flo_pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-flo8-spom-act1-norm.svg', plot=flo_act, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-avt2-input-norm.svg', plot=avt_input, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-avt2-spom-pma1-norm.svg', plot=avt_pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-avt2-spom-act1-norm.svg', plot=avt_act, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-ypt52-input-norm.svg', plot=ypt_input, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-ypt52-spom-pma1-norm.svg', plot=ypt_pma, width=16, height=10, units="cm")
ggsave('tfiib-chip-qpcr_scer-ypt52-spom-act1-norm.svg', plot=ypt_act, width=16, height=10, units="cm")
