library(tidyverse)
library(forcats)
library(ggthemes)

df = read_csv("TFIIB_DST1.csv") %>%
    select(c(1,2,15,16,17,18,19)) %>%
    magrittr::set_colnames(c("sample", "factor_id", "gray_min", "gray_max", "gray_mean", "gray_median", "auc")) %>%
    select(sample, factor_id, auc) %>%
    spread(key=factor_id, value=auc) %>%
    magrittr::set_colnames(c("sample", "background", "spikein", "tfiib")) %>%
    separate(col=sample, into=c("group", "replicate"), sep="-", remove=FALSE, convert=TRUE) %>%
    mutate(group=fct_recode(group, WT="wt", `spt6-1004`="spt")) %>%
    mutate_at(vars(spikein, tfiib), funs(.-background)) %>%
    mutate(spikenorm=tfiib/spikein,
           group=fct_inorder(group, ordered=TRUE) %>% fct_rev()) %>%
    group_by(group) %>%
    mutate(group_mean = mean(spikenorm))


#rescale data so that mean of WT group is 1
wt_og_mean = df %>% filter(group=="WT") %>% distinct(group_mean) %>% pull(group_mean)

df = df %>%
    mutate(spikenorm_scaled=spikenorm+1-wt_og_mean)

summary_df = df %>%
    summarise(group_mean_scaled = mean(spikenorm_scaled),
              group_sem = sd(spikenorm_scaled)/sqrt(n()))

max_y = summary_df %>%
    mutate(max_y = (group_mean_scaled+1.96*group_sem)*1.05) %>%
    pull(max_y) %>% max()

barplot = ggplot() +
    geom_col(data = summary_df, aes(x=group, y=group_mean_scaled, fill=group)) +
    geom_errorbar(data = summary_df, width=0.2,
                  aes(x=group, ymin=group_mean_scaled-1.96*group_sem,
                      ymax=group_mean_scaled+1.96*group_sem)) +
    geom_jitter(data = df, aes(x=group, y=spikenorm_scaled),
                width=0.2, size=1) +
    geom_text(data = summary_df, parse=TRUE,
              aes(x=group, y=(group_mean_scaled-1.96*group_sem)-0.05,
                  label = paste(round(group_mean_scaled,2), "%+-%", round(1.96*group_sem, 2))),
              size=7/72*25.4) +
    scale_fill_ptol(guide=FALSE) +
    scale_y_continuous(limits = c(0, max_y),
                       expand=c(0,0),
                       name="relative signal (a.u.)") +
    ggtitle("TFIIB levels by Western blot",
            subtitle = "normalized to Dst1 spike-in") +
    theme_light() +
    theme(text = element_text(size=9, color="black", face="plain"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=9, color="black"),
          axis.text.y = element_text(size=7, color="black"),
          plot.title = element_text(size=9),
          plot.subtitle = element_text(size=9))

ggsave("tfiib_western.svg", plot=barplot, width=6, height=11, units="cm")
