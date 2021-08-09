source("function-library.R")

library(tidyverse)
library(cowplot)

H <- 1.2
W <- 1.5

class_conversion <- c(
  normal = "Normal",
  mds = "MDS",
  anemia = "Anaemia"
)

fine_class_conversion <- c(
  `normal finding` = "Normal",
  `MDS SF3B1 mutated` = "SF3B1-mutant",
  `MDS RUNX1 mutated` = "RUNX1-mutant",
  `MDS SRSF2 mutated` = "SRSF2-mutant",
  `MDS U2AF1 mutated` = "U2AF1-mutant",
  `Iron deficiency anaemia` = "Iron deficiency",
  `Megaloblastic anemia` = "Megaloblastic"
)

fine_colours <- c(
  "grey50",
  "#3B4992FF","#008280FF","#008B45FF","#631879FF",
  "#fc5f64","#BB0021FF")
names(fine_colours) <- fine_class_conversion

all_conditions <- read_csv("data/all_classes.csv",progress = F) %>%
  mutate(Condition = factor(class_conversion[coarse_class],levels = class_conversion)) %>%
  mutate(fine_class = factor(fine_class_conversion[fine_class],levels = fine_class_conversion)) %>% 
  merge(read_csv("data/blood_count_data.csv"),by = "id",progress = F,all.x = T) %>%
  merge(read_csv("data/demographic_data.csv"),by = "id",progress = F,all.x = T) %>%
  group_by(Condition) %>%
  mutate(N_coarse = length(unique(id))) %>%
  group_by(fine_class) %>%
  mutate(N_fine = length(unique(id))) %>%
  ungroup %>%
  subset(sex != "n.k") %>%
  subset(!grepl("[a-zA-Z]",wbc_ul))  %>%
  subset(!grepl("[a-zA-Z]",hb_g_dl)) %>%
  subset(!grepl("[a-zA-Z]",plt_ul))

n_condition_coarse <- all_conditions %>%
  arrange(Condition) %>% 
  select(Condition,N = N_coarse) %>%
  distinct

n_condition_fine <- all_conditions %>%
  arrange(fine_class) %>%
  select(fine_class,N = N_fine) %>%
  distinct

age_plot <- all_conditions %>%
  ggplot(aes(x = Condition,y = age,fill = fine_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL,
                    labels = sprintf("%s (n=%s)",n_condition_fine$fine_class,n_condition_fine$N)) + 
  ylab("Age") + 
  scale_x_discrete(labels = sprintf("%s\n(n=%s)",n_condition_coarse$Condition,n_condition_coarse$N)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.2,"cm")) + 
  ggsave("figures/cohort-age-condition.pdf",height = H,width = W,useDingbats = F)

wbc_plot <- all_conditions %>%
  ggplot(aes(x = Condition,y = as.numeric(wbc_ul),fill = fine_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL,
                    labels = sprintf("%s (n=%s)",n_condition_fine$fine_class,n_condition_fine$N)) + 
    ylab("WBC count (/uL)") + 
  scale_x_discrete(labels = sprintf("%s\n(n=%s)",n_condition_coarse$Condition,n_condition_coarse$N)) +
  scale_y_continuous(trans = 'log10') +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.2,"cm")) + 
  ggsave("figures/cohort-wbc-condition.pdf",height = H,width = W,useDingbats = F)

hb_plot <- all_conditions %>%
  ggplot(aes(x = Condition,y = as.numeric(hb_g_dl),fill = fine_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL,
                    labels = sprintf("%s (n=%s)",n_condition_fine$fine_class,n_condition_fine$N)) + 
  ylab("Haemoglobin concentration (g/dL)") + 
  scale_x_discrete(labels = sprintf("%s\n(n=%s)",n_condition_coarse$Condition,n_condition_coarse$N)) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.2,"cm")) + 
  ggsave("figures/cohort-hg-condition.pdf",height = H,width = W,useDingbats = F)

plt_plot <- all_conditions %>%
  ggplot(aes(x = Condition,y = as.numeric(plt_ul),fill = fine_class)) + 
  geom_boxplot(size = 0.25,outlier.size = 0.5) + 
  theme_pretty(base_size = 6) + 
  scale_fill_manual(values = fine_colours,name = NULL,
                    labels = sprintf("%s (n=%s)",n_condition_fine$fine_class,n_condition_fine$N)) + 
  ylab("Platelet count (/uL)") + 
  scale_x_discrete(labels = sprintf("%s\n(n=%s)",n_condition_coarse$Condition,n_condition_coarse$N)) +
  scale_y_continuous(trans = 'log10',labels = scientific) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.2,"cm")) + 
  ggsave("figures/cohort-plt-condition.pdf",height = H,width = W,useDingbats = F)

plot_grid(
  age_plot,wbc_plot,hb_plot,plt_plot,
  align = "hv"
) + 
  ggsave(filename = "figures/cohort-age-wbc-hb-plt-condition.pdf",height = H*2,width = W*2,useDingbats = F)

(all_conditions %>%
    ggplot(aes(x = Condition,y = age,fill = fine_class)) + 
    geom_boxplot(size = 0.25,outlier.size = 0.5) + 
    theme_pretty(base_size = 6) + 
    scale_fill_manual(values = fine_colours,name = NULL,
                      labels = sprintf("%s (n=%s)",n_condition_fine$fine_class,n_condition_fine$N)) + 
    ylab("Age") + 
    scale_x_discrete(labels = sprintf("%s\n(n=%s)",n_condition_coarse$Condition,n_condition_coarse$N)) +
    theme(axis.title.x = element_blank(),
          legend.key.height = unit(0.2,"cm"),
          legend.key.width = unit(0.2,"cm"))) %>%
  get_legend() %>%
  ggsave(filename = "figures/fine-legend.pdf",height = 0.7,width = 1.1,useDingbats = F)

all_conditions %>%
  ggplot(aes(x = Condition,fill = sex)) + 
  geom_bar(position = "dodge") + 
  theme_pretty(base_size = 6) + 
  scale_fill_aaas(labels = c("Female","Male"),name = NULL) +
  ylab("Number of individuals") + 
  scale_x_discrete(labels = sprintf("%s\n(n=%s)",n_condition_coarse$Condition,n_condition_coarse$N)) +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.2,"cm"),
        legend.key.width = unit(0.2,"cm")) + 
  scale_y_continuous(expand = c(0,0)) +
  ggsave("figures/cohort-sex-condition.pdf",height = H*1.2,width = W*0.9,useDingbats = F)

glm(Condition == "MDS" ~ age + sex, data = all_conditions,
    family = stats::binomial) %>% 
  summary

glm(as.numeric(wbc_ul) ~ Condition, data = all_conditions) %>% 
  summary

t.test(as.numeric(wbc_ul) ~ fine_class,data = all_conditions[grepl("Megaloblastic|Normal",all_conditions$fine_class),])

glm(as.numeric(hb_g_dl) ~ Condition, data = all_conditions) %>% 
  summary

glm(as.numeric(plt_ul) ~ Condition, data = all_conditions) %>% 
  summary

t.test(as.numeric(plt_ul) ~ fine_class,
       data = all_conditions[grepl("Megaloblastic|Normal",all_conditions$fine_class),])

t.test(as.numeric(wbc_ul) ~ fine_class == "SF3B1-mutant",
       data = all_conditions[grepl("MDS",all_conditions$Condition),])

t.test(as.numeric(hb_g_dl) ~ fine_class == "SF3B1-mutant",
       data = all_conditions[grepl("MDS",all_conditions$Condition),])

t.test(as.numeric(plt_ul) ~ fine_class == "SF3B1-mutant",
       data = all_conditions[grepl("MDS",all_conditions$Condition),])

t.test(as.numeric(wbc_ul) ~ fine_class == "SF3B1-mutant",
       data = all_conditions[grepl("SF3B1|Normal",all_conditions$fine_class),])

t.test(as.numeric(plt_ul) ~ fine_class == "SF3B1-mutant",
       data = all_conditions[grepl("SF3B1|Normal",all_conditions$fine_class),])
