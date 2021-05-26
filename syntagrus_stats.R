library('lme4')   # glm
library(optimx)
library(dplyr)    # tibbles, piping


library('partykit') # partition trees
library("rcompanion")  # Cramer's V, better chi-squared

library("treemap")
library(ggplot2)
library(Cairo)    # for saving plots with non-ascii chars
library(scales)   # used here for scale fill on graph

LANG = 'ru'
w <- 400    # width and height for plotting
h <- 400

stgr_full = read.csv('./syntagrus/result/syntagrus_5.csv', encoding='UTF-8')


pronominality_en <- c('both pronominal', 'pronominal DO',
                      'pronominal IO', 'both non-pronominal')
pronominality_ru <- c('оба местоим.', 'местоим. DO',
                      'местоим. IO', 'оба неместоим.')
pronominality <- if(LANG=='ru') pronominality_ru else pronominality_en

animacy_en <- c('both animate', 'animate DO', 'animate IO', 'both inanimate')
animacy_ru <- c('оба одушевл.', 'одушевл. DO',
                'одушевл. IO', 'оба неодушевл.')
animacy <- if(LANG=='ru') animacy_ru else animacy_en


proper_nounness_en <- c('both proper', 'proper DO', 'proper IO', 'both common')
proper_nounness_ru <- c('оба собственные', 'собственное DO',
                'собственное IO', 'оба нарицательные')
proper_nounness <- if(LANG=='ru') proper_nounness_ru else proper_nounness_en

stgr_full$pronominality = with(
  stgr_full,
  ifelse(
    stgr_full$DO_pos=='PRON',
    ifelse(stgr_full$IO_pos=='PRON', pronominality[1], pronominality[2]),
    ifelse(stgr_full$IO_pos=='PRON', pronominality[3], pronominality[4])
  )
)

stgr_full$animacy = with(
  stgr_full,
  ifelse(
    stgr_full$DO_animacy=='од',
    ifelse(stgr_full$IO_animacy=='од', animacy[1], animacy[2]),
    ifelse(stgr_full$IO_animacy=='од', animacy[3], animacy[4])
  )
)

stgr_full$proper_nounness = with(
  stgr_full,
  ifelse(
    stgr_full$DO_is_proper_name==1,
    ifelse(stgr_full$IO_is_proper_name==1, 'both proper', 'proper DO'),
    ifelse(stgr_full$IO_is_proper_name==1, 'proper IO', 'both common')
  )
)

# убрать местоименные
stgr_no_prons <- with(stgr_full, stgr_full[DO_pos!='PRON' & IO_pos!='PRON',])

# исключить вхождения с размеченными эллидированными глаголами (42 из оставшихся)
# исключить вхождения с пустыми одушевлённостями (около 120 из оставшихся)
stgr_no_prons_elidedV <- with(stgr_no_prons,
                              stgr_no_prons[verb !='',])
stgr_no_prons_elidedV_empty <- with(
    stgr_no_prons_elidedV,
    stgr_no_prons_elidedV[DO_animacy!='' & IO_animacy!='',]
)

# фильтровать порядки и ДАТ + Предложные IO
orders = c('V IO DO', 'V DO IO')
order_marked_as_1 = 'V IO DO'
io_head_pos_marked_as_1 = "V"
stopifnot(length(orders)==2)

stgr_dat_adp <- with(stgr_no_prons_elidedV_empty, stgr_no_prons_elidedV_empty[
    (order %in% orders)
    & ((IO_head_pos=='V' & IO_case=='дат') | (IO_head_pos=='PR')),
  ])

stgr_base <- stgr_dat_adp
stgr <- transform(stgr_base,
                  io_length_cat=as.factor(ifelse(
                    IO_length>=10, '>=10', IO_length
                  )),
                  do_length_cat=as.factor(ifelse(
                    DO_length>=10, '>=10', DO_length
                  )),
                  io_dat = ifelse(
                    stgr_base$IO_head_pos==io_head_pos_marked_as_1,
                    1, 0
                  ),
                  verb_inf=as.factor(stgr_base$verb_inf),
                  IO_animacy=as.factor(stgr_base$IO_animacy),
                  DO_animacy=as.factor(stgr_base$DO_animacy),
                  DO_number=as.factor(stgr_base$DO_number),
                  IO_number=as.factor(stgr_base$IO_number),
                  animacy=as.factor(stgr_base$animacy),
                  proper_nounness=as.factor(stgr_base$proper_nounness),
                  order_bin = ifelse(order==order_marked_as_1, 1, 0),
                  order = as.factor(stgr_base$order)
)


stgr_dat = stgr[stgr$IO_head_pos=='V',]  # только дативы
stgr_adp = stgr[stgr$IO_head_pos=='PR',]  # только предлоги

stgr_dat <- transform(stgr_dat,
                      io_depth_sc=scale(IO_depth),
                      do_depth_sc=scale(DO_depth),
                      io_length_sc=scale(IO_length),
                      do_length_sc=scale(DO_length)
                      )

stgr_adp <- transform(stgr_adp,
                      io_depth_sc=scale(IO_depth),
                      do_depth_sc=scale(DO_depth),
                      io_length_sc=scale(IO_length),
                      do_length_sc=scale(DO_length)
                      )


# all orders full
## (full is stgr without elided stuff)
with(stgr_full,
     stgr_full[
       ((IO_head_pos=='V' & IO_case=='дат') | (IO_head_pos=='PR'))
       & (DO_animacy!='' & IO_animacy!='')
       & verb !='',
       ]
    ) %>%
  group_by(order) %>%
  summarise(group_total=n()) %>%
  mutate(percent = round(group_total / sum(group_total) * 100, digits=3)) %>%
  arrange(percent)

## all orders dative
with(stgr_full,
     stgr_full[
       (IO_head_pos=='V' & IO_case=='дат')
       & (DO_animacy!='' & IO_animacy!='')
       & verb !='',
     ]
) %>%
  group_by(order) %>%
  summarise(group_total=n()) %>%
  mutate(percent = round(group_total / sum(group_total) * 100, digits=3)) %>%
  arrange(percent)

# ## dative chi (against what?)
# with(stgr_full,
#      stgr_full[
#        (IO_head_pos=='V' & IO_case=='дат')
#        & order %in% orders
#        & (DO_animacy!='' & IO_animacy!='')
#        & verb !='',
#      ]
# ) %>%


# bar plot pronominality
##  full sample (full is stgr without elided stuff)
pron_full <- with(stgr_full,
     stgr_full[
       ((IO_head_pos=='V' & IO_case=='дат') | (IO_head_pos=='PR'))
       & (DO_animacy!='' & IO_animacy!='')
       & verb !=''
       & (DO_pos=='PRON' | IO_pos=='PRON')
       & (order %in% orders),
     ]) 

plot_bar_pronom_order <- function(data){
  plt <- data %>%
    group_by(pronominality, order) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    mutate(percent = round(count / sum(count) * 100, digits=3)) %>%
    ggplot(aes(x=pronominality, y=count, fill=order, group=order,
               label=paste(count, " (", round(percent,2), "%)", sep="")),
    ) +
    geom_bar(stat = "identity", width=0.8, position=position_dodge(1)) +
    geom_text(size = 3, vjust=-.5, position = position_dodge(1.2)) +
    scale_fill_discrete(name = ifelse(LANG=="ru", "Порядок", "Order")) +
    labs(x = ifelse(LANG=="ru", "Местоимённость", "Pronominality"),
         y = ifelse(LANG=="ru", "Количество", "Count"))
  plt
}

plt <- plot_bar_pronom_order(pron_full)
# plt + theme(legend.position = "none"
# ggsave("./syntagrus/result/graphs/fig_01_1_pronominality_full_nolegend.png", width=5.55, height = 3.29)
# ggsave("./syntagrus/result/graphs/fig_01_1_pronominality_full_nolegend.pdf", width=5.55, height=3.29)

ggsave("./syntagrus/result/graphs/fig_01_1_pronominality_full.png", width=6.05, height = 3.29)
ggsave("./syntagrus/result/graphs/fig_01_1_pronominality_full.pdf", width=6.05, height=3.29)

##  dative sample
pron_dat <- with(stgr_full,
      stgr_full[
        (IO_head_pos=='V' & IO_case=='дат')
        & (DO_animacy!='' & IO_animacy!='')
        & verb !=''
        & (DO_pos=='PRON' | IO_pos=='PRON')
        & (order %in% orders),
      ]) 

plt_dat <- plot_bar_pronom_order(pron_dat)
plt_dat

ggsave("./syntagrus/result/graphs/fig_01_2_pronominality_dat.png", width=6.05, height = 3.29)
ggsave("./syntagrus/result/graphs/fig_01_2_pronominality_dat.pdf", width=6.05, height = 3.29)



# length - depth corr
cor.test(rbind(stgr$DO_length, stgr$IO_length),
         rbind(stgr$DO_depth, stgr$IO_depth),
         alternative="greater")



# animacy and lengths
## num animate
anim_IO_dat_counts = table(stgr_dat[, c("IO_animacy")])
rbind(anim_IO_dat_counts, round(anim_IO_dat_counts / sum(anim_IO_dat_counts), 4))

anim_DO_dat_counts = table(stgr_dat[, c("DO_animacy")])
rbind(anim_DO_dat_counts, round(anim_DO_dat_counts / sum(anim_DO_dat_counts), 4))


## chi-square, mosaic plot
animacies_axes_ru = c("одушевлённость", "глубина дополнения")
animacies_axes_en = c("animacy", "object group length")
animacies_axes = if(LANG=="ru") animacies_axes_ru else animacies_axes_en

animacy_values_ru = c("неодушевлённое", "одушевлённое")
animacy_values_en = c("inanimate", "animate")
animacy_values = if(LANG=="ru") animacy_values_ru else animacy_values_en

make_anim_length_table <- function(stgr_sample){
  stgr_chi_do <- stgr_sample[, c('DO_animacy', 'DO_depth')]
  stgr_chi_io <- stgr_sample[, c('IO_animacy', 'IO_depth')]
  
  names(stgr_chi_do) <- c("animacy", "depth")
  names(stgr_chi_io) <- c("animacy", "depth")
  
  # putting IO and DO together to look at animacy independently of role
  stgr_chi <- rbind(stgr_chi_do, stgr_chi_io)
  
  stgr_chi$length_cat <- as.factor(with(
    stgr_chi,
    ifelse(depth<1, "<1", ifelse(depth>2, ">2", "1-2")))
  )
  stgr_chi_tbl = table(stgr_chi$animacy, stgr_chi$length_cat)
  rownames(stgr_chi_tbl) <- animacy_values
  stgr_chi_tbl <- stgr_chi_tbl[,c(">2", "1-2", "<1")]
  
  stgr_chi_tbl
}

stgr_chi_tbl_dat <- make_anim_length_table(stgr_dat)

stgr_chi_anim_lencat <- chisq.test(stgr_chi_tbl_dat)
stgr_chi_anim_lencat

### Хи-квадрат и V Крамера
cramerV(stgr_chi_tbl_dat, R=1000, bias.correct = TRUE, ci=TRUE, histogram = TRUE, verbose = TRUE, digits = 7)


plot_anim_length_mosaic <- function(stgr_tbl){
  par(mar=c(2,2,1,1)+0.1)  # shrink margins after deleting title
  mosaicplot(
    stgr_tbl, color=T, shade=T,
    main=NULL,
    # main="DO IO animacy length correlation (full sample)",
    xlab=animacies_axes[1],
    ylab=animacies_axes[2]
  )
}

cairo_pdf("./syntagrus/result/graphs/fig_02_1DAT_do_io_anim_depth_ru.pdf",
          width=5, height=4.2, pointsize=16, 
          # width=w, height=h
)
plot_anim_length_mosaic(stgr_chi_tbl_dat)
dev.off()

png("./syntagrus/result/graphs/fig_02_1DAT_do_io_anim_length_1_ru.png",
    res=300, width=w*3.5, height=h*3,
    pointsize = 20)
plot_anim_length_mosaic(stgr_chi_tbl_dat)
dev.off()


# одушевлённость (дативы)
animacy_order_dat = table(stgr_dat[, c("animacy", "order")])
animacy_order_dat

chisq.test(animacy_order_dat)
cramerV(animacy_order_dat, R=1000, bias.correct = TRUE, ci=TRUE, histogram = TRUE, verbose = TRUE, digits = 7)


plot_bar_animacy_order <- function(data){
  plt <- data %>%
    group_by(animacy, order) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    mutate(percent = round(count / sum(count) * 100, digits=3)) %>%
    ggplot(aes(x=animacy, y=count, fill=order, group=order,
               label=paste(count, " (", round(percent,2), "%)", sep="")),
    ) +
    geom_bar(stat = "identity", width=0.8, position=position_dodge(1)) +
    geom_text(size = 3, vjust=-.5, position = position_dodge(1.2)) +
    scale_fill_discrete(name = ifelse(LANG=="ru", "Порядок", "Order")) +
    labs(x = ifelse(LANG=="ru", "Одушевлённость", "Animacy"),
         y = ifelse(LANG=="ru", "Количество", "Count")) +
    theme(legend.position = c(0.5, 0.8))
  plt
}

cairo_pdf("./syntagrus/result/graphs/fig_03_nums_order_animacy_dat.pdf",
          width=6.75, height=3.25, pointsize=16, 
          # width=w, height=h
)
plot_bar_animacy_order(stgr_dat)
dev.off()

png("./syntagrus/result/graphs/fig_03_nums_order_animacy_dat.png",
    res=300, width=w*6.75, height=h*3.25,
    pointsize = 20)
plot_bar_animacy_order(stgr_dat)
dev.off()


# длина
## `animacy_len_order.png`
# do_anim_length_ord <- stgr_dat[, c('DO_animacy', 'DO_length', 'order')]
# io_anim_length_ord <- stgr_dat[, c('IO_animacy', 'IO_length', 'order')]
# 
# names(do_anim_length_ord) <- c("animacy", "length", "order")
# names(io_anim_length_ord) <- c("animacy", "length", "order")
# 
# ### putting IO and DO together to look at animacy independently of role
# anim_length_ord <- rbind(do_anim_length_ord, io_anim_length_ord)
# 
# anim_length_ord$length_cat <- as.factor(with(
#   anim_length_ord,
#   ifelse(length<2, "<1", ifelse(length>4, ">4", "2-4")))
# )
# 
# 
# mosaicplot(
#   table(anim_length_ord[, c("animacy", "length_cat", "order")]), color=T, shade=T,
#   main=NULL,
#   # main="DO IO animacy length correlation (full sample)",
#   xlab=animacies_axes[1],
#   ylab=animacies_axes[2]
# )

## горизонтальный длины
### DO
labels_len <- c("Длина", "Length")
labels_dep <- c("Глубина", "Depth")

length_breaks <- c(1,2,3,4,5,6,8,12,42)  # cut-offs
depth_breaks <- c(0,1,2,3,4,5,6,8,12,42)
length_tags <- c("1","2","3","4","5","[6-8)", "[8-12)","[12-42]")
depth_tags <- c("0","1","2","3","4","5","[6-8)", "[8-12)","[12-42]")

hist(stgr_dat$do_length_sc,# breaks = length_breaks,
     freq=TRUE,
     xlab="Длина DO", ylab="Частота", labels=TRUE, main=NULL,
     col="light gray")
hist(stgr_dat$IO_length, breaks = length_breaks, freq=TRUE,
     xlab="Длина IO", ylab="Частота", labels=TRUE)

get_obj_stgr <- function(data, type="DO", param="length"){
  breaks <- if(param=="length") length_breaks else depth_breaks
  tags <- if(param=="length") length_tags else depth_tags
  obj_stgr <- data[, c("order", paste(type, "_", param, sep=""))]
  obj_stgr$param_cat <- cut(obj_stgr[, paste(type, "_", param, sep="")], 
                             breaks=breaks, 
                             include.lowest=TRUE, 
                             right=FALSE, 
                             labels=tags)
  obj_stgr
}

plot_bar_length_order <- function(data, color="grey", param="length"){
  x_labels = if(param=="length") labels_len else labels_dep
  plt <- data %>%
    group_by(param_cat, order) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>%
    mutate(percent = round(count / sum(count) * 100, digits=3)) %>%
    filter(order=="V IO DO") %>%
    ggplot(aes(x=param_cat, y=percent,# fill=order, group=order,
               label=paste(count, " (", round(percent,2), "%)", sep="")),    ) +
    geom_bar(stat = "identity", width=0.8, position=position_dodge(1),
             fill=color) +
    geom_text(size = 4, hjust=1.1, position = position_dodge(1)) +
    scale_fill_discrete(name = ifelse(LANG=="ru", "Порядок", "Order")) +
    scale_y_continuous(labels=percent_format(scale=1)) +
    labs(x = ifelse(LANG=="ru", x_labels[1], x_labels[2]),
         y = ifelse(LANG=="ru", "Доля порядка IO DO", "Ratio of IO DO order")) +
    coord_flip()
  plt
}

params <- c("length", "depth")
types <- c("DO", "IO")

i <- 6
for (param in params){
  for (type in types){
    obj_param <- get_obj_stgr(stgr_dat, type=type, param=param)
    filename = paste("fig_0", i, "_", type, "_", param, "_", "R", sep="")
    
    cairo_pdf(paste("./syntagrus/result/graphs/", filename, ".pdf", sep=""),
              width=6.5, height=4, pointsize=20, 
    )
    print(plot_bar_length_order(obj_param, "grey", param=param))
    dev.off()
    i = i+1
    
    png(paste("./syntagrus/result/graphs/", filename, ".png", sep=""),
        res=300, width=w*6.5, height=h*4,
        pointsize = 20)
    print(plot_bar_length_order(obj_param, "grey", param=param))
    dev.off()
  }
}

### неумещающийся график
param="length"
type="IO"

obj_param <- get_obj_stgr(stgr_dat, type=type, param=param)
filename = paste("fig_0", 7, "_", type, "_", param, "_", "R", sep="")

op <- par(mar = c(4,10,4,2) + 0.1)
cairo_pdf(paste("./syntagrus/result/graphs/", filename, ".pdf", sep=""),
          width=6.75, height=4, pointsize=20, 
)
print(plot_bar_length_order(obj_param, "grey", param=param))
dev.off()
par(op)



# порядки и одушевлённость: дативы
prog_doioanimdat = table(stgr[stgr$IO_head_pos=='V', c('animacy', 'order')])

xi_doioanimprep <- chisq.test(prog_doioanimdat)
xi_doioanimprep
cramerV(prog_doioanimdat, R=1000, bias.correct = TRUE, ci=TRUE, histogram = TRUE, verbose = TRUE, digits = 7)

cairo_pdf("./syntagrus/result/graphs/fig_05_do_io_anim_order_dat_ru.pdf",
    width=6.5, height=4, pointsize = 12)
# png("./syntagrus/result/graphs/do_io_anim_order_dat.png",
#     res=300, width=400*5, height=400*4,
#     pointsize = 14)
op <- par(mar=c(3,3,1,1)+0.1)  # shrink margins after deleting title
mosaicplot(
  prog_doioanimdat[c(2,1,3,4),], color=T, shade=T,
  # main="DO IO animacy order correlation (dative sample)",
  main=NULL,
  xlab = ifelse(LANG=="ru", "Одушевлённость", "Animacy"),
  ylab = "Порядок")
par(op)
dev.off()




# линейная регрессия
### (упомянуть про шкалирование и центрирование)
### (гистограмму длин см. выше)


fm_length <- order_bin ~ io_length_sc + do_length_sc + (1 | verb_inf)
fm_len_animcat_add  <- order_bin ~ animacy + io_length_sc + do_length_sc + (1 | verb_inf)
fm_len_anim_int <- 
  order_bin ~ io_length_sc * IO_animacy + do_length_sc * DO_animacy + (1 | verb_inf)

fm_len_properness_add <-
  order_bin ~ proper_nounness + io_length_sc + do_length_sc + (1 | verb_inf)
fm_len_properness_int <-
  order_bin ~ io_length_sc * IO_is_proper_name + do_length_sc * DO_is_proper_name + (1 | verb_inf)

fm_len_properness_anim_add <-
  order_bin ~ proper_nounness + animacy + io_length_sc + do_length_sc + (1 | verb_inf)
fm_len_properness_anim_int <- 
  order_bin ~ proper_nounness * animacy + io_length_sc + do_length_sc + (1 | verb_inf)


method='nlminb'
nAGQ=20

gm_length <- glmer(fm_length, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                   control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_length)

gm_len_animcat_add <- glmer(fm_len_animcat_add, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                   control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_len_animcat_add)
anova(gm_length, gm_len_animcat_add)

### BEST MODEL!!!
gm_len_anim_int <- glmer(fm_len_anim_int, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                            control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_len_anim_int)
anova(gm_length, gm_len_anim_int)
anova(gm_len_animcat_add, gm_len_anim_int)


gm_len_properness_add <- glmer(fm_len_properness_add, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                         control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_len_properness_add)
anova(gm_length, gm_len_properness_add)

gm_len_properness_int <- glmer(fm_len_properness_int, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                               control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_len_properness_int)
anova(gm_length, gm_len_properness_int)


gm_len_properness_anim_add <- glmer(fm_len_properness_anim_add, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                               control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_len_properness_anim_add)
anova(gm_len_anim_int, gm_len_properness_anim_add)
anova(gm_len_properness_add, gm_len_properness_anim_add)

gm_len_properness_anim_int <- glmer(fm_len_properness_anim_int, data=stgr_dat, family=binomial, nAGQ = nAGQ,
                                    control = glmerControl(optimizer ='optimx', optCtrl=list(method=method)))
summary(gm_len_properness_anim_int)
anova(gm_len_anim_int, gm_len_properness_anim_int)
anova(gm_len_properness_add, gm_len_properness_anim_int)
anova(gm_len_properness_anim_add, gm_len_properness_anim_int)

### regression check
#### при неодушевлённости IO_length almost -2 * DO_length
#### значит ли это что среди неод. длинных IO доля DO IO будет больше
#### чем доля IO DO среди неод. длинных DO?
#### (нужно перекодирование единичных порядков для DO это должно быть DO IO)

len_stgr_dat <- data.frame(stgr_dat)
len_stgr_dat$DO_len_cat <- cut(len_stgr_dat[, "DO_length"], 
                          breaks=length_breaks, 
                          include.lowest=TRUE, 
                          right=FALSE, 
                          labels=length_tags)
len_stgr_dat$IO_len_cat <- cut(len_stgr_dat[, "IO_length"], 
                            breaks=length_breaks, 
                            include.lowest=TRUE, 
                            right=FALSE, 
                            labels=length_tags)

len_stgr_dat %>%
  filter(IO_animacy=="неод") %>%
  group_by(IO_len_cat) %>%
  summarise(
    IO_first_count=sum(order_bin),
    IO_first_ratio=mean(order_bin)) %>%
  ggplot(aes(x=IO_len_cat, y=IO_first_ratio,# fill=order, group=order,
             label=paste(IO_first_count, " (", round(IO_first_ratio,2), "%)", sep=""))) +
  geom_bar(stat = "identity", width=0.8, position=position_dodge(1)) +
  geom_text(size = 4, vjust=-0.5, position = position_dodge(1)) +
  # scale_fill_discrete(name = ifelse(LANG=="ru", "Порядок", "Order")) +
  # scale_y_continuous(labels=percent_format(scale=1)) +
  # labs(x = ifelse(LANG=="ru", x_labels[1], x_labels[2]),
  #      y = ifelse(LANG=="ru", "Доля порядка IO DO", "Ratio of IO DO order")) +
  ylim(0, 1) 

len_stgr_dat %>%
  filter(DO_animacy=="неод") %>%
  mutate(DO_order_bin = as.numeric(!order_bin)) %>%
  group_by(DO_len_cat) %>%
  summarise(
    DO_first_count=sum(DO_order_bin),
    DO_first_ratio=mean(DO_order_bin)) %>%
  ggplot(aes(x=DO_len_cat, y=DO_first_ratio,# fill=order, group=order,
             label=paste(DO_first_count, " (", round(DO_first_ratio,2), "%)", sep=""))) +
  geom_bar(stat = "identity", width=0.8, position=position_dodge(1)) +
  geom_text(size = 4, vjust=-0.5, position = position_dodge(1)) +
  ylim(0,1)

### Кажется правда, что быть длинным для IO хуже!
### доля порядков с IO первым убывает монотонно и совсем низкая для больших значений



# individual verbs
require(tidyverse)

regr_labels_ru = c("незначимо", "значимо")
regr_labels_en = c("insignificant", "significant")
regr_labels = if(LANG=="ru") regr_labels_ru else regr_labels_en

regr_xlabel_ru = "Влияние глагола (\u03b1=%s)"
regr_xlabel_en = "Verb influence (\u03b1=%s)"
regr_xlabel = if(LANG=="ru") regr_xlabel_ru else regr_xlabel_en

regr_xaxis_en = "\u03b2 model coefficient"
regr_xaxis_ru = "\u03b2 коэффициент модели"
regr_xaxis = if(LANG=="ru") regr_xaxis_ru else regr_xaxis_en

regr_yaxis_en = "Verb"
regr_yaxis_ru = "Глагол"
regr_yaxis = if(LANG=="ru") regr_yaxis_ru else regr_yaxis_en


min_freq = 10
freq_verbs <- c(as.character(
  stgr_dat %>%
    count(verb_inf, sort = TRUE) %>%
    filter(n >= min_freq) %>%
    pull(var=verb_inf))
  )


ci95 <- qnorm(0.975)
p_ <- .05
# fm_graph = order_bin ~ verb_inf
fm_graph = order_bin ~ do_length_sc + io_length_sc + verb_inf
fm_anim_graph = order_bin ~ verb_inf + animacy + do_length_sc + io_length_sc
fm_anim_int_graph = order_bin ~ verb_inf + io_length_sc * IO_animacy + do_length_sc * DO_animacy

coefs_to_keep <- c("animacyоба одушевл.", "animacyодушевл. DO", "animacyодушевл. IO",
                   "do_length_sc", "io_length_sc") 
coefs_to_keep <- c(coefs_to_keep, freq_verbs)


stgr_dat %>%
  group_by(verb_inf) %>%
  filter(n() >= min_freq) %>%
  ungroup() %>%
  glm(fm_anim_int_graph, ., family = binomial(link = "logit")) %>%
  broom::tidy() %>%
  mutate(term = str_remove(term, "verb_inf"),
         term = forcats::fct_rev(term),
         p.adjust = p.adjust(p.value, method = "BH"),
         line_colour = p.value < p_) %>%
  ggplot(aes(x = estimate, y = term, colour = line_colour)) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_point() + 
  geom_segment(aes(x = estimate - std.error*ci95, 
                   xend = estimate + std.error*ci95,
                   yend = term)) +
  scale_colour_manual(values = c("black", "red"), 
                      labels=regr_labels,
                      name = sprintf(regr_xlabel, p_)) +
  theme_minimal() +
  xlim(-5, 5) +
  labs(
    title=NULL,
    # title = "Verb influence on object order",
    x = regr_xaxis,
    y = regr_yaxis) +
  theme(legend.position = "bottom")

# ggsave("verbs_influence_position_MORECOEFF.pdf", scale=1.5)
# ggsave("verbs_influence_position_MORECOEFF.png", scale = 1.5)

ggsave("./syntagrus/result/graphs/fig_10_DAT_verbs_influence_position_MORECOEFF_3.pdf",
       width=10, height = 8, units="cm", scale = 1.5, device=cairo_pdf)
ggsave("./syntagrus/result/graphs/fig_10_DAT_verbs_influence_position_MORECOEFF_3.png",
       width=10, height = 8, units="cm", scale = 1.5)

