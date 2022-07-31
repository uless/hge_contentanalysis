#LSC 875 Content Analysis Final Project
#@Anqi Shao
#Dec, 2021

library(tidyverse)
library(dplyr)
library(stm)
library(data.table)
library(tidytext)
library(plyr)
library(syuzhet)


##### Data Prep #####
#load("hge_news.Rda")
hge_newsdf$PD <- gsub(' January ', '-01-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' February ', '-02-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' March ', '-03-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' April ', '-04-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' May ', '-05-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' June ', '-06-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' July ', '-07-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' August ', '-08-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' September ', '-09-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' October ', '-10-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' November ', '-11-', hge_newsdf$PD)
hge_newsdf$PD <- gsub(' December ', '-12-', hge_newsdf$PD)
hge_newsdf$PD_new <- as.Date(hge_newsdf$PD, format = "%d-%m-%Y")

hge_newsdf$TD <- gsub("License this article from Dow Jones Reprint Service.*", "", hge_newsdf$TD)

hge_newsdf$day <- as.numeric(hge_newsdf$PD_new)
min(hge_newsdf$day) #15344
hge_newsdf$Day <- hge_newsdf$day - 15343
hge_newsdf <- hge_newsdf %>% select(doc_id,TD,SN,PD_new,Day)

#Nov.25 and 24, nothing
#Nov.26, 8 pieces; Nov. 27, 11 pieces
t <- hge_newsdf%>%filter(PD_new=="2018-11-23")

hge_newsdf <- hge_newsdf %>%
  mutate(baby = ifelse(Day > 2517, "After","Before"))

newsrank<-hge_newsdf %>% count(Day)
i <- ggplot(newsrank, aes(Day, n))
i + geom_step(direction = "hv")


##### STM total#####
processed <- textProcessor(hge_newsdf$TD, 
                           metadata = hge_newsdf,
                           removestopwords = TRUE)
out <- prepDocuments(processed$documents, 
                     processed$vocab, 
                     processed$meta,
                     lower.thresh = 15) 
topicnum35 <- searchK(out$documents, out$vocab, K = c(5:35), 
                      prevalence = ~ SN + baby, 
                      data = out$meta)

plot(topicnum35)

hge_topics <- stm(out$documents,out$vocab, 
                  K=10, prevalence = ~ SN + baby,
                  #max.em.its = 500,
                  #gamma.prior='L1',
                  seed = 53705,
                  data=out$meta, 
                  init.type = "Spectral")

plot(hge_topics,labeltype = c("frex"),frexw = 0.5,n=5,
     main = "10 Prevalent Topics on Human Gene Editing 2012-2019")


hge_effects <- estimateEffect(1:10 ~ SN + baby, hge_topics, meta = out$meta,
                         uncertainty = "Global")

summary(hge_topics)

summary(hge_effects)


plot(
  hge_effects,
  covariate = "baby",
  topics = c(2,5,6,7,8,9),
  model = hge_effects,
  method = "difference",
  cov.value1 = "After",
  cov.value2 = "Before",
  xlab = "Before the CRISPR baby         ...        On and after the CRISPR baby",
  main = "Effect of CRISPR baby Time (2018-11-25) on topic prevalence for news coverage on human gene editing, 
  with mean and 95% confidence intervals",
  xlim = c(-0.25, 0.25),
  labeltype = "custom",
  custom.labels = c("Topic 2: mosquito, crop, farmer plant, label",
                    "Topic 5: crispr, embryo, crisprca, patent, edit",
                    "Topic 6: trump, trump', budget, republican, senat",
                    "Topic 7: tumor, immun, sickl, cancer, blood",
                    "Topic 8: biogen, novarti, price, spark, drug",
                    "Topic 9: mice, sequenc, mutat, brain, chromosom"))

likeinfer <-as.data.frame(hge_topics$theta)
colnames(likeinfer)<-c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", 
                       "Topic6", "Topic7", "Topic8", "Topic9", "Topic10")

likeinfer<-as.data.frame(likeinfer)
row_handler <- function(row.data){  
  index <- which(row.data == max(row.data)) 
  out <- names(row.data[index])
  return(out)}
inferred<-likeinfer %>%
  by_row(..f = row_handler, .collate = "rows", .to = "InferTopic")
df_infer<-bind_cols(out$meta, inferred$InferTopic)
names(df_infer)[7] <- "Topic"

#word count
df_infer$wcount <- sapply(df_infer$TD, 
                               function(x) length(unlist(strsplit(as.character(x), "\\W+"))))

df_infer <- df_infer%>%
  mutate(article_type = ifelse(wcount > 800, "in-depth reports",
                               ifelse(wcount < 200, "notices","news")))

table(df_infer$article_type)

df_plot <- df_infer%>%select(doc_id,Topic,article_type,PD_new)

##### STM by the date #####
hge_before <- hge_newsdf%>%filter(baby == "Before")
hge_after <- hge_newsdf%>%filter(baby == "After")

##### [STM before baby] #####
processed_before <- textProcessor(hge_before$TD, 
                           metadata = hge_before,
                           removestopwords = TRUE)
out_before <- prepDocuments(processed_before$documents, 
                     processed_before$vocab, 
                     processed_before$meta,
                     lower.thresh = 15) 
topicnum_before <- searchK(out_before$documents, out_before$vocab, K = c(5:20), 
                      prevalence = ~ SN + Day, 
                      data = out_before$meta)

plot(topicnum_before)
#11 topics most optimal

topics_before <- stm(out_before$documents,out_before$vocab, 
                  K=11, prevalence = ~ SN + Day,
                  #max.em.its = 500,
                  #gamma.prior='L1',
                  seed = 53705,
                  data=out_before$meta, 
                  init.type = "Spectral")

plot(topics_before,labeltype = c("frex"),frexw = 0.5,n=5,
     main = "11 Prevalent Topics on Human Gene Editing Before the CRISPR baby scandal")

before_effects <- estimateEffect(1:11 ~ SN + Day, topics_before, meta = out_before$meta,
                              uncertainty = "Global")
summary(before_effects)

#topic 2 -
#topic 5 +
#topic 8 +
#topic 9 -
#topic 11 +

likeinfer <-as.data.frame(topics_before$theta)
colnames(likeinfer)<-c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", 
                       "Topic6", "Topic7", "Topic8", "Topic9", "Topic10","Topic11")

likeinfer<-as.data.frame(likeinfer)
row_handler <- function(row.data){  
  index <- which(row.data == max(row.data)) 
  out <- names(row.data[index])
  return(out)}
inferred<-likeinfer %>%
  by_row(..f = row_handler, .collate = "rows", .to = "InferTopic")
dfbeforeinfer<-bind_cols(out_before$meta, inferred$InferTopic)
names(dfbeforeinfer)[7] <- "Topic"

#word count
dfbeforeinfer$wcount <- sapply(dfbeforeinfer$TD, 
                              function(x) length(unlist(strsplit(as.character(x), "\\W+"))))

dfbeforeinfer <- dfbeforeinfer%>%
  mutate(article_type = ifelse(wcount > 800, "in-depth reports",
                               ifelse(wcount < 200, "notices","news")))

table(dfbeforeinfer$article_type)

dfbefore_plot <- dfbeforeinfer%>%select(doc_id,Topic,article_type,PD_new)
                   

##### [STM after baby] #####
processed_after <- textProcessor(hge_after$TD, 
                                  metadata = hge_after,
                                  removestopwords = TRUE)
out_after <- prepDocuments(processed_after$documents, 
                            processed_after$vocab, 
                            processed_after$meta,
                            lower.thresh = 15) 
topicnum_after <- searchK(out_after$documents, out_after$vocab, K = c(3:15), 
                             prevalence = ~ SN + Day, 
                             data = out_after$meta)

plot(topicnum_after)
#8 topics most optimal

unnested <- NULL
nrc_lexicon <- NULL

topics_after <- stm(out_after$documents,out_after$vocab, 
                     K=8, prevalence = ~ SN + Day,
                     #max.em.its = 500,
                     #gamma.prior='L1',
                     seed = 53705,
                     data=out_after$meta, 
                     init.type = "Spectral")

plot(topics_after,labeltype = c("frex"),frexw = 0.5,n=5,
     main = "8 Prevalent Topics on Human Gene Editing Before the CRISPR Baby Scandal")

after_effects <- estimateEffect(1:6 ~ SN + Day, topics_after, meta = out_after$meta,
                                 uncertainty = "Global")
summary(after_effects)
#topic 4-

likeinfer <-as.data.frame(topics_after$theta)
colnames(likeinfer)<-c("Topic1", "Topic2", "Topic3", "Topic4", "Topic5", 
                       "Topic6", "Topic7", "Topic8", "Topic9", "Topic10","Topic11")

likeinfer<-as.data.frame(likeinfer)
row_handler <- function(row.data){  
  index <- which(row.data == max(row.data)) 
  out <- names(row.data[index])
  return(out)}
inferred<-likeinfer %>%
  by_row(..f = row_handler, .collate = "rows", .to = "InferTopic")
dfafterinfer<-bind_cols(out_after$meta, inferred$InferTopic)
names(dfafterinfer)[7] <- "Topic"

#word count
dfafterinfer$wcount <- sapply(dfafterinfer$TD, 
                               function(x) length(unlist(strsplit(as.character(x), "\\W+"))))

dfafterinfer <- dfafterinfer%>%
  mutate(article_type = ifelse(wcount > 800, "in-depth reports",
                               ifelse(wcount < 200, "notices","news")))

table(dfafterinfer$article_type)

dfafter_plot <- dfafterinfer%>%select(doc_id,Topic,article_type,PD_new)


##### Sentiment #####
nrc_lexicon <- get_sentiments("nrc")

unnested <- hge_newsdf %>%
  unnest_tokens(word,TD) %>%  # unnest the words
  left_join(nrc_lexicon)%>%     # join with the lexicon to have sentiments
  left_join(hge_newsdf)

table_sentiment <- as.data.frame.matrix(table(unnested$doc_id, unnested$sentiment))
table_sentiment <- tibble::rowid_to_column(table_sentiment, "index")

names(table_sentiment)[1] <- 'doc_id'

hge_sentiment <- dplyr::full_join(hge_newsdf, table_sentiment, by = "doc_id")
t<- tail(hge_sentiment)

mean_sentiment <- hge_sentiment%>%
  select(PD_new,Day,c(anger:trust))%>%
  group_by(month = lubridate::floor_date(PD_new, "month")) %>%
  summarise_all("mean")

mean_sentiment<-mean_sentiment %>% mutate(across(where(is.numeric), scale))

ggplot(mean_sentiment, aes(PD_new)) + 
  geom_line(aes(y = anger, color = "anger")) + 
  geom_line(aes(y = anticipation, color = "anticipation"))+
  geom_line(aes(y = disgust, color = "disgust")) + 
  geom_line(aes(y = fear, color = "fear")) + 
  geom_line(aes(y = joy, color = "joy")) + 
  geom_line(aes(y = sadness, color = "sadness")) + 
  geom_line(aes(y = surprise, color = "surprise")) + 
  geom_line(aes(y = trust, color = "trust"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "Sentiment evolving in media coverage of human gene editing") + 
  labs(color='Sentiment') +
  ylab("sentiment scores") + xlab("Date")

senti_reg <- full_join(mean_sentiment,hge_newsdf,by = "PD_new")

summary(aov(mean_sentiment$trust~mean_sentiment$anticipation))

mean(mean_sentiment$anticipation)

summary(aov(senti_reg$trust ~ senti_reg$baby.y))

patterns <- "believ*"

hge_newsdf$believe <- str_count(hge_newsdf$TD, "believ*")

table(hge_newsdf$believe)
table(hge_sentiment$trust)

##### moral foundations #####
moral_raw <- read.csv("all-sent.csv")
moral_raw = moral_raw[-1,]

hge_moral <- cbind(hge_newsdf,moral_raw)
hge_moral$PD_new <- substr(hge_moral$PD_new,1,nchar(hge_moral$PD_new)-1)

mean_moral <- hge_moral%>%
  select(PD_new,Day,c(care_p:sanctity_p))%>%
  group_by(month = lubridate::floor_date(PD_new, "month")) %>%
  summarise_all("mean")

ggplot(mean_moral, aes(PD_new)) + 
  geom_line(aes(y = care_p, color = "care_p")) + 
  geom_line(aes(y = fairness_p, color = "fairness_p"))+
  geom_line(aes(y = loyalty_p, color = "loyalty_p")) + 
  geom_line(aes(y = authority_p, color = "authority_p")) + 
  geom_line(aes(y = sanctity_p, color = "sanctity_p")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(title = "Moral foundation evolving in media coverage of human gene editing") + 
  labs(color='Moral foundation') +
  ylab("Moral foundation probability") + xlab("Date")

mean(mean_moral$fairness_p)

summary(aov(mean_moral$care_p ~ mean_moral$fairness_p))
  
##### find examples #####
findexamples <- function(model,text,number, topic){
  findThoughts(model,texts = as.character(text), n = number, topics = topic)$docs[[1]]
}
findexamples(topics_after,hge_after$TD,3,topic=8)

