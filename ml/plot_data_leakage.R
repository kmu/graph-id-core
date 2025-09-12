library(ggplot2)
library(ggplot2)
# library(wesanderson)



theme.paper <- function(plt, base_size=11){
    plt <- plt + theme_bw(base_size) # + theme(legend.position = "none")
    plt <- plt + theme(

        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.ticks.length = unit(0.05, "cm"),
        # # axis.ticks.margin = unit(0.05, "cm"),
        # axis.ticks = element_line(size=0.2),
        # plot.margin = unit(c(0.3,0.3,0.3,0.3),"cm"),
        axis.text = element_text(size = base_size, colour = "black")
    )
    return(plt)
}

df <- read.csv("train3_df8.csv")
df_test <- df[df$kind != 'train',]

max_v = max(max(df$actual) , max(df$pred))
min_v = min(min(df$actual) , min(df$pred))

print(head(df_test))

plt <- ggplot()+ geom_abline(slope=1, col="grey", size=0.5) + geom_point(data=df_test[df_test$kind!='observed',], aes(x=actual, y=pred), col="#446DF6" , size=1, stroke =0) + geom_point(data=df_test[df_test$kind=='observed',], aes(y=pred, x=actual, ),shape=4, col="#01172F", size=0.6, stroke =0.2)  + xlim(min_v, max_v) + ylim(min_v, max_v)
theme.paper(plt) + ylab("Predicted ")
ggsave("data_leakagae.pdf", w=2.5, h=2.5, dpi=600)
ggsave("data_leakagae.png", w=2.5, h=2.5, dpi=600)

df_train <- df[df$kind == 'train',]

plt <- ggplot()+ geom_abline(slope=1, col="grey", size=0.5) + geom_point(data=df_train, aes(y=pred, x=actual, ),shape=4, col="#01172F", size=0.6, stroke =0.2)  + xlim(min_v, max_v) + ylim(min_v, max_v)
theme.paper(plt) + ylab("Predicted ")
ggsave("data_train.pdf", w=2.5, h=2.5, dpi=600)
ggsave("data_train.png", w=2.5, h=2.5, dpi=600)




plt <- ggplot()+ geom_abline(slope=1, col="black", size=0.5) + geom_point(data=df_test[df_test$kind!='observed',], aes(x=actual, y=pred), col="#888888" , size=1, stroke =0, alpha=0.5)  + xlim(min_v, max_v) + ylim(min_v, max_v)
theme.paper(plt) + ylab("Predicted")
ggsave("data_leakagae_observed.pdf", w=2.5, h=2.5, dpi=600)
ggsave("data_leakagae_observed.png", w=2.5, h=2.5, dpi=600)


plt <- ggplot()+ geom_abline(slope=1, col="black", size=0.5) + geom_point(data=df_test[df_test$kind=='observed',], aes(x=actual, y=pred), col="#4FADF0" , size=1, stroke =0, alpha=0.5)  + xlim(min_v, max_v) + ylim(min_v, max_v)
theme.paper(plt) + ylab("Predicted")
ggsave("data_leakagae_unique.pdf", w=2.5, h=2.5, dpi=600)
ggsave("data_leakagae_unique.png", w=2.5, h=2.5, dpi=600)
