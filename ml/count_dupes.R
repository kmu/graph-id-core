df = read.csv("train3_df.csv")

print(head(df))

length(df$graph_id)
length(unique(df$graph_id))