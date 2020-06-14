df_format <- read_csv ("summary_table.csv")

formattable (df_format,align=c("c","c","c","c","c","c","c","c","l"),
             list(
  Gene_Name= formatter("span",style = ~ style(color = "black",font.weight = "bold")),            
  Start_Sequence1=color_tile ("#E4F9E3","#E4F9E3"),
  Stop_Sequence1=color_tile ("#E4F9E3","#E4F9E3"),
  Gene_Length1=color_tile ("#E4F9E3","#E4F9E3"),
  Start_Sequence2=color_tile ("#D1E6FA","#D1E6FA"),
  Stop_Sequence2=color_tile ("#D1E6FA","#D1E6FA"),
  Gene_Length2=color_tile ("#D1E6FA","#D1E6FA"),
  
  Percentage_Similarity =color_bar ("#FA614B")
))
