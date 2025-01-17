library(RColorBrewer)

# define colors
color_aknno <- c(brewer.pal(8, "Set2"), brewer.pal(12, "Set3"), brewer.pal(8, "Pastel2"),
                 brewer.pal(8, "Accent"))[-c(6, 8, 10, 17, 20)]
color_aknno[c(22, 26, 27, 31, 28, 23)] <- c("#2171b5", "#FF7F00", "#238b45", "#E41A1C", "#dd3497", "#984EA3")
color_aknno[c(3, 12, 11, 24, 2, 30)] <- c(color_aknno[c(11, 3)], "#6a51a3", "#cc4c02", "#E6AB02", "#525252")
