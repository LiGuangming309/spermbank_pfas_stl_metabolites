
library(circlize)

################################绘制LOD的图形#################################


library(showtext)
library(sysfonts)


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

# windowsFonts(A = windowsFont("Calibri")) 
###绘制剂量反应图
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/result_lod_picture.pdf",width =80, height =80, pointsize = 80,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(80,80), mar = c(0,0,0,0))  ##,mgp=c(3,2,0)
###修改分面中每个标题的颜色

circos.initialize(letters[1:16], xlim = c(0, 1))
# 创建带自定义边框的 track

circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = "#9ecae1",   ###设置透明度  track.height =cm_h(4) 设置扇形的高度 
             bg.border = adjustcolor("#9ecae1", alpha.f =0.8))
set_track_gap(mm_h(8)) # 2mm 设置间距的宽度
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#FCCDE5", alpha.f = 0.8),
             bg.border = adjustcolor("#FCCDE5", alpha.f = 0.8))
set_track_gap(mm_h(8)) # 0.5cm
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#FB8072", alpha.f = 0.8),
             bg.border = adjustcolor("#FB8072", alpha.f = 0.8))





sectors_1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_1 <- c("3.27 (1.06, 7.22)", "0.5 (0.3, 0.84)",
              "0.08 (0.05, 0.15)","0.24 (0.12, 0.42)",
              "0.02 (0.01, 0.02)","0.79 (0.45, 1.46)",
              "0.004 (0.004, 0.06)","0.2 (0.2, 0.2)",
              "0.001 (0.001, 0.002)","1.82 (1.35, 2.53)",
              "0.45 (0.29, 0.66)","0.33 (0.2, 0.57)",
              "0.02 (0.01, 0.04)","0.1 (0.02, 0.63)",
              "0.004 (0.004, 0.06)","0.08 (0.02, 0.16)")

for (i in seq_along(sectors_1)) {
  highlight.sector(
    sectors_1[i],
    track.index = 1,
    text = labels_1[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}




sectors_2 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_2 <- c("0.2 (0.1, 0.43)", "0.03 (0.02, 0.04)",
              "0.01 (0.01, 0.01)","0.02 (0.01, 0.03)",
              "0.07 (0.04, 0.1)#","0.02 (0.01, 0.04)",     ###* 代表×e-10  ,# 代表×e-100
              "0.004 (0.004, 0.01)","<LOD",
              "<LOD","1.13 (0.08, 0.21)",
              "0.03 (0.03, 0.03)","0.01 (0.004, 0.02)",
              "0.02 (0.02, 0.02)","0.02 (0.02, 0.02)",
              "<LOD","<LOD")

for (i in seq_along(sectors_2)) {
  highlight.sector(
    sectors_2[i],
    track.index = 2,
    text = labels_2[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}




sectors_3 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_3 <- c("0.0516", "0.0034", "0.0092", "0.0035", "0.0005", "0.0005", "0.0049", "0.0031", 
              "0.0015", "0.0068", "0.045", "0.004", "0.0026", "0.0319", "0.0059", "0.0022")

for (i in seq_along(sectors_3)) {
  highlight.sector(
    sectors_3[i],
    track.index = 3,
    text = labels_3[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



sectors_4 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i",
               "j", "k", "l", "m", "n", "o", 
               "p")
labels_4 <- c("PFOS", "PFHxS", "PFHpS", "PFUdA", "8:2 Cl-PFESA", "6:2 Cl-PFESA", "PFBS", "PFDS", "PFOSA", 
              "PFOA", "PFNA", "PFDA", "PFDoA", "PFBA", "PFTrDA", 
              "6:2 diPAP")

for (i in seq_along(sectors_4)) {
  highlight.sector(
    sectors_4[i],
    track.index = 1,
    text = labels_4[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "75mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



highlight.sector("p", col = "#00FF0040",border =adjustcolor("#00FF0040", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(0.05, 0.01, 0.25, 0.01))   ####中间、上、外、下

highlight.sector(c("j", "k", "l", "m", "n", "o"), col = adjustcolor( "#4292c6", alpha.f = 0.4),border = adjustcolor("#4292c6", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(-1, 0, 0.3, 0))   ####中间、上、外、下

highlight.sector(c("a", "b", "c", "d", "e", "f", "g", "h", "i"), col = adjustcolor("#c6dbef", alpha.f = 0.4),border = adjustcolor("#c6dbef", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(-1, 0, 0.3, 0))   ####中间、上、外、下
# highlight.sector("d", col = NA, border = "red", lwd = 2)
# highlight.sector("e", col = "#0000FF40", track.index = c(2, 3))
# highlight.sector(c("f", "g"), col = NA, border = "green", 
#                  lwd = 2, track.index = c(2, 3), padding = c(0.1, 0.2, 0.1, 0.1))



circos.clear()     ####重置图形内部

par(opar)
dev.off()













# 加载必要的包
library(grid)

# 定义颜色
colors <- c("#8DD3C7", "#FFED6F", "#BEBADA", "#FB8072", "#80B1D3")

# 设置图形区域
grid.newpage()

# 计算正方形的数量和位置
num_squares <- length(colors)
n_col <- 3  # 设置列数为3列
n_row <- ceiling(num_squares / n_col)  # 根据列数计算行数

# 设置正方形的大小
square_size <- unit(1, "inches")  # 每个正方形的边长，单位为英寸
spacing <- unit(0.1, "inches")  # 正方形之间的间隔

# 计算整体的宽度和高度
total_width <- n_col * square_size + (n_col - 1) * spacing
total_height <- n_row * square_size + (n_row - 1) * spacing

# 计算偏移量，使正方形在图形中居中
x_offset <- (unit(1, "npc") - total_width) / 2
y_offset <- (unit(1, "npc") - total_height) / 2

# 绘制正方形
for (i in 1:num_squares) {
  # 计算当前正方形的列和行位置
  row <- floor((i - 1) / n_col)
  col <- (i - 1) %% n_col
  
  # 计算每个正方形的横纵坐标，并进行居中调整
  x_pos <- x_offset + col * (square_size + spacing) + square_size / 2
  y_pos <- y_offset + row * (square_size + spacing) + square_size / 2
  
  # 绘制正方形
  grid.rect(x = unit(x_pos, "npc"),
            y = unit(y_pos, "npc"),
            width = square_size,
            height = square_size,
            gp = gpar(fill = colors[i], col = "white"))
}

















################################绘制LOD的图形--最终使用的代码#################################



library(showtext)
library(sysfonts)


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

# windowsFonts(A = windowsFont("Calibri")) 
###绘制剂量反应图
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/result_lod_picture_final.pdf",width =80, height =80, pointsize = 80,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(80,80), mar = c(0,0,0,0))  ##,mgp=c(3,2,0)
###修改分面中每个标题的颜色

circos.initialize(letters[1:16], xlim = c(0, 1))
# 创建带自定义边框的 track

circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col ="#9ecae1",   ###设置透明度  track.height =cm_h(4) 设置扇形的高度 
             bg.border = "#9ecae1")
set_track_gap(mm_h(8)) # 2mm 设置间距的宽度
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#FCCDE5", alpha.f = 0.8),
             bg.border = "#FCCDE5")
set_track_gap(mm_h(12)) # 0.5cm
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#FB8072", alpha.f = 0.8),
             bg.border = "#FB8072")





sectors_1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_1 <- c("3.27 (1.06, 7.22)", "0.5 (0.3, 0.84)",
              "0.08 (0.05, 0.15)","0.24 (0.12, 0.42)",
              "0.02 (0.01, 0.02)","0.79 (0.45, 1.46)",
              "0.004 (0.004, 0.06)","0.2 (0.2, 0.2)",
              "0.001 (0.001, 0.002)","1.82 (1.35, 2.53)",
              "0.45 (0.29, 0.66)","0.33 (0.2, 0.57)",
              "0.02 (0.01, 0.04)","0.1 (0.02, 0.63)",
              "0.004 (0.004, 0.06)","0.08 (0.02, 0.16)")

for (i in seq_along(sectors_1)) {
  highlight.sector(
    sectors_1[i],
    track.index = 1,
    text = labels_1[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}




sectors_2 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_2 <- c("0.2 (0.1, 0.43)", "0.03 (0.02, 0.04)",
              "0.01 (0.01, 0.01)","0.02 (0.01, 0.03)",
              "0.07 (0.04, 0.1)#","0.02 (0.01, 0.04)",     ###* 代表×e-10  ,# 代表×e-100
              "0.004 (0.004, 0.01)","<LOD",
              "<LOD","1.13 (0.08, 0.21)",
              "0.03 (0.03, 0.03)","0.01 (0.004, 0.02)",
              "0.02 (0.02, 0.02)","0.02 (0.02, 0.02)",
              "<LOD","<LOD")

for (i in seq_along(sectors_2)) {
  highlight.sector(
    sectors_2[i],
    track.index = 2,
    text = labels_2[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}




sectors_3 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_3 <- c("0.0516", "0.0034", "0.0092", "0.0035", "0.0005", "0.0005", "0.0049", "0.0031", 
              "0.0015", "0.0068", "0.045", "0.004", "0.0026", "0.0319", "0.0059", "0.0022")

for (i in seq_along(sectors_3)) {
  highlight.sector(
    sectors_3[i],
    track.index = 3,
    text = labels_3[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



sectors_4 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i",
               "j", "k", "l", "m", "n", "o", 
               "p")
labels_4 <- c("PFOS", "PFHxS", "PFHpS", "PFUdA", "8:2 Cl-PFESA", "6:2 Cl-PFESA", "PFBS", "PFDS", "PFOSA", 
              "PFOA", "PFNA", "PFDA", "PFDoA", "PFBA", "PFTrDA", 
              "6:2 diPAP")

for (i in seq_along(sectors_4)) {
  highlight.sector(
    sectors_4[i],
    track.index = 1,
    text = labels_4[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "75mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



highlight.sector("p", col = "#00FF0040",border =adjustcolor("#00FF0040", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(0.05, 0.01, 0.25, 0.01))   ####中间、上、外、下

highlight.sector(c("j", "k", "l", "m", "n", "o"), col = adjustcolor( "#08519c", alpha.f = 0.4),border = adjustcolor("#08519c", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(-1, 0, 0.3, 0))   ####中间、上、外、下

highlight.sector(c("a", "b", "c", "d", "e", "f", "g", "h", "i"), col = adjustcolor("#6baed6", alpha.f = 0.4),border = adjustcolor("#6baed6", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(-1, 0, 0.3, 0))   ####中间、上、外、下
# highlight.sector("d", col = NA, border = "red", lwd = 2)
# highlight.sector("e", col = "#0000FF40", track.index = c(2, 3))
# highlight.sector(c("f", "g"), col = NA, border = "green", 
#                  lwd = 2, track.index = c(2, 3), padding = c(0.1, 0.2, 0.1, 0.1))



circos.clear()     ####重置图形内部

par(opar)
dev.off()

















font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

# windowsFonts(A = windowsFont("Calibri")) 
###绘制剂量反应图
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/result_lod_picture_final_lengend.pdf",width =18, height =18, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(18,18), mar = c(0,0,0,0))  ##,mgp=c(3,2,0)
# 加载必要的包
library(grid)

# 定义颜色
colors <- c("#9ecae1",
            adjustcolor("#FCCDE5", alpha.f = 0.8),
            adjustcolor("#FB8072", alpha.f = 0.8),
            "#00FF0040",
            adjustcolor( "#08519c", alpha.f = 0.4),
            adjustcolor("#6baed6", alpha.f = 0.4))

# 设置图形区域
grid.newpage()

# 计算正方形的数量和位置
num_squares <- length(colors)
n_col <- 3  # 设置列数为3列
n_row <- ceiling(num_squares / n_col)  # 根据列数计算行数

# 设置正方形的大小
square_size <- unit(1, "inches")  # 每个正方形的边长，单位为英寸
spacing <- unit(0.1, "inches")  # 正方形之间的间隔

# 计算整体的宽度和高度
total_width <- n_col * square_size + (n_col - 1) * spacing
total_height <- n_row * square_size + (n_row - 1) * spacing

# 计算偏移量，使正方形在图形中居中
x_offset <- (unit(1, "npc") - total_width) / 2
y_offset <- (unit(1, "npc") - total_height) / 2

# 绘制正方形
for (i in 1:num_squares) {
  # 计算当前正方形的列和行位置
  row <- floor((i - 1) / n_col)
  col <- (i - 1) %% n_col
  
  # 计算每个正方形的横纵坐标，并进行居中调整
  x_pos <- x_offset + col * (square_size + spacing) + square_size / 2
  y_pos <- y_offset + row * (square_size + spacing) + square_size / 2
  
  # 绘制正方形
  grid.rect(x = unit(x_pos, "npc"),
            y = unit(y_pos, "npc"),
            width = square_size,
            height = square_size,
            gp = gpar(fill = colors[i], col = "white"))
}


par(opar)
dev.off()












################################绘制LOD的图形--最终使用的代码添加LOD百分比#################################



library(showtext)
library(sysfonts)

# 
font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

# windowsFonts(A = windowsFont("Calibri"))
###绘制剂量反应图
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/result_lod_picture_final_ALL.pdf",width =80, height =80, pointsize = 80,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(80,80), mar = c(0,0,0,0))  ##,mgp=c(3,2,0)
##修改分面中每个标题的颜色

circos.initialize(letters[1:16], xlim = c(0, 1))
# 创建带自定义边框的 track

circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col =adjustcolor("#631879FF", alpha.f = 0.3),   ###设置透明度  track.height =cm_h(4) 设置扇形的高度 
             bg.border = "#631879FF")
set_track_gap(mm_h(12)) # 2mm 设置间距的宽度
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#9ecae1", alpha.f = 0.6),
             bg.border = "#9ecae1")
set_track_gap(mm_h(10)) # 0.5cm
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#FB8072", alpha.f = 0.6),
             bg.border = "#FB8072")
set_track_gap(mm_h(8)) # 2mm 设置间距的宽度
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col = adjustcolor("#FCCDE5", alpha.f = 0.8),
             bg.border = "#FCCDE5")
set_track_gap(mm_h(6)) # 0.5cm
circos.track(ylim = c(0, 1),track.height =cm_h(8),bg.col =adjustcolor("#FF410DFF", alpha.f = 0.8),   ###设置透明度  track.height =cm_h(4) 设置扇形的高度 
             bg.border = "#FF410DFF")

# circos.clear()     ####重置图形内部




sectors_1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_1 <- c("3.27 (1.06, 7.22)", "0.5 (0.3, 0.84)",
              "0.08 (0.05, 0.15)","0.24 (0.12, 0.42)",
              "0.02 (0.01, 0.02)","0.79 (0.45, 1.46)",
              "0.004 (0.004, 0.06)","0.2 (0.2, 0.2)",
              "0.001 (0.001, 0.002)","1.82 (1.35, 2.53)",
              "0.45 (0.29, 0.66)","0.33 (0.2, 0.57)",
              "0.02 (0.01, 0.04)","0.1 (0.02, 0.63)",
              "0.004 (0.004, 0.06)","0.08 (0.02, 0.16)")

for (i in seq_along(sectors_1)) {
  highlight.sector(
    sectors_1[i],
    track.index = 1,
    text = labels_1[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}




sectors_2 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_2 <- c("0.2 (0.1, 0.43)", "0.03 (0.02, 0.04)",
              "0.01 (0.01, 0.01)","0.02 (0.01, 0.03)",
              "0.07 (0.04, 0.1)#","0.02 (0.01, 0.04)",     ###* 代表×e-10  ,# 代表×e-100
              "0.004 (0.004, 0.01)","<LOD",
              "<LOD","1.13 (0.08, 0.21)",
              "0.03 (0.03, 0.03)","0.01 (0.004, 0.02)",
              "0.02 (0.02, 0.02)","0.02 (0.02, 0.02)",
              "<LOD","<LOD")

for (i in seq_along(sectors_2)) {
  highlight.sector(
    sectors_2[i],
    track.index = 2,
    text = labels_2[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



sectors_1_1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_1_1 <- c("90.51", "100", "98.53", "94.22", "99.12", "100", "34.57", "13.03", 
                "25.95", "99.8", "99.12", "99.61", "88.25", "52.5", "43.78", "87.07")

for (i in seq_along(sectors_1_1)) {
  highlight.sector(
    sectors_1_1[i],
    track.index = 3,
    text = labels_1_1[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



sectors_2_1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_2_1 <- c("90.19", "98.21", "9.57", "84.57", "57.89", "99.4", "27.15", "<LOD", 
                "<LOD", "99.28", "19.74", "76.08", "18.42", "18.42", "<LOD", "<LOD")

for (i in seq_along(sectors_2_1)) {
  highlight.sector(
    sectors_2_1[i],
    track.index = 4,
    text = labels_2_1[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



sectors_3 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p")
labels_3 <- c("0.0516", "0.0034", "0.0092", "0.0035", "0.0005", "0.0005", "0.0049", "0.0031", 
              "0.0015", "0.0068", "0.045", "0.004", "0.0026", "0.0319", "0.0059", "0.0022")

for (i in seq_along(sectors_3)) {
  highlight.sector(
    sectors_3[i],
    track.index = 5,
    text = labels_3[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "-2mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



sectors_4 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i",
               "j", "k", "l", "m", "n", "o", 
               "p")
labels_4 <- c("PFOS", "PFHxS", "PFHpS", "PFUdA", "8:2 Cl-PFESA", "6:2 Cl-PFESA", "PFBS", "PFDS", "PFOSA", 
              "PFOA", "PFNA", "PFDA", "PFDoA", "PFBA", "PFTrDA", 
              "6:2 diPAP")

for (i in seq_along(sectors_4)) {
  highlight.sector(
    sectors_4[i],
    track.index = 1,
    text = labels_4[i],
    col = NA,  # 不填充背景
    facing = "bending.inside",
    niceFacing = TRUE,
    text.vjust = "75mm",
    cex = 1.4,
    font = 2,
    font.family = "Arial"
  )
}



highlight.sector("p", col = "#00FF0040",border = adjustcolor("#00FF0040", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2,3,4,5), padding = c(0.05, 0.01, 0.25, 0.01))   ####中间、上、外、下

highlight.sector(c("j", "k", "l", "m", "n", "o"), col = adjustcolor( "#08519c", alpha.f = 0.4),border = adjustcolor("#08519c", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(-1, 0, 0.3, 0))   ####中间、上、外、下

highlight.sector(c("a", "b", "c", "d", "e", "f", "g", "h", "i"), col = adjustcolor("#91D1C2FF", alpha.f = 0.4),border = adjustcolor("#91D1C2FF", alpha.f = 0.6),
                 lwd = 2, track.index = c(1,2, 3), padding = c(-1, 0, 0.3, 0))   ####中间、上、外、下
# highlight.sector("d", col = NA, border = "red", lwd = 2)
# highlight.sector("e", col = "#0000FF40", track.index = c(2, 3))
# highlight.sector(c("f", "g"), col = NA, border = "green", 
#                  lwd = 2, track.index = c(2, 3), padding = c(0.1, 0.2, 0.1, 0.1))



circos.clear()     ####重置图形内部
# 
par(opar)
dev.off()
# 
# 


#显示配色方案，以npg为例
# install.packages("scales")
library("ggsci")
library("scales")
pal= pal_npg("nrc")(10)

show_col(pal)


show_col(pal_aaas("default", alpha = 0.6)(10))

show_col(pal_bmj("default", alpha = 0.6)(10))

show_col(pal_jama("default", alpha = 0.6)(10))

show_col(pal_nejm("default", alpha = 0.6)(10))


# install.packages("ggsci")






##########将用到的图像全部重新生成小方块，导出


font_add("A","C:/Users/zhimo/AppData/Local/Microsoft/Windows/Fonts/Arial.ttf")   ####字体存放位置

# windowsFonts(A = windowsFont("Calibri"))
###绘制剂量反应图
cairo_pdf("result/20240407_PFOA_sperm/Multi_linear_regression_20240407/describe/result_lod_picture_final_lengend.pdf",width =18, height =18, pointsize = 20,family = "Arial")
opar <- par(no.readonly = TRUE)
par(pin = c(18,18), mar = c(0,0,0,0))  ##,mgp=c(3,2,0)
加载必要的包
library(grid)

# 定义颜色
colors <- c(adjustcolor("#631879FF", alpha.f = 0.3),
            adjustcolor("#9ecae1", alpha.f = 0.6),
            adjustcolor("#FB8072", alpha.f = 0.6),
            adjustcolor("#FCCDE5", alpha.f = 0.8),
            adjustcolor("#FF410DFF", alpha.f = 0.8),
            adjustcolor("#00FF0040", alpha.f = 0.6),
            adjustcolor( "#08519c", alpha.f = 0.4),
            adjustcolor("#91D1C2FF", alpha.f = 0.4))

# 设置图形区域
grid.newpage()

# 计算正方形的数量和位置
num_squares <- length(colors)
n_col <- 3  # 设置列数为3列
n_row <- ceiling(num_squares / n_col)  # 根据列数计算行数

# 设置正方形的大小
square_size <- unit(1, "inches")  # 每个正方形的边长，单位为英寸
spacing <- unit(0.1, "inches")  # 正方形之间的间隔

# 计算整体的宽度和高度
total_width <- n_col * square_size + (n_col - 1) * spacing
total_height <- n_row * square_size + (n_row - 1) * spacing

# 计算偏移量，使正方形在图形中居中
x_offset <- (unit(1, "npc") - total_width) / 2
y_offset <- (unit(1, "npc") - total_height) / 2

# 绘制正方形
for (i in 1:num_squares) {
  # 计算当前正方形的列和行位置
  row <- floor((i - 1) / n_col)
  col <- (i - 1) %% n_col
  
  # 计算每个正方形的横纵坐标，并进行居中调整
  x_pos <- x_offset + col * (square_size + spacing) + square_size / 2
  y_pos <- y_offset + row * (square_size + spacing) + square_size / 2
  
  # 绘制正方形
  grid.rect(x = unit(x_pos, "npc"),
            y = unit(y_pos, "npc"),
            width = square_size,
            height = square_size,
            gp = gpar(fill = colors[i], col = "white"))
}


par(opar)
dev.off()







