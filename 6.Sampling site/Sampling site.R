china_shp <- "中国省级地图GS（2019）1719号.geojson"
nine <- "九段线GS（2019）1719号.geojson"
china <- sf::read_sf(china_shp)
nine_line <- sf::read_sf(nine)

library(ggspatial)
library(cowplot)
library(ggplot2)
map <- ggplot() + 
    geom_sf(data = china,fill=NA) + 
    geom_sf(data = nine_line,color='gray50',size=.8)+
    coord_sf(ylim = c(-2387082,1654989),crs="+proj=laea +lat_0=40 +lon_0=104")+
    # spatial-aware automagic north arrow
    annotation_north_arrow(location = "tl", which_north = "false",
                           style = north_arrow_fancy_orienteering,
    )+
    #theme_bw()+
    theme(
        text = element_text(size = 18,face = "bold"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey80",size=.2)
        )

nine_map <- ggplot() +
    geom_sf(data = china,fill='NA') + 
    geom_sf(data = nine_line,color='gray70',size=1.)+
    #geom_sf(data = scatter_df_tro,aes(fill=class,size=data),shape=21,colour='black',stroke=.25)+
    coord_sf(ylim = c(-4028017,-1877844),xlim = c(117131.4,2115095),crs="+proj=laea +lat_0=40 +lon_0=104")+
    theme(
        #aspect.ratio = 1.25, #调节长宽比
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="grey10",linetype=1,size=1.),
        plot.margin=unit(c(0,0,0,0),"mm"))


gg_inset_map = ggdraw() +
    draw_plot(map) +
    draw_plot(nine_map, x = 0.89, y = 0.07, width = 0.1, height = 0.3)


pdf(file = "采样位点.pdf",width = 6,height = 6)
gg_inset_map
dev.off()





