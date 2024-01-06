#!-usr-bin-env Rscript

## TFM - Máster universitario de Bioinformática y Bioestadística (UOC, UB)
## Assessing the properties of asymptotic PERMANOVA test through comprehensive simulations in the context of genetic studies
## Aitor Invernón de Campos

## Cretion of TFM's Gantt chart
## Based on Giorgio Comai's ganttrify R project (https:--github.com-giocomai-ganttrify)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Se cargan los paquetes, instalándolos si es necesario
package_list <- c("knitr", "tidyverse", "RColorBrewer", "kableExtra", "scales", "ganttrify", "tidyverse", "egg", "ggpp", "latex2exp")
# Info: ?knitr, ?tidyverse, ?RColorBrewer, ?kableExtra, ?scales, ?ganttrify, ?tidyverse, ?egg, ?ggpp

for (i in package_list){
  if(! i %in% installed.packages()){
    install.packages(i)
  }
  require(i, character.only = TRUE)
  library(i, character.only = TRUE)
}

## Se limpian todos los objetos y la consola
rm(list = ls(all.names = TRUE))
cat("\014")

## Creación del DF para el diagrama de Gantt
wp <- c("PEC1 - Definición y plan de trabajo",
        "PEC1 - Definición y plan de trabajo",
        "PEC1 - Definición y plan de trabajo",
        "PEC1 - Definición y plan de trabajo",
        "PEC1 - Definición y plan de trabajo",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC2 - Desarrollo del trabajo (Fase 1)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC3 - Desarrollo del trabajo (Fase 2)",
        "PEC4 - Cierre de la memoria y presentación",
        "PEC4 - Cierre de la memoria y presentación",
        "PEC4 - Cierre de la memoria y presentación",
        "PEC4 - Cierre de la memoria y presentación",
        "Defensa pública del TFM",
        "Defensa pública del TFM")

activity <- c("Primera tutoría: propuesta oficial de TFM (título, objetivos principales, etc.).", # "PEC1 - Definición y plan de trabajo"
              "Redacción de la PEC1 - Definición y plan de trabajo.", # "PEC1 - Definición y plan de trabajo"
              "Redacción de la memoria: adaptación del contenido de la PEC1 a los correspondientes apartados.", # "PEC1 - Definición y plan de trabajo"
              "Enmiendas basadas en las sugerencias realizadas por los tutores.", # "PEC1 - Definición y plan de trabajo"
              "Entrega de la PEC1.", # "PEC1 - Definición y plan de trabajo"
              "Redacción de la memoria: inicio de la redacción del capítulo Estado del arte.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Abordar las pruebas en R, primer objetivo.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Redacción de la memoria: inclusión de los resultados obtenidos para el primer objetivo.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Continuar con las pruebas en R, segundo objetivo.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Redacción de la memoria: inclusión de los resultados obtenidos para el segundo objetivo.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Valorar, tras los resultados obtenidos para el segundo objetivo, la posibilidad de implementar mejoras en MANTA.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Entrega de la PEC2.", # "PEC2 - Desarrollo del trabajo (Fase 1)"
              "Redacción de la memoria: acabar, si es necesario, la escritura del capítulo Estado del arte.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Seguir con las pruebas en R, tercer objetivo.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Redacción de la memoria: inclusión de los resultados obtenidos para el tercer objetivo.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Pruebas secundarias en R, otros objetivos.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Redacción de la memoria: inclusión, si cabe, de los resultados secundarios obtenidos.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Redacción de la memoria: inicio de los capítulos Discusión y Conclusiones y trabajos futuros.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Entrega de la PEC3.", # "PEC3 - Desarrollo del trabajo (Fase 2)"
              "Finalizar la redacción de la memoria.", # "PEC4 - Cierre de la memoria y presentación"
              "Creación  de la presentación basada en la memoria final.", # "PEC4 - Cierre de la memoria y presentación"
              "Grabación en vídeo de la presentación.", # "PEC4 - Cierre de la memoria y presentación"
              "Entrega de la memoria, la presentación y el vídeo final obtenido.", # "PEC4 - Cierre de la memoria y presentación"
              "Preparación de la defensa", # "Defensa pública del TFM"
              "Fecha límite de la defensa") # "Defensa pública del TFM"

start_date <- as.Date(c('2023-09-26', # "PEC1 - Definición y plan de trabajo"
                        '2023-10-02', # "PEC1 - Definición y plan de trabajo"
                        '2023-10-10', # "PEC1 - Definición y plan de trabajo"
                        '2023-10-15', # "PEC1 - Definición y plan de trabajo"
                        '2023-10-17', # "PEC1 - Definición y plan de trabajo"
                        '2023-10-18', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-10-23', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-11-06', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-11-06', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-11-13', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-11-17', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-11-21', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                        '2023-11-22', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-11-27', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-12-04', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-12-11', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-12-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-12-17', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-12-19', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                        '2023-12-20', # "PEC4 - Cierre de la memoria y presentación"
                        '2024-01-08', # "PEC4 - Cierre de la memoria y presentación"
                        '2024-01-13', # "PEC4 - Cierre de la memoria y presentación"
                        '2024-01-16', # "PEC4 - Cierre de la memoria y presentación"
                        '2024-01-17', # "Defensa pública del TFM"
                        '2024-02-02')) # "Defensa pública del TFM"

end_date <- as.Date(c('2023-10-01', # "PEC1 - Definición y plan de trabajo"
                      '2023-10-10', # "PEC1 - Definición y plan de trabajo"
                      '2023-10-15', # "PEC1 - Definición y plan de trabajo"
                      '2023-10-17', # "PEC1 - Definición y plan de trabajo"
                      '2023-10-17', # "PEC1 - Definición y plan de trabajo"
                      '2023-11-05', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2023-11-05', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2023-11-12', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2023-11-19', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2023-11-19', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2023-11-20', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2023-11-21', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                      '2024-01-14', # "PEC4 - Cierre de la memoria y presentación"
                      '2024-01-14', # "PEC4 - Cierre de la memoria y presentación"
                      '2024-01-15', # "PEC4 - Cierre de la memoria y presentación"
                      '2024-01-16', # "PEC4 - Cierre de la memoria y presentación"
                      '2024-02-02', # "Defensa pública del TFM"
                      '2024-02-02')) # "Defensa pública del TFM"

status <- c('C', # "PEC1 - Definición y plan de trabajo"
            'C', # "PEC1 - Definición y plan de trabajo"
            'C', # "PEC1 - Definición y plan de trabajo"
            'C', # "PEC1 - Definición y plan de trabajo"
            'C', # "PEC1 - Definición y plan de trabajo"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'C', # "PEC2 - Desarrollo del trabajo (Fase 1)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC3 - Desarrollo del trabajo (Fase 2)"
            'I', # "PEC4 - Cierre de la memoria y presentación"
            'I', # "PEC4 - Cierre de la memoria y presentación"
            'I', # "PEC4 - Cierre de la memoria y presentación"
            'I', # "PEC4 - Cierre de la memoria y presentación"
            "I", # "Defensa pública del TFM"
            "I") # "Defensa pública del TFM"

spot_type <- c("T", # "PEC1 - Definición y plan de trabajo"
               "R", # "PEC1 - Definición y plan de trabajo"
               "R", # "PEC1 - Definición y plan de trabajo"
               "V-C", # "PEC1 - Definición y plan de trabajo"
               "E1", # "PEC1 - Definición y plan de trabajo"
               "R", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "PR", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "R", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "PR", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "R", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "V-C", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "E2", # "PEC2 - Desarrollo del trabajo (Fase 1)"
               "R", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "PR", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "R", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "PR", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "R", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "R", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "E3", # "PEC3 - Desarrollo del trabajo (Fase 2)"
               "R", # "PEC4 - Cierre de la memoria y presentación"
               "PF", # "PEC4 - Cierre de la memoria y presentación"
               "GP", # "PEC4 - Cierre de la memoria y presentación"
               "E4", # "PEC4 - Cierre de la memoria y presentación"
               "P", # "Defensa pública del TFM"
               "D") # "Defensa pública del TFM"

               # "Preparación", # "Defensa pública del TFM"
               # "Defensa")) # "Defensa pública del TFM"

spot_T <- "T => Tutoría\n"
spot_R <- "R => Redacción\n"
spot_PR <- "PR => Pruebas en R\n"
spot_E1 <- "E1 => Entrega - PEC1\n"
spot_E2 <- "E2 => Entrega - PEC2\n"
spot_E3 <- "E3 => Entrega - PEC3\n"
spot_EF <- "EF => Entrega Final - PEC4\n"
spot_VC <- "V-C => Valoraciones - Correcciones\n"
spot_PF <- "PF => Presentación Final\n"
spot_GP <- "GP => Grabación Presentación\n"
spot_PD <- "PD => Preparación Defensa\n"
spot_D <- "D => Defensa"
spot_legend <- paste0(spot_T, spot_R, spot_PR, spot_E1, spot_E2, spot_E3, spot_EF, spot_VC, spot_PF, spot_GP, spot_PD, spot_D )

spot_date <- as.Date(c('2023-09-28', # "PEC1 - Definición y plan de trabajo"
                       '2023-10-06', # "PEC1 - Definición y plan de trabajo"
                       '2023-10-12', # "PEC1 - Definición y plan de trabajo"
                       '2023-10-16', # "PEC1 - Definición y plan de trabajo"
                       '2023-10-17', # "PEC1 - Definición y plan de trabajo"
                       '2023-10-27', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-10-29', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-11-09', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-11-12', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-11-16', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-11-18', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-11-21', # "PEC2 - Desarrollo del trabajo (Fase 1)"
                       '2023-12-20', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-20', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-20', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-30', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-30', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-30', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-30', # "PEC3 - Desarrollo del trabajo (Fase 2)"
                       '2023-12-30', # "PEC4 - Cierre de la memoria y presentación"
                       '2024-01-11', # "PEC4 - Cierre de la memoria y presentación"
                       '2024-01-14', # "PEC4 - Cierre de la memoria y presentación"
                       '2024-01-16', # "PEC4 - Cierre de la memoria y presentación"
                       '2024-01-24', # "Defensa pública del TFM"                       
                       '2024-02-02')) # "Defensa pública del TFM"

                       # "Preparación", # "Defensa pública del TFM"
                       # "Defensa")) # "Defensa pública del TFM"

cat("\014") # Se limpia la consola

df <- tibble(data.frame(wp, activity, start_date, end_date, status))
dfspots <- tibble(data.frame(activity, spot_type, spot_date))

cat("\014") # Se limpia la consola


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

## Limpiar todos los gráficos
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

## Aplicación de la función ganttrify()
DateLabelLimits <- c(min(sort(as.Date(union(union(start_date, spot_date), end_date)))),
                     max(sort(as.Date(union(union(start_date, spot_date), end_date)))))

ganttrify(project = df,
          spots = dfspots,
          by_date = TRUE,
          exact_date = TRUE,
          project_start_date = "2023-09-26",
          axis_text_align = "left",
          size_text_relative = 0.75,
          label_wrap = 40,
          font_family = "serif", # families <- c("sans", "serif", "mono", "symbol"),
          alpha_wp = 0.9,
          alpha_activity = 0.6,
          line_end_wp = "round", # alternative values: "butt" or "square"
          line_end_activity = "round", # alternative values: "butt" or "square",
          spot_size_text_relative = 0.75,
          spot_fill = ggplot2::alpha(c("white"), 0.25),
          spot_padding = ggplot2::unit(0.05, "lines"),
          spot_border = NA,
          month_number_label = TRUE,
          month_date_label = FALSE,
          colour_palette = c("#CD6155", "#AED6F1", "#73C6B6", "#C39BD3", "#F8C471"), # Tantos como wp diferentes
          mark_quarters = FALSE) +
  ggplot2::labs(title = expression(bold("Diagrama de Gantt del proyecto")),
                subtitle = expression(atop("TFM - Máster universitario de Bioinformática y Bioestadística (UOC, UB)",
                                           paste(italic(paste("Assessing the properties of asymptotic PERMANOVA test ",
                                           "through comprehensive simulations in the context of genetic studies")))))) +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(size=6,angle=60, hjust=.85),
        axis.title = element_text(size=8, face="bold"),
        plot.margin = margin(.5, .5, .5, .5, "cm")) +
  labs(x = "Fecha") + scale_x_date(limits = DateLabelLimits, date_breaks = "2 days", labels = date_format("%d-%m-%Y")) +
  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label), data = data.frame(x = 0.94, y = 0.97, label = spot_legend))

cat("\014") # Will clear the console.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


p <- ganttrify(project = df,
               spots = dfspots,
               by_date = TRUE,
               exact_date = TRUE,
               project_start_date = "2023-09-26",
               axis_text_align = "left",
               size_text_relative = 0.75,
               label_wrap = 40,
               font_family = "serif", # families <- c("sans", "serif", "mono", "symbol"),
               alpha_wp = 0.9,
               alpha_activity = 0.6,
               line_end_wp = "round", # alternative values: "butt" or "square"
               line_end_activity = "round", # alternative values: "butt" or "square",
               spot_size_text_relative = 0.75,
               spot_fill = ggplot2::alpha(c("white"), 0.25),
               spot_padding = ggplot2::unit(0.05, "lines"),
               spot_border = NA,
               month_number_label = TRUE,
               month_date_label = FALSE,
               colour_palette = c("#CD6155", "#AED6F1", "#73C6B6", "#C39BD3", "#F8C471"), # Tantos como wp diferentes
               mark_quarters = FALSE) +
  ggplot2::labs(title = expression(bold("Diagrama de Gantt del proyecto")),
                subtitle = expression(atop("TFM - Máster universitario de Bioinformática y Bioestadística (UOC, UB)",
                                           paste(italic(paste("Assessing the properties of asymptotic PERMANOVA test ",
                                                              "through comprehensive simulations in the context of genetic studies")))))) +
  theme(plot.title = element_text(hjust=0.5),
        plot.subtitle = element_text(hjust=0.5),
        axis.text.x = element_text(size=6,angle=60, hjust=.85),
        axis.title = element_text(size=8, face="bold"),
        plot.margin = margin(.5, .5, .5, .5, "cm")) +
  labs(x = "Fecha") + scale_x_date(limits = DateLabelLimits, date_breaks = "2 days", labels = date_format("%d-%m-%Y")) +
  ggpp::geom_text_npc(aes(npcx = x, npcy = y, label = label), data = data.frame(x = 0.94, y = 0.97, label = spot_legend))

p_escalado <- set_panel_size(p, width  = unit(24.7, "cm"), height = unit(37, "cm"))
ggsave(p_escalado, filename = "TFMGantt.pdf", width = 29.7, height = 42, units = "cm")

## Limpiar todos los gráficos
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

cat("\014") # Will clear the console.
