#!/usr/bin/env Rscript

## TFM - Máster universitario de Bioinformática y Bioestadística (UOC, UB)
## Assessing the properties of asymptotic PERMANOVA test through comprehensive simulations in the context of genetic studies
## Aitor Invernón de Campos

## Script for the evaluation of asymptotic PERMANOVA in complex models
## Model: Y ~ A + B + AB
## Based on Diego Garrido-Martín's script Y~A+B+AB.2.R (https://github.com/dgarrimar/manta-sim/blob/sim0/bin/Y~A%2BB%2BAB.2.R)


# # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # #
# #                           # #
# #  S I M U L A C I O N E S  # #
# #                           # #
# # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # #


# Objetivos
# ---------

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PRIMER OBJETIVO #                                                                               #
# # # # # # # # # #                                                                               #
#                                                                                                 #
# Estudiar la pérdida de potencia de la versión asintótica de PERMANOVA (MANTA) con respecto      #
# a MANOVA y otros métodos, profundizando en la afectación de la variación del nivel α de         #
# significación  considerado sobre la potencia de cada uno.                                       #
#                                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# a) Se comparará la potencia de ambos métodos según la variación de la correlación (c o Cor) tanto
#    para el conjunto de datos simulados sin transformar como bajo las transformaciones consideradas
#    (log, sqrt, escalado de datos por desviación típica y normalización Min- Max) para los tres
#    modelos considerados ("mvnorm", "simplex" y "multinom").

# b) Se repite (a) pero dejando la Cor fija.

# c) Ver cómo varía la potencia  para (a) y (b) para diferentes valores de α ∈ [0.05, 0.01, 0.001].

# <=> "Potencia del método" ≡ para las S simulaciones consideradas, representa la fracción de los
#     p-valores del factor a estudio que se encuentran por debajo del nivel α definido.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# SEGUNDO OBJETIVO  #                                                                             #
# # # # # # # # # # #                                                                             #
#                                                                                                 #
# Comparación MANTA-MANOVA de la potencia para un conjunto de datos no transformado bajo          #
# el modelo "mvnorm".                                                                             #                                     
#                                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# a) Analizar todas las combinaciones de valores posibles de las variables ∆, Var y Cor.

# b) Como en (a) pero forzando la matriz de covarianzas y variando también q.

# <=> "Potencia del método" ≡ para las S simulaciones consideradas, representa la fracción de los
#     p-valores del factor a estudio que se encuentran por debajo del nivel α definido.

# <=> Forzado de la "Matriz de covarianzas":

#     i)   x ≡ matrix f(q) .𝜏. rnorm(q)
#     ii)  eigen(tcrossproduct(x))$values
#     iii) Crear "Matriz de covarianzas" que sustituye la de la simulación de datos de Y .𝜏. f(q ≡ nº de respuestas)
#     iv)  Aplicar MANTA y MANOVA para diferentes valores de q ∈ [3, 4, 5, > 5]
#          <=> (q>5) no debería subir el tiempo de computación en exceso pero sí llevar a "sitios extraños"!!!

# ---------




# # # # # # # # # # # # # # # # # # # # # # # #
#  Multivariate normal distribution "mvnorm"  #
# # # # # # # # # # # # # # # # # # # # # # # #
#       Matriz de correlación homogénea       #
#    (mismo valor Cor fuera de la diagonal)   #
# # # # # # # # # # # # # # # # # # # # # # # #

# Inicialización
# --------------

# Se eliminan todos los objetos y variables del entorno.
rm(list = ls(all.names = TRUE))

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Se indican los archivos de soporte que incluyen las diferentes funciones a utilizar.
# Funciones que se usarán de fx.R: dM, dH, dE, interDist, geodesic, step2distance,
# label, step2h1, sim.simplex, Sim.simplex, sim.mvnorm, Sim.mvnorm.
f = as.character(getwd())
source(sprintf("%s/fx.R", f)) # Funciones creadas por el tutor: Diego Garrido-Martín.
source(sprintf("%s/fx_Aitor.R", f)) # Funciones creadas por el autor del TFM: Aitor Invernón de Campos.

# --------------

# Características y definición de variables
# -----------------------------------------

# => Se usará la función de simulación de datos "Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor)".
# => Se simularán todas las combinaciones Δ - Cor, donde Δ = 0 evalua H0 y Δ > 0 evaluará H1.
# => La opción dejar fijo Δ pero variar n eleva demasiado el tiempo de computación. <= HACER LA PRUEBA MIDIÉNDOLO!!!

# Listados de las variables necesarias y sus equivalencias:
list(# FIJAS: a priori, no deberían influir en los resultados.
  "modelSim" <- "mvnorm", # m_values ∈ c("mvnorm", "simplex", "multinom")
  "a" <- 2, # a_values ∈ c(1:6)
  "b" <- 3, # b_values ∈ c(1:6)
  "u" <- 1, # u_values ∈ c(seq(0.1, 2, 0.1))
  "w" <- "B", # w_values ∈ c("A", "B", "AB")
  "hk" <- 1,
  "k" <- 0, # k_values ∈ c(0, 1)
  "plotheatmap" <- F, # plotheatmap_values ∈ as.logical(c("F", "T"))
  "D" <- "unif-0-1",
  "l" <- 1000,
  "s" <- 0.1,
  "p" <- 1,
  "pdist" <- "norm", # pdist_values ∈ c("norm", "gamma", "beta")
  "x" <- 10,
  # VARIABLES: a estudiar su potencial influencia en los resultados.
  "Alpha_values" <- c(0.05, 0.01, 0.001), # Niveles de significación (α) más comunes.
  "S_values" <- c(1E2, 1E3, 1E4), # Número de simulaciones, aumentar S debería aumentar la precisión.
                                  # Máx. 1E4 si se quiere aumentar la precisión sin elevar el tiempo de computación.
  "n_values" <- c(seq(200, 400, 100)), # Tamaño de la muestra.
                                       # Valor que equilibra la representatividad de los resultados con el tiempo de computación.
  "Var_values" <- c("equal", "unequal", "unequalALT1", "unequalALT2"), # Varianza de las variables de Y.
                                                                       # Ver "Sim.mvnorm" para detalles del cálculo de "vars".
  "q_values" <- c(3), # Número de respuestas .𝜏. para q > 5 "SITIOS EXTRAÑOS"!
                      # Variable <=> "simplex" .𝜏. c(3, 5, 8, 10)
                      #              ∀ q ∈ "qlocstdev.norm.tsv" [c(2, 3, 5, 8, 10, 12, 15, 20, 25)]
                      # No aplica <=> "mvnorm" .𝜏. c(3)
  "loc_values" <- c(NA), # Ubicación del modelo generador "simplex" (p / position / loc).
                         # e.g. Si (loc = 1/3) ≡> "Centro del simplex" | (loc -> 1) ≡> "Vértice del simplex"
                         # Variable <=> "simplex" .𝜏. c(1, 2, 3, 5, 8, 10)
                         #              ∀ loc ∈ "qlocstdev.norm.tsv" [c(1, 2, 3, 5, 8, 10)]
                         # No aplica <=> "mvnorm" .𝜏. c(NA)
  "delta_values" <- c(seq(0, 0.35, 0.0175)), # delta (∆) => Δ = 0 evalúa H0 y Δ > 0 evaluará H1.
                                             # El "step" usado depende de la prueba-error del método para diferentes Δ ≠ 0,
                                             # teniendo en cuenta cómo la [Potencia -> 1] en cada método.
  "Cor_values" <- c(seq(0, 0.8, 0.2)), # Correlación de las variables del conjunto de datos Y (Cor).
                                       # Siempre: 0 <= Cor < 1
                                       # Variable <=> "mvnorm" .𝜏. c(seq(0, 0.8, 0.2))
                                       # No aplica <=> "simplex" .𝜏. c(NA)
  "lambda_values" <- c(NA) # Parámetro lambda (distrib. Poisson).
                           # Variable <=> "multinom" .𝜏. c(seq(200, 1200, 200))
                           # No aplica <=> "mvnorm" y en "simplex" .𝜏. c(NA)
)

list(# Equivalencias de algunas variables:
  "chunk" <- k,
  "DistDef" <- dd <- D,
  "p_dist" <- pdist,
  "heterosk" <- H <- hk,
  "cores" <- x,
  "fx" <- f)

# Número total de simulaciones:
Nmax_simul <- length(Alpha_values)*length(S_values)*length(n_values)*
  length(Var_values)*length(q_values)*length(loc_values)*
  length(delta_values)*length(Cor_values)*length(lambda_values)

# DFs y vectores que almacenarán todos los resultados:
DF_Results <- data.frame()

cat("\014") # Clean console

# Combinación de valores de variables actual:
writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values),
                  "}\n  S = {", toString(S_values),
                  "}\n  n = {", toString(n_values),
                  "}\n  Var = {", toString(Var_values),
                  "}\n  q = {", toString(q_values),
                  "}\n  loc = {", toString(loc_values),                  
                  "}\n  delta = {", toString(delta_values),
                  "}\n  Cor = {", toString(Cor_values),
                  "}\n  lambda = {", toString(lambda_values),
                  "}\n\n Se simularán ",
                  length(Alpha_values)*length(S_values)*length(n_values)*
                    length(Var_values)*length(q_values)*length(loc_values)*
                    length(delta_values)*length(Cor_values)*length(lambda_values),
                  " escenarios bajo el modelo ", modelSim, ".", "\n"))

# -----------------------------------------

# Simulación
# ----------

if (modelSim != ""){
  
  "m" <- modelSim
  
  # Restringimos los valores a simular:
  Alpha_values_sim <- Alpha_values[1]
  S_values_sim <- S_values[2]
  n_values_sim <- n_values[2]
  Var_values_sim <- Var_values[]
  q_values_sim <- q_values[]
  loc_values_sim <- loc_values[]
  delta_values_sim <- delta_values[]
  Cor_values_sim <- Cor_values[]
  lambda_values_sim <- lambda_values[]
  
  cat("\014") # Clean console
  
  writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values_sim),
                    "}\n  S = {", toString(S_values_sim),
                    "}\n  n = {", toString(n_values_sim),
                    "}\n  Var = {", toString(Var_values_sim),
                    "}\n  q = {", toString(q_values_sim),
                    "}\n  loc = {", toString(loc_values_sim),
                    "}\n  delta = {", toString(delta_values_sim),
                    "}\n  Cor = {", toString(Cor_values_sim),
                    "}\n  lambda = {", toString(lambda_values_sim),
                    "}\n\n Se simularán ",
                    length(Alpha_values_sim)*length(S_values_sim)*length(n_values_sim)*
                      length(Var_values_sim)*length(q_values_sim)*length(loc_values_sim)*
                      length(delta_values_sim)*length(Cor_values_sim)*length(lambda_values_sim),
                    " escenarios bajo el modelo ", modelSim, ".", "\n"))
  
  # Se crea el directorio de la simulación según convenga:
  
  results_dir <- "Resultados"
  model_dir <- paste("modelSim_", modelSim)
  
  Simul_count <- 0
  t_S_0 <- Sys.time()
  
  for(Alpha in Alpha_values_sim){
    for(S in S_values_sim){
      for(n in n_values_sim){
        for(Var in Var_values_sim){
          "v" <- Var
          for(q in q_values_sim){
            # Variable dependiente de q:
            "mu" <- rep(10, q)
            for(loc in loc_values_sim){
              "position" <- p <- loc
              for(delta in delta_values_sim){
                "d" <- delta
                for(Cor in Cor_values_sim){
                  "c" <- Cor
                  for(lambda in lambda_values_sim){
                    "Lambda" <- lambda <- l
                    
                    Simul_count <- Simul_count + 1
                    
                    if (modelSim == "mvnorm"){
                      
                      stdev <- NA
                      CompMantaManova_mvnorm()
                      
                    } else if (modelSim == "simplex") {
                      
                      CompMantaManova_simplex()
                      
                    } else if (modelSim == "multinom") {
                      
                      CompMantaManova_multinom()
                      
                    }
                    
                    DF_Results <- rbind(DF_Results, DF_CompPot_res)
                    DF_Results_byDatos <- DF_Results[order(DF_Results$Datos), ]
                    
                    # # # # # # # # # # # # # # # # # # # # # # # # # 
                    #   Almacenaje de los datos de cada simulación  #
                    # # # # # # # # # # # # # # # # # # # # # # # # #
                    
                    sim_path = file.path(getwd(), paste0("Resultados/Modelo ", modelSim),
                                         paste0("Sim. ", format(t_S_0, '%d-%m-%Y %H h %M min')))
                    dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path", sim_path, envir = .GlobalEnv)
                    
                    sim_path_DFsSim = file.path(sim_path, paste0("DFs Simulaciones"))
                    dir.create(sim_path_DFsSim, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path_DFsSim", sim_path_DFsSim, envir = .GlobalEnv)
                    
                    # Guardamos la información necesaria hasta este punto:
                    # => deparse(substitute(df)) solventa un problema a la hora de guardar el DF como CSV.
                    #    (https://stackoverflow.com/questions/37998967)
                    
                    # Guardamos en formato "csv" el DF con los resultados de la combinación actual de variables:
                    write.csv(DF_CompPot_res, file = file.path(file.path(sim_path_DFsSim),
                                                               paste0(deparse(substitute(DF_CompPot_res)),
                                                                      " [m = ", toString(modelSim),
                                                                      ", alpha = ", toString(Alpha),
                                                                      ", S = ", toString(S),
                                                                      ", n = ", toString(n),
                                                                      ", Var = ", toString(Var),
                                                                      ", q = ", toString(q),
                                                                      ", loc = ", toString(loc),
                                                                      ", stdev = ", toString(stdev),
                                                                      ", delta = ", toString(sprintf("%1.3f", delta)),
                                                                      ", Cor = ", toString(Cor),
                                                                      ", lambda = ", toString(lambda),
                                                                      "]", ".csv")), row.names = FALSE)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # DF para realizar las gráficas: se unen las potencias y se añade una columna de método (color de la gráfica).
  MANTA_col <- c(rep("MANTA", dim(DF_Results)[1]))
  MANOVA_col <- c(rep("MANOVA", dim(DF_Results)[1]))
  DF_Results <- add_column(DF_Results, MANTA_col, .after = 12)
  DF_Results <- add_column(DF_Results, MANOVA_col, .after = 15)
  DF_Results_MANTA <- DF_Results[1:15]
  DF_Results_MANOVA <- DF_Results[c(1:12, 16:18)]
  colnames(DF_Results_MANOVA) <- colnames(DF_Results_MANTA)
  DF_Results_Graph <- rbind(DF_Results_MANTA, DF_Results_MANOVA)
  colnames(DF_Results_Graph) <- c("Datos", "Modelo", "alpha", "S", "n", "Var", "q", "loc", "stdev",
                                  "Cor", "delta", "lambda", "Método", "Potencia", "tcomp")
  
  # Cálculo del porcentaje de diferencia simple y del RPD ("Relative Percentage Difference") entre
  # la Potencia de MANTA y MANOVA en cada caso. Se añaden ambas columnas al DF de resultados para su
  # posterior tratamiento.
  Column_pctgDif <- c(abs(DF_Results$Pot_MANTA-DF_Results$Pot_MANOVA)*100)
  Column_RPD <- c((abs(DF_Results$Pot_MANTA-DF_Results$Pot_MANOVA)/
                     ((DF_Results$Pot_MANTA+DF_Results$Pot_MANOVA)/2))*100)
  DF_Results_pctgDif_RPD <- add_column(DF_Results, Column_pctgDif, .after = dim(DF_Results)[2])
  DF_Results_pctgDif_RPD <- add_column(DF_Results_pctgDif_RPD, Column_RPD, .after = dim(DF_Results_pctgDif_RPD)[2])
  colnames(DF_Results_pctgDif_RPD) <- c(colnames(DF_Results), "δ Potencia (%)", "RPD Potencia (%)")
  
  # Guardamos en formato "csv" los DF con los resultados finales:
  write.csv(DF_Results, file = file.path(file.path(sim_path),
                                         paste0(deparse(substitute(DF_Results)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_MANTA, file = file.path(file.path(sim_path),
                                               paste0(deparse(substitute(DF_Results_MANTA)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_MANOVA, file = file.path(file.path(sim_path),
                                                paste0(deparse(substitute(DF_Results_MANOVA)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_Graph, file = file.path(file.path(sim_path),
                                               paste0(deparse(substitute(DF_Results_Graph)), ".csv"))
            , row.names = FALSE)

  
  write.csv(DF_Results_pctgDif_RPD, file = file.path(file.path(sim_path),
                                                paste0(deparse(substitute(DF_Results_pctgDif_RPD)), ".csv"))
            , row.names = FALSE)

}

# ----------

# Resultados
# ----------

# # # # # # # # # # # # # # # # # # # # # # #
# Simulación: "Sim. 31-12-2023 20 h 10 min" #
# # # # # # # # # # # # # # # # # # # # # # #
#     Matriz de correlación homogénea       #
#  (mismo valor Cor fuera de la diagonal)   #
# # # # # # # # # # # # # # # # # # # # # # #

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Directorio de la simulación a representar:
setwd(file.path(getwd(), "Resultados/Modelo mvnorm/Sim. 31-12-2023 20 h 10 min"))
DF_Results_2graph <- read.csv("DF_Results.csv")
DF_Results_Graph_2graph <- read.csv("DF_Results_Graph.csv")

## Para mvnorm forzamos que todos los valores de loc y lambda son NA.
DF_Results_2graph$loc <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$loc <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))

## Para mvnorm no se aplicará ni Yscaled ni YnormMinMax ya que se genera otro estadístico diferente al que se está estudiando
## (p_valor) y que, a parte de que no tiene porque ser invariante a transformaciones, NO SIRVE PARA COMPARAR MANTA-MANOVA.
## Para simplex solo se estudiará MANTA, así que no es necesario corregir ni Ylog ni Yclr.

# => Filtramos las transformaciones no necesarias:
DF_Results_Graph_2graph <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos ==
                                                     c(unique(DF_Results_Graph_2graph$Datos)[1:4]), ]
# => Otros arreglos para la gráfica: etiquetas de Var + formato de Cor.
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "equal"] <- "Equal"
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "unequal"] <- "Unequal Type I"
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "unequalALT1"] <- "Unequal Type II"
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "unequalALT2"] <- "Unequal Type III"
DF_Results_Graph_2graph[,'Cor'] = format(round(DF_Results_Graph_2graph[,'Cor'], 1), nsmall = 1)

# (a) + (b): Gráficas "∆-Potencia" comparativas entre MANTA-MANOVA (Datos sin transformar + Datos transformados).

## => Gráficas tipo Scatter Plot "∆-Potencia" entre MANTA y MANOVA agrupadas según "Cor - Var" manteniendo el
##    "Tipo de datos - S - α" fijos, usando una matriz de correlación homogénea (mismo valor Cor fuera de la diagonal):
Count_graph = 0

for (Datos_graph in unique(DF_Results_Graph_2graph$Datos)){
  for (S_graph in unique(DF_Results_Graph_2graph$S)){
    for (alpha_graph in unique(DF_Results_Graph_2graph$alpha)){
      
      Count_graph = Count_graph + 1
      
      DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos == Datos_graph
                                         & DF_Results_Graph_2graph$S == S_graph
                                         & DF_Results_Graph_2graph$alpha == alpha_graph, ]
      
      # PDF
      graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/∆-Potencia MANTA-MANOVA"))
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      pdf(file = file.path(graph_path,
                           paste0("Grid Var-Cor [Datos = ", toString(Datos_graph),
                                  ", S = ", toString(S_graph),
                                  ", alpha = ", toString(alpha_graph), "].pdf")),
          width = 16, height = 12)
      
      p <- DeltaPotScatPlot_Method_mvnorm_1(DF_temp)
      print(p)
      
      dev.off()
      
      # PNG
      graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/∆-Potencia MANTA-MANOVA"))
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      ppi = 600
      png(file = file.path(graph_path,
                           paste0("Grid Var-Cor [Datos = ", toString(Datos_graph),
                                  ", S = ", toString(S_graph),
                                  ", alpha = ", toString(alpha_graph), "].png")),
          width = 16*ppi, height = 12*ppi, res = ppi)
      
      p <- DeltaPotScatPlot_Method_mvnorm_1(DF_temp)
      print(p)
      
      dev.off()
      
    }
  }
}

cat("\014") # Clean console

# (c): ¿Es MANOVA invariante a las transformaciones, y MANTA?

## => Estudio de la posible invarianza a la transformación de datos de cada método:
Count_graph = 0

for (Metodo_graph in unique(DF_Results_Graph_2graph$Método)){
  for (S_graph in unique(DF_Results_Graph_2graph$S)){
    for (alpha_graph in unique(DF_Results_Graph_2graph$alpha)){
      
      Count_graph = Count_graph + 1
      
      DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Método == Metodo_graph
                                         & DF_Results_Graph_2graph$S == S_graph
                                         & DF_Results_Graph_2graph$alpha == alpha_graph, ]
      
      # PDF
      graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/Estudio invarianza"))
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      pdf(file = file.path(graph_path,
                           paste0("Grid Var-Cor [Método = ", toString(Metodo_graph),
                                  ", S = ", toString(S_graph),
                                  ", alpha = ", toString(alpha_graph), "].pdf")),
          width = 16, height = 12)
      
      p <- DeltaPotScatPlot_Method_mvnorm_3(DF_temp, Metodo_graph)
      print(p)
      
      dev.off()
      
      # PNG
      graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/Estudio invarianza"))
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      ppi = 600
      png(file = file.path(graph_path,
                           paste0("Grid Var-Cor [Método = ", toString(Metodo_graph),
                                  ", S = ", toString(S_graph),
                                  ", alpha = ", toString(alpha_graph), "].png")),
          width = 16*ppi, height = 12*ppi, res = ppi)
      
      p <- DeltaPotScatPlot_Method_mvnorm_3(DF_temp, Metodo_graph)
      print(p)
      
      dev.off()
      
    }
  }
}

cat("\014") # Clean console


### ??? ###

# (d): Estudio de la similitud de resultados entre MANTA y MANOVA (∆Potencia = P_MANTA - P_MANOVA).

## => Mediante la "Diferencia en valor absoluto" entre ambas medidas para todas las simulaciones:
DF_AbsDiffPotMANTAMANOVA <- data.frame(DF_Results_2graph[1:12], abs(DF_Results_2graph[14] - DF_Results_2graph[17]))
DF_AbsDiffPotMANTAMANOVA <- DF_AbsDiffPotMANTAMANOVA[DF_AbsDiffPotMANTAMANOVA$alpha == 0.05 &
                                                       DF_AbsDiffPotMANTAMANOVA$Var == "equal", ]
colnames(DF_AbsDiffPotMANTAMANOVA) <- c(colnames(DF_Results_2graph[1:12]), "∆P MANTA-MANOVA")

## => Mediante la "Diferencia de porcentaje relativa" (RPD) entre ambas medidas para todas las simulaciones.

# ----------


# # # # # # # # # # # # # # # # # # # # # # # #
#  Multivariate normal distribution "mvnorm"  #
# # # # # # # # # # # # # # # # # # # # # # # #
#     Matriz de correlación inhomogénea       #
#  (valores aleatorios fuera de la diagonal)  #
# # # # # # # # # # # # # # # # # # # # # # # #

# Inicialización
# --------------

# Se eliminan todos los objetos y variables del entorno.
rm(list = ls(all.names = TRUE))

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Se indican los archivos de soporte que incluyen las diferentes funciones a utilizar.
# Funciones que se usarán de fx.R: dM, dH, dE, interDist, geodesic, step2distance,
# label, step2h1, sim.simplex, Sim.simplex, sim.mvnorm, Sim.mvnorm.
f = as.character(getwd())
source(sprintf("%s/fx.R", f)) # Funciones creadas por el tutor: Diego Garrido-Martín.
source(sprintf("%s/fx_Aitor.R", f)) # Funciones creadas por el autor del TFM: Aitor Invernón de Campos.

# --------------

# Características y definición de variables
# -----------------------------------------

# => Se usará la función de simulación de datos "Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor)".
# => Se simularán todas las combinaciones Δ - Cor, donde Δ = 0 evalua H0 y Δ > 0 evaluará H1.
# => La opción dejar fijo Δ pero variar n eleva demasiado el tiempo de computación. <= HACER LA PRUEBA MIDIÉNDOLO!!!

# Listados de las variables necesarias y sus equivalencias:
list(# FIJAS: a priori, no deberían influir en los resultados.
  "modelSim" <- "mvnorm", # m_values ∈ c("mvnorm", "simplex", "multinom")
  "a" <- 2, # a_values ∈ c(1:6)
  "b" <- 3, # b_values ∈ c(1:6)
  "u" <- 1, # u_values ∈ c(seq(0.1, 2, 0.1))
  "w" <- "B", # w_values ∈ c("A", "B", "AB")
  "hk" <- 1,
  "k" <- 0, # k_values ∈ c(0, 1)
  "plotheatmap" <- F, # plotheatmap_values ∈ as.logical(c("F", "T"))
  "D" <- "unif-0-1",
  "l" <- 1000,
  "s" <- 0.1,
  "p" <- 1,
  "pdist" <- "norm", # pdist_values ∈ c("norm", "gamma", "beta")
  "x" <- 10,
  # VARIABLES: a estudiar su potencial influencia en los resultados.
  "Alpha_values" <- c(0.05, 0.01, 0.001), # Niveles de significación (α) más comunes.
  "S_values" <- c(1E2, 1E3, 1E4), # Número de simulaciones, aumentar S debería aumentar la precisión.
                                  # Máx. 1E4 si se quiere aumentar la precisión sin elevar el tiempo de computación.
  "n_values" <- c(seq(200, 400, 100)), # Tamaño de la muestra.
                                       # Valor que equilibra la representatividad de los resultados con el tiempo de computación.
  "Var_values" <- c("equal", "unequal", "unequalALT1", "unequalALT2"), # Varianza de las variables de Y.
                                                                       # Ver "Sim.mvnorm" para detalles del cálculo de "vars".
  "q_values" <- c(3), # Número de respuestas .𝜏. para q > 5 "SITIOS EXTRAÑOS"!
                      # Variable <=> "simplex" .𝜏. c(3, 5, 8, 10)
                      #              ∀ q ∈ "qlocstdev.norm.tsv" [c(2, 3, 5, 8, 10, 12, 15, 20, 25)]
                      # No aplica <=> "mvnorm" .𝜏. c(3)
  "loc_values" <- c(NA), # Ubicación del modelo generador "simplex" (p / position / loc).
                         # e.g. Si (loc = 1/3) ≡> "Centro del simplex" | (loc -> 1) ≡> "Vértice del simplex"
                         # Variable <=> "simplex" .𝜏. c(1, 2, 3, 5, 8, 10)
                         #              ∀ loc ∈ "qlocstdev.norm.tsv" [c(1, 2, 3, 5, 8, 10)]
                         # No aplica <=> "mvnorm" .𝜏. c(NA)
  "delta_values" <- c(seq(0, 0.35, 0.0175)), # delta (∆) => Δ = 0 evalúa H0 y Δ > 0 evaluará H1.
                                             # El "step" usado depende de la prueba-error del método para diferentes Δ ≠ 0,
                                             # teniendo en cuenta cómo la [Potencia -> 1] en cada método.
  "Cor_values" <- c(0), # Correlación de las variables del conjunto de datos Y (Cor).
                        # Siempre: 0 <= Cor < 1
                        # Variable <=> "mvnorm" .𝜏. c(seq(0, 0.8, 0.2))
                        # No aplica <=> "simplex" .𝜏. c(NA)
  "lambda_values" <- c(NA) # Parámetro lambda (distrib. Poisson).
                           # Variable <=> "multinom" .𝜏. c(seq(200, 1200, 200))
                           # No aplica <=> "mvnorm" y en "simplex" .𝜏. c(NA)
)

list(# Equivalencias de algunas variables:
  "chunk" <- k,
  "DistDef" <- dd <- D,
  "p_dist" <- pdist,
  "heterosk" <- H <- hk,
  "cores" <- x,
  "fx" <- f)

# Número total de simulaciones:
Nmax_simul <- length(Alpha_values)*length(S_values)*length(n_values)*
  length(Var_values)*length(q_values)*length(loc_values)*
  length(delta_values)*length(Cor_values)*length(lambda_values)

# DFs y vectores que almacenarán todos los resultados:
DF_Results <- data.frame()

cat("\014") # Clean console

# Combinación de valores de variables actual:
writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values),
                  "}\n  S = {", toString(S_values),
                  "}\n  n = {", toString(n_values),
                  "}\n  Var = {", toString(Var_values),
                  "}\n  q = {", toString(q_values),
                  "}\n  loc = {", toString(loc_values),                  
                  "}\n  delta = {", toString(delta_values),
                  "}\n  Cor = {", toString(Cor_values),
                  "}\n  lambda = {", toString(lambda_values),
                  "}\n\n Se simularán ",
                  length(Alpha_values)*length(S_values)*length(n_values)*
                    length(Var_values)*length(q_values)*length(loc_values)*
                    length(delta_values)*length(Cor_values)*length(lambda_values),
                  " escenarios bajo el modelo ", modelSim, ".", "\n"))

# -----------------------------------------

# Simulación
# ----------

if (modelSim != ""){
  
  "m" <- modelSim
  
  # Restringimos los valores a simular:
  Alpha_values_sim <- Alpha_values[1]
  S_values_sim <- S_values[2]
  n_values_sim <- n_values[2]
  Var_values_sim <- Var_values[]
  q_values_sim <- q_values[]
  loc_values_sim <- loc_values[]
  delta_values_sim <- delta_values[]
  Cor_values_sim <- Cor_values[]
  lambda_values_sim <- lambda_values[]
  
  cat("\014") # Clean console
  
  writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values_sim),
                    "}\n  S = {", toString(S_values_sim),
                    "}\n  n = {", toString(n_values_sim),
                    "}\n  Var = {", toString(Var_values_sim),
                    "}\n  q = {", toString(q_values_sim),
                    "}\n  loc = {", toString(loc_values_sim),
                    "}\n  delta = {", toString(delta_values_sim),
                    "}\n  Cor = {", toString(Cor_values_sim),
                    "}\n  lambda = {", toString(lambda_values_sim),
                    "}\n\n Se simularán ",
                    length(Alpha_values_sim)*length(S_values_sim)*length(n_values_sim)*
                      length(Var_values_sim)*length(q_values_sim)*length(loc_values_sim)*
                      length(delta_values_sim)*length(Cor_values_sim)*length(lambda_values_sim),
                    " escenarios bajo el modelo ", modelSim, ".", "\n"))
  
  # Se crea el directorio de la simulación según convenga:
  
  results_dir <- "Resultados"
  model_dir <- paste("modelSim_", modelSim)
  
  Simul_count <- 0
  t_S_0 <- Sys.time()
  
  for(Alpha in Alpha_values_sim){
    for(S in S_values_sim){
      for(n in n_values_sim){
        for(Var in Var_values_sim){
          "v" <- Var
          for(q in q_values_sim){
            # Variable dependiente de q:
            "mu" <- rep(10, q)
            for(loc in loc_values_sim){
              "position" <- p <- loc
              for(delta in delta_values_sim){
                "d" <- delta
                for(Cor in Cor_values_sim){
                  "c" <- Cor
                  for(lambda in lambda_values_sim){
                    "Lambda" <- lambda <- l
                    
                    Simul_count <- Simul_count + 1
                    
                    if (modelSim == "mvnorm"){
                      
                      stdev <- NA
                      CompMantaManova_mvnorm()
                      
                    } else if (modelSim == "simplex") {
                      
                      CompMantaManova_simplex()
                      
                    } else if (modelSim == "multinom") {
                      
                      CompMantaManova_multinom()
                      
                    }
                    
                    DF_Results <- rbind(DF_Results, DF_CompPot_res)
                    DF_Results_byDatos <- DF_Results[order(DF_Results$Datos), ]
                    
                    # # # # # # # # # # # # # # # # # # # # # # # # # 
                    #   Almacenaje de los datos de cada simulación  #
                    # # # # # # # # # # # # # # # # # # # # # # # # #
                    
                    sim_path = file.path(getwd(), paste0("Resultados/Modelo ", modelSim),
                                         paste0("Sim. ", format(t_S_0, '%d-%m-%Y %H h %M min')))
                    dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path", sim_path, envir = .GlobalEnv)
                    
                    sim_path_DFsSim = file.path(sim_path, paste0("DFs Simulaciones"))
                    dir.create(sim_path_DFsSim, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path_DFsSim", sim_path_DFsSim, envir = .GlobalEnv)
                    
                    # Guardamos la información necesaria hasta este punto:
                    # => deparse(substitute(df)) solventa un problema a la hora de guardar el DF como CSV.
                    #    (https://stackoverflow.com/questions/37998967)
                    
                    # Guardamos en formato "csv" el DF con los resultados de la combinación actual de variables:
                    write.csv(DF_CompPot_res, file = file.path(file.path(sim_path_DFsSim),
                                                               paste0(deparse(substitute(DF_CompPot_res)),
                                                                      " [m = ", toString(modelSim),
                                                                      ", alpha = ", toString(Alpha),
                                                                      ", S = ", toString(S),
                                                                      ", n = ", toString(n),
                                                                      ", Var = ", toString(Var),
                                                                      ", q = ", toString(q),
                                                                      ", loc = ", toString(loc),
                                                                      ", stdev = ", toString(stdev),
                                                                      ", delta = ", toString(sprintf("%1.3f", delta)),
                                                                      ", Cor = ", toString(Cor),
                                                                      ", lambda = ", toString(lambda),
                                                                      "]", ".csv")), row.names = FALSE)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # DF para realizar las gráficas: se unen las potencias y se añade una columna de método (color de la gráfica).
  MANTA_col <- c(rep("MANTA", dim(DF_Results)[1]))
  MANOVA_col <- c(rep("MANOVA", dim(DF_Results)[1]))
  DF_Results <- add_column(DF_Results, MANTA_col, .after = 12)
  DF_Results <- add_column(DF_Results, MANOVA_col, .after = 15)
  DF_Results_MANTA <- DF_Results[1:15]
  DF_Results_MANOVA <- DF_Results[c(1:12, 16:18)]
  colnames(DF_Results_MANOVA) <- colnames(DF_Results_MANTA)
  DF_Results_Graph <- rbind(DF_Results_MANTA, DF_Results_MANOVA)
  colnames(DF_Results_Graph) <- c("Datos", "Modelo", "alpha", "S", "n", "Var", "q", "loc", "stdev",
                                  "Cor", "delta", "lambda", "Método", "Potencia", "tcomp")
  
  # Cálculo del porcentaje de diferencia simple y del RPD ("Relative Percentage Difference") entre
  # la Potencia de MANTA y MANOVA en cada caso. Se añaden ambas columnas al DF de resultados para su
  # posterior tratamiento.
  Column_pctgDif <- c(abs(DF_Results$Pot_MANTA-DF_Results$Pot_MANOVA)*100)
  Column_RPD <- c((abs(DF_Results$Pot_MANTA-DF_Results$Pot_MANOVA)/
                     ((DF_Results$Pot_MANTA+DF_Results$Pot_MANOVA)/2))*100)
  DF_Results_pctgDif_RPD <- add_column(DF_Results, Column_pctgDif, .after = dim(DF_Results)[2])
  DF_Results_pctgDif_RPD <- add_column(DF_Results_pctgDif_RPD, Column_RPD, .after = dim(DF_Results_pctgDif_RPD)[2])
  colnames(DF_Results_pctgDif_RPD) <- c(colnames(DF_Results), "δ Potencia (%)", "RPD Potencia (%)")
  
  # Guardamos en formato "csv" los DF con los resultados finales:
  write.csv(DF_Results, file = file.path(file.path(sim_path),
                                         paste0(deparse(substitute(DF_Results)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_MANTA, file = file.path(file.path(sim_path),
                                               paste0(deparse(substitute(DF_Results_MANTA)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_MANOVA, file = file.path(file.path(sim_path),
                                                paste0(deparse(substitute(DF_Results_MANOVA)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_Graph, file = file.path(file.path(sim_path),
                                               paste0(deparse(substitute(DF_Results_Graph)), ".csv"))
            , row.names = FALSE)
  
  
  write.csv(DF_Results_pctgDif_RPD, file = file.path(file.path(sim_path),
                                                     paste0(deparse(substitute(DF_Results_pctgDif_RPD)), ".csv"))
            , row.names = FALSE)
  
}

# ----------

# Resultados
# ----------

# # # # # # # # # # # # # # # # # # # # # # #
# Simulación: "Sim. 01-01-2024 22 h 33 min" #
# # # # # # # # # # # # # # # # # # # # # # #
#    Matriz de correlación inhomogénea      #
# (valores aleatorios fuera de la diagonal) #
# # # # # # # # # # # # # # # # # # # # # # #

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Directorio de la simulación a representar:
setwd(file.path(getwd(), "Resultados/Modelo mvnorm/Sim. 01-01-2024 22 h 33 min"))
DF_Results_2graph <- read.csv("DF_Results.csv")
DF_Results_Graph_2graph <- read.csv("DF_Results_Graph.csv")

## Para mvnorm forzamos que todos los valores de loc y lambda son NA.
DF_Results_2graph$loc <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$loc <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))

## Para mvnorm no se aplicará ni Yscaled ni YnormMinMax ya que se genera otro estadístico diferente al que se está estudiando
## (p_valor) y que, a parte de que no tiene porque ser invariante a transformaciones, NO SIRVE PARA COMPARAR MANTA-MANOVA.
## Para simplex solo se estudiará MANTA, así que no es necesario corregir ni Ylog ni Yclr.

# => Filtramos las transformaciones no necesarias:
DF_Results_Graph_2graph <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos ==
                                                     c(unique(DF_Results_Graph_2graph$Datos)[1:4]), ]
# => Otros arreglos para la gráfica: etiquetas de Var + formato de Cor.
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "equal"] <- "Equal"
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "unequal"] <- "Unequal Type I"
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "unequalALT1"] <- "Unequal Type II"
DF_Results_Graph_2graph$Var[DF_Results_Graph_2graph$Var == "unequalALT2"] <- "Unequal Type III"
DF_Results_Graph_2graph[,'Cor'] = format(round(DF_Results_Graph_2graph[,'Cor'], 1), nsmall = 1)

# (a) + (b): Gráficas "∆-Potencia" comparativas entre MANTA-MANOVA (Datos sin transformar + Datos transformados).

## => Gráficas tipo Scatter Plot "∆-Potencia" entre MANTA y MANOVA agrupadas según la "Var" usada, y manteniendo el
##    "Tipo de datos - S - α" fijos, forzando una matriz de correlación inhomogénea (diferentes valores aleatorios
##    fuera de la diagonal):
Count_graph = 0

for (S_graph in unique(DF_Results_Graph_2graph$S)){
  for (alpha_graph in unique(DF_Results_Graph_2graph$alpha)){
    
    Count_graph = Count_graph + 1
    
    DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$S == S_graph
                                       & DF_Results_Graph_2graph$alpha == alpha_graph, ]
    
    # PDF
    graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/∆-Potencia MANTA-MANOVA (corr aleatoria)"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    pdf(file = file.path(graph_path,
                         paste0("Grid Var-Datos [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].pdf")),
        width = 16, height = 12)
    
    p <- DeltaPotScatPlot_Method_mvnorm_2(DF_temp)
    print(p)
    
    dev.off()
    
    # PNG
    graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/∆-Potencia MANTA-MANOVA (corr aleatoria)"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    ppi = 600
    png(file = file.path(graph_path,
                         paste0("Grid Var-Datos [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].png")),
        width = 16*ppi, height = 12*ppi, res = ppi)
    
    p <- DeltaPotScatPlot_Method_mvnorm_2(DF_temp)
    print(p)
    
    dev.off()
  }
}

cat("\014") # Clean console

# (c): ¿Es en este caso MANOVA invariante a las transformaciones, y MANTA?

## => Estudio de la posible invarianza a la transformación de datos de cada método:
Count_graph = 0

for (S_graph in unique(DF_Results_Graph_2graph$S)){
  for (alpha_graph in unique(DF_Results_Graph_2graph$alpha)){
    
    Count_graph = Count_graph + 1
    
    DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$S == S_graph
                                       & DF_Results_Graph_2graph$alpha == alpha_graph, ]
    
    # PDF
    graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/Estudio invarianza"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    pdf(file = file.path(graph_path,
                         paste0("Grid Var-Método [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].pdf")),
        width = 16, height = 12)
    
    p <- DeltaPotScatPlot_Method_mvnorm_4(DF_temp)
    print(p)
    
    dev.off()
    
    # PNG
    graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/Estudio invarianza"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    ppi = 600
    png(file = file.path(graph_path,
                         paste0("Grid Var-Método [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].png")),
        width = 16*ppi, height = 12*ppi, res = ppi)
    
    p <- DeltaPotScatPlot_Method_mvnorm_4(DF_temp)
    print(p)
    
    dev.off()
    
  }
}

cat("\014") # Clean console

# ----------




# # # # # # # # # # # # # # # # #
#  Simplex algorithm "simplex"  #
# # # # # # # # # # # # # # # # # 
#         ∆ max = 0.35          #
# # # # # # # # # # # # # # # # #

# Inicialización
# --------------

# Se eliminan todos los objetos y variables del entorno.
rm(list = ls(all.names = TRUE))

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Se indican los archivos de soporte que incluyen las diferentes funciones a utilizar.
# Funciones que se usarán de fx.R: dM, dH, dE, interDist, geodesic, step2distance,
# label, step2h1, sim.simplex, Sim.simplex, sim.mvnorm, Sim.mvnorm.
f = as.character(getwd())
source(sprintf("%s/fx.R", f)) # Funciones creadas por el tutor: Diego Garrido-Martín.
source(sprintf("%s/fx_Aitor.R", f)) # Funciones creadas por el autor del TFM: Aitor Invernón de Campos.

# --------------

# Características y definición de variables
# -----------------------------------------

# => Se usará la función de simulación de datos "Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = F, pdist)".
# => Se simularán todas las combinaciones Δ - q - Loc ≡ f(q)  (Cor NO APLICA en este método) .𝜏. primero se hará
#    [q fijo <-> Δ variable] y luego [Δ fijo <-> q variable] ∀ Loc ≡ f(q).
#    Donde Δ = 0 evalua H0 y Δ > 0 evaluará H1. La opción dejar fijo Δ pero variar n eleva demasiado
#    el tiempo de computación. <= HACER LA PRUEBA MIDIÉNDOLO!!!
# => Donde Loc indica la ubicación del modelo generador "simplex" (p / position / loc).
#    e.g. Si (loc = 1/3) ≡> "Centro del simplex" | (loc -> 1) ≡> "Vértice del simplex"

# Listados de las variables necesarias y sus equivalencias:
list(# FIJAS: a priori, no deberían influir en los resultados.
  "modelSim" <- "simplex", # m_values ∈ c("mvnorm", "simplex", "multinom")
  "a" <- 2, # a_values ∈ c(1:6)
  "b" <- 3, # b_values ∈ c(1:6)
  "u" <- 1, # u_values ∈ c(seq(0.1, 2, 0.1))
  "w" <- "B", # w_values ∈ c("A", "B", "AB")
  "hk" <- 1,
  "k" <- 0, # k_values ∈ c(0, 1)
  "plotheatmap" <- F, # plotheatmap_values ∈ as.logical(c("F", "T"))
  "D" <- "unif-0-1",
  "l" <- 1000,
  "s" <- 0.1,
  "p" <- 1,
  "pdist" <- "norm", # pdist_values ∈ c("norm", "gamma", "beta")
  "x" <- 10,
  # VARIABLES: a estudiar su potencial influencia en los resultados.
  "Alpha_values" <- c(0.05, 0.01, 0.001), # Niveles de significación (α) más comunes.
  "S_values" <- c(1E2, 1E3, 1E4), # Número de simulaciones, aumentar S debería aumentar la precisión.
                                  # Máx. 1E4 si se quiere aumentar la precisión sin elevar el tiempo de computación.
  "n_values" <- c(seq(200, 400, 100)), # Tamaño de la muestra.
                                       # Valor que equilibra la representatividad de los resultados con el tiempo de computación.
  "Var_values" <- c("equal", "unequal", "unequalALT1", "unequalALT2"), # Varianza de las variables de Y.
                                                                       # Ver "Sim.mvnorm" para detalles del cálculo de "vars".
  "q_values" <- c(3, 5, 8, 10), # Número de respuestas .𝜏. para q > 5 "SITIOS EXTRAÑOS"!
                                # Variable <=> "simplex" .𝜏. c(3, 5, 8, 10)
                                #              ∀ q ∈ "qlocstdev.norm.tsv" [c(2, 3, 5, 8, 10, 12, 15, 20, 25)]
                                # No aplica <=> "mvnorm" .𝜏. c(3)
  "loc_values" <- c(1, 2, 3, 5, 8, 10), # Ubicación del modelo generador "simplex" (p / position / loc).
                                        # e.g. Si (loc = 1/3) ≡> "Centro del simplex" | (loc -> 1) ≡> "Vértice del simplex"
                                        # Variable <=> "simplex" .𝜏. c(1, 2, 3, 5, 8, 10)
                                        #              ∀ loc ∈ "qlocstdev.norm.tsv" [c(1, 2, 3, 5, 8, 10)]
                                        # No aplica <=> "mvnorm" .𝜏. c(NA)
  "delta_values" <- c(seq(0, 0.025, 0.025/50)), # delta (∆) => Δ = 0 evalúa H0 y Δ > 0 evaluará H1.
                                                # El "step" usado depende de la prueba-error del método para diferentes Δ ≠ 0,
                                                # teniendo en cuenta cómo la [Potencia -> 1] en cada método.
  "Cor_values" <- c(NA), # Correlación de las variables del conjunto de datos Y (Cor).
                         # Siempre: 0 <= Cor < 1
                         # Variable <=> "mvnorm" .𝜏. c(seq(0, 0.8, 0.2))
                         # No aplica <=> "simplex" .𝜏. c(NA)
  "lambda_values" <- c(NA) # Parámetro lambda (distrib. Poisson).
                           # Variable <=> "multinom" .𝜏. c(seq(200, 1200, 200))
                           # No aplica <=> "mvnorm" y en "simplex" .𝜏. c(NA)
)

list(# Equivalencias de algunas variables:
  "chunk" <- k,
  "DistDef" <- dd <- D,
  "p_dist" <- pdist,
  "heterosk" <- H <- hk,
  "cores" <- x,
  "fx" <- f)

# Número total de simulaciones:
Nmax_simul <- length(Alpha_values)*length(S_values)*length(n_values)*
  length(Var_values)*length(q_values)*length(loc_values)*
  length(delta_values)*length(Cor_values)*length(lambda_values)

# DFs y vectores que almacenarán todos los resultados:
DF_Results <- data.frame()

cat("\014") # Clean console

# Combinación de valores de variables actual:
writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values),
                  "}\n  S = {", toString(S_values),
                  "}\n  n = {", toString(n_values),
                  "}\n  Var = {", toString(Var_values),
                  "}\n  q = {", toString(q_values),
                  "}\n  loc = {", toString(loc_values),                  
                  "}\n  delta = {", toString(delta_values),
                  "}\n  Cor = {", toString(Cor_values),
                  "}\n  lambda = {", toString(lambda_values),
                  "}\n\n Se simularán ",
                  length(Alpha_values)*length(S_values)*length(n_values)*
                    length(Var_values)*length(q_values)*length(loc_values)*
                    length(delta_values)*length(Cor_values)*length(lambda_values),
                  " escenarios bajo el modelo ", modelSim, ".", "\n"))

# -----------------------------------------

# Simulación
# ----------

if (modelSim != ""){
  
  "m" <- modelSim
  
  # Restringimos los valores a simular:
  Alpha_values_sim <- Alpha_values[1]
  S_values_sim <- S_values[2]
  n_values_sim <- n_values[2]
  Var_values_sim <- Var_values[1]
  q_values_sim <- q_values[1:2] # c(3, 5, 8, 10) ∀ q ∈ "qlocstdev.norm.tsv" [c(2, 3, 5, 8, 10, 12, 15, 20, 25)]
  loc_values_sim <- loc_values[] # ∀ loc ∈ "qlocstdev.norm.tsv" [c(1, 2, 3, 5, 8, 10)]
                                 # loc_values_q_3_sim <- NULL
                                 # loc_values_q_5_sim <- NULL
  loc_values_q_3_sim <- loc_values[1:4]
  loc_values_q_5_sim <- loc_values[1:3]
  delta_values_sim <- delta_values[]
  Cor_values_sim <- Cor_values[]
  lambda_values_sim <- lambda_values[]
  
  NumMaxSim <- length(Alpha_values_sim)*length(S_values_sim)*length(n_values_sim)*
    length(Var_values_sim)*length(q_values_sim)*length(loc_values_sim)*
    length(delta_values_sim)*length(Cor_values_sim)*length(lambda_values_sim)
  
  # Cuando no se simula el mismo número de loc para cada q:
  if(is.null(loc_values_q_3_sim) & is.null(loc_values_q_5_sim)){
    NumMaxSim <- NumMaxSim
    cat("\014") # Clean console
    writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values_sim),
                      "}\n  S = {", toString(S_values_sim),
                      "}\n  n = {", toString(n_values_sim),
                      "}\n  Var = {", toString(Var_values_sim),
                      "}\n  q = {", toString(q_values_sim),
                      "}\n  loc = {", toString(loc_values_sim),
                      "}\n  delta = {", toString(delta_values_sim),
                      "}\n  Cor = {", toString(Cor_values_sim),
                      "}\n  lambda = {", toString(lambda_values_sim),
                      "}\n\n Se simularán ", NumMaxSim,
                      " escenarios bajo el modelo ", modelSim, ".", "\n"))
  } else{
    NumMaxSim <- NumMaxSim - (abs(length(loc_values_q_3_sim) -
                                    length(loc_values_q_5_sim)))*
      length(delta_values_sim)
    cat("\014") # Clean console
    writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values_sim),
                      "}\n  S = {", toString(S_values_sim),
                      "}\n  n = {", toString(n_values_sim),
                      "}\n  Var = {", toString(Var_values_sim),
                      "}\n  Para q = {", toString(q_values_sim[1]),
                      "} => loc = {", toString(loc_values_q_3_sim),
                      "}\n  Para q = {", toString(q_values_sim[2]),
                      "} => loc = {", toString(loc_values_q_5_sim),
                      "}\n  delta = {", toString(delta_values_sim),
                      "}\n  Cor = {", toString(Cor_values_sim),
                      "}\n  lambda = {", toString(lambda_values_sim),
                      "}\n\n Se simularán ", NumMaxSim,
                      " escenarios bajo el modelo ", modelSim, ".", "\n"))
    }
  
  # Se crea el directorio de la simulación según convenga:
  
  results_dir <- "Resultados"
  model_dir <- paste("modelSim_", modelSim)
  
  Simul_count <- 0
  t_S_0 <- Sys.time()
  
  for(Alpha in Alpha_values_sim){
    for(S in S_values_sim){
      for(n in n_values_sim){
        for(Var in Var_values_sim){
          "v" <- Var
          for(q in q_values_sim){
            # Variable dependiente de q:
            "mu" <- rep(10, q)
                 if(q == 3){
                   loc_values_sim <- loc_values[1:4]
                 } else if(q == 5){
                     loc_values_sim <- loc_values[1:3]
                     }
            for(loc in loc_values_sim){
              "position" <- p <- loc
              for(delta in delta_values_sim){
                "d" <- delta
                for(Cor in Cor_values_sim){
                  "c" <- Cor
                  for(lambda in lambda_values_sim){
                    "Lambda" <- lambda <- l
                    
                    Simul_count <- Simul_count + 1
                    
                    if (modelSim == "mvnorm"){
                      
                      CompMantaManova_mvnorm()
                      
                    } else if (modelSim == "simplex") {
                      
                      CompMantaManova_simplex()
                      
                    } else if (modelSim == "multinom") {
                      
                      CompMantaManova_multinom()
                      
                    }
                    
                    DF_Results <- rbind(DF_Results, DF_CompPot_res)
                    DF_Results_byDatos <- DF_Results[order(DF_Results$Datos), ]
                    
                    # # # # # # # # # # # # # # # # # # # # # # # # # 
                    #   Almacenaje de los datos de cada simulación  #
                    # # # # # # # # # # # # # # # # # # # # # # # # #
                    
                    sim_path = file.path(getwd(), paste0("Resultados/Modelo ", modelSim),
                                         paste0("Sim. ", format(t_S_0, '%d-%m-%Y %H h %M min')))
                    dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path", sim_path, envir = .GlobalEnv)
                    
                    sim_path_DFsSim = file.path(sim_path, paste0("DFs Simulaciones"))
                    dir.create(sim_path_DFsSim, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path_DFsSim", sim_path_DFsSim, envir = .GlobalEnv)
                    
                    # Guardamos la información necesaria hasta este punto:
                    # => deparse(substitute(df)) solventa un problema a la hora de guardar el DF como CSV.
                    #    (https://stackoverflow.com/questions/37998967)
                    
                    # Guardamos en formato "csv" el DF con los resultados de la combinación actual de variables:
                    write.csv(DF_CompPot_res, file = file.path(file.path(sim_path_DFsSim),
                                                               paste0(deparse(substitute(DF_CompPot_res)),
                                                                      " [m = ", toString(modelSim),
                                                                      ", alpha = ", toString(Alpha),
                                                                      ", S = ", toString(S),
                                                                      ", n = ", toString(n),
                                                                      ", Var = ", toString(Var),
                                                                      ", q = ", toString(q),
                                                                      ", loc = ", toString(loc),
                                                                      ", stdev = ", toString(stdev),
                                                                      ", delta = ", toString(sprintf("%1.3f", delta)),
                                                                      ", Cor = ", toString(Cor),
                                                                      ", lambda = ", toString(lambda),
                                                                      "]", ".csv")), row.names = FALSE)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  MANTA_col <- c(rep("MANTA", dim(DF_Results)[1]))
  MANOVA_col <- c(rep("MANOVA", dim(DF_Results)[1]))
  DF_Results <- add_column(DF_Results, MANTA_col, .after = 12)
  DF_Results <- add_column(DF_Results, MANOVA_col, .after = 15)
  DF_Results_MANTA <- DF_Results[1:15]
  DF_Results_MANOVA <- DF_Results[c(1:12, 16:18)]
  colnames(DF_Results_MANOVA) <- colnames(DF_Results_MANTA)
  DF_Results_Graph <- rbind(DF_Results_MANTA, DF_Results_MANOVA)
  colnames(DF_Results_Graph) <- c("Datos", "Modelo", "alpha", "S", "n", "Var", "q", "loc", "stdev",
                                  "Cor", "delta", "lambda", "Método", "Potencia", "tcomp")

  # Guardamos en formato "csv" los DF con los resultados finales:

  write.csv(DF_Results, file = file.path(file.path(sim_path),
                                         paste0(deparse(substitute(DF_Results)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_Graph, file = file.path(file.path(sim_path),
                                         paste0(deparse(substitute(DF_Results_Graph)), ".csv"))
            , row.names = FALSE)
  
}

# ----------

# Resultados
# ----------

# # # # # # # # # # # # # # # # # # # # # # #
# Simulación: "Sim. 21-12-2023 21 h 11 min" #
# # # # # # # # # # # # # # # # # # # # # # #
#              ∆ max = 0.35                 #
# # # # # # # # # # # # # # # # # # # # # # #

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

## Directorio de la simulación a representar:
setwd(file.path(getwd(), "Resultados/Modelo simplex/Sim. 21-12-2023 21 h 11 min"))
# setwd(file.path(getwd(), "Resultados/Modelo simplex/Sim. 22-12-2023 21 h 20 min"))
DF_Results_2graph <- read.csv("DF_Results.csv")
DF_Results_Graph_2graph <- read.csv("DF_Results_Graph.csv")

## Para simplex forzamos que todos los valores de Cor y lambda son NA.
DF_Results_2graph$Cor <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$Cor <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))

## Se estudiará MANTA, teniendo en cuenta los datos sin transformar y las transformaciones log-ratio, sqrt y clr:
DF_Results_2graph <- DF_Results_2graph[DF_Results_2graph$Datos == unique(DF_Results_2graph$Datos)[1:4], ]
DF_Results_Graph_2graph <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos ==
                                                     unique(DF_Results_Graph_2graph$Datos)[1:4], ]
DF_Results_Graph_2graph <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Método ==
                                                     unique(DF_Results_Graph_2graph$Método[1]), ]

# (a): ¿Es MANTA invariante a las transformaciones bajo una distribución de datos simulada
#       mediante el modelo simplex?

## => Grid "q - loc" del estudio de la posible invarianza a la transformación de datos de MANTA:
for (S_graph in unique(DF_Results_Graph_2graph$S)){
  for (alpha_graph in unique(DF_Results_Graph_2graph$alpha)){
    
    Count_graph = Count_graph + 1
    
    DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$S == S_graph
                                       & DF_Results_Graph_2graph$alpha == alpha_graph, ]
    
    # PDF
    graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/Estudio Invarianza MANTA"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    pdf(file = file.path(graph_path,
                         paste0("Grid q-loc [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].pdf")),
        width = 16, height = 12)
    
    p <- DeltaPotScatPlot_Method_simplex_1(DF_temp)
    print(p)
    
    dev.off()
    
    # PNG
    graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/Estudio Invarianza MANTA"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    ppi = 600
    png(file = file.path(graph_path,
                         paste0("Grid q-loc [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].png")),
        width = 16*ppi, height = 12*ppi, res = ppi)
    
    p <- DeltaPotScatPlot_Method_simplex_1(DF_temp)
    print(p)
    
    dev.off()
  }
}

cat("\014") # Clean console

# ----------


# # # # # # # # # # # # # # # # #
#  Simplex algorithm "simplex"  #
# # # # # # # # # # # # # # # # # 
#        ∆ max = 0.025          #
# # # # # # # # # # # # # # # # #

# Inicialización
# --------------

# Se eliminan todos los objetos y variables del entorno.
rm(list = ls(all.names = TRUE))

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Se indican los archivos de soporte que incluyen las diferentes funciones a utilizar.
# Funciones que se usarán de fx.R: dM, dH, dE, interDist, geodesic, step2distance,
# label, step2h1, sim.simplex, Sim.simplex, sim.mvnorm, Sim.mvnorm.
f = as.character(getwd())
source(sprintf("%s/fx.R", f)) # Funciones creadas por el tutor: Diego Garrido-Martín.
source(sprintf("%s/fx_Aitor.R", f)) # Funciones creadas por el autor del TFM: Aitor Invernón de Campos.

# --------------

# Características y definición de variables
# -----------------------------------------

# => Se usará la función de simulación de datos "Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = F, pdist)".
# => Se simularán todas las combinaciones Δ - q - Loc ≡ f(q)  (Cor NO APLICA en este método) .𝜏. primero se hará
#    [q fijo <-> Δ variable] y luego [Δ fijo <-> q variable] ∀ Loc ≡ f(q).
#    Donde Δ = 0 evalua H0 y Δ > 0 evaluará H1. La opción dejar fijo Δ pero variar n eleva demasiado
#    el tiempo de computación. <= HACER LA PRUEBA MIDIÉNDOLO!!!
# => Donde Loc indica la ubicación del modelo generador "simplex" (p / position / loc).
#    e.g. Si (loc = 1/3) ≡> "Centro del simplex" | (loc -> 1) ≡> "Vértice del simplex"

# Listados de las variables necesarias y sus equivalencias:
list(# FIJAS: a priori, no deberían influir en los resultados.
  "modelSim" <- "simplex", # m_values ∈ c("mvnorm", "simplex", "multinom")
  "a" <- 2, # a_values ∈ c(1:6)
  "b" <- 3, # b_values ∈ c(1:6)
  "u" <- 1, # u_values ∈ c(seq(0.1, 2, 0.1))
  "w" <- "B", # w_values ∈ c("A", "B", "AB")
  "hk" <- 1,
  "k" <- 0, # k_values ∈ c(0, 1)
  "plotheatmap" <- F, # plotheatmap_values ∈ as.logical(c("F", "T"))
  "D" <- "unif-0-1",
  "l" <- 1000,
  "s" <- 0.1,
  "p" <- 1,
  "pdist" <- "norm", # pdist_values ∈ c("norm", "gamma", "beta")
  "x" <- 10,
  # VARIABLES: a estudiar su potencial influencia en los resultados.
  "Alpha_values" <- c(0.05, 0.01, 0.001), # Niveles de significación (α) más comunes.
  "S_values" <- c(1E2, 1E3, 1E4), # Número de simulaciones, aumentar S debería aumentar la precisión.
                                  # Máx. 1E4 si se quiere aumentar la precisión sin elevar el tiempo de computación.
  "n_values" <- c(seq(200, 400, 100)), # Tamaño de la muestra.
                                       # Valor que equilibra la representatividad de los resultados con el tiempo de computación.
  "Var_values" <- c("equal", "unequal", "unequalALT1", "unequalALT2"), # Varianza de las variables de Y.
                                                                       # Ver "Sim.mvnorm" para detalles del cálculo de "vars".
  "q_values" <- c(3, 5, 8, 10), # Número de respuestas .𝜏. para q > 5 "SITIOS EXTRAÑOS"!
                                # Variable <=> "simplex" .𝜏. c(3, 5, 8, 10)
                                #              ∀ q ∈ "qlocstdev.norm.tsv" [c(2, 3, 5, 8, 10, 12, 15, 20, 25)]
                                # No aplica <=> "mvnorm" .𝜏. c(3)
  "loc_values" <- c(1, 2, 3, 5, 8, 10), # Ubicación del modelo generador "simplex" (p / position / loc).
                                        # e.g. Si (loc = 1/3) ≡> "Centro del simplex" | (loc -> 1) ≡> "Vértice del simplex"
                                        # Variable <=> "simplex" .𝜏. c(1, 2, 3, 5, 8, 10)
                                        #              ∀ loc ∈ "qlocstdev.norm.tsv" [c(1, 2, 3, 5, 8, 10)]
                                        # No aplica <=> "mvnorm" .𝜏. c(NA)
  "delta_values" <- c(seq(0, 0.025, 0.025/50)), # delta (∆) => Δ = 0 evalúa H0 y Δ > 0 evaluará H1.
                                                # El "step" usado depende de la prueba-error del método para diferentes Δ ≠ 0,
                                                # teniendo en cuenta cómo la [Potencia -> 1] en cada método.
  "Cor_values" <- c(NA), # Correlación de las variables del conjunto de datos Y (Cor).
                         # Siempre: 0 <= Cor < 1
                         # Variable <=> "mvnorm" .𝜏. c(seq(0, 0.8, 0.2))
                         # No aplica <=> "simplex" .𝜏. c(NA)
  "lambda_values" <- c(NA) # Parámetro lambda (distrib. Poisson).
                           # Variable <=> "multinom" .𝜏. c(seq(200, 1200, 200))
                           # No aplica <=> "mvnorm" y en "simplex" .𝜏. c(NA)
)

list(# Equivalencias de algunas variables:
  "chunk" <- k,
  "DistDef" <- dd <- D,
  "p_dist" <- pdist,
  "heterosk" <- H <- hk,
  "cores" <- x,
  "fx" <- f)

# Número total de simulaciones:
Nmax_simul <- length(Alpha_values)*length(S_values)*length(n_values)*
  length(Var_values)*length(q_values)*length(loc_values)*
  length(delta_values)*length(Cor_values)*length(lambda_values)

# DFs y vectores que almacenarán todos los resultados:
DF_Results <- data.frame()

cat("\014") # Clean console

# Combinación de valores de variables actual:
writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values),
                  "}\n  S = {", toString(S_values),
                  "}\n  n = {", toString(n_values),
                  "}\n  Var = {", toString(Var_values),
                  "}\n  q = {", toString(q_values),
                  "}\n  loc = {", toString(loc_values),                  
                  "}\n  delta = {", toString(delta_values),
                  "}\n  Cor = {", toString(Cor_values),
                  "}\n  lambda = {", toString(lambda_values),
                  "}\n\n Se simularán ",
                  length(Alpha_values)*length(S_values)*length(n_values)*
                    length(Var_values)*length(q_values)*length(loc_values)*
                    length(delta_values)*length(Cor_values)*length(lambda_values),
                  " escenarios bajo el modelo ", modelSim, ".", "\n"))

# -----------------------------------------

# Simulación
# ----------

if (modelSim != ""){
  
  "m" <- modelSim
  
  # Restringimos los valores a simular:
  Alpha_values_sim <- Alpha_values[1]
  S_values_sim <- S_values[2]
  n_values_sim <- n_values[2]
  Var_values_sim <- Var_values[1]
  q_values_sim <- q_values[1:2] # c(3, 5, 8, 10) ∀ q ∈ "qlocstdev.norm.tsv" [c(2, 3, 5, 8, 10, 12, 15, 20, 25)]
  loc_values_sim <- loc_values[] # ∀ loc ∈ "qlocstdev.norm.tsv" [c(1, 2, 3, 5, 8, 10)]
                                 # loc_values_q_3_sim <- NULL
                                 # loc_values_q_5_sim <- NULL
  loc_values_q_3_sim <- loc_values[1:4]
  loc_values_q_5_sim <- loc_values[1:3]
  delta_values_sim <- delta_values[]
  Cor_values_sim <- Cor_values[]
  lambda_values_sim <- lambda_values[]
  
  NumMaxSim <- length(Alpha_values_sim)*length(S_values_sim)*length(n_values_sim)*
    length(Var_values_sim)*length(q_values_sim)*length(loc_values_sim)*
    length(delta_values_sim)*length(Cor_values_sim)*length(lambda_values_sim)
  
  # Cuando no se simula el mismo número de loc para cada q:
  if(is.null(loc_values_q_3_sim) & is.null(loc_values_q_5_sim)){
    NumMaxSim <- NumMaxSim
    cat("\014") # Clean console
    writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values_sim),
                      "}\n  S = {", toString(S_values_sim),
                      "}\n  n = {", toString(n_values_sim),
                      "}\n  Var = {", toString(Var_values_sim),
                      "}\n  q = {", toString(q_values_sim),
                      "}\n  loc = {", toString(loc_values_sim),
                      "}\n  delta = {", toString(delta_values_sim),
                      "}\n  Cor = {", toString(Cor_values_sim),
                      "}\n  lambda = {", toString(lambda_values_sim),
                      "}\n\n Se simularán ", NumMaxSim,
                      " escenarios bajo el modelo ", modelSim, ".", "\n"))
  } else{
    NumMaxSim <- NumMaxSim - (abs(length(loc_values_q_3_sim) -
                                    length(loc_values_q_5_sim)))*
      length(delta_values_sim)
    cat("\014") # Clean console
    writeLines(paste0("\nCon la combinación de valores actual:\n\n  alpha = {", toString(Alpha_values_sim),
                      "}\n  S = {", toString(S_values_sim),
                      "}\n  n = {", toString(n_values_sim),
                      "}\n  Var = {", toString(Var_values_sim),
                      "}\n  Para q = {", toString(q_values_sim[1]),
                      "} => loc = {", toString(loc_values_q_3_sim),
                      "}\n  Para q = {", toString(q_values_sim[2]),
                      "} => loc = {", toString(loc_values_q_5_sim),
                      "}\n  delta = {", toString(delta_values_sim),
                      "}\n  Cor = {", toString(Cor_values_sim),
                      "}\n  lambda = {", toString(lambda_values_sim),
                      "}\n\n Se simularán ", NumMaxSim,
                      " escenarios bajo el modelo ", modelSim, ".", "\n"))
  }
  
  # Se crea el directorio de la simulación según convenga:
  
  results_dir <- "Resultados"
  model_dir <- paste("modelSim_", modelSim)
  
  Simul_count <- 0
  t_S_0 <- Sys.time()
  
  for(Alpha in Alpha_values_sim){
    for(S in S_values_sim){
      for(n in n_values_sim){
        for(Var in Var_values_sim){
          "v" <- Var
          for(q in q_values_sim){
            # Variable dependiente de q:
            "mu" <- rep(10, q)
            if(q == 3){
              loc_values_sim <- loc_values[1:4]
            } else if(q == 5){
              loc_values_sim <- loc_values[1:3]
            }
            for(loc in loc_values_sim){
              "position" <- p <- loc
              for(delta in delta_values_sim){
                "d" <- delta
                for(Cor in Cor_values_sim){
                  "c" <- Cor
                  for(lambda in lambda_values_sim){
                    "Lambda" <- lambda <- l
                    
                    Simul_count <- Simul_count + 1
                    
                    if (modelSim == "mvnorm"){
                      
                      CompMantaManova_mvnorm()
                      
                    } else if (modelSim == "simplex") {
                      
                      CompMantaManova_simplex()
                      
                    } else if (modelSim == "multinom") {
                      
                      CompMantaManova_multinom()
                      
                    }
                    
                    DF_Results <- rbind(DF_Results, DF_CompPot_res)
                    DF_Results_byDatos <- DF_Results[order(DF_Results$Datos), ]
                    
                    # # # # # # # # # # # # # # # # # # # # # # # # # 
                    #   Almacenaje de los datos de cada simulación  #
                    # # # # # # # # # # # # # # # # # # # # # # # # #
                    
                    sim_path = file.path(getwd(), paste0("Resultados/Modelo ", modelSim),
                                         paste0("Sim. ", format(t_S_0, '%d-%m-%Y %H h %M min')))
                    dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path", sim_path, envir = .GlobalEnv)
                    
                    sim_path_DFsSim = file.path(sim_path, paste0("DFs Simulaciones"))
                    dir.create(sim_path_DFsSim, recursive = TRUE, showWarnings = FALSE)
                    assign("sim_path_DFsSim", sim_path_DFsSim, envir = .GlobalEnv)
                    
                    # Guardamos la información necesaria hasta este punto:
                    # => deparse(substitute(df)) solventa un problema a la hora de guardar el DF como CSV.
                    #    (https://stackoverflow.com/questions/37998967)
                    
                    # Guardamos en formato "csv" el DF con los resultados de la combinación actual de variables:
                    write.csv(DF_CompPot_res, file = file.path(file.path(sim_path_DFsSim),
                                                               paste0(deparse(substitute(DF_CompPot_res)),
                                                                      " [m = ", toString(modelSim),
                                                                      ", alpha = ", toString(Alpha),
                                                                      ", S = ", toString(S),
                                                                      ", n = ", toString(n),
                                                                      ", Var = ", toString(Var),
                                                                      ", q = ", toString(q),
                                                                      ", loc = ", toString(loc),
                                                                      ", stdev = ", toString(stdev),
                                                                      ", delta = ", toString(sprintf("%1.3f", delta)),
                                                                      ", Cor = ", toString(Cor),
                                                                      ", lambda = ", toString(lambda),
                                                                      "]", ".csv")), row.names = FALSE)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  MANTA_col <- c(rep("MANTA", dim(DF_Results)[1]))
  MANOVA_col <- c(rep("MANOVA", dim(DF_Results)[1]))
  DF_Results <- add_column(DF_Results, MANTA_col, .after = 12)
  DF_Results <- add_column(DF_Results, MANOVA_col, .after = 15)
  DF_Results_MANTA <- DF_Results[1:15]
  DF_Results_MANOVA <- DF_Results[c(1:12, 16:18)]
  colnames(DF_Results_MANOVA) <- colnames(DF_Results_MANTA)
  DF_Results_Graph <- rbind(DF_Results_MANTA, DF_Results_MANOVA)
  colnames(DF_Results_Graph) <- c("Datos", "Modelo", "alpha", "S", "n", "Var", "q", "loc", "stdev",
                                  "Cor", "delta", "lambda", "Método", "Potencia", "tcomp")
  
  # Guardamos en formato "csv" los DF con los resultados finales:
  
  write.csv(DF_Results, file = file.path(file.path(sim_path),
                                         paste0(deparse(substitute(DF_Results)), ".csv"))
            , row.names = FALSE)
  
  write.csv(DF_Results_Graph, file = file.path(file.path(sim_path),
                                               paste0(deparse(substitute(DF_Results_Graph)), ".csv"))
            , row.names = FALSE)
  
}

# ----------

# Resultados
# ----------

# # # # # # # # # # # # # # # # # # # # # # #
# Simulación: "Sim. 02-01-2024 01 h 56 min" #
# # # # # # # # # # # # # # # # # # # # # # #
#              ∆ max = 0.025                #
# # # # # # # # # # # # # # # # # # # # # # #

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

## Directorio de la simulación a representar:
setwd(file.path(getwd(), "Resultados/Modelo simplex/Sim. 02-01-2024 01 h 56 min"))
DF_Results_2graph <- read.csv("DF_Results.csv")
DF_Results_Graph_2graph <- read.csv("DF_Results_Graph.csv")

## Para simplex forzamos que todos los valores de Cor y lambda son NA.
DF_Results_2graph$Cor <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$Cor <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))

## Se estudiará MANTA, teniendo en cuenta los datos sin transformar y las transformaciones log-ratio, sqrt y clr:
DF_Results_2graph <- DF_Results_2graph[DF_Results_2graph$Datos == unique(DF_Results_2graph$Datos)[1:4], ]
DF_Results_Graph_2graph <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos ==
                                                     unique(DF_Results_Graph_2graph$Datos)[1:4], ]
DF_Results_Graph_2graph <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Método ==
                                                     unique(DF_Results_Graph_2graph$Método[1]), ]

# (a): ¿Es MANTA invariante a las transformaciones bajo una distribución de datos simulada
#       mediante el modelo simplex?

## => Grid "q - loc" del estudio de la posible invarianza a la transformación de datos de MANTA:
Count_graph = 0

for (S_graph in unique(DF_Results_Graph_2graph$S)){
  for (alpha_graph in unique(DF_Results_Graph_2graph$alpha)){
    
    Count_graph = Count_graph + 1
    
    DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$S == S_graph
                                       & DF_Results_Graph_2graph$alpha == alpha_graph, ]
    
    # PDF
    graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/Estudio Invarianza MANTA"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    pdf(file = file.path(graph_path,
                         paste0("Grid q-loc [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].pdf")),
        width = 16, height = 12)
    
    p <- DeltaPotScatPlot_Method_simplex_1(DF_temp)
    print(p)
    
    dev.off()
    
    # PNG
    graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/Estudio Invarianza MANTA"))
    dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
    
    ppi = 600
    png(file = file.path(graph_path,
                         paste0("Grid q-loc [S = ", toString(S_graph),
                                ", alpha = ", toString(alpha_graph), "].png")),
        width = 16*ppi, height = 12*ppi, res = ppi)
    
    p <- DeltaPotScatPlot_Method_simplex_1(DF_temp)
    print(p)
    
    dev.off()
    
    # ZOOMED
    DeltaMin_vect <- c(0.000, 0.010)
    PotMin_vect <- c(0.05, 0.65)
    DeltaMax_vect <- c(0.010, 0.025)
    PotMax_vect <- c(0.65, 1)
    
    for (i in 1:length(DeltaMin_vect)){
      
      DeltaMin = DeltaMin_vect[i]
      PotMin = PotMin_vect[i]
      DeltaMax = DeltaMax_vect[i]
      PotMax = PotMax_vect[i]
      
      # PDF
      graph_path = file.path(paste0(getwd(), "/Gráficas/PDF/Estudio Invarianza MANTA"))
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      pdf(file = file.path(graph_path,
                           paste0("Grid q-loc [Zoom ∆ = (", sprintf("%1.3f", DeltaMin),
                                  ", ", sprintf("%1.3f", DeltaMax),
                                  "), Zoom P = (", sprintf("%1.3f", PotMin),
                                  ", ", sprintf("%1.3f", PotMax),
                                  "), S = ", toString(S_graph),
                                  ", alpha = ", toString(alpha_graph), "].pdf")),
          width = 16, height = 12)
      
      p <- DeltaPotScatPlot_Method_simplex_1_ZoomDelta(DF_temp, DeltaMin, DeltaMax, PotMin, PotMax)
      print(p)
      
      dev.off()
      
      # PNG
      graph_path = file.path(paste0(getwd(), "/Gráficas/PNG/Estudio Invarianza MANTA"))
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      ppi = 600
      png(file = file.path(graph_path,
                           paste0("Grid q-loc [Zoom ∆ = (", sprintf("%1.3f", DeltaMin),
                                  ", ", sprintf("%1.3f", DeltaMax),
                                  "| Zoom P = (", sprintf("%1.3f", PotMin),
                                  ", ", sprintf("%1.3f", PotMax),
                                  "), S = ", toString(S_graph),
                                  ", alpha = ", toString(alpha_graph), "].png")),
          width = 16*ppi, height = 12*ppi, res = ppi)
      
      p <- DeltaPotScatPlot_Method_simplex_1_ZoomDelta(DF_temp, DeltaMin, DeltaMax, PotMin, PotMax)
      print(p)
      
      dev.off()
      
    }
  }
}

cat("\014") # Clean console

# ----------






























################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# # # # # # # # # # # # # # # # # # # # # # #
# Simulación: "Sim. 21-12-2023 18 h 22 min" #
# # # # # # # # # # # # # # # # # # # # # # #

# Directorio de trabajo principal.
work_dir <- "/Users/aitor/Desktop/RProjects/TFMgit"
setwd(work_dir)

# Directorio de la simulación a representar:
setwd(file.path(getwd(), "Resultados/Modelo mvnorm/Sim. 21-12-2023 18 h 22 min"))
DF_Results_2graph <- read.csv("DF_Results.csv")
DF_Results_Graph_2graph <- read.csv("DF_Results_Graph.csv")

## Para mvnorm forzamos que todos los valores de loc y lambda son NA.
DF_Results_2graph$loc <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$loc <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))
DF_Results_Graph_2graph$lambda <- c(rep(NA, dim(DF_Results_2graph)[1]))

# (a) + (b): Gráficas "∆-Potencia" comparativas entre MANTA-MANOVA (Datos sin transformar + Datos transformados).

## => Gráficas de todas las combinaciones posibles de "Tipo de datos - S - α - Var - Cor":
Count_graph = 0

for (Datos_graph in unique(DF_Results_2graph$Datos)){
  for (S_graph in unique(DF_Results_2graph$S)){
    for (alpha_graph in unique(DF_Results_2graph$alpha)){
      for (Var_graph in unique(DF_Results_2graph$Var)){
        for (Cor_graph in unique(DF_Results_2graph$Cor)){
          
          Count_graph = Count_graph + 1
          
          DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos == Datos_graph
                                             & DF_Results_Graph_2graph$S == S_graph
                                             & DF_Results_Graph_2graph$alpha == alpha_graph
                                             & DF_Results_Graph_2graph$Var == Var_graph
                                             & DF_Results_Graph_2graph$Cor == Cor_graph, ]
          
          # PDF
          graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Gráficas/PDF"),
                                 Datos_graph, 
                                 paste0("Var = ", Var_graph), 
                                 paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
          dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
          
          pdf(file = file.path(graph_path,
                               paste0("Gráfica [Datos = (", toString(Datos_graph),
                                      "), S = (", toString(S_graph),
                                      "), alpha = (", toString(alpha_graph),
                                      "), Var = (", toString(Var_graph),
                                      "), Cor = (", sprintf("%1.1f", Cor_graph), ")]", ".pdf")),
              width = 16, height = 12)
          
          p <- DeltaPotGraph_CorColor(DF_temp)
          print(p)
          
          dev.off()
          
          # PNG
          graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Gráficas/PNG"),
                                 Datos_graph, 
                                 paste0("Var = ", Var_graph), 
                                 paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
          dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
          
          ppi = 600
          png(file = file.path(graph_path,
                               paste0("Gráfica [Datos = (", toString(Datos_graph),
                                      "), S = (", toString(S_graph),
                                      "), alpha = (", toString(alpha_graph),
                                      "), Var = (", toString(Var_graph),
                                      "), Cor = (", sprintf("%1.1f", Cor_graph), ")]", ".png")),
              width = 16*ppi, height = 12*ppi, res = ppi)
          
          p <- DeltaPotGraph_CorColor(DF_temp)
          print(p)
          
          dev.off()
          
        }
      }
    }
  }
}

cat("\014") # Clean console

# (c): ¿Es MANOVA invariante a las transformaciones, y MANTA?

## => Gráficas de todas las combinaciones posibles de "S - α - Var - Cor" para el estudio
##    de la posible invarianza a la transformación de datos de cada método:
Count_graph = 0
Metodo_graph_list = c("MANTA", "MANOVA")

for(Metodo_graph in Metodo_graph_list){
  for (S_graph in unique(DF_Results_2graph$S)){
    for (alpha_graph in unique(DF_Results_2graph$alpha)){
      for (Var_graph in unique(DF_Results_2graph$Var)){
        for (Cor_graph in unique(DF_Results_2graph$Cor)){
          
          Count_graph = Count_graph + 1
          
          DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Método == Metodo_graph
                                             & DF_Results_Graph_2graph$S == S_graph
                                             & DF_Results_Graph_2graph$alpha == alpha_graph
                                             & DF_Results_Graph_2graph$Var == Var_graph
                                             & DF_Results_Graph_2graph$Cor == Cor_graph, ]
          
          # PDF
          graph_path = file.path(paste0(getwd(), "/Estudio Invarianza ", Metodo_graph, "/Gráficas/PDF"),
                                 paste0("Var = ", Var_graph), 
                                 paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
          dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
          
          pdf(file = file.path(graph_path,
                               paste0("Gráfica [S = (", toString(S_graph),
                                      "), alpha = (", toString(alpha_graph),
                                      "), Var = (", toString(Var_graph),
                                      "), Cor = (", sprintf("%1.1f", Cor_graph), ")]", ".pdf")),
              width = 16, height = 12)
          
          p <- DeltaPotGraph_DataTypeColor_mvnorm(DF_temp, Metodo_graph)
          print(p)
          
          dev.off()
          
          # PNG
          graph_path = file.path(paste0(getwd(), "/Estudio Invarianza ", Metodo_graph, "/Gráficas/PNG"),
                                 paste0("Var = ", Var_graph), 
                                 paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
          dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
          
          ppi = 600
          png(file = file.path(graph_path,
                               paste0("S = (", toString(S_graph),
                                      "), alpha = (", toString(alpha_graph),
                                      "), Var = (", toString(Var_graph),
                                      "), Cor = (", sprintf("%1.1f", Cor_graph), ")]", ".png")),
              width = 16*ppi, height = 12*ppi, res = ppi)
          
          p <- DeltaPotGraph_DataTypeColor_mvnorm(DF_temp, Metodo_graph)
          print(p)
          
          dev.off()
          
        }
      }
    }
  }
}

cat("\014") # Clean console

# (d): Estudio de la similitud de resultados entre MANTA y MANOVA (∆Potencia = P_MANTA - P_MANOVA).

## Comparativa mediante boxplot:

### Opción 1
Count_graph = 0

for (S_graph in unique(DF_Results_2graph$S)){
  for (alpha_graph in unique(DF_Results_2graph$alpha)){
    for (Var_graph in unique(DF_Results_2graph$Var)){
      for (Cor_graph in unique(DF_Results_2graph$Cor)){
        for (Datos_graph in unique(DF_Results_2graph$Datos)){
          
          Count_graph = Count_graph + 1
          
          DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos == Datos_graph
                                             & DF_Results_Graph_2graph$S == S_graph
                                             & DF_Results_Graph_2graph$alpha == alpha_graph
                                             & DF_Results_Graph_2graph$Var == Var_graph
                                             & DF_Results_Graph_2graph$Cor == Cor_graph, ]
          
          # PDF
          graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (Opción 1)/PDF"),
                                 paste0("Cor = ", Cor_graph), 
                                 paste0("Var = ", Var_graph),
                                 paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
          dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
          
          pdf(file = file.path(graph_path,
                               paste0("Boxplot [Datos = (", toString(Datos_graph),
                                      "), S = (", toString(S_graph),
                                      "), alpha = (", toString(alpha_graph),
                                      "), Var = (", toString(Var_graph),
                                      "), Cor = (", sprintf("%1.1f", Cor_graph),
                                      ")]", ".pdf")),
              width = 16, height = 12)
          
          p <- DeltaPotBoxPlot_Method_mvnorm_1(DF_temp)
          print(p)
          
          dev.off()
          
          # Guardamos en formato "csv" el DF utilizado para la gráfica actual:
          write.csv(DF_graph, file = file.path(graph_path, paste0(deparse(substitute(DF_graph)),
                                                                  Count_graph, ".csv")), row.names = FALSE)
          
          # # PNG
          # graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (Opción 1)/PNG"),
          #                        paste0("Cor = ", Cor_graph), 
          #                        paste0("Var = ", Var_graph),
          #                        paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
          # dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
          # 
          # ppi = 600
          # png(file = file.path(graph_path,
          #                      paste0("Boxplot [Datos = (", toString(Datos_graph),
          #                             "), S = (", toString(S_graph),
          #                             "), alpha = (", toString(alpha_graph),
          #                             "), Var = (", toString(Var_graph),
          #                             "), Cor = (", sprintf("%1.1f", Cor_graph),
          #                             ")]", ".png")),
          #     width = 16*ppi, height = 12*ppi, res = ppi)
          # 
          # p <- DeltaPotBoxPlot_Method_mvnorm_1(DF_temp)
          # print(p)
          # 
          # dev.off()
        }
      }
    }
  }
}

cat("\014") # Clean console

### Opción 2
Count_graph = 0

for (S_graph in unique(DF_Results_2graph$S)){
  for (alpha_graph in unique(DF_Results_2graph$alpha)){
    for (Var_graph in unique(DF_Results_2graph$Var)){
      for (Cor_graph in unique(DF_Results_2graph$Cor)){
        
        Count_graph = Count_graph + 1
        
        DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos == Datos_graph
                                           & DF_Results_Graph_2graph$S == S_graph
                                           & DF_Results_Graph_2graph$alpha == alpha_graph
                                           & DF_Results_Graph_2graph$Var == Var_graph
                                           & DF_Results_Graph_2graph$Cor == Cor_graph, ]
        
        # PDF
        graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (Opción 2)/PDF"),
                               paste0("Cor = ", Cor_graph), 
                               paste0("Var = ", Var_graph),
                               paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
        dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
        
        pdf(file = file.path(graph_path,
                             paste0("Boxplot [Datos = (", toString(Datos_graph),
                                    "), S = (", toString(S_graph),
                                    "), alpha = (", toString(alpha_graph),
                                    "), Var = (", toString(Var_graph),
                                    "), Cor = (", sprintf("%1.1f", Cor_graph),
                                    ")]", ".pdf")),
            width = 16, height = 12)
        
        p <- DeltaPotBoxPlot_Method_mvnorm_2(DF_temp)
        print(p)
        
        dev.off()
        
        # Guardamos en formato "csv" el DF utilizado para la gráfica actual:
        write.csv(DF_graph, file = file.path(graph_path, paste0(deparse(substitute(DF_graph)),
                                                                Count_graph, ".csv")), row.names = FALSE)
        
        # # PNG
        # graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (Opción 2)/PNG"),
        #                        paste0("Cor = ", Cor_graph), 
        #                        paste0("Var = ", Var_graph),
        #                        paste0("alpha = ", sprintf("%1.3f", alpha_graph)))
        # dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
        # 
        # ppi = 600
        # png(file = file.path(graph_path,
        #                      paste0("Boxplot [Datos = (", toString(Datos_graph),
        #                             "), S = (", toString(S_graph),
        #                             "), alpha = (", toString(alpha_graph),
        #                             "), Var = (", toString(Var_graph),
        #                             "), Cor = (", sprintf("%1.1f", Cor_graph),
        #                             ")]", ".png")),
        #     width = 16*ppi, height = 12*ppi, res = ppi)
        # 
        # p <- DeltaPotBoxPlot_Method_mvnorm_2(DF_temp)
        # print(p)
        # 
        # dev.off()
        
      }
    }
  }
}

cat("\014") # Clean console

### Opción 3: boxplots "∆ vs. Potencia" .𝜏. "facet_grid(alpha ~ Datos)" para "S - Var - Cor" fijos.
Count_graph = 0

for (S_graph in unique(DF_Results_2graph$S)){
  for (Var_graph in unique(DF_Results_2graph$Var)){
    for (Cor_graph in unique(DF_Results_2graph$Cor)){
      
      Count_graph = Count_graph + 1
      
      DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$S == S_graph
                                         & DF_Results_Graph_2graph$Var == Var_graph
                                         & DF_Results_Graph_2graph$Cor == Cor_graph, ]
      
      # PDF
      graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (alpha - Datos - Método)"),
                             paste0("S = ", S_graph), Var_graph, "PDF")
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      pdf(file = file.path(graph_path,
                           paste0("Boxplot (alpha - Datos - Método) [S = (", toString(S_graph),
                                  "), Var = (", toString(Var_graph),
                                  "), Cor = (", sprintf("%1.1f", Cor_graph),
                                  ")]", ".pdf")),
          width = 16, height = 12)
      
      ## Boxplots comparativos "∆ vs. Potencia" .𝜏. "facet_grid(alpha ~ Datos)" para "S - Var - Cor"
      ## fijo, y agrupados por "método".
      p <- DeltaPotBoxPlot_Method_mvnorm_2(DF_temp)
      print(p)
      
      dev.off()
      
      # PNG
      graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (alpha - Datos - Método)"),
                             paste0("S = ", S_graph), Var_graph, "PNG")
      dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
      
      ppi = 600
      png(file = file.path(graph_path,
                           paste0("Boxplot (alpha - Datos - Método) [S = (", toString(S_graph),
                                  "), Var = (", toString(Var_graph),
                                  "), Cor = (", sprintf("%1.1f", Cor_graph),
                                  ")]", ".png")),
          width = 16*ppi, height = 12*ppi, res = ppi)
      
      ## Boxplots comparativos "∆ vs. Potencia" .𝜏. "facet_grid(alpha ~ Datos)" para "S - Var - Cor"
      ## fijo, y agrupados por "método".
      p <- DeltaPotBoxPlot_Method_mvnorm_2(DF_temp)
      print(p)
      
      dev.off()
    }
  }
}

cat("\014") # Clean console

### Opción 4: boxplots "∆ vs. Potencia" .𝜏. "facet_wrap(~ Cor)" para "Tipo de datos - S -
### α - Var" fijos.
Count_graph = 0

for (Datos_graph in unique(DF_Results_2graph$Datos)){
  for (S_graph in unique(DF_Results_2graph$S)){
    for (alpha_graph in unique(DF_Results_2graph$alpha)){
      for (Var_graph in unique(DF_Results_2graph$Var)){
        
        Count_graph = Count_graph + 1
        
        DF_temp <- DF_Results_Graph_2graph[DF_Results_Graph_2graph$Datos == Datos_graph
                                           & DF_Results_Graph_2graph$S == S_graph
                                           & DF_Results_Graph_2graph$alpha == alpha_graph
                                           & DF_Results_Graph_2graph$Var == Var_graph, ]
        
        # PDF
        graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (Cor - Método)"),
                               paste0("S = ", S_graph), Datos_graph, Var_graph, "PDF")
        dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
        
        pdf(file = file.path(graph_path,
                             paste0("Boxplot (Cor - Método) [Datos = (", Datos_graph,
                                    "), S = (", toString(S_graph),
                                    "), alpha = (", sprintf("%1.3f", alpha_graph),
                                    "), Var = (", toString(Var_graph),
                                    ")]", ".pdf")),
            width = 16, height = 12)
        
        ## Boxplots comparativos "∆ vs. Potencia" .𝜏. "facet_grid(alpha ~ Datos)" para "S - Var - Cor"
        ## fijo, y agrupados por "método".
        p <- DeltaPotBoxPlot_Method_mvnorm_1(DF_temp)
        print(p)
        
        dev.off()
        
        # PNG
        graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Boxplot (Cor - Método)"),
                               paste0("S = ", S_graph), Datos_graph, Var_graph, "PNG")
        dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
        
        ppi = 600
        png(file = file.path(graph_path,
                             paste0("Boxplot (Cor - Método) [Datos = (", Datos_graph,
                                    "), S = (", toString(S_graph),
                                    "), alpha = (", sprintf("%1.3f", alpha_graph),
                                    "), Var = (", toString(Var_graph),
                                    ")]", ".png")),
            width = 16*ppi, height = 12*ppi, res = ppi)
        
        ## Boxplots comparativos "∆ vs. Potencia" .𝜏. "facet_grid(alpha ~ Datos)" para "S - Var - Cor"
        ## fijo, y agrupados por "método".
        p <- DeltaPotBoxPlot_Method_mvnorm_1(DF_temp)
        print(p)
        
        dev.off()
        
      }
    }
  }
}

cat("\014") # Clean console

## Mediante la "Diferencia en valor absoluto" entre ambas medidas para todas las simulaciones:
DF_AbsDiffPotMANTAMANOVA <- data.frame(DF_Results_2graph[1:12], abs(DF_Results_2graph[14] - DF_Results_2graph[17]))
DF_AbsDiffPotMANTAMANOVA <- DF_AbsDiffPotMANTAMANOVA[DF_AbsDiffPotMANTAMANOVA$alpha == 0.05 &
                                                       DF_AbsDiffPotMANTAMANOVA$Var == "equal", ]
colnames(DF_AbsDiffPotMANTAMANOVA) <- c(colnames(DF_Results_2graph[1:12]), "∆P MANTA-MANOVA")

for (Cor in unique(DF_AbsDiffPotMANTAMANOVA$Cor)){
  
  # PDF
  graph_path = file.path(paste0(getwd(), "/Comparación MANTA-MANOVA/Gráficas/PDF"))
  dir.create(graph_path, recursive = TRUE, showWarnings = FALSE)
  
  pdf(file = file.path(graph_path,
                       paste0("Gráfica ∆P MANTA-MANOVA [Cor = ", sprintf("%1.1f", Cor), "]", ".pdf")),
      width = 16, height = 12)
  
  p <- ggplot(DF_AbsDiffPotMANTAMANOVA[DF_AbsDiffPotMANTAMANOVA$Cor == Cor, ]
              , aes(x = delta, y = `∆P MANTA-MANOVA`, color = Datos)) +
    geom_point()
  
  print(p)
  
  dev.off()
}

cat("\014") # Clean console

## Mediante la "Diferencia de porcentaje relativa" (RPD) entre ambas medidas para todas las simulaciones:


cat("\014") # Clean console

## Multinomial distribution
# ------------------------

# => Se usará la función de simulación de datos "Sim.multinom(ch, q, n, Lambda, delta, loc)".
# => Se simularán todas las combinaciones Δ - Loc (¿también Lambda, Cor?).
# => Donde Δ = 0 evalua H0 y Δ > 0 evaluará H1.
# => La opción dejar fijo Δ pero variar n eleva demasiado el tiempo de computación. <= HACER LA PRUEBA MIDIÉNDOLO!!!
# => Donde Loc indica la ubicación del modelo generador "simplex" (p / position / loc).
# => (loc = 0) ≡ "Centro del simplex" | (loc = 1) ≡ "Vértice"

# Listados de las variables necesarias y sus equivalencias:
list(# FIJAS: a priori, no deberían influir en los resultados.
  "modelSim" <- "multinom", # m_values <- c("mvnorm", "simplex", "multinom")
  "m" <- "multinom", # Equivalencia
  "S" <- 1E3, # Máx. 1E4 si se quiere aumentar la precisión sin elevar el tiempo de computación.
  "S_values" <- c(seq(1E2, 1E4, 1E2)),
  "n" <- 300, # Valor que equilibra la representatividad de los resultados con el tiempo de computación.
  "n_values" <- c(seq(200, 400, 50)),
  "q" <- 3, # q_values <- c(1:4)
  "loc" <- 1, # Ubicación del modelo generador "simplex" (p / position / loc)
  # => (loc = 0) ≡ "Centro del simplex" | (loc = 1) ≡ "Vértice"
  "a" <- 2, # a_values <- c(1:6)
  "b" <- 3, # b_values <- c(1:6)
  "u" <- 1, # u_values <- c(seq(0.1, 2, 0.1))
  "w" <- "B", # c("A", "B", "AB")
  "hk" <- 1,
  "Var" <- "equal", # c("equal", "unequal")
  "k" <- 0, # c(0, 1)
  "plotheatmap" <- F, # as.logical(c("F", "T"))
  "D" <- "unif-0-1",
  "l" <- 1000, # Parámetro lambda (distrib. Poisson): indica el tamaño del modelo generador 'multinom'
  "s" <- 0.1,
  "p" <- 1,
  "pdist" <- "norm", # c("norm", "gamma", "beta")
  "p_dist" <- "norm", # Equivalencia
  "x" <- 10,
  # VARIABLES: a estudiar su potencial influencia en los resultados.
  "S_values" <- c(NA, 1E2, 1E3, 1E4), # Número de simulaciones, aumentar S debería aumentar la precisión
  "n_values" <- c(NA, seq(200, 400, 100)), # Tamaño de la muestra
  "q_values" <- c(NA, seq(3, 10, 1)), # Número de respuestas .𝜏. para q > 5 "SITIOS EXTRAÑOS"!
  "Var_values" <- c(NA, "equal", "unequal"), # Varianza de las variables de Y
  "Alpha_values" <- c(NA, 0.05, 0.01, 0.001), # Niveles de significación (α) más comunes
  "delta_values" <- c(NA, seq(0, 1, 0.1)), # delta (∆) => Δ = 0 evalua H0 y Δ > 0 evaluará H1
  # POSIBLES PROBLEMAS CON LOS VALORES DE Δ para MANOVA en "simplex" y "multinom"!!! <= PQ???
  "Cor_values" <- c(NA, seq(0, 0.8, 0.2)), # Correlación de las variables del conjunto de datos Y (Cor) => Siempre: 0 <= Cor < 1
  "loc_values" <- c(NA, seq(0.0, 1.0, 0.1)), # Ubicación del modelo generador "simplex" (p / position / loc)
  "lambda_values" <- c(NA, seq(200, 1000, 200)) # Parámetro lambda (distrib. Poisson)
  # => (loc = 0) ≡ "Centro del simplex" | (loc = 1) ≡ "Vértice"
)


# Variables dependientes.

# f(q):
"mu" <- rep(10, q)

# f(loc, q):
# Para evitar que el valor de Δ haga que H1 esté fuera del "simplex",
# teniendo en cuenta que [0 ≤ Δ ≤ 1], se acotan sus valores posibles: 
deltamax = 1 - c(c(loc, rep(1, q-1))/sum(c(loc, rep(1, q-1))))[[1]]
delta_values <- c(seq(0.0, 1.0, 0.1))
delta_values <- c(delta_values[delta_values < deltamax], deltamax)

# ------------------------

## Simulación (original)
# ---------------------

# DFs que almacenarán todos los resultados:
DF_OBJ1 = data.frame()
DF_OBJ1MANTA = data.frame()
DF_OBJ1MANOVA = data.frame()


# Lista de almacenaje de gráficas.
graph_list = list()

cat("\014") # Clean console


if (modelSim != ""){
  
  cat("\014") # Clean console
  
  writeLines(paste("\nCon la combinación de valores actuales [alpha = (", toString(Alpha_values),
                   "), delta = (", toString(delta_values), "), Cor = (",
                   toString(Cor_values), ")] se simularán ", length(Alpha_values)*length(delta_values)*length(Cor_values),
                   "escenarios bajo el modelo ", modelSim, ".", "\n"))
  
  for(Alpha in Alpha_values){
    
    # Se crea el directorio de la simulación vigente según el valor de α considerado:
    
    work_dir <- getwd()
    results_dir <- "Resultados"
    Alpha_dir <- paste("Nivel Significación (alpha = ", Alpha, ")")
    
    for(delta in delta_values){
      
      for(Cor in Cor_values){
        
        list(# Equivalencias de variables:
          "m" <- modelSim,
          "c" <- Cor,
          "v" <- Var,
          "chunk" <- k,
          "DistDef" <- dd <- D,
          "Lambda" <- l,
          "stdev" <- s,
          "position" <- loc <- p,
          "p_dist" <- pdist,
          "heterosk" <- H <- hk,
          "d" <- delta,
          "cores" <- x,
          "fx" <- f)
        
        if (modelSim == "mvnorm"){
          
          CompMantaManova_mvnorm()
          
          } else if (modelSim == "simplex") {
            
            CompMantaManova_simplex()
            
            } else if (modelSim == "multinom") {
              
              CompMantaManova_multinom()
              
              }
        
        DF_OBJ1 <- rbind(DF_OBJ1, DF_CompPot_res)
        DF_OBJ1 <- DF_OBJ1[order(DF_OBJ1$`Tipo de datos`),]
        DF_OBJ1MANTA <- rbind(DF_OBJ1MANTA, DF_MANTA_Results)
        DF_OBJ1MANOVA <- rbind(DF_OBJ1MANOVA, DF_MANOVA_Results)
        
        # # # # # # # # # # # # # # # # # # # # # #
        #  Almacenaje de los datos de simulación  #
        # # # # # # # # # # # # # # # # # # # # # #
        
        if (modelSim == "simplex") {
          sim_dir_name = paste("[model = ", m, ", S = ",  S, ", n = ",  n, ", a = ",  a, ", b = ",  b,
                               ", q (Q) = ", q, ", loc (L) = ",  loc,", stdev = ", stdev,", pdist = ", pdist,
                               ", u = ",  u, ", w = ",  w, ", delta = ",  delta, ", Cor = ",  Cor, "]", sep = "")
        } else if (modelSim == "mvnorm") {
          sim_dir_name = paste("[model = ", m, ", S = ",  S, ", n = ",  n, ", a = ",  a, ", b = ",  b,
                               ", q = ", q, ", u = ",  u, ", w = ",  w, ", delta = ",  delta, ", Cor = ",  Cor, "]", sep = "")
        } else if (modelSim == "multinom") {
          sim_dir_name =  paste("[model = ", m, ", S = ",  S, ", n = ",  n, ", a = ",  a, ", b = ",  b,
                                ", q = ", q, ", u = ",  u, ", w = ",  w, ", delta = ",  delta, ", Cor = ",  Cor, "]", sep = "")
        }
        
        sim_path = file.path(work_dir, results_dir, Alpha_dir, sim_dir_name)
        dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
        assign("sim_path", sim_path, envir = .GlobalEnv)
        
        # Guardamos la información necesaria hasta este punto:
        # => deparse(substitute(df)) solventa un problema a la hora de guardar el DF como CSV.
        #    (https://stackoverflow.com/questions/37998967)
        
        if (modelSim == "simplex") {
          # Guardamos en formato "csv" el DF con los valores de "stdev" para la combinación actual de
          # α - ∆ - Cor:
          write.csv(stdevDF, file = file.path(sim_path, paste0(deparse(substitute(stdevDF)), ".csv")),  row.names = FALSE)
          
          # Guardamos en formato "csv" el DF con los resultados de la combinación actual de α - ∆ - Cor:
          write.csv(DF_CompPot_res,
                    file = file.path(sim_path, paste0(deparse(substitute(DF_CompPot_res)), " (delta = ", d, " | Cor = ", Cor, ")", ".csv")),
                    row.names = FALSE)
          
        } else {
          write.csv(DF_CompPot_res,
                    file = file.path(sim_path, paste0(deparse(substitute(DF_CompPot_res)), " (delta = ", d, " | Cor = ", Cor, ")", ".csv")),
                    row.names = FALSE)
        }
        
        # # Generamos la gráfica "∆ vs. Potencia" para las diferentes Y consideradas en MANTA:
        # Method_name = "MANTA"
        # Method_col = 5
        # 
        # plot_MANTA = DeltaPotGraph()
        # graph_list = c(graph_list, list(plot_MANTA))
        # 
        # # Generamos la gráfica "∆ vs. Potencia" para las diferentes Y consideradas en MANOVA:
        # Method_name = "MANOVA"
        # Method_col = 7
        # 
        # plot_MANOVA = DeltaPotGraph()
        # graph_list = c(graph_list, list(plot_MANOVA))
        
      }
      
    }
    
    # Guardamos en formato "csv" los DF con los resultados de todas las combinaciones α - ∆ - Cor:
    write.csv(DF_Results,
              file = file.path(file.path(work_dir, results_dir),
                               paste0(deparse(substitute(DF_Results)),
                                      " [alpha = (", toString(Alpha_values),
                                      "), S = (", toString(S_values),
                                      "), q = (", toString(q_values), 
                                      "), Cor = (", toString(Cor_values),
                                      "), delta = (", toString(delta_values),
                                      "), loc = (", toString(loc_values),
                                      "), lambda = (", toString(lambda_values),
                                      ")]", ".csv")),
              row.names = FALSE)
    
    # DFs de los resultados según el "Tipo de datos" usado:
    DF_OBJ1_Y <- DF_OBJ1[DF_OBJ1$`Tipo de datos` == 'Datos sin transformar', 2:8]
    DF_OBJ1_Y_Log <- DF_OBJ1[DF_OBJ1$`Tipo de datos` == 'Transformación logarítmica', 2:8]
    DF_OBJ1_Y_sqrt <- DF_OBJ1[DF_OBJ1$`Tipo de datos` == 'Transformación raíz cuadrada', 2:8]
    DF_OBJ1_Y_EscDesvTip <- DF_OBJ1[DF_OBJ1$`Tipo de datos` == 'Escalado de datos por desviación típica', 2:8]
    DF_OBJ1_Y_MinMax <- DF_OBJ1[DF_OBJ1$`Tipo de datos` == 'Normalización Min- Max', 2:8]
    
    DF_OBJ1_MANTA_TypeY <- cbind(DF_OBJ1_Y[1:4], DF_OBJ1_Y_Log[4], DF_OBJ1_Y_sqrt[4],
                                 DF_OBJ1_Y_EscDesvTip[4], DF_OBJ1_Y_MinMax[4])
    colnames(DF_OBJ1_MANTA_TypeY) <- c("d", "c", "alpha", "Y", "Y_Log", "Y_sqrt", "Y_EscDesvTip", "Y_MinMax")
    DF_OBJ1_MANOVA_TypeY <- cbind(DF_OBJ1_Y[1:3], DF_OBJ1_Y[6], DF_OBJ1_Y_Log[6], DF_OBJ1_Y_sqrt[6],
                                  DF_OBJ1_Y_EscDesvTip[6], DF_OBJ1_Y_MinMax[6])
    
    DF_OBJ1_Y_TYPE <- unique(DF_OBJ1$`Tipo de datos`)
    DF_OBJ1_delta_vals <- unique(DF_OBJ1$`Delta (d)`)
    DF_OBJ1_Cor_vals <- unique(DF_OBJ1$`Correlation (c)`)
    DF_OBJ1_Alpha_vals <- unique(DF_OBJ1$`Nivel de significación (alpha)`)
    DF_OBJ1_Methods <- c("MANTA", "MANOVA")
    
    # # Se imprimen las gráficas resultantes:
    # pdf(file = file.path(file.path(work_dir, results_dir),
    #                      paste0("Prueba.pdf")), width = 16, height = 12)
    # 
    # # Prueba gráfica comparación métodos "∆ - Potencia" con (α - Cor)_fijo para los diferentes Y:
    # Graph_Data1 <- DF_OBJ1_Y[DF_OBJ1_Y$`Correlation (c)` == 0.0,]
    # colnames(Graph_Data1) <- c("d", "c", "alpha", "MANTA", "tMANTA", "MANOVA", "tMANOVA")
    # Graph_Data2 <- DF_OBJ1_Y_Log[DF_OBJ1_Y_Log$`Correlation (c)` == 0.0,]
    # colnames(Graph_Data2) <- c("d", "c", "alpha", "MANTA", "tMANTA", "MANOVA", "tMANOVA")
    # 
    # ggplot(data = Graph_Data1,
    #        mapping = aes(x = Graph_Data1$d, y = Graph_Data1$MANTA, color = DF_OBJ1_Methods[[1]])) +
    #   geom_line(linetype = "dashed") +
    #   geom_line(data = Graph_Data2,
    #             mapping = aes(x = Graph_Data2$d, y = Graph_Data2$MANTA, color = DF_OBJ1_Methods[[2]])) +
    #   labs(title = TeX(r"($\Delta$ vs. Potencia)", italic = TRUE),
    #        subtitle = paste0("Método: ", Method_name, "\n
    #                          model = ", m, ", alpha = ", Alpha, ", S = ",  S,
    #                          ", n = ",  n, ", a = ",  a, ", b = ",  b, ", q = ", q,
    #                          ", u = ",  u, ", w = ",  w, "]                     "),
    #        x = TeX(r"(Potencia del método)", italic = TRUE),
    #        y = TeX(r"(Parámetro de generación de la H1 ($\Delta$))", italic = TRUE)) +
    #   theme_light() +
    #   theme(plot.title = element_text(hjust = 0.5, size = 14),
    #         plot.subtitle = element_text(hjust = 0.5, size = 8),
    #         axis.title.x = element_text(size = 10),
    #         axis.title.y = element_text(size = 10),
    #         legend.title = element_text(size = 14),
    #         legend.text = element_text(size = 6),
    #         legend.position = c(0.15, 0.85),
    #         legend.background = element_rect(fill=alpha('lightgray', 0.75)))
    # 
    # pMANTA <- ggplot(data = Graph_Data1,
    #             mapping = aes(x = Graph_Data1$d)) +
    #   geom_line(aes(y = DF_OBJ1_Y[DF_OBJ1_Y$`Correlation (c)` == 0.0,][,4]),
    #             size = .3, color = "black", linetype = "solid") +
    #   geom_line(aes(y = DF_OBJ1_Y_Log[DF_OBJ1_Y_Log$`Correlation (c)` == 0.0,][,4]),
    #             size = .6, color = "red", linetype = "dashed") +
    #   geom_line(aes(y = DF_OBJ1_Y_sqrt[DF_OBJ1_Y_sqrt$`Correlation (c)` == 0.0,][,4]),
    #             size = .6, color = "blue", linetype = "dashed") +
    #   geom_line(aes(y = DF_OBJ1_Y_EscDesvTip[DF_OBJ1_Y_EscDesvTip$`Correlation (c)` == 0.0,][,4]),
    #             size = .6, color = "green", linetype = "dashed") +
    #   geom_line(aes(y = DF_OBJ1_Y_MinMax[DF_OBJ1_Y_MinMax$`Correlation (c)` == 0.0,][,4]),
    #             size = .6, color = "orange", linetype = "dashed") +
    #   labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
    #        subtitle = paste0("Método: ", DF_OBJ1_Methods[[1]], "\n
    #                          [model = ", m, ", alpha = ", Alpha, ", S = ",  S,
    #                          ", n = ",  n, ", a = ",  a, ", b = ",  b, ", q = ", q,
    #                          ", u = ",  u, ", w = ",  w, "]                     "),
    #        x = TeX(paste0("Potencia del método"), italic = TRUE),
    #        y = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE)) +
    #   theme_light() +
    #   theme(plot.title = element_text(hjust = 0.5, size = 14),
    #         plot.subtitle = element_text(hjust = 0.5, size = 8),
    #         axis.title.x = element_text(size = 10),
    #         axis.title.y = element_text(size = 10),
    #         legend.title = element_text(size = 14),
    #         legend.text = element_text(size = 6),
    #         legend.position = c(0.08, 0.90),
    #         legend.background = element_rect(fill = alpha('lightgray', 0.75)))
    # 
    # pMANOVA <- ggplot(data = Graph_Data1,
    #                  mapping = aes(x = Graph_Data1$d)) +
    #   geom_line(aes(y = DF_OBJ1_Y[DF_OBJ1_Y$`Correlation (c)` == 0.0,][,6]),
    #             size = .3, color = "black", linetype = "solid") +
    #   geom_line(aes(y = DF_OBJ1_Y_Log[DF_OBJ1_Y_Log$`Correlation (c)` == 0.0,][,6]),
    #             size = .6, color = "red", linetype = "dashed") +
    #   geom_line(aes(y = DF_OBJ1_Y_sqrt[DF_OBJ1_Y_sqrt$`Correlation (c)` == 0.0,][,6]),
    #             size = .6, color = "blue", linetype = "dashed") +
    #   geom_line(aes(y = DF_OBJ1_Y_EscDesvTip[DF_OBJ1_Y_EscDesvTip$`Correlation (c)` == 0.0,][,6]),
    #             size = .6, color = "green", linetype = "dashed") +
    #   geom_line(aes(y = DF_OBJ1_Y_MinMax[DF_OBJ1_Y_MinMax$`Correlation (c)` == 0.0,][,6]),
    #             size = .6, color = "orange", linetype = "dashed") +
    #   labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
    #        subtitle = paste0("Método: ", DF_OBJ1_Methods[[2]], "\n
    #                          [model = ", m, ", alpha = ", Alpha, ", S = ",  S,
    #                          ", n = ",  n, ", a = ",  a, ", b = ",  b, ", q = ", q,
    #                          ", u = ",  u, ", w = ",  w, "]                     "),
    #        x = TeX(paste0("Potencia del método"), italic = TRUE),
    #        y = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE)) +
    #   theme_light() +
    #   theme(plot.title = element_text(hjust = 0.5, size = 14),
    #         plot.subtitle = element_text(hjust = 0.5, size = 8),
    #         axis.title.x = element_text(size = 10),
    #         axis.title.y = element_text(size = 10),
    #         legend.title = element_text(size = 14),
    #         legend.text = element_text(size = 6),
    #         legend.position = c(0.08, 0.90),
    #         legend.background = element_rect(fill = alpha('lightgray', 0.75)))
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    # while (!is.null(dev.list()))  dev.off()
    # 
    # cat("\014") # Clean console
    # 
    # 
    # ############
    # # GRÁFICAS #
    # ############
    # 
    # # Guardamos el "mapa de calor" de los diferentes niveles de los factores A, B y A:B generados en
    # # la simulación del conjunto de datos Y:
    # plot_heatmap = pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
    # graph_list = c(graph_list, list(plot_heatmap))
    # 
    # # Se imprimen las gráficas resultantes:
    # pdf(file = file.path(file.path(work_dir, results_dir, Alpha_dir),
    #                      paste0("Gráficas [alpha = (", toString(Alpha_values),
    #                             "), delta = (", toString(delta_values), "), Cor = (",
    #                             toString(Cor_values), ")]", ".pdf")), width = 16, height = 12)
    # for (p in graph_list) {
    #   print(p)
    #   }
    # 
    # dev.off()
    
  }
}

# ---------------------

# UTILIDADES
# ----------

## Cuando hay problemas al guardar las imágenes:
dev.set(dev.next()) # Opción 1: hasta ver en consola "quartz_off_screen; 3".
dev.set(dev.prev()) # Opción 2: hasta ver en consola "quartz_off_screen; 3".
while (!is.null(dev.list()))  dev.off() # Opción 3
while (dev.cur()[[1]] != 1) { # Opción 4: Se asegura el cierre de todos los dispositivos gráficos.
  dev.off()
}

## Prueba rápida:
Alpha <- 0.05
delta <- 0.2
Cor <- 0.2

cat("\014") # Clean console

## Benchmark:
microbenchmark(
  "expression",
  list = NULL, times = 10L, unit = NULL, check = NULL, control = list(), setup = NULL
)

## Sacar columnas de un DF:
# DF <- subset(DF, select = -col#)

# ----------

# Código tras la reunión con Diego (17/11/23)
# -------------------------------------------

## Caso particular interesante <=> Forzar matrices de covarianza específicas:

X = matrix(rnorm(q^2), q, q)
crossprod(X)

#### Código para tema p-valores (3r objetivo) <=> Revisar y ver qué hace:

library(CompQuadForm)
library(car)
library(MASS)
library(MCMCpack)
library(data.table)
library(microbenchmark)
library(ggplot2)
library(survey)

# Comparison methods tail

set.seed(1)
q <- 3
e <- runif(q)
h <- rep(1, q)

d <- function(acc){
  pv <- davies(ss, lambda = e, h = h, acc = acc)
  if(pv$Qq < 0 || pv$Qq > 1){
    return(pv)
  } else {
    return(pv$Qq)
  }
}
e2 <- e[e/sum(e)>1e-3]
df <- c()
for (ss in 1:50){
  print(ss)
  acc <- 1e-14
  D <- d(acc)
  while (length(D) > 1) {
    acc <- acc * 10
    D <- d(acc)
  }
  if(D < 1e-14) D <- 1e-14
  
  I <- imhof(ss, lambda = e, h = h, epsabs  = 1e-14, epsrel = 1e-14)$Qq
  if(I < 1e-14) I <- 1e-14
  
  acc <- 1e-14
  F <- farebrother(ss, lambda = e2, h = h, eps = acc)$Qq
  while (F <= 0) {
    acc <- acc * 10
    F <- farebrother(ss, lambda = e2, h = h, eps = acc)$Qq
    print(acc)
  }
  if(F < 1e-14) F <- 1e-14
  
  St <- pchisqsum(ss, df = h, a = e, lower.tail = FALSE, method = "satterthwaite")
  if(St < 1e-14) St <- 1e-14
  
  Sd <- pchisqsum(ss, df = h, a = e, lower.tail = FALSE, method = "saddlepoint")
  if(Sd < 1e-14) Sd <- 1e-14
  
  L <- liu(ss, lambda = e, h = h)
  if(L < 1e-14) L <- 1e-14
  
  
  df <- rbind(df, data.frame("SS" = ss,
                             "Davies" = D,
                             "Farebrother" = F,
                             "Imhof" = I,
                             "Liu" = L,
                             "Satterthwaite" = St,
                             "Saddlepoint" = Sd))
}

df.melt <- melt(df, id = "SS")
colnames(df.melt) = c("SS", "Method", "p")
df.melt[df.melt$p < 0, "p" ] <- 1
df.melt[df.melt$p == 0, "p" ] <- 1e-15

ggplot(df.melt) +
  geom_point(aes(x = SS , y = -log10(p), col = Method, shape = Method)) +
  scale_shape() +
  ylab(expression(-log[10]*(italic(p)))) +
  xlab("Test statistic") +
  theme_classic(base_size = 25) +
  ylim(0,15)

barplot((df$Farebrother-df$Saddlepoint)/df$Farebrother, 
        names.arg = sprintf("%.2e", df$Farebrother), las = 2)

set.seed(1)
ss <- 1
mb_1 <- microbenchmark(Davies = davies(ss, lambda = e, h = h, acc = 1e-8)$Qq,
                       Farebrother = farebrother(ss, lambda = e, h = h, eps = 1e-14)$Qq,
                       Imhof = imhof(ss, lambda = e, h = h, epsabs  = 1e-14, epsrel = 1e-14)$Qq,
                       Saddle = pchisqsum(ss, a = e, df = h, lower.tail = FALSE, method = "saddlepoint"), 
                       unit = "ms")
mb_1 <- data.frame(Method=mb_1$expr, Time=mb_1$time)

ss <- 20
mb_20 <- microbenchmark(Davies = davies(ss, lambda = e, h = h, acc = 1e-8)$Qq,
                        Farebrother = farebrother(ss, lambda = e, h = h, eps = 1e-14)$Qq,
                        Imhof = imhof(ss, lambda = e, h = h, epsabs  = 1e-14, epsrel = 1e-14)$Qq,
                        Saddle = pchisqsum(ss, a = e, df = h, lower.tail = FALSE, method = "saddlepoint"), 
                        unit = "ms")
mb_20 <- data.frame(Method=mb_20$expr, Time=mb_20$time)

ss <- 50
mb_50 <- microbenchmark(Davies = davies(ss, lambda = e, h = h, acc = 1e-8)$Qq,
                        Farebrother = farebrother(ss, lambda = e, h = h, eps = 1e-14)$Qq,
                        Imhof = imhof(ss, lambda = e, h = h, epsabs  = 1e-14, epsrel = 1e-14)$Qq,
                        Saddle = pchisqsum(ss, a = e, df = h, lower.tail = FALSE, method = "saddlepoint"), 
                        unit = "ms")
mb_50 <- data.frame(Method=mb_50$expr, Time=mb_50$time)

mb_1$Time <- mb_1$Time/1e3

ggplot(mb_1) +
  geom_boxplot(aes(x = Method, y = Time/1e3, fill = Method)) +
  ylab(expression(Time~(ms))) +
  theme_classic(base_size = 20)

ggplot(subset(mb_50, Method != "Imhof")) +
  geom_boxplot(aes(x = Method, y = Time, fill = Method)) +
  ylab(expression(Time~(mu*s))) +
  theme_classic(base_size = 20) +
  theme(legend.position = "none")

# -------------------------------------------