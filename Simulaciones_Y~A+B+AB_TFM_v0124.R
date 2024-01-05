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












# UTILIDADES
# ----------

## Cuando hay problemas al guardar las imágenes:
dev.set(dev.next()) # Opción 1: hasta ver en consola "quartz_off_screen; 3".
dev.set(dev.prev()) # Opción 2: hasta ver en consola "quartz_off_screen; 3".
while (!is.null(dev.list()))  dev.off() # Opción 3
while (dev.cur()[[1]] != 1) { # Opción 4: Se asegura el cierre de todos los dispositivos gráficos.
  dev.off()
}

# ----------

# Benchmark
# ---------

microbenchmark(
  "expression",
  list = NULL, times = 10L, unit = NULL, check = NULL, control = list(), setup = NULL
)

# ---------

# Sacar columnas de un DF
# -----------------------

# DF <- subset(DF, select = -col#)

# -----------------------