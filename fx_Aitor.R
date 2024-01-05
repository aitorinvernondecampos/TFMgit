########################### 
# Funciones de simulación #
###########################

# Función que simula el conjunto de datos Y, lo transforma, y aplica los diferentes
# métodos (MANTA, MANOVA,...) para el modelo "mvnorm":
CompMantaManova_mvnorm <- function() {

  # Instalación y carga de los paquetes necesarios
  # ----------------------------------------------
  
  # Se cargan los paquetes, instalándolos si es necesario.
  
  #    require() => a warning if a package is not installed and then continue to execute the code
  #    library() => an error and stop the execution of the code
  
  if (!require("devtools")) install.packages("devtools")
  require(devtools)
  library(devtools)
  
  if (!require("versions")) install.packages("versions")
  require(versions)
  library(versions)
  
  if (!require("manta")) devtools::install_github("dgarrimar/manta")
  require(manta)
  library(manta)
  
  package_list <- c("CompQuadForm", "car", "Hmisc", "MASS", "copula", "vegan",
                    "randomcoloR", "plyr", "svMisc", "reshape2", "hms", "progress",
                    "ggplot2", "arm", "pheatmap", "latex2exp", "tibble", "compositions",
                    "glue")
  
  # Hmisc => Necesario para la función label
  # randomcoloR => La función randomcoloR() genera un color aleatóriamente
  # plyr => Necesario para la función mapvalues()
  # svMisc => Permite usar progress() para indicar visualmente el progreso de la simulación
  # reshape2 => La función melt() permite convertir un objeto "dist" en "DF"
  # hms => La función as_hms() permite convertir "difftime" en DD:HH:MM:SS
  # progress => Para mostrar una barra de progreso de la simulación
  # latex2exp => Toma una cadena LaTeX, la analiza y devuelve la expresión plotmath más cercana
  # glue => Permite pasar expresiones (símbolos griegos, equaciones, etc.) a las etiquetas de facet.grid
  
  for (i in package_list){
    if(! i %in% installed.packages()){
      install.packages(i)
    }
    require(i, character.only = TRUE)
    library(i, character.only = TRUE)
  }
  
  # ----------------------------------------------
  
  
  # DF de los parámetros del modelo
  # -------------------------------
  
  # Estos son: q, a, b, n, u, m, c, v, w, S, plotheatmap, k, D, l, s, p, p_dist, H, d, t, x, f.
  
  # Dependent variables: the parameter q (numerical type) is equivalent to the number of dependent variables.
  # (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 4)
  
  # Factor A Levels: the parameter a (numerical type) is equivalent to the number of levels of factor A.
  # (Value by default: 2 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
  
  # Factor B Levels: the parameter b (numerical type) is equivalent to the number of levels of factor B."
  # (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
  
  # Number of samples: the parameter n (numerical type) is equivalent to the number of samples.
  # (Value by default: 300 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 400)
  
  # Unbalance 1st B Level: the parameter u (numerical type) is equivalent to the unbalance degree of the first level regarding the rest of the levels of the factor B.
  # (Value by default: 1 / Possible values: all positive rational numbers / Max. recommended to avoid an excess of interactions: < 10)
  
  # Model H0/H1: the parameter m (character type) is equivalent to the model that generates H0/H1.
  # (Value by default: mvnorm / Possible values: mvnorm, simplex or multinom)
  
  # Y correlation: the parameter c (numerical type) is equivalent to the correlation of the Y variables.
  # (Value by default: 0 / Possible values: all rational numbers between -1.0 and 1.0)
  
  # Y variance: the parameter v (character type) is equivalent to the variance of the Y variables.
  # (Value by default: equal / Possible values: equal or unequal)
  
  # Changing Factor: the parameter w (character type) determines which factor changes.
  # (Value by default: B / Possible values: A, B or AB)
  
  # Number of simulations: the parameter S (numerical type) is equivalent to the number of simulations.
  # (Value by default: 1E3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 1E6)
  
  # Factors heatmap: the parameter plot (character type) determine whether or not to show the heatmap of the factor distribution.
  # (Value by default: F / Possible values: F or T)
  
  # Chunk number or k: the parameter k or Chunk number (numeric type) determine whether to compute tstats or ado.
  # (Value by default: 0 / Possible values: 0 (ado) or not 0 (tstats))
  
  # D, dd or DistDef: Multivariate non-normal distribution definition.
  # (Value by default: unif-0-1 / Possible values: unif-0-1)
  
  # l or lambda: lambda parameter (Poisson distribution) to generate size for 'multinom' generator model.
  # (Value by default: 1000 / Possible values: 1000)
  
  # s or stdev: stdev for the 'simplex' generator model.
  # (Value by default: 0.1 / Possible values: all rational numbers >0)
  
  # p, loc or position: location of the 'simplex' generator model.
  # (Value by default: 1 / Possible values: ?)
  
  # p_dist: Distribution for 'simplex' generator model.
  # (Value by default: norm / Possible values: norm, gamma or beta)
  
  # hk, H or heterosk: Heteoskedasticity degree (level 1).
  # (Value by default: 1 / Possible values: 1)
  
  # d or delta: H1 generation parameter.
  # (Value by default: 0 / Possible values: 0 for H0, > 0 for H1)
  
  # t or transf: Data transformation type.
  # (Value by default: none / Possible values: none, sqrt, log,...)
  # t_values = c("none", "sqrt", "log", "...")
  # transf_values <- t # Equivalencia
  
  # o or output: Output file name.
  # (Value by default: outputfiletemp / Possible values: NULL or outputfiletemp)
  # o = "outputfiletemp"
  
  # x or cpu: Number of cores.
  # (Value by default: 10 / Possible values: ?)
  
  # f or fx: Path to helper functions.
  # f = as.character(getwd())
  
  # Se crea una lista con estos parámetros y se asignan al entorno global.
  Parameter_List <- list("q" = q, "a" = a, "b" = b, "n" = n,
                         "u" = u, "m" = m, "c" = c, "v" = v,
                         "w" = w, "S" = S, "plotheatmap" = plotheatmap,
                         "k" = k, "D" = D, "l" = l, "s" = s, "p" = p, "p_dist" = p_dist,
                         "H" = H, "d" = d) # "t" = t, "o" = 0, "x" = x, "f" = f)
  list2env(Parameter_List, .GlobalEnv)
  
  Parameter_defin <- c("Dependent variables",
                       "Factor A Levels",
                       "Factor B Levels",
                       "Number of samples",
                       "Unbalance 1st B Level",
                       "Model H0/H1",
                       "Y correlation",
                       "Y variance",
                       "Changing Factor",
                       "Number of simulations",
                       "Factors heatmap",
                       "Chunk number",
                       "Multivariate non-normal distribution definition",
                       "Lambda parameter of the Poisson distribution",
                       "Stdev for the 'simplex' generator model",
                       "Location of the 'simplex' generator model",
                       "Distribution for 'simplex' generator model: norm, gamma or beta",
                       "Heteroskedasticity degree (level 1)",
                       "Delta (H1 generation parameter)")
                       # "Data transformation: sqrt, log or none",
                       # "Output file name",
                       # "Number of cores",
                       # "Path to helper functions")
  
  Parameter_names <- c("q",
                       "a",
                       "b",
                       "n",
                       "u",
                       "m / modelSim",
                       "c / Cor",
                       "v / Var",
                       "w",
                       "S",
                       "plotheatmap",
                       "k / chunk",
                       "D / DistDef / dd",
                       "l / Lambda",
                       "s /stdev",
                       "p / position / loc",
                       "p_dist / pdist",
                       "H / hk",
                       "d / delta")
                       # "t / transf",
                       # "Output file name",
                       # "x / cores",
                       # "f / fx")
  
  Parameter_values <- c()
  for (i in 1:length((Parameter_List))){
    Parameter_values <- append(Parameter_values, toString(Parameter_List[[i]]))}
  
  Parameter_class <- c()
  for (i in 1:length((Parameter_List))){
    Parameter_class <- append(Parameter_class, class(Parameter_List[[i]]))}
  
  Parameter_DF = data.frame(Parameter_defin, Parameter_names, Parameter_values, Parameter_class)
  colnames(Parameter_DF) = c("Parameter definition", "Parameter names", "Parameter value", "Parameter class type")
  
  # Se asignan los DF generados al entorno global.
  Parameter_DF_List <- list("Parameter_DF" = Parameter_DF, "Parameter_values" = Parameter_values,
                            "Parameter_names" = Parameter_names,"Parameter_class" = Parameter_class,
                            "Parameter definition" = Parameter_defin)
  list2env(Parameter_DF_List, .GlobalEnv)
  
  # -------------------------------
  
  
  # Niveles de los factores y su interacción
  # ----------------------------------------
  
  # Se obtienen los niveles de los factores y su interacción según los parámetros: n, a, b, u, w. 
  # Especificando si se ha de mostrar el mapa de calor con el valor lógico "plotheatmap".
  
  label <- function(a, b, n, u, w, plot = plotheatmap, seed = 1325) {
    
    set.seed(seed)
    
    if (w == "A") {
      ua <- u
      ub <- 1
    } else if (w == "B") {
      ua <- 1
      ub <- u
    } else if (w == "AB") {
      ua <- ub <- u
    }
    
    xa <- c(ua, rep(1, a-1))
    pa <- round(xa/sum(xa)*n)
    r <- n - sum(pa)
    sel <- sample(1:length(pa), size = abs(r))
    if (r > 0) {
      pa[sel] <- pa[sel] + 1
    } else if (r < 0) {
      pa[sel] <- pa[sel] - 1
    }
    A <- rep(1:a, times = pa)
    # A <- gl(a, n/a, length = n)
    
    B <- c()
    xb <- c(ub, rep(1, b-1))
    for (i in table(A)) {
      pb <- round(xb/sum(xb)*i)
      r <- i - sum(pb)
      sel <- sample(1:length(pb), size = abs(r))
      if (r > 0) {
        pb[sel] <- pb[sel] + 1
      } else if (r < 0) {
        pb[sel] <- pb[sel] - 1
      } 
      B <- c(B, rep(1:b, times = pb))
    }
    if (plotheatmap) {
      pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
    }
    return(list(factor(A), factor(B)))
  }
  
  labs <- label(a, b, n, u, w)
  A <- labs[[1]]
  B <- labs[[2]]
  
  if (w == "A") {
    ch <- A
  } else if (w == "B") {
    ch <- B
  } else if (w == "AB") {
    ch <- A:B
  } else {
    stop(sprintf("Unknown factor: '%s'.", w))
  }
  
  # Se asignan al entorno los parámetros generados en esta sección:
  Parameter_gen2_List <- list("ch" = ch, "A" = A, "B" = B)
  list2env(Parameter_gen2_List, .GlobalEnv)
  
  # ----------------------------------------
  
  
  # Simulación del conjunto de datos
  # y aplicación de MANTA y MANOVA
  # --------------------------------
  
  # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
  # las p del factor B para cada método y conjunto de datos considerado:
  Manta_pB = rep(NA, S)
  Manova_pB = rep(NA, S)
  Manta_Ylog_pB = rep(NA, S)
  Manova_Ylog_pB = rep(NA, S)
  Manta_Yclr_pB = rep(NA, S)
  Manova_Yclr_pB = rep(NA, S)
  Manta_Ysqrt_pB = rep(NA, S)
  Manova_Ysqrt_pB = rep(NA, S)
  Manta_Yscaled_pB = rep(NA, S)
  Manova_Yscaled_pB = rep(NA, S)
  Manta_YnormMinMax_pB = rep(NA, S)
  Manova_YnormMinMax_pB = rep(NA, S)
  
  # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
  # lostiempos de simulación para cada método y conjunto de datos considerado: 
  t_Manta_pB = rep(NA, S)
  t_Manova_pB = rep(NA, S)
  t_Manta_Ylog_pB = rep(NA, S)
  t_Manova_Ylog_pB = rep(NA, S)
  t_Manta_Yclr_pB = rep(NA, S)
  t_Manova_Yclr_pB = rep(NA, S)
  t_Manta_Ysqrt_pB = rep(NA, S)
  t_Manova_Ysqrt_pB = rep(NA, S)
  t_Manta_Yscaled_pB = rep(NA, S)
  t_Manova_Yscaled_pB = rep(NA, S)
  t_Manta_YnormMinMax_pB = rep(NA, S)
  t_Manova_YnormMinMax_pB = rep(NA, S)
  
  t_start_simulation <- Sys.time()
  
  writeLines(paste0("SIMULACIÓN ", Simul_count, "\n",
                    "\n- Variables simuladas: ",
                    "alpha = ", Alpha, " | S = ", S, " | n = ", n,  " | Var = ", Var,
                    " | Cor = ", Cor, " | delta = ", delta,
                    # " | q = ", q, " | loc = ", loc, " | stdev = ", stdev, " | lambda = ", lambda,
                    "\n- Inicio de la simulación: ", format(t_start_simulation, '%d/%m/%Y a las %H:%M:%S'),
                    "\n- Progreso:\n"))
  
  pb <- txtProgressBar(min = 0, max = S, style = 3, width = 50, char = "=") 
  
  for (i in 1:S) {
    
    Y <- Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor) # Cambiar mu = 0 <-> núm. grande (para transformaciones)
    
    # Transformaciones posibles de los datos simulados (Y).
    # => Para mvnorm no se aplicará ni Yscaled ni YnormMinMax ya que se genera otro estadístico diferente al que se está estudiando
    #    (p_valor) y que, a parte de que no tiene porque ser invariante a transformaciones, NO SIRVEN PARA COMPARAR MANTA-MANOVA.
    # => Para simplex solo se estudiará MANTA, así que no es necesario corregir ni Ylog ni Yclr.
    Ylog <- log(Y+1) # Transformación logarítmica. Se suma "1" para evitar "log(0) = -Inf" y valores < 0.
    Yclr <- (clr(Y))^2 # Transformación Centered Log Ratio (clr): clr(x) = (ln(x_i) - 1/D ∑_(j=1 -> D) ln(x_j))_i
                       # Resultados < 0 dan error al aplicar MANOVA y extraer "summary(manova(Yclr ~ A + B + A:B))$stats[2,6]":
                       # Error => Error in summary.manova(manova(Yclr ~ A + B + A:B)) : 
                       #          residuals have rank 2 < 3
    Ysqrt <- sqrt(Y) # Transformación sqrt.
    # Yscaled <- cbind(Y[,1]/sd(Y[,1]), Y[,2]/sd(Y[,2]), Y[,3]/sd(Y[,3])) # Escalado de datos: se divide cada valor por la desviación 
    #                                                                     # típica del conjunto de datos.
    # YnormMinMax <-  cbind((Y[,1]-min(Y[,1]))/(max(Y[,1])-min(Y[,1])), # Normalización Min- Max: si queremos que cada variable 
    #                                                                   # contribuya por igual al análisis.
    #                       (Y[,2]-min(Y[,2]))/(max(Y[,2])-min(Y[,2])),
    #                       (Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3])))
    
    # Verificación de NAs en los conjuntos de datos creados:
    # if (sum(is.na(Y)) != 0) {
    #   cat("\nEl conjunto de datos Y contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    #   }
    # if (sum(is.na(Ylog)) != 0) {
    #   cat("\nEl conjunto de datos Ylog contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    #   }
    # if (sum(is.na(Yclr)) != 0) {
    #   cat("\nEl conjunto de datos Ylog contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    #   }
    # if (sum(is.na(Ysqrt)) != 0) {
    #   cat("\nEl conjunto de datos Ysqrt contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    #   }
    # if (sum(is.na(Yscaled)) != 0) {
    #   cat("\nEl conjunto de datos Yscaled contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    #   }
    # if (sum(is.na(YnormMinMax)) != 0) {
    #   cat("\nEl conjunto de datos YnormMinMax contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    #   }
    
    # Aplicación de Manta y Manova para cada conjunto de datos (Y, Ylog, Ysqrt,...):
    
    # Y:
    t_Sim_start <- Sys.time()
    Manta_pB[i] <- manta(Y ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    t_Sim_start <- Sys.time()
    Manova_pB[i] <- summary(manova(Y ~ A + B + A:B))$stats[2,6]
    t_Manova_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # YLog:
    t_Sim_start <- Sys.time()
    Manta_Ylog_pB[i] <- manta(Ylog ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_Ylog_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    t_Sim_start <- Sys.time()
    Manova_Ylog_pB[i] <- summary(manova(Ylog ~ A + B + A:B))$stats[2,6]
    t_Manova_Ylog_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # Yclr:
    t_Sim_start <- Sys.time()
    Manta_Yclr_pB[i] <- manta(Yclr ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_Yclr_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    t_Sim_start <- Sys.time()
    Manova_Yclr_pB[i] <- summary(manova(Yclr ~ A + B + A:B))$stats[2,6]
    t_Manova_Yclr_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # Ysqrt:
    t_Sim_start <- Sys.time()
    Manta_Ysqrt_pB[i] <- manta(Ysqrt ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_Ysqrt_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    t_Sim_start <- Sys.time()
    Manova_Ysqrt_pB[i] <- summary(manova(Ysqrt ~ A + B + A:B))$stats[2,6]
    t_Manova_Ysqrt_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # # Yscaled:
    # t_Sim_start <- Sys.time()
    # Manta_Yscaled_pB[i] <- manta(Yscaled ~ A + B + A:B)$aov.tab[2,6]
    # t_Manta_Yscaled_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # t_Sim_start <- Sys.time()
    # Manova_Yscaled_pB[i] <- summary(manova(Yscaled ~ A + B + A:B))$stats[2,6]
    # t_Manova_Yscaled_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # # YnormMinMax:
    # t_Sim_start <- Sys.time()
    # Manta_YnormMinMax_pB[i] <- manta(YnormMinMax ~ A + B + A:B)$aov.tab[2,6]
    # t_Manta_YnormMinMax_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # t_Sim_start <- Sys.time()
    # Manova_YnormMinMax_pB[i] <- summary(manova(YnormMinMax ~ A + B + A:B))$stats[2,6]
    # t_Manova_YnormMinMax_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  t_end_simulation <- Sys.time()
  
  writeLines(paste0("\n- Tiempo de simulación: ",
                    round(as.numeric(t_end_simulation -
                                       t_start_simulation, units = "secs"), 0),
                    " s.\n\n"))
  
  # --------------------------------
  
  
  # Cálculo de la potencia y los tiempos de computación
  # ---------------------------------------------------
  
  # Se calculan los tiempos de computación para cada método:
  T_Manta_pB = sum(t_Manta_pB)
  T_Manova_pB = sum(t_Manova_pB)
  T_Manta_Ylog_pB = sum(t_Manta_Ylog_pB)
  T_Manova_Ylog_pB = sum(t_Manova_Ylog_pB)
  T_Manta_Yclr_pB = sum(t_Manta_Yclr_pB)
  T_Manova_Yclr_pB = sum(t_Manova_Yclr_pB)
  T_Manta_Ysqrt_pB = sum(t_Manta_Ysqrt_pB)
  T_Manova_Ysqrt_pB = sum(t_Manova_Ysqrt_pB)
  T_Manta_Yscaled_pB = sum(t_Manta_Yscaled_pB)
  T_Manova_Yscaled_pB = sum(t_Manova_Yscaled_pB)
  T_Manta_YnormMinMax_pB = sum(t_Manta_YnormMinMax_pB)
  T_Manova_YnormMinMax_pB = sum(t_Manova_YnormMinMax_pB)
  
  # Se calcula la potencia de cada método para el α dado:
  t_Sim_start <- Sys.time()
  Manta_pot <- mean(Manta_pB < Alpha)
  t_Manta_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_pB
  t_Sim_start <- Sys.time()
  Manova_pot <- mean(Manova_pB < Alpha)
  t_Manova_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_pB
  Y_Manta_Manova_pot <- cbind(Manta_pot, t_Manta_pot, Manova_pot, t_Manova_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Ylog_pot <- mean(Manta_Ylog_pB < Alpha)
  t_Manta_Ylog_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Ylog_pB
  t_Sim_start <- Sys.time()
  Manova_Ylog_pot <- mean(Manova_Ylog_pB < Alpha)
  t_Manova_Ylog_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Ylog_pB
  Ylog_Manta_Manova_pot <- cbind(Manta_Ylog_pot, t_Manta_Ylog_pot, Manova_Ylog_pot, t_Manova_Ylog_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Yclr_pot <- mean(Manta_Yclr_pB < Alpha)
  t_Manta_Yclr_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Yclr_pB
  t_Sim_start <- Sys.time()
  Manova_Yclr_pot <- mean(Manova_Yclr_pB < Alpha)
  t_Manova_Yclr_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Yclr_pB
  Yclr_Manta_Manova_pot <- cbind(Manta_Yclr_pot, t_Manta_Yclr_pot, Manova_Yclr_pot, t_Manova_Yclr_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Ysqrt_pot <- mean(Manta_Ysqrt_pB < Alpha)
  t_Manta_Ysqrt_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Ysqrt_pB
  t_Sim_start <- Sys.time()
  Manova_Ysqrt_pot <- mean(Manova_Ysqrt_pB < Alpha)
  t_Manova_Ysqrt_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Ysqrt_pB
  Ysqrt_Manta_Manova_pot <- cbind(Manta_Ysqrt_pot, t_Manta_Ysqrt_pot, Manova_Ysqrt_pot, t_Manova_Ysqrt_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Yscaled_pot <- mean(Manta_Yscaled_pB < Alpha)
  t_Manta_Yscaled_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Yscaled_pB
  t_Sim_start <- Sys.time()
  Manova_Yscaled_pot <- mean(Manova_Yscaled_pB < Alpha)
  t_Manova_Yscaled_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Yscaled_pB
  Yscaled_Manta_Manova_pot <- cbind(Manta_Yscaled_pot, t_Manta_Yscaled_pot, Manova_Yscaled_pot, t_Manova_Yscaled_pot)
  
  t_Sim_start <- Sys.time()
  Manta_YnormMinMax_pot <- mean(Manta_YnormMinMax_pB < Alpha)
  t_Manta_YnormMinMax_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_YnormMinMax_pB
  t_Sim_start <- Sys.time()
  Manova_YnormMinMax_pot <- mean(Manova_YnormMinMax_pB < Alpha)
  t_Manova_YnormMinMax_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_YnormMinMax_pB
  YnormMinMax_Manta_Manova_pot <- cbind(Manta_YnormMinMax_pot, t_Manta_YnormMinMax_pot, Manova_YnormMinMax_pot, t_Manova_YnormMinMax_pot)
  
  # ---------------------------------------------------
  
  
  # DF finales
  # ----------
  
  # DF que almacena la potencia de cada método:
  DF_CompPot_res_Ytype <- c("Datos sin transformar", "Transformación logarítmica", "Transformación Centered Log Ratio (clr)",
                    "Transformación raíz cuadrada", "Escalado de datos por desviación típica", "Normalización Min- Max")
  DF_CompPot_res_modelSim_iter <- rep(modelSim, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_Alpha_iter <- rep(Alpha, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_S_iter <- rep(S, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_n_iter <- rep(n, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_Var_iter <- rep(Var, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_q_iter <- rep(q, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_loc_iter <- rep(loc, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_stdev_iter <- rep(stdev, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_Cor_iter <- rep(Cor, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_delta_iter <- rep(delta, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_lambda_iter <- rep(lambda, length(DF_CompPot_res_Ytype))
  DF_CompPot_res <- data.frame(DF_CompPot_res_Ytype, DF_CompPot_res_modelSim_iter, DF_CompPot_res_Alpha_iter,
                               DF_CompPot_res_S_iter, DF_CompPot_res_n_iter, DF_CompPot_res_Var_iter,
                               DF_CompPot_res_q_iter, DF_CompPot_res_loc_iter, DF_CompPot_res_stdev_iter, 
                               DF_CompPot_res_Cor_iter, DF_CompPot_res_delta_iter, DF_CompPot_res_lambda_iter,
                               rbind(Y_Manta_Manova_pot, Ylog_Manta_Manova_pot, Yclr_Manta_Manova_pot, Ysqrt_Manta_Manova_pot,
                                     Yscaled_Manta_Manova_pot, YnormMinMax_Manta_Manova_pot))
  colnames(DF_CompPot_res) <- c("Datos", "Modelo", "alpha", "S", "n", "Var", "q", "loc", "stdev",
                                "Cor", "delta", "lambda", "Pot_MANTA", "tcomp_MANTA", "Pot_MANOVA", "tcomp_MANOVA")
  
  # Se asignan los parámetros generados en esta sección:
  DF_CompPot_res_List <- list("DF_CompPot_res" = DF_CompPot_res)
  list2env(DF_CompPot_res_List, .GlobalEnv)
  
  # ----------
  
  }

# Función de asignación de la "stdev" para el modelo "simplex:
funcStdev_Simplex <- function(q, loc) {
  
  # Se carga la tabla correspondiente de "stdev":
  tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", f, pdist), h = T)
  colnames(tbl) <- c("Q", "L", "S")
  tbl_name <- "simplex_table_Q_L_S"
  write.csv(tbl, file = paste(tbl_name, ".csv", sep = ""), row.names = FALSE)
  
  if (! q %in% unique(tbl$Q)) {
    stop(sprintf("stdev not precomputed for q = %s", q))}
  
  # Se busca la "stdev" que corresponde a los valores [q, loc] seleccionados.
  stdev <- subset(tbl, Q == q & L == loc)$S
  
  # Se asignan los parámetros generados en esta sección:
  funcStdev_Simplex_ResList <- list("stdev" = stdev)
  list2env(funcStdev_Simplex_ResList, .GlobalEnv)
  #return(list(funcStdev_Simplex_ResList, list2env(funcStdev_Simplex_ResList, .GlobalEnv)))
  
}

# Función que calcula el "loc" necesario para tener una distribución "simplex" dada
# (de tamaño f(q)) .τ. Σsimplex = 1:
funcWhatLoc_Simplex <- function(x) {
  
  # x ≡ valor buscado para la primera posición del vector "simplex"
  #   => y0 = {x, ...} .τ. y0[1] = x
   
  loc_calc <- (x * (q-1)) / (1 - x)
  
  # Se asignan los parámetros generados en esta sección:
  funcWhatLoc_Simplex_ResList <- list("loc_calc" = loc_calc)
  list2env(funcWhatLoc_Simplex_ResList, .GlobalEnv)
  
}

# Función que calcula el "loc" necesario para tener una distribución "simplex" dada
# (de tamaño f(q)) .τ. Σsimplex = 1:
funcWhatY0_Simplex <- function(q, loc) {
  
  # Distribución de "simplex" según q y loc:
  
  x <- c(loc, rep(1, q-1))
  Y0_calc <- x/sum(x)
  
  # Se asignan los parámetros generados en esta sección:
  funcWhatY0_Simplex_ResList <- list("Y0_1_calc" = Y0_calc)
  list2env(funcWhatY0_Simplex_ResList, .GlobalEnv)
  
}

# Función que simula el conjunto de datos Y y aplica el método MANTA
# para el modelo "simplex":
CompMantaManova_simplex <- function() {
  
  # Instalación y carga de los paquetes necesarios
  # ----------------------------------------------
  
  # Se cargan los paquetes, instalándolos si es necesario.
  
  #    require() => a warning if a package is not installed and then continue to execute the code
  #    library() => an error and stop the execution of the code
  
  if (!require("devtools")) install.packages("devtools")
  require(devtools)
  library(devtools)
  
  if (!require("versions")) install.packages("versions")
  require(versions)
  library(versions)
  
  if (!require("manta")) devtools::install_github("dgarrimar/manta")
  require(manta)
  library(manta)
  
  package_list <- c("CompQuadForm", "car", "Hmisc", "MASS", "copula", "vegan",
                    "randomcoloR", "plyr", "svMisc", "reshape2", "hms", "progress",
                    "ggplot2", "arm", "pheatmap", "latex2exp", "tibble", "compositions",
                    "glue")
  
  # Hmisc => Necesario para la función label
  # randomcoloR => La función randomcoloR() genera un color aleatóriamente
  # plyr => Necesario para la función mapvalues()
  # svMisc => Permite usar progress() para indicar visualmente el progreso de la simulación
  # reshape2 => La función melt() permite convertir un objeto "dist" en "DF"
  # hms => La función as_hms() permite convertir "difftime" en DD:HH:MM:SS
  # progress => Para mostrar una barra de progreso de la simulación
  # latex2exp => Toma una cadena LaTeX, la analiza y devuelve la expresión plotmath más cercana
  # glue => Permite pasar expresiones (símbolos griegos, equaciones, etc.) a las etiquetas de facet.grid
  
  for (i in package_list){
    if(! i %in% installed.packages()){
      install.packages(i)
    }
    require(i, character.only = TRUE)
    library(i, character.only = TRUE)
  }
  
  # ----------------------------------------------
  
  
  # DF de los parámetros del modelo
  # -------------------------------
  
  # Estos son: q, a, b, n, u, m, c, v, w, S, plotheatmap, k, D, l, s, p, p_dist, H, d, t, x, f.
  
  # Dependent variables: the parameter q (numerical type) is equivalent to the number of dependent variables.
  # (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 4)
  
  # Factor A Levels: the parameter a (numerical type) is equivalent to the number of levels of factor A.
  # (Value by default: 2 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
  
  # Factor B Levels: the parameter b (numerical type) is equivalent to the number of levels of factor B."
  # (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
  
  # Number of samples: the parameter n (numerical type) is equivalent to the number of samples.
  # (Value by default: 300 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 400)
  
  # Unbalance 1st B Level: the parameter u (numerical type) is equivalent to the unbalance degree of the first level regarding the rest of the levels of the factor B.
  # (Value by default: 1 / Possible values: all positive rational numbers / Max. recommended to avoid an excess of interactions: < 10)
  
  # Model H0/H1: the parameter m (character type) is equivalent to the model that generates H0/H1.
  # (Value by default: mvnorm / Possible values: mvnorm, simplex or multinom)
  
  # Y correlation: the parameter c (numerical type) is equivalent to the correlation of the Y variables.
  # (Value by default: 0 / Possible values: all rational numbers between -1.0 and 1.0)
  
  # Y variance: the parameter v (character type) is equivalent to the variance of the Y variables.
  # (Value by default: equal / Possible values: equal or unequal)
  
  # Changing Factor: the parameter w (character type) determines which factor changes.
  # (Value by default: B / Possible values: A, B or AB)
  
  # Number of simulations: the parameter S (numerical type) is equivalent to the number of simulations.
  # (Value by default: 1E3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 1E6)
  
  # Factors heatmap: the parameter plot (character type) determine whether or not to show the heatmap of the factor distribution.
  # (Value by default: F / Possible values: F or T)
  
  # Chunk number or k: the parameter k or Chunk number (numeric type) determine whether to compute tstats or ado.
  # (Value by default: 0 / Possible values: 0 (ado) or not 0 (tstats))
  
  # D, dd or DistDef: Multivariate non-normal distribution definition.
  # (Value by default: unif-0-1 / Possible values: unif-0-1)
  
  # l or lambda: lambda parameter (Poisson distribution) to generate size for 'multinom' generator model.
  # (Value by default: 1000 / Possible values: 1000)
  
  # s or stdev: stdev for the 'simplex' generator model.
  # (Value by default: 0.1 / Possible values: all rational numbers >0)
  
  # p, loc or position: location of the 'simplex' generator model.
  # (Value by default: 1 / Possible values: ?)
  
  # p_dist: Distribution for 'simplex' generator model.
  # (Value by default: norm / Possible values: norm, gamma or beta)
  
  # hk, H or heterosk: Heteoskedasticity degree (level 1).
  # (Value by default: 1 / Possible values: 1)
  
  # d or delta: H1 generation parameter.
  # (Value by default: 0 / Possible values: 0 for H0, > 0 for H1)
  
  # t or transf: Data transformation type.
  # (Value by default: none / Possible values: none, sqrt, log,...)
  # t_values = c("none", "sqrt", "log", "...")
  # transf_values <- t # Equivalencia
  
  # o or output: Output file name.
  # (Value by default: outputfiletemp / Possible values: NULL or outputfiletemp)
  # o = "outputfiletemp"
  
  # x or cpu: Number of cores.
  # (Value by default: 10 / Possible values: ?)
  
  # f or fx: Path to helper functions.
  # f = as.character(getwd())
  
  # Se crea una lista con estos parámetros y se asignan al entorno global.
  Parameter_List <- list("q" = q, "a" = a, "b" = b, "n" = n,
                         "u" = u, "m" = m, "c" = c, "v" = v,
                         "w" = w, "S" = S, "plotheatmap" = plotheatmap,
                         "k" = k, "D" = D, "l" = l, "s" = s, "p" = p, "p_dist" = p_dist,
                         "H" = H, "d" = d) # "t" = t, "o" = 0, "x" = x, "f" = f)
  list2env(Parameter_List, .GlobalEnv)
  
  Parameter_defin <- c("Dependent variables",
                       "Factor A Levels",
                       "Factor B Levels",
                       "Number of samples",
                       "Unbalance 1st B Level",
                       "Model H0/H1",
                       "Y correlation",
                       "Y variance",
                       "Changing Factor",
                       "Number of simulations",
                       "Factors heatmap",
                       "Chunk number",
                       "Multivariate non-normal distribution definition",
                       "Lambda parameter of the Poisson distribution",
                       "Stdev for the 'simplex' generator model",
                       "Location of the 'simplex' generator model",
                       "Distribution for 'simplex' generator model: norm, gamma or beta",
                       "Heteroskedasticity degree (level 1)",
                       "Delta (H1 generation parameter)")
  # "Data transformation: sqrt, log or none",
  # "Output file name",
  # "Number of cores",
  # "Path to helper functions")
  
  Parameter_names <- c("q",
                       "a",
                       "b",
                       "n",
                       "u",
                       "m / modelSim",
                       "c / Cor",
                       "v / Var",
                       "w",
                       "S",
                       "plotheatmap",
                       "k / chunk",
                       "D / DistDef / dd",
                       "l / Lambda",
                       "s /stdev",
                       "p / position / loc",
                       "p_dist / pdist",
                       "H / hk",
                       "d / delta")
  # "t / transf",
  # "Output file name",
  # "x / cores",
  # "f / fx")
  
  Parameter_values <- c()
  for (i in 1:length((Parameter_List))){
    Parameter_values <- append(Parameter_values, toString(Parameter_List[[i]]))}
  
  Parameter_class <- c()
  for (i in 1:length((Parameter_List))){
    Parameter_class <- append(Parameter_class, class(Parameter_List[[i]]))}
  
  Parameter_DF = data.frame(Parameter_defin, Parameter_names, Parameter_values, Parameter_class)
  colnames(Parameter_DF) = c("Parameter definition", "Parameter names", "Parameter value", "Parameter class type")
  
  # Se asignan los DF generados al entorno global.
  Parameter_DF_List <- list("Parameter_DF" = Parameter_DF, "Parameter_values" = Parameter_values,
                            "Parameter_names" = Parameter_names,"Parameter_class" = Parameter_class,
                            "Parameter definition" = Parameter_defin)
  list2env(Parameter_DF_List, .GlobalEnv)
  
  # -------------------------------
  
  
  # Niveles de los factores y su interacción
  # ----------------------------------------
  
  # Se obtienen los niveles de los factores y su interacción según los parámetros: n, a, b, u, w. 
  # Especificando si se ha de mostrar el mapa de calor con el valor lógico "plotheatmap".
  
  label <- function(a, b, n, u, w, plot = plotheatmap, seed = 1325) {
    
    set.seed(seed)
    
    if (w == "A") {
      ua <- u
      ub <- 1
    } else if (w == "B") {
      ua <- 1
      ub <- u
    } else if (w == "AB") {
      ua <- ub <- u
    }
    
    xa <- c(ua, rep(1, a-1))
    pa <- round(xa/sum(xa)*n)
    r <- n - sum(pa)
    sel <- sample(1:length(pa), size = abs(r))
    if (r > 0) {
      pa[sel] <- pa[sel] + 1
    } else if (r < 0) {
      pa[sel] <- pa[sel] - 1
    }
    A <- rep(1:a, times = pa)
    # A <- gl(a, n/a, length = n)
    
    B <- c()
    xb <- c(ub, rep(1, b-1))
    for (i in table(A)) {
      pb <- round(xb/sum(xb)*i)
      r <- i - sum(pb)
      sel <- sample(1:length(pb), size = abs(r))
      if (r > 0) {
        pb[sel] <- pb[sel] + 1
      } else if (r < 0) {
        pb[sel] <- pb[sel] - 1
      } 
      B <- c(B, rep(1:b, times = pb))
    }
    if (plotheatmap) {
      pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
    }
    return(list(factor(A), factor(B)))
  }
  
  labs <- label(a, b, n, u, w)
  A <- labs[[1]]
  B <- labs[[2]]
  
  if (w == "A") {
    ch <- A
  } else if (w == "B") {
    ch <- B
  } else if (w == "AB") {
    ch <- A:B
  } else {
    stop(sprintf("Unknown factor: '%s'.", w))
  }
  
  # Se asignan al entorno los parámetros generados en esta sección:
  Parameter_gen2_List <- list("ch" = ch, "A" = A, "B" = B)
  list2env(Parameter_gen2_List, .GlobalEnv)
  
  # ----------------------------------------
  
  
  # Asignación de la "stdev" para el modelo "simplex"
  # -------------------------------------------------
  
  funcStdev_Simplex(q, loc)
  s <- stdev
  
  # -------------------------------------------------
  
  
  # Simulación del conjunto de datos
  # y aplicación de MANTA y MANOVA
  # --------------------------------
  
  # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
  # las p del factor B para cada método y conjunto de datos considerado:
  Manta_pB = rep(NA, S)
  Manova_pB = rep(NA, S)
  Manta_Ylog_pB = rep(NA, S)
  Manova_Ylog_pB = rep(NA, S)
  Manta_Yclr_pB = rep(NA, S)
  Manova_Yclr_pB = rep(NA, S)
  Manta_Ysqrt_pB = rep(NA, S)
  Manova_Ysqrt_pB = rep(NA, S)
  Manta_Yscaled_pB = rep(NA, S)
  Manova_Yscaled_pB = rep(NA, S)
  Manta_YnormMinMax_pB = rep(NA, S)
  Manova_YnormMinMax_pB = rep(NA, S)
  
  # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
  # lostiempos de simulación para cada método y conjunto de datos considerado: 
  t_Manta_pB = rep(NA, S)
  t_Manova_pB = rep(NA, S)
  t_Manta_Ylog_pB = rep(NA, S)
  t_Manova_Ylog_pB = rep(NA, S)
  t_Manta_Yclr_pB = rep(NA, S)
  t_Manova_Yclr_pB = rep(NA, S)
  t_Manta_Ysqrt_pB = rep(NA, S)
  t_Manova_Ysqrt_pB = rep(NA, S)
  t_Manta_Yscaled_pB = rep(NA, S)
  t_Manova_Yscaled_pB = rep(NA, S)
  t_Manta_YnormMinMax_pB = rep(NA, S)
  t_Manova_YnormMinMax_pB = rep(NA, S)
  
  t_start_simulation <- Sys.time()
  
  writeLines(paste0("SIMULACIÓN ", Simul_count, "\n",
                    "\n- Variables simuladas: ",
                    "alpha = ", Alpha, " | S = ", S, " | n = ", n,  " | Var = ", Var,
                    " | q = ", q, " | loc = ", loc, " | stdev = ", stdev,
                    " | delta = ", delta, # " | Cor = ", Cor, " | lambda = ", lambda,
                    "\n- Inicio de la simulación: ", format(t_start_simulation, '%d/%m/%Y a las %H:%M:%S'),
                    "\n- Progreso:\n"))
  
  pb <- txtProgressBar(min = 0, max = S, style = 3, width = 50, char = "=")
  
  for (i in 1:S) {
    
    # if (i == 1) { # Sanity check
    #   
    #   tol <- 0.01 # "tol <- 0.1" no arregla el problema "Error in -Y[i, k] + tol : non-numeric argument to binary operator"
    #               # Por defecto: 0.01
    #   check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, pdist)
    #   wm <- which.max(check$exp)
    #   if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
    #     stop("Deviation from expected centroid greater than tolerance.")
    #   }
    #   
    # }
    
      Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = F, pdist)
    
    # Transformaciones posibles de los datos simulados (Y).
    # => Para mvnorm no se aplicará ni Yscaled ni YnormMinMax ya que se genera otro estadístico diferente al que se está estudiando
    #    (p_valor) y que, a parte de que no tiene porque ser invariante a transformaciones, NO SIRVEN PARA COMPARAR MANTA-MANOVA.
    # => Para simplex solo se estudiará MANTA, así que no es necesario corregir ni Ylog ni Yclr.
    Ylog <- log(Y+1) # Transformación logarítmica. Se suma "1" para evitar "log(0) = -Inf" y valores < 0.
    Yclr <- clr(Y) # Transformación Centered Log Ratio (clr): clr(x) = (ln(x_i) - 1/D ∑_(j=1 -> D) ln(x_j))_i
                   # Resultados < 0 dan error al aplicar MANOVA y extraer "summary(manova(Yclr ~ A + B + A:B))$stats[2,6]":
                   # Error => Error in summary.manova(manova(Yclr ~ A + B + A:B)) : 
                   #          residuals have rank 2 < 3
    Ysqrt <- sqrt(Y) # Transformación sqrt.
    # Yscaled <- cbind(Y[,1]/sd(Y[,1]), Y[,2]/sd(Y[,2]), Y[,3]/sd(Y[,3])) # Escalado de datos: se divide cada valor por la desviación 
    #                                                                     # típica del conjunto de datos.
    # YnormMinMax <-  cbind((Y[,1]-min(Y[,1]))/(max(Y[,1])-min(Y[,1])), # Normalización Min- Max: si queremos que cada variable 
    #                                                                   # contribuya por igual al análisis.
    #                       (Y[,2]-min(Y[,2]))/(max(Y[,2])-min(Y[,2])),
    #                       (Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3])))
    
    # # Verificación de NAs en los conjuntos de datos creados:
    # if (sum(is.na(Y)) != 0) {
    #   cat("\nEl conjunto de datos Y contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    # }
    # if (sum(is.na(Ylog)) != 0) {
    #   cat("\nEl conjunto de datos Ylog contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    # }
    # if (sum(is.na(Yclr)) != 0) {
    #   cat("\nEl conjunto de datos Ylog contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    # }
    # if (sum(is.na(Ysqrt)) != 0) {
    #   cat("\nEl conjunto de datos Ysqrt contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    # }
    # if (sum(is.na(Yscaled)) != 0) {
    #   cat("\nEl conjunto de datos Yscaled contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    # }
    # if (sum(is.na(YnormMinMax)) != 0) {
    #   cat("\nEl conjunto de datos YnormMinMax contiene al menos un valor NA.\n\n")
    #   #next # O "stop".
    # }
    
    # Aplicación de Manta y Manova para cada conjunto de datos (Y, Ylog, Ysqrt,...):
    
    # Y:
    t_Sim_start <- Sys.time()
    Manta_pB[i] <- manta(Y ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # t_Sim_start <- Sys.time()
    # Manova_pB[i] <- summary(manova(Y ~ A + B + A:B))$stats[2,6]
    # t_Manova_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # YLog:
    t_Sim_start <- Sys.time()
    Manta_Ylog_pB[i] <- manta(Ylog ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_Ylog_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # t_Sim_start <- Sys.time()
    # Manova_Ylog_pB[i] <- summary(manova(Ylog ~ A + B + A:B))$stats[2,6]
    # t_Manova_Ylog_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # Yclr:
    t_Sim_start <- Sys.time()
    Manta_Yclr_pB[i] <- manta(Yclr ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_Yclr_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # t_Sim_start <- Sys.time()
    # Manova_Yclr_pB[i] <- summary(manova(Yclr ~ A + B + A:B))$stats[2,6]
    # t_Manova_Yclr_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # Ysqrt:
    t_Sim_start <- Sys.time()
    Manta_Ysqrt_pB[i] <- manta(Ysqrt ~ A + B + A:B)$aov.tab[2,6]
    t_Manta_Ysqrt_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # t_Sim_start <- Sys.time()
    # Manova_Ysqrt_pB[i] <- summary(manova(Ysqrt ~ A + B + A:B))$stats[2,6]
    # t_Manova_Ysqrt_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    # # Yscaled:
    # t_Sim_start <- Sys.time()
    # Manta_Yscaled_pB[i] <- manta(Yscaled ~ A + B + A:B)$aov.tab[2,6]
    # t_Manta_Yscaled_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # # t_Sim_start <- Sys.time()
    # # Manova_Yscaled_pB[i] <- summary(manova(Yscaled ~ A + B + A:B))$stats[2,6]
    # # t_Manova_Yscaled_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # 
    # # YnormMinMax:
    # t_Sim_start <- Sys.time()
    # Manta_YnormMinMax_pB[i] <- manta(YnormMinMax ~ A + B + A:B)$aov.tab[2,6]
    # t_Manta_YnormMinMax_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    # # t_Sim_start <- Sys.time()
    # # Manova_YnormMinMax_pB[i] <- summary(manova(YnormMinMax ~ A + B + A:B))$stats[2,6]
    # # t_Manova_YnormMinMax_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  t_end_simulation <- Sys.time()
  
  writeLines(paste0("\n- Tiempo de simulación: ",
                    round(as.numeric(t_end_simulation -
                                       t_start_simulation, units = "secs"), 0),
                    " s.\n\n"))
  
  # --------------------------------
  
  
  # Cálculo de la potencia y los tiempos de computación
  # ---------------------------------------------------
  
  # Se calculan los tiempos de computación para cada método:
  T_Manta_pB = sum(t_Manta_pB)
  T_Manova_pB = sum(t_Manova_pB)
  T_Manta_Ylog_pB = sum(t_Manta_Ylog_pB)
  T_Manova_Ylog_pB = sum(t_Manova_Ylog_pB)
  T_Manta_Yclr_pB = sum(t_Manta_Yclr_pB)
  T_Manova_Yclr_pB = sum(t_Manova_Yclr_pB)
  T_Manta_Ysqrt_pB = sum(t_Manta_Ysqrt_pB)
  T_Manova_Ysqrt_pB = sum(t_Manova_Ysqrt_pB)
  T_Manta_Yscaled_pB = sum(t_Manta_Yscaled_pB)
  T_Manova_Yscaled_pB = sum(t_Manova_Yscaled_pB)
  T_Manta_YnormMinMax_pB = sum(t_Manta_YnormMinMax_pB)
  T_Manova_YnormMinMax_pB = sum(t_Manova_YnormMinMax_pB)
  
  # Se calcula la potencia de cada método para el α dado:
  t_Sim_start <- Sys.time()
  Manta_pot <- mean(Manta_pB < Alpha)
  t_Manta_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_pB
  t_Sim_start <- Sys.time()
  Manova_pot <- mean(Manova_pB < Alpha)
  t_Manova_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_pB
  Y_Manta_Manova_pot <- cbind(Manta_pot, t_Manta_pot, Manova_pot, t_Manova_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Ylog_pot <- mean(Manta_Ylog_pB < Alpha)
  t_Manta_Ylog_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Ylog_pB
  t_Sim_start <- Sys.time()
  Manova_Ylog_pot <- mean(Manova_Ylog_pB < Alpha)
  t_Manova_Ylog_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Ylog_pB
  Ylog_Manta_Manova_pot <- cbind(Manta_Ylog_pot, t_Manta_Ylog_pot, Manova_Ylog_pot, t_Manova_Ylog_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Yclr_pot <- mean(Manta_Yclr_pB < Alpha)
  t_Manta_Yclr_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Yclr_pB
  t_Sim_start <- Sys.time()
  Manova_Yclr_pot <- mean(Manova_Yclr_pB < Alpha)
  t_Manova_Yclr_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Yclr_pB
  Yclr_Manta_Manova_pot <- cbind(Manta_Yclr_pot, t_Manta_Yclr_pot, Manova_Yclr_pot, t_Manova_Yclr_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Ysqrt_pot <- mean(Manta_Ysqrt_pB < Alpha)
  t_Manta_Ysqrt_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Ysqrt_pB
  t_Sim_start <- Sys.time()
  Manova_Ysqrt_pot <- mean(Manova_Ysqrt_pB < Alpha)
  t_Manova_Ysqrt_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Ysqrt_pB
  Ysqrt_Manta_Manova_pot <- cbind(Manta_Ysqrt_pot, t_Manta_Ysqrt_pot, Manova_Ysqrt_pot, t_Manova_Ysqrt_pot)
  
  t_Sim_start <- Sys.time()
  Manta_Yscaled_pot <- mean(Manta_Yscaled_pB < Alpha)
  t_Manta_Yscaled_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Yscaled_pB
  t_Sim_start <- Sys.time()
  Manova_Yscaled_pot <- mean(Manova_Yscaled_pB < Alpha)
  t_Manova_Yscaled_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Yscaled_pB
  Yscaled_Manta_Manova_pot <- cbind(Manta_Yscaled_pot, t_Manta_Yscaled_pot, Manova_Yscaled_pot, t_Manova_Yscaled_pot)
  
  t_Sim_start <- Sys.time()
  Manta_YnormMinMax_pot <- mean(Manta_YnormMinMax_pB < Alpha)
  t_Manta_YnormMinMax_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_YnormMinMax_pB
  t_Sim_start <- Sys.time()
  Manova_YnormMinMax_pot <- mean(Manova_YnormMinMax_pB < Alpha)
  t_Manova_YnormMinMax_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_YnormMinMax_pB
  YnormMinMax_Manta_Manova_pot <- cbind(Manta_YnormMinMax_pot, t_Manta_YnormMinMax_pot, Manova_YnormMinMax_pot, t_Manova_YnormMinMax_pot)
  
  # ---------------------------------------------------
  
  
  # DF finales
  # ----------
  
  # DF que almacena la potencia de cada método:
  DF_CompPot_res_Ytype <- c("Datos sin transformar", "Transformación logarítmica", "Transformación Centered Log Ratio (clr)",
                            "Transformación raíz cuadrada", "Escalado de datos por desviación típica", "Normalización Min- Max")
  DF_CompPot_res_modelSim_iter <- rep(modelSim, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_Alpha_iter <- rep(Alpha, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_S_iter <- rep(S, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_n_iter <- rep(n, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_Var_iter <- rep(Var, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_q_iter <- rep(q, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_loc_iter <- rep(loc, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_stdev_iter <- rep(stdev, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_Cor_iter <- rep(Cor, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_delta_iter <- rep(delta, length(DF_CompPot_res_Ytype))
  DF_CompPot_res_lambda_iter <- rep(lambda, length(DF_CompPot_res_Ytype))
  DF_CompPot_res <- data.frame(DF_CompPot_res_Ytype, DF_CompPot_res_modelSim_iter, DF_CompPot_res_Alpha_iter,
                               DF_CompPot_res_S_iter, DF_CompPot_res_n_iter, DF_CompPot_res_Var_iter,
                               DF_CompPot_res_q_iter, DF_CompPot_res_loc_iter, DF_CompPot_res_stdev_iter, 
                               DF_CompPot_res_Cor_iter, DF_CompPot_res_delta_iter, DF_CompPot_res_lambda_iter,
                               rbind(Y_Manta_Manova_pot, Ylog_Manta_Manova_pot, Yclr_Manta_Manova_pot, Ysqrt_Manta_Manova_pot,
                                     Yscaled_Manta_Manova_pot, YnormMinMax_Manta_Manova_pot))
  colnames(DF_CompPot_res) <- c("Datos", "Modelo", "alpha", "S", "n", "Var", "q", "loc", "stdev",
                                "Cor", "delta", "lambda", "Pot_MANTA", "tcomp_MANTA", "Pot_MANOVA", "tcomp_MANOVA")
  
  # Se asignan los parámetros generados en esta sección:
  DF_CompPot_res_List <- list("DF_CompPot_res" = DF_CompPot_res)
  list2env(DF_CompPot_res_List, .GlobalEnv)
  
  # ----------
  
}

# Función que simula el conjunto de datos Y y aplica el método MANTA
# para el modelo "multinom":
# ¡¡¡POR AHORA NO!!!


############################# 
# Funciones gráficas mvnorm #
#############################

# Función para crear las gráficas "∆ vs. Potencia" según el método usado.
# Agrupación de "color" ≡ Método - Cor
DeltaPotGraph_CorColor <- function(DF_Results) {
  
  # Genera la gráfica "∆ vs. Potencia" para los diferentes Y considerados en el método dado.
  
  DF_graph <- DF_Results[DF_Results$Datos == Datos_graph
                         & DF_Results$S == S_graph
                         & DF_Results$alpha == alpha_graph
                         & DF_Results$Var == Var_graph
                         & DF_Results$Cor > Cor_graph - 0.1 & DF_Results$Cor < Cor_graph + 0.1, ]
  # Color:
  DF_graph$color <- paste0(DF_graph$Método, " (Cor = ", DF_graph$Cor, ")") # Opción 1
  # DF_graph$color <- paste0(DF_graph$Método, " (Cor = ", DF_graph$Cor,
  #                          " | Var = ", DF_graph$Var, ")") # Opción 2
  
  p <- ggplot(data = DF_graph, mapping = aes(x = DF_graph$delta, y = DF_graph$Potencia,
                                             color = color)) +
    geom_line(size = .6, linetype = "solid") + geom_point() +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
         subtitle = paste0("Comparación MANTA - MANOVA", "\n\n",
                           "Tipo de datos = ", unique(DF_graph$Datos),
                           " | Cor = ", unique(DF_graph$Cor),"\n\n",
                           "[modelo = ", DF_graph$Modelo, " | alpha = ", DF_graph$alpha,
                           " | S = ", DF_graph$S, " | n = ", DF_graph$n, " | a = ",  a,
                           " | b = ",  b, " | q = ", DF_graph$q, " | Var = ", DF_graph$Var,
                           " | u = ",  u, " | w = ",  w, " | loc = ", DF_graph$loc, " | lambda = ", DF_graph$lambda,
                           "]                     "),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.position = c(0.09, 0.72),
          legend.background = element_rect(fill = alpha('lightgray', 0.75))) +
    guides(color = guide_legend(title = "Método (Cor)"))
  
  return(list(p))
  
}


# Función para crear las gráficas "∆ vs. Potencia" para un solo método, comparando
# los resultados para las diferentes transformaciones, bajo la combinación de variables
# "S - α - Var - Cor" dada.
# Agrupación de "color" ≡ Tipo de datos
DeltaPotGraph_DataTypeColor_mvnorm <- function(DF_Results, Metodo_graph) {
  
  # Genera la gráfica "∆ vs. Potencia" para los diferentes Y considerados en el método dado.
  
  DF_graph <- DF_Results[DF_Results$S == S_graph
                         & DF_Results$alpha == alpha_graph
                         & DF_Results$Var == Var_graph
                         & DF_Results$Cor > Cor_graph - 0.1 & DF_Results$Cor < Cor_graph + 0.1, ]
  # Color:
  DF_graph$color <- paste0(DF_graph$Método, " (", DF_graph$Datos, ")")
  
  p <- ggplot(data = DF_graph, mapping = aes(x = DF_graph$delta, y = DF_graph$Potencia,
                                             color = color)) +
    geom_line(size = .6, linetype = "solid") + geom_point() +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
         subtitle = paste0("Estudio de la posible invarianza a la transformación de datos de ",
                           Metodo_graph, "\n\n",
                           "Cor = ", unique(DF_graph$Cor),"\n\n",
                           "[modelo = ", DF_graph$Modelo, " | alpha = ", DF_graph$alpha,
                           " | S = ", DF_graph$S, " | n = ", DF_graph$n, " | a = ",  a,
                           " | b = ",  b, " | u = ",  u, " | w = ",  w,
                           " | Var = ", DF_graph$Var, " | q = ", DF_graph$q,
                           " | loc = ", DF_graph$loc, " | stdev = ", unique(DF_graph$stdev),
                           " | lambda = ", DF_graph$lambda, "]                     "),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.position = c(0.09, 0.72),
          legend.background = element_rect(fill = alpha('lightgray', 0.75))) +
    guides(color = guide_legend(title = "Método (Tipo de datos)"))
  
  return(list(p))
  
}


# Funciones para estudiar la similitud de resultados entre MANTA y MANOVA teniendo en cuenta
# los resultados para las diferentes transformaciones, bajo la combinación de variables
# "S - α - Var - Cor" dada.

## Boxplots comparativos "∆ vs. Potencia" .𝜏.:

### "facet_wrap(~ Cor)" para "Tipo de datos - S -α - Var" fijo, y agrupados por "método".
DeltaPotBoxPlot_Method_mvnorm_1 <- function(DF_Results) {
  
  # Genera los boxplots comparativos "∆ vs. Potencia" entre MANTA y MANOVA mostrando
  # los resultados según el "Cor" utilizado, manteniendo "Tipo de datos - S - α - Var"
  # fijos, y agrupando por "Método".
  
  Cor.labs <- list("0" = "0.0",
                   "0.2" = "0.2",
                   "0.4" = "0.4",
                   "0.6" = "0.6",
                   "0.8" = "0.8")
  # levels(DF_Results$Cor) <- levels(factor(c(sprintf("%1.1f", DF_Results$Cor))))
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, fill = Método)) +
    geom_boxplot() +
    facet_wrap(~ Cor, labeller = labeller(Cor.labs)) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
         subtitle = paste0("Comparación MANTA - MANOVA", "\n\n",
                           "Tipo de datos = ", DF_Results$Datos, " | ",
                           paste0(expression(TeX(r'($\alpha$)')), sprintf(" = %1.3f", DF_Results$alpha)), "\n\n",
                           "[modelo = ", DF_Results$Modelo, " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w, " | Var = ", DF_Results$Var,
                           " | q = ", DF_Results$q, " | loc = ", DF_Results$loc, " | lambda = ", DF_Results$lambda,
                           "]                     "),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

### "facet_grid(alpha ~ Datos)" para "S - Var - Cor" fijo, y agrupados por "método".
DeltaPotBoxPlot_Method_mvnorm_2 <- function(DF_Results) {
  
  # Genera los boxplots comparativos "∆ vs. Potencia" entre MANTA y MANOVA mostrando
  # los resultados según el "tipo de datos" considerado y el "α" utilizado, manteniendo
  # "Var" fijo, y agrupando por "Método".
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, fill = Método)) +
    geom_boxplot() +
    facet_grid((glue('alpha*" = {alpha}"') ~ glue("'{Datos}'")), labeller = label_parsed) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
         subtitle = paste0("Comparación MANTA - MANOVA", "\n\n",
                           "Cor = ", sprintf("%1.1f", DF_Results$Cor), "\n\n",
                           "[modelo = ", DF_Results$Modelo, " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w, " | Var = ", DF_Results$Var,
                           " | q = ", DF_Results$q, " | loc = ", DF_Results$loc, " | lambda = ", DF_Results$lambda,
                           "]                     "),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

# Como los Boxplot no son muy útiles para estudiar el estadístico p-valor de interés, y su relación
# con la potencia del método dado, se usarán mejor gráficos de puntos ("scatter plots"): 

## Scatter plots comparativos "∆ vs. Potencia" .𝜏.:

### "facet_grid(Var ~ Cor)" para "Tipo de datos - S - α" fijo, y agrupados por "método".
DeltaPotScatPlot_Method_mvnorm_1 <- function(DF_Results) {
  
  # Genera los scatter plots comparativos "∆ vs. Potencia" entre MANTA y MANOVA mostrando
  # los resultados en una matriz de gráficos con el tipo de varianza "Var" en las filas,
  # y el valor de "Cor" utilizado en las columnas. Se ha mantenido el "Tipo de datos - S
  # - α" fijos.
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, color = Método)) +
    geom_line(size = .4, linetype = "solid") + geom_point() +
    facet_grid(Var ~ Cor, labeller = label_both) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia (", r'($\alpha$)', " = ",
                            sprintf("%1.3f", DF_Results$alpha), ")"), italic = TRUE),
         subtitle = paste0("Comparación MANTA - MANOVA", "\n",
                           "Tipo de datos = ", DF_Results$Datos, "\n",
                           "modelo = ", DF_Results$Modelo, " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

### "facet_grid(Var ~ Datos)" para "S - α" fijo, y agrupados por "método", usando una 
### matriz de correlación ("corr") aleatoria => Cor no relevante.
DeltaPotScatPlot_Method_mvnorm_2 <- function(DF_Results) {
  
  # Genera los scatter plots comparativos "∆ vs. Potencia" entre MANTA y MANOVA mostrando
  # los resultados en una matriz de gráficos con el tipo de varianza "Var" en las filas,
  # y el "Tipo de datos" utilizado en las columnas. Se ha mantenido "S - α" fijos, y 
  # utilizado una matriz de correlación ("corr") aleatoria.
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, color = Método)) +
    geom_line(size = .4, linetype = "solid") + geom_point() +
    facet_grid(Var ~ Datos, labeller = label_both) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia (", r'($\alpha$)', " = ",
                            sprintf("%1.3f", DF_Results$alpha), ")"), italic = TRUE),
         subtitle = paste0("Comparación MANTA - MANOVA (Matriz de correlación aleatoria)", "\n",
                           "modelo = ", DF_Results$Modelo, " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

### "facet_grid(Var ~ Cor)" para "Método - S - α" fijo, y agrupados por "Tipo de datos".
DeltaPotScatPlot_Method_mvnorm_3 <- function(DF_Results, Graph_Method) {
  
  # Genera los scatter plots comparativos "∆ vs. Potencia" del método escogido, mostrando
  # los resultados en una matriz de gráficos con el tipo de varianza "Var" en las filas,
  # y el valor de "Cor" utilizado en las columnas. Se ha mantenido el "Método - S - α" fijos.
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, color = Datos)) +
    geom_line(size = .4, linetype = "solid") + geom_point() +
    facet_grid(Var ~ Cor, labeller = label_both) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia (", r'($\alpha$)', " = ",
                            sprintf("%1.3f", DF_Results$alpha), ")"), italic = TRUE),
         subtitle = paste0("Estudio de la posible invarianza a la transformación de datos de ",
                           Graph_Method, "\n",
                           "modelo = ", DF_Results$Modelo, " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

### "facet_grid(Var ~ Método)" para "S - α" fijo, y agrupados por "Tipo de datos",
### usando una matriz de correlación ("corr") aleatoria => Cor no relevante.
DeltaPotScatPlot_Method_mvnorm_4 <- function(DF_Results) {
  
  # Genera los scatter plots comparativos "∆ vs. Potencia" mostrando los resultados
  # en una matriz de gráficos con el tipo de varianza "Var" en las filas, y el método
  # utilizado en las columnas. Se ha mantenido el "S - α" fijos.
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, color = Datos)) +
    geom_line(size = .4, linetype = "solid") + geom_point() +
    facet_grid(Var ~ Método, labeller = label_both) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia (", r'($\alpha$)', " = ",
                            sprintf("%1.3f", DF_Results$alpha), ")"), italic = TRUE),
         subtitle = paste0("Estudio de la posible invarianza a la transformación ",
                           "de datos de los métodos MANTA y MANOVA", "\n",
                           "modelo = ", DF_Results$Modelo, " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

############################## 
# Funciones gráficas simplex #
##############################

# Función para crear las gráficas "∆ vs. Potencia" para un solo método, comparando
# los resultados para las diferentes transformaciones, bajo la combinación de variables
# "S - α - Var - q - loc - stdev" dada.
# Agrupación de "color" ≡ Tipo de datos
DeltaPotGraph_DataTypeColor_simplex <- function(DF_Results, Metodo_graph) {
  
  # Genera la gráfica "∆ vs. Potencia" para los diferentes Y considerados en el método dado.
  
  DF_graph <- DF_Results[DF_Results$S == S_graph
                         & DF_Results$alpha == alpha_graph
                         & DF_Results$Var == Var_graph
                         & DF_Results$q == q_graph
                         & DF_Results$loc == loc_graph, ]
  # Color:
  DF_graph$color <- paste0(DF_graph$Método, " (", DF_graph$Datos, ")")
  
  p <- ggplot(data = DF_graph, mapping = aes(x = DF_graph$delta, y = DF_graph$Potencia,
                                             color = color)) +
    geom_line(size = .6, linetype = "solid") + geom_point() +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia"), italic = TRUE),
         subtitle = paste0("Estudio de la posible invarianza a la transformación de datos de ",
                           Metodo_graph, "\n\n",
                           "q = ", unique(DF_graph$q), " | loc = ", unique(DF_graph$loc),
                           " | stdev = ", unique(DF_graph$stdev), "\n\n",
                           "[modelo = ", DF_graph$Modelo, " | alpha = ", DF_graph$alpha,
                           " | S = ", DF_graph$S, " | n = ", DF_graph$n, " | a = ",  a,
                           " | b = ",  b, " | u = ",  u, " | w = ",  w,
                           " | Var = ", DF_graph$Var, " | Cor = ", DF_graph$Cor,
                           " | lambda = ", DF_graph$lambda, "]                     "),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.position = c(0.09, 0.72),
          legend.background = element_rect(fill = alpha('lightgray', 0.75))) +
    guides(color = guide_legend(title = "Método (Tipo de datos)"))
  
  return(list(p))
  
}

## Scatter plots comparativos "∆ vs. Potencia" .𝜏.:

### "facet_grid(q ~ loc)" para una "stdev" dada, manteniendo "S - α" fijos.
DeltaPotScatPlot_Method_simplex_1 <- function(DF_Results) {
  
  # Genera los scatter plots comparativos "∆ vs. Potencia" para el método MANTA, teniendo
  # en cuenta las diferentes transformaciones aplicadas, y mostrando los resultados en una
  # matriz de gráficos con el número de respuestas "q" en las filas, y "loc" en las columnas.
  # Se ha mantenido la combinación "S - α" fija.
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, color = Datos)) +
    geom_line(size = .4, linetype = "solid") + geom_point() +
    facet_grid(q ~ loc, labeller = label_both) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia (", r'($\alpha$)', " = ",
                            sprintf("%1.3f", DF_Results$alpha), ")"), italic = TRUE),
         subtitle = paste0("Estudio de la posible invarianza a la transformación de datos de MANTA",
                           "\n", "modelo = ", DF_Results$Modelo, " | pdist = ", pdist,
                           " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}

### "facet_grid(q ~ loc)" para una "stdev" dada, manteniendo "S - α" fijos, haciendo zoom
### en una zona específica del eje de abcisas "∆".
DeltaPotScatPlot_Method_simplex_1_ZoomDelta <- function(DF_Results, DeltaMin = NA, DeltaMax = NA,
                                                        PotMin = NA, PotMax = NA) {
  
  # Genera los scatter plots comparativos "∆ vs. Potencia" para el método MANTA, teniendo
  # en cuenta las diferentes transformaciones aplicadas, y mostrando los resultados en una
  # matriz de gráficos con el número de respuestas "q" en las filas, y "loc" en las columnas.
  # Se ha mantenido la combinación "S - α" fija.
  
  p <- ggplot(DF_Results, aes(x = delta, y = Potencia, color = Datos)) +
    geom_line(size = .4, linetype = "solid") + geom_point() +
    facet_grid(q ~ loc, labeller = label_both) +
    scale_x_continuous(limits = c(DeltaMin, DeltaMax)) +
    scale_y_continuous(limits = c(PotMin, PotMax)) +
    labs(title = TeX(paste0(r'($\Delta$)', " vs. Potencia (", r'($\alpha$)', " = ",
                            sprintf("%1.3f", DF_Results$alpha), ")"), italic = TRUE),
         subtitle = paste0("Estudio de la posible invarianza a la transformación de datos de MANTA",
                           "\n", "modelo = ", DF_Results$Modelo, " | pdist = ", pdist,
                           " | S = ", DF_Results$S, " | n = ", DF_Results$n,
                           " | a = ",  a, " | b = ",  b, " | u = ",  u, " | w = ",  w),
         x = TeX(paste0(" Parámetro de generación de la H1 (", r'($\Delta$)', ") "), italic = TRUE),
         y = TeX(paste0("Potencia del método"), italic = TRUE)) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 6),
          legend.background = element_rect(fill = alpha('lightgray', 0.75)))
  
  return(list(p))
  
}







########## 
# BACKUP #
##########
# -------


# # Sim.simplex
# pdist <- p_dist
# x <- c(loc, rep(1, q-1))
# y0 <- x/sum(x)
# 
# if (check) {
#   e <- c(1, rep(0, q-1))
#   y <- step2h1(y0, e, delta)
#   Y <- sim.simplex(q, n, y, stdev*hk, pdist)
#   return(data.frame(exp = y, obs = colMeans(Y)))
# }
# 
# # sim.simplex
# p0 <- y0
# tol <- 1e-10
# Y <- NULL
# 
# if (is.null(Y)) {
#   Y <- t(matrix(p0, nrow = q, ncol = n))
# }
# 
# elist <- list()
# for (i in 1:q) {
#   elist[[i]] <- rep(0, q)
#   elist[[i]][i] <- 1 
# }
# 
# for (i in 1:n) {
#   
#   if (pdist == "norm") {
#     steps <- rnorm(q, mean = 0, sd = stdev) 
#   } else if (pdist == "gamma") {
#     steps <- (rgamma(q, shape = 1, scale = 100) - 1)/gammasd(1, 100)*stdev
#   } else if (pdist == "beta") {
#     steps <- (rbeta(q, shape1 = 0.5, shape2 = 0.5) - 0.5)/betasd(0.5, 0.5)*stdev
#   } else {
#     stop (sprintf("Unknown pdist '%s'.", pdist))
#   }
#   
#   elist <- elist[sample(1:q)]
#   for (j in 1:q) {
#     k <- which(elist[[j]] == 1)
#     if (steps[j] < -Y[i, k]) {
#       steps[j] <- -Y[i, k] + tol
#       warning("One observation out of the simplex was corrected.")
#     } else if (steps[j] > (1 - Y[i,k]) ) {
#       steps[j] <- (1 - Y[i, k]) - tol
#       warning("One observation out of the simplex was corrected.")
#     }
#     Y[i, ] <- step2h1(Y[i, ], elist[[j]], steps[j])
#   }
# }
# 
# # Sim.simplex
# if (delta != 0) {
#   
#   # Generate e and H1 (+/- delta)
#   e <- c(1, rep(0, q-1))
#   y <- step2h1(y0, e, delta)
#   yprime <- step2h1(y0, e, -delta)
#   
#   if (any(grepl(":", B))) { # Then we are changing the interaction
#     levs <- levels(B)
#     b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]   # Recover b
#     B <- mapvalues(B, from = levs, to = 1:length(levs))             # Relevel
#     Y[B == (1+b), ] <- sim.simplex(q, sum(B == (1+b)), yprime, stdev, pdist)
#     Y[B == (2+b), ] <- sim.simplex(q, sum(B == (2+b)), y, stdev, pdist) 
#   } 
#   
#   Y[B == 1, ] <- sim.simplex(q, sum(B == 1), y, stdev*hk, pdist)   
#   Y[B == 2, ] <- sim.simplex(q, sum(B == 2), yprime, stdev, pdist)
#   
# } else {
#   
#   if (any(grepl(":", B))) { # Then we are changing the interaction
#     levs <- levels(B)
#     b <- as.numeric(unlist(strsplit(levs[length(levs)], ":")))[2]  # Recover b
#     B <- mapvalues(B, from = levs, to = 1:length(levs))            # Relevel
#   } 
#   
#   Y[B == 1, ] <- sim.simplex(q, sum(B == 1), y0 , stdev*hk, pdist)
# }


# CompMantaManovafunc_mvnorm <- function() {
#   
#   # # # # # # # # # # # # # # # # # # # # # # # # # #
#   #  Instalación y carga de los paquetes necesarios #
#   # # # # # # # # # # # # # # # # # # # # # # # # # #
#   
#   # Se cargan los paquetes, instalándolos si es necesario.
#   
#   #    require() => a warning if a package is not installed and then continue to execute the code
#   #    library() => an error and stop the execution of the code
#   
#   if (!require("devtools")) install.packages("devtools")
#   require(devtools)
#   library(devtools)
#   
#   if (!require("versions")) install.packages("versions")
#   require(versions)
#   library(versions)
#   
#   if (!require("manta")) devtools::install_github("dgarrimar/manta")
#   require(manta)
#   library(manta)
#   
#   package_list <- c("CompQuadForm", "car", "Hmisc", "MASS", "copula", "vegan",
#                     "randomcoloR", "plyr", "svMisc", "reshape2", "hms", "progress",
#                     "ggplot2", "arm", "pheatmap", "latex2exp")
#   
#   # Hmisc => Necesario para la función label
#   # randomcoloR => La función randomcoloR() genera un color aleatóriamente
#   # plyr => Necesario para la función mapvalues()
#   # svMisc => Permite usar progress() para indicar visualmente el progreso de la simulación
#   # reshape2 => La función melt() permite convertir un objeto "dist" en "DF"
#   # hms => La función as_hms() permite convertir "difftime" en DD:HH:MM:SS
#   # progress => Para mostrar una barra de progreso de la simulación
#   # latex2exp => Toma una cadena LaTeX, la analiza y devuelve la expresión plotmath más cercana
#   
#   for (i in package_list){
#     if(! i %in% installed.packages()){
#       install.packages(i)
#     }
#     require(i, character.only = TRUE)
#     library(i, character.only = TRUE)
#   }
#   
#   # # # # # # # # # # # # # # # # # #
#   # DF de los parámetros del modelo #
#   # # # # # # # # # # # # # # # # # #
#   
#   # Estos son: q, a, b, n, u, m, c, v, w, S, plotheatmap, k, D, l, s, p, p_dist, H, d, t, x, f.
#   
#   # Dependent variables: the parameter q (numerical type) is equivalent to the number of dependent variables.
#   # (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 4)
#   
#   # Factor A Levels: the parameter a (numerical type) is equivalent to the number of levels of factor A.
#   # (Value by default: 2 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
#   
#   # Factor B Levels: the parameter b (numerical type) is equivalent to the number of levels of factor B."
#   # (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
#   
#   # Number of samples: the parameter n (numerical type) is equivalent to the number of samples.
#   # (Value by default: 300 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 400)
#   
#   # Unbalance 1st B Level: the parameter u (numerical type) is equivalent to the unbalance degree of the first level regarding the rest of the levels of the factor B.
#   # (Value by default: 1 / Possible values: all positive rational numbers / Max. recommended to avoid an excess of interactions: < 10)
#   
#   # Model H0/H1: the parameter m (character type) is equivalent to the model that generates H0/H1.
#   # (Value by default: mvnorm / Possible values: mvnorm, simplex or multinom)
#   m <- modelSim # Equivalencia
#   
#   # Y correlation: the parameter c (numerical type) is equivalent to the correlation of the Y variables.
#   # (Value by default: 0 / Possible values: all rational numbers between -1.0 and 1.0)
#   c <- Cor # Equivalencia
#   
#   # Y variance: the parameter v (character type) is equivalent to the variance of the Y variables.
#   # (Value by default: equal / Possible values: equal or unequal)
#   v <- Var # Equivalencia
#   
#   # Changing Factor: the parameter w (character type) determines which factor changes.
#   # (Value by default: B / Possible values: A, B or AB)
#   
#   # Number of simulations: the parameter S (numerical type) is equivalent to the number of simulations.
#   # (Value by default: 1E3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 1E6)
#   
#   # Factors heatmap: the parameter plot (character type) determine whether or not to show the heatmap of the factor distribution.
#   # (Value by default: F / Possible values: F or T)
#   
#   # Chunk number or k: the parameter k or Chunk number (numeric type) determine whether to compute tstats or ado.
#   # (Value by default: 0 / Possible values: 0 (ado) or not 0 (tstats))
#   chunk <- k # Equivalencia
#   
#   # D, dd or DistDef: Multivariate non-normal distribution definition.
#   # (Value by default: unif-0-1 / Possible values: unif-0-1)
#   DistDef <- dd <- D # Equivalencia
#   
#   # l or lambda: lambda parameter (Poisson distribution) to generate size for 'multinom' generator model.
#   # (Value by default: 1000 / Possible values: 1000)
#   Lambda <- l # Equivalencia
#   
#   # s or stdev: stdev for the 'simplex' generator model.
#   # (Value by default: 0.1 / Possible values: all rational numbers >0)
#   stdev <- s # Equivalencia
#   
#   # p, loc or position: location of the 'simplex' generator model.
#   # (Value by default: 1 / Possible values: ?)
#   position <- loc <- p # Equivalencia
#   
#   # p_dist: Distribution for 'simplex' generator model.
#   # (Value by default: norm / Possible values: norm, gamma or beta)
#   p_dist <- pdist # Equivalencia
#   
#   # hk, H or heterosk: Heteoskedasticity degree (level 1).
#   # (Value by default: 1 / Possible values: 1)
#   hk = as.numeric(1)
#   heterosk <- H <- hk # Equivalencia
#   
#   # d or delta: H1 generation parameter.
#   # (Value by default: 0 / Possible values: 0 for H0, > 0 for H1)
#   d <- delta # Equivalencia
#   
#   # t or transf: Data transformation type.
#   # (Value by default: none / Possible values: none, sqrt, log,...)
#   # t_values = c("none", "sqrt", "log", "...")
#   # transf_values <- t # Equivalencia
#   
#   # o or output: Output file name.
#   # (Value by default: outputfiletemp / Possible values: NULL or outputfiletemp)
#   # o = "outputfiletemp"
#   
#   # x or cpu: Number of cores.
#   # (Value by default: 10 / Possible values: ?)
#   cores <- x # Equivalencia
#   
#   # f or fx: Path to helper functions.
#   # f = as.character(getwd())
#   fx <- f # Equivalencia
#   
#   # Se crea una lista con estos parámetros y se asignan al entorno global.
#   Parameter_List <- list("q" = q, "a" = a, "b" = b, "n" = n,
#                          "u" = u, "m" = m, "c" = c, "v" = v,
#                          "w" = w, "S" = S, "plotheatmap" = plotheatmap,
#                          "k" = k, "D" = D, "l" = l, "s" = s, "p" = p, "p_dist" = p_dist,
#                          "H" = H, "d" = d) # "t" = t, "o" = 0, "x" = x, "f" = f)
#   list2env(Parameter_List, .GlobalEnv)
#   
#   Parameter_defin <- c("Dependent variables",
#                        "Factor A Levels",
#                        "Factor B Levels",
#                        "Number of samples",
#                        "Unbalance 1st B Level",
#                        "Model H0/H1",
#                        "Y correlation",
#                        "Y variance",
#                        "Changing Factor",
#                        "Number of simulations",
#                        "Factors heatmap",
#                        "Chunk number",
#                        "Multivariate non-normal distribution definition",
#                        "Lambda parameter of the Poisson distribution",
#                        "Stdev for the 'simplex' generator model",
#                        "Location of the 'simplex' generator model",
#                        "Distribution for 'simplex' generator model: norm, gamma or beta",
#                        "Heteroskedasticity degree (level 1)",
#                        "Delta (H1 generation parameter)")
#   # "Data transformation: sqrt, log or none",
#   # "Output file name",
#   # "Number of cores",
#   # "Path to helper functions")
#   
#   Parameter_names <- c("q",
#                        "a",
#                        "b",
#                        "n",
#                        "u",
#                        "m / modelSim",
#                        "c / Cor",
#                        "v / Var",
#                        "w",
#                        "S",
#                        "plotheatmap",
#                        "k / chunk",
#                        "D / DistDef / dd",
#                        "l / Lambda",
#                        "s /stdev",
#                        "p / position / loc",
#                        "p_dist / pdist",
#                        "H / hk",
#                        "d / delta")
#   # "t / transf",
#   # "Output file name",
#   # "x / cores",
#   # "f / fx")
#   
#   Parameter_values <- c()
#   for (i in 1:length((Parameter_List))){
#     Parameter_values <- append(Parameter_values, toString(Parameter_List[[i]]))}
#   
#   Parameter_class <- c()
#   for (i in 1:length((Parameter_List))){
#     Parameter_class <- append(Parameter_class, class(Parameter_List[[i]]))}
#   
#   Parameter_DF = data.frame(Parameter_defin, Parameter_names, Parameter_values, Parameter_class)
#   colnames(Parameter_DF) = c("Parameter definition", "Parameter names", "Parameter value", "Parameter class type")
#   
#   # Se asignan los DF generados al entorno global.
#   Parameter_DF_List <- list("Parameter_DF" = Parameter_DF, "Parameter_values" = Parameter_values,
#                             "Parameter_names" = Parameter_names,"Parameter_class" = Parameter_class,
#                             "Parameter definition" = Parameter_defin)
#   list2env(Parameter_DF_List, .GlobalEnv)
#   
#   # # # # # # # # # # # # # # # # # # # # # # # #
#   #   Niveles de los factores y su interacción  #
#   # # # # # # # # # # # # # # # # # # # # # # # #
#   
#   # Se obtienen los niveles de los factores y su interacción según los parámetros: n, a, b, u, w. 
#   # Especificando si se ha de mostrar el mapa de calor con el valor lógico "plotheatmap".
#   
#   label <- function(a, b, n, u, w, plot = plotheatmap, seed = 1325) {
#     
#     set.seed(seed)
#     
#     if (w == "A") {
#       ua <- u
#       ub <- 1
#     } else if (w == "B") {
#       ua <- 1
#       ub <- u
#     } else if (w == "AB") {
#       ua <- ub <- u
#     }
#     
#     xa <- c(ua, rep(1, a-1))
#     pa <- round(xa/sum(xa)*n)
#     r <- n - sum(pa)
#     sel <- sample(1:length(pa), size = abs(r))
#     if (r > 0) {
#       pa[sel] <- pa[sel] + 1
#     } else if (r < 0) {
#       pa[sel] <- pa[sel] - 1
#     }
#     A <- rep(1:a, times = pa)
#     # A <- gl(a, n/a, length = n)
#     
#     B <- c()
#     xb <- c(ub, rep(1, b-1))
#     for (i in table(A)) {
#       pb <- round(xb/sum(xb)*i)
#       r <- i - sum(pb)
#       sel <- sample(1:length(pb), size = abs(r))
#       if (r > 0) {
#         pb[sel] <- pb[sel] + 1
#       } else if (r < 0) {
#         pb[sel] <- pb[sel] - 1
#       } 
#       B <- c(B, rep(1:b, times = pb))
#     }
#     if (plotheatmap) {
#       pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
#     }
#     return(list(factor(A), factor(B)))
#   }
#   
#   labs <- label(a, b, n, u, w)
#   A <- labs[[1]]
#   B <- labs[[2]]
#   
#   if (w == "A") {
#     ch <- A
#   } else if (w == "B") {
#     ch <- B
#   } else if (w == "AB") {
#     ch <- A:B
#   } else {
#     stop(sprintf("Unknown factor: '%s'.", w))
#   }
#   
#   # Se asignan al entorno los parámetros generados en esta sección:
#   Parameter_gen2_List <- list("ch" = ch, "A" = A, "B" = B)
#   list2env(Parameter_gen2_List, .GlobalEnv)
#   
#   # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   # #  Asignación de la "stdev" para el modelo "simplex"  #
#   # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   # 
#   # # Para: m = modelSim = "simplex" 
#   # # Se carga la tabla correspondiente de "stdev":
#   # if (modelSim == "simplex") {
#   #   
#   #   tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", f, pdist), h = T)
#   #   colnames(tbl) <- c("Q", "L", "S")
#   #   tbl_name <- "simplex_table_Q_L_S"
#   #   write.csv(tbl, file = paste(tbl_name, ".csv", sep = ""), row.names = FALSE)
#   #   
#   #   if (! q %in% unique(tbl$Q)) {
#   #     stop(sprintf("stdev not precomputed for q = %s", q))} 
#   #   
#   #   # Se busca la "stdev" que corresponde a los valores [q, loc] seleccionados.
#   #   stdev <- subset(tbl, Q == q & L == loc)$S
#   #   
#   #   # Se crea y guarda un DF con los valores [Q, L, S] que cumplen la condición,
#   #   # asignándole un nombre que refleja los parámetros seleccionados:
#   #   stdevDF <- data.frame(subset(tbl, Q == q & L == loc)$Q,
#   #                         subset(tbl, Q == q & L == loc)$L,
#   #                         subset(tbl, Q == q & L == loc)$S)
#   #   colnames(stdevDF) <- c("q (Q)", "loc (L)", "stdev (S)")
#   #   stdev_name = paste("stdev_", modelSim, sep = "")
#   #   assign(stdev_name, stdevDF[[3]], envir = .GlobalEnv)
#   #   list2env(stdevDF, .GlobalEnv)
#   #   
#   #   # Se asignan los parámetros generados en esta sección:
#   #   Parameter_gen3_List <- list("tbl" = tbl)
#   #   list2env(Parameter_gen3_List, .GlobalEnv)}
#   
#   # # # # # # # # # # # # # # # # # # # #
#   #  Simulación del conjunto de datos   #
#   #  y aplicación de MANTA y MANOVA     #
#   # # # # # # # # # # # # # # # # # # # #
#   
#   # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
#   # las p del factor B para cada método y conjunto de datos considerado:
#   Manta_pB = rep(NA, S)
#   Manova_pB = rep(NA, S)
#   Manta_Ylog_pB = rep(NA, S)
#   Manova_Ylog_pB = rep(NA, S)
#   Manta_Ysqrt_pB = rep(NA, S)
#   Manova_Ysqrt_pB = rep(NA, S)
#   Manta_Yscaled_pB = rep(NA, S)
#   Manova_Yscaled_pB = rep(NA, S)
#   Manta_YnormMinMax_pB = rep(NA, S)
#   Manova_YnormMinMax_pB = rep(NA, S)
#   
#   # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
#   # lostiempos de simulación para cada método y conjunto de datos considerado: 
#   t_Manta_pB = rep(NA, S)
#   t_Manova_pB = rep(NA, S)
#   t_Manta_Ylog_pB = rep(NA, S)
#   t_Manova_Ylog_pB = rep(NA, S)
#   t_Manta_Ysqrt_pB = rep(NA, S)
#   t_Manova_Ysqrt_pB = rep(NA, S)
#   t_Manta_Yscaled_pB = rep(NA, S)
#   t_Manova_Yscaled_pB = rep(NA, S)
#   t_Manta_YnormMinMax_pB = rep(NA, S)
#   t_Manova_YnormMinMax_pB = rep(NA, S)
#   
#   t_start_simulation <- Sys.time()
#   
#   writeLines(paste("\n", "La simulación de datos para la combinación de variables [alpha = ", Alpha, ", delta = ",
#                    delta, ", Cor = ", Cor, "] se ha iniciado:", format(t_start_simulation), "\n"))
#   
#   pb <- txtProgressBar(min = 0, max = S, style = 3, width = 50, char = "=") 
#   
#   for (i in 1:S) {
#     # if (modelSim == "simplex") {
#     #   if (i == 1) { # Sanity check
#     #     tol <- 0.01
#     #     check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, pdist)
#     #     wm <- which.max(check$exp)
#     #     if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
#     #       stop("Deviation from expected centroid greater than tolerance.")
#     #     }
#     #   }
#     #   Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = F, pdist)
#     # } else 
#     if (modelSim == "mvnorm") {
#       Y <- Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor) # Cambiar 0 <-> núm grande (para transformaciones)
#     } else if (modelSim == "multinom") {
#       Y <- Sim.multinom(ch, q, n, Lambda, delta, loc)
#     }
#     
#     # Transformaciones posibles de los datos simulados (Y):
#     Ylog <- log(Y+1) # Transformación logarítmica. Se suma "1" para evitar "log(0) = -Inf" y valores < 0.
#     Ysqrt <- sqrt(Y) # Transformación sqrt.
#     Yscaled <- cbind(Y[,1]/sd(Y[,1]), Y[,2]/sd(Y[,2]), Y[,3]/sd(Y[,3])) # Escalado de datos: se divide cada valor por la desviación 
#     # típica del conjunto de datos.
#     YnormMinMax <-  cbind((Y[,1]-min(Y[,1]))/(max(Y[,1])-min(Y[,1])), # Normalización Min- Max: si queremos que cada variable 
#                           # contribuya por igual al análisis.
#                           (Y[,2]-min(Y[,2]))/(max(Y[,2])-min(Y[,2])),
#                           (Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3])))
#     
#     # Verificación de NAs en los conjuntos de datos creados:
#     if (sum(is.na(Y)) != 0) {
#       cat("\nEl conjunto de datos Y contiene al menos un valor NA.\n\n")
#       #next # O "stop".
#     }
#     if (sum(is.na(Ylog)) != 0) {
#       cat("\nEl conjunto de datos Ylog contiene al menos un valor NA.\n\n")
#       #next # O "stop".
#     }
#     if (sum(is.na(Ysqrt)) != 0) {
#       cat("\nEl conjunto de datos Ysqrt contiene al menos un valor NA.\n\n")
#       #next # O "stop".
#     }
#     if (sum(is.na(Yscaled)) != 0) {
#       cat("\nEl conjunto de datos Yscaled contiene al menos un valor NA.\n\n")
#       #next # O "stop".
#     }
#     if (sum(is.na(YnormMinMax)) != 0) {
#       cat("\nEl conjunto de datos YnormMinMax contiene al menos un valor NA.\n\n")
#       #next # O "stop".
#     }
#     
#     # Aplicación de Manta y Manova para cada conjunto de datos (Y, Ylog, Ysqrt,...):
#     
#     # Y: MANOVA no funciona para "simplex".
#     # ==> Si hago "summary(manova(..., tol = 0))" SI FUNCIONA PERO MUESTRA NAs!!!
#     t_Sim_start <- Sys.time()
#     Manta_pB[i] <- manta(Y ~ A + B + A:B)$aov.tab[2,6]
#     t_Manta_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     t_Sim_start <- Sys.time()
#     Manova_pB[i] <- summary(manova(Y ~ A + B + A:B))$stats[2,6]
#     t_Manova_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     
#     # YLog:
#     t_Sim_start <- Sys.time()
#     Manta_Ylog_pB[i] <- manta(Ylog ~ A + B + A:B)$aov.tab[2,6]
#     t_Manta_Ylog_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     t_Sim_start <- Sys.time()
#     Manova_Ylog_pB[i] <- summary(manova(Ylog ~ A + B + A:B))$stats[2,6]
#     t_Manova_Ylog_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     
#     # Ysqrt:
#     t_Sim_start <- Sys.time()
#     Manta_Ysqrt_pB[i] <- manta(Ysqrt ~ A + B + A:B)$aov.tab[2,6]
#     t_Manta_Ysqrt_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     t_Sim_start <- Sys.time()
#     Manova_Ysqrt_pB[i] <- summary(manova(Ysqrt ~ A + B + A:B))$stats[2,6]
#     t_Manova_Ysqrt_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     
#     # Yscaled: MANOVA no funciona para "simplex".
#     # ==> Si hago "summary(manova(..., tol = 0))" SI FUNCIONA PERO MUESTRA NAs!!!
#     t_Sim_start <- Sys.time()
#     Manta_Yscaled_pB[i] <- manta(Yscaled ~ A + B + A:B)$aov.tab[2,6]
#     t_Manta_Yscaled_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     t_Sim_start <- Sys.time()
#     Manova_Yscaled_pB[i] <- summary(manova(Yscaled ~ A + B + A:B))$stats[2,6]
#     t_Manova_Yscaled_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     
#     # YnormMinMax: MANOVA no funciona para "simplex".
#     # ==> Si hago "summary(manova(..., tol = 0))" SI FUNCIONA PERO MUESTRA NAs!!!
#     t_Sim_start <- Sys.time()
#     Manta_YnormMinMax_pB[i] <- manta(YnormMinMax ~ A + B + A:B)$aov.tab[2,6]
#     t_Manta_YnormMinMax_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     t_Sim_start <- Sys.time()
#     Manova_YnormMinMax_pB[i] <- summary(manova(YnormMinMax ~ A + B + A:B))$stats[2,6]
#     t_Manova_YnormMinMax_pB[i] <- difftime(Sys.time(), t_Sim_start, units="auto")
#     
#     setTxtProgressBar(pb, i)
#   }
#   
#   close(pb)
#   
#   # Se calculan los tiempos de computación para cada método:
#   T_Manta_pB = sum(t_Manta_pB)
#   T_Manova_pB = sum(t_Manova_pB)
#   T_Manta_Ylog_pB = sum(t_Manta_Ylog_pB)
#   T_Manova_Ylog_pB = sum(t_Manova_Ylog_pB)
#   T_Manta_Ysqrt_pB = sum(t_Manta_Ysqrt_pB)
#   T_Manova_Ysqrt_pB = sum(t_Manova_Ysqrt_pB)
#   T_Manta_Yscaled_pB = sum(t_Manta_Yscaled_pB)
#   T_Manova_Yscaled_pB = sum(t_Manova_Yscaled_pB)
#   T_Manta_YnormMinMax_pB = sum(t_Manta_YnormMinMax_pB)
#   T_Manova_YnormMinMax_pB = sum(t_Manova_YnormMinMax_pB)
#   
#   # Se calcula la potencia de cada método para el α dado:
#   t_Sim_start <- Sys.time()
#   Manta_pot <- mean(Manta_pB < Alpha)
#   t_Manta_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_pB
#   t_Sim_start <- Sys.time()
#   Manova_pot <- mean(Manova_pB < Alpha)
#   t_Manova_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_pB
#   Y_Manta_Manova_pot <- cbind(Manta_pot, t_Manta_pot, Manova_pot, t_Manova_pot)
#   
#   t_Sim_start <- Sys.time()
#   Manta_Ylog_pot <- mean(Manta_Ylog_pB < Alpha)
#   t_Manta_Ylog_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Ylog_pB
#   t_Sim_start <- Sys.time()
#   Manova_Ylog_pot <- mean(Manova_Ylog_pB < Alpha)
#   t_Manova_Ylog_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Ylog_pB
#   Ylog_Manta_Manova_pot <- cbind(Manta_Ylog_pot, t_Manta_Ylog_pot, Manova_Ylog_pot, t_Manova_Ylog_pot)
#   
#   t_Sim_start <- Sys.time()
#   Manta_Ysqrt_pot <- mean(Manta_Ysqrt_pB < Alpha)
#   t_Manta_Ysqrt_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Ysqrt_pB
#   t_Sim_start <- Sys.time()
#   Manova_Ysqrt_pot <- mean(Manova_Ysqrt_pB < Alpha)
#   t_Manova_Ysqrt_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Ysqrt_pB
#   Ysqrt_Manta_Manova_pot <- cbind(Manta_Ysqrt_pot, t_Manta_Ysqrt_pot, Manova_Ysqrt_pot, t_Manova_Ysqrt_pot)
#   
#   t_Sim_start <- Sys.time()
#   Manta_Yscaled_pot <- mean(Manta_Yscaled_pB < Alpha)
#   t_Manta_Yscaled_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_Yscaled_pB
#   t_Sim_start <- Sys.time()
#   Manova_Yscaled_pot <- mean(Manova_Yscaled_pB < Alpha)
#   t_Manova_Yscaled_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_Yscaled_pB
#   Yscaled_Manta_Manova_pot <- cbind(Manta_Yscaled_pot, t_Manta_Yscaled_pot, Manova_Yscaled_pot, t_Manova_Yscaled_pot)
#   
#   t_Sim_start <- Sys.time()
#   Manta_YnormMinMax_pot <- mean(Manta_YnormMinMax_pB < Alpha)
#   t_Manta_YnormMinMax_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manta_YnormMinMax_pB
#   t_Sim_start <- Sys.time()
#   Manova_YnormMinMax_pot <- mean(Manova_YnormMinMax_pB < Alpha)
#   t_Manova_YnormMinMax_pot <- difftime(Sys.time(), t_Sim_start, units="auto") + T_Manova_YnormMinMax_pB
#   YnormMinMax_Manta_Manova_pot <- cbind(Manta_YnormMinMax_pot, t_Manta_YnormMinMax_pot, Manova_YnormMinMax_pot, t_Manova_YnormMinMax_pot)
#   
#   # DF1 que almacena la potencia de cada método:
#   DF_CompPot_res_Ytype <- c("Datos sin transformar", "Transformación logarítmica",
#                      "Transformación raíz cuadrada", "Escalado de datos por desviación típica",
#                      "Normalización Min- Max")
#   DF_CompPot_res_delta_iter <- rep(delta, length(DF_CompPot_res_Ytype))
#   DF_CompPot_res_Cor_iter <- rep(Cor, length(DF_CompPot_res_Ytype))
#   DF_CompPot_res_Alpha_iter <- rep(Alpha, length(DF_CompPot_res_Ytype))
#   DF_CompPot_res <- data.frame(DF_CompPot_res_Ytype, DF_CompPot_res_delta_iter, DF_CompPot_res_Cor_iter, DF_CompPot_res_Alpha_iter,
#                         rbind(Y_Manta_Manova_pot, Ylog_Manta_Manova_pot, Ysqrt_Manta_Manova_pot,
#                               Yscaled_Manta_Manova_pot, YnormMinMax_Manta_Manova_pot))
#   colnames(DF_CompPot_res) <- c("Tipo de datos", "Delta (d)", "Correlation (c)", "Nivel de significación (alpha)",
#                          "Potencia de MANTA", "t computación MANTA", "Potencia de MANOVA", "t computación MANOVA")
#   
#   # Se asignan los parámetros generados en esta sección:
#   DF_CompPot_res_List <- list("DF_CompPot_res" = DF_CompPot_res)
#   list2env(DF_CompPot_res_List, .GlobalEnv)
#   
#   # DF2 que almacena la potencia de cada método:
#   
#   DF_MANTA_Results <- cbind(rep(m, length(Manta_pot)), rep(Alpha, length(Manta_pot)), rep(Cor, length(Manta_pot)),
#                             rep(delta, length(Manta_pot)), rep("MANTA", length(Manta_pot)),
#                             Manta_pot, t_Manta_pot, Manta_Ylog_pot, t_Manta_Ylog_pot,
#                             Manta_Ysqrt_pot, t_Manta_Ysqrt_pot, Manta_Yscaled_pot, t_Manta_Yscaled_pot,
#                             Manta_YnormMinMax_pot, t_Manta_YnormMinMax_pot)
#   colnames(DF_MANTA_Results) <- c("m", "alpha", "Cor", "delta", "Método",
#                                   "PotY", "tPotY", "PotYLog", "tPotYLog", "PotYsqrt", "tPotYsqrt",
#                                   "PotYscaled", "tPotYscaled", "PotYnormMinMax", "tPotYnormMinMax")
#   
#   DF_MANOVA_Results <- cbind(rep(m, length(Manova_pot)), rep(Alpha, length(Manova_pot)), rep(Cor, length(Manova_pot)),
#                              rep(delta, length(Manova_pot)), rep("MANTA", length(Manova_pot)),
#                              Manova_pot, t_Manova_pot, Manova_Ylog_pot, t_Manova_Ylog_pot,
#                              Manova_Ysqrt_pot, t_Manova_Ysqrt_pot, Manova_Yscaled_pot, t_Manova_Yscaled_pot,
#                              Manova_YnormMinMax_pot, t_Manova_YnormMinMax_pot)
#   colnames(DF_MANOVA_Results) <- c("m", "alpha", "Cor", "delta", "Método",
#                                    "PotY", "tPotY", "PotYLog", "tPotYLog", "PotYsqrt", "tPotYsqrt",
#                                    "PotYscaled", "tPotYscaled", "PotYnormMinMax", "tPotYnormMinMax")
#   
#   # Se asignan los parámetros generados en esta sección:
#   DF2_pot_List <- list("DF_MANTA_Results" = DF_MANTA_Results, "DF_MANOVA_Results" = DF_MANOVA_Results)
#   list2env(DF2_pot_List, .GlobalEnv)
#   
# }