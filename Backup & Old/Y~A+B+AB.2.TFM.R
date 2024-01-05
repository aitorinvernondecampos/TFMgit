#!/usr/bin/env Rscript

## TFM - Máster universitario de Bioinformática y Bioestadística (UOC, UB)
## Assessing the properties of asymptotic PERMANOVA test through comprehensive simulations in the context of genetic studies
## Aitor Invernón de Campos

## Script for the evaluation of asymptotic PERMANOVA in complex models
## Model: Y ~ A + B + AB
## Based on Diego Garrido-Martín's script Y~A+B+AB.2.R (https://github.com/dgarrimar/manta-sim/blob/sim0/bin/Y~A%2BB%2BAB.2.R)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


## Modificación manual de los parámetros de entrada:

### Función que permite modificar los parámetros.

rm(list = ls(all.names = TRUE)) # Will clear all objects includes hidden objects.

Parameter_mod_DF <- function() {
  
  mod_param <- as.character(readline(writeLines(c("Do you want to select the parameters manually?"," ",
                                                "   Possible values: y or n", " "))))
  
  cat("\014") # Clean console
  
  Parameter_mod_vect_vals <- c() # Vector que almacena los valores de los parámetros
  Parameter_mod_vect_names <- c() # Vector que almacena los nombres de los parámetros
  Parameter_mod_vect_class <- c() # Vector que almacena la clase de cada parámetro
  
  Parameter_no_mod_vect_vals <- c() # Vector que almacena los valores de los parámetros no modificables
  Parameter_no_mod_vect_names <- c() # Vector que almacena los nombres de los parámetros no modificables
  Parameter_no_mod_vect_class <- c() # Vector que almacena la clase de cada parámetro no modificable
  
  if (mod_param == "y") {
    q <- as.numeric(readline(writeLines(c("The parameter q (numerical type) is equivalent to the number of dependent variables."," ",
                                          "   Value by default: 3", "   Possible values: all natural numbers",
                                          "   Max. recommended to avoid an excess of interactions: < 4", " ",
                                          "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, q)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Dependent variables")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(q))
  } else {
    q <- 3
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, q)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Dependent variables")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(q))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    a <- as.numeric(readline(writeLines(c("The parameter a (numerical type) is equivalent to the number of levels of factor A."," ",
                                          "   Value by default: 2", "   Possible values: all natural numbers",
                                          "   Max. recommended to avoid an excess of interactions: < 6", " ",
                                          "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, a)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Factor A Levels")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(a))
  } else {
    a <- 2
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, a)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Factor A Levels")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(a))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    b <- as.numeric(readline(writeLines(c("The parameter b (numerical type) is equivalent to the number of levels of factor B."," ",
                                          "   Value by default: 3", "   Possible values: all natural numbers",
                                          "   Max. recommended to avoid an excess of interactions: < 6", " ",
                                          "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, b)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Factor B Levels")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(a))
  } else {
    b <- 3
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, b)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Factor B Levels")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(a))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    n <- as.numeric(readline(writeLines(c("The parameter n (numerical type) is equivalent to the number of samples."," ",
                                          "   Value by default: 100", "   Possible values: all natural numbers",
                                          "   Max. recommended to avoid an excess of computational time: < 200", " ",
                                          "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, n)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Number of samples")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(n))
  } else {
    n <- 100
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, n)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Number of samples")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(n))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    u <- as.numeric(readline(writeLines(c("The parameter u (numerical type) is equivalent to the unbalance degree of the first level",
                                          "regarding the rest of the levels of the factor B.", " ",
                                          "   Value by default: 1", "   Possible values: all positive rational numbers",
                                          "   Max. recommended to avoid an excess of unbalancing: < 10", " ",
                                          "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, u)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Unbalance 1st B Level")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(u))
  } else {
    u <- 1
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, u)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Unbalance 1st B Level")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(u))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    m <- as.character(readline(writeLines(c("The parameter m (character type) is equivalent to the model that generates H0/H1."," ",
                                            "   Value by default: mvnorm", "   Possible values: mvnorm, simplex or multinom",
                                            " ", "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, m)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Model H0/H1")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(m))
  } else {
    m <- "mvnorm"
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, m)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Model H0/H1")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(m))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    c <- as.numeric(readline(writeLines(c("The parameter c (numerical type) is equivalent to the correlation of the Y variables."," ",
                                          "   Value by default: 0", "   Possible values: all rational numbers between -1.0 and 1.0",
                                          " ", "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, c)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Y correlation")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(c))
  } else {
    c <- 0
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, c)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Y correlation")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(c))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    v <- as.character(readline(writeLines(c("The parameter v (character type) is equivalent to the variance of the Y variables."," ",
                                            "   Value by default: equal", "   Possible values: equal or unequal",
                                            " ", "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, v)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Y variance")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(v))
  } else {
    v <- "equal"
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, v)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Y variance")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(v))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    w <- as.character(readline(writeLines(c("The parameter w (character type) determines which factor changes."," ",
                                            "   Value by default: B", "   Possible values: A, B or AB",
                                            " ", "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, w)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Changing Factor")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(w))
  } else {
    w <- "B"
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, w)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Changing Factor")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(w))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    simulations <- as.numeric(readline(writeLines(c("The parameter S (numerical type) is equivalent to the number of simulations."," ",
                                                    "   Value by default: 1E3", "   Possible values: all natural numbers",
                                                    "   Max. recommended to avoid an excess of computational time: < 1E6", " ",
                                                    "What is its value?", " "))))
    S <- round(simulations)
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, S)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Number of simulations")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(S))
  } else {
    S <- 1E3
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, S)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Number of simulations")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(S))
  }
  
  cat("\014") # Clean console
  
  if (mod_param == "y") {
    plotheatmap <- as.logical(readline(
      writeLines(c("The parameter plot (character type) determine whether or not to show the heatmap of the factor distribution.",
                   " ", "   Value by default: F", "   Possible values: F or T", " ", "What is its value?", " "))))
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, plotheatmap)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Factors heatmap")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(plotheatmap))
  } else {
    plotheatmap <- F
    Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, plotheatmap)
    Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Factors heatmap")
    Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(plotheatmap))
  }
  
  cat("\014") # Clean console
  
  # if (mod_param == "y") {
  #   k <- as.numeric(readline(writeLines(c("Chunk number (numeric type) determine whether to compute tstats or ado.",
  #                                         " ", "   Value by default: 0", "   Possible values: 0 (ado) or not 0 (tstats)",
  #                                         " ", "What is its value?", " "))))
  #   Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, k)
  #   Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "k")
  #   Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(k))
  # } else {
  #   k <- 0
  #   Parameter_mod_vect_vals <- append(Parameter_mod_vect_vals, k)
  #   Parameter_mod_vect_names <- append(Parameter_mod_vect_names, "Chunk number")
  #   Parameter_mod_vect_class <- append(Parameter_mod_vect_class, class(k))
  # }
  # 
  # cat("\014") # Clean console
  
  
  Parameter_mod_DF = data.frame(Parameter_mod_vect_vals, Parameter_mod_vect_class)
  row.names(Parameter_mod_DF) = Parameter_mod_vect_names
  colnames(Parameter_mod_DF) = c("Parameter value", "Parameter type")
  
  D = "unif-0-1"
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(D))
  l = 1000
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(l))
  s = 0.1
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(s))
  p = 1
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(p))
  p_dist = "norm"
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(p_dist))
  H = 1
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(H))
  d = 0
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(d))
  t = "none"
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(t))
  # o = "outputfiletemp"
  # Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(o))
  # k = 1
  # Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(k))
  x = 10
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(x))
  f = getwd()
  Parameter_no_mod_vect_class <- append(Parameter_no_mod_vect_class, class(f))

  
  Parameter_no_mod_vect_names <- c("Multivariate non-normal distribution definition", "Lambda parameter of the Poisson distribution",
                                   "Stdev for the 'simplex' generator model", "Location of the 'simplex' generator model",
                                   "Distribution for 'simplex' generator model: norm, gamma or beta", "Heteroskedasticity degree (level 1)",
                                   "Delta (H1 generation parameter)", "Data transformation: sqrt, log or none", # "Output file name",
                                   "Number of cores", "Path to helper functions")
  Parameter_no_mod_vect_vals <- c(D, l, s, p, p_dist, H, d, t, x, f)

  Parameter_no_mod_DF = data.frame(Parameter_no_mod_vect_vals, Parameter_no_mod_vect_class)
  row.names(Parameter_no_mod_DF) = Parameter_no_mod_vect_names
  colnames(Parameter_no_mod_DF) = c("Parameter value", "Parameter type")

  
  cat("\014") # Clean console
  
  # Se asignan los parámetros modificados o no al entorno global:
  Parameter_mod_List <- list("q" = q, "a" = a, "b" = b, "n" = n, "u" = u,
                             "m" = m, "c" = c, "v" = v, "w" = w, "S" = S,
                             "plotheatmap" = plotheatmap,
                             "Parameter_DF" = Parameter_mod_DF,
                             "Parameter_mod_vect_vals" = Parameter_mod_vect_vals,
                             "Parameter_mod_vect_names" = Parameter_mod_vect_names,
                             "Parameter_mod_vect_class" = Parameter_mod_vect_class)
  list2env(Parameter_mod_List, .GlobalEnv)
  
  Parameter_no_mod_List <- list("D" = D, "l" = l, "s" = s, "p" = p, "p_dist" = p_dist,
                                "H" = H, "d" = d, "t" = t, "x" = x, "f" = f)
  list2env(Parameter_no_mod_List, .GlobalEnv)

  message1 = "The following parameters have been modified and assigned to the global environment: "
  writeLines(paste(message1, "\n"))
  print(Parameter_mod_DF)
  writeLines("\n")
  message2 = "While the values of these other parameters are not modifiable: "
  writeLines(paste(message2, "\n"))
  print(Parameter_no_mod_DF)
  writeLines("\n")
  
  
  # Se cargan los paquetes, instalándolos si es necesario:

  #    require() => a warning if a package is not installed and then continue to execute the code
  #    library() => an error and stop the execution of the code

  if (!require("devtools")) install.packages("devtools")
  require(devtools)
  library(devtools)
  
  if (!require("versions")) install.packages("versions")
  require(versions)
  library(devtools)
  
  if (!require("manta")) devtools::install_github("dgarrimar/manta")
  require(manta)
  library(devtools)
  
  package_list <- c("CompQuadForm", "car", "Hmisc", "MASS", "copula", "vegan",
                    "randomcoloR", "plyr", "svMisc", "reshape2", "hms", "progress",
                    "ggplot2", "arm", "pheatmap")
  # Hmisc => Necesario para la función label
  # randomcoloR => La función randomcoloR() genera un color aleatóriamente
  # plyr => Necesario para la función mapvalues()
  # svMisc => Permite usar progress() para indicar visualmente el progreso de la simulación
  # reshape2 => La función melt() permite convertir un objeto "dist" en "DF"
  # hms => La función as_hms() permite convertir "difftime" en DD:HH:MM:SS
  # progress => Para mostrar una barra de progreso de la simulación
  
  for (i in package_list){
    if(! i %in% installed.packages()){
      install.packages(i)
    }
    require(i, character.only = TRUE)
    library(i, character.only = TRUE)
  }

  
  # Se indica el archivo de soporte que incluye las funciones necesarias.
  source(sprintf("%s/fx.R", f))
  
  
  # Equivalencias de parámetros.
  modelSim <- m
  fx <- f
  stdev <- s
  pdist <- p_dist
  delta <- d
  Lambda <- l
  Cor <- c
  Var <- v
  # chunk <- k
  hk <- H
  transf <- t
  cores <- x
  # output <- o
  DistDef <- dd <- D 
  position <- loc <- p
  
  # Se asignan los parámetros generados en esta sección:
  Parameter_gen1_List <- list("modelSim" = m, "fx" = f, "stdev" = s, "pdist" = p_dist, "delta" = d,
                              "Lambda" = l, "Cor" = c, "Var" = v, # "chunk" = k,
                              "hk" = H, "transf" = t, "cores" = x, # "output" = o,
                              "DistDef" = D, "dd" = D, "position" = p, "loc" = p)
  list2env(Parameter_gen1_List, .GlobalEnv)
  
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  #   Niveles de los factores y su interacción  #
  # # # # # # # # # # # # # # # # # # # # # # # #
  
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
  
  # Se asignan los parámetros generados en esta sección:
  Parameter_gen2_List <- list("ch" = ch, "A" = A, "B" = B)
  list2env(Parameter_gen2_List, .GlobalEnv)
  
  
  # Para: m = modelSim = "simplex" 
  # Se carga la tabla correspondiente de "stdev":
  if (modelSim == "simplex") {
    
    tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", fx, p_dist), h = T)
    colnames(tbl) <- c("Q", "L", "S")
    tbl_name <- "simplex_table_Q_L_S"
    write.csv(tbl, file = paste(tbl_name, ".csv", sep = ""), row.names=FALSE)
    assign(tbl, tbl, envir = .GlobalEnv)
    
    if (! q %in% unique(tbl$Q)) {
      stop(sprintf("stdev not precomputed for q = %s", q))
    } 
    
    # Se busca la "stdev" que corresponde a los valores [q, loc] seleccionados.
    stdev <- subset(tbl, Q == q & L == loc)$S
    
    # Se crea y guarda un DF con los valores [Q, L, S] que cumplen la condición,
    # asignándole un nombre que refleja los parámetros seleccionados:
    stdevDF <- data.frame(subset(tbl, Q == q & L == loc)$Q,
                          subset(tbl, Q == q & L == loc)$L,
                          subset(tbl, Q == q & L == loc)$S)
    colnames(stdevDF) <- c("q (Q)", "loc (L)", "stdev (S)")
    stdev_name = paste("stdev_", modelSim, sep = "")
    assign(stdev_name, stdevDF[[3]], envir = .GlobalEnv)
    
    # Se asignan al entorno los parámetros generados en esta sección:
    Parameter_gen3_List <- list("tbl" = tbl)
    list2env(Parameter_gen3_List, .GlobalEnv)
    
  }
  
  
  # Se crea el directorio de la simulación vigente:
  work_dir <- getwd()
  results_dir <- "Resultados"
  sim_dir_name = paste("[model = ", m, ", a = ",  a, ", b = ",  b, ", n = ",  n,
                        ", u = ",  u, ", w = ",  w, ", S = ",  S, ", q (Q) = ", q,
                        ", loc (L) = ",  loc, "]", sep = "")
  sim_path = file.path(work_dir, results_dir, sim_dir_name)
  dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
  assign("sim_path", sim_path, envir = .GlobalEnv)
  
  
  # Guardamos la información necesaria hasta este punto:
  if (modelSim == "simplex") {
    
    write.csv(stdevDF, file = file.path(sim_path, paste(stdev_name, ".csv", sep = "")), row.names=FALSE)
    
  }

  pdf(file = file.path(sim_path, "FactorLevelsHeatmap.pdf"), width = 8, height = 8)
  pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
  dev.off()

  
  # # # # # # # # # # # # # # # # # # # # # #
  #   Simulación según el modelo aplicado   #
  # # # # # # # # # # # # # # # # # # # # # #
  
  t_start_simulation <- Sys.time()
  writeLines(paste("\n", "Data simulation is running at", format(t_start_simulation), "\n"))
  
  pb <- txtProgressBar(min = 0, max = S, style = 3, width = 50, char = "=") 
  
  tstats <- c()
  set.seed(1325)
  
    
    for (i in 1:S) {
      
      # Para un modelo Simplex:
      if (modelSim == "simplex") {
        
        if (i == 1) { # Sanity check
          tol <- 0.01
          check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, "norm")
          wm <- which.max(check$exp)
          if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
            stop("Deviation from expected centroid greater than tolerance.")
          }
        }
        
        Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, pdist)
        sd <- mean(apply(Y, 2, sd))
        
        assign("Y", Y, envir = .GlobalEnv)
        assign("sd", sd, envir = .GlobalEnv)
        
        # Transformaciones:
        Y_sqrt <- sqrt(Y)
        Y_log <- log(Y + 1)
        
      } else if (modelSim == "mvnorm") {
        
        Y <- Sim.mvnorm(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor) # Cambiar 0 <-> núm grande (para transformaciones)
        
        assign("Y", Y, envir = .GlobalEnv)
        
      } else if (modelSim == "multinom") {
        
        # N <- rpois(1, Lambda)
        Y <- Sim.multinom(ch, q, n, Lambda, delta, loc)
        
        assign("Y", Y, envir = .GlobalEnv)
        
      } else if (modelSim == "copula") {
        
        Y <- Sim.copula(ch, q, n, mu = rep(0, q), delta, hk, Var, Cor, dd)
        
        assign("Y", Y, envir = .GlobalEnv)
        
      } else {
        stop(sprintf("Unknown option: modelSim = '%s'.", modelSim))
      }
      
      # # Transformaciones sobre las variables simuladas: 
      # if (transf == "sqrt") {
      #   Y <- sqrt(Y)
      # } else if (transf == "log") {
      #   Y <- log(Y+1)
      # } else if (transf == "none") {
      #   
      # } else {
      #   stop(sprintf("Unknown option: transf = '%s'.", transf))
      # }
      
      # # lm and residuals
      # Y <- scale(Y, center = T, scale = F)
      # fit <- lm(Y ~ A + B + A:B)
      # R <- fit$residuals
      # 
      # # Sums of squares
      # UU <- Anova(fit, type = "II") # SS type II 
      # SS <- lapply(UU$SSP, function(x){sum(diag(x))})
      # SSe <- sum(diag(UU$SSPE))
      # 
      # # Df
      # df.e <- fit$df.residual       # df.e <- (n-1) - sum(Df)
      # Df <- table(fit$assign)[-1]
      # names(Df) <- attributes(fit$terms)$term.labels
      # 
      # # Statistic 
      # f <- unlist(lapply(SS, function(x){x/SSe})) 
      # # We divide by SSe for numerical stability when running CompQuadForm::davies
      # Fs <- f/Df*df.e 
      # 
      # # Asymptotic test statistic (McArtor) 
      # G <- crossprod(Y)
      # eG <- eigen(G/n, only.values = T)$values # n or df.e ?
      # 
      # # Asymptotic test statistic (Our proposal)
      # e <- eigen(cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
      
      # # Se indica visualmente el progreso de la simulación:
      # progress(i, progress.bar = TRUE)
      
      # tstats_func(Y)
      # tstats <- rbind(tstats, c(tstats_func(Y)$Fs, unlist(tstats_func(Y)$SS),
      #                           tstats_func(Y)$SSe, tstats_func(Y)$e, tstats_func(Y)$eG))
      
      manta(Y ~ A + B + A:B)
      manova(Y ~ ., data = data.frame(covariates, snp = X[, snp]))
      summary(manova(Y ~ A + B + A:B))
      
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
  
    # Tiempo de ejecución de la simulación
    t_simulation <- difftime(Sys.time(), t_start_simulation, units="auto")
    assign("t_simulation", t_simulation, envir = .GlobalEnv)
    writeLines(paste("\n", "  It was performed in ", format(t_simulation), sep =""))
    
    
    colnames(tstats) <- c("Fs (A)", "Fs (B)", "Fs (A:B)",
                          "SS (A)", "SS (B)", "SS (A:B)", "SSe",
                          "e (A)", "e (B)", "e (A:B)",
                          "eG (A)", "eG (B)", "eG (A:B)")
    write.csv(tstats, file = file.path(sim_path, paste("tstats.csv", sep = "")), row.names = F, quote = F)

    
    if (transf == "sqrt") {
      Y_transf <- sqrt(Y)
      assign("Y_transf", Y_transf, envir = .GlobalEnv)
    } else if (transf == "log") {
      Y_transf <- log(Y+1)
      assign("Y_transf", Y_transf, envir = .GlobalEnv)
    } else if (transf == "none") {
      Y_transf <- Y
      assign("Y_transf", Y_transf, envir = .GlobalEnv)
    } else {
      stop(sprintf("Unknown data transformation option: transf = '%s'.", transf))
    }
    
    Y_transf_DF <- data.frame("A" = Y_transf[,1], "B" = Y_transf[,2], "AB" = Y_transf[,3])
    assign("Y_transf_DF", Y_transf_DF, envir = .GlobalEnv)
       
    Y_transf_scaled <- scale(Y_transf, center = T, scale = F)
    assign("Y_transf_scaled", Y_transf_scaled, envir = .GlobalEnv)
    
    d_Y_transf_scaled <- as.dist(interDist(Y_transf_scaled))
    assign("d_Y_transf_scaled", d_Y_transf_scaled, envir = .GlobalEnv)
    
    d_Y_transf_scaled_DF <- melt(as.matrix(d_Y_transf_scaled), varnames = c("row", "col"), value.name = "distance")
    assign("d_Y_transf_scaled_DF", d_Y_transf_scaled_DF, envir = .GlobalEnv)
    
    # Gráfico de la matriz de distancias:
    png(file = file.path(sim_path, "DistanceMatrix.png"), width = 600, height = 600)
    # pdf(file = file.path(sim_path, "DistanceMatrix.pdf"), width = 8, height = 8)
    
    ggplot(d_Y_transf_scaled_DF, aes(row, col)) +
      geom_tile(aes(fill = distance), colour = "white") +
      scale_fill_gradient(low = "white", high = randomColor()) +
      geom_text(aes(label = round(distance, 1)), size = 1)
    
    dev.off()
    
    
    # # # # # #
    #  MANTA  #
    # # # # # #
    
    ## MANTA with Y not scaled
    t_start_MANTA_not_scaled <- Sys.time()
    
    MANTA_t_end_prev <- t_start_MANTA_not_scaled + system.time({manta(Y_transf ~ A + B + A:B, data = Y_transf_DF)})[3]
    writeLines(paste("\n", "\n", "\n", "MANTA is running at", format(t_start_MANTA_not_scaled),
                     "\n", "\n", "It will finish approximately at", MANTA_t_end_prev))
    
    MANTA_results_not_scaled <- manta(Y_transf ~ A + B + A:B, data = Y_transf_DF)
    
    # Tiempo de ejecución del cálculo mediante MANTA scaled
    t_MANTA_not_scaled <- difftime(Sys.time(), t_start_MANTA_not_scaled, units="auto")
    assign("t_MANTA_not_scaled", t_MANTA_not_scaled, envir = .GlobalEnv)  
    writeLines(paste(" It was performed in ", format(t_MANTA_not_scaled), "\n", "\n", " Resulting in:", sep =""))
    print(MANTA_results_not_scaled)
    
    capture.output(MANTA_results_not_scaled, file = file.path(sim_path, paste("MANTA_results_not_scaled.txt", sep = "")))
    capture.output(MANTA_results_not_scaled, file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    assign("MANTA_results_not_scaled", MANTA_results_not_scaled, envir = .GlobalEnv)
    capture.output(writeLines(paste("\n", "\n", strrep("-", 88), "\n", "\n")),
                   file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)

    
    # # # # #
    #  ADO  #
    # # # # #
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #   WARNING: Adonis uses Type I SS, while we use Type II SS!  #
    #                                                             #
    #            => Only valid to study A:B under H0              #
    #            => Always with Y scaled                          #
    #                                                             #
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    ## Adonis permutations (raw data + type I SS)
    t_start_ado_raw_type_I_SS <- Sys.time()
    
    ado_raw_type_I_SS_t_end_prev <- t_start_ado_raw_type_I_SS + 
      system.time({adonis2(d_Y_transf_scaled ~ A + B + A:B, permutations = S, parallel = cores)})[3]
    writeLines(paste("\n", "\n", "\n", "ADO (raw data + type I SS) is running at", format(t_start_ado_raw_type_I_SS),
                     "\n", "\n", "It will finish approximately at", ado_raw_type_I_SS_t_end_prev))

    ado_raw_type_I_SS <- adonis2(d_Y_transf_scaled ~ A + B + A:B, permutations = S, parallel = cores)
    
    # Tiempo de ejecución del cálculo mediante Adonis permutations (raw data + type I SS)
    t_ado_raw_type_I_SS <- difftime(Sys.time(), t_start_ado_raw_type_I_SS, units="auto")
    assign("t_ado_raw_type_I_SS", t_ado_raw_type_I_SS, envir = .GlobalEnv)
    writeLines(paste(" It was performed in ", format(t_ado_raw_type_I_SS), "\n", "\n", " Resulting in:", "\n", sep =""))
    print(ado_raw_type_I_SS)
    
    capture.output(ado_raw_type_I_SS, file = file.path(sim_path, paste("ado_raw_type_I_SS.txt", sep = "")))
    capture.output(ado_raw_type_I_SS, file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    assign("ado_raw_type_I_SS", ado_raw_type_I_SS, envir = .GlobalEnv)
    capture.output(writeLines(paste("\n", "\n", strrep("-", 88), "\n", "\n")),
                   file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)

    
    ## Adonis permutations (strata + type I SS)
    perm <- how(nperm = S)
    setBlocks(perm) <- ch
    
    t_start_ado_strata_type_I_SS <- Sys.time()
    
    ado_strata_type_I_SS_t_end_prev <- t_start_ado_strata_type_I_SS +
      system.time({adonis2(d_Y_transf_scaled ~ A + B + A:B, permutations = perm, parallel = cores)})[3]
    writeLines(paste("\n", "\n", "\n", "ADO (strata + type I SS) is running at", format(t_start_ado_strata_type_I_SS),
                     "\n", "\n", "It will finish approximately at", ado_strata_type_I_SS_t_end_prev))
    
    ado_strata_type_I_SS_exec_prev <- as.numeric()
    ado_strata_type_I_SS <- adonis2(d_Y_transf_scaled ~ A + B + A:B, permutations = perm, parallel = cores)
    
    # Tiempo de ejecución del cálculo mediante Adonis permutations (strata + type I SS)
    t_ado_strata_type_I_SS <- difftime(Sys.time(), t_start_ado_strata_type_I_SS, units="auto")
    assign("t_ado_strata_type_I_SS", t_ado_strata_type_I_SS, envir = .GlobalEnv)
    writeLines(paste(" It was performed in ", format(t_ado_strata_type_I_SS), "\n", "\n", " Resulting in:", "\n", sep =""))
    print(ado_strata_type_I_SS)
    
    capture.output(ado_strata_type_I_SS, file = file.path(sim_path, paste("ado_strata_type_I_SS.txt", sep = "")))
    capture.output(ado_strata_type_I_SS, file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    assign("ado_strata_type_I_SS", ado_strata_type_I_SS, envir = .GlobalEnv)
    capture.output(writeLines(paste("\n", "\n", strrep("-", 88), "\n", "\n")),
                   file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    
    # ado <- list()
    # ado$f.perms <- matrix(NA, nrow = 1, ncol = 7)
    
    # Yp <- Y
    # for (p in 1:S) {
    #     print(p)
    #     for (i in 1:b) {Yp[B == i, ] <- Y[sample(which(B == i)),]}
    #     fit <- lm(Yp ~ A + B + A:B)
    #     R <- fit$residuals
    #     UU <- Anova(fit, type = "II") # SS type II
    #     SS <- lapply(UU$SSP, function(x){sum(diag(x))})
    #     SSe <- sum(diag(UU$SSPE))
    #     df.e <- fit$df.residual
    #     Df <- table(fit$assign)[-1]
    #     f <- unlist(lapply(SS, function(x){x/SSe}))
    #     ado$f.perms <- rbind(ado$f.perms, cbind(t(f/Df*df.e), t(as.numeric(SS)), t(SSe)))
    # }
    # ado$f.perms <- ado$f.perms[-1,]
    
    
    # # # # # # # # # # # # # #
    #  Tiempos de ejecución   #
    # # # # # # # # # # # # # #
    
    ExecutionTime_DF <- data.frame(as_hms(c(t_simulation, t_MANTA_not_scaled,
                                            t_ado_raw_type_I_SS, t_ado_strata_type_I_SS)))
    row.names(ExecutionTime_DF) = c("Simulation", "MANTA with Y not scaled",
                                    "Adonis with permutations and Y scaled (raw data + type I SS)",
                                    "Adonis with permutations and Y scaled (strata + type I SS)")
    colnames(ExecutionTime_DF) = c("Execution time (DD:HH:MM:SS)")
    capture.output(ExecutionTime_DF, file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    assign("ExecutionTime_DF", ExecutionTime_DF, envir = .GlobalEnv)
    capture.output(writeLines(paste("\n", "\n", strrep("-", 88), "\n", "\n")),
                   file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    
    
    # # # # # # # # # # # # # #
    #  Resumen de resultados  #
    # # # # # # # # # # # # # #

    Pr_gt_F_DF <- data.frame(c(MANTA_results_not_scaled$aov.tab[1,6],
                               ado_raw_type_I_SS$`Pr(>F)`[1],
                               ado_strata_type_I_SS$`Pr(>F)`[1]),
                             c(MANTA_results_not_scaled$aov.tab[2,6],
                               ado_raw_type_I_SS$`Pr(>F)`[2],
                               ado_strata_type_I_SS$`Pr(>F)`[2]),
                             c(MANTA_results_not_scaled$aov.tab[3,6],
                               ado_raw_type_I_SS$`Pr(>F)`[3],
                               ado_strata_type_I_SS$`Pr(>F)`[3]),
                             c(as.character(as_hms(t_MANTA_not_scaled)),
                               as.character(as_hms(t_ado_raw_type_I_SS)),
                               as.character(as_hms(t_ado_strata_type_I_SS))))
    row.names(Pr_gt_F_DF) = c("MANTA with Y not scaled",
                              "ADO + Permut (Y scaled, raw data, type I, SS)",
                              "ADO + Permut (Y scaled, strata, type I, SS)")
    colnames(Pr_gt_F_DF) <- c("Pr(>F) A", "Pr(>F) B", "Pr(>F) A:B", "Exec. Time (h:m:s)")
    capture.output(Pr_gt_F_DF, file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    assign("Pr_gt_F_DF", Pr_gt_F_DF, envir = .GlobalEnv)
    capture.output(writeLines(paste("\n", "\n", strrep("-", 88), "\n", "\n")),
                   file = file.path(sim_path, paste("Results.txt", sep = "")), append = TRUE)
    
    writeLines(paste("\n", "\n", "\n", "Results summary:", "\n"))
    print(Pr_gt_F_DF)
}

cat("\014") # Clean console


### Aplicación de la función:

# "y" => Se deberá introducir el valor deseado para todos los parámetros.
# "n" => Se considerarán todos los parámetros con su valor por defecto.

Parameter_mod_DF()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# #                                                                                   # # #
# #   Loop de simulación y aplicación de los diferentes métodos para diversos casos   # # #
# #                                                                                   # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Se eliminan todos los objetos y variables del entorno.
rm(list = ls(all.names = TRUE))


# # # # # # # # # # # # # # # # # # # # # # # # # #
#  Instalación y carga de los paquetes necesarios #
# # # # # # # # # # # # # # # # # # # # # # # # # #

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
                  "ggplot2", "arm", "pheatmap")

# Hmisc => Necesario para la función label
# randomcoloR => La función randomcoloR() genera un color aleatóriamente
# plyr => Necesario para la función mapvalues()
# svMisc => Permite usar progress() para indicar visualmente el progreso de la simulación
# reshape2 => La función melt() permite convertir un objeto "dist" en "DF"
# hms => La función as_hms() permite convertir "difftime" en DD:HH:MM:SS
# progress => Para mostrar una barra de progreso de la simulación

for (i in package_list){
  if(! i %in% installed.packages()){
    install.packages(i)
  }
  require(i, character.only = TRUE)
  library(i, character.only = TRUE)
}

cat("\014") # Clean console


# # # # # # # # # # # # # # # # # # # #
#  Asignación y gestión de variables  #
# # # # # # # # # # # # # # # # # # # #

# Asignación y gestión de las diferentes variables (modificables o no).

# Dependent variables: the parameter q (numerical type) is equivalent to the number of dependent variables.
# (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 4)
q_values <- c(1:4)

# Factor A Levels: the parameter a (numerical type) is equivalent to the number of levels of factor A.
# (Value by default: 2 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
a_values <- c(1:6)

# Factor B Levels: the parameter b (numerical type) is equivalent to the number of levels of factor B."
# (Value by default: 3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 6)
b_values <- c(1:6)

# Number of samples: the parameter n (numerical type) is equivalent to the number of samples.
# (Value by default: 300 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 400)
n_values <- c(seq(200, 400, 100))

# Unbalance 1st B Level: the parameter u (numerical type) is equivalent to the unbalance degree of the first level regarding the rest of the levels of the factor B.
# (Value by default: 1 / Possible values: all positive rational numbers / Max. recommended to avoid an excess of interactions: < 10)
u_values <- c(seq(0.1, 10, 0.1))

# Model H0/H1: the parameter m (character type) is equivalent to the model that generates H0/H1.
# (Value by default: mvnorm / Possible values: mvnorm, simplex or multinom)
m_values <- c("mvnorm", "simplex", "multinom")

# Y correlation: the parameter c (numerical type) is equivalent to the correlation of the Y variables.
# (Value by default: 0 / Possible values: all rational numbers between -1.0 and 1.0)
c_values <- c(seq(0, .8, 0.2)) # Positivo y < 1 siempre.

# Y variance: the parameter v (character type) is equivalent to the variance of the Y variables.
# (Value by default: equal / Possible values: equal or unequal)
v_values <- c("equal", "unequal")

# Changing Factor: the parameter w (character type) determines which factor changes.
# (Value by default: B / Possible values: A, B or AB)
w_values <- c("A", "B", "AB")

# Number of simulations: the parameter S (numerical type) is equivalent to the number of simulations.
# (Value by default: 1E3 / Possible values: all natural numbers / Max. recommended to avoid an excess of interactions: < 1E6)
S_values <- c(seq(1E2, 1E4, 1E2))

# Factors heatmap: the parameter plot (character type) determine whether or not to show the heatmap of the factor distribution.
# (Value by default: F / Possible values: F or T)
plotheatmap_values <- c("F", "T")

# Chunk number: the parameter k or Chunk number (numeric type) determine whether to compute tstats or ado.
# (Value by default: 0 / Possible values: 0 (ado) or not 0 (tstats))
k_values <- c(0, 1)

# Se asignan las variables generadas a los diferentes parámetros de simulación:
# q, a, b, n, u, m, c, v, w, S, plotheatmap, k, D, l, s, p, p_dist, H, d, t, x, f.
q = as.numeric(q_values[3])
a = as.numeric(a_values[2])
b = as.numeric(b_values[3])
n = as.numeric(n_values[1])
u = as.numeric(u_values[10])
m = as.character(m_values[1])
c = as.numeric(c_values[1])
v = as.character(v_values[1])
w = as.character(w_values[2])
S = as.numeric(S_values[1])
plotheatmap = as.logical(plotheatmap_values[2])
k = as.numeric(k_values[1])
D = as.character("unif-0-1")
l = as.numeric(1000)
s = as.numeric(0.1)
p = as.numeric(1)
p_dist = as.character("norm")
H = as.numeric(1)
d = as.numeric(0)
t = as.character("none")
# o = "outputfiletemp"
# k = 1
x = as.numeric(10)
f = as.character(getwd())

# Se indica el archivo de soporte que incluye las funciones.
# Funciones que se usarán: dM, dH, dE, interDist, geodesic, step2distance,
# label, step2h1, sim.simplex, Sim.simplex, sim.mvnorm, Sim.mvnorm.
source(sprintf("%s/fx.R", f))

# Se crea una lista con los parámetros modificados o no, y se asignan al entorno global.
Parameter_List <- list("q" = q, "a" = a, "b" = b, "n" = n, "u" = u, "m" = m, "c" = c,
                       "v" = v, "w" = w, "S" = S, "plotheatmap" = plotheatmap, "k" = k,
                       "D" = D, "l" = l, "s" = s, "p" = p, "p_dist" = p_dist, "H" = H,
                       "d" = d, "t" = t, "x" = x, "f" = f)
list2env(Parameter_List, .GlobalEnv)

# Equivalencias de parámetros para las diferentes funciones de "fx.R".
modelSim <- m
fx <- f
stdev <- s
pdist <- p_dist
delta <- d
Lambda <- l
Cor <- c
Var <- v
# chunk <- k
hk <- H
transf <- t
cores <- x
# output <- o
DistDef <- dd <- D 
position <- loc <- p

# Se asignan al entorno los nuevos parámetros generados en esta sección:
Parameter_gen1_List <- list("modelSim" = m, "fx" = f, "stdev" = s, "pdist" = p_dist, "delta" = d,
                            "Lambda" = l, "Cor" = c, "Var" = v, "chunk" = k,
                            "hk" = H, "transf" = t, "cores" = x, # "output" = o,
                            "DistDef" = D, "dd" = D, "position" = p, "loc" = p)
list2env(Parameter_gen1_List, .GlobalEnv)


# # # # # # # # # # # # # # # # # #
# DF de los parámetros del modelo #
# # # # # # # # # # # # # # # # # #

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
                     "Delta (H1 generation parameter)",
                     "Data transformation: sqrt, log or none",
                     # "Output file name",
                     "Number of cores",
                     "Path to helper functions")

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
                     "d / delta",
                     "t / transf",
                     # "Output file name",
                     "x / cores",
                     "f / fx")

Parameter_values <- c()
for (i in 1:length((Parameter_List))){
  Parameter_values <- append(Parameter_values, Parameter_List[[i]])}

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

cat("\014") # Clean console


# # # # # # # # # # # # # # # # # # # # # # # #
#   Niveles de los factores y su interacción  #
# # # # # # # # # # # # # # # # # # # # # # # #

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

# Función de simulación:
CompMantaManovafunc <- function() {
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #  Asignación de la "stdev" para el modelo "simplex"  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Para: m = modelSim = "simplex" 
  # Se carga la tabla correspondiente de "stdev":
  if (modelSim == "simplex") {
    
    tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", fx, p_dist), h = T)
    colnames(tbl) <- c("Q", "L", "S")
    tbl_name <- "simplex_table_Q_L_S"
    write.csv(tbl, file = paste(tbl_name, ".csv", sep = ""), row.names=FALSE)
    assign(tbl, tbl, envir = .GlobalEnv)
    
    if (! q %in% unique(tbl$Q)) {
      stop(sprintf("stdev not precomputed for q = %s", q))} 
    
    # Se busca la "stdev" que corresponde a los valores [q, loc] seleccionados.
    stdev <- subset(tbl, Q == q & L == loc)$S
    
    # Se crea y guarda un DF con los valores [Q, L, S] que cumplen la condición,
    # asignándole un nombre que refleja los parámetros seleccionados:
    stdevDF <- data.frame(subset(tbl, Q == q & L == loc)$Q,
                          subset(tbl, Q == q & L == loc)$L,
                          subset(tbl, Q == q & L == loc)$S)
    colnames(stdevDF) <- c("q (Q)", "loc (L)", "stdev (S)")
    stdev_name = paste("stdev_", modelSim, sep = "")
    assign(stdev_name, stdevDF[[3]], envir = .GlobalEnv)
    
    # Se asignan los parámetros generados en esta sección:
    Parameter_gen3_List <- list("tbl" = tbl)
    list2env(Parameter_gen3_List, .GlobalEnv)}
  
  
  # # # # # # # # # # # # # # # # # # # # # #
  #  Almacenaje de los datos de simulación  #
  # # # # # # # # # # # # # # # # # # # # # #
  
  # Se crea el directorio de la simulación vigente:
  work_dir <- getwd()
  results_dir <- "Resultados"
  sim_dir_name = paste("[model = ", m, ", a = ",  a, ", b = ",  b, ", n = ",  n,
                       ", u = ",  u, ", w = ",  w, ", S = ",  S, ", q (Q) = ", q,
                       ", loc (L) = ",  loc, "]", sep = "")
  sim_path = file.path(work_dir, results_dir, sim_dir_name)
  dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
  assign("sim_path", sim_path, envir = .GlobalEnv)
  
  
  # Guardamos la información necesaria hasta este punto:
  if (modelSim == "simplex") {
    
    write.csv(stdevDF, file = file.path(sim_path, paste(stdev_name, ".csv", sep = "")), row.names=FALSE)
    
  }
  
  pdf(file = file.path(sim_path, "FactorLevelsHeatmap.pdf"), width = 8, height = 8)
  pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
  dev.off()
  
  
  # # # # # # # # # # # # # # # # # # # #
  #  Simulación del conjunto de datos   #
  #  y aplicación de MANTA y MANOVA     #
  # # # # # # # # # # # # # # # # # # # #
  
  # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
  # las p del factor B para cada método y conjunto de datos considerado:
  Manta_pB = rep(NA, S)
  Manova_pB = rep(NA, S)
  Manta_Ylog_pB = rep(NA, S)
  Manova_Ylog_pB = rep(NA, S)
  Manta_Ysqrt_pB = rep(NA, S)
  Manova_Ysqrt_pB = rep(NA, S)
  Manta_Yscaled_pB = rep(NA, S)
  Manova_Yscaled_pB = rep(NA, S)
  Manta_YnormMinMax_pB = rep(NA, S)
  Manova_YnormMinMax_pB = rep(NA, S)
  
  for (i in 1:S) {
    if (modelSim == "simplex") {
      if (i == 1) { # Sanity check
        tol <- 0.01
        check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, pdist)
        wm <- which.max(check$exp)
        if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
          stop("Deviation from expected centroid greater than tolerance.")
        }
      }
      Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = F, pdist)
    } else if (modelSim == "mvnorm") {
      Y <- Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor) # Cambiar 0 <-> núm grande (para transformaciones)
    } else if (modelSim == "multinom") {
      Y <- Sim.multinom(ch, q, n, Lambda, delta, loc)
    }
    
    # Transformaciones posibles de los datos simulados (Y):
    Ylog <- log(Y) # Transformación logarítmica.
    Ysqrt <- sqrt(Y) # Transformación sqrt.
    Yscaled <- cbind(Y[,1]/sd(Y[,1]), Y[,2]/sd(Y[,2]), Y[,3]/sd(Y[,3])) # Escalado de datos: se divide cada valor por la desviación típica del conjunto de datos.
    YnormMinMax <-  cbind((Y[,1]-min(Y[,1]))/(max(Y[,1])-min(Y[,1])), # Normalización Min- Max: si queremos que cada variable contribuya por igual al análisis.
                          (Y[,2]-min(Y[,2]))/(max(Y[,2])-min(Y[,2])),
                          (Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3])))
    
    # Aplicación de Manta y Manova para cada conjunto de datos (Y, Ylog, Ysqrt):
    Manta_pB[i] <- manta(Y ~ A + B + A:B)$aov.tab[2,6]
    Manova_pB[i] <- summary(manova(Y ~ A + B + A:B))$stats[2,6]
    
    Manta_Ylog_pB[i] <- manta(Ylog ~ A + B + A:B)$aov.tab[2,6]
    Manova_Ylog_pB[i] <- summary(manova(Ylog ~ A + B + A:B))$stats[2,6]
    
    Manta_Ysqrt_pB[i] <- manta(Ysqrt ~ A + B + A:B)$aov.tab[2,6]
    Manova_Ysqrt_pB[i] <- summary(manova(Ysqrt ~ A + B + A:B))$stats[2,6]
    
    Manta_Yscaled_pB[i] <- manta(Yscaled ~ A + B + A:B)$aov.tab[2,6]
    Manova_Yscaled_pB[i] <- summary(manova(Yscaled ~ A + B + A:B))$stats[2,6]
    
    Manta_YnormMinMax_pB[i] <- manta(YnormMinMax ~ A + B + A:B)$aov.tab[2,6]
    Manova_YnormMinMax_pB[i] <- summary(manova(YnormMinMax ~ A + B + A:B))$stats[2,6]
  }
  
  # Establecemos el nivel de significación α:
  alpha = 0.05
  
  # Se calcula la potencia de cada método para el α dado:
  Manta_pot <- mean(Manta_pB < alpha)
  Manova_pot <- mean(Manova_pB < alpha)
  Y_Manta_Manova_pot <- cbind(Manta_pot, Manova_pot)
  
  Manta_Ylog_pot <- mean(Manta_Ylog_pB < alpha)
  Manova_Ylog_pot <- mean(Manova_Ylog_pB < alpha)
  Ylog_Manta_Manova_pot <- cbind(Manta_Ylog_pot, Manova_Ylog_pot)
  
  Manta_Ysqrt_pot <- mean(Manta_Ysqrt_pB < alpha)
  Manova_Ysqrt_pot <- mean(Manova_Ysqrt_pB < alpha)
  Ysqrt_Manta_Manova_pot <- cbind(Manta_Ysqrt_pot, Manova_Ysqrt_pot)
  
  Manta_Yscaled_pot <- mean(Manta_Yscaled_pB < alpha)
  Manova_Yscaled_pot <- mean(Manova_Yscaled_pB < alpha)
  Yscaled_Manta_Manova_pot <- cbind(Manta_Yscaled_pot, Manova_Yscaled_pot)
  
  Manta_YnormMinMax_pot <- mean(Manta_YnormMinMax_pB < alpha)
  Manova_YnormMinMax_pot <- mean(Manova_YnormMinMax_pB < alpha)
  YnormMinMax_Manta_Manova_pot <- cbind(Manta_YnormMinMax_pot, Manova_YnormMinMax_pot)
  
  # DF que almacena la potencia de cada método:
  DF_pot <- data.frame(rbind(Y_Manta_Manova_pot, Ylog_Manta_Manova_pot, Ysqrt_Manta_Manova_pot,
                             Yscaled_Manta_Manova_pot, YnormMinMax_Manta_Manova_pot))
  colnames(DF_pot) <- c("MANTA", "MANOVA")
  row.names(DF_pot) <- c("Datos sin transformar", "Transformación logarítmica",
                         "Transformación raíz cuadrada", "Escalado de datos por desviación típica",
                         "Normalización Min- Max")
}


# # # # # # # # # # # #
# # # # # # # # # # # #
# #  SIMULACIONES   # #
# # # # # # # # # # # #
# # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PRIMER OBJETIVO #                                                                               #
# # # # # # # # # #                                                                               #
#                                                                                                 #
# Estudiar la pérdida de potencia de la versión asintótica de PERMANOVA (MANTA) con respecto      #
# a MANOVA y otros métodos, profundizando en la afectación de la variación del nivel α de         #
# significación  considerado sobre la potencia de cada uno.                                       #
#                                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# => Se comparará la potencia de ambos métodos según la variación de la correlación (c) tanto para
#    el conjunto de datos simulados sin transformar como bajo las transformaciones consideradas
#    (log, sqrt, escalado de datos por desviación típica y normalización Min- Max).

# => Potencia del método ≡ para las S simulaciones consideradas, representa la fracción de los
#    p-valores del factor a estudio que se encuentran por debajo del nivel α definido.


# Multivariate normal distribution
# --------------------------------

modelSim = "mvnorm" # No Multinomial distribution

# => Se usará la función de simulación de datos "Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor)".
# => Se simularán todas las combinaciones Δ - Cor.
# => Donde Δ = 0 evalua H0 y Δ > 0 evaluará H1.
# => La opción dejar fijo Δ pero variar n eleva demasiado el tiempo de computación. <= HACER LA PRUEBA MIDIÉNDOLO!!!

# Variables fijadas, a priori, no deberían influir en los resultados:
S = 1E3 # Máx. 1E4 si se quiere aumentar la precisión sin elevar el tiempo de computación.
n = 300 # Equilibra la representatividad de los resultados con el tiempo de computación.
q = 3
a = 2
b = 3
u = 1
hk = 1
Var = "equal"

# Variables a estudio, potencialmente influyentes en los resultados:
delta_values <- c(seq(0, 1, 0.1)) # delta
Cor_values <- c(seq(0, .8, 0.2)) # Cor


for (n in N){
  for (q in Q){
    CompMantaManovafunc()
  }
  }





## Código tras la reunión con Diego (17/11/23):

#### Código principal para simular:

# Multivariate normal distribution
# --------------------------------

modelSim = "mvnorm" # No Multinomial distribution

# Variables fijadas, a priori, no deberían influir en los resultados:
q = 3
n = 300
a = 2
b = 3
u = 1
hk = 1

# Variables a estudio, potencialmente influyentes en los resultados:
S = 1000 
delta = 0.1
Cor = 0
Var = "equal"


# Simplex
# -------

modelSim = "simplex"

# Variables fijadas, a priori, no deberían influir en los resultados:
q = 3
n = 300
a = 2
b = 3
u = 1
hk = 1
pdist = "norm"

# Variables a estudio, potencialmente influyentes en los resultados:
S = 1000 
delta = 0.1
Cor = 0
Var = "equal"
loc = 1
stdev = 0.025 # Varía según q, loc y pdist


# Multinomial distribution
# ------------------------

modelSim = "multinom"

# Variables fijadas, a priori, no deberían influir en los resultados:
q = 3
n = 300
a = 2
b = 3
u = 1
hk = 1
pdist = "norm"

# Variables a estudio, potencialmente influyentes en los resultados:
S = 1000 
delta = 0.1
Cor = 0
Var = "equal"
loc = 1
stdev = 0.025 # Varía según q, loc y pdist
lambda = 1000


CompMantaManovafunc <- function() {
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  #  Asignación de la "stdev" para el modelo "simplex"  #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  # Para: m = modelSim = "simplex" 
  # Se carga la tabla correspondiente de "stdev":
  if (modelSim == "simplex") {
    
    tbl <- read.table(sprintf("%s/qlocstdev.%s.tsv", fx, p_dist), h = T)
    colnames(tbl) <- c("Q", "L", "S")
    tbl_name <- "simplex_table_Q_L_S"
    write.csv(tbl, file = paste(tbl_name, ".csv", sep = ""), row.names=FALSE)
    assign(tbl, tbl, envir = .GlobalEnv)
    
    if (! q %in% unique(tbl$Q)) {
      stop(sprintf("stdev not precomputed for q = %s", q))} 
    
    # Se busca la "stdev" que corresponde a los valores [q, loc] seleccionados.
    stdev <- subset(tbl, Q == q & L == loc)$S
    
    # Se crea y guarda un DF con los valores [Q, L, S] que cumplen la condición,
    # asignándole un nombre que refleja los parámetros seleccionados:
    stdevDF <- data.frame(subset(tbl, Q == q & L == loc)$Q,
                          subset(tbl, Q == q & L == loc)$L,
                          subset(tbl, Q == q & L == loc)$S)
    colnames(stdevDF) <- c("q (Q)", "loc (L)", "stdev (S)")
    stdev_name = paste("stdev_", modelSim, sep = "")
    assign(stdev_name, stdevDF[[3]], envir = .GlobalEnv)
    
    # Se asignan los parámetros generados en esta sección:
    Parameter_gen3_List <- list("tbl" = tbl)
    list2env(Parameter_gen3_List, .GlobalEnv)}
  
  
  # # # # # # # # # # # # # # # # # # # # # #
  #  Almacenaje de los datos de simulación  #
  # # # # # # # # # # # # # # # # # # # # # #
  
  # Se crea el directorio de la simulación vigente:
  work_dir <- getwd()
  results_dir <- "Resultados"
  sim_dir_name = paste("[model = ", m, ", a = ",  a, ", b = ",  b, ", n = ",  n,
                       ", u = ",  u, ", w = ",  w, ", S = ",  S, ", q (Q) = ", q,
                       ", loc (L) = ",  loc, "]", sep = "")
  sim_path = file.path(work_dir, results_dir, sim_dir_name)
  dir.create(sim_path, recursive = TRUE, showWarnings = FALSE)
  assign("sim_path", sim_path, envir = .GlobalEnv)
  
  
  # Guardamos la información necesaria hasta este punto:
  if (modelSim == "simplex") {
    
    write.csv(stdevDF, file = file.path(sim_path, paste(stdev_name, ".csv", sep = "")), row.names=FALSE)
    
  }
  
  pdf(file = file.path(sim_path, "FactorLevelsHeatmap.pdf"), width = 8, height = 8)
  pheatmap::pheatmap(cbind(A,B), cluster_cols = F, cluster_rows = F)
  dev.off()
  
  
  # # # # # # # # # # # # # # # # # # # #
  #  Simulación del conjunto de datos   #
  #  y aplicación de MANTA y MANOVA     #
  # # # # # # # # # # # # # # # # # # # #
  
  # Se inicializan los vectores, de tamaño igual al número de simulaciones S, que acumularán
  # las p del factor B para cada método y conjunto de datos considerado:
  Manta_pB = rep(NA, S)
  Manova_pB = rep(NA, S)
  Manta_Ylog_pB = rep(NA, S)
  Manova_Ylog_pB = rep(NA, S)
  Manta_Ysqrt_pB = rep(NA, S)
  Manova_Ysqrt_pB = rep(NA, S)
  Manta_Yscaled_pB = rep(NA, S)
  Manova_Yscaled_pB = rep(NA, S)
  Manta_YnormMinMax_pB = rep(NA, S)
  Manova_YnormMinMax_pB = rep(NA, S)
  
  for (i in 1:S) {
    if (modelSim == "simplex") {
      if (i == 1) { # Sanity check
        tol <- 0.01
        check <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = T, pdist)
        wm <- which.max(check$exp)
        if (abs(check[wm, "obs"] - check[wm, "exp"]) > tol) {
          stop("Deviation from expected centroid greater than tolerance.")
        }
      }
      Y <- Sim.simplex(ch, q, n, loc, delta, hk, stdev, check = F, pdist)
    } else if (modelSim == "mvnorm") {
      Y <- Sim.mvnorm(ch, q, n, mu = rep(10, q), delta, hk, Var, Cor) # Cambiar 0 <-> núm grande (para transformaciones)
    } else if (modelSim == "multinom") {
      Y <- Sim.multinom(ch, q, n, Lambda, delta, loc)
    }
    
    # Transformaciones posibles de los datos simulados (Y):
    Ylog <- log(Y) # Transformación logarítmica.
    Ysqrt <- sqrt(Y) # Transformación sqrt.
    Yscaled <- cbind(Y[,1]/sd(Y[,1]), Y[,2]/sd(Y[,2]), Y[,3]/sd(Y[,3])) # Escalado de datos: se divide cada valor por la desviación típica del conjunto de datos.
    YnormMinMax <-  cbind((Y[,1]-min(Y[,1]))/(max(Y[,1])-min(Y[,1])), # Normalización Min- Max: si queremos que cada variable contribuya por igual al análisis.
                          (Y[,2]-min(Y[,2]))/(max(Y[,2])-min(Y[,2])),
                          (Y[,3]-min(Y[,3]))/(max(Y[,3])-min(Y[,3])))
    
    # Aplicación de Manta y Manova para cada conjunto de datos (Y, Ylog, Ysqrt):
    Manta_pB[i] <- manta(Y ~ A + B + A:B)$aov.tab[2,6]
    Manova_pB[i] <- summary(manova(Y ~ A + B + A:B))$stats[2,6]
    
    Manta_Ylog_pB[i] <- manta(Ylog ~ A + B + A:B)$aov.tab[2,6]
    Manova_Ylog_pB[i] <- summary(manova(Ylog ~ A + B + A:B))$stats[2,6]
    
    Manta_Ysqrt_pB[i] <- manta(Ysqrt ~ A + B + A:B)$aov.tab[2,6]
    Manova_Ysqrt_pB[i] <- summary(manova(Ysqrt ~ A + B + A:B))$stats[2,6]
    
    Manta_Yscaled_pB[i] <- manta(Yscaled ~ A + B + A:B)$aov.tab[2,6]
    Manova_Yscaled_pB[i] <- summary(manova(Yscaled ~ A + B + A:B))$stats[2,6]
    
    Manta_YnormMinMax_pB[i] <- manta(YnormMinMax ~ A + B + A:B)$aov.tab[2,6]
    Manova_YnormMinMax_pB[i] <- summary(manova(YnormMinMax ~ A + B + A:B))$stats[2,6]
  }
  
  # Establecemos el nivel de significación α:
  alpha = 0.05
  
  # Se calcula la potencia de cada método para el α dado:
  Manta_pot <- mean(Manta_pB < alpha)
  Manova_pot <- mean(Manova_pB < alpha)
  Y_Manta_Manova_pot <- cbind(Manta_pot, Manova_pot)
  
  Manta_Ylog_pot <- mean(Manta_Ylog_pB < alpha)
  Manova_Ylog_pot <- mean(Manova_Ylog_pB < alpha)
  Ylog_Manta_Manova_pot <- cbind(Manta_Ylog_pot, Manova_Ylog_pot)
  
  Manta_Ysqrt_pot <- mean(Manta_Ysqrt_pB < alpha)
  Manova_Ysqrt_pot <- mean(Manova_Ysqrt_pB < alpha)
  Ysqrt_Manta_Manova_pot <- cbind(Manta_Ysqrt_pot, Manova_Ysqrt_pot)
  
  Manta_Yscaled_pot <- mean(Manta_Yscaled_pB < alpha)
  Manova_Yscaled_pot <- mean(Manova_Yscaled_pB < alpha)
  Yscaled_Manta_Manova_pot <- cbind(Manta_Yscaled_pot, Manova_Yscaled_pot)
  
  Manta_YnormMinMax_pot <- mean(Manta_YnormMinMax_pB < alpha)
  Manova_YnormMinMax_pot <- mean(Manova_YnormMinMax_pB < alpha)
  YnormMinMax_Manta_Manova_pot <- cbind(Manta_YnormMinMax_pot, Manova_YnormMinMax_pot)
  
  # DF que almacena la potencia de cada método:
  DF_pot <- data.frame(rbind(Y_Manta_Manova_pot, Ylog_Manta_Manova_pot, Ysqrt_Manta_Manova_pot,
                             Yscaled_Manta_Manova_pot, YnormMinMax_Manta_Manova_pot))
  colnames(DF_pot) <- c("MANTA", "MANOVA")
  row.names(DF_pot) <- c("Datos sin transformar", "Transformación logarítmica",
                         "Transformación raíz cuadrada", "Escalado de datos por desviación típica",
                         "Normalización Min- Max")
}




























df = rbind(c("n", "100"), c("q", "2,5"), c("a", "2:3"))
rownames(df) <- df[,1]
df = df[,2, drop = F]
for (n in N){
  for (q in Q){

  }
}

#### Caso particular interesante <=> Forzar matrices de covarianza específicas:

X = matrix(rnorm(q^2), q, q)
crossprod(X)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


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
