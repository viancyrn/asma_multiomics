
# Proyecto: Asma multiomics
# Dataset: GSE201872
# Objetivo= descargar y explorar los metadatos de GEO

# 1. Instalar GEOquery si todavía no está instalado
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("GEOquery")
}

# 2. Cargar la librería
library(GEOquery)

# 3. Descargar los metadatos del dataset GSE201872
gse <- getGEO(
  "GSE201872",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# 4. Comprobar cuántas plataformas contiene el dataset
length(gse)

# 5. Identificar las plataformas disponibles
names(gse)

# 6. Comprobar el número de muestras por plataforma
sapply(gse, ncol)

# 7. Extraer los metadatos de ambas plataformas
metadata_list <- lapply(gse, Biobase::pData)

# 8. Añadir el identificador de plataforma
metadata_list <- Map(
  function(df, platform) {
    df$platform <- platform
    df
  },
  metadata_list,
  c("GPL13534", "GPL21145")
)

# 9. Comprobar el número de muestras por plataforma
sapply(metadata_list, nrow)

# 10. Identificar las columnas clínicas disponibles
lapply(metadata_list, function(df) {
  grep("characteristics", colnames(df), value = TRUE)
})

# 11. Explorar los nombres de las variables
colnames(metadata_list[[1]])

# 12. Identificar las columnas que contienen características clínicas
clinical_columns <- lapply(metadata_list, function(df) {
  grep("characteristics", colnames(df), value = TRUE)
})

# Mostrar las columnas clínicas de ambas plataformas
clinical_columns


# 13. Mostrar las características clínicas de las primeras cinco muestras
head(
  metadata_list[[1]][, clinical_columns[[1]], drop = FALSE],
  5
)

# 14. Convertir las características clínicas en columnas separadas
parse_metadata <- function(df) {
  
  # Seleccionar las columnas de características clínicas
  char_columns <- grep(
    "^characteristics_ch1",
    colnames(df),
    value = TRUE
  )
  
  # Extraer los nombres de las variables
  variable_names <- sub(
    ":.*$",
    "",
    as.character(df[1, char_columns])
  )
  
  # Extraer los valores después de los dos puntos
  clinical_data <- lapply(df[char_columns], function(x) {
    trimws(sub("^[^:]*:", "", as.character(x)))
  })
  
  # Convertir los valores en un dataframe
  clinical_data <- as.data.frame(
    clinical_data,
    stringsAsFactors = FALSE
  )
  
  # Asignar los nombres de las variables
  colnames(clinical_data) <- variable_names
  
  # Conservar los identificadores y la plataforma
  clinical_data$geo_accession <- df$geo_accession
  clinical_data$platform <- df$platform
  
  return(clinical_data)
}

# 15. Aplicar la función a ambas plataformas
clinical_list <- lapply(metadata_list, parse_metadata)

# 16. Comprobar las dimensiones
lapply(clinical_list, dim)

# 17. Mostrar las primeras muestras
head(clinical_list[[1]])

# 18. Comprobar si ambas plataformas tienen las mismas variables
identical(
  names(clinical_list[[1]]),
  names(clinical_list[[2]])
)

# 19. Contar los grupos de asma en cada plataforma
lapply(clinical_list, function(df) {
  table(df$asthma, useNA = "ifany")
})

# 21. Combinar los metadatos clínicos de ambas plataformas
clinical_metadata <- dplyr::bind_rows(clinical_list)

# 22. Comprobar las dimensiones de la tabla combinada
dim(clinical_metadata)

# 23. Verificar la distribución de asma por plataforma
table(
  clinical_metadata$platform,
  clinical_metadata$asthma
)

# 24. Comprobar si existen identificadores GSM duplicados
sum(duplicated(clinical_metadata$geo_accession))

# 25. Definir las variables clínicas numéricas
numeric_columns <- c(
  "age", "step", "ige", "bmi",
  "fvc_abs", "fev1_abs", "fev1_pct",
  "feno", "bal_neut", "bal_eo",
  "blood_eo_abs",
  "ancestry_pc1", "ancestry_pc2", "ancestry_pc3"
)

# 26. Convertir las variables numéricas y reconocer los valores faltantes
clinical_metadata[numeric_columns] <- lapply(
  clinical_metadata[numeric_columns],
  function(x) {
    as.numeric(dplyr::na_if(x, "NA"))
  }
)

# 27. Comprobar los tipos de datos
sapply(clinical_metadata[numeric_columns], class)

# 28. Comprobar los valores faltantes
colSums(is.na(clinical_metadata[numeric_columns]))

# 29. Recuperar los valores originales de ambas plataformas
clinical_raw <- dplyr::bind_rows(clinical_list)

# 30. Identificar los valores que no se pueden convertir a números
for (column in numeric_columns) {
  
  original <- clinical_raw[[column]]
  converted <- suppressWarnings(as.numeric(original))
  
  invalid_values <- unique(
    original[is.na(converted) & !is.na(original)]
  )
  
  if (length(invalid_values) > 0) {
    cat("\nVariable:", column, "\n")
    print(invalid_values)
  }
}

# 31. Resumir la edad y el IMC por plataforma y grupo de asma
clinical_metadata |>
  dplyr::group_by(platform, asthma) |>
  dplyr::summarise(
    n = dplyr::n(),
    age_median = median(age, na.rm = TRUE),
    age_min = min(age, na.rm = TRUE),
    age_max = max(age, na.rm = TRUE),
    bmi_median = median(bmi, na.rm = TRUE),
    .groups = "drop"
  )

# 32. Examinar la distribución del sexo por plataforma y asma
table(
  clinical_metadata$platform,
  clinical_metadata$asthma,
  clinical_metadata$Sex
)

# 33. Examinar la distribución del tabaquismo actual
table(
  clinical_metadata$platform,
  clinical_metadata$asthma,
  clinical_metadata$current_smok
)

# 34. Examinar la distribución de los lotes técnicos
table(
  clinical_metadata$platform,
  clinical_metadata$asthma,
  clinical_metadata$batch
)

# 35. Examinar la procedencia de los participantes
table(
  clinical_metadata$platform,
  clinical_metadata$asthma,
  clinical_metadata[["recruitment source"]]
)

# 36. Filtrar candidatos sin obesidad y sin tabaquismo actual
eligible_candidates <- clinical_metadata |>
  dplyr::filter(
    platform == "GPL13534",
    bmi < 30,
    current_smok == "No"
  )

# 37. Contar los candidatos según el estado de asma
table(eligible_candidates$asthma)

# 38. Examinar los candidatos según su cohorte de origen
table(
  eligible_candidates[["recruitment source"]],
  eligible_candidates$asthma
)

# 39. Examinar la distribución por sexo
table(
  eligible_candidates$asthma,
  eligible_candidates$Sex
)

# 40. Examinar la distribución por lote técnico
table(
  eligible_candidates$asthma,
  eligible_candidates$batch
)

# 41. Resumir la edad y el IMC de los candidatos elegibles
eligible_candidates |>
  dplyr::group_by(asthma) |>
  dplyr::summarise(
    n = dplyr::n(),
    age_median = median(age),
    age_min = min(age),
    age_max = max(age),
    bmi_median = median(bmi),
    bmi_min = min(bmi),
    bmi_max = max(bmi),
    .groups = "drop"
  )

# 42. Examinar la edad y el IMC según la cohorte de origen
eligible_candidates |>
  dplyr::group_by(`recruitment source`, asthma) |>
  dplyr::summarise(
    n = dplyr::n(),
    age_median = median(age),
    bmi_median = median(bmi),
    .groups = "drop"
  )

# 43. Crear una tabla de candidatos de ABRIDGE y CABFIR
candidate_overview <- eligible_candidates |>
  dplyr::filter(
    `recruitment source` %in% c("ABRIDGE", "CABFIR")
  ) |>
  dplyr::select(
    geo_accession,
    geo_id,
    `recruitment source`,
    asthma,
    age,
    Sex,
    bmi,
    batch,
    ChIP
  ) |>
  dplyr::arrange(
    `recruitment source`,
    asthma,
    age
  )

# 44. Abrir la tabla de candidatos en RStudio
View(candidate_overview)

# 45. Guardar la tabla de candidatos
dir.create("results", showWarnings = FALSE)

write.csv(
  candidate_overview,
  "results/GSE201872_candidate_overview.csv",
  row.names = FALSE
)

# Crear figuros de análisis exploratorio de datos
# 1. Instalar ggplot2 si todavía no está instalado
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# 2. Preparar los datos para la primera figura
platform_counts <- clinical_metadata |>
  dplyr::count(platform, asthma)

# 3. Crear la figura de distribución por plataforma
fig_platform <- ggplot2::ggplot(
  platform_counts,
  ggplot2::aes(
    x = platform,
    y = n,
    fill = asthma
  )
) +
  ggplot2::geom_col(
    position = "dodge",
    width = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_dodge(width = 0.65),
    vjust = -0.4,
    size = 3
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asthma" = "#0072B2",
      "Control" = "#E69F00"
    )
  ) +
  ggplot2::labs(
    title = "DNA methylation samples by platform",
    subtitle = "GSE201872 | Before eligibility filtering",
    x = "Illumina array platform",
    y = "Number of samples",
    fill = "Clinical group"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

# 4. Mostrar la figura
print(fig_platform)

# 5. Guardar la figura en formato PNG
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

ggplot2::ggsave(
  filename = "results/figures/01_platform_distribution.png",
  plot = fig_platform,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 6. Crear la carpeta para las figuras de metilación
dir.create(
  "results/methylation/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# 7. Guardar la primera figura en alta resolución
ggplot2::ggsave(
  filename = "results/methylation/figures/01_platform_distribution.png",
  plot = fig_platform,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 8. Preparar los datos antes y después de aplicar los criterios
filter_counts <- dplyr::bind_rows(
  clinical_metadata |>
    dplyr::filter(platform == "GPL13534") |>
    dplyr::count(asthma) |>
    dplyr::mutate(stage = "Before filtering"),
  eligible_candidates |>
    dplyr::count(asthma) |>
    dplyr::mutate(stage = "After filtering")
)

# 9. Definir el orden de las etapas
filter_counts$stage <- factor(
  filter_counts$stage,
  levels = c("Before filtering", "After filtering")
)

# 10. Crear la figura de selección de candidatos
fig_filtering <- ggplot2::ggplot(
  filter_counts,
  ggplot2::aes(x = stage, y = n, fill = asthma)
) +
  ggplot2::geom_col(
    position = "dodge",
    width = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_dodge(width = 0.65),
    vjust = -0.4,
    size = 3
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asthma" = "#0072B2",
      "Control" = "#E69F00"
    )
  ) +
  ggplot2::labs(
    title = "Eligibility filtering of methylation samples",
    subtitle = "GPL13534 | BMI < 30 and current non-smoker",
    x = NULL,
    y = "Number of samples",
    fill = "Clinical group"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

# 11. Guardar la segunda figura
ggplot2::ggsave(
  "results/methylation/figures/02_eligibility_filtering.png",
  plot = fig_filtering,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 12. Contar los candidatos e incluir combinaciones con cero muestras
recruitment_counts <- eligible_candidates |>
  dplyr::count(`recruitment source`, asthma) |>
  tidyr::complete(
    `recruitment source`,
    asthma = c("Asthma", "Control"),
    fill = list(n = 0)
  )
 
# 13. Crear la figura de distribución por cohorte
fig_recruitment <- ggplot2::ggplot(
  recruitment_counts,
  ggplot2::aes(
    x = `recruitment source`,
    y = n,
    fill = asthma
  )
) +
  ggplot2::geom_col(
    position = "dodge",
    width = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_dodge(width = 0.65),
    vjust = -0.4,
    size = 2
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asthma" = "#0072B2",
      "Control" = "#E69F00"
    )
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title = "Eligible candidates by recruitment source",
    subtitle = "GPL13534 | BMI < 30 and current non-smoker",
    x = "Recruitment source",
    y = "Number of candidates",
    fill = "Clinical group"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

# 15. Guardar la tercera figura
ggplot2::ggsave(
  filename = "results/methylation/figures/03_recruitment_source.png",
  plot = fig_recruitment,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 16. Seleccionar los candidatos de ABRIDGE y CABFIR
age_candidates <- eligible_candidates |>
  dplyr::filter(
    `recruitment source` %in% c("ABRIDGE", "CABFIR")
  )

# 17. Crear la figura de distribución de edad
fig_age <- ggplot2::ggplot(
  age_candidates,
  ggplot2::aes(
    x = asthma,
    y = age,
    fill = asthma
  )
) +
  ggplot2::stat_boxplot(
    geom = "errorbar",
    width = 0.2
  ) +
  ggplot2::geom_boxplot(
    width = 0.5,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    size = 1,
    alpha = 0.8
  ) +
  ggplot2::facet_wrap(
    ~ `recruitment source`
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asthma" = "#0072B2",
      "Control" = "#E69F00"
    )
  ) +
  ggplot2::labs(
    title = "Age distribution of eligible candidates",
    subtitle = "GPL13534 | BMI < 30 and current non-smoker",
    x = "Clinical group",
    y = "Age (years)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "none"
  )

# 18. Guardar la figura de edad con los extremos de los bigotes
ggplot2::ggsave(
  filename = "results/methylation/figures/04_age_distribution.png",
  plot = fig_age,
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 19. Crear la figura de distribución del IMC
fig_bmi <- ggplot2::ggplot(
  age_candidates,
  ggplot2::aes(
    x = asthma,
    y = bmi,
    fill = asthma
  )
) +
  ggplot2::stat_boxplot(
    geom = "errorbar",
    width = 0.2
  ) +
  ggplot2::geom_boxplot(
    width = 0.5,
    alpha = 0.5,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  ggplot2::facet_wrap(
    ~ `recruitment source`
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asthma" = "#0072B2",
      "Control" = "#E69F00"
    )
  ) +
  ggplot2::labs(
    title = "BMI distribution of eligible candidates",
    subtitle = "GPL13534 | BMI < 30 and current non-smoker",
    x = "Clinical group",
    y = "BMI (kg/m²)"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "none"
  )

# 20. Guardar la figura de distribución del IMC
ggplot2::ggsave(
  filename = "results/methylation/figures/05_bmi_distribution.png",
  plot = fig_bmi,
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)

# 21. Contar los candidatos según cohorte, asma y sexo
sex_counts <- age_candidates |>
  dplyr::count(`recruitment source`, asthma, Sex)

# 22. Crear la figura de distribución por sexo
fig_sex <- ggplot2::ggplot(
  sex_counts,
  ggplot2::aes(
    x = Sex,
    y = n,
    fill = asthma
  )
) +
  ggplot2::geom_col(
    position = "dodge",
    width = 0.65
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = n),
    position = ggplot2::position_dodge(width = 0.65),
    vjust = -0.4,
    size = 3
  ) +
  ggplot2::facet_wrap(
    ~ `recruitment source`
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Asthma" = "#0072B2",
      "Control" = "#E69F00"
    )
  ) +
  ggplot2::scale_y_continuous(
    breaks = 0:12,
    expand = ggplot2::expansion(mult = c(0, 0.15))
  ) +
  ggplot2::labs(
    title = "Sex distribution of eligible candidates",
    subtitle = "GPL13534 | BMI < 30 and current non-smoker",
    x = "Sex",
    y = "Number of candidates",
    fill = "Clinical group"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top",
    panel.grid.major.x = ggplot2::element_blank()
  )

# 26. Guardar la figura de distribución por sexo
ggplot2::ggsave(
  filename = "results/methylation/figures/06_sex_distribution.png",
  plot = fig_sex,
  width = 9,
  height = 5,
  dpi = 300,
  bg = "white"
)

