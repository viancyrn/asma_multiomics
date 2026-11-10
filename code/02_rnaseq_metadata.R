# 1. Cargar los paquetes necesarios
library(GEOquery)
library(dplyr)

# 2. Descargar los metadatos de RNA-seq
gse_rnaseq <- getGEO(
  "GSE201955",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# 3. Extraer los metadatos de las muestras
rnaseq_metadata <- bind_rows(
  lapply(gse_rnaseq, Biobase::pData)
)

# 4. Examinar la estructura de los metadatos
dim(rnaseq_metadata)
colnames(rnaseq_metadata)

# 5. Examinar los identificadores de las muestras
head(
  rnaseq_metadata[, c(
    "title",
    "geo_accession",
    "description",
    "relation"
  )],
  10
)

# 6. Identificar los participantes presentes en ambas plataformas
shared_ids <- intersect(
  eligible_candidates$geo_id,
  rnaseq_metadata$description
)

# 7. Contar los participantes compartidos
length(shared_ids)

# 8. Mostrar los primeros identificadores compartidos
head(shared_ids, 10)

# 9. Seleccionar los candidatos con datos de ambas ómicas
multiomics_candidates <- eligible_candidates %>%
  filter(geo_id %in% shared_ids)

# 10. Contar los participantes según su condición de asma
multiomics_candidates %>%
  count(asthma)

# 11. Extraer los identificadores y el estado de asma de RNA-seq
rnaseq_asthma <- rnaseq_metadata %>%
  transmute(
    geo_id = description,
    asthma_rnaseq = `asthma:ch1`
  )

# 12. Comparar el estado de asma entre ambas plataformas
asthma_comparison <- multiomics_candidates %>%
  select(geo_id, asthma) %>%
  left_join(rnaseq_asthma, by = "geo_id")

# 13. Resumir los participantes compartidos según el estado de asma
asthma_comparison %>%
  count(
    asthma,
    asthma_rnaseq,
    name = "n_participants"
  ) %>%
  mutate(
    status_match = asthma == asthma_rnaseq
  )

# 14. Comprobar identificadores duplicados en los candidatos de metilación
multiomics_candidates %>%
  count(geo_id) %>%
  filter(n > 1)

# 15. Comprobar identificadores duplicados en los metadatos de RNA-seq
rnaseq_metadata %>%
  count(description) %>%
  filter(n > 1)

# 16. Examinar la distribución de los participantes compartidos
# según la fuente de reclutamiento y el estado de asma

multiomics_candidates %>%
  count(`recruitment source`, asthma) %>%
  arrange(`recruitment source`, asthma)

# 17. Examinar la relación entre la fuente de reclutamiento y el lote
multiomics_candidates %>%
  count(`recruitment source`, batch, asthma) %>%
  arrange(`recruitment source`, batch, asthma)

# 18. Resumir la edad y el sexo por grupo y fuente de reclutamiento
multiomics_candidates %>%
  filter(`recruitment source` %in% c("ABRIDGE", "CABFIR")) %>%
  group_by(`recruitment source`, asthma) %>%
  summarise(
    n = n(),
    mean_age = round(mean(as.numeric(age), na.rm = TRUE), 1),
    sd_age = round(sd(as.numeric(age), na.rm = TRUE), 1),
    female = sum(Sex == "Female", na.rm = TRUE),
    male = sum(Sex == "Male", na.rm = TRUE),
    .groups = "drop"
  )

# 19. Examinar las características individuales de los candidatos
selection_overview <- multiomics_candidates %>%
  filter(`recruitment source` %in% c("ABRIDGE", "CABFIR")) %>%
  select(
    geo_id,
    `recruitment source`,
    batch,
    asthma,
    age,
    Sex,
    bmi
  ) %>%
  arrange(`recruitment source`, asthma, age)

View(selection_overview)

# 20. Visualizar las edades disponibles por diagnóstico y cohorte
library(ggplot2)

fig_age_overlap <- ggplot(
  selection_overview,
  aes(
    x = `recruitment source`,
    y = as.numeric(age),
    color = asthma
  )
) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(~ asthma) +
  labs(
    title = "Age distribution of multi-omics candidates",
    subtitle = "Eligible participants with matched DNA methylation and RNA-seq data",
    x = "Recruitment source",
    y = "Age (years)",
    color = "Diagnosis"
  ) +
  theme_minimal(base_size = 12)

print(fig_age_overlap)

# 21. Comparar las edades disponibles por diagnóstico y cohorte
selection_overview %>%
  group_by(asthma, `recruitment source`) %>%
  summarise(
    ages = paste(sort(as.numeric(age)), collapse = ", "),
    .groups = "drop"
  )

# 22. Verificar los criterios de elegibilidad antes de la selección
selection_overview_check <- multiomics_candidates %>%
  filter(`recruitment source` %in% c("ABRIDGE", "CABFIR")) %>%
  summarise(
    total = n(),
    bmi_below_30 = sum(as.numeric(bmi) < 30, na.rm = TRUE),
    no_current_smoking = sum(current_smok == "No", na.rm = TRUE),
    both_criteria = sum(
      as.numeric(bmi) < 30 & current_smok == "No",
      na.rm = TRUE
    )
  )

print(selection_overview_check)

# 23. Seleccionar cinco participantes por diagnóstico y cohorte
# Utilizar una semilla fija para garantizar la reproducibilidad
set.seed(2026)

selected_samples <- multiomics_candidates %>%
  filter(`recruitment source` %in% c("ABRIDGE", "CABFIR")) %>%
  group_by(`recruitment source`, asthma) %>%
  slice_sample(n = 5) %>%
  ungroup()

# 24. Verificar la distribución de los participantes seleccionados
selected_samples %>%
  count(`recruitment source`, asthma)

# 25. Resumir las características de los participantes seleccionados
selection_summary <- selected_samples %>%
  group_by(`recruitment source`, asthma) %>%
  summarise(
    n = n(),
    mean_age = round(mean(as.numeric(age)), 1),
    sd_age = round(sd(as.numeric(age)), 1),
    female = sum(Sex == "Female"),
    male = sum(Sex == "Male"),
    mean_bmi = round(mean(as.numeric(bmi)), 1),
    .groups = "drop"
  )

print(selection_summary)

# 26. Comparar las características globales entre asma y controles
selected_samples %>%
  group_by(asthma) %>%
  summarise(
    n = n(),
    mean_age = round(mean(as.numeric(age)), 1),
    sd_age = round(sd(as.numeric(age)), 1),
    mean_bmi = round(mean(as.numeric(bmi)), 1),
    female = sum(Sex == "Female"),
    male = sum(Sex == "Male"),
    .groups = "drop"
  )

# 27. Mostrar los participantes seleccionados y sus características
selected_samples %>%
  select(
    geo_id,
    `recruitment source`,
    batch,
    asthma,
    age,
    Sex,
    bmi,
    current_smok
  ) %>%
  arrange(`recruitment source`, asthma, age) %>%
  print()

# 28. Comparar las características clínicas entre ambas ómicas
# Armonizar la codificación del tabaquismo: N = No, Y = Yes

rnaseq_clinical <- rnaseq_metadata %>%
  transmute(
    geo_id = description,
    age_rnaseq = as.numeric(`age:ch1`),
    bmi_rnaseq = as.numeric(`bmi:ch1`),
    smoking_rnaseq = recode(
      `current_smok:ch1`,
      "N" = "No",
      "Y" = "Yes"
    )
  )

# Verificar la concordancia clínica entre las plataformas
clinical_comparison <- selected_samples %>%
  select(geo_id, age, bmi, current_smok) %>%
  left_join(rnaseq_clinical, by = "geo_id") %>%
  summarise(
    n = n(),
    age_matches = sum(age == age_rnaseq, na.rm = TRUE),
    bmi_matches = sum(
      near(as.numeric(bmi), bmi_rnaseq),
      na.rm = TRUE
    ),
    smoking_matches = sum(
      current_smok == smoking_rnaseq,
      na.rm = TRUE
    )
  )

print(clinical_comparison)

# 29. Vincular los identificadores GEO de ambas ómicas

rnaseq_accessions <- rnaseq_metadata %>%
  transmute(
    geo_id = description,
    rnaseq_gsm = geo_accession
  )

sample_manifest <- selected_samples %>%
  transmute(
    geo_id,
    methylation_gsm = geo_accession,
    rnaseq_group = asthma,
    recruitment_source = `recruitment source`,
    batch,
    age,
    sex = Sex,
    bmi,
    current_smok
  ) %>%
  left_join(rnaseq_accessions, by = "geo_id") %>%
  relocate(rnaseq_gsm, .after = methylation_gsm)

# Verificar que las veinte muestras tengan ambos identificadores GEO
sample_manifest %>%
  summarise(
    participants = n(),
    methylation_gsm_available = sum(!is.na(methylation_gsm)),
    rnaseq_gsm_available = sum(!is.na(rnaseq_gsm))
  )

# 30. Guardar la selección reproducible de las veinte muestras

dir.create("results", showWarnings = FALSE)

readr::write_csv(
  sample_manifest,
  "results/multiomics_sample_manifest.csv"
)

# Verificar que el archivo se haya guardado correctamente
file.exists("results/multiomics_sample_manifest.csv")

# 31. Leer el archivo conservando el lote como variable categórica
saved_manifest <- readr::read_csv(
  "results/multiomics_sample_manifest.csv",
  col_types = readr::cols(
    batch = readr::col_character(),
    .default = readr::col_guess()
  )
)

# 32. Verificar la integridad de los datos guardados
dim(saved_manifest)

dplyr::n_distinct(saved_manifest$geo_id)

all.equal(
  as.data.frame(saved_manifest),
  as.data.frame(sample_manifest),
  check.attributes = FALSE
)