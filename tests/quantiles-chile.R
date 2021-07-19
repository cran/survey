library(survey)
datos <- readRDS("datos_ejemplo.rds")

design <- svydesign(id = ~id_directorio, strata = ~estrato, weights = ~f_pers, check.strata = TRUE, data = datos)
set.seed(234262762)
repdesign <- as.svrepdesign(design, type = "subbootstrap", replicates=20)
options(survey.lonely.psu="remove")

values<-datos$ing_t_p[datos$CL_GRUPO_OCU_08=="ISCO08_6"]

suppressWarnings({
f0<-coef(svyquantile(~ing_t_p, subset(design,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="math"))
f0.5<-coef(svyquantile(~ing_t_p, subset(design,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="school"))
})
all.equal(c(values[1],mean(values)), as.vector(c(f0,f0.5)))

suppressWarnings({
f0<-coef(svyquantile(~ing_t_p, subset(repdesign,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="math"))
f0.5<-coef(svyquantile(~ing_t_p, subset(repdesign,CL_GRUPO_OCU_08=="ISCO08_6"),quantiles=c(0.5), qrule="school"))
})
all.equal(c(values[1],mean(values)), as.vector(c(f0,f0.5)))
