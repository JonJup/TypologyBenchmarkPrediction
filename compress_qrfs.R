taxa <- c("diatoms", "fishes", "macroinvertebrates", "macrophytes")
files <- list.files("shiny/models/diatoms")

for (i in taxa){
        cat(sprintf("starting with %s\n", i))
        for (j in files){
                
                cond = file.exists(paste0("models_f_upload/",i, "_",j))
                
                if (cond) {
                        cat(sprintf("     skipping %s\n", j))
                        next
                }
                
                
                cat(sprintf("        compressing %s\n", j))
                x <- readRDS(paste0(c("shiny", "models", i,j), collapse = "/"))
                saveRDS(x, paste0("models_f_upload/",i, "_",j), compress = "xz")
        }
}

