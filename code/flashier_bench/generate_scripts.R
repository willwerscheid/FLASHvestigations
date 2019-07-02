for (data.name in c("gtex", "trachea", "sp_trachea")) {
  for (prior.family in c("point.normal", "normal.scale.mix")) {
    for (method in c("flash", "flashier")) {
      script <- paste0("data.name <- \"", data.name, "\"")
      script <- c(script, paste0("prior.family <- \"", prior.family, "\""))
      script <- c(script, paste0("source(\"./bench_", method, ".R\")"))
      script <- c(script, "t.init; t.greedy; t.backfit")

      if (prior.family == "point.normal") {
        abbr.prior <- "pn"
      } else {
        abbr.prior <- "smn"
      }
      fname <- paste0("./code/flashier_bench/bench_",
                      data.name, "_", abbr.prior, "_", method, ".R")

      fconn <- file(fname)
      writeLines(script, fconn)
      close(fconn)
    }
  }
}
