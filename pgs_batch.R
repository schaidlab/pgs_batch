"
Usage:
  pgs_batch.R batch (--n_batches=<n_batches> | --n_per_batch=<n_per_batch>) [--dir=<dir> --force]
  pgs_batch.R download --batch_id=<batch_id> [--dir=<dir> --target_build=<target_build> --resume]
  pgs_batch.R create_samplesheet --id=<id> --genos_path_prefix=<genos_path_prefix> --format=<format> [--dir=<dir> --genos_single_file]
  pgs_batch.R calc --id=<id> --batch_id=<batch_id> [--dir=<dir> --profile=<profile> --target_build=<target_build> --min_overlap=<min_overlap> --max_cpus=<max_cpus> --max_memory=<max_memory> --custom_config --resume --extra_args=<extra_args> --offline]
  pgs_batch.R (-h | --help)
  pgs_batch.R --version

Options:
  -h --help
  --version
  --dir=<dir>                             Working directory. If null, then current.
  --id=<id>                               Analysis ID, f.x., name of cohort.
  --n_batches=<n_batches>                 Split scoring files into n_batches number of batches.
  --n_per_batch=<n_per_batch>             Split into batches that have n_per_batch scoring files.
  --force
  --batch_id=<batch_id>                   Run for specific batch.
  --target_build=<target_build>           Genome build [default: GRCh38].
  --resume                                Resume if something fails.
  --genos_path_prefix=<genos_path_prefix> Genotype path prefix. Assumes one file per chromosome ending on the chromosome number. Otherwise, use genos_single_file flag (not recommended, slow)
  --format=<format>                       Genotype format: vcf, bfile (plink 1), or pfile (plink 2).
  --genos_single_file
  --profile=<profile>                     docker, singularity or conda [default: singularity].
  --min_overlap=<min_overlap>             [default: 0].
  --max_cpus=<max_cpus>                   [default: 4].
  --max_memory=<max_memory>               [default: 16.GB].
  --extra_args=<extra_args>               Specify arbitrary pgsc_calc parameters. Details: https://pgsc-calc.readthedocs.io/en/latest/reference/params.html.
  --offline                               Use if working in an offline environment. Make sure to download containers first (see docs).

" -> doc

library(docopt)
a <- docopt(doc)

if (is.null(a$dir)) {
    a$dir <- getwd()
}

nextflow <- stringr::str_glue("{a$dir}/nextflow")
pgsc_calc <- "pgscatalog/pgsc_calc -r v2.0.0-alpha.5"
if (a$offline) {
    if (!fs::dir_exists(stringr::str_glue("{a$dir}/pgsc_calc-2.0.0-alpha.5"))) {
        stop("pgsc_calc-2.0.0-alpha.5 directory doesn't exist. Follow the docs (under 'Offline')")
    }
    nxf_sc_export <- stringr::str_glue("export NXF_SINGULARITY_CACHEDIR={a$dir}/pgsc_calc-2.0.0-alpha.5/nxf_sc")
    pgsc_calc <- stringr::str_glue("{a$dir}/pgsc_calc-2.0.0-alpha.5/main.nf")
}
pgsc_calc_metadata <- stringr::str_glue("{a$dir}/pgs_all_metadata_scores_20240510.csv")

run_batch <- function(a) {
    metadata <- data.table::fread(pgsc_calc_metadata)
    pgs_ids <- metadata[[1]]

    dir_batches <- stringr::str_glue("{a$dir}/batches")

    if (length(fs::dir_ls(dir_batches)) > 0) {
        if (a$force) {
            fs::dir_delete(dir_batches)
            fs::dir_create(dir_batches)
        } else {
            stop(stringr::str_glue("Do you want to overwrite the files in {dir_batches}? Use --force to overwrite"))
        }
    }

    cat(stringr::str_glue("Batching... Outputting in {dir_batches}"), "\n", sep = "")

    if (!is.null(a$n_per_batch)) {
        batches <- split(pgs_ids, ceiling(seq_along(pgs_ids) / as.numeric(a$n_per_batch)))
    } else {
        batches <- split(pgs_ids, cut(seq_along(pgs_ids), as.numeric(a$n_batches), labels = FALSE))
    }

    for (i in 1:length(batches)) {
        readr::write_lines(batches[[i]], stringr::str_glue("{dir_batches}/batch{i}"))
    }
}

get_api_paths <- function(pgs_ids, target_build) {
    pgs_paths <- stringr::str_glue("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_ids}/ScoringFiles/Harmonized/{pgs_ids}_hmPOS_{a$target_build}.txt.gz")
}

run_download <- function(a) {
    allowed_target_builds <- c("GRCh37", "GRCh38")
    if (!a$target_build %in% allowed_target_builds) {
        stop("--target_build has to be paste(allowed_target_builds, collapse = '/')")
    }
    batch_id <- a$batch_id
    resume <- a$resume

    dir_scoringfiles <- stringr::str_glue("{a$dir}/scoringfiles/batch{a$batch_id}")
    fs::dir_create(dir_scoringfiles)

    batch_path <- stringr::str_glue("{a$dir}/batches/batch{batch_id}")
    batch <- readr::read_lines(batch_path)
    pgs_paths <- get_api_paths(batch, a$target_build)
    destfiles <- stringr::str_glue("{dir_scoringfiles}/{fs::path_file(pgs_paths)}")
    curl::multi_download(pgs_paths, destfiles = destfiles, resume = resume)
}

create_samplesheet <- function(a) {
    allowed_formats <- c("vcf", "bfile", "pfile")
    if (!a$format %in% allowed_formats) {
        stop("--format has to be paste(allowed_formats, collapse = '/')")
    }
    if (a$genos_single_file) {
        path_prefix <- a$genos_path_prefix
        chrom <- NA
    } else {
        path_prefix <- stringr::str_glue("{a$genos_path_prefix}{1:22}")
        chrom <- 1:22
    }

    samplesheet <- data.frame(sampleset = a$id, path_prefix = path_prefix, chrom = chrom, format = a$format)

    out_path <- stringr::str_glue("{a$dir}/samplesheet_{a$id}.csv")
    cat(stringr::str_glue("Writing samplesheet file to {out_path}"), "\n", sep = "")
    data.table::fwrite(samplesheet, out_path, sep = ",")
}

run_calc <- function(a) {
    pgsc_calc_input <- stringr::str_glue("{a$dir}/samplesheet_{a$id}.csv")
    pgsc_calc_scores <- stringr::str_glue("--scorefile \"{a$dir}/scoringfiles/batch{a$batch_id}/*{a$target_build}.txt.gz\"")

    pgsc_calc_custom_config <- stringr::str_glue(" -c {a$dir}/custom.config")
    pgsc_calc_resume <- ifelse(a$resume, " -resume", "")
    pgsc_calc_extra_args <- ifelse(!is.null(a$extra_args), stringr::str_glue(" {a$extra_args}"), "")

    dir_runs <- stringr::str_glue("{a$dir}/runs/{a$id}/batch{a$batch_id}")
    dir_results <- stringr::str_glue("{a$dir}/results/{a$id}/batch{a$batch_id}")
    fs::dir_create(dir_runs)

    setwd(dir_runs)

    if (a$offline) {
        nxf_plugins_pin <- stringr::str_glue(" -plugins nf-validation@1.1.3")
    }
    cmd_pgsc_calc <- stringr::str_glue("{nextflow} run {pgsc_calc}{nxf_plugins_pin} -profile {a$profile} --input {pgsc_calc_input} {pgsc_calc_scores} --target_build {a$target_build} --outdir {dir_results} --min_overlap {a$min_overlap} --fast_match --parallel --max_cpus {a$max_cpus} --max_memory {a$max_memory}{pgsc_calc_custom_config}{pgsc_calc_resume}{pgsc_calc_extra_args}")
    if (a$offline) {
        nxf_home_export <- stringr::str_glue("export NXF_HOME={a$dir}/.nextflow")
        cmd_pgsc_calc <- stringr::str_glue("{nxf_home_export}; {nxf_sc_export}; {cmd_pgsc_calc}")
    }

    system(cmd_pgsc_calc)
}

if (a$batch) {
    run_batch(a)
}
if (a$download) {
    run_download(a)
}
if (a$create_samplesheet) {
    create_samplesheet(a)
}
if (a$calc) {
    run_calc(a)
}
