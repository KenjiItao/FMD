#' ---
#' title: "FMD modelling"
#' author: "F.Omata"
#' date: "2017/10/06"
#' ---
#'
#'
#'knitr::opts_chunk$set(echo = TRUE)
#'
#'
#' # Plots for FMD Modelling
#'
#' - Analysis and plot for FMD
#' - Please run this file with the needed libraries installed and a script "multiplot.R"
#' - To install package run: `install.packages("package_name")`
#' - This script will make a directory called "dir" to save all plots
#'
#'
rm(list=ls())
dir.create(file.path("./figs"), showWarnings = FALSE)

require(deSolve)
require(ggplot2)
require(dplyr)
require(tidyverse)
library(latex2exp)

source("./multiplot.R")
#'
#'
#' ## Initial parameters
#'
#' latent period cattle: 3.59
#' subclinical period cattle: 2.04
#' clinical period cattle: 2.35
#'
#' latent period swine: 3.07
#' subclinical period swine: 2.27
#' clinical period swine: 3.42
#'
#' R_0: 38.4 *Need verification for parameter*
#' cattle to cattle: cattle to swine: swine to cattle: swine to swine = 1:0.77:2.45:3.01
#'
#' from: http://www.maff.go.jp/j/wpaper/w_maff/h22/pdf/z_1_1_5.pdf
#' 平成22年度 食料・農業・農村白書（平成23年5月31日公表）
#' 発生農場は 292 農場、発生自治体数は 11 市町、家畜への被害は牛 69,454 頭、
#' 豚 227,949 頭、その他（山羊、羊、イノシシ、水牛等）405 頭
#'
#'
## For cattle
latent_period_cow <- 3.59
subclinical_period_cow <- 2.04
clinical_period_cow <- 2.35 * 20  # most interersting parameter

## For swine
latent_period_swine <- 3.07
subclinical_period_swine <- 2.27
clinical_period_swine <- 3.42 * 20  # most interersting parameter

time_undiscovered <- 22.0

## image_type = "pdf"
image_type = "png"
## image_type = "eps"
image_font_size = 10
image_line_size = 0.5

ode_solve <- function(y, times, func, parms) {
    return(rk(y = y, times = times, func = func, parms = parms, method = rkMethod("rk4"), hini = 0.01))
}

## day_depop_start is days between first confirmation and start of depopulation
params_func <- function(beta., day_depop_start) {
    params <- list(
        beta_cowtocow = beta.,
        beta_cowtoswine = 0.77 * beta.,
        beta_swinetocow = 2.45 * beta.,
        beta_swinetoswine = 3.01 * beta.,
        lambda_cow = 1 / latent_period_cow,
        lambda_swine = 1 / latent_period_swine,
        gamma_cow = 1 / subclinical_period_cow,
        gamma_swine = 1 / subclinical_period_swine,
        ## recover_cow = 1 / clinical_period_cow,
        ## recover_swine = 1 / clinical_period_swine,
        depopulation_capacity = 7,
        depopulation_initial_capacity = 1,
        day_depopulation_start = day_depop_start)  ## need confirmation
    return(params)
}

save_image <- function(path, plt) {
    ## ggsave(filename = path, plot = plt, dpi = 300, width = 10.0, height = 8.0)
    ggsave(filename = path, plot = plt, dpi = 300, width = 85, height = 60, units = "mm")
}
#'
#'
#' ### Initial conditions
#' - From Itao san's email from 2017/10/03
#' estimated onset: 04/17
#' depopulation completion date: 04/21
#' - onset-day:=4+subclinical period  -> Depopulation starts 4 days after the emergence of I1.
#' - Initial values: I1_cow=1, S_cow=8,000, S_swine=500
#' - Run with different $\beta$, finding the one with similar outcome
#'
#'
initial <- c(S_c = 11030, E_c = 0, I1_c = 0, I2_c = 0, S_s = 654, E_s = 0, I1_s = 1, I2_s = 0)

#'
#'
#' ### equation of model
#'
#'
diff <- function(time, state, params) {
    with(as.list(c(state, params)), {
        ## initial の I1_c = 1 なので最初の感染は time = 0.
        ## time = time_undiscovered になると depopulation_initial_capacity で屠殺開始
        ## time = day_depop_start + time_undiscovered になると depopulation_capacity で本格的に屠殺開始
        if (exists("depopulation_data")) {
            day = floor(time)
            val <- depopulation_data[depopulation_data["days"] == day, "completion.of.depopulation"]
            if (length(val) == 0) {
                depopulation_capacity <- 0
            } else {
                depopulation_capacity <- val
            }
        }
        if (time > day_depopulation_start + time_undiscovered) {
            delta_swine <- min(I2_s, depopulation_capacity)
            delta_cow <- min(I2_c, depopulation_capacity - delta_swine)
        } else if (time > time_undiscovered) {
            delta_swine <- min(I2_s, depopulation_initial_capacity)
            delta_cow <- min(I2_c, depopulation_initial_capacity - delta_swine)
        } else {
            delta_swine <- 0
            delta_cow <- 0
        }
        cow_to_cow <- beta_cowtocow * (I1_c + I2_c) * S_c
        swine_to_cow <- beta_swinetocow * (I1_s + I2_s) * S_c
        ec_to_i1c <- lambda_cow * E_c
        i1c_to_i2c <- gamma_cow * I1_c
        cow_to_swine <- beta_cowtoswine * (I1_c + I2_c) * S_s
        swine_to_swine <- beta_swinetoswine * (I1_s + I2_s) * S_s
        es_to_i1s <- lambda_swine * E_s
        i1s_to_i2s <- gamma_swine * I1_s

        S_c_new <- -cow_to_cow - swine_to_cow
        E_c_new <- cow_to_cow + swine_to_cow - ec_to_i1c
        I1_c_new <- ec_to_i1c - i1c_to_i2c
        I2_c_new <- i1c_to_i2c - delta_cow
        S_s_new <- -cow_to_swine - swine_to_swine
        E_s_new <- cow_to_swine + swine_to_swine - es_to_i1s
        I1_s_new <- es_to_i1s - i1s_to_i2s
        I2_s_new <- i1s_to_i2s - delta_swine

        val <- c(S_c = S_c_new, E_c = E_c_new, I1_c = I1_c_new, I2_c = I2_c_new,
                 S_s = S_s_new, E_s = E_s_new, I1_s = I1_s_new, I2_s = I2_s_new)
        if (exists("I2_c_incidence")) {
            val <- c(val, c(I2_c_incidence = i1c_to_i2c))
        }
        if (exists("I2_s_incidence")) {
            val <- c(val, c(I2_s_incidence = i1s_to_i2s))
        }
        if (exists("R_c")) {
            val <- c(val, c(R_c = delta_cow))
        }
        if (exists("R_s")) {
            val <- c(val, c(R_s = delta_swine))
        }
        return(list(val))
    })
}

#'
#'
#' ## Analysis 1: Estimation of $\beta$
#' - How does the difference of $\beta$ contribute to final size?
#'
#' from: http://lin.alic.go.jp/alic/month/domefore/2011/aug/spe-02.htm
#' 宮崎県における口蹄疫からの復興への取り組み
#' "4月20日の発生から130日を要した口蹄疫の防疫措置は完了し、"
#'
#' - Thus max_day=130 may be appropriate.
#' - max_day=150 from Itao-san's advice
#'
#'
final_size_est <- function(params. = params_tmp, initial. = initial, max_day = 150, func. = diff) {
    params. <- c(params.)
    days <- 0:max_day
    initial_state <- c(initial., c(R_c = 0.0, R_s = 0.0))
    column_names <- names(initial_state)
    final_size <- NULL
    result_last_value <- NULL
    while (TRUE) {
        result <- rbind(result_last_value, ode_solve(y = initial_state, times = days, func = func., parms = params.))
        nrow <- nrow(result)
        variation <- result[2:nrow, "R_c"] + result[2:nrow, "R_s"] - (result[1:(nrow - 1), "R_c"] + result[1:(nrow - 1), "R_s"])
        if (tail(variation, 1) < 1.0e-6) {
            for (i in rev(1:length(variation))) {
                if (variation[i] >= 1.0e-6) {
                    max_day <- days[i]
                    final_size <- as.numeric(result[i, "R_c"] + result[i, "R_s"])
                    break
                }
            }
            break
        }
        days <- days + max_day
        ## print(sprintf("beta_cowtocow=%.8f variation=%.8f", params.[["beta_cowtocow"]], variation))
        initial_state <- as.numeric(result[nrow, 2:ncol(result)])
        names(initial_state) <- column_names
        result_last_value <- result[ceiling(nrow / 2):nrow,]
    }
    if (is.null(final_size)) {
        print("Can not calculate final size")
        print(result)
        quit()
    }
    return(final_size)
}

plot_epicurve <- function(result, name_prefix, depop = NULL) {
    max_y <- ceiling(max(result["I2_c"]) / 10.0) * 10
    graph_min_y <- -max_y * 0.02
    graph_max_y <- max_y * 1.02
    plt <- ggplot(result, aes_string("time"))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.89, 0.8),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(1, 4, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    if (!is.null(depop)) {
        plt <- plt + geom_segment(aes(x = depop, xend = depop, y = graph_min_y, yend = graph_max_y), size = image_line_size)
        plt <- plt + annotate("text", label = "Depopulation start", x = depop + 27, y = max_y * 0.95, size = 3)
    }
    ## plt <- plt + geom_line(aes_string(y = "S_c", linetype = shQuote("S_c")), size = image_line_size)
    plt <- plt + geom_line(aes_string(y = "E_c", linetype = shQuote("E_c")), size = image_line_size)
    plt <- plt + geom_line(aes_string(y = "I1_c", linetype = shQuote("I1_c")), size = image_line_size)
    plt <- plt + geom_line(aes_string(y = "I2_c", linetype = shQuote("I2_c")), size = image_line_size)
    plt <- plt + scale_linetype_manual(values = c("twodash", "dotted", "dashed"), labels = lapply(c("$E_c$", "$I_{1c}$", "$I_{2c}$"), TeX))
    plt <- plt + scale_x_continuous(limits = c(0, 120), expand = c(0, 0))
    plt <- plt + scale_y_continuous(expand = c(0, 0), limits = c(graph_min_y, graph_max_y))
    plt <- plt + xlab("Day") + ylab("Number")
    save_image(paste0("./figs/", name_prefix, "_cattle.", image_type), plt)

    max_y <- ceiling(max(result["I2_s"]) / 5.0) * 5
    graph_min_y <- -max_y * 0.02
    graph_max_y <- max_y * 1.02
    plt <- ggplot(result, aes_string("time"))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.89, 0.8),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(1, 4, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    if (!is.null(depop)) {
        plt <- plt + geom_segment(aes(x = depop, xend = depop, y = graph_min_y, yend = graph_max_y), size = image_line_size)
        plt <- plt + annotate("text", label = "Depopulation start", x = depop + 27, y = max_y * 0.95, size = 3)
    }
    ## plt <- plt + geom_line(aes_string(y = "S_s", linetype = shQuote("S_s")), size = image_line_size)
    plt <- plt + geom_line(aes_string(y = "E_s", linetype = shQuote("E_s")), size = image_line_size)
    plt <- plt + geom_line(aes_string(y = "I1_s", linetype = shQuote("I1_s")), size = image_line_size)
    plt <- plt + geom_line(aes_string(y = "I2_s", linetype = shQuote("I2_s")), size = image_line_size)
    plt <- plt + scale_linetype_manual(values = c("twodash", "dotted", "dashed"), labels = lapply(c("$E_s$", "$I_{1s}$", "$I_{2s}$"), TeX))
    plt <- plt + scale_x_continuous(limits = c(0, 120), expand = c(0, 0))
    plt <- plt + scale_y_continuous(expand = c(0, 0), limits = c(graph_min_y, graph_max_y))
    plt <- plt + xlab("Day") + ylab("Number")
    save_image(paste0("./figs/", name_prefix, "_swine.", image_type), plt)
}

#'
#'
#' ### power of $\beta$
#' $\beta * 10^{-3}$ - $\beta * 10^3$
#'
#'
run_analysis1_1 <- function() {
    pows <- seq(-6, -3, 0.02)
    final_size_list <- c()
    for (p in pows) {
        params_tmp <- params_func(10^p, 21)
        final_size_list <- append(final_size_list, final_size_est(params_tmp, initial, 150, diff))
    }

    df_final_size <- cbind.data.frame(pows, final_size_list)

    plt <- ggplot(aes(x = pows, y = final_size_list), data = df_final_size)
    plt <- plt + geom_line(size = image_line_size)
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       plot.margin = unit(c(5, 5, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + scale_x_continuous(expand = c(0, 0))
    plt <- plt + xlab("Power of beta: log10") + ylab("Final size")
    ## show(plt)

    save_image(paste0("./figs/fig1-1_final_size_by_power_beta.", image_type), plt)
}

#'
#'
#' - $\beta * 10^p$
#' - $\beta*10^{-0.20}$: `r beta * 10^(-0.20)` -> `r dplyr::filter(df_final_size, pows == "-0.2")$final_size_list`
#' - $\beta*10^{-0.18}$: `r beta * 10^(-0.18)` -> `r dplyr::filter(df_final_size, pows == "-0.18")$final_size_list`
#'
#' ###  Change of final size when $\beta$ is proportionaly changed
#'
#'
#'
run_analysis1_2 <- function() {
    beta_low <- 1.2e-5
    beta_high <- 1.5e-5
    betas <- seq(beta_low, beta_high, length.out=1000)
    final_size_list <- c()

    for (b in betas) {
        params_tmp <- params_func(b, 21)
        final_size_list <- append(final_size_list, final_size_est(params_tmp, initial, 150, diff))
    }

    df_final_size <- cbind.data.frame(betas, final_size_list)
    print(df_final_size)

    plt <- ggplot(aes(x = betas, y = final_size_list), data = df_final_size)
    plt <- plt + geom_line(size = image_line_size)
    ## final size of actual outbreak: needs confirmation
    plt <- plt + geom_segment(x = 0, xend = 40, y = 291, yend = 16, size = image_line_size)
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       plot.margin = unit(c(5, 5, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + scale_x_continuous(expand = c(0, 0))
    ## plt <- plt + scale_y_continuous(expand = c(0, 0))
    plt <- plt + xlab("beta") + ylab("Final size")
    ## show(plt)

    save_image(paste0("./figs/fig1-2_linear_beta_vs_finalsize2.", image_type), plt)

    #'
    #'
    #' - Result
    #'
    #'
    print(df_final_size[800:900, ])
}

#'
#'
#' ## Analysis 2: Final size simulation
#' - Initiation day of depopulationを連続的に変化させて、final sizeの振る舞いを調べた。
#'
#' from: http://lin.alic.go.jp/alic/month/domefore/2011/aug/spe-02.htm
#' 宮崎県における口蹄疫からの復興への取り組み
#' "4月20日の発生から130日を要した口蹄疫の防疫措置は完了し、"
#'
#'
run_analysis2 <- function(beta) {
    max_depop_date <- 40
    path_data <- "result_FigFinalSizeVSd.csv"
    if(file.exists(path_data)) {
        df_final_size <- read.csv(path_data)
    } else {
        depop_date_without_undiscovered <- 0:max_depop_date
        final_size_list <- c()
        for (d in depop_date_without_undiscovered) {
            params_tmp <- params_func(beta, d)
            final_size_list <- append(final_size_list, final_size_est(params_tmp, initial, 300, diff))
        }
        df_final_size <- cbind.data.frame(depop_date_without_undiscovered, final_size_list)
        write.csv(df_final_size, path_data)
    }

    max_y <- ceiling(max(df_final_size) / 1000.0) * 1000
    ## graph_min_y <- -max_y * 0.02
    ## graph_max_y <- max_y * 1.02
    graph_min_y <- 0.0
    graph_max_y <- 12000.0

    plt <- ggplot(aes(x = depop_date_without_undiscovered, y = final_size_list), data = df_final_size)
    plt <- plt + geom_line(size = image_line_size)
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       plot.margin = unit(c(2, 3, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + scale_x_continuous(limits = c(0, max_depop_date), expand = c(0, 0))
    plt <- plt + scale_y_continuous(limits = c(graph_min_y, graph_max_y), expand = c(0, 0), breaks = c(0, 4000, 8000, 12000))
    plt <- plt + xlab("Initiation of depopulation\n(d, days after index case)") + ylab("Final size")
    image_path <- paste0("./figs/FigFinalSizeVSd.", image_type)
    save_image(image_path, plt)
    print(df_final_size[20:25, ])
}
#'
#'
#' ## Analysis 3: Sensitivity analysis of beta
#' - Final size simulation was done under different beta values.
#' - $0.6-1.4 *\beta$ was used for simulation.
#'
#'
## This process may take some time
run_analysis3 <- function(beta) {
    mults <- seq(0.6, 1.4, 0.05)
    max_depop_date <- 100
    depop_date_without_undiscovered <- 0:max_depop_date
    df_depop_mult <- data.frame()
    for (m in mults) {
        beta_new <- m * beta
        final_size_list <- c()
        for (d in depop_date_without_undiscovered) {
            params_tmp <- params_func(beta_new, d)
            final_size_list <- append(final_size_list, final_size_est(params_tmp, initial, 150, diff))
        }
        mult <- rep(m, length(final_size_list))
        df_depop <- cbind.data.frame(depop_date_without_undiscovered, final_size_list, mult)
        df_depop_mult <- rbind(df_depop_mult, df_depop)
    }

    plt <- ggplot(aes(x = depop_date_without_undiscovered, y = final_size_list, group = factor(mult)), data = df_depop_mult)
    plt <- plt + geom_line(aes(colour = mult), size = image_line_size)
    plt <- plt + scale_colour_gradient2(low="blue", mid="gray", high="red", midpoint=1.0)
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       plot.margin = unit(c(5, 5, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + xlab("Initiation day of depopulation") + ylab("Final size")
    image_path <- paste0("./figs/fig3_final_size_by_beta_and_initiation_day.", image_type)
    save_image(image_path, plt)
}

run_analysis3_2 <- function(beta) {
    mults <- seq(0.7, 1.3, 0.1)
    beta_labels <- sprintf("beta%02d", seq(1, length(mults)))
    max_depop_date <- 160
    depop_date_without_undiscovered <- 0:max_depop_date
    path_data <- "result_FigFinalSizeVariousBeta.csv"
    if(file.exists(path_data)) {
        df_depop_mult <- read.csv(path_data)
    } else {
        df_depop_mult <- data.frame(depop_date = depop_date_without_undiscovered)
        for (m in mults) {
            beta_new <- m * beta
            final_size_list <- c()
            for (d in depop_date_without_undiscovered) {
                params_tmp <- params_func(beta_new, d)
                final_size_list <- append(final_size_list, final_size_est(params_tmp, initial, 150, diff))
            }
            df_depop <- data.frame(final_size_list)
            df_depop_mult <- cbind(df_depop_mult, df_depop)
        }
        colnames(df_depop_mult) <- c("depop_date_without_undiscovered", beta_labels)
        write.csv(df_depop_mult, path_data)
    }


    max_y <- max(df_depop_mult)
    ## graph_min_y <- -max_y * 0.02
    ## graph_max_y <- max_y * 1.02
    graph_min_y <- -500.0
    graph_max_y <- 12000.0
    graph_max_x <- 100.0
    image_font_size <- 8

    plt <- ggplot(aes(depop_date_without_undiscovered), data = df_depop_mult)
    for (name in beta_labels) {
        plt <- plt + geom_line(aes_string(y = name, linetype = shQuote(name)), size = image_line_size)
    }
    plt <- plt + scale_linetype_discrete(labels = lapply(sprintf("$\\beta = %.2e$", mults * beta), TeX))
    plt <- plt + scale_x_continuous(limits = c(0.0, graph_max_x), expand = c(0, 0))
    plt <- plt + scale_y_continuous(limits = c(graph_min_y, graph_max_y), expand = c(0, 0), breaks = c(0, 4000, 8000, 12000))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.box.margin = unit(c(0, 0, 0, 0), "mm"),
                       legend.box.spacing = unit(c(0, 0, 0, 0), "mm"),
                       legend.key.size = grid::unit(2.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(2, 1, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + xlab("Initiation of depopulation\n(d, days after index case)") + ylab("Final size")
    image_path <- paste0("./figs/FigFinalSizeVariousBeta.", image_type)
    save_image(image_path, plt)
}

#'
#' ## Analysis 4: Equational interpretation of this model
#'
run_analysis4 <- function(beta, output_name) {
    path_data <- paste0("result_", output_name, ".csv")
    max_day <- 300
    days <- 0:(max_day + time_undiscovered)
    params = params_func(beta, max_day) # Consider case so that we don't depopulate
    if (file.exists(path_data)) {
        result_ed <- read.csv(path_data)
    } else {
        result <- ode_solve(y = initial, times = days, func = diff, parms = params)
        result <- as.data.frame(result)
        result[, "time"] <- result[, "time"] - time_undiscovered

        result_ed <- result %>%
            dplyr::mutate(FOI = params$beta_cowtocow * (I1_c + I2_c) * S_c +
                              params$beta_swinetocow * (I1_s + I2_s) * S_c +
                              params$beta_cowtoswine * (I1_c + I2_c) * S_s +
                              params$beta_swinetoswine * (I1_s + I2_s) * S_s,
                          capacity = rep(params[["depopulation_capacity"]], nrow(result)),
                          init_capacity = rep(1.0, nrow(result)))
        write.csv(result_ed, path_data)
    }

    graph_max_day <- max_day
    max_y <- NA
    for (i in 1:nrow(result_ed)) {
        row <- result_ed[i,]
        if (row[["FOI"]] > params[["depopulation_capacity"]]) {
            graph_max_day <- ceiling((row[["time"]] * 1.5) / 10.0) * 10
            val <- result_ed[min(nrow(result_ed), ceiling(i + graph_max_day - row[["time"]])), "FOI"]
            max_y <- ceiling(val * 1.1 / 10.0) * 10
            break
        }
    }

    graph_max_y <- 15.0 * 1.02
    plt <- ggplot(aes(time), data = result_ed)
    plt <- plt + geom_line(aes(y = FOI, linetype = "1"), size = image_line_size)
    plt <- plt + geom_line(aes(y = capacity, linetype = "2"), size = image_line_size)
    plt <- plt + geom_line(aes(y = init_capacity, linetype = "3"), size = image_line_size)
    plt <- plt + scale_linetype_manual(values = c("solid", "twodash", "dashed"), labels = c("Increasing rate of infection", "Capacity of full-fledged depopulatation", "Capacity of initial depopulatation"))
    plt <- plt + scale_x_continuous(limits = c(-time_undiscovered, 30.0), expand = c(0, 0))
    plt <- plt + scale_y_continuous(limits = c(0, graph_max_y), expand = c(0, 0), breaks = c(0, 1, 7, 15))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = ceiling(image_font_size * 0.9)),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = ceiling(image_font_size * 0.9)),
                       legend.title = element_blank(),
                       legend.position = c(0.48, 0.85),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.key.width = unit(2.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(1, 3, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + xlab("Day") + ylab("Rate")
    save_image(paste0("./figs/", output_name, ".", image_type), plt)

    print(result_ed[0:(graph_max_day + time_undiscovered), c("time", "FOI")])
}

run_analysis4_2 <- function(beta) {
    max_day <- 300
    days <- 0:(max_day + time_undiscovered)
    params = params_func(beta, max_day) # Consider case so that we don't depopulate
    initial_val <- c(initial, c(I2_c_incidence = 0.0, I2_s_incidence = 0.0))
    result <- ode_solve(y = initial_val, times = days, func = diff, parms = params)
    result <- as.data.frame(result)
    result[, "time"] <- result[, "time"] - time_undiscovered
    n <- nrow(result)
    s <- as.vector(result[, "I2_c_incidence"] + result[, "I2_s_incidence"])
    new_cases <- c(0.0, s[2:n] - s[1:(n - 1)])
    result <- cbind(result, data.frame(FOI = new_cases, capacity = rep(params[["depopulation_capacity"]], n)))

    graph_max_day <- max_day
    max_y <- NA
    for (i in 1:nrow(result)) {
        row <- result[i,]
        if (row[["FOI"]] > params[["depopulation_capacity"]]) {
            graph_max_day <- ceiling((row[["time"]] * 1.5) / 10.0) * 10
            val <- result[min(nrow(result), ceiling(i + graph_max_day - row[["time"]])), "FOI"]
            max_y <- ceiling(val * 1.1 / 10.0) * 10
            break
        }
    }

    graph_max_y <- max_y * 1.02
    plt <- ggplot(aes(time), data = result)
    plt <- plt + geom_line(aes(y = FOI, linetype = "Daily new case"), size = image_line_size)
    plt <- plt + geom_line(aes(y = capacity, linetype = "Depopulation capacity"), size = image_line_size)
    plt <- plt + scale_x_continuous(limits = c(0, graph_max_day), expand = c(0, 0))
    plt <- plt + scale_y_continuous(limits = c(0, graph_max_y), expand = c(0, 0))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.32, 0.88),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(1, 1, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + xlab("Day") + ylab("Cases")
    save_image(paste0("./figs/fig4_actual_new_cases.", image_type), plt)
}
#'
#'
#' Analysis 4: Sensitivity analysis of 牧場間の感染のパラメータ
#' - Final size simulation was done under different 牧場間の感染のパラメータ values.
#' - How to implement the change parameters?
#' - No need to implement.

## Create epi curve
run_analysis5 <- function(beta) {
    max_day <- 120
    depop_start <- 21
    days <- 0:(max_day + time_undiscovered)

    params = params_func(beta, depop_start)
    initial_val <- c(initial, c(R_c = 0.0, R_s = 0.0))
    result <- ode_solve(y = initial_val, times = days, func = diff, parms = params)
    result[, "time"] <- result[, "time"] - time_undiscovered
    plot_epicurve(as.data.frame(result), "epicurve_depop_21", params[["day_depopulation_start"]])
}

## 新規感染数 (I2 の増分)
run_analysis6 <- function(beta) {
    max_day <- 120
    depop_start <- 21
    days <- 0:(max_day + time_undiscovered + 1)
    params = params_func(beta, depop_start)
    initial_val <- c(initial, c(I2_c_incidence = 0.0, I2_s_incidence = 0.0))
    result <- ode_solve(y=initial_val, times=days, func=diff, parms=params)
    ind <- 1:(nrow(result) - 1)
    I2_c_incidence_cum <- result[, "I2_c_incidence"]
    I2_s_incidence_cum <- result[, "I2_s_incidence"]
    I2_c_new <- I2_c_incidence_cum[ind + 1] - I2_c_incidence_cum[ind]
    I2_s_new <- I2_s_incidence_cum[ind + 1] - I2_s_incidence_cum[ind]
    incidence <- data.frame(I2_c_new = I2_c_new, I2_s_new = I2_s_new, I2_total_new = I2_c_new + I2_s_new)
    result_ed <- cbind(result[ind, ], incidence)
    result_ed[, "time"] <- result_ed[, "time"] - time_undiscovered

    max_y <- ceiling(max(incidence[, "I2_total_new"]) / 5.0) * 5
    graph_min_y <- -max_y * 0.02
    graph_max_y <- max_y * 1.02

    ## plt <- ggplot(aes(time, FOI), data = result_ed)
    plt <- ggplot(aes(time), data = result_ed)
    plt <- plt + geom_line(aes(y = I2_total_new, linetype = "Daily new I2 total"), size = image_line_size)
    plt <- plt + geom_line(aes(y = I2_c_new, linetype = "Daily new I2_c"), size = image_line_size)
    plt <- plt + geom_line(aes(y = I2_s_new, linetype = "Daily new I2_s"), size = image_line_size)
    plt <- plt + scale_linetype_manual(values = c("twodash", "dotted", "dashed"), labels = lapply(c("Daily new $I_{2c} + I_{2s}$", "Daily new $I_{2c}", "Daily new $I_{2s}$"), TeX))

    plt <- plt + scale_x_continuous(limits = c(0, max_day), expand = c(0, 0))
    plt <- plt + scale_y_continuous(limits = c(graph_min_y, graph_max_y), expand = c(0, 0))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.73, 0.85),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(1, 4, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + xlab("Day") + ylab("Number")
    save_image(paste0("./figs/fig6_new_infection.", image_type), plt)
}

calc_incidence <- function(beta, depop_start, max_day, data_start, depop_data = NULL) {
    params <- params_func(beta, depop_start)
    if (!is.null(depop_data)) {
        params[["depopulation_data"]] <- depop_data
    }
    initial_mle <- c(initial, c(I2_c_incidence = 0.0, I2_s_incidence = 0.0))
    result <- ode_solve(y = initial_mle, times = 0:max_day, func = diff, parms = params)
    ind <- data_start:(max_day - 1)
    sim_reported <- (result[ind + 1, "I2_c_incidence"] - result[ind, "I2_c_incidence"]) + (result[ind + 1, "I2_s_incidence"] - result[ind, "I2_s_incidence"])
    return(sim_reported)
}

## Estimate beta when depopulation is given by step functions
nll_estimate_beta <- function(beta) {
    print(beta)
    depop_start <- 21
    ## データの day がシミュレーションの day のいつに対応するのか決める必要がある
    data_start <- 22
    data <- read.csv("epi_depop_curve.csv")
    data_reported <- data[, "reported"]
    max_day <- length(data_reported) + data_start + 1
    sim_reported <- calc_incidence(beta, depop_start, max_day, data_start)
    val <- -sum(dpois(data_reported, sim_reported, log = TRUE))
    print(val)
    return(val)
}

## Estimate beta when depopulation is given by actual data
nll_estimate_beta_actual_depop <- function(beta) {
    print(beta)
    depop_start <- 21
    data_start <- 22
    data <- read.csv("epi_depop_curve.csv")
    data[, "days"] <- data[, "days"] + data_start
    data_reported <- data[, "reported"]
    max_day <- length(data_reported) + data_start + 1
    sim_reported <- calc_incidence(beta, depop_start, max_day, data_start, data)
    val <- -sum(dpois(data_reported, sim_reported, log = TRUE))
    print(val)
    return(val)
}

run_analysis7 <- function(beta_initial) {
    estimation <- optim(beta_initial, nll_estimate_beta, gr = NULL, method = "BFGS", control = list(parscale = beta_initial))
    print(estimation)
    beta_estimate <- estimation$par
    print(calc_incidence(beta_estimate, 21, 100, 22))
    params <- params_func(beta_estimate, 21)
    print(final_size_est(params))
}

run_analysis8 <- function(beta_initial) {
    estimation <- optim(beta_initial, nll_estimate_beta_actual_depop, gr = NULL, method = "BFGS", control = list(parscale = beta_initial))
    print(estimation)
    beta_estimate <- estimation$par
    print(calc_incidence(beta_estimate, 21, 100, 22))
    params <- params_func(beta_estimate, 21)
    print(final_size_est(params))
}

create_actual_epicurve <- function() {
    data <- read.csv("epi_depop_curve_orig.csv")
    depop_model <- c(rep(1.0, 22), rep(7.0, 55))
    data <- cbind(data, data.frame(depop_model = depop_model))

    plt <- ggplot(data, aes_string("days"))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.8, 0.9),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(2, 1, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + geom_line(aes_string(y = "completion.of.depopulation", linetype = shQuote("Depopulated")), size = image_line_size, color = "grey")
    plt <- plt + geom_line(aes_string(y = "reported", linetype = shQuote("Reported")), size = image_line_size)
    plt <- plt + scale_x_continuous(expand = c(0, 0))
    plt <- plt + scale_y_continuous(expand = c(0, 0), limits = c(0, 16), breaks = c(0.0, 5.0, 10.0, 15.0))
    plt <- plt + xlab("Day") + ylab("Number of farms")
    save_image(paste0("./figs/FigDepopEpicurve.", image_type), plt)


    plt <- ggplot(data, aes_string("days"))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.6, 0.9),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(2, 1, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + geom_line(aes_string(y = "completion.of.depopulation", linetype = shQuote("Actual depopulation")), size = image_line_size, color = "grey")
    plt <- plt + geom_line(aes_string(y = "depop_model", linetype = shQuote("Depopulation capacity in model")), size = image_line_size)
    plt <- plt + scale_x_continuous(expand = c(0, 0))
    plt <- plt + scale_y_continuous(expand = c(0, 0), limits = c(0, 16), breaks = c(0.0, 5.0, 10.0, 15.0))
    plt <- plt + xlab("Day") + ylab("Number of farms")
    save_image(paste0("./figs/FigDepopModel.", image_type), plt)
}

create_graph_model_depopulation <- function() {
    data <- read.csv("epi_depop_curve_orig.csv")
    depop_model <- c(rep(1.0, 22), rep(7.0, 55))
    data <- cbind(data, data.frame(depop_model = depop_model))

    plt <- ggplot()
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = ceiling(image_font_size * 0.9)),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = ceiling(image_font_size * 0.7)),
                       legend.title = element_blank(),
                       legend.direction = "horizontal",
                       legend.position = c(0.5, 1.06),
                       legend.key.size = unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       legend.box.margin = unit(c(0, 0, 0, 0), "mm"),
                       legend.box.spacing = unit(c(0, 0, 0, 0), "mm"),
                       plot.margin = unit(c(5, 1, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())

    xmax <- 76
    ymin <- 0.0
    ymax <- 14.0

    init_depop_start <- 0
    depop_start <- 21.0
    init_depop_cap <- 1.0
    depop_cap <- 7.0

    plt <- plt + scale_x_continuous(limits = c(-time_undiscovered, xmax), expand = c(0, 0), breaks = c(-time_undiscovered, 0.0, init_depop_start, depop_start), labels = c(-time_undiscovered, 0.0, init_depop_start, sprintf("d = %d", depop_start)))
    plt <- plt + scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax), breaks = c(0.0, 1.0, 7.0, 14.0))
    plt <- plt + xlab("Day") + ylab("Depopulation capacity")

    plt <- plt + geom_segment(aes(x = 0.0, xend = 0.0, y = ymin, yend = ymax), size = image_line_size, color = "#cccccc", linetype = "dashed")
    plt <- plt + geom_segment(aes(x = init_depop_start, xend = init_depop_start, y = ymin, yend = ymax), size = image_line_size, color = "grey", linetype = "dashed")
    plt <- plt + geom_segment(aes(x = depop_start, xend = depop_start, y = ymin, yend = ymax), size = image_line_size, color = "grey", linetype = "dashed")

    plt <- plt + annotate("text", label = "Undiscovered", x = -time_undiscovered / 2.0, y = ymax * 0.92, size = 2.5)
    plt <- plt + annotate("text", label = "Initial\ndepopulation", x = (init_depop_start + depop_start) / 2.0, y = ymax * 0.9, size = 2.5)
    plt <- plt + annotate("text", label = "Full-fledged depopulation", x = (depop_start + xmax) / 2.0, y = ymax * 0.92, size = 2.5)

    plt <- plt + geom_line(aes_string(x = data[, "days"], y = data[, "completion.of.depopulation"], linetype = shQuote("Actual depopulation")), size = image_line_size, color = "grey")
    plt <- plt + geom_line(aes_string(x = c(0.0, init_depop_start, init_depop_start, depop_start, depop_start, xmax), y = c(0.0, 0.0, init_depop_cap, init_depop_cap, depop_cap, depop_cap), linetype = shQuote("Depopulation capacity in model")), size = image_line_size)

    save_image(paste0("./figs/FigDepopSetting.", image_type), plt)
}

create_data_fitting <- function(beta) {
    data <- read.csv("epi_depop_curve_orig.csv")[c("days", "reported")]

    depop_start <- 21
    data_start <- 22
    max_day <- nrow(data) + data_start + 1
    sim_reported <- calc_incidence(beta, depop_start, max_day, 1)

    data <- rbind(data.frame(days = -data_start:-1, reported = rep(0, data_start)), data)
    data <- cbind(data, data.frame(model = sim_reported))

    plt <- ggplot(data, aes_string("days"))
    plt <- plt + theme_bw()
    plt <- plt + theme(axis.text = element_text(size = image_font_size),
                       axis.title = element_text(size = image_font_size, face = "bold"),
                       text = element_text(size = image_font_size),
                       legend.text = element_text(size = image_font_size),
                       legend.title = element_blank(),
                       legend.position = c(0.85, 0.88),
                       legend.key.size = grid::unit(1.0, "lines"),
                       legend.key.height = unit(1.0, "lines"),
                       legend.background = element_rect(fill = alpha("white", 0.0)),
                       plot.margin = unit(c(2, 2, 1, 1), "mm"),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    plt <- plt + geom_line(aes_string(y = "reported", linetype = shQuote("Reported")), size = image_line_size, color = "grey")
    plt <- plt + geom_line(aes_string(y = "model", linetype = shQuote("Model")), size = image_line_size)
    plt <- plt + scale_x_continuous(expand = c(0, 0))
    plt <- plt + scale_y_continuous(expand = c(0, 0), limits = c(0, 16), breaks = c(0.0, 5.0, 10.0, 15.0))
    plt <- plt + xlab("Day") + ylab("Number of farms")
    save_image(paste0("./figs/FigEpicurveDataFitting.", image_type), plt)
}

## These two functions do not support new settings
## run_analysis1_1()
## run_analysis1_2()

## run_analysis3(beta)
## run_analysis4_2(beta)
## run_analysis5(beta)
## run_analysis6(beta)

beta <- 1.057085e-05
run_analysis2(beta)
run_analysis3_2(beta)
run_analysis4(beta, "FigIncreasingRateOfInfection")
run_analysis4(beta * 0.7, "FigIncreasingRateOfInfectionSmallBeta")
create_actual_epicurve()
create_graph_model_depopulation()
create_data_fitting(beta)

## Parameter estimation
## beta <- 1.054251e-05
## run_analysis7(beta)                     # 1.057085e-05
## run_analysis8(beta)                     # 1.027725e-05
