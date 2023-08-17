#' Calculate fish osmoregulation index according to fitted cubic function
#'
#' @import ggplot2
#' @param kpa vector of kpa (x)
#' @param MO2 vector of MO2 (y)
#' @param period vector coded to indicate measurements under "SMR"
#' @param plot should a plot of the fitted function be produced?
#'
#' @return list containing fitted model coefficients, Efron's pseudo R2, regulation index, and p50
#' @export
#' @examples
#' \dontrun{
#' # with plot
#' regIndex(data$kpa, data$MO2, data$period)
#' # without plot
#' data.frame(t(unlist(regIndex(data$kpa, data$MO2, data$period, plot = FALSE))))
#' }

regIndex = function(kpa, MO2, period, plot = TRUE) {
  x = func1 = func_min = func_max = area = NULL

  data = data.frame(kpa = kpa, MO2 = MO2, period = period)
  ### nonlinear least squares
  nlsmod = stats::nls(MO2 ~ y0 + a * kpa + b * kpa^2 + c * kpa^3,
               start = list(y0 = 0, a = 0, b = 0, c = 0),
               algorithm = "port", upper = c(0, Inf, Inf, Inf),
               data = data)

  # fitted coefficients
  fitted_coefs = summary(nlsmod)$coefficients[, "Estimate"]

  # find root of fitted equation
  realRoots = function(roots) {Re(roots)[abs(Im(roots)) < 1e-8]}
  lower_bound = realRoots(polyroot(fitted_coefs))
  stopifnot("fitted equation has more than one root. consider using a different function" = length(lower_bound) == 1)

  # Efron's pseudo R-squared
  R2 = round(1 - sum(stats::residuals(nlsmod)^2) /
               sum((data[, as.character(stats::formula(nlsmod)[[2]])] -
                      mean(data[, as.character(stats::formula(nlsmod)[[2]])]))^2), 3)

  # means at SMR
  conf_means = colMeans(data[data$period %in% c(1, "SMR"), c("kpa", "MO2")])

  # p50
  MO2_50 = conf_means[["MO2"]] / 2
  p50 = realRoots(polyroot(fitted_coefs - c(MO2_50, rep(0, length(fitted_coefs) - 1))))[1]

  # fitted cubic function
  rootfun = function(x) {
    summary(nlsmod)$parameters["y0", "Estimate"] +
      summary(nlsmod)$parameters["a", "Estimate"] * x +
      summary(nlsmod)$parameters["b", "Estimate"] * x^2 +
      summary(nlsmod)$parameters["c", "Estimate"] * x^3
  }

  # slope of perfect conformity
  slope = conf_means[["MO2"]] / conf_means[["kpa"]]

  # area between curves
  intfun = function(x) {rootfun(x) - slope * x}
  reg_area = stats::integrate(intfun, lower_bound, conf_means[["kpa"]])$value
  # triangle area
  tri_area = conf_means[["MO2"]] * conf_means[["kpa"]] / 2
  # regulation_index
  reg_ind = round(reg_area / tri_area, 3)

  ### plot
  if (isTRUE(plot)) {
    func_area = data.frame(x = seq(lower_bound, conf_means[["kpa"]], by = 0.01),
                           func1 = sapply(seq(lower_bound, conf_means[["kpa"]], by = 0.01), rootfun),
                           func2 = sapply(seq(lower_bound, conf_means[["kpa"]], by = 0.01), function (x) {slope * x}))
    func_area$func_min = pmin(func_area$func1, func_area$func2)
    func_area$func_max = pmax(func_area$func1, func_area$func2)
    func_area$area = ifelse(func_area$func1 < func_area$func2, "a", "b")

    p = ggplot(func_area, aes(x = x, y = func1)) +
      geom_ribbon(aes(ymin = func_min, ymax = func_max, fill = area), alpha = 0.5,
                  show.legend = FALSE) +
      geom_point(data = data, aes(x = kpa, y = MO2)) +
      geom_smooth(data = data, aes(x = kpa, y = MO2),
                  method = "nls",
                  formula = y ~ y0 + a * x + b * x^2 + c * x^3,
                  method.args = list(start = list(y0 = 0, a = 0, b = 0, c = 0),
                                     algorithm = "port", upper = c(0, Inf, Inf, Inf)),
                  se = FALSE,
                  fullrange = TRUE,
                  na.rm = TRUE,
                  colour = "blue") +
      geom_segment(x = 0, y = 0,
                   xend = conf_means[["kpa"]], yend = conf_means[["MO2"]],
                   linewidth = 1, lty = 2) +
      geom_segment(x = 0, y = conf_means[["MO2"]],
                   xend = conf_means[["kpa"]], yend = conf_means[["MO2"]],
                   lty = 2, linewidth = 1) +
      geom_segment(x = 0, y = MO2_50, xend = p50, yend = MO2_50, colour = "red", linewidth = 1) +
      geom_segment(x = p50, y = 0, xend = p50, yend = MO2_50, colour = "red", linewidth = 1) +
      scale_x_continuous(expand = c(0, 0), labels = scales::comma) +
      scale_y_continuous(expand = c(0, 0), labels = scales::comma) +
      coord_cartesian(xlim = c(0, pretty(max(data$kpa), n = 10)[pretty(max(data$kpa), n = 10) > max(data$kpa)][1]),
                      ylim = c(0, pretty(max(data$MO2), n = 10)[pretty(max(data$MO2), n = 10) > max(data$MO2)][1])) +
      labs(x = "kpa", y = expression(ring(M)*O[2])) +
      theme_classic(base_size = 12)

    list(`coefficients` = round(fitted_coefs, 3), `pseudo R2` = R2,
         `regulation index` = reg_ind, p50 = round(p50, 3), plot = p)
  } else {
    list(`coefficients` = round(fitted_coefs, 3), `pseudo R2` = R2,
         `regulation index` = reg_ind, p50 = round(p50, 3))
  }
}
