pwr <- function(output, sig.level, power, delta) {
  stats::power.t.test(n = NULL, sd = max(na.omit(output$sd)), 
                      sig.level = sig.level, power = power, delta = delta,
                      type = "paired", alternative = "two.sided")
}
