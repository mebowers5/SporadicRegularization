pwr <- function(output) {
  stats::power.t.test(n = NULL, sd = max(na.omit(output$sd)), 
                      sig.level = 0.05, power = 0.95, delta = 1,
                      type = "paired", alternative = "two.sided")
}