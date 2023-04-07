output <- function(filename) {
  results <- read.csv(file = filename) %>%
    select(-X) %>%
    group_by(gid) %>%
    summarise(mean = mean(dif),
              sd = sd(dif),
              se = plotrix::std.error(dif)) %>%
    na.omit()
}