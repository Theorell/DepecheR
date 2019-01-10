# Function to give the right cluster annotation to each event
# for each cluster investigation
dViolinsCoFunction1 <- function(event, n) {
    ifelse(event == n, return(n), return("All clusters"))
}
