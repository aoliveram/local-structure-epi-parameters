# Generating scale-free networks using igraph
library(igraph)
library(netplot)
library(ggplot2)

nnets <- 1000L
nagents <- 1000L
set.seed(1231333)

# Generating the networks
networks <- replicate(nnets, {
    g <- sample_pa(nagents, m = 2, directed = FALSE)
    list(source = as.integer(ends(g, E(g))[, 1]),
         target = as.integer(ends(g, E(g))[, 2]))
}, simplify = FALSE)


if (interactive()) {
# Reading one to check the result
net <- graph_from_edgelist(as.matrix(as.data.frame(networks[[1]])))

nplot(net, skip.vertex = TRUE)

# Plot the degree sequence of net
# op <- par(logy = TRUE)
table(degree(net)) |> 
    as.data.frame() |>
    ggplot(aes(x = Var1, y = Freq)) +
        geom_point() +
        scale_y_log10() +
        labs(x = "Degree", y = "Frequency")

# Saving the figure
ggsave("data/degree-distribution.pdf", width = 5, height = 5)
}
# hist(, xlab = "Degree", breaks = 100)
# par(op)

# Saving the networks in the data folder
saveRDS(networks, "data/networks-scalefree.rds", compress = FALSE)
