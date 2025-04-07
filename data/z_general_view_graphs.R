library(igraph)

networks <- list.files(
  "data/graphs", pattern = "[0-9]+-(ergm|sf|swp[0-9]{2}|degseq|er)\\.rds$",full.names = TRUE) 


# Assuming 'networks' contains your list of network file paths
# Select one network file (e.g., the first one)
network_file <- networks[1]

# --------------------------------------------------

library(network)

# Load your network object from an .rds file
network_data <- readRDS(network_file)


num_vertices <- network.size(network_data)
is_directed <- is.directed(network_data)
has_hyperedges <- get.network.attribute(network_data, "hyper")
has_loops <- has.loops(network_data)
has_multiple_edges <- get.network.attribute(network_data, "multiple")
is_bipartite <- get.network.attribute(network_data, "bipartite")

cat("Network Attributes:\n")
cat("Vertices:", num_vertices, "\n")
cat("Directed:", is_directed, "\n")
cat("Hyperedges:", has_hyperedges, "\n")
cat("Loops:", has_loops, "\n")
cat("Multiple Edges:", has_multiple_edges, "\n")
cat("Bipartite:", is_bipartite, "\n")

num_edges <- network.edgecount(network_data)
cat("Total Edges:", num_edges, "\n")

vertex_attributes <- list.vertex.attributes(network_data)
edge_attributes <- list.edge.attributes(network_data)

cat("Vertex Attributes:\n", paste(vertex_attributes, collapse=", "), "\n")
cat("Edge Attributes:\n", paste(edge_attributes, collapse=", "), "\n")


# Extract distinct values for each specified vertex attribute
distinct_initialsNum <- unique(get.vertex.attribute(network_data, "initialsNum"))
distinct_lunch <- unique(get.vertex.attribute(network_data, "lunch"))
distinct_na <- unique(get.vertex.attribute(network_data, "na"))
distinct_unique <- unique(get.vertex.attribute(network_data, "unique"))

# Print the distinct values for each attribute
cat("Distinct values for initialsNum:", distinct_initialsNum, "\n")
cat("Distinct values for lunch:", distinct_lunch, "\n")
cat("Distinct values for na:", distinct_na, "\n")
cat("Distinct values for unique:", distinct_unique, "\n")

cat("Distinct values for edge attribute 'na':", unique(get.edge.attribute(network_data, "na", unlist = TRUE)), "\n")
