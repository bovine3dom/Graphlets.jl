# Motifs

Interface:

`Motif.getmotifs(G,k; norm=true, verbose=false)`

Gets the frequencies of unique subgraphs of size k appearing in G.

`norm` currently just divides the counts of each motif by the total number of motifs appearing in the graph.

Coloured graph support is available in [ColoredGraphs.jl](https://github.com/bovine3dom/ColoredGraphs.jl).

# Backends

We currently just use a Julian reimplementation of [Kavosh](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-318), but we are investigating the use of g-tries as we have heard good things about it.

[Nauty.jl](https://github.com/bovine3dom/Nauty.jl) is used for determining isomorphism.
