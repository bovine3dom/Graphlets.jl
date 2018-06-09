# Motifs

Interface:

`Motif.getmotifs(g,k; norm=true, verbose=false)`

Gets the frequencies of unique subgraphs of size k appearing in G.

`norm` currently just divides the counts of each motif by the total number of motifs appearing in the graph.

Coloured graph support is a work in progress.

# Backends

We currently just use a Julian reimplementation of Kavosh<!-- cite -->, but we are investigating the use of g-tries as we have heard good things about it.

Nauty.jl is used for determining isomorphism.
