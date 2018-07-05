module Motifs
    # TODO:
    #   - give variables sensible names
    #   - give functions sensible names
    #       - follow Julia convention
    #   - remove superfluous type annotaitons
    #   - Multithreading...
    #       - Each root node is independent

    # Notes
    # canonical labelling -> graph: (where nodes is number of nodes in subgraph)
    # a = BitArray(64,[nodes])
    # a.chunks = Nauty.canonical_labelling(g)[1]
    # Array{Int64,2}(a[end-nodes+1:end,:])

    module kavosh
        import LightGraphs
        const lg = LightGraphs
        import ColoredGraphs
        import Nauty
        import IterTools
        const it = IterTools

        struct MotifSig
            canong::Array{UInt64,1}
            partition::Array{Int32,1}
        end
        function Base.hash(m::MotifSig)
            return hash((m.canong,m.partition))
        end
        
        function Base.:(==)(a::MotifSig, b::MotifSig)
            return a.canong == b.canong && a.partition == b.partition
        end
        # Find frequencies of all unique connected subgraphs of size k in G
        function getmotifs(G,k; norm=true, verbose=false, colored=false)::Dict{MotifSig,Float64}
            # No speedup compared to []
            answers = Dict{MotifSig,Int64}()
            Visited = zeros(Bool,lg.nv(G))
            # For each node u
            for u in lg.vertices(G)
                if verbose
                    print("\r", round(u/lg.nv(G)*100), "% done")
                end
                S = Dict{Int64,Array{Int64,1}}()
                Visited .= false
                Visited[u] = true
                # S are parents?
                S[1] = [u]
                Enumerate_Vertex(G,u,S,k-1,1,Visited,answers; colored=colored)
            end
            norm ? normalise(answers) : answers
        end

        function normalise(answers)
            normalisation = sum(values(answers))
            d = Dict()
            for (key,v) in answers
                d[key] = v/normalisation
            end
            d
        end

        # Take graph, root vertex, "Selection": the vertices that are part of the current motif, the number of nodes left to choose
        function Enumerate_Vertex(G::GraphType,u::Int64,S,Remainder::Int64,i::Int64,Visited::Array{Bool,1},answers; colored=false)::Void where GraphType <: lg.AbstractGraph
            # If there are no more nodes to choose, terminate
            s = copy(S) # Stops shorter trees from accidentally sharing data. Must be a neater way of doing this.
            if Remainder == 0
                temp = vcat(values(s)...)
                # Most of memory and time usage is in the G[temp] call
                # Quicker if we bake in options
                #= k = Nauty.canonical_form(G[temp]).canong =#

                # Adding colouring - need to include partition with this somehow;
                # could probably just append it.
                
                # When we include support for complicated colours, we'll need to add a key of coloured nodes to each subgraph so that Nauty doesn't think that colours can be swapped
                canonfunc = colored ? ColoredGraphs.nauty : Nauty.baked_canonical_form
                k = canonfunc(G[temp])
                # Human readable alternative
                #k = Nauty.label_to_adj(Nauty.canonical_form(G[temp])[1],3)
                answers[MotifSig(k.canong,k.partition)] = get(answers,k,0) + 1
                return
            else
                # Find vertices that could be part of unique motifs
                ValList = Validate(G,s[i],u,Visited)
                # Make sure we don't try to pick more than we can
                n = min(length(ValList),Remainder)
                # Pick k nodes from our current depth
                for k = 1:n
                    for combination in it.subsets(ValList,k)
                        # set the current selection at the current depth to them nodes
                        # Might get a performance boost if s was just a list of numbers, and Parents was stored separately.
                        s[i+1] = combination
                        # and Remainder-k from other depths
                        Enumerate_Vertex(G,u,s,Remainder-k,i+1,Visited,answers; colored=colored)
                        # Repeat for all combinations of k nodes at current depth
                    end
                end
                # Tidy up
                for v in ValList
                    Visited[v] = false
                end
            end
        end

        # Take graph, selected vertices of previous layer, and the root vertex, return vertices that could form unique motifs
        # This is the bit where the labels are considered. Only labels bigger than root are considered. This is to stop double counting.
        function Validate(G::GraphType,Parents::Array{Int64,1},u::Int64,Visited::Array{Bool,1})::Array{Int64,1} where GraphType <: lg.AbstractGraph
            ValList = Array{Int64,1}()
            # For all of the immediate neighbours of the parents
            for v in Parents
                # Original paper has this as the neighbours of U. 100% sure that it's a typo.
                for w in lg.all_neighbors(G,v)
                    # If the label of the neighbour is greater than the parent, and the neighbour has not yet been visited,
                    # I think perhaps it should be the root node, actually
                    if (u < w) && !Visited[w]
                        # Mark it as visited and add it to our list of candidates
                        Visited[w] = true
                        push!(ValList, w)
                    end
                end
            end
            return ValList
        end

        # Questions: what's the difference between ValList and Visited?
    end

"""
    getmotifs(G,k;norm=true, verbose=false)::Dict{Array{UInt64,1},Float64}

Compute the frequencies of the unique connected subgraphs of size `k` in `G`.

The dict labels are adjacency matrices in the nauty format. Use `Nauty.label_to_adj` to convert a label to a LightGraphs compatible adjacency matrix.
"""
getmotifs = kavosh.getmotifs

# TODO:
#       - move normalisation out of subgraphs, support different methods
#       - support different backends
#       - support coloured graphs

end
