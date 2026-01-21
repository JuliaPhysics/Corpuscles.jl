function split_and_keep(input_string::AbstractString, regex::Regex)
    tokens = String[]
    last_idx = firstindex(input_string)
    for match in eachmatch(regex, input_string)
        start_idx = match.offset
        # push non-matching segment if non-empty
        if start_idx > last_idx
            push!(tokens, String(input_string[last_idx:start_idx-1]))
        end
        # push the matched delimiter token itself
        push!(tokens, String(match.match))
        last_idx = start_idx + ncodeunits(match.match)
    end
    if last_idx <= ncodeunits(input_string)
        push!(tokens, String(input_string[last_idx:end]))
    end
    return tokens
end

"""
    tokenize(s::String) -> Vector{String}

Tokenize a string into pieces separated by non-alphanumeric characters,
keeping delimiters like `+`, `/`, `-`, and `_` as part of the tokens.
This better respects particle-name structure while still being robust
to small typos.
"""
function tokenize(s::String)::Vector{String}
    return split_and_keep(s, r"[^A-Za-z0-9+/\-_]")
end

"""
    jaccard_similarity(A::Vector{String}, B::Vector{String}) -> Float64

Compute a Jaccard-like similarity between two token vectors, treating them
as sets of unique tokens. This favors overlaps in the set of characters
and ignores multiplicity and order.
"""
function jaccard_similarity(A::Vector{String}, B::Vector{String})::Float64
    union_size = length(join(union(A, B)))
    intersection_size = length(join(intersect(A, B)))
    if union_size == 0
        return 1.0
    end
    return intersection_size / union_size
end

"""
    closest_key_token_based(user_input::String, keys::Vector{String}) -> Vector{String}

Return up to the top 5 keys ordered by decreasing token-based Jaccard
similarity to the `user_input`. This is used only in error paths to
provide helpful "did you mean" suggestions.
"""
function closest_key_token_based(user_input::String, keys::Vector{String})
    isempty(keys) && return String[]
    user_tokens = tokenize(user_input)
    similarities = [(key, jaccard_similarity(user_tokens, tokenize(key))) for key in keys]
    sort!(similarities, by = x -> x[2], rev = true)
    n = min(length(similarities), 5)
    return getindex.(similarities[1:n], 1)
end
