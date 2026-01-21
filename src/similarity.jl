function split_and_keep(input_string, regex)
    tokens = []
    last_idx = 1
    for match in eachmatch(regex, input_string)
        push!(tokens, input_string[last_idx:match.offset-1])
        push!(tokens, match.match)
        last_idx = match.offset + length(match) + 1
    end
    push!(tokens, input_string[last_idx:end])
    return tokens
end

function tokenize(s::String)::Vector{String}
    # tokenize as individual characters to better capture small typos
    return [string(c) for c in collect(s)]
end

function jaccard_similarity(A::Vector{String}, B::Vector{String})::Float64
    union_size = length(join(union(A, B)))
    intersection_size = length(join(intersect(A, B)))
    if union_size == 0
        return 1.0
    end
    return intersection_size / union_size
end

function closest_key_token_based(user_input::String, keys::Vector{String})
    user_tokens = tokenize(user_input)
    similarities = [(key, jaccard_similarity(user_tokens, tokenize(key))) for key in keys]
    sort!(similarities, by=x -> x[2], rev=true)
    return getindex.(similarities[1:5], 1)
end
