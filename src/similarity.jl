function split_and_keep(input_string::AbstractString, regex::Regex)
    tokens = String[]
    last_idx = firstindex(input_string)
    for match in eachmatch(regex, input_string)
        start_idx = match.offset
        # push non-matching segment if non-empty
        if start_idx > last_idx
            push!(tokens, String(input_string[last_idx:(start_idx - 1)]))
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

Tokenize a string into alphanumeric chunks and delimiter tokens.

The characters `+`, `/`, `-`, `_`, `(`, `)`, and `^` are treated as
**separators** and are kept as *stand‑alone* tokens. This way
input such as `"D_s"` and `"D(s)+"` or `"rho^+"` and `"rho(770)+"`
get broken into very similar token sequences, which improves fuzzy matching.
"""
function tokenize(s::String)::Vector{String}
    return split_and_keep(s, r"[+/\-_()^]")
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

    # Score keys but report suggestions in terms of "physical" particle names
    # (e.g. "D(s)+", "rho(770)+") rather than internal alias keys like "D_s_plus".
    scored = Dict{String, Float64}()
    for key in keys
        score = jaccard_similarity(user_tokens, tokenize(key))
        # Map key to a display name if possible
        display = try
            Particle(key).name
        catch
            key
        end
        # Keep the best score we have seen for a given display name
        if !haskey(scored, display) || score > scored[display]
            scored[display] = score
        end
    end

    similarities = collect(scored)  # Vector of Pair(display => score)
    sort!(similarities, by = x -> last(x), rev = true)
    n = min(length(similarities), 5)
    return first.(similarities[1:n])
end
