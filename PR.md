# Add similarity matching for particle name lookups

## Summary
Implements similarity-based suggestions when particle names are not found, improving user experience with typo-tolerant lookups.

## Changes
- **New file**: `src/similarity.jl` - Adds tokenization, Jaccard similarity, and closest key matching functions
- **Updated**: `src/literals.jl` - Uses similarity matching to suggest similar particle names on lookup failures
- **Code formatting**: Consistent spacing around operators and type parameters across codebase
- **Tests**: Added similarity fallback tests

## Implementation
Uses character-level tokenization and Jaccard similarity to find the 5 closest matching particle names when an exact match is not found, providing helpful suggestions in error messages.
