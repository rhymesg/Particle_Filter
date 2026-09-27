function weight = auxiliary_weights(posterior_likelihood, lookahead, ancestors)
% Correct APF second-stage weights for the selected ancestor proposal.
% See docs/particle-filtering.md for the two-stage importance ratio.
denominator = lookahead(ancestors);
assert(all(isfinite(denominator)) && all(denominator > 0), ...
    'APF:InvalidLookahead', 'Selected lookahead likelihoods must be positive.');
assert(all(isfinite(posterior_likelihood)) && all(posterior_likelihood >= 0), ...
    'APF:InvalidLikelihood', 'Likelihoods must be finite and nonnegative.');
log_weight = log(posterior_likelihood) - log(denominator);
assert(any(isfinite(log_weight)), 'APF:ZeroLikelihood', ...
    'All second-stage likelihoods underflowed or were zero.');
weight = exp(log_weight-max(log_weight));
weight = weight/sum(weight);
end
