function s_mle = SPADE_MLE(Q,q,sigma)
    % estimate the separation
    s_mle = 4 * sigma * sqrt(sum(q .* Q / sum(Q,'all')));
end
