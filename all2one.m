function [toOne] = all2one(X)
        X_min = min(X);
        X_max = max(X);
        toOne = (X - X_min) ./ (X_max - X_min);
end