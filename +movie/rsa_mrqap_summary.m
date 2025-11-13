function rsa_mrqap_summary(stats)
% RSA_MRQAP_SUMMARY  Print an R-like summary for rsa_mrqap output.

names = stats.names; if isempty(names), names = arrayfun(@(i)sprintf('b%d',i), 1:size(stats.beta,1), 'uni',0); end
T = size(stats.beta,2);
fprintf('MRQAP (%s), nPerm=%d, scheme=%s', stats.info.test, stats.info.nPerm, stats.scheme);
if isfield(stats.info,'strata') && ~isempty(stats.info.strata), fprintf(' [stratified]'); end
fprintf('\n');

if T>1
    fprintf('Time-resolved: T=%d, maxT=%d (FWER)\n', T, stats.info.maxT);
end

% Report OLS (use time-average for scalars)
R2   = mean(stats.ols.R2);
adjR2= mean(stats.ols.adjR2);
sigma= mean(stats.ols.sigma);
AIC  = mean(stats.ols.AIC);
BIC  = mean(stats.ols.BIC);
dof  = stats.ols.dof;
fprintf('Residual SE: %.4f on %d df | R^2: %.4f  Adj R^2: %.4f  AIC: %.3f  BIC: %.3f\n', ...
        sigma, dof, R2, adjR2, AIC, BIC);

% Table header
if T==1
    fprintf('\nCoefficients:\n');
    fprintf('%-14s %10s %10s %10s %10s\n', '', 'Estimate', 'Std.Err', 't', 'p(perm)');
    for j=1:numel(names)
        pj = NaN;
        if j==2 && isfield(stats,'p_main') && ~isempty(stats.p_main), pj = stats.p_main(1); end
        if isfield(stats,'p_all') && ~isempty(stats.p_all) && size(stats.p_all,1)>=j
            if ~isnan(stats.p_all(j,1)), pj = stats.p_all(j,1); end
        end
        fprintf('%-14s %10.6f %10.6f %10.3f %10s\n', names{j}, stats.beta(j,1), stats.se(j,1), stats.tval(j,1), pstr(pj));
    end
else
    % For T>1, print main predictor row with min/median/max across time
    j = 2;  % main
    est = stats.beta(j,:); se = stats.se(j,:); t = stats.tval(j,:);
    if isfield(stats,'p_main') && ~isempty(stats.p_main), p = stats.p_main(:)'; else, p = NaN(1,T); end
    fprintf('\n%s (time-resolved):\n', names{j});
    fprintf('  Estimate  min/med/max:  %.6f  %.6f  %.6f\n', min(est), median(est), max(est));
    fprintf('  Std.Err   min/med/max:  %.6f  %.6f  %.6f\n', min(se),  median(se),  max(se));
    fprintf('  t         min/med/max:  %.3f    %.3f    %.3f\n', min(t),   median(t),   max(t));
    if any(isfinite(p)), fprintf('  p(perm)  min/med/max:  %s  %s  %s\n', pstr(min(p)), pstr(median(p)), pstr(max(p))); end
end
fprintf('\n');

function s = pstr(x)
    if ~isfinite(x), s = 'NA'; return; end
    if x < 1e-4, s = sprintf('%.1e', x); else, s = sprintf('%.4f', x); end
end
end
