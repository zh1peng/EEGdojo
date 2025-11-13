function [S, meta] = build_similarity_single(B, varargin)
% Univariate-only behavioral similarity (or distance) for IS-RSA/MRQAP.
% Models: 'nn','absdiff','rank-absdiff','rbf','divergence','convergence'/'annak-min','annak-mean',
%         and for binary/nominal: 'same' (similarity 1/0) or Hamming distance.
%
% OPTIONS:
%   'return'   : 'similarity' (default) | 'distance'
%   'uniModel' : see above (default 'nn' for numeric; 'same' for nominal/binary)
%   'rank'     : true/false (default true for numeric/ordinal)
%   'rbfSigma' : [] (auto median heuristic) or >0

ip = inputParser;
ip.addParameter('return','similarity');
ip.addParameter('uniModel','nn');
ip.addParameter('rank',[]);
ip.addParameter('rbfSigma',[]);
ip.parse(varargin{:});
opt = ip.Results;

assert(isvector(B),'Provide a univariate vector.');
x = B(:);

% type dispatch
if islogical(x) || (isnumeric(x) && all(ismember(x(~isnan(x)),[0 1])))
    dtype = 'binary';
elseif iscategorical(x) || isstring(x) || iscellstr(x)
    if ~iscategorical(x), x = categorical(x); end
    dtype = 'nominal';
else
    dtype = 'numeric';
end

N = numel(x);
S = nan(N);

% helpers
mapD = @(D) (strcmpi(opt.return,'similarity') * (1 - (D / max(D(:)+eps)))) + ...
            (strcmpi(opt.return,'distance')   * D);
rbfSim = @(D,sig) exp(-(D.^2) ./ (2*(sig^2)));

switch dtype
    case 'binary'
        xb = double(x);
        D  = abs(xb - xb.');              % Hamming distance
        if any(strcmpi(opt.uniModel, {'same'}))
            S = 1 - D;                    % similarity
            if strcmpi(opt.return,'distance'), S = 1 - S; end
        else
            % default to Hamming distance/sim mapping
            S = mapD(D);
        end

    case 'nominal'
        S = double(x == x.');
        if strcmpi(opt.return,'distance'), S = 1 - S; end

    otherwise % numeric / ordinal
        if isempty(opt.rank), opt.rank = true; end
        if opt.rank
            r  = tiedrank(x);
            r01 = (N==1) * 1 + (N>1) * ((r-1)/(N-1));
        end

        m = lower(string(opt.uniModel));
        switch m
            case {'nn','nearest-neighbor','nn-absdiff'}
                r = tiedrank(x);                         % keep rank-based by definition
                S = 1 - abs(r - r.')/max(N-1,1);
                if strcmpi(opt.return,'distance'), S = 1 - S; end

            case {'absdiff','rank-absdiff'}
                if strcmpi(m,'rank-absdiff'), z = tiedrank(x); D = abs(z - z.');
                else,                              D = abs(x - x.'); end
                S = mapD(D);

            case {'rbf','gaussian'}
                D = abs(x - x.');
                if isempty(opt.rbfSigma)
                    tri = D(triu(true(N),1)); sig = median(tri(tri>0));
                    if isempty(sig) || ~isfinite(sig) || sig==0, sig = 1; end
                else
                    sig = opt.rbfSigma;
                end
                S = rbfSim(D, sig);
                if strcmpi(opt.return,'distance'), S = 1 - S; end

            case {'divergence'}
                if ~exist('r01','var')
                    r  = tiedrank(x); r01 = (N==1)*1 + (N>1)*((r-1)/(N-1));
                end
                S = 1 - (r01 + r01.')/2;
                if strcmpi(opt.return,'distance'), S = 1 - S; end

            case {'convergence','annak-min'}
                if ~exist('r01','var')
                    r  = tiedrank(x); r01 = (N==1)*1 + (N>1)*((r-1)/(N-1));
                end
                S = min(r01, r01.');
                if strcmpi(opt.return,'distance'), S = 1 - S; end

            case {'annak-mean'}
                if ~exist('r01','var')
                    r  = tiedrank(x); r01 = (N==1)*1 + (N>1)*((r-1)/(N-1));
                end
                S = (r01 + r01.')/2;
                if strcmpi(opt.return,'distance'), S = 1 - S; end

            otherwise
                error('Unknown uniModel "%s".', m);
        end
end

% symmetry + diagonal
S = 0.5*(S + S.');
if strcmpi(opt.return,'similarity'), S(1:N+1:end) = 1; else, S(1:N+1:end) = 0; end

meta = struct('N',N,'P',1,'mode','univariate','uniModel',char(opt.uniModel),...
              'return',lower(opt.return));
end
