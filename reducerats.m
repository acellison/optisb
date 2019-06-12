function results = reducerats(cfs,nfs,nds,tol,maxdigits,numapprox,isbtarget,accurateisb)
    % Approximate the free weights with nearby rationals
    cfsets = cell(length(cfs),1);
    for i=1:length(cfsets)
        c = cfs(i);
        capprox = ratapprox(c,tol,maxdigits);
        dd = maxdigits;
        while isempty(capprox)
            dd = dd+1;
            capprox = ratapprox(c,tol,dd);
        end
        halflen = floor(numapprox/2);
        if length(capprox)>numapprox
            first = floor(length(capprox)/2)-halflen;
            capprox = capprox(first:first+numapprox-1);
        end
        cfsets{i} = capprox;
    end
    combs = makecombs(cfsets);

    nfs = sym(nfs);
    nds = sym(nds);
    Vd = rot90(vander(1./nds.^2));
    Vdi = inv(Vd);
    Vf = rot90(vandern(1./nfs.^2,length(nds)));
    b = sym([1;zeros(length(nds)-1,1)]);

    Nmax = double(max(max(nfs),max(nds)));
    ncores = corecount([double(nfs) double(nds)]);
    y = gbspolys(Nmax,'sym');

	% Loop over the combos until we find a match            
    results = {};
    for i=1:length(combs)
        cfs = combs(i,:)';
        cds = Vdi*(b-Vf*cfs);

        [~,d] = numden(cds);
        ndigits_dep = max(ceil(log10(double(d))));
        [~,d] = numden(cfs);
        ndigits_free = max(ceil(log10(double(d))));

        optpol = cfs'*y(nfs,:) + cds'*y(nds,:);
        optpol = fliplr(optpol);        

        % First approximate isb with fast algorithm
        isbn = isbn_est(fliplr(double(optpol)));
        if isbn >= isbtarget
            good = 1;
            if accurateisb
                % Now compute it with high precision
                isbn = double(isb_vpa(optpol)/(Nmax+1));
                good = isbn >= isbtarget;
            end
            if good
                config = struct('p',length(nds)*2,'N',Nmax, ...
                                'nds',nds,'nfs',nfs,'ncores',ncores, ...
                                'ndigits_dep',ndigits_dep,'ndigits_free',ndigits_free, ...
                                'cd',cds,'cf',cfs,'optpol',optpol,'isbn',isbn);
                results{end+1} = config;
            end
        end
    end
end

function combs = makecombs(vectors)
    n = numel(vectors); % number of vectors
    combs = cell(1,n); % pre-define to generate comma-separated list
    [combs{end:-1:1}] = ndgrid(vectors{end:-1:1}); % the reverse order in these two
    % comma-separated lists is needed to produce the rows of the result matrix 
    combs = cat(n+1, combs{:}); %concat the n n-dim arrays along dimension n+1
    combs = reshape(combs,[],n); %reshape to obtain desired matrix    
end

