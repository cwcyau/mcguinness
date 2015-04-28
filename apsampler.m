function options = gibbsampler(infile, outfile, options)

if nargin < 3
	options.burnin = 10000;
	options.thin = 10;
	options.maxIters = 100000;
	options.a = 1;
	options.b = 10;
end
if nargin == 3
	if isfield(options, 'u') & isfield(options, 'v')
		lam = (1-options.u)/options.u;
		options.a = (1/(1+lam))*( lam/(options.v*(1+lam)^2) - 1 );
		options.b = lam*options.a;
	else
		options.a = 1;
		options.b = 10;	
	end
end

rand('state', 1);
randn('state', 1);

col{1} = 'k';
col{2} = 'r';
col{3} = 'b';

% read in data
if exist(infile, 'file') 
	[sample, cellNo, y] = textread(infile, '%s %n %n', 'headerlines', 1, 'delimiter', '\t');
else
	disp(['Error! Data file: ' infile ' does not exist Please check that the correct directory/file name has been entered.']);
	return;
end

sampleNames = [];
sampleNames{1} = 'Control';
j = 2;
for i = 1 : length(sample)
	loc = strmatch( sample{i}, sampleNames, 'exact' );
	if isempty(loc)
		sampleNames{j} = sample{i};
		j = j + 1;
	end
end
nSamples = length(sampleNames);


cellId = [];
cellIdx = [];
groups = [];
Y = [];

% for each sample
groupId = 1;
for s = [ 1 : nSamples ]

    % find cells for this sample
    sampleLoc = strmatch(sampleNames(s), sample, 'exact');

    % get number of cells
    cellNumbers = unique( cellNo(sampleLoc) );
    nCells = length(cellNumbers);

    % for each cell
    groups = [ groups groupId*ones(1, nCells) ];
    cellId = [ cellId; cellNumbers([1:nCells]) ];
	cellIdx = [ cellIdx [1:nCells] ];    

    for c = 1 : nCells
        cellLoc = find( cellNo(sampleLoc) == cellNumbers(c) );
        Y = [ Y (y(sampleLoc(cellLoc))) ];
    end
    
    groupId = groupId + 1;

end

nGroups 		= length(unique(groups));
[T, N] 			= size(Y);

cellNames = unique(cellId);
nCells = length(cellNames);


% define hyperparameters
for i = 1 : nCells
	loc = find( cellNames == cellId(i) );
    Yi = Y(:, loc);
	Ymean(i) = mean(Yi(:));
    s_a(i) = 1;
    s_b(i) = 1/var(Yi(:));
end

lambda_u = mean(Y(:));
r_u = 1/(0.01*var(Y(:)));

lambda_m = 0;
r_m = 1/(var(Ymean(:)));
    
lambda_d = 0;
r_d = 1/var(Y(:));

p_a = options.a*ones(1, nGroups);
p_b = options.b*ones(1, nGroups);


% initialise parameters
u = lambda_u*ones(1, nGroups);
m = lambda_m*ones(1, nCells);
d = lambda_d*ones(1, nCells);
s = s_a.*s_b;
p = p_a./(p_a + p_b);
Z = ceil(rand(T, N));


% preallocate arrays
res_vec = zeros(options.maxIters, 	T*N, 			'single');
p_vec 	= zeros(options.maxIters, 	nGroups, 		'single');
u_vec 	= zeros(options.maxIters, 	nGroups, 		'single');
m_vec 	= zeros(options.maxIters, 	nCells, 		'single');
z_vec 	= zeros(T, 			N, options.maxIters, 	'uint8');
d_vec	= zeros(options.maxIters, 	nCells, 		'single');
s_vec 	= zeros(options.maxIters, 	nCells, 		'single');

p_a_vec = zeros(options.maxIters, 	nGroups,		'single');
p_b_vec = zeros(options.maxIters, 	nGroups,		'single');

lambda_d_vec = zeros(options.maxIters, 	1,		'single');
lambda_u_vec = zeros(options.maxIters, 	1,		'single');
lambda_m_vec = zeros(options.maxIters, 	1,		'single');

r_d_vec 	= zeros(options.maxIters, 	1,		'single');
r_u_vec 	= zeros(options.maxIters, 	1,		'single');
r_m_vec 	= zeros(options.maxIters, 	1,		'single');

for it = 1 : options.maxIters

    if mod(it, 1000) == 0
        disp(it);
    end

    % update allocations
    for i = 1 : N

        j = groups(i);
        q = find( cellNames == cellId(i) );

        pr(1, :) = lognormpdf( Y(:, i)', m(q) + u(j), 1/sqrt(s(q)) ) + log(1-p(j));
        pr(2, :) = lognormpdf( Y(:, i)', m(q) + d(q) + u(j), 1/sqrt(s(q)) ) + log(p(j));

        pr = exp( pr - repmat(logsumexp(pr, 1), [2 1]) );

        rndNo = rand(1, T);

        loc = find( rndNo >= pr(1, :) );
        
        Z(:, i) = 0;
        Z(loc, i) = 1;

    end

    % update mean levels for each treatment group
    for j = 1 : nGroups

        groupLoc 	= find( groups == j );
        nCell_j 	= length(groupLoc);
        
        S 	= 0;
        Sx 	= 0;

        for i = 1 : nCell_j
            
            cellNo = find( cellNames == cellId(i) );
            
            d_i = d(cellNo);
            m_i = m(cellNo);
            s_i = s(cellNo);
            
            Y_i = Y(:, groupLoc(i))';
            Z_i = Z(:, groupLoc(i))';

            S  = S  + s_i*T;
            Sx = Sx + sum( s_i.*( Y_i - m_i - d_i.*Z_i ) );
            
        end

        u_mean  = (Sx + r_u*lambda_u)/(S + r_u);
        u_std   = 1/sqrt(S + r_u);

        u(j) = u_mean + u_std*randn;

    end

	% update hyperparameters for u
	lambda_u_mu = ( (1/var(Ymean(:)))*mean(Ymean(:)) + r_u*sum(u) )/( 1/(var(Ymean(:))) + nGroups*r_u );
	lambda_u_std = 1/sqrt(var(Ymean(:)) + nGroups*r_u);
	
	a = nGroups + 1;
	b = var(Ymean(:)) + sum( (u - lambda_u).^2 );

    lambda_u = lambda_u_std*randn + lambda_u_mu;
	r_u = gamrnd( a, 1/b );
    
    % update mean levels and variance for each cell
    for i = 1 : nCells
        
		loc = find( cellId == cellNames(i) );
        
		S 	= 0;
		Sx 	= 0;
		for k = 1 : length(loc)
            
            j = groups(loc(k));

			Y_i = Y(:, loc(k));
			Z_i = Z(:, loc(k));
			s_i = s(i);
			d_i = d(i);
            u_j = u(j);
	
			S   = S + T*s_i;
			Sx  = Sx + sum( s_i*( Y_i - u_j - d_i.*Z_i ) );

		end

        m_mean 	= (Sx + r_m*lambda_m)/(S + r_m);
        m_std 	= 1/sqrt(S + r_m);        

		m(i) = m_mean + m_std*randn;

    end
    	
	% update hyperparameters for m
	lambda_m_mu = r_m*sum(m) /( 1/(var(Ymean(:))) + nCells*r_m );
	lambda_m_std = 1/sqrt(var(Ymean(:)) + nCells*r_m);
	
	a = nCells + 1;
	b = var(Ymean(:)) + sum( (m - lambda_m).^2 );

    lambda_m = lambda_m_std*randn + lambda_m_mu;
	r_m = gamrnd( a, 1/b );
    
    % update variance
    for i = 1 : nCells

        cellLoc = find( cellId == cellNames(i) );
        
        S = 0;
        a = 0;
        b = 0;
        for q = 1 : length(cellLoc)

            j = groups(cellLoc(q));
            
            m_i = m(i);
            d_i = d(i);
            u_j = u(j);
            Y_i = Y(:, cellLoc(q));
            Z_i = Z(:, cellLoc(q));
            
            S = sum( ( Y_i - u_j - m_i - d_i.*Z_i ).^2 );

            a = a + T;
            b = b + S;

        end

        a = a/2 + s_a(i);
        b = b/2 + s_b(i);

        s(i) = gamrnd( a, 1/b );
        
    end


    % update d
    for i = 1 : nCells

        cellLoc = find( cellId == cellNames(i) );
              
    	a = 0;
    	b = 0;
  
        for q = 1 : length(cellLoc)

            j = groups(cellLoc(q));
            
            s_i = s(q);
            m_i = m(i);
            u_j = u(j);
            Y_i = Y(:, cellLoc(q));
            Z_i = Z(:, cellLoc(q));
            
            activeLoc = find( Z_i == 1 );
            nActive = length(activeLoc);
            
            S   = nActive*s_i;
            Sx  = sum( s_i*( Y_i(activeLoc) - u_j - m_i ) );

            a = a + S;
            b = b + Sx;

        end
       
		% sample from truncated normal using auxilary variable method
		d_mean = (b + r_d*lambda_d)/(a + r_d);
		d_std = 1/sqrt(a + r_d);
			
		d_y = rand*exp(-0.5*((d(i)-d_mean)/d_std)^2);
		d_x = max(-d_mean/d_std, -sqrt(-2*log(d_y))) + rand*( sqrt(-2*log(d_y)) - max(-d_mean/d_std, -sqrt(-2*log(d_y))) );
		
		d(i) = d_std*d_x + d_mean;
 
    end	


	% update hyperparameters for d
	lambda_d_mu = r_d*sum(d)/( 1/var(Y(:)) + nCells*r_d );
	lambda_d_std = 1/sqrt(var(Y(:)) + nCells*r_d);

	a = nCells + 1;
	b = var(Y(:)) + sum( (d - lambda_d).^2 );
    
	lambda_d = lambda_d_std*randn + lambda_d_mu;
	r_d = gamrnd( a, 1/b );

    
    % update weights
    for j = 1 : nGroups

        loc = find( groups == j );
        Z_j = Z(:, loc);

        n_j = length( find( Z_j == 1 ) );
        n   = length( Z_j(:) );

        a = n_j + p_a(j);
        b = n - n_j + p_b(j);

        p(j) = betarnd(a, b);

    end
    
    % residuals
    for i = 1 : N

        q = find( cellNames == cellId(i) );

        j = groups(i);

        Y_i = Y(:, i);
        Z_i = Z(:, i);
        s_i = s(q);
        d_i = d(q);
        
        res(:, i) = sqrt(s_i)*( Y_i - m(q) - u(j) - d_i*Z_i );
        
    end
    
    res_vec(it, :) = res(:);

    p_vec(it, :) = p;
    u_vec(it, :) = u;
    m_vec(it, :) = m;
    z_vec(:, :, it) = Z;
    d_vec(it, :) = d;
	s_vec(it, :) = s;

	lambda_d_vec(it, :) = lambda_d;
	lambda_m_vec(it, :) = lambda_m;
	lambda_u_vec(it, :) = lambda_u;

	r_d_vec(it, :) = r_d;
	r_m_vec(it, :) = r_m;
	r_u_vec(it, :) = r_u;

	p_a_vec(it, :) = p_a;
	p_b_vec(it, :) = p_b;

end

range = (options.burnin+1) : options.thin : options.maxIters;
res = res_vec(range, :);



%hnd = figure(1);
%% set(hnd, 'Position', [1 1 1024 768]);

%mplot = 2;
%nplot = 3;

%subplot(mplot, nplot, 1);
%hold on;
%for j = 1 : nSamples
%	loc = find( groups == j );
%   	plot(loc, Y(:, loc), 'o', 'color', col{j});
%end
%axis square;
%xlim([ 0 size(Y, 2)+1]);
%ylim([ 0 1.1*max(Y(:)) ]);
%for i = 1 : nCells
%	cellStr{i} = num2str(cellId(i));
%end
%set(gca, 'Box', 'On', 'XTick', [1:1:size(Y, 2)], 'XTickLabel', cellStr);
%xlabel('Cell');
%ylabel('Intensity');
%%  title('(a) Data');

%subplot(mplot, nplot, 2);
%boxplot(Y);
%axis square;
%xlim([ 0 size(Y, 2)+1]);
%ylim([ 0 1.1*max(Y(:)) ]);
%set(gca, 'Box', 'On', 'XTick', [1:1:size(Y, 2)], 'XTickLabel', cellStr);
%xlabel('Cell')
%ylabel('Intensity');
%%  title('(b) Boxplot');

%subplot(mplot, nplot, 3);
%hold on;
%bins = linspace(0, 1, 100);
%for j = 1 : nSamples
%	nvals = hist( p_vec(range, j), bins );
%	plot(bins, nvals/sum(nvals), 'color', col{j}, 'LineStyle', '-');
%end
%legend(sampleNames);
%axis square;
%set(gca, 'Box', 'On');
%xlabel('Probability of large event');
%ylabel('Posterior Probability');
%%  title('(c) Event probabilities');

%subplot(mplot, nplot, 4);
%hold on;
%line([-5 5], [-5 5], 'Color', [1 1 1]*0.5);
%hnd = qqplot( randn(1, 10000), (res(:)-mean(res(:)))/std(res(:)) );
%set(hnd(1), 'Color', 'k', 'Marker', 'o');
%set(hnd(2), 'Visible', 'Off');
%set(hnd(3), 'Visible', 'Off');
%set(gca, 'Box', 'On');
%axis square;
%axis([-5 5 -5 5]);
%xlabel('Normal');
%ylabel('Residuals');
%%  title('(d) Residuals analysis');

%subplot(mplot, nplot, 5);
%hold on;
%u2 = u_vec(range, :);
%bins = linspace(min(u2(:)), max(u2(:)), 100);
%for j = 1 : nSamples
% 	loc = find( groups == j );
% 	nCellGroup = length(loc);
%	nvals = hist( u_vec(range, j), bins );
% 	plot(bins, nvals/sum(nvals), 'color', col{j}, 'LineStyle', '-');
%end
%legend(sampleNames);
%axis square;
%set(gca, 'Box', 'On');
%xlabel('Group amplitude');
%ylabel('Posterior Probability');
%%  title('(e) Amplitude analysis');

%subplot(mplot, nplot, 6);
%hold on;
%imagesc(mean(z_vec(:, :, range), 3));
%colorbar;
%axis square;
%xlim([ 1 size(Y, 2) ]);
%ylim([ 1 T ]);
%set(gca, 'Box', 'On', 'XTick', [1:1:size(Y, 2)], 'XTickLabel', cellStr);
%xlabel('Cell')
%ylabel('Data point');
%%  title('(f) Events');

d_vec = d_vec(range, :);
s_vec = s_vec(range, :);
m_vec = m_vec(range, :);
u_vec = u_vec(range, :);
res_vec = res_vec(range, :);
z_vec = z_vec(:, :, range);
p_vec = p_vec(range, :);
lambda_d_vec = lambda_d_vec(range, :);

save(outfile, 'lambda_d_vec', 'd_vec', 'm_vec', 'u_vec', 'res_vec', 'z_vec', 'p_vec', 'sampleNames', 'groups', 'cellId', 'nSamples', 'Y', 's_vec', 'options');

%print(outfileplot, '-depsc', '-r300');
