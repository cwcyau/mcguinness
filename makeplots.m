function makeplots(datfile, resultsfile, outfile)

rand('state', 1);
randn('state', 1);

col{1} = 'k';
col{2} = 'r';
col{3} = 'b';

load(resultsfile);

[sample, cellNo, y] = textread(datfile, '%s %n %n', 'headerlines', 1, 'delimiter', '\t');

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
groups = [];
Y = [];

groupId = 1;
for s = [ 1 : nSamples ]

    %         find cells for this sample
    sampleLoc = strmatch(sampleNames(s), sample, 'exact');

    %         get number of cells
    cellNumbers = unique( cellNo(sampleLoc) );
    nCells = length(cellNumbers);

    %         for each cell
    groups = [ groups groupId*ones(1, nCells) ];
    cellId = [ cellId [1:nCells] ];

    for c = 1 : nCells
        cellLoc = find( cellNo(sampleLoc) == cellNumbers(c) );
        Y = [ Y (y(sampleLoc(cellLoc))) ];
    end

    groupId = groupId + 1;

end

nGroups = length(unique(groups));
[T, N] = size(Y);

cellNames = unique(cellId);
nCells = length(cellNames);




range = 1 : size(u_vec, 1);

hnd = figure(1); clf;
set(hnd, 'Position', [1 1 1024 768]);

mplot = 2;
nplot = 3;
lineSz = 2;

figure(1);

subplot(2, 2, 1);
hold on;
for j = 1 : nSamples
    loc = find( groups == j );
    plot(loc, Y(:, loc), 'x', 'color', col{j}, 'MarkerSize', 6, 'LineWidth', 2);
end
axis square;
xlim([ 0 size(Y, 2)+1]);
ylim([ 0 1.1*max(Y(:)) ]);
cellStr = [];
for i = 1 : nCells
    cellStr{i} = num2str(cellId(i));
end
set(gca, 'Box', 'On', 'XTick', [1:1:size(Y, 2)], 'XTickLabel', cellStr);
xlabel('Cell');
ylabel('Amplitude');



subplot(2, 2, 3);
hold on;
bins = linspace(0, 1, 100);
for j = 1 : nSamples
    nvals = hist( p_vec(range, j), bins );
    plot(bins, nvals/sum(nvals), 'color', col{j}, 'LineStyle', '-', 'LineWidth', 2);
end
axis square;
set(gca, 'Box', 'On');
xlabel('Probability of large event');
ylabel('Posterior Probability');

subplot(2, 2, 2);
hold on;
u2 = u_vec(range, :);
bins = linspace(min(u2(:)), max(u2(:)), 50);
for j = 1 : nSamples

    loc = find( groups == j );
    nCellGroup = length(loc);

    ampdiff = repmat(u_vec(range, 1), [1 nCellGroup]) + m_vec(range, cellId(loc));

    nvals = hist( u_vec(range, j), bins );

    plot(bins, nvals/sum(nvals), 'color', col{j}, 'LineStyle', '-', 'LineWidth', 2);

end
ax = axis;
%ax = [ 0 120 ];

xlim([ax(1) ax(2)])
axis square;
set(gca, 'Box', 'On');
xlabel('Group amplitude');
ylabel('Posterior Probability');

%
% sample predictive distribution
%
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

[nIters, nGroups] = size(u_vec);
[T, N] = size(Y);

cellNo = 1;
ypred = zeros(nIters, nGroups);
for it = 1 : nIters
    
    sigma = 1./sqrt(s_vec(it, cellNo));
    d = d_vec(it, cellNo);
    
    for j = 1 : nGroups
    
        u = u_vec(it, j);
        
        z = randsample(2, 1, 'true', [1-p_vec(it, j) p_vec(it, j)] );
        
        zpred(it, j) = z;
        ypred(it, j) = u + sigma*randn + (z-1)*d;
    
    end
    
end
    
subplot(2, 2, 4);
hold on;

bins = linspace( min(ypred(:)), max(ypred(:)), 30);

for j = 1 : nGroups
    
    nvals = hist(ypred(:, j), bins);
%         plot(bins, nvals/sum(nvals), '-', 'color', col{j}, 'LineWidth', 3);
    
    loc = find( zpred(:, j) == 2 );
    nvals2 = hist(ypred(loc, j), bins);
    plot(bins, nvals2/sum(nvals), '--', 'color', col{j}, 'LineWidth', lineSz);
    
    loc = find( zpred(:, j) == 1 );
    nvals1 = hist(ypred(loc, j), bins);
    plot(bins, nvals1/sum(nvals), '-', 'color', col{j}, 'LineWidth', lineSz);
    
end
%xlim([0 200]);
axis square;
set(gca, 'Box', 'On');
xlabel('Amplitude');
ylabel('Predictive Probability');

print(outfile, '-r600', '-dpdf');

    
