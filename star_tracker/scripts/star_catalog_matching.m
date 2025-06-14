% star_catalog_matching.m
% - loads measured LOS vectors (e) from previous step
% - loads a user‐provided catalog of unit vectors
% - computes all pairwise dot‐products and absolute errors
% - writes sorted results to CSV

%% --- USER SETTINGS ---
dataDir = fullfile(pwd,'..','output');   % where `line_of_sight_vectors.mat` lives
outDir  = fullfile(pwd,'..','output');
if ~exist(outDir,'dir'), mkdir(outDir); end

% load measured vectors
S = load(fullfile(dataDir,'line_of_sight_vectors.mat'));
e_meas = S.e;           % 3×N_meas
N_meas = size(e_meas,2);

% define your catalog here: field = name, value = [3×1] unit vector
catalog.alpha = [ 0.0262; -0.0312; 0.0224 ];
catalog.Beta  = [ 0.0011;  0.0599;  -0.0759 ];
catalog.gamma = [ 0.9997;  0.9977;   0.9969 ];
% catalog.delta = ...      % add more as needed

catNames = fieldnames(catalog);
e_cat    = cellfun(@(n) catalog.(n), catNames,'Uni',0);
N_cat    = numel(e_cat);

%% --- compute all pairwise matches ---
measPairs = nchoosek(1:N_meas,2);
catPairs  = nchoosek(1:N_cat,2);

results = [];
for i = 1:size(measPairs,1)
    a = measPairs(i,1);  b = measPairs(i,2);
    dp_m = dot(e_meas(:,a), e_meas(:,b));
    for j = 1:size(catPairs,1)
        c = catPairs(j,1);  d = catPairs(j,2);
        dp_c = dot(e_cat{c}, e_cat{d});
        err  = abs(dp_m - dp_c);
        results(end+1,:) = [a,b,c,d,dp_m,dp_c,err]; 
    end
end

% sort by error
[~,I] = sort(results(:,7));
R = results(I,:);

% write to table
T = array2table(R, ...
    'VariableNames', {
        'meas_i','meas_j','cat_i','cat_j', ...
        'dot_meas','dot_cat','abs_error'});
writetable(T, fullfile(outDir,'star_match_results.csv'));

disp('Matching complete. See star_match_results.csv in output/');
