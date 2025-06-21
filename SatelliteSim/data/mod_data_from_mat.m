% load the struct
s = load('igrfcoefs.mat','coefs');
coefs = s.coefs;                   
epochs = numel(coefs);

% pull out the years
years = [coefs.year]';
writematrix(years,'years.csv');

% find the max length of any gh vector
lens = arrayfun(@(c) numel(c.gh), coefs);
maxlen = max(lens);

% build a padded gh matrix
ghall = zeros(epochs, maxlen);
for i = 1:epochs
    v = double(coefs(i).gh(:));   % make sure it’s a column
    ghall(i,1:numel(v)) = v';     % fill as a row
end

% write it out
writematrix(ghall,'gh.csv');
