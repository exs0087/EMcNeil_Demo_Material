load('igrfcoefs.mat','coefs');
N = numel(coefs);
lens = arrayfun(@(c)numel(c.gh),coefs);
maxL = max(lens);
G = zeros(N, maxL);
for i=1:N
  ghi = coefs(i).gh;
  G(i,1:numel(ghi)) = ghi;
end
writematrix(G, 'gh.csv');
writematrix([coefs.year]', 'years.csv');
