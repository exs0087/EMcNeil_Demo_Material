% star_image_processing.m
%
% Generic pipeline for star-image centroiding & LOS vector computation,
% loading the image from a .mat file, without requiring the Image Processing Toolbox.

function star_image_processing(matFilename, imageVarName, threshold, focalLength, principalPoint)
    if nargin<1, matFilename    = 'star_image.mat';    end
    if nargin<2, imageVarName   = 'image';    end
    if nargin<3, threshold      = 30;            end
    if nargin<4, focalLength    = 60;            end
    if nargin<5, principalPoint = [6, 6];        end

    %% ─── PART 0: LOAD IMAGE ────────────────────────────────────────────────
    data = load(matFilename, imageVarName);
    I_raw = data.(imageVarName);
    I = double(I_raw);

    %% ─── PART 1: THRESHOLD & MANUAL CC LABELING ───────────────────────────
    bw = I >= threshold;
    [H,W] = size(bw);
    label = zeros(H,W);
    nextLabel = 0;
    neighbors = [-1 0; 1 0; 0 -1; 0 1];  % 4-conn

    for r = 1:H
        for c = 1:W
            if bw(r,c) && label(r,c)==0
                % new component
                nextLabel = nextLabel + 1;
                queue = [r c];
                label(r,c) = nextLabel;
                qhead = 1;
                while qhead <= size(queue,1)
                    rr = queue(qhead,1);
                    cc = queue(qhead,2);
                    qhead = qhead + 1;
                    for k = 1:4
                        nr = rr + neighbors(k,1);
                        nc = cc + neighbors(k,2);
                        if nr>=1 && nr<=H && nc>=1 && nc<=W ...
                           && bw(nr,nc) && label(nr,nc)==0
                            label(nr,nc) = nextLabel;
                            queue(end+1,:) = [nr nc];  %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
    numStars = nextLabel;

    %% ─── PART 2: CENTROIDS & OVERLAY ──────────────────────────────────────
    figure('Color','w');
    imagesc(flipud(I .* bw));
    colormap(parula);
    colorbar; axis image off; hold on;
    title('Thresholded Star Image with Centroids');

    centroids = zeros(numStars,2);
    for lbl = 1:numStars
        [rows, cols] = find(label==lbl);
        vals = I(sub2ind([H W], rows, cols));
        % convert pixel‐coords for plotting
        U = cols - 0.5;
        V = H - (rows - 0.5);
        C = sum(vals);
        cu = (U' * vals) / C;
        cv = (V' * vals) / C;
        centroids(lbl,:) = [cu, cv];
        plot(cu, cv, 'r*', 'MarkerSize', 12, 'LineWidth', 1.5);
        text(cu+1, cv+1, sprintf('%d',lbl), ...
             'Color','w','FontSize',12,'FontWeight','bold');
    end

    %% ─── PART 3: COMPUTE UNIT LOS VECTORS ────────────────────────────────
    u0 = principalPoint(1);
    v0 = principalPoint(2);
    e_meas = zeros(numStars,3);
    for k = 1:numStars
        u = centroids(k,1);
        v = centroids(k,2);
        vc = [ (u - u0);
               (v - v0);
                focalLength ];
        e_meas(k,:) = vc' / norm(vc);
    end

    %% ─── PART 4: STAR‐PAIR TABLE ─────────────────────────────────────────
    pairs    = nchoosek(1:numStars,2);
    numPairs = size(pairs,1);
    starPairs = table('Size',[numPairs,4], ...
                      'VariableTypes',{'string','string','double','double'}, ...
                      'VariableNames',{'Star_i','Star_j','Dot','Angle_deg'});

    for idx = 1:numPairs
        i = pairs(idx,1);
        j = pairs(idx,2);
        dp = dot(e_meas(i,:), e_meas(j,:));
        ang = acosd(min(max(dp,-1),1));
        starPairs(idx,:) = {sprintf('star%d',i), ...
                            sprintf('star%d',j), ...
                            dp, ang};
    end

    %% ─── SAVE RESULTS ────────────────────────────────────────────────────
    save('detected_centroids.mat','centroids','e_meas','starPairs');
    saveas(gcf,'centroid_overlay.png');

    fprintf('Found %d star clusters.\n', numStars);
    disp(starPairs);
end
