clear, clc

disp('Step 2: IFMs')

mu(:, :, 2) = [0.75 0.85 0.50
    0.65 0.50 0.65
    0.50 0.85 0.75
    0.75 0.50 0.65];

nu(:, :, 2) = [0.15 0.10 0.40
    0.25 0.40 0.25
    0.40 0.10 0.15
    0.15 0.40 0.25];

mu(:, :, 1) = [0.35 0.50 0.25
    0.25 0.05 0.50
    0.75 0.65 0.35
    0.50 0.50 0.65];

nu(:, :, 1) = [0.55 0.40 0.65
    0.65 0.95 0.40
    0.15 0.25 0.55
    0.40 0.40 0.25];

mu(:, :, 4) = [0.85 0.75 0.75
    0.95 0.65 0.75
    0.75 0.85 0.65
    0.65 0.65 0.75];

nu(:, :, 4) = [0.10 0.15 0.15
    0.05 0.25 0.15
    0.15 0.10 0.25
    0.25 0.25 0.15];

mu(:, :, 3) = [0.25 0.35 0.50
    0.50 0.65 0.35
    0.15 0.50 0.75
    0.75 0.50 0.65];

nu(:, :, 3) = [0.65 0.55 0.40
    0.40 0.25 0.55
    0.80 0.40 0.15
    0.15 0.40 0.25];

[m, n] = size(mu(:, :, 1));
t = 4;

for k = 1:t
    temp = [reshape(mu(:, :, k)', 1, []); reshape(nu(:, :, k)', 1, [])];
    eval(['Y', num2str(k), '=  sprintf(''(%4.2f, %4.2f) (%4.2f, %4.2f) (%4.2f, %4.2f) \n'', temp)'])
end
clear Y1 Y2 Y3 Y4

mm = median(mu, 3);
mn = median(nu, 3);

disp('Step 3: Median matrix of F:')

mod = [];
mod = [mod, mm(:, 1), mn(:, 1), mm(:, 2), mn(:, 2), mm(:, 3), mn(:, 3)];

Median_m = sprintf('(%4.2f, %4.2f) (%4.2f, %4.2f) (%4.2f, %4.2f)  \n', mod')

%%% Step 3: Compute inversions of Fk

for k = 1:t
    Zm(:, :, k) = reshape(mu(:, :, k)', 1, 12);   % reshape to 1x12 (row-wise)
    Tm{k} = Zm(:, :, k);
end
Zm;

%% Step 5: Calculate the number of inversions of Y(k).

disp('Inversion counts for each element of mu(k):')
for k = 1:t
    sm{k} = 0;
    for j = 2:length(Tm{k})
        um{k} = sum((Tm{k}(1:j-1) > Tm{k}(j)));
        sm{k} = [sm{k}, um{k}];
    end
    inv_num_m(k) = sum(sm{k});
end
inv_num_m

for k = 1:t
    Zn(:, :, k) = reshape(nu(:, :, k)', 1, 12);   % reshape to 1x12 (row-wise)
    Tn{k} = Zn(:, :, k);
end
Zn;

disp('Inversion counts for each element of nu(k):')
for k = 1:t
    sn{k} = 0;
    for j = 2:length(Tn{k})
        un{k} = sum((Tn{k}(1:j-1) > Tn{k}(j)));
        sn{k} = [sn{k}, un{k}];
    end
    inv_num_n(k) = sum(sn{k});
end

inv_num_n

disp('Total inversion counts:')
inv_num = inv_num_m + inv_num_n

% Step 6: Calculate the number of inversions of F(median).

Y_m = reshape(mm', 1, 12);

smo = 0;
for j = 2:length(Y_m)
    um = sum((Y_m(1:j-1) > Y_m(j)));
    smo = [smo, um];
end
disp('Inversion counts for each element of mu in Ymode:')
smo_matrix = reshape(smo', 4, 3);
disp('The smo matrix is:');
smo_matrix
inv_num_mo = sum(smo)

Y_n = reshape(mn', 1, 12);

sno = 0;
for j = 2:length(Y_n)
    un = sum((Y_n(1:j-1) > Y_n(j)));
    sno = [sno, un];
end
disp('Inversion counts for each element of nu in Ymode:')
sno_matrix = reshape(sno', 4, 3);
disp('The sno matrix is:');
sno_matrix
inv_num_no = sum(sno)

inv_num_mode = inv_num_mo + inv_num_no

disp('Step 3: Entropy values of Yk:')

for k = 1:t
    Pm(:, :, k) = mu(:, :, k) ./ sum(sum(mu(:, :, k)));
    Em(k) = -sum(sum((Pm(:, :, k) .* log(Pm(:, :, k))))) ./ log(m * n);

    Pn(:, :, k) = nu(:, :, k) ./ sum(sum(nu(:, :, k)));
    En(k) = -sum(sum((Pn(:, :, k) .* log(Pn(:, :, k))))) ./ log(m * n);
end

Em
En

disp('Step 3: Total entropy of Yk:')
entropyk = Em + En

Pme = mm ./ sum(sum(mm));
Eme = -sum(sum((Pme .* log(Pme)))) ./ log(m * n)

Pne = mn ./ sum(sum(mn));
Ene = -sum(sum((Pne .* log(Pne)))) ./ log(m * n)

disp('Total entropy of E*:')
entropy_median = Eme + Ene

disp('Inversion-based CCs:')
IC = inv_num_mode ./ (inv_num_mode + abs(inv_num - inv_num_mode))

disp('Entropy-based CCs:')
EC = entropy_median ./ (entropy_median + abs(entropyk - entropy_median))

%% Step 7: Calculate the relative closeness of inversions.

disp('Step 6: Integrated CC:')
CC = (IC + EC) / 2

%% Step 8: Calculate the inversion-based weights of DMs

disp('Step 7: Weights of DMs:')
lambda = CC / sum(CC)

%%%% Step 9: Compute weighted decision matrix

disp('Step 9: Weighted to decision makers')

for k = 1:4
    tau(:, :, k) = 1 - (1 - mu(:, :, k)).^lambda(k);
    upsilon(:, :, k) = nu(:, :, k).^lambda(k);
    tmp = [reshape(tau(:, :, k)', 1, []); reshape(upsilon(:, :, k)', 1, [])];
    eval(['F', num2str(k), '=  sprintf(''(%4.2f, %4.2f)  (%4.2f, %4.2f)  (%4.2f, %4.2f)\n'', tmp)']);
end

%%%% Convert to alternative matrix

disp('Step 10: Group decision matrix')

xii = permute(tau, [3, 2, 1]);
oi = permute(upsilon, [3, 2, 1]);
for k = 1:4
    temp = [reshape(xii(:, :, k)', 1, []); reshape(oi(:, :, k)', 1, [])];
    eval(['H', num2str(k), '=  sprintf(''(%4.2f, %4.2f)  (%4.2f, %4.2f)  (%4.2f, %4.2f)  \n'', temp)']);
end

disp('Step 11: Weighted to attributes')

weights = [0.4, 0.4, 0.2];

for k = 1:4
    for j = 1:3
        xi(:, j, k) = 1 - (1 - xii(:, j, k)).^weights(j);
        o(:, j, k) = oi(:, j, k).^weights(j);
    end
    temp = [reshape(xi(:, :, k)', 1, []); reshape(o(:, :, k)', 1, [])];
    eval(['G', num2str(k), '=  sprintf(''(%4.2f, %4.2f) & (%4.2f, %4.2f) & (%4.2f, %4.2f)  \n'', temp)']);
end

%% Step 12: Establish reference matrices (PIS and NIS)
% xi and o dimensions: m ¡Á n ¡Á t (alternatives ¡Á attributes ¡Á experts)

% Positive ideal solution: max for xi, min for o (t ¡Á n)
xi_posi = max(xi, [], 3);
xi_nega = min(xi, [], 3);
o_posi = min(o, [], 3);
o_nega = max(o, [], 3);

%% Step 13: Compute GU+ and GU- for each alternative (matrix projection)

GU_posi = zeros(m, 1);
GU_nega = zeros(m, 1);

% Norm squares of G+ and G- (scalars)
G_plus_norm2 = sum(xi_posi(:).^2 + o_posi(:).^2 + (1 - xi_posi(:) - o_posi(:)).^2);
G_nega_norm2 = sum(xi_nega(:).^2 + o_nega(:).^2 + (1 - xi_nega(:) - o_nega(:)).^2);

for i = 1:m
    % Extract data for the i-th alternative (t ¡Á n)
    Gi_xi = squeeze(xi(i, :, :))';  % t ¡Á n
    Gi_o = squeeze(o(i, :, :))';    % t ¡Á n

    % Norm square of G_i
    Gi_norm2 = sum(Gi_xi(:).^2 + Gi_o(:).^2 + (1 - Gi_xi(:) - Gi_o(:)).^2);

    % Inner product of G_i and G+
    pi_Gi = 1 - Gi_xi(:) - Gi_o(:);
    pi_posi = 1 - xi_posi(:) - o_posi(:);
    inner_posi = sum(Gi_xi(:) .* xi_posi(:) + Gi_o(:) .* o_posi(:) + pi_Gi .* pi_posi);

    % Inner product of G_i and G-
    pi_nega = 1 - xi_nega(:) - o_nega(:);
    inner_nega = sum(Gi_xi(:) .* xi_nega(:) + Gi_o(:) .* o_nega(:) + pi_Gi .* pi_nega);

    % Normalized projection using formula (19)
    % NProj_Y(X) = |1 - XY/|X|^2 - XY/|Y|^2| / (|1 - XY/|X|^2 - XY/|Y|^2| + |1 - XY/|Y|^2|)

    % GU+: projection of G_i onto G+
    A = 1 - inner_posi / Gi_norm2 - inner_posi / G_plus_norm2;
    B = 1 - inner_posi / G_plus_norm2;
    GU_posi(i) = abs(A) / (abs(A) + abs(B));

    % GU-: projection of G_i onto G-
    C = 1 - inner_nega / Gi_norm2 - inner_nega / G_nega_norm2;
    D = 1 - inner_nega / G_nega_norm2;
    GU_nega(i) = abs(C) / (abs(C) + abs(D));
end

% Compute integrated correlation coefficient GC using formula (28)
GC = GU_posi ./ (GU_posi + GU_nega);

% Display results
disp('GU+ (normalized projection of each alternative onto PIS):');
disp(GU_posi);
disp('GU- (normalized projection of each alternative onto NIS):');
disp(GU_nega);
disp('GC (integrated correlation coefficient, i.e., URC):');
disp(GC);

% Ranking
[sorted_GC, idx] = sort(GC, 'descend');
disp('Alternative ranking (best to worst):');
for i = 1:m
    fprintf('Alternative %d: %.4f\n', idx(i), sorted_GC(i));
end