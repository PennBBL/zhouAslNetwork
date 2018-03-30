%%%%%%%%%%%%%
%% Glasser %%
%%%%%%%%%%%%%

%% Define numbers and subjects
nsub=1373;
nreg=360;
nedge=64620;


%% Read in FA connectivity matrices
cd('/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/deterministic_dec2016/FA/GlasserPNC/')
fa_network_files = dir('*FA_GlasserPNC.mat');
nfiles = length(fa_network_files);
fa_sq = zeros(nsub, nedge);

for k = 1:nfiles
    fa_net = load(fa_network_files(k).name);
    fa_net = fa_net.connectivity;
    fa_net = fa_net - diag(diag(fa_net));
    fa_sq(k,:) = squareform(fa_net);
end

%% Apply community detection across subjects

q_mat = zeros(nsub,362);

for s=1:nsub

    A_fa=squareform(fa_sq(s,:));
    A = A_fa - diag(diag(A_fa));
    
    k = full(sum(A));
    twom = sum(k); 
    B = A - k'*k/twom;
    [S,Q] = genlouvain(B); 
    Q = Q/twom;
    
    subj_name = fa_network_files(s).name
    q_mat(s,1) = str2num(strtok(subj_name, '_'));
    q_mat(s,2) = Q;
    q_mat(s,3:362) = S';
end

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('modularity.txt',q_mat, ' ')

%%%%%%%%%%%%%%
%% Schaefer %%
%%%%%%%%%%%%%%

%% Define numbers and subjects
nsub=1689;
nreg=200;
nedge=19900;

%% Read in FA connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/faNetwork/')
fa_network_files = dir('*.mat');
nfiles = length(fa_network_files);
fa_sq = zeros(nsub, nedge);

for k = 1:10 %nfiles
    fa_net = load(fa_network_files(k).name);
    fa_net = fa_net.connectivity;
    fa_net = fa_net - diag(diag(fa_net));
    fa_sq(k,:) = squareform(fa_net);
end


%% Apply community detection across subjects

compath = '/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityAffiliation.1D';
outpath = '/data/joy/BBL/projects/zhouCbfNetworks/results/';

for s=1:10
    A_fa=squareform(fa_sq(s,:));
    A = A_fa - diag(diag(A_fa));
    quality(A, compath, outpath)
end

