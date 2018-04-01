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
nsub=nfiles;
nreg=200;
nedge=19900;

%% Read in FA connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/faNetwork/')
fa_network_files = dir('*.mat');
nfiles = length(fa_network_files);
fa_sq = zeros(nsub, nedge);

for k = 1:nfiles
    fa_net = load(fa_network_files(k).name);
    fa_net = fa_net.connectivity;
    fa_net = fa_net - diag(diag(fa_net));
    fa_sq(k,:) = squareform(fa_net);
end


%% Apply community detection across subjects

Yeo_part=dlmread('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityAffiliation.1D');
q_mat = zeros(nsub,2);

S=Yeo_part;
gamma=1;

for i = 1:nsub
	A=squareform(fa_sq(i,:));
	N = size(A,1);
	twomu = 0;
	for s=1
    	k=sum(A(:,:,s));
    	twom=sum(k);
    	twomu=twomu+twom;
    	indx=[1:N]+(s-1)*N;
    	B(indx,indx)=A(:,:,s)-gamma*k'*k/twom;
	end
	q_mat(i,2) = sum(B(bsxfun(@eq,S,S.'))) ./ twomu;
    
    subj_name = fa_network_files(i).name;
    q_mat(i,1) = str2num(strtok(subj_name, '_'));
end

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('modularitySchaefer.txt',q_mat, ' ')
