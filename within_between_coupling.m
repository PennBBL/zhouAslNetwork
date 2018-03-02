clear all 
close all

%% Define numbers and subjects
nsub=26;
nreg=200;
nedge=19900;

%% Read in BOLD connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/boldNetwork/prelim_data_n30/')
bold_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/boldNetwork/prelim_data_n30/*network.txt');
nfiles = length(bold_network_files);
bold_files = cell(1, nfiles);
bold_sq = zeros(nsub, nedge);

for k = 1:nfiles
    bold_net = dlmread(bold_network_files(k).name);
    bold_net = bold_net - diag(diag(bold_net));
    bold_sq(k,:) = squareform(bold_net);
end

%% Read in CBF connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/cbfProc/prelim_data_n30/')
cbf_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/cbfProc/prelim_data_n30/*network.txt');
nfiles = length(cbf_network_files);
cbf_files = cell(1, nfiles);
cbf_sq = zeros(nsub, nedge);

for k = 1:nfiles
    cbf_net = dlmread(cbf_network_files(k).name);
    cbf_net = squareform(cbf_net);
    cbf_net = cbf_net - diag(diag(cbf_net));
    cbf_sq(k,:) = squareform(cbf_net);
end

%% Read in FA connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/')
fa_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/*FA_matrixts.csv');
nfiles = length(fa_network_files);
fa_files = cell(1, nfiles);
fa_sq = zeros(nsub, nedge);

for k = 1:nfiles
    fa_net = csvread(fa_network_files(k).name, 1, 0);
    fa_net = fa_net - diag(diag(fa_net));
    fa_sq(k,:) = squareform(fa_net);
end

%% Read in ODI connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/')
odi_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/*ODI_matrixts.csv');
nfiles = length(odi_network_files);
odi_files = cell(1, nfiles);
odi_sq = zeros(nsub, nedge);

for k = 1:nfiles
    odi_net = csvread(odi_network_files(k).name, 1, 0);
    odi_net = odi_net - diag(diag(odi_net));
    odi_sq(k,:) = squareform(odi_net);
end

%% Read in ICVF connectivity matrices
cd('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/')
icvf_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/*ICVF_matrixts.csv');
nfiles = length(icvf_network_files);
icvf_files = cell(1, nfiles);
icvf_sq = zeros(nsub,nedge);

for k = 1:nfiles
    icvf_net = csvread(icvf_network_files(k).name, 1, 0);
    icvf_net = icvf_net - diag(diag(icvf_net));
    icvf_sq(k,:) = squareform(icvf_net);
end

%% For running single subjects
%bold_net=dlmread('/data/joy/BBL/projects/zhouCbfNetworks/data/boldNetwork/105176/20170510x10571/net/SchaeferPNC_200/105176_20170510x10571_SchaeferPNC_200_network.txt');

%cbf_net=dlmread('/data/joy/BBL/projects/zhouCbfNetworks/data/cbfProc/105176/20170510x10571/fcon/schaefer200/105176_20170510x10571_schaefer200_network.txt');
%cbf_net = squareform(cbf_net);

%fa_net=csvread('/data/jux/BBL/projects/multishell_diffusion/processedData/Connectivity/93787_20160524x10163_FA_matrixsc.csv', 1, 0)

%% Read in community index (denotes which module each region is assigned to)
Yeo_part=dlmread('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityAffiliation.1D');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global-level coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial correlation
%for k = 1:length(bold_sq);
%     [r p] = corr(bold_sq(k,:), cbf_sq(k,:))
% end

%% For running single subjects
%bold_net = bold_net - diag(diag(bold_net));
%cbf_net = cbf_net - diag(diag(cbf_net));
%fa_net = fa_net - diag(diag(fa_net));

%sq_bold=squareform(bold_net)';
%sq_cbf=squareform(cbf_net)';
%sq_fa=squareform(fa_net)';

% [r p]=corr(sq_bold, sq_cbf)
% [r p]=corr(sq_fa, sq_cbf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Module-level Coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define regional community index (ci)
ci=Yeo_part;
numComm=length(unique(ci));

%%%%%%%%%%%%%%%%
%% BOLD vs CBF %
%%%%%%%%%%%%%%%%

%% Allocate empty matrices
corr_boldCbf_within_mat = zeros(nsub,numComm);
boldCbf_commComm_mat= zeros(nsub, numComm, numComm);

%% Loop through subjects
for s=1:nsub

	% Allocate empty community by community matrix to fill in
	comm_comm_mat=zeros(numComm,numComm);
	
	% Define connectivity matrices
	A_bold=squareform(bold_sq(s,:));
	A_cbf=squareform(cbf_sq(s,:));

	%% Define Modules and Nodes in network
	unique_S=unique(ci);
	numNodes=length(A_bold);

	%% Number of communities 
	numComm=length(unique_S);

	%% Set diagonal of adjacency matrix to nan
	A_bold = A_bold - diag(diag(A_bold));
	A_cbf = A_cbf - diag(diag(A_cbf));

	%% Iterate through each module to calculate correlation of edge weights (within and between)
	com1 = 1;

	for i=unique_S'
		com2 = 1;
		% Define index for nodes in each community
		comidx = find(ci==i);
		not_comidx = find(ci~=i);
	
		for j = unique_S'
			comidx_2= find(ci==j);
			% Pair-wise Between-module coupling
			current_edges_bold=A_bold(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_bold=triu(current_edges_bold,1)';
			end
		
			current_edges_cbf=A_cbf(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_cbf=triu(current_edges_cbf,1)';
            end
            
			current_edges_cbf=current_edges_cbf(current_edges_bold~=0);
			current_edges_bold=current_edges_bold(current_edges_bold~=0);
			
			% Define a community X community matrix where elements represent within/between coupling
			comm_comm_mat(com1,com2)=corr(current_edges_bold, current_edges_cbf, 'type', 'Spearman');
			com2= com2 + 1;
		end
	
		%% Within module connectivity
		current_edges_bold_within = squareform(A_bold(comidx,comidx));
		within_thresh_idx=find(current_edges_bold_within==0); % Define index for removing disconnect edges
		% current_edges_bold_within([within_thresh_idx])=[]; % Remove disconnected edges

		current_edges_cbf_within = squareform(A_cbf(comidx,comidx));
		% current_edges_cbf_within([within_thresh_idx])=[]; % Remove disconnected edges

		%% Correlation between edges within/between modules
		corr_boldCbf_within_mat (s,com1) = corr(current_edges_bold_within', current_edges_cbf_within', 'type', 'Spearman');
		com1 = com1 + 1;
	end
	boldCbf_commComm_mat(s,:,:)=comm_comm_mat;
end

figure; imagesc(squeeze(boldCbf_commComm_mat(1,:,:)))

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('boldCbf_commComm.txt',boldCbf_commComm_mat, ' ')
dlmwrite('boldCbf_within.txt',corr_boldCbf_within_mat, ' ')

%%%%%%%%%%%%%%%
%% CBF vs FA %%
%%%%%%%%%%%%%%%

corr_faCbf_within_mat = zeros(nsub,numComm);
faCbf_commComm_mat= zeros(nsub, numComm, numComm);

%% Loop through subjects
for s=2:nsub

	% Allocate empty community by community matrix to fill in
	comm_comm_mat=zeros(numComm,numComm);
	
	% Define connectivity matrices
	A_fa=squareform(fa_sq(s,:));
	A_cbf=squareform(cbf_sq(s,:));

	%% Define Modules and Nodes in network
	unique_S=unique(ci);
	numNodes=length(A_fa);

	%% Number of communities 
	numComm=length(unique_S);

	%% Set diagonal of adjacency matrix to nan
	A_fa = A_fa - diag(diag(A_fa));
	A_cbf = A_cbf - diag(diag(A_cbf));

	%% Iterate through each module to calculate correlation of edge weights (within and between)
	com1 = 1;

	for i=unique_S'
		com2 = 1;
		% Define index for nodes in each community
		comidx = find(ci==i);
		not_comidx = find(ci~=i);
	
		for j = unique_S'
			comidx_2= find(ci==j);
			% Pair-wise Between-module coupling
			current_edges_fa=A_fa(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_fa=triu(current_edges_fa,1)';
			end
		
			current_edges_cbf=A_cbf(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_cbf=triu(current_edges_cbf,1)';
            end
            
			current_edges_cbf=current_edges_cbf(current_edges_fa~=0);
			current_edges_fa=current_edges_fa(current_edges_fa~=0);
			
			% Define a community X community matrix where elements represent within/between coupling
			comm_comm_mat(com1,com2)=corr(current_edges_fa, current_edges_cbf, 'type', 'Spearman');
			com2= com2 + 1;
		end
	
		%% Within module connectivity
		current_edges_fa_within = squareform(A_fa(comidx,comidx));
		within_thresh_idx=find(current_edges_fa_within==0); % Define index for removing disconnect edges
		% current_edges_bold_within([within_thresh_idx])=[]; % Remove disconnected edges

		current_edges_cbf_within = squareform(A_cbf(comidx,comidx));
		% current_edges_cbf_within([within_thresh_idx])=[]; % Remove disconnected edges

		%% Correlation between edges within/between modules
		corr_faCbf_within_mat (s,com1) = corr(current_edges_fa_within', current_edges_cbf_within', 'type', 'Spearman');
		com1 = com1 + 1;
	end
	faCbf_commComm_mat(s,:,:)=comm_comm_mat;
end

figure; imagesc(squeeze(faCbf_commComm_mat(2,:,:)))

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('faCbf_commComm_ts.txt',faCbf_commComm_mat, ' ')
dlmwrite('faCbf_within_ts.txt',corr_faCbf_within_mat, ' ')

%%%%%%%%%%%%%%%%
%% CBF vs ODI %%
%%%%%%%%%%%%%%%%

corr_odiCbf_within_mat = zeros(nsub,numComm);
odiCbf_commComm_mat= zeros(nsub, numComm, numComm);

%% Loop through subjects

% First ODI sub has all 0s, omitted for now
for s=2:nsub

	% Allocate empty community by community matrix to fill in
	comm_comm_mat=zeros(numComm,numComm);
	
	% Define connectivity matrices
	A_odi=squareform(odi_sq(s,:));
	A_cbf=squareform(cbf_sq(s,:));

	%% Define Modules and Nodes in network
	unique_S=unique(ci);
	numNodes=length(A_odi);

	%% Number of communities 
	numComm=length(unique_S);

	%% Set diagonal of adjacency matrix to nan
	A_odi = A_odi - diag(diag(A_odi));
	A_cbf = A_cbf - diag(diag(A_cbf));

	%% Iterate through each module to calculate correlation of edge weights (within and between)
	com1 = 1;

	for i=unique_S'
		com2 = 1;
		% Define index for nodes in each community
		comidx = find(ci==i);
		not_comidx = find(ci~=i);
	
		for j = unique_S'
			comidx_2= find(ci==j);
			% Pair-wise Between-module coupling
			current_edges_odi=A_odi(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_odi=triu(current_edges_odi,1)';
			end
		
			current_edges_cbf=A_cbf(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_cbf=triu(current_edges_cbf,1)';
            end
            
			current_edges_cbf=current_edges_cbf(current_edges_odi~=0);
			current_edges_odi=current_edges_odi(current_edges_odi~=0);
			
			% Define a community X community matrix where elements represent within/between coupling
			comm_comm_mat(com1,com2)=corr(current_edges_odi, current_edges_cbf, 'type', 'Spearman');
			com2= com2 + 1;
		end
	
		%% Within module connectivity
		current_edges_odi_within = squareform(A_odi(comidx,comidx));
		within_thresh_idx=find(current_edges_odi_within==0); % Define index for removing disconnect edges
		% current_edges_bold_within([within_thresh_idx])=[]; % Remove disconnected edges

		current_edges_cbf_within = squareform(A_cbf(comidx,comidx));
		% current_edges_cbf_within([within_thresh_idx])=[]; % Remove disconnected edges

		%% Correlation between edges within/between modules
		corr_odiCbf_within_mat (s,com1) = corr(current_edges_odi_within', current_edges_cbf_within', 'type', 'Spearman');
		com1 = com1 + 1;
	end
	odiCbf_commComm_mat(s,:,:)=comm_comm_mat;
end

figure; imagesc(squeeze(odiCbf_commComm_mat(2,:,:)))

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('odiCbf_commComm_ts.txt',odiCbf_commComm_mat, ' ')
dlmwrite('odiCbf_within_ts.txt',corr_odiCbf_within_mat, ' ')

%%%%%%%%%%%%%%%%%
%% CBF vs ICVF %%
%%%%%%%%%%%%%%%%%

corr_icvfCbf_within_mat = zeros(nsub,numComm);
icvfCbf_commComm_mat= zeros(nsub, numComm, numComm);

%% Loop through subjects

% First ODI sub has all 0s, omitted for now
for s=2:nsub

	% Allocate empty community by community matrix to fill in
	comm_comm_mat=zeros(numComm,numComm);
	
	% Define connectivity matrices
	A_icvf=squareform(icvf_sq(s,:));
	A_cbf=squareform(cbf_sq(s,:));

	%% Define Modules and Nodes in network
	unique_S=unique(ci);
	numNodes=length(A_icvf);

	%% Number of communities 
	numComm=length(unique_S);

	%% Set diagonal of adjacency matrix to nan
	A_icvf = A_icvf - diag(diag(A_icvf));
	A_cbf = A_cbf - diag(diag(A_cbf));

	%% Iterate through each module to calculate correlation of edge weights (within and between)
	com1 = 1;

	for i=unique_S'
		com2 = 1;
		% Define index for nodes in each community
		comidx = find(ci==i);
		not_comidx = find(ci~=i);
	
		for j = unique_S'
			comidx_2= find(ci==j);
			% Pair-wise Between-module coupling
			current_edges_icvf=A_icvf(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_icvf=triu(current_edges_icvf,1)';
			end
		
			current_edges_cbf=A_cbf(comidx,comidx_2);
			if i==j
				% current_edges_bold([nan_idx])=0; % Set NaNs along diagonal to 0
				current_edges_cbf=triu(current_edges_cbf,1)';
            end
            
			current_edges_cbf=current_edges_cbf(current_edges_icvf~=0);
			current_edges_icvf=current_edges_icvf(current_edges_icvf~=0);
			
			% Define a community X community matrix where elements represent within/between coupling
			comm_comm_mat(com1,com2)=corr(current_edges_icvf, current_edges_cbf, 'type', 'Spearman');
			com2= com2 + 1;
		end
	
		%% Within module connectivity
		current_edges_icvf_within = squareform(A_icvf(comidx,comidx));
		within_thresh_idx=find(current_edges_icvf_within==0); % Define index for removing disconnect edges
		% current_edges_bold_within([within_thresh_idx])=[]; % Remove disconnected edges

		current_edges_cbf_within = squareform(A_cbf(comidx,comidx));
		% current_edges_cbf_within([within_thresh_idx])=[]; % Remove disconnected edges

		%% Correlation between edges within/between modules
		corr_icvfCbf_within_mat (s,com1) = corr(current_edges_icvf_within', current_edges_cbf_within', 'type', 'Spearman');
		com1 = com1 + 1;
	end
	icvfCbf_commComm_mat(s,:,:)=comm_comm_mat;
end

figure; imagesc(squeeze(icvfCbf_commComm_mat(2,:,:)))

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('icvfCbf_commComm_ts.txt',icvfCbf_commComm_mat, ' ')
dlmwrite('icvfCbf_within_ts.txt',corr_icvfCbf_within_mat, ' ')
