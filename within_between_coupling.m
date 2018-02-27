clear all 
close all

%% Read in connectivity matrices
bold_net=dlmread('/data/joy/BBL/projects/zhouCbfNetworks/results/BOLDnetwork/105176/20170510x10571/net/SchaeferPNC_200/105176_20170510x10571_SchaeferPNC_200_network.txt');

cbf_net=dlmread('/data/joy/BBL/projects/zhouCbfNetworks/data/105176/20170510x10571/fcon/schaefer200/105176_20170510x10571_schaefer200_network.txt');
cbf_net = squareform(cbf_net);

%% Read in community index (denotes which module each region is assigned to)
Yeo_part=dlmread('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityAffiliation.1D');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Global-level coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set diagonal of adjacency matrix to nan
bold_net = bold_net - diag(diag(bold_net));
cbf_net = cbf_net - diag(diag(cbf_net));

sq_bold=squareform(bold_net)';
sq_cbf=squareform(cbf_net)';

[r p]=corr(sq_bold, sq_cbf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Module-level Coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define regional community index (ci)
ci=Yeo_part;

%% Define number of communities and subjects
numComm=length(unique(ci));
nsub=1;

%% Allocate empty matrices
corr_boldCbf_within_mat = zeros(nsub,numComm);
boldCbf_commComm_mat= zeros(nsub, numComm, numComm);

%% Loop through subjects
for s=1:nsub

	% Allocate empty community by community matrix to fill in
	comm_comm_mat=zeros(numComm,numComm);
	
	% Define connectivity matrices
	A_bold=bold_net;
	A_cbf=cbf_net;

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

figure; imagesc(squeeze(boldCbf_commComm_mat))