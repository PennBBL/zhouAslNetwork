%% Define numbers and subjects
nsub=26;
nreg=200;
nedge=19900;

% %% Read in CBF connectivity matrices
% cd('/data/joy/BBL/projects/zhouCbfNetworks/data/cbfProc/prelim_data_n30/')
% cbf_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/cbfProc/prelim_data_n30/*network.txt');
% nfiles = length(cbf_network_files);
% cbf_files = cell(1, nfiles);
% cbf_sq = zeros(nsub, nedge);
% 
% for k = 1:nfiles
%     cbf_net = dlmread(cbf_network_files(k).name);
%     cbf_net = squareform(cbf_net);
%     cbf_net = cbf_net - diag(diag(cbf_net));
%     cbf_sq(k,:) = squareform(cbf_net);
% end
% 
% %% Read in FA connectivity matrices
% cd('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/')
% fa_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/*FA_matrixts.csv');
% nfiles = length(fa_network_files);
% fa_files = cell(1, nfiles);
% fa_sq = zeros(nsub, nedge);
% 
% for k = 1:nfiles
%     fa_net = csvread(fa_network_files(k).name, 1, 0);
%     fa_net = fa_net - diag(diag(fa_net));
%     fa_sq(k,:) = squareform(fa_net);
% end

%% Read in community index (denotes which module each region is assigned to)
Yeo_part=dlmread('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityAffiliation.1D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Module-level Coupling %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define regional community index (ci)
ci=Yeo_part;
numComm=length(unique(ci));

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

            current_nodes_fa=sum(A_fa(comidx,:));
            current_nodes_fa=current_nodes_fa';
		
            current_nodes_cbf=sum(A_cbf(comidx,:));
            current_nodes_cbf=current_nodes_cbf';
            
            current_nodes_cbf=current_nodes_cbf(current_nodes_fa~=0);
            current_nodes_fa=current_nodes_fa(current_nodes_fa~=0);
			
			% Define a community X community matrix where elements represent within/between coupling
            comm_comm_mat(com1,com2)=corr(current_nodes_fa, current_nodes_cbf, 'type', 'Spearman');
			com2= com2 + 1;
		end
	
		%% Within module connectivity
        current_nodes_fa_within = squareform(A_fa(comidx,comidx));
		within_thresh_idx=find(current_edges_fa_within==0); % Define index for removing disconnect edges
		
%         current_edges_bold_within([within_thresh_idx])=[]; % Remove disconnected edges

        current_nodes_cbf_within = squareform(A_cbf(comidx,comidx));
		
%         current_edges_cbf_within([within_thresh_idx])=[]; % Remove disconnected edges

		% Correlation between nodes within/between modules
		corr_faCbf_within_mat (s,com1) = corr(current_nodes_fa_within', current_nodes_cbf_within', 'type', 'Spearman');
		com1 = com1 + 1;
	end
	faCbf_commComm_mat(s,:,:)=comm_comm_mat;
end

figure; imagesc(squeeze(faCbf_commComm_mat(2,:,:)));

%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('faCbf_commComm_nodeStrength.txt',faCbf_commComm_mat, ' ')
dlmwrite('faCbf_within_nodeStrength.txt',corr_faCbf_within_mat, ' ')
