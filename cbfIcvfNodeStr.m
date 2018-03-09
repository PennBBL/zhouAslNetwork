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
%% Read in ICVF connectivity matrices
% cd('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/')
% icvf_network_files = dir('/data/joy/BBL/projects/zhouCbfNetworks/data/noddiProc/prelim_data_n30/*ICVF_matrixts.csv');
% nfiles = length(icvf_network_files);
% icvf_files = cell(1, nfiles);
% icvf_sq = zeros(nsub,nedge);
% 
% for k = 1:nfiles
%     icvf_net = csvread(icvf_network_files(k).name, 1, 0);
%     icvf_net = icvf_net - diag(diag(icvf_net));
%     icvf_sq(k,:) = squareform(icvf_net);
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
%% CBF vs icvf %%
%%%%%%%%%%%%%%%

corr_icvfCbf_within_mat = zeros(nsub,numComm);
icvfCbf_commComm_mat= zeros(nsub, numComm, numComm);

%% Loop through subjects
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
	

            current_nodes_icvf=sum(A_icvf(comidx,:));
            current_nodes_icvf=current_nodes_icvf';
		
            current_nodes_cbf=sum(A_cbf(comidx,:));
            current_nodes_cbf=current_nodes_cbf';
            
%             current_nodes_cbf=current_nodes_cbf(current_nodes_icvf~=0);
%             current_nodes_icvf=current_nodes_icvf(current_nodes_icvf~=0);
			
			% Define a community X community matrix where elements represent within/between coupling
            comm_comm_mat(com1,1)=corr(current_nodes_icvf, current_nodes_cbf, 'type', 'Spearman');
			com2= com2 + 1;
        
		%% Within module connectivity
        %current_nodes_icvf_within = squareform(A_icvf(comidx,comidx));
%		within_thresh_idx=find(current_edges_icvf_within==0); % Define index for removing disconnect edges
		
%         current_edges_bold_within([within_thresh_idx])=[]; % Remove disconnected edges

        %current_nodes_cbf_within = squareform(A_cbf(comidx,comidx));
		
%         current_edges_cbf_within([within_thresh_idx])=[]; % Remove disconnected edges

		% Correlation between nodes within/between modules
		%corr_icvfCbf_within_mat (s,com1) = corr(current_nodes_icvf_within', current_nodes_cbf_within', 'type', 'Spearman');
		com1 = com1 + 1;
	end
	icvfCbf_commComm_mat(s,:,:)=comm_comm_mat;
end

figure; imagesc(squeeze(icvfCbf_commComm_mat(2,:,:)));



%% Write matrices in results directory
cd('/data/joy/BBL/projects/zhouCbfNetworks/results/')
dlmwrite('icvfCbf_commComm_nodeStrength.txt',icvfCbf_commComm_mat, ' ')
%dlmwrite('icvfCbf_within_nodeStrength.txt',corr_icvfCbf_within_mat, ' ')
