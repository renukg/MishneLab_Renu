clc;
clear;
close all;

% Directory containing the files
dataDir = 'D:\UCSD_Acads\ProfGal_Research\data32\fMRIData\REST';

% Results diectory
resultDir = 'D:\UCSD_Acads\ProfGal_Research\test_hypsearch\test1_SLC10';

% All sessions in the directory
sesList = dir(fullfile(dataDir, 'session-*'));

RUN_LSSC = 1;
RUN_DICE_SIMILARITY = 1;
%SUBJECTS = ['SLC01', 'SLC03', 'SLC04', 'SLC05', 'SLC06', 'SLC07', 'SLC08', 'SLC09', 'SLC10'];
SUBJECTS = {'SLC10'};
cfg.thrcluster=[0.9]; % @renu - handle this well 

h_normalize = 0:1;
h_pca = 0:1;
h_kNN = 8:8:40;
h_self_tune = 4:8:36;

NumHypParms = length(h_normalize) * length(h_pca) * length(h_kNN) * length(h_self_tune);
% HYPERPARAMETER SEARCH @renu 

count_runs = 0;
StoreDice = struct();
for s = 1:length(SUBJECTS)
    StoreDice.(SUBJECTS{s}) = struct();
    StoreDice.(SUBJECTS{s}).normalize = [];
    StoreDice.(SUBJECTS{s}).pca = [];
    StoreDice.(SUBJECTS{s}).kNN = [];
    StoreDice.(SUBJECTS{s}).self_tune = [];
    StoreDice.(SUBJECTS{s}).Dice = {};
    StoreDice.(SUBJECTS{s}).avgDice = [];
end

for pca = h_pca
    for normalize = h_normalize
        for kNN = h_kNN
            for self_tune = h_self_tune
                if kNN < self_tune
                    continue;
                end

                count_runs = count_runs + 1;

                fprintf('\rProgress: %d/%d \n', count_runs, NumHypParms);

                % Configuring segmentation
                cfg.preProcess=false;
                cfg.N_TRIALS=1;
                cfg.n_clust = [100 ];
                cfg.makePlots = false;
                %cfg.thrcluster=[0.9:0.03:0.99];
                cfg.thrcluster=[0.9];
                cfg.isoverlap = false;
                cfg.min_filled_area = 0.98;
                cfg.computeTempcorr = false; %@renu - default false
                cfg.outputfilepath = resultDir;
                cfg.N_EIG = 51; % @renu  - 51 is deafult value
                cfg.normalize = normalize; %@renu - bool 1 is default value
                cfg.pca = pca;  %@renu - bool 1 is default value
                
                configParams = [];
                configParams.self_tune = self_tune; %32; % @renu make 16
                configParams.kNN = kNN; %50; % @renu make 32
                configParams.normalization = 'markov';
                configParams.maxInd = cfg.N_EIG;
                configParams.doNormalize = cfg.normalize;
                configParams.doPCA = cfg.pca;
                
                if (RUN_LSSC)
                    % Run LSSC after concatenating all runs of a session for each subject
                    for ses = 1:length(sesList) % processing each session
                        sesname = sesList(ses).name;
                        sesDir = [dataDir, '\', sesname];
                        fileList = dir(fullfile(sesDir, '*.nii.gz'));
                        
                        % Extract the names from fileList
                        fileNames = {fileList.name};  % Get all file names as a cell array
                        
                        % Initialize a logical array to keep track of matches
                        isMatch = false(size(fileNames));
                        
                        % Check if each filename contains any of the prescribed SUBJECTS
                        for i = 1:length(SUBJECTS)
                            isMatch = isMatch | contains(fileNames, SUBJECTS{i});
                        end
                        
                        % Filter the fileList based on the matching indices
                        filteredFileList = fileList(isMatch);

                        V = [];
                        prev_sub = '';
                        prev_ses = '';
                        for i = 1:length(filteredFileList) % processing each .nii file
                            fname = filteredFileList(i).name;
                            filePath = fullfile(sesDir, fname);
                    
                            sub_match = regexp(fname, 'sub-(\w{5})', 'tokens');
                            ses_match = regexp(fname, 'ses-(\d+)', 'tokens');
                            run_match = regexp(fname, 'run-(\d+)', 'tokens');
                            sub_value = sub_match{1}{1};
                            ses_value = str2double(ses_match{1}{1});
                            run_value = str2double(run_match{1}{1});
                            
                            % if i > 1 && ~ismember(prev_sub, SUBJECTS)  % edit 'SUBJECTS' if you want to run for selective subjects
                            %     prev_sub = sub_value;
                            %     prev_ses = ses_value;
                            %     if (i == 2) % Not the 1st subject as per the data sequence
                            %         V = [];
                            %     end
                            %     continue;
                            % end
                
                            V0 = niftiread(filePath);
                            V0 = squeeze(V0);
                    
                            if ((i == 1 || strcmp(prev_sub, sub_value)) && i ~= length(filteredFileList))    % runs belong to same subject
                                if isempty(V)
                                    V = V0;
                                else
                                    V = cat(3, V, V0); % concatenate runs
                                end 
                            else                % New subject - perform LSSC to V matrix and initialize to new subject's data
                                if (i == length(filteredFileList)) % if last file
                                    prev_sub = sub_value;
                                    prev_ses = ses_value;
                                end
                                % Run LSSC for the current data (V)
                                num_time_samples = size(V, 3);
                        
                                % Create brain mask
                                %maxV1 = (max(V(:,:,1,:),[],3)); % @Gal doubt
                                maxV = max(V, [], 3);
                                brain_mask = maxV ~= 0;
                                brain_mask(29:32,:) = 0; 
                                [R,C] = size(brain_mask);
                            
                                % Flatten the spatial dimensions
                                V_flat = reshape(V, [], size(V, 3));
                                
                                % Extract relevant pixels
                                allregionspix = find(brain_mask);
                                dFoF_masked = V_flat(allregionspix, :);
                                
                               
                                % Configuring segmentation @renu
                                cfg.NROWS = R;
                                cfg.NCOLS = C;
                                cfg.title_str = ['sub_' prev_sub 'session_' num2str(prev_ses) '_preproc_0_normal_1_pca_1_neig_51_nclust_100'];
                               
                                % Run segmentation   
                                %disp(configParams)
                                [labels_all, mergedA_all] = runROI_meso_nlm_new_hypsearch(cfg, configParams, ...
                                    dFoF_masked, allregionspix, brain_mask);  
                    
                                % save outputs
                                outfname = ['sub_' prev_sub, 'ses_' num2str(prev_ses),...
                                    '_nclust_',num2str(cfg.n_clust(1)),...
                                    '_thrcluster_',num2str(cfg.thrcluster(1)),...
                                    '_selftune_', num2str(configParams.self_tune),...
                                    '_kNN_', num2str(configParams.kNN), ...
                                    '_doNormalize_', num2str(configParams.doNormalize),...
                                    '_doPCA_', num2str(configParams.doPCA),...
                                    '_norm_', configParams.normalization, '_maxInd_', num2str(configParams.maxInd),'_lssc_out.mat'];
                                storepath = fullfile(resultDir, 'run_fmri_sessions');
                                fulloutpath = fullfile(storepath, outfname);
                                results = struct();
                                results.filename = outfname;
                                results.labels = labels_all;
                                results.mergedA = mergedA_all;
                                save(fulloutpath, '-struct', 'results');
                    
                                V = V0; % initialize to the new subject
                            end
                    
                            prev_sub = sub_value;
                            prev_ses = ses_value;
                        end
                    end
                end
                
                if (RUN_DICE_SIMILARITY)
                    processed_dir = fullfile(resultDir, 'run_fmri_sessions');
                    
                    % All processed files
                    pFileList = dir(fullfile(processed_dir, '*.mat'));

                    % Compute Dice similarity between all the session-pairs for each subject
                    NUM_SES_PAIRS = 3;
                    displayDice = zeros(length(SUBJECTS), NUM_SES_PAIRS);
                    for i = 1:length(SUBJECTS)
                        sub = SUBJECTS{i};
                
                        dice_pairwise = cell(1, NUM_SES_PAIRS);
                        pairs = cell(1, NUM_SES_PAIRS);
                        for j = 1:NUM_SES_PAIRS % 3 session pairs
                            r1 = j;
                            r2 = rem(j, NUM_SES_PAIRS) + 1;
                            k1 = sprintf('ses_%d', r1);
                            k2 = sprintf('ses_%d', r2);

                            resfname1 = ['sub_' sub, 'ses_' num2str(r1),...
                                    '_nclust_',num2str(cfg.n_clust(1)),...
                                    '_thrcluster_',num2str(cfg.thrcluster(1)),...
                                    '_selftune_', num2str(configParams.self_tune),...
                                    '_kNN_', num2str(configParams.kNN), ...
                                    '_doNormalize_', num2str(configParams.doNormalize),...
                                    '_doPCA_', num2str(configParams.doPCA),...
                                    '_norm_', configParams.normalization, '_maxInd_', num2str(configParams.maxInd),'_lssc_out.mat'];

                            resfname2 = ['sub_' sub, 'ses_' num2str(r2),...
                                    '_nclust_',num2str(cfg.n_clust(1)),...
                                    '_thrcluster_',num2str(cfg.thrcluster(1)),...
                                    '_selftune_', num2str(configParams.self_tune),...
                                    '_kNN_', num2str(configParams.kNN), ...
                                    '_doNormalize_', num2str(configParams.doNormalize),...
                                    '_doPCA_', num2str(configParams.doPCA),...
                                    '_norm_', configParams.normalization, '_maxInd_', num2str(configParams.maxInd),'_lssc_out.mat'];
                
                            LSSC_out_pair = cell(1, 2);
                            LSSC_out_pair{1} = load(fullfile(processed_dir, resfname1));
                            LSSC_out_pair{2} = load(fullfile(processed_dir, resfname2));
                            
                            check_unique_flag = 1;
                            % Compute dice coefficient for each cluster threshold value
                            dice_values = zeros(1, length(cfg.thrcluster));
                            for thr_id = 1:length(cfg.thrcluster)
                                % Pairing similar cluster labels 
                                pairing = pairComponents(LSSC_out_pair{1}.mergedA{thr_id}, LSSC_out_pair{2}.mergedA{thr_id});
                                p_1 = pairing.p1;
                                p_2 = pairing.p2;
                                up_1 = unique(p_1);
                                up_2 = unique(p_2);

                                if (length(p_1) ~= length(up_1)) || (length(p_2) ~= length(up_2))
                                    check_unique_flag = 0;
                                    break;
                                end
                            
                                labels_1 = LSSC_out_pair{1}.labels{thr_id};
                                labels_2 = LSSC_out_pair{2}.labels{thr_id};
                            
                                % Re-labeling labels for 2nd one
                                for k1=1:size(labels_2, 1)
                                    for k2=1:size(labels_2, 2)
                                        if (labels_2(k1, k2) ~= 0)
                                            tp = find(p_2 == labels_2(k1, k2));
                                            if (isempty(tp))
                                                labels_2(k1, k2) = 0;
                                            else
                                                labels_2(k1, k2) = p_1(tp);
                                            end
                                            
                                        end
                                    end
                                end
                            
                                dice_similarity = multiclass_dice_coefficient(labels_1, labels_2);
                                dice_values(thr_id) = dice_similarity;
                            end

                            if check_unique_flag == 0
                                continue;
                            end

                            dice_pairwise{j} = dice_values;
                            pairs{j} = sprintf('ses_%d-%d', r1, r2);
                            displayDice(i, j) = dice_values;
                        end
                        pFileGrp.(sub).dice = dice_pairwise;
                        pFileGrp.(sub).pairs = pairs;
                        
                        StoreDice.(sub).normalize(end+1) = normalize;
                        StoreDice.(sub).pca(end+1) = pca;
                        StoreDice.(sub).kNN(end+1) = kNN;
                        StoreDice.(sub).self_tune(end+1) = self_tune;
                        StoreDice.(sub).Dice{end+1} = displayDice(i, :);
                        StoreDice.(sub).avgDice(end+1) = mean(displayDice(i, :));

                    end
                    
                    data_dice = displayDice;
                    data_dice(:, 4) = mean(displayDice, 2);
                    subLabels = SUBJECTS;
                    groupLabels = {'Session 1-2', 'Session 2-3', 'Session 3-1'};
                    % figure;
                    % %bar(data_dice);
                    % bar(displayDice);
                    % set(gca, 'XTickLabel', subLabels);
                    % title('Comparison between sessions');
                    % xlabel('Subjects');
                    % ylabel('Dice value');
                    % legend('Session 1-2', 'Session 2-3', 'Session 3-1');
                    % ylim([0 1]);

                    % figure;
                    % bar(mean(displayDice, 2));
                    % set(gca, 'XTickLabel', subLabels);
                    % title('Comparison between sessions');
                    % xlabel('Subjects');
                    % ylabel('Dice value');
                    % legend('Average');
                    % ylim([0 1]);
                end
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%