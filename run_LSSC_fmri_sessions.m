clc;
clear;
close all;

% Directory containing the files
dataDir = 'D:\UCSD_Acads\ProfGal_Research\data32\fMRIData\REST';

% Results diectory
resultDir = 'D:\UCSD_Acads\ProfGal_Research\test_run4';

% All sessions in the directory
sesList = dir(fullfile(dataDir, 'session-*'));

RUN_LSSC = 0;
RUN_DICE_SIMILARITY = 1;
cfg.thrcluster=[0.9]; % @renu - handle this well 

if (RUN_LSSC)
    % Run LSSC after concatenating all runs of a session for each subject
    for ses = 1:length(sesList) % processing each session
        sesname = sesList(ses).name;
        sesDir = [dataDir, '\', sesname];
        fileList = dir(fullfile(sesDir, '*.nii.gz'));
        
        V = [];
        prev_sub = '';
        prev_ses = '';
        for i = 1:length(fileList) % processing each .nii file
            fname = fileList(i).name;
            filePath = fullfile(sesDir, fname);
    
            sub_match = regexp(fname, 'sub-(\w{5})', 'tokens');
            ses_match = regexp(fname, 'ses-(\d+)', 'tokens');
            run_match = regexp(fname, 'run-(\d+)', 'tokens');
            sub_value = sub_match{1}{1};
            ses_value = str2double(ses_match{1}{1});
            run_value = str2double(run_match{1}{1});
            
            V0 = niftiread(filePath);
            V0 = squeeze(V0);
    
            if ((i == 1 || strcmp(prev_sub, sub_value)) && i ~= length(fileList))    % runs belong to same subject
                if isempty(V)
                    V = V0;
                else
                    V = cat(3, V, V0); % concatenate runs
                end 
            else                % New subject - perform LSSC to V matrix and initialize to new subject's data
                if (i == length(fileList)) % if last file
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
                
                % Configuring segmentation
                cfg.preProcess=false;
                cfg.N_TRIALS=1;
                cfg.n_clust = [100 ];
                cfg.makePlots = false;
                %cfg.thrcluster=[0.9:0.03:0.99];
                cfg.thrcluster=[0.9];
                cfg.NROWS = R;
                cfg.NCOLS = C;
                cfg.isoverlap = false;
                cfg.min_filled_area = 0.98;
                cfg.computeTempcorr = true; %@renu
                cfg.title_str = ['sub_' prev_sub 'session_' num2str(prev_ses) '_preproc_0_normal_1_pca_1_neig_51_nclust_100'];
                cfg.outputfilepath = resultDir;
               
                % Run segmentation          
                [labels_all, mergedA_all] = runROI_meso_nlm_new_v1(cfg, dFoF_masked, allregionspix, brain_mask);  
    
                % save outputs
                outfname = ['sub_' prev_sub 'ses_' num2str(prev_ses) '_lssc_out.mat'];
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

    % Group the files subject-wise
    pFileGrp = struct();
    for i = 1:length(pFileList)
        fname = pFileList(i).name;
        fpath = fullfile(processed_dir, fname);

        sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
        ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
        sub_value = sub_match{1}{1};
        ses_value = str2double(ses_match{1}{1});

        % Add the file to the subject group
        ses_key = sprintf('ses_%d', ses_value);
        if ~isfield(pFileGrp, sub_value)
            pFileGrp.(sub_value) = struct();
        end
        pFileGrp.(sub_value).(ses_key) = fpath;
    end

    % Compute Dice similarity between all the session-pairs for each subject
    subNames = fieldnames(pFileGrp);
    NUM_SES_PAIRS = 3;
    displayDice = zeros(length(subNames), NUM_SES_PAIRS);
    for i = 1:length(subNames)
        sub = subNames{i};

        dice_pairwise = cell(1, NUM_SES_PAIRS);
        pairs = cell(1, NUM_SES_PAIRS);
        for j = 1:NUM_SES_PAIRS % 3 session pairs
            r1 = j;
            r2 = rem(j, NUM_SES_PAIRS) + 1;
            k1 = sprintf('ses_%d', r1);
            k2 = sprintf('ses_%d', r2);

            LSSC_out_pair = cell(1, 2);
            LSSC_out_pair{1} = load(pFileGrp.(sub).(k1));
            LSSC_out_pair{2} = load(pFileGrp.(sub).(k2));
            
            % Compute dice coefficient for each cluster threshold value
            
            dice_values = zeros(1, length(cfg.thrcluster));
            for thr_id = 1:length(cfg.thrcluster)
                % Pairing similar cluster labels 
                pairing = pairComponents(LSSC_out_pair{1}.mergedA{thr_id}, LSSC_out_pair{2}.mergedA{thr_id});
                p_1 = pairing.p1;
                p_2 = pairing.p2;
            
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
            dice_pairwise{j} = dice_values;
            pairs{j} = sprintf('ses_%d-%d', r1, r2);
            displayDice(i, j) = dice_values;
        end
        pFileGrp.(sub).dice = dice_pairwise;
        pFileGrp.(sub).pairs = pairs;
    end
    
    data_dice = displayDice;
    data_dice(:, 4) = mean(displayDice, 2);
    subLabels = {'SLC01', 'SLC03', 'SLC04', 'SLC05', 'SLC06', 'SLC07', 'SLC08', 'SLC09', 'SLC10'};
    groupLabels = {'Session 1-2', 'Session 2-3', 'Session 3-1'};
    figure;
    %bar(data_dice);
    bar(displayDice);
    set(gca, 'XTickLabel', subLabels);
    title('Comparison between sessions');
    xlabel('Subjects');
    ylabel('Dice value');
    legend('Session 1-2', 'Session 2-3', 'Session 3-1');
    ylim([0 1]);

    figure;
    bar(mean(displayDice, 2));
    set(gca, 'XTickLabel', subLabels);
    title('Comparison between sessions');
    xlabel('Subjects');
    ylabel('Dice value');
    legend('Average');
    ylim([0 1]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dice_overall = multiclass_dice_coefficient(P1, P2) % renu
    % Get unique classes from both P1 and P2
    classes = unique([P1(:); P2(:)]);
    classes(classes == 0) = []; % Remove the background class if it's represented by 0

    % Initialize sums for intersections and total pixels
    intersection_sum = 0;
    total_p1_sum = 0;
    total_p2_sum = 0;

    % Loop over each class to calculate intersection and union
    for i = 1:length(classes)
        class = classes(i);
        p1_class = (P1 == class);
        p2_class = (P2 == class);

        intersection_sum = intersection_sum + sum(p1_class(:) & p2_class(:));
        total_p1_sum = total_p1_sum + sum(p1_class(:));
        total_p2_sum = total_p2_sum + sum(p2_class(:));
    end

    % Calculate the overall Dice coefficient
    dice_overall = (2 * intersection_sum) / (total_p1_sum + total_p2_sum);
end