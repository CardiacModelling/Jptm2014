% Code to examine the results of the variability study
close all
clear all

apd90_col_idx = 14;

models = {'ohara_rudy_2011','tentusscher_model_2006_epi','grandi_pasqualini_bers_2010'};
variants = {'AZ_orig_data','GSK_data','AZ_Gintant_hERG'};
% 
% % Multipliers of TQT conc to use.
concentration_cut_offs = [1 10 100 1000];

positive_effect_cut_off = 5.0; % ms

drug_name_file = ['../../input_data/AZ_HTS_data.txt'];

fileID = fopen(drug_name_file);
C = textscan(fileID, '%s %f %f %f %f %f %f %f %f %f %f','HeaderLines',1);
drug_names = C{1};
for i=1:length(drug_names)
   fprintf('%i: %s\n',i,drug_names{i}) 
end

[exp_names, exp_effect, exp_conc] = ...
    import_experimental_data(['clinical_tqt_results' filesep 'experimental_data.txt']);

% We're going to populate this structure with a list of positive /
% negative for contingency tables for each drug compound.
results = cell(length(models),length(variants),length(concentration_cut_offs),2); 

for v = 1:length(variants)
    variant = variants{v};
    
    for m = 1:length(models)
        model = models{m};    
        
        
        
        %% To load from this folder
        data_folder = ['..' filesep model filesep variant filesep];
        
        %% Loop over all the drugs, load data and plot.
        drugs_of_interest = 1:length(drug_names);
        %drugs_of_interest = 1;
        
        doses_for_cut_off = zeros(length(models),length(drugs_of_interest)) + 1000;
        doses_95_for_cut_off = zeros(length(models),length(drugs_of_interest)) + 1000;
 
        for drug_index = drugs_of_interest
            drug_name = drug_names(drug_index);
            
            % Find the row index for this drug in the experimental data file.
            for i=1:length(exp_names)
                if strcmp(exp_names{i},drug_name{1})
                    exp_data_row_idx = i;
                    break;
                end
            end
            
            figure(drug_index)
            subplot(1,3,v)
            if (v==2)
                th = title([drug_name{1} ' ' variant],'FontSize',16);
            else
                th = title(variant,'FontSize',16);
            end
            set(th,'interpreter','none')
            
            try
                filename = [data_folder 'voltage_results_' drug_name{1} '.dat'];
                importfile(filename);
            catch
                continue
            end
            
            parameter_captions = textdata(2:11);
            num_random_runs = max(data(:,1));
            temporary_experiment_params = zeros(num_random_runs,10);

            for run=(num_random_runs+1):-1:1
                this_run_rows = find(data(:,1)==run-1);
                temporary_experiment_params(run,:) = data(this_run_rows(1),2:11);
                
                doses = data(this_run_rows,12);
                apd90s = data(this_run_rows,apd90_col_idx);
                control_apd = apd90s(1); % Should be the same for each run, random or not.
                
                if run > 1
                    % Uncomment to overlay the individual traces
                    %semilogx(doses,apd90s-control_apd,'k-')
                    %hold on
                else
                    % This is the original 'no noise' trace.
                    if m==1
                        colour_to_use = 'b'; % O'Hara
                    elseif m==2
                        colour_to_use = 'r'; % Ten Tusscher 2006
                    else
                        colour_to_use = 'g'; % Grandi 2010
                    end
                    plot(log10(doses),apd90s-control_apd,[colour_to_use '-'],'Linewidth',2)
                    hold on
                end
            end

            set(gca,'FontSize',16)
            xlabh = xlabel('Conc (uM)','FontSize',18);
            ylabel('Change in APD90 (ms)','FontSize',18)
            ylim([-10 20])
            xlim(log10([min(doses) max(doses)+0.25*max(doses)]))
            set(gca,'XTick',[-3:1:2]);
            %format_ticks(gca,{'10^{-3}','10^{-2}','10^{-1}','10^0','10^1','10^2'});
            
            % Plot the magic cut off line
            plot(log10([doses(2) doses(end)]),[positive_effect_cut_off positive_effect_cut_off],'b--')
            
            % store the concentration at which we go above cut off (if ever).
            cut_off_idx = find(apd90s-control_apd>=positive_effect_cut_off,1,'first');
            if (~isempty(cut_off_idx))
                doses_for_cut_off_pos(m, drug_index) = doses(cut_off_idx);
            end
            
            % Look for shorteners instead of prolongers.
            cut_off_idx = find(apd90s-control_apd<=-positive_effect_cut_off,1,'first');
            if (~isempty(cut_off_idx))
                doses_for_cut_off_neg(m, drug_index) = doses(cut_off_idx);
            end
            
            if (num_random_runs > 0)
                
                %% Calculate two-tailed confidence interval
                concs = unique(data(:,12));
                confidence_wanted = 95;
                
                % We want
                lower_idx = floor(((100.0 - confidence_wanted)./2.0) * (num_random_runs/100.0));
                % This will throw an error if the number doesn't divide nicely - we may
                % then need to round up?
                upper_idx = num_random_runs - lower_idx;
                
                upper_confidence_interval = zeros(length(concs),1);
                lower_confidence_interval = zeros(length(concs),1);
                for conc_index = 1:length(concs)
                    this_conc_rows = find(data(:,12)==concs(conc_index));
                    
                    sorted_apd90s = sort(data(this_conc_rows, apd90_col_idx));
                    
                    upper_confidence_interval(conc_index) = sorted_apd90s(upper_idx)-control_apd;
                    lower_confidence_interval(conc_index) = sorted_apd90s(lower_idx)-control_apd;
                end
                
                patch_x = log10(doses(2:end))';
                patch_y = [lower_confidence_interval(2:end)'; upper_confidence_interval(2:end)'];
                px=[patch_x,fliplr(patch_x)];
                py=[patch_y(1,:), fliplr(patch_y(2,:))];
                p = patch(px,py,1,'FaceColor',colour_to_use,'EdgeColor','none');
                alpha(p, 0.2)
                               
                cut_off_idx = find(upper_confidence_interval>=positive_effect_cut_off,1,'first');
                
                if (~isempty(cut_off_idx))
                    doses_95_for_cut_off_pos(m, drug_index) = doses(cut_off_idx);
                end
                
                cut_off_idx = find(lower_confidence_interval<=-positive_effect_cut_off,1,'first');
                if (~isempty(cut_off_idx))
                    doses_95_for_cut_off_neg(m, drug_index) = doses(cut_off_idx);
                end
                
            end
            
            % Plot experimental data too
            plot(log10([doses(2) doses(end)]),[exp_effect(exp_data_row_idx)  exp_effect(exp_data_row_idx)],'k--')
           
            for c = 1:length(concentration_cut_offs)
                for prolong_or_short = 1:2
                    
                    if prolong_or_short==1 && c == 2
                        disp([variant ' ' model ' ' num2str(concentration_cut_offs(c))])
                    end
                                        
                    % If we have clinical TQT concentration data for this compound...
                    if (exp_conc(exp_data_row_idx) > 0)

                        % Mark the experimental data point on the graph.
                        plot(log10(exp_conc(exp_data_row_idx)), exp_effect(exp_data_row_idx), 'yx')
                        plot(log10(exp_conc(exp_data_row_idx)), exp_effect(exp_data_row_idx), 'ro')
                        plot(log10([exp_conc(exp_data_row_idx) exp_conc(exp_data_row_idx)]), ylim, 'k--')

                        %% Work out some measures of success                       
                        
                        % Generate a range of concentrations based on
                        % conc_cut_offs spread on log scale.
                        max_conc = concentration_cut_offs(c)*exp_conc(exp_data_row_idx);
                        if (max_conc > doses(end))
                            % Don't do wild extrapolation!
                            max_conc = doses(end);
                        end
                        min_conc = (1.0/concentration_cut_offs(c))*exp_conc(exp_data_row_idx);
                        concs_to_test = 10.0.^((linspace(log10(min_conc), log10(max_conc), 100)));
                                                
                        apd_prediction_at_exp_concs = interp1(doses,apd90s,concs_to_test);
                                                                              
                        % First - binary
                        sim_positive = false;
                        exp_positive = false;
                        if (prolong_or_short==1) 
                            % Prolongation is the 'positive effect'.
                            threshold_for_positive = control_apd + positive_effect_cut_off;
                            
                            % See how many of tested concs are positive.
                            temp_idxes = find(apd_prediction_at_exp_concs >= threshold_for_positive);
                            
                            if m==1 && v==1 && c==2
                               apd_prediction_at_exp_concs(temp_idxes)
                               threshold_for_positive
                            end
                            
                            if (exp_effect(exp_data_row_idx) >= positive_effect_cut_off)
                                exp_positive = true;
                                % If TQT is positive, then see if any of
                                % the tested concs are positive
                                if (~isempty(temp_idxes))
                                    sim_positive = true;
                                    if m==1 && v==1 && c==2
                                        disp('TQT +ve, SIM +ve')
                                    end
                                else
                                    if m==1 && v==1 && c==2
                                        disp('TQT +ve, SIM -ve')
                                    end 
                                end
                            else
                                % If TQT is negative we have agreement
                                % unless ALL of the tested concs are
                                % positive
                                if (length(temp_idxes)==length(apd_prediction_at_exp_concs))
                                    sim_positive = true;
                                    if m==1 && v==1 && c==2
                                        disp('TQT -ve, SIM +ve')
                                    end
                                else
                                    if m==1 && v==1 && c==2
                                        disp('TQT -ve, SIM -ve')
                                    end
                                end
                            end

                        else % Shortening is a 'positive' effect
                            threshold_for_positive = control_apd - positive_effect_cut_off;
                            % See how many of tested concs are positive.
                            temp_idxes = find(apd_prediction_at_exp_concs <= threshold_for_positive);
                            
                            if (exp_effect(exp_data_row_idx) <= -positive_effect_cut_off)
                                exp_positive = true;
                                % If TQT is positive, then see if any of
                                % the tested concs are positive
                                if (~isempty(temp_idxes))
                                    sim_positive = true;
                                end
                            else
                                % If TQT is negative we have agreement
                                % unless ALL of the tested concs are
                                % positive
                                if (length(temp_idxes)==length(apd_prediction_at_exp_concs))
                                    sim_positive = true;
                                end
                            end                            
                        end                        

                        % Draw up a 'confusion matrix' of drugs that showed an effect
                        % and those that didn't.
                        results{m,v,c,prolong_or_short} = [results{m,v,c,prolong_or_short}; sim_positive exp_positive];

                        %% Get the simulated percent change and plot against the real percent
                        % change at the estimated concentration.

                        if m==1 && v==1 && c==2
                            % Interpolate the simulated APD change                    
                            fprintf('%i: Drug: %s, TQT %5.3f ms,\t Sim %5.3f ms\n',...
                                exp_data_row_idx,drug_name{1},exp_effect(exp_data_row_idx),interp1(doses,apd90s-control_apd,exp_conc(exp_data_row_idx)))
                        end
                        %% Interpolate the simulated APD change                    
                        %fprintf('Drug: %s, TQT %5.3f ms,\t Sim %5.3f ms\n',...
                        %    drug_name{1},exp_effect(exp_data_row_idx),interp1(doses,apd90s-control_apd,exp_conc(exp_data_row_idx)))
                    end
                end            
            end % concentration cut offs
        end % Drug        
    end % models
end % variants

save('results.mat','results','models','variants','positive_effect_cut_off','concentration_cut_offs','-mat')

% Print summary to screen
for v = 1:length(variants)
    for m = 1:length(models)
        for c = 1:length(concentration_cut_offs)
            %for prolong_or_short = 1:2
            for prolong_or_short = 1

                % Work out totals for each box of a contingency table.
                true_pos = 0;
                false_pos = 0;
                false_neg = 0;
                true_neg = 0;
                % Data stored in this cell as sim_positive, TQT_positive
                for i=1:size(results{m,v,c,prolong_or_short},1)
                    if results{m,v,c,prolong_or_short}(i,1) && results{m,v,c,prolong_or_short}(i,2)
                        true_pos = true_pos +1;
                    end
                    if results{m,v,c,prolong_or_short}(i,1) && ~results{m,v,c,prolong_or_short}(i,2)
                        false_pos = false_pos +1;
                    end
                    if ~results{m,v,c,prolong_or_short}(i,1) && results{m,v,c,prolong_or_short}(i,2)
                        false_neg = false_neg +1;
                    end
                    if ~results{m,v,c,prolong_or_short}(i,1) && ~results{m,v,c,prolong_or_short}(i,2)
                        true_neg = true_neg +1;
                    end
                end

                % Sim_Expt names
                if (prolong_or_short==1)
                    fprintf('\nProlongation beyond %g ms:\n',positive_effect_cut_off)
                else
                    fprintf('\nShortening beyond %g ms:\n',positive_effect_cut_off)
                end
                fprintf('Model: %s\nData: %s\nSim correct within %i times TQT conc\n\tTQT +\t TQT -\n', models{m}, variants{v}, concentration_cut_offs(c))
                fprintf('Sim +\t%i\t%i\t%i\n',true_pos, false_pos, true_pos+false_pos)
                fprintf('Sim -\t%i\t%i\t%i\n',false_neg, true_neg, false_neg+true_neg)
                fprintf('\t%i\t%i\t%i\n',true_pos+false_neg,false_pos+true_neg,true_pos+false_pos+false_neg+true_neg)
                fprintf('Accuracy = %4.3g%%\n',100.0.*(true_pos+true_neg)./(true_pos+false_pos+false_neg+true_neg))
                fprintf('Sensitivity = %4.3g%%\n',100.0.*true_pos/(true_pos+false_neg))
                fprintf('Specificity = %4.3g%%\n',100.0.*true_neg/(false_pos+true_neg))
                fprintf('Positive predictive value = %4.3g%%\n',100.0.*true_pos/(true_pos+false_pos))
                fprintf('Negative predictive value = %4.3g%%\n',100.0.*true_neg/(false_neg+true_neg))
            end % Prolonger or shortener
        end % conc cut offs
    end % models
end % variants

