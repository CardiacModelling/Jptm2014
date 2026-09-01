close all
clear all

models = {'tt06','grandi','ohara'};
model_names = {'ten Tusscher ''06','Grandi ''10','O''Hara ''11'};
channels = {'Kr','Ks','CaL','Na','to','K1'};

d = importdata('outputs_detailed_time.csv');
time = d.data(:,1);

num_models = length(models);
num_channels = length(channels);

counter = 1;
for c = 1:num_channels
    for m=1:num_models
        d = importdata([models{m} '_i' channels{c} '_block.csv']);
        subplot(num_channels,num_models,counter)                
        plot(time, d(:,2:end-1),'k-')
        hold on
        plot(time, d(:,end), 'k-','Linewidth', 2.0)
        
        xlim([-30 600])
        
        ylim([-90 60])
        set(gca,'FontSize',14)
        
        if c==1 % Top row
            title(model_names{m},'FontSize',16)
        end
        
        if c==num_channels
            xlabel('Time (ms)','FontSize',14)
        end
        
        if m==1 % First column
            ylabel('Voltage (mV)','FontSize',14)
        end
        counter = counter + 1;
    end
end
