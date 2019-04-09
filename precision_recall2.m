%UNIVERSIDADE FEDERAL DO RIO DE JANEIRO
%LABORATÓRIO DE SINAIS, MULTIMÍDIA E TELECOMUNICAÇÕES
%AUTOR: TIAGO SIMÕES C. RODRIGUES
%SCRIPT PARA CÁLCULO DE PRECISION E RECALL DO ONSET_GD

clear all;
clc;
%Get data from result and ground_truth directories
cd result;
result_list = dir('*.mat');%use regex to exclude the '.' and '..' directories
cd ../ground_truth;
gt_list = dir('*.mat');
cd ..;

for index = 1:length(result_list)% for each instrument
    
    %Load data of result and ground truth to matlab workspace
    cd result;
    instrument_name = result_list(index).name;
    instrument_gt_name = gt_list(index).name;
    instrument = struct2cell(load(instrument_name));% to get the first value of the first instrument 'instrument{1}(1)', of the second: 'instrument{2}(1)'
    cd ../ground_truth;
    instrument_gt = struct2cell(load(instrument_gt_name));
    cd ..;
    true_positive = 0;
    false_positive = 0;
    false_negative = 0;
    for current_sample = 1:length(instrument_gt)%for each sample in the current instrument
        
        size_instr = length(instrument{current_sample});
        size_gt = length(instrument_gt{current_sample});
        %take the lowest value of sample size to avoid overlimit error
        %if length(instrument{current_sample})< length(instrument_gt{current_sample})
        %    low_size = length(instrument{current_sample});
        %else
        %    low_size = length(instrument_gt{current_sample});
        %end
        
        for data = 1:size_gt
           
           find = false;
           data_count = 1;% to reduce the number of loops
           
           while ~find && (data_count < size_instr) 
               
               if (abs(instrument_gt{current_sample}(data)-instrument{current_sample}(data_count))) <= 0.020
                   true_positive = true_positive + 1;
                   find = true;
               end
               data_count = data_count+1;       
           end
           
           if ~find
              false_negative = false_negative + 1; 
           end
           
        end
        false_positive = false_positive +abs(((size_gt - false_negative) - size_instr));
    end
    
    X = sprintf('Instrument: %s\nPrecision: %f\nRecall: %f\n',instrument_name(1:end-4),(true_positive)/(true_positive+false_positive),(true_positive)/(true_positive+false_negative));%In the article, the form how calculate Recall is different than this one
    disp(X);
end   
