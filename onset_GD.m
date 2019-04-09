%Copyright (C) 2017 Speech and Music Technology Lab,
%Indian Institute of Technology Madras
                
%This file is part of GD based onset detection(onset_GD).                         
%onset_GD is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or   
%(at your option) any later version.
                
%This software is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.                             

%You should have received a copy of the GNU General Public License          
%.  If not, see <http://www.gnu.org/licenses/>

% AM Demoduation and Grp delay based onset detection 
% Initial Author :  P A Manoj Kumar
% Edited by : Jilt Sebastian
% Note: The placement of the file is essential. The binary WordSegmentWithSielnceRemoval and the script extrema.m must be in the
% same folder. This is the single resolution program, hence wsF is a parameter!



clear all;
close all;
clc;

%txt_file = 'result_1.txt';

warning('off','all');
warning;

% Parameter definitions
downsampling_rate = 10;
smoothening_factor = 44;    
winScaleFactor = 30;

thres = 0.005;

cd test_here; %place the test files here
listing = dir(pwd);
cd ..;

gt_name = {};
gt_strokes = [];
names = {};
strokes = [];
gt = dir('./test_here_gt/');
% Begin of paramterization


for indexA = 1:length(downsampling_rate)
    for indexB = 1:length(smoothening_factor)
            for indexD = 1:length(thres)
                fprintf('Paramterization details : dwnsmpl_rate = %d\n smth_factor = %d, thres = %f\n',downsampling_rate(indexA),smoothening_factor(indexB),thres(indexD));
                %fid3 = fopen(mat_file,'w');
                %fprintf(fid3,'#!MLF!#\n');
                
                for ii = 3:1:length(listing)
                    %======================================================================
                    % Part1: Using Amplitude Demodulation, and applying Group delay on it
                    filename = listing(ii).name;
                    orig_filename = filename;
                    names{ii-2} = orig_filename;
                    cd test_here;
                    [Y,Fs] = audioread(orig_filename);%change wavread to audioread
                    % Plot the signal in time domain
                    %dt = 1/Fs;
                    %t = 0:dt:(length(Y)*dt)-dt;
                    %plot(t,Y,'r');
                    %legend('Agogo signal');
                    %xlabel('seconds');
                    %ylabel('Amplitude');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    cd ..;    
                    %already_write = false;
                    DF=diff(Y);   % Differentiate it to emphasize the frequency components in amplitude
                    
                    %Plot diff signal in time domain
                    %figure;
                    %t1 = 0:dt:(length(DF)*dt)-dt;
                    %plot(t1,DF,'r');
                    %legend('Derivative of Music signal');
                    %xlabel('seconds');
                    %ylabel('Amplitude');
                    fprintf('Working on %s\n',filename);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    a = hilbert(DF);
                    Z = abs(a+DF);
                    D = downsample(Z,downsampling_rate(indexA));%reduce calculations
                    S = smooth(D,smoothening_factor(indexB),'moving');%yy = smooth(y) smooths the data in the column vector y using a moving average filter. Results are returned in the column vector yy. The default span for the moving average is 5.
                        %yy = smooth(y,span,method) sets the span of method to span. For the loess and lowess methods, span is a percentage of the total number of data points, less than or equal to 1 
                    % Plot the envelope estimated using Hilbert transform
                    %figure;
                    %t2 = 0:dt:(length(S)*dt)-dt;
                    %plot(t2,S,'r');
                    %legend('Envelope of Music signal');
                    %xlabel('seconds');
                    %ylabel('Amplitude');
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   
                    %X = ifft(S.^2);% It's not necessary i think
                    
                    fprintf('Working on %s\n',filename);
                    assignin('base','S',S);% assignin(ws, 'var', val) assigns the value val to the variable var in the workspace ws
                    assignin('base','Y',Y);%ws can have a value of 'base' or 'caller' to denote the MATLAB base workspace or the workspace of the caller function.
                    grp_delay = ones(length(S),1);
                    gd_sum = ones(length(S),1);
                    
                    for wsfIndex = 1:length(winScaleFactor)%pode tirar o par, vale a pena?
                        
                        tempDir = sprintf('temp_%d',wsfIndex);
                        mkdir(tempDir); 
                        cd(tempDir);
                        energy_file_name = strcat(filename(1:end-4),'.en');
                        dlmwrite(energy_file_name,S*1000,'\n');%Write matrix S*1000 to 'energy...' file, delimited by the '\n' character
                        spec_file_name = strcat(energy_file_name(1:end-2),'spec');
                        copyfile(energy_file_name,spec_file_name);% assim foi criado o arquivo spec
                        % Invoking the binary - ele so copiou o config file para um temp, cade o binario?                       
                        copyfile('../fe-words.base_ref','fe-words.base');
                        ctrl_file = 'fe-words.base';
                        temp_ctrl_file = strcat('temp.base');% NAO HOUVE CONCATENACAO, PQ DO USO?
                        % Changing the winscalefactor parameter in config file
                        a = importdata(ctrl_file);
                        a = struct2cell(a);
                        a{1}(3) = winScaleFactor(wsfIndex);
                        
                        fprintf('Window scale factor is %d\n',winScaleFactor(wsfIndex));
                        fid0 = fopen(temp_ctrl_file,'w');
                        for i = 1:length(a{1})
                            fprintf(fid0,'%s %s %f\n',char(a{2}(i,1)),char(a{2}(i,2)),a{1}(i));
                        end
                        copyfile(temp_ctrl_file,ctrl_file);
                        delete(temp_ctrl_file);
                        fclose(fid0);

                        dummy1 = 'b';
                        dummy2 = 'c';
                        dummy3 = 'd';
                        dummy4 = 'e';
                        dump = 'dump.txt';

                        system(sprintf('../WordSegmentWithSilenceRemoval %s %s %s %s %s %s %s > %s 2>&1',ctrl_file,energy_file_name,spec_file_name,dummy1,dummy2,dummy3,dummy4,dump));

                        delete(energy_file_name);
                        temp = load(spec_file_name);
                        delete(spec_file_name);
                        temp = temp(:,1);%change 5 to 1
                        temp(length(S)+1:end) = [];%essa linha
                        grp_delay = grp_delay.*temp;
                        temp = temp - mean(temp);
                        gd_sum = gd_sum + cumsum(temp);
                        cd ..; 
                    end
                    
                    
                    grp_delay = -diff(gd_sum);%missing '-'
                    %figure;
                    %t3 = 0:dt:(length(grp_delay)*dt)-dt;
                    %plot(t3,-grp_delay,'r');
                    %legend('Minimum phase group delay computed on the envelope');
                    %xlabel('seconds');
                    %ylabel('Amplitude');
                    
                    grp_delay = smooth(grp_delay,2*smoothening_factor(indexB),'moving');   % A moving average with 1 ms interval
                    
                    grp_delay = grp_delay/max(grp_delay);
                    %figure;
                    %t4 = 0:dt:(length(grp_delay)*dt)-dt;
                    %plot(t4,-grp_delay,'r');
                    %legend('Minimum phase group delay computed on the envelope smoothed');
                    %xlabel('seconds');
                    %ylabel('Amplitude');
                    %======================================================================



                    %======================================================================
                    % Part2: Reading the contents of group delay file, and getting the onsets
                    threshold = thres(indexD);

                    stroke_loc = zeros(1,length(grp_delay));
                    % Go to each minima, and calculate height till next maxima. Keep a threshold on this to decide if stroke!
                    t = 1:length(grp_delay);
                    [ymax,imax,ymin,imin] = extrema(grp_delay);

                    % sort the minimas and maximas;
                    temp_min = sortrows([imin ymin]);
                    imin = temp_min(:,1)';
                    ymin = temp_min(:,2)';
                    clear temp_min;

                    temp_max = sortrows([imax ymax]);
                    imax = temp_max(:,1)';
                    ymax = temp_max(:,2)';
                    clear temp_max;

                    if (imin(1) < imax(1) )  % fine, just truncate the maximum
                        imin(1) = []; ymin(1) = [];

                        if (length(imin) > length(imax) )
                            imin(length(imax)+1:end) = [];
                            ymin(length(imax)+1:end) = [];
                        elseif (length(imin) < length(imax) )
                            imax(length(imin)+1:end) = [];
                            ymax(length(imin)+1:end) = [];
                        end
                    else                                                    

%                             imax(1) = []; ymax(1) = [];
                        if (length(imin) > length(imax) )
                            disp('this shouldnt have come');
                            imin(length(imax)+1:end) = [];
                            ymin(length(imax)+1:end) = [];
                        elseif (length(imin) < length(imax) )

                            imax(length(imin)+1:end) = [];
                            ymax(length(imin)+1:end) = [];
                        end
                    end

                    assignin('base','ymax',ymax);
                    assignin('base','imax',imax);
                    assignin('base','ymin',ymin);
                    assignin('base','imin',imin);
                    assignin('base','grp_delay',grp_delay);

                    %==================================================================
                     % Algorithm1  for stroke location
                     index_stroke = 1;
                     peak_valley_heights = ymax - ymin;

                     for index = 1:1:length(peak_valley_heights)
                         if (peak_valley_heights(index) > threshold)%para resolver o problema do tanta e tamborim so diminuir esse limite?
                     %         fprintf('Interest point at %d\n',ceil((imin(index) + imax(index))/2)); 
                             stroke_loc(index_stroke) = ceil((imin(index) + imax(index))/2);
                             index_stroke = index_stroke + 1;
                         end
                     end
                     %==================================================================
                    stroke_loc(stroke_loc==0) = [];% substitute the zeros for an empty list

                    assignin('base','stroke_loc',stroke_loc);
                    assignin('base','peaks',peaks);

                    %======================================================================



                    %======================================================================
                    % Printing in standard MLF format
                    dangerflag = 0;
                    cd test_here;
                    [X,Fs] = audioread(filename);
                    cd ..;
                    fid3 = fopen('result_1.txt','a');
                    length_wav_file = length(X)*1/Fs;
                    stroke_loc = stroke_loc*downsampling_rate(indexA)/Fs;  % Converting into seconds
                    if (isempty(stroke_loc))        % Provision for null strokes
                        filename = filename(1:end-4); 
                        fprintf(fid3,'\"*/%s.lab\"\n0\t%f\n%f\t%f',filename,length_wav_file-0.0100000,length_wav_file-0.0100000,length_wav_file);
                    else                   
                        if (length_wav_file - 0.01 < stroke_loc(end) )  % Artificial stroke at end due to compuation of group delay function
                            stroke_loc(end) = [];
                            if (isempty(stroke_loc))        % Provision for only one stroke, that too at end of file
                                filename = filename(1:end-4); 
                                fprintf(fid3,'\"*/%s.lab\"\n0\t%f\n%f\t%f',filename,length_wav_file-0.0100000,length_wav_file-0.0100000,length_wav_file);
                                dangerflag = 1;
                            end
                        end
                        if (dangerflag~=1)
                            clear X;
                            filename = filename(1:end-4);
                            assignin('base','filename',filename);
                            assignin('base','stroke_loc',stroke_loc);
                            %fprintf(fid3,'\"*/%s.lab\"\n',filename);
                            index3 = 1;
                            while (index3 <= length(stroke_loc) )
                                %fprintf(fid3,'%f\n',stroke_loc(index3));
                                strokes{ii-2}(index3) = stroke_loc(index3);%somehow  two extra slots appears in stroke's vector that are unexpected
                                
                                index3 = index3 + 1;
                                
                            end
                            %strokes(ii-2) = transpose(strokes(ii-2));
                        end
                    end
                    fclose(fid3);
                    result = table(names,strokes);

                end
                % predefinitions to calculate precision and recall
                cd test_here_gt;
                for current_gt = 3:1:length(gt)
                    gt_name{current_gt-2} = gt(current_gt).name;
                    gt_file = importdata(gt(current_gt).name);
                    gt_file = gt_file(:,1);%just the first column
                    
                    index4 = 1;
                    while(index4 <= length(gt_file))
                        gt_strokes{current_gt-2}(index4) = gt_file(index4);
                        index4 = index4+1;
                    end
                      
                end
                cd ..;
                ground_truth = table(gt_name,gt_strokes);
                
                
            end
    end
    save result result;
    save ground_truth ground_truth;
end
