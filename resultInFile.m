function [] = resultInFile( name_file,column_vector,save_folder )
%resultInFile writes the data to as file and save in the folder that you
%specify
%   Detailed explanation goes here

mkdir(save_folder);
test = sprintf('./%s/%s',string(save_folder),string(name_file));
fid = fopen(test,'w');

for data = 1:length(column_vector)
    fprintf(fid,sprintf('%d\n',column_vector{data})); 
end
fclose(fid);
end

