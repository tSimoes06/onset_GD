function [precision, recall] = precisionRecall( ground_truth_strokes,result_strokes)
% precisionRecall calculates precision and recall for the result_strokes of
% one instrument,

true_positive = 0;
false_positive = 0;
false_negative = 0;

total_of_tracks = length(result_strokes);

for current_track = 1:total_of_tracks
    
    size_GT_strokes = length(ground_truth_strokes{current_track});
    size_inst_strokes = length(result_strokes{current_track});
    
    for data = 1:size_GT_strokes
        
        find = false;
        data_count = 1; %to reduce the number of loops
        
        while ~find && (data_count < size_inst_strokes)
            
            if (abs(ground_truth_strokes{current_track}(data)-result_strokes{current_track}(data_count))) <= 0.02
                true_positive = true_positive +1;
                find = true;
            end
            data_count = data_count + 1;
        end
        
        if ~find
            false_negative = false_negative + 1;
        end    
    end        
    false_positive = false_positive + abs(((size_GT_strokes - false_negative) - size_inst_strokes));
end
precision = true_positive/(true_positive + false_positive);
recall = true_positive/(true_positive + false_negative);


end

