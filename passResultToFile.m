%Generate txt files of the result, just run it if you have result
warning off;
for index5 = 1: length(result.strokes)
    resultInFile(strcat(names{index5}(1:end-4),'.txt'),strokes(index5),'tamborim_result');
end