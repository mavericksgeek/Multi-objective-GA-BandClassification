clc;
clear all;
%%
region = 'salinasA_corrected';
class_namelist = {'class1', 'class2', 'class3', 'class4', 'class5', 'class6', 'noClass'};
colormap_table = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];
% region = 'paviaU';
% class_namelist = {'Asphalt', 'BareSoil', 'Bitumen', 'Gravel', 'Meadows', 'PaintedMetalSheets', 'SelfBlockingBricks', 'Shadows', 'Trees'};
% colormap_table = [0.5 0.5 0.5; 0.4 0.2 0; 0 0 1; 0 1 1; 0 0.7 0; 1 0 0; 0.4 0 0; 0 0 0; 0 0.4 0];

%%
no_of_class = length(class_namelist) - 1;
load(strcat(region, '.mat'));
original_image = salinasA_corrected;
% original_image = paviaU;
[no_of_row, no_of_col, no_of_band] = size(original_image);
no_of_pixel = no_of_col* no_of_row;
TestSet = double(reshape(original_image, no_of_pixel, no_of_band));

%%
for i = 1 : no_of_class
    trainingData{i} = xlsread(strcat(class_namelist{i}, 'DataTraining.xlsx'));
    trainingLabel{i} = repmat(i, [1 size(trainingData{i},1)]);
end
result = zeros(length(TestSet(:,1)),1);

%% build models
k = 0;
for i = 1 : no_of_class-1
    for j = i+1 : no_of_class
        k = k+1;
        TrainingSet = [trainingData{i}; trainingData{j}];
        GroupTrain = [trainingLabel{i}'; trainingLabel{j}'];
        models(k) = svmtrain(TrainingSet, GroupTrain, 'kernel_function','rbf','rbf_sigma', 10);
    end
end

% classify test cases
parfor j=1:size(TestSet,1)
    for i=1:k
        int_result(j,i) = svmclassify(models(i),TestSet(j,:));
    end
end
[result, f] = mode(int_result, 2);
%%
classified_image = reshape(result, no_of_row, no_of_col);

h = figure;
imshow(classified_image, colormap_table);
save('classified_image.mat', 'classified_image');
