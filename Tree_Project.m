%Implementing and Statistically Analyzing Protein Trees - by Mehri Haji Abdolghafar

close all;
clear all;
clc;

%Read tables from excel files
data1 = readtable('Datasets\Basket_Zelle_Axon_Infos_Links.xlsx','sheet', 'NodeAnalysis');
data2 = readtable('Datasets\Basket_Zelle_Axon_Infos_Links.xlsx','sheet', 'Nearest Terminal (1)');


%Build a table
numOfData = size(data1, 1);

sorted_node = table('Size', [numOfData 5], 'VariableTypes', {'string', 'double', 'string', 'double', 'double'}, 'VariableNames',...
    {'Node', 'node_index', 'Coordinate', 'node_index_in_tree', 'Distance_from_Root'});

for i=1 : numOfData
    sorted_node{i, 2}= i;
end
sorted_node{:, 1} = data1{:, 2};
sorted_node = sortrows(sorted_node);

%Build the tree, calculate the distances from root & edges length
root_crd = data2{sorted_node{1, 2}, 2};
sorted_node{1, 3} = root_crd;
t = tree(root_crd);
parentNode = 1;
sorted_node{1, 4} = 1;
sorted_node{1, 5} = 0;
char_root_crd = char(root_crd);
root_crd_num = str2double(strsplit(char_root_crd(1, 2:end-1), ','));

edgel = table('Size', [numOfData-1, 7], 'VariableTypes', {'double', 'double', 'double', 'string', 'string', 'double', 'double'},....
    'VariableNames', {'Distance_between_node_and_its_parent', 'node_index_in_tree', 'parent_index_in_tree', ...
    'Node_Coordinate', 'Parent_Coordinate', 'Node_index', 'Parent_index'});

for j=2 : numOfData
    length = strlength(sorted_node{j, 1});
    
    sn = char(sorted_node{j, 1});
    pr = sn(1,1:length - 2);
    
    for s=1 : numOfData
        ch = sorted_node{s, 1};
        if strcmp(ch, pr)
            break;
        end
    end
    
    
    parentNode = sorted_node{s, 4};
    crd = data2{sorted_node{j, 2}, 2};
    sorted_node{j, 3} = crd;
    [t sorted_node{j, 4}] = t.addnode(parentNode, crd);
    
    %Calculate the distances from root
    char_crd = char(crd);
    crd_num = str2double(strsplit(char_crd(1, 2:end-1), ','));
    
    dis_from_root = sqrt (sum ((minus(crd_num, root_crd_num)) .^ 2));
    sorted_node{j, 5} = dis_from_root;
    
    %Calculate edges length
    prt_char_crd = char(sorted_node{s, 3});
    prt_crd_num = str2double(strsplit(prt_char_crd(1, 2:end-1), ','));
    
    edgel{j-1, 1} = sqrt (sum ((minus(crd_num, prt_crd_num)) .^ 2));
    edgel{j-1, 2} = sorted_node{j, 4};
    edgel{j-1, 3} = parentNode;
    edgel{j-1, 4} = sorted_node{j, 3};
    edgel{j-1, 5} = sorted_node{s, 3};
    edgel{j-1, 6} = sorted_node{j, 2};
    edgel{j-1, 7} = sorted_node{s, 2};
end
disp(t.tostring)

%Count the number of edges at each level
levelContent = flatten(t);
level = size(levelContent, 1) - 1;
levelNum = zeros(level, 1);
for k=1 : level
    levelNum(k) = size(levelContent{k+1}, 2);
end

%Find Distances between leaves and their parents
n = 0;
dis_leaf_parent = table('Size', [numOfData 7], 'VariableTypes', {'double', 'double', 'double', 'string', 'string', 'double', 'double'},....
    'VariableNames', {'Distance_between_leaf_and_its_parent', 'leaf_index_in_tree', 'parent_index_in_tree', ...
    'Leaf_Coordinate', 'Parent_Coordinate', 'Leaf_index', 'Parent_index'});
for m=1 : numOfData
    if(t.isleaf(m))
        n = n + 1;
        
       ind = find(edgel{:, 2} == m);
       
       dis_leaf_parent{n, 1} = edgel{ind, 1};
       dis_leaf_parent{n, 2} = m;
       dis_leaf_parent{n, 3} = edgel{ind, 3};
       dis_leaf_parent{n, 4} = edgel{ind, 4};
       dis_leaf_parent{n, 5} = edgel{ind, 5};
       dis_leaf_parent{n, 6} = edgel{ind, 6};
       dis_leaf_parent{n, 7} = edgel{ind, 7};
       
    end
end

dis_leaf_parent(n+1:end, :) = [];

%********Figures*********
%Distance from Root Histogram
figure('Name', 'Distance from Root Histogram and Kernel Density Function')
dis_from_root_hist = histfit(sorted_node{:,5}, 25, 'kernel');
xlabel('Distance from Root')
ylabel('Frequency')
dis_from_root_hist_XRange = dis_from_root_hist(1).XData;
dis_from_root_hist_YRange = dis_from_root_hist(1).YData;

%Histogram of the number of edges at each level
figure('Name', 'Histogram of the Number of Edges at each level')
levelNum_hist = bar(1:level, levelNum);
xlabel('Level')
ylabel('Number of Edges')

%Discrete probability of the number of edges at each level
np = levelNum ./ (numOfData - 1);
figure('Name', 'Discrete probability of the number of edges at each level')
levelNum_p = plot(1:level, np, '*');
xlabel('Level')
ylabel('Edge Distribution of each level')

%Edges Length Histogram
figure('Name', 'Edges Length Histogram and Kernel Density Function')
edgel_hist = histfit(edgel{:, 1}, 25, 'kernel');
xlabel('Edges Length')
ylabel('Frequency')
edgel_hist_XRange = edgel_hist(1).XData;
edgel_hist_YRange = edgel_hist(1).YData;

%Distances between leaves and their parents Histogram
figure('Name', 'Distances between leaves and their parents Histogram and Kernel Density Function')
dis_leaf_parent_hist = histfit(dis_leaf_parent{:, 1}, 10, 'kernel');
xlabel('Distance between leaf and its parent')
ylabel('Frequency')
dis_leaf_parent_hist_XRange = dis_leaf_parent_hist(1).XData;
dis_leaf_parent_hist_YRange = dis_leaf_parent_hist(1).YData;
