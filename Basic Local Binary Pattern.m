%==========================================================================
%==============================Introduction================================
%==========================================================================
% Subject:
%    Face Recognition using Local Binary Pattern
% Author:
%    Pouya Hosseini <970131010>
%    hosseini.pouya7279@gmail.com
%    University of Yasuj
%    Summer 2022
%==========================================================================

% Close all figures, clear command window & workspace
close all;
clear;
clc;

%==========================================================================
%============================Database Images===============================
%==========================================================================
% Folder of image database
Folder = "Image Datasets\Yale\Database Images\";

% Get list of all jpg files in the folder
Filelist = dir(fullfile(Folder,'*.gif'));

% Get number of images in the folder
Num_of_files = numel(Filelist);

% Store images into a cell
Cell_folder = "Image Datasets\Yale\Database Images\";
Cell_filelist = dir(fullfile(Cell_folder,'*.gif'));
Cell_num_of_files = numel(Cell_filelist);

% Preallocate array & cell
Euclidean_distance = zeros(Num_of_files,1);
Cell_read = cell(Cell_num_of_files,1);
Cell_LBP = cell(Cell_num_of_files,1);
Cell_histogram = cell(Cell_num_of_files,1);
Cell_euclidean_distance = cell(Cell_num_of_files,1);
Cell_query = cell(Cell_num_of_files,1);
Cell_location = cell(Cell_num_of_files,1);

% Read images in  the cell
for k = 1 : Cell_num_of_files
    Cell_read{k} = imread(fullfile(Cell_folder,Cell_filelist(k).name));
end

% Read images in the database
for k = 1 : Num_of_files
	Base_file_name = Filelist(k).name;
	Full_file_name = fullfile(Filelist(k).folder,Base_file_name);
	fprintf('Processing image #%d of %d : "%s"\n', ...
             k,Num_of_files,Base_file_name);
	Base_image = imread(Full_file_name);

	% Convert image to double
    Image = im2double(Base_image);

    % Convert image to grayscale if RGB
    if size(Image,3) == 3
        Image = rgb2gray(Image);
    end

    % Get the dimensions of the image
    [Rows_1,Columns_1] = size(Image);

    % Preallocate array
    LBP = zeros(size(Image),'uint8');
    Eight_bit_number = zeros(1,Num_of_files);
    
    % Local binary patterning algorithm 
    % Size of window ====> 3*3
    for i_2 = 2 : Rows_1 - 1
        for j_2 = 2 : Columns_1 - 1

            % Find out center pixel
            Center_pixel = Image(i_2,j_2);

            % Pixel quantization
		    % Upper left
		    Pixel_7 = Image(i_2-1,j_2-1) >= Center_pixel;
            % Upper middle
		    Pixel_6 = Image(i_2-1,j_2) >= Center_pixel;
            % Upper right
		    Pixel_5 = Image(i_2-1,j_2+1) >= Center_pixel;
            % Middle right
		    Pixel_4 = Image(i_2,j_2+1) >= Center_pixel;
            % Lower right
		    Pixel_3 = Image(i_2+1,j_2+1) >= Center_pixel;
            % Lower middle
		    Pixel_2 = Image(i_2+1,j_2) >= Center_pixel;
            % Lower left
		    Pixel_1 = Image(i_2+1,j_2-1) >= Center_pixel;
            % Middle left
		    Pixel_0 = Image(i_2,j_2-1) >= Center_pixel;
		
            % Find out eight bit number
		    Eight_bit_number(k) = uint8(...
			Pixel_7 * 2^7 + Pixel_6 * 2^6 + ...
			Pixel_5 * 2^5 + Pixel_4 * 2^4 + ...
			Pixel_3 * 2^3 + Pixel_2 * 2^2 + ...
			Pixel_1 * 2 + Pixel_0);

            % LBP image
		    LBP(i_2,j_2) = Eight_bit_number(k);
            Cell_LBP{k} = LBP;
        end
    end

%--------------------------------------------------------------------------
% Divide database LBP
% Figure out where to divide it
[Rows_2,Columns_2] = size(LBP);
Rows_3 = round(Rows_2 / 3);
Column_3 = round(Columns_2 / 3);

% Divide image to 9 blocks (3*3)
% Upper left
Block_1 = LBP(1:Rows_3,1:Column_3);

% Upper middle
Block_2 = LBP(1:Rows_3,Column_3+1:2*Column_3);

% Upper right
Block_3 = LBP(1:Rows_3,2*Column_3+1:end);

% Middle left
Block_4 = LBP(Rows_3+1:2*Rows_3,1:Column_3);

% Center
Block_5 = LBP(Rows_3+1:2*Rows_3,Column_3+1:2*Column_3);

% Middle right
Block_6 = LBP(Rows_3+1:2*Rows_3,2*Column_3+1:end);

% Lower left
Block_7 = LBP(2*Rows_3+1:end,1:Column_3);

% Lower middle
Block_8 = LBP(2*Rows_3+1:end,Column_3+1:2*Column_3);

% Lower right
Block_9 = LBP(2*Rows_3+1:end,2*Column_3+1:end);

% Histogram of each block
Histogram_1 = imhist(uint8(Block_1));
Histogram_2 = imhist(uint8(Block_2));
Histogram_3 = imhist(uint8(Block_3));
Histogram_4 = imhist(uint8(Block_4));
Histogram_5 = imhist(uint8(Block_5));
Histogram_6 = imhist(uint8(Block_6));
Histogram_7 = imhist(uint8(Block_7));
Histogram_8 = imhist(uint8(Block_8));
Histogram_9 = imhist(uint8(Block_9));

% Total histogram
Total_histogram = [Histogram_1(:);Histogram_2(:);Histogram_3(:);
                   Histogram_4(:);Histogram_5(:);Histogram_6(:);
                   Histogram_7(:);Histogram_8(:);Histogram_9(:)];

Cell_histogram{k} = Total_histogram';

% Concatenate arrays within a loop
All_total_histogram = vertcat(Cell_histogram{:});

% Show images (if you want to show images uncomment it!)
% figure("Name",['Image Number ',num2str(k)]);
% set(gcf,'position',[365,110,800,600]);
% 
% subplot(2,2,1);
% imshow(Base_image)
% title(['Image Number ',num2str(k)]);
% 
% subplot(2,2,2);
% imshow(uint8(LBP))
% title('LBP Image');
% 
% subplot(2,2,(3:4));
% bar(Total_histogram,'BarWidth',1,'EdgeColor','none');
% title('LBP Histogram');

end

%==========================================================================
%===========================Matching Process===============================
%==========================================================================
fprintf('\n')
Num_of_top = input('>>Please enter the number of top images: ');

% Prelocate array & cell
Location_of_top = cell(1,Num_of_files);
Class = zeros(1,Num_of_files);
Cell_counter = zeros(1,Num_of_files);
Percentage_of_each_image = zeros(1,Num_of_files);

% Making a class for comparison
for c = 1:165
    Class(c) = c;
end

Class = repelem(Class,11);

% Find euclidean distance & accuracy
for m = 1 : Num_of_files
    Query = All_total_histogram(m,:);
    Cell_query{m,1} = Query;

    for n = 1 : Num_of_files
        Euclidean_distance(n,1) = sqrt(sum((Query- ...
                                  All_total_histogram(n,:)).^2));
    end

    Cell_euclidean_distance{m,1} = Euclidean_distance;
    [Sorted_values,Location] = sort(Euclidean_distance);
    Cell_location{m,1} = Location;

% Find & show matching images (if you want to show images uncomment it!)
%      for i = 1 : Num_of_top
%          figure("Name",['Matching Database Image Number ',num2str(i)]);
%          set(gcf,'position',[600,300,800,600]);
%          imshow(Cell{Location(i,1)})
%          title(['Matching Database Image Number ',num2str(i)]);
%      end

%==========================================================================
%========================Accuracy of Algorithm=============================
%==========================================================================
     % Put location of top images in a cell
     Location_of_top{m}=Cell_location{m,1}(1:Num_of_top)';
     
     % Set a counter for accuracy
     Counter = 0;

     for z = 1:Num_of_top

         if Class(Location_of_top{1,m}(z)) == Class(m)
             Counter = Counter + 1;
         end

         Cell_counter(1,m) = Counter;
     end

     Percentage_of_each_image(m) = (Counter/Num_of_top)*100;
end

Total_accuracy_percentage = mean(Percentage_of_each_image);
fprintf('The Total Percentage is %0.1f \n',Total_accuracy_percentage)