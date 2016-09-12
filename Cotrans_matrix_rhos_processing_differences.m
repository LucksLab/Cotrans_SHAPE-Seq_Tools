%Cotrans_matrix_rhos_processing_differences.m takes two tabular formatting SHAPE-Seq
%reactivities and alignment numbers and plots their differences on a heatmap-type plot

%Written by Kyle E. Watters, 2016
%Version 0.0.2
%Last documentation update: 9/12/2016

%Copyright (C) 2016  Julius B. Lucks, Angela M Yu, and Kyle E. Watters.
%All rights reserved.
%Distributed under the terms of the GNU General Public License, see 'LICENSE'.


clear 

%Choose the plot limits (even)
lims = 5;

%Set a block of figure save options. Switch to '0' to turn off saving.
%Saves automatically one directory level above data location (if rho data selected).
save_plots_fig = 0;
save_plots_tif = 0;
save_plots_svg = 0;

%Get file save base name (will be appended for filenames)
base_filename = 'empty';   %By default this is determined using the reactivity directory name, fill in to choose

%Locate where the first tabular reactivities file is
[filename_rho,path_rho,indexx] = uigetfile('.txt', 'Select the Location of the first reactivity data set (?1 of ?1-?2)');
%Build the base name for file saving if not provided
    if base_filename == 'empty'
        base_filename = '/';
        directories = strsplit(path_rho,'/');
        directories = directories(~cellfun('isempty',directories)); 
        for path = 1:size(directories,2)
            base_filename = strjoin([base_filename, directories(path)],'/');
        end
    end
    
if (filename_rho ~= 0)  %skip plotting the cotrans table if hit cancel
    matrix_file = strcat(path_rho,filename_rho);
    %Load the matrix data into arrays
    cotrans_matrix_1 = textread(matrix_file,'','delimiter','\t','emptyvalue',NaN);
    [rows1,columns1] = size(cotrans_matrix_1);
    next_path = strcat(path_rho,'.txt');
    %Locate where the second tabular reactivities file is
    [filename_rho,path_rho,indexx] = uigetfile(next_path, 'Select the Location of the second reactivity data set (?2 of ?1-?2)');
    if (filename_rho ~= 0 )  %skip plotting the cotrans table if hit cancel
        matrix_file = strcat(path_rho,filename_rho);
        %Load the matrix data into arrays
        cotrans_matrix_2 = textread(matrix_file,'','delimiter','\t','emptyvalue',NaN);
        [rows,columns] = size(cotrans_matrix_2);
        if rows1 ~= rows || columns1 ~= columns
            error('The two selected matrices are not the same size. Please rerun and select two matrices that are the same size.')
        end
        cotrans_matrix_diff = cotrans_matrix_1 - cotrans_matrix_2;
        %Pad an extra row and column of 0s so that the pcolor vertices are
        %plotted in the proper location (first/last dropped otherwise)
        pad1 = zeros(1,columns);
        pad1(pad1 == 0) = NaN;
        pad2 = zeros(rows+1,1);
        pad2(pad2 == 0) = NaN;
        cotrans_matrix_diff = [pad1 ; cotrans_matrix_diff];
        cotrans_matrix_diff = horzcat(cotrans_matrix_diff,pad2);
        
        % create plot
        load('red_blue_hot_cold.mat');
        figure(1)
        fig1 = pcolor(flipud(cotrans_matrix_diff));
        colormap(red_blue_hot_cold);
        shortest_len = sum(isfinite(cotrans_matrix_1(1,:)));    %Count the number of columns in 1st row of the original matrix to see how long the shortest length is
        xtickmarks = 0.5:10:columns+0.5;  %determine how to do tick marks, centering
        xticklabels = 0:10:columns;  %determine how to do labels 
        first_ylabel = shortest_len;
        if shortest_len < 10   %If very small, start numbering at ten
            first_ylabel = 10;
        elseif mod(shortest_len,10) < 5 && mod(shortest_len,10) ~= 0  %Always round up to the nearest ten
            first_ylabel = first_ylabel + 5;
        end
        first_ylabel = round(first_ylabel/10)*10;
        yticklabels = fliplr(first_ylabel:10:(rows+shortest_len));
        ytickmarks = (rows+shortest_len-yticklabels(size(yticklabels))+0.5):10:rows+0.5;  %determine how to do tick marks, starting labels with shortest length, centered
    
        set(gca,'YTick',ytickmarks)
        set(gca,'YTickLabel',yticklabels)
        set(gca,'XTick',xtickmarks)
        set(gca,'XTickLabel',xticklabels)
        set(fig1, 'EdgeColor', 'none')
        set(gcf, 'Position', [300 300 800 650]);  %Sets the size of the matrix
        xlabel('Nucleotide Position')
        ylabel('Transcript Length')
        caxis([-lims,lims]);              %Change this value to adjust the colormap limits/scaling
        h = colorbar();
        ylabel(h,'\Delta Reactivity (rho)')
    else
        
    end
    
    %Save plot if requested
    save_files = save_plots_fig + save_plots_tif + save_plots_svg;
    if save_files > 0
        if save_plots_fig == 1
            saveas(fig1,strcat(base_filename,'_diff.fig'))
        end
        if save_plots_tif == 1
            saveas(fig1,strcat(base_filename,'_diff.tif'))
        end
        if save_plots_svg == 1
            saveas(fig1,strcat(base_filename,'_diff.svg'))
        end
    end
end
