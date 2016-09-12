%Cotrans_matrix_rhos_processing_2D.m takes tabular formatting SHAPE-Seq
%reactivities and alignment numbers and plots them on a heatmap-type plot

%Written by Kyle E. Watters, 2016
%Version 0.0.2
%Last documentation update: 9/12/2016

%Copyright (C) 2016  Julius B. Lucks, Angela M Yu, and Kyle E. Watters.
%All rights reserved.
%Distributed under the terms of the GNU General Public License, see 'LICENSE'.

clear
load('jet_hires_white_bckgd')

%Choose upper limit of rho plot
upper_lim = 4;

%Set a block of figure save options. Switch to '0' to turn off saving.
%Saves automatically one directory level above data location (if rho data selected).
save_plots_fig = 0;
save_plots_tif = 0;
save_plots_svg = 0;

%Get file save base name (will be appended for filenames)
base_filename = 'empty';   %By default this is determined using the reactivity directory name, fill in to choose

%Locate where the tabular reactivities file is
[filename_rho,path_rho,indexx] = uigetfile('.txt', 'Select the Location of the Reactivities Data (Rhos)');

file_skips = 0;
if (filename_rho ~= 0 )  %skip plotting the cotrans table if hit cancel
    %Build the base name for file saving if not provided
    if base_filename == 'empty'
        base_filename = '/';
        directories = strsplit(path_rho,'/');
        directories = directories(~cellfun('isempty',directories)); 
        for path = 1:size(directories,2)
            base_filename = strjoin([base_filename, directories(path)],'/');
        end
    end
    
    matrix_file = strcat(path_rho,filename_rho);
    
    %Load the matrix data into arrays
    raw_matrix = textread(matrix_file,'','delimiter','\t','emptyvalue',NaN);
    [rows,columns] = size(raw_matrix);
    %Pad an extra row and column of 0s so that the pcolor vertices are
    %plotted in the proper location (first/last dropped otherwise)
    pad1 = zeros(1,columns);
    pad1(pad1 == 0) = NaN;
    pad2 = zeros(rows+1,1);
    pad2(pad2 == 0) = NaN;
    cotrans_matrix = [pad1 ; raw_matrix];
    cotrans_matrix = horzcat(cotrans_matrix,pad2);
    
    % create plot  
    figure(1)
    fig1 = pcolor(flipud(cotrans_matrix));
    %colormap(jet(300))
    colormap(parula(300))
    shortest_len = sum(isfinite(raw_matrix(1,:)));    %Count the number of columns in 1st row of the original matrix to see how long the shortest length is
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
    caxis([0,upper_lim]);              %Change this value to adjust the colormap limits/scaling
    h = colorbar();
    ylabel(h,'Reactivity (rho)')
    next_path = strcat(path_rho,'.txt');

    
else
    next_path = '.txt';
    file_skips = file_skips + 1;
end

%Locate where the tabular treated aligned reads file is
[filename_1,path_1,indexx] = uigetfile(next_path, 'Select the Location of the Treated Aligned Reads Data');
if (filename_1 ~= 0 )  %skip plotting the cotrans table if hit cancel
    treated_file = strcat(path_1,filename_1);
    %Load the matrix data into arrays
    treated_matrix = textread(treated_file,'','delimiter','\t','emptyvalue',NaN);
    log_treated_matrix = log10(treated_matrix);
    [rows,columns] = size(treated_matrix);
    %Pad an extra row and column of 0s so that the pcolor vertices are
    %plotted in the proper location (first/last dropped otherwise)
    pad1 = zeros(1,columns);
    pad1(pad1 == 0) = NaN;
    pad2 = zeros(rows+1,1);
    pad2(pad2 == 0) = NaN;
    log_treated_matrix = [pad1 ; log_treated_matrix];
    log_treated_matrix = horzcat(log_treated_matrix,pad2);
    
    % create plot
    figure(2)
    fig2 = pcolor(flipud(log_treated_matrix));
    colormap(jet_hires_white_bckgd)  %Need this colormap so that the spots with no reads appear white
    %colormap(jet(300))
    %colormap(parula(300))
    shortest_len = sum(isfinite(treated_matrix(1,:)))-1;    %Count the number of columns in 1st row to see how long the shortest length is (-1 from RT offset)
    xtickmarks = horzcat(2.5,11.5:10:columns);  %determine how to do tick marks, centering (shifted 1 to right for RT offset)
    xticklabels = horzcat(1, 10:10:columns);
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
    set(fig2, 'EdgeColor', 'none')
    set(gcf, 'Position', [300 300 800 650]);   %Sets the size of the matrix
    xlabel('Nucleotide Position')
    ylabel('Transcript Length')
    caxis([0,6]);                %Change this value to adjust the colormap limits/scaling
    colorbar()
    h = colorbar;
    ylabel(h,'Log(Number of Treated Reads Aligned)')
    
    next_path = strcat(path_1,'.txt');
else
    file_skips = file_skips + 1;
end

%Locate where the tabular untreated aligned reads file is
[filename_2,path_2,indexx] = uigetfile(next_path, 'Select the Location of the Untreated Aligned Reads Data');
if (filename_2 ~= 0 )  %skip plotting the cotrans table if hit cancel
    untreated_file = strcat(path_2,filename_2);
    %Load the matrix data into arrays
    untreated_matrix = textread(untreated_file,'','delimiter','\t','emptyvalue',NaN);
    log_untreated_matrix = log10(untreated_matrix);
    [rows,columns] = size(untreated_matrix);
    %Pad an extra row and column of 0s so that the pcolor vertices are
    %plotted in the proper location (first/last dropped otherwise)
    pad1 = zeros(1,columns);
    pad1(pad1 == 0) = NaN;
    pad2 = zeros(rows+1,1);
    pad2(pad2 == 0) = NaN;
    log_untreated_matrix = [pad1 ; log_untreated_matrix];
    log_untreated_matrix = horzcat(log_untreated_matrix,pad2);  
    
    % create plot
    figure(3)
    fig3 = pcolor(flipud(log_untreated_matrix));
    colormap(jet_hires_white_bckgd)  %Need this colormap so that the spots with no reads appear white
    %colormap(jet(300))
    %colormap(parula(300))
    [rows,columns] = size(untreated_matrix);
    shortest_len = sum(isfinite(untreated_matrix(1,:)))-1;    %Count the number of columns in 1st row to see how long the shortest length is (-1 from RT offset)
    xtickmarks = horzcat(2.5,11.5:10:columns);  %determine how to do tick marks, centering (shifted 1 to right for RT offset)
    xticklabels = horzcat(1, 10:10:columns);  %determine how to do labels 
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
    set(fig3, 'EdgeColor', 'none')
    set(gcf, 'Position', [300 300 800 650]);  %Sets the size of the matrix
    xlabel('Nucleotide Position')
    ylabel('Transcript Length')
    caxis([0,6]);             %Change this value to adjust the colormap limits/scaling
    colorbar()
    h = colorbar;
    ylabel(h,'Log(Number of Untreated Reads Aligned)')
else
    file_skips = file_skips + 1;
end

%Save plot if requested
save_files = save_plots_fig + save_plots_tif + save_plots_svg;
if save_files > 0 && file_skips < 3
    if save_plots_fig == 1
       saveas(fig1,strcat(base_filename,'_rhos.fig'))
       saveas(fig2,strcat(base_filename,'_treated_reads.fig'))
       saveas(fig3,strcat(base_filename,'_untreated_reads.fig'))
   end
   if save_plots_tif == 1
       saveas(fig1,strcat(base_filename,'_rhos.tif'))
       saveas(fig2,strcat(base_filename,'_treated_reads.tif'))
       saveas(fig3,strcat(base_filename,'_untreated_reads.tif'))
   end
   if save_plots_svg == 1
       saveas(fig1,strcat(base_filename,'_rhos.svg'))
       saveas(fig2,strcat(base_filename,'_treated_reads.svg'))
       saveas(fig3,strcat(base_filename,'_untreated_reads.svg'))    
   end
end
