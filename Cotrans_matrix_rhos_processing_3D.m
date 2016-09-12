%Cotrans_matrix_rhos_processing_3D.m takes tabular formatting SHAPE-Seq
%reactivities and alignment numbers and plots them in a 3D bar chart
%colored by reactivity

%Written by Kyle E. Watters, 2016
%Version 0.0.2
%Last documentation update: 9/12/2016

%Copyright (C) 2016  Julius B. Lucks, Angela M Yu, and Kyle E. Watters.
%All rights reserved.
%Distributed under the terms of the GNU General Public License, see 'LICENSE'.

clear

%CHOOSE ONLY ONE OF THE FOLLOWING COLORING OPTIONS (set desired to 1)
scaled_solid_colors = 0;   %Color each bar a solid color based on it's 2D color (rho value)
scaled_gradient_colors = 1;   %Color bar height in a gradient fashion based on rho
one_color_only = 0; %Color all the bars the same color (light gray default, change near bottom)
bar_edges_on = 1; %If '0', the black edges on the bars are removed

%Set the upper reactivity limit
upper_limit = 4;

%Set a block of figure save options. Switch to '0' to turn off saving.
%Saves automatically one directory level above data location (if rho data selected).
save_plots_fig = 1;
save_plots_tif = 0;
save_plots_svg = 0;

%Get file save base name (will be appended for filenames)
base_filename = 'empty';   %By default this is determined using the reactivity directory name, fill in to choose

%Choose whether to record a movie (a simple zaxis rotation by default)
record_movie = 0;

%Locate where the tabular reactivities file is
[filename_rho,path_rho,indexx] = uigetfile('.txt', 'Select the Location of the Reactivities Data (Rhos)');
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
    cotrans_matrix = textread(matrix_file,'','delimiter','\t','emptyvalue',NaN);
    
    % create plot  
    figure(1)
    fig1 = bar3(cotrans_matrix);
    [rows,columns] = size(cotrans_matrix);
    cmap = colormap(gcf,parula(300));   %Change color map here
    shortest_len = sum(isfinite(cotrans_matrix(1,:)));    %Count the number of columns in 1st row to see how long the shortest length is
    xtickmarks = 0.5:10:columns;  %determine how to do tick marks, centering
    xticklabels = 0:10:columns;  %determine how to do labels 
    first_ylabel = shortest_len;
    if shortest_len < 10   %If very small, start numbering at ten
        first_ylabel = 10;
    elseif mod(shortest_len,10) < 5 && mod(shortest_len,10) ~= 0  %Always round up to the nearest ten
        first_ylabel = first_ylabel + 5;
    end
    first_ylabel = round(first_ylabel/10)*10;
    ytickmarks = 1:10:columns;  %determine how to do tick marks, starting with 20
    yticklabels = first_ylabel:10:columns;
    set(gca,'YTick',ytickmarks)
    set(gca,'YTickLabel',yticklabels)
    set(gca,'XTick',xtickmarks)
    set(gca,'XTickLabel',xticklabels)
    set(gcf, 'Position', [300 300 800 650]);  %Sets the size of the matrix
    xlabel('Nucleotide Position')
    ylabel('Transcript Length')
    zlabel('Reactivity (\rho)')
    
    caxis([0,upper_limit])
    for i = 1:length(fig1)
        zdata = fig1(i).ZData;
        fig1(i).CData = zdata;
        if scaled_solid_colors == 1
            fig1(i).FaceColor = 'flat';
        elseif scaled_gradient_colors == 1
            fig1(i).FaceColor = 'interp'; 
        elseif one_color_only == 1
            fig1(i).FaceColor = [0.9 0.9 0.9];   %This color can be changed to make the whole plot the same color
        end
    end

    if bar_edges_on == 0
        set(fig1, 'EdgeColor', 'none')    %removes all edge colors
    end
    colorbar();                        
    h = colorbar;
    ylabel(h,'Reactivity (\rho)')     

    %Remove NaN values from the upper triangle
    for i = 1:numel(fig1)
        index = logical(kron(isnan(cotrans_matrix(:,i)),ones(6,1)));
        zData = get(fig1(i),'ZData');
        zData(index,:) = nan;
        set(fig1(i),'ZData',zData);
    end
    
    %Create a movie commands for the 3D bars
    if record_movie == 1
        delete(h)               %Remove the colorbar, as it's distracting
        v = VideoWriter('saved_cotrans_movie.avi');  %Make a video writer object 
        elevation = 55;        %Starting elevation (degrees above x-axis)
        azimuth = -10;        %Starting azimuth (degrees about z-axis)
        view(azimuth,elevation)         %Starting view orientation
        camdolly(-.2,0,0)     %Starting camera location (moved to)
        step_factor = 10;
        frames_per_slide = 24*step_factor;
        open(v)
        for i = 1:frames_per_slide          %Zoom and pan to upper part
            camdolly(-.015/step_factor,0.005/step_factor,0); zoom(1+.05/step_factor)
            frame = getframe(gcf);
            writeVideo(v,frame);    %Sequentially add movie frames
        end
        for i = (2*frames_per_slide):(3*frames_per_slide)  %Slide down the 5' half of the plot
            camdolly(.0125/step_factor,-0.05/step_factor,0)  
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        for i = (3*frames_per_slide):(4.2*frames_per_slide)    %Slide down the 5' half of the plot
            camdolly(0.05/step_factor,0.005/step_factor,0)  
            frame = getframe(gcf);
            writeVideo(v,frame);
        end
        close(v)  %Close the movie
    end
end


%Save plot if requested
save_files = save_plots_fig + save_plots_tif + save_plots_svg;
if save_files > 0 
    if save_plots_fig == 1
       saveas(gcf,strcat(base_filename,'_rhos_3D.fig'),'fig')
   end
   if save_plots_tif == 1
       saveas(gcf,strcat(base_filename,'_rhos_3D.tif'),'tiff')
   end
   if save_plots_svg == 1
       saveas(gcf,strcat(base_filename,'_rhos_3D.svg'),'svg')   
   end
end

