function fre_loc_vis(Split,Names,w,J_w,f)
% This file is to visualize the location of recovered frequencies.
% 
% Inputs
% - Split: 0 or 1. If 1, split the window into several pieces, and plot each piece in a
% separate window; if 0, plot the whole domain in one window
% - Names: a string or a structure containing names for splitted windows
% - w: frequency grids
% - J_w: iamge function values
% - f: location of true frequencies
%
% Return
% - figures: save the plotted figures as pdf
% 
% Created by JYI, 10/24/2019

%%

if Split==1

    figure; 
    subplot(2,2,1); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.00 0.125]);
    
    subplot(2,2,2); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.125 0.250]);
    
    subplot(2,2,3); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.250 0.375]);
    
    subplot(2,2,4); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.375 0.5]);
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3),fig_pos(4)];
    print(fig,Names.F1,'-dpdf'); % 'Noiseless_Hankel_FreIden_P1'
    close all;

    figure; 
    subplot(2,2,1); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.50 0.625]);
    
    subplot(2,2,2); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.625 0.750]);
    
    subplot(2,2,3); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.750 0.875]);
    
    subplot(2,2,4); hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0.875 1.00]);
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3),fig_pos(4)];
    print(fig,Names.F2,'-dpdf'); % 'Noiseless_Hankel_FreIden_P2'
    close all;
else
    figure; hold on
    plot(w,J_w,'-');
    xlabel('Frequency f');
    ylabel('Imaging function J(f)');
    plot(f,zeros(size(f)),'+');
    xlim([0 1.00]);
    
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3),fig_pos(4)];
    print(fig,Names,'-dpdf'); % 'Noiseless_Hankel_FreIden_P2'
    close all;
    
end
% 
% if Zoom==1
% 
% 
% end


end