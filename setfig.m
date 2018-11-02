% This script sets the property of figure.
% Prease write after making all the texts and diagrams.
% 
% 
%   background color
%   set font
%   % colormap
% 

% figure handle
fig = gcf;

% background color
fig.Color = 'black';

% set font
% list = listfonts      % available fonts
myFontName = 'ÉÅÉCÉäÉI';
myFontSize = 16;
myAxisFontSize = 16;
myFontColor = 'white';
set(findall(fig,'type','axes'),'FontName', myFontName,'fontsize',myFontSize,'Color',myFontColor)
set(findall(fig,'type','text'),'FontName', myFontName,'fontSize',myAxisFontSize,'Color',myFontColor) 

% colormap
% myColormap = jet;
% colormap(myColormap);

% end