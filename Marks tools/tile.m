function tile(varargin)

% TILE tiles figure windows on screen
%
% TILE by itself lays out all the open MatLab figure windows into a grid, starting from the top-left corner
%
% TILE(M,A) places the figures on monitor M (M=1,2), set to [] to omit; 
% tiles just the figure windows specified by A
%
% Mark Humphries 6/07/2006

figure_h = sort(get(0,'Children'));
% screen = get(0,'ScreenSize');
monitors = get(0,'MonitorPositions');

if isempty(figure_h)
    error('No figure windows are open');
end

targets = figure_h;
monitor = 1;
if nargin >= 1 & isnumeric(varargin{1}) monitor = varargin{1}; end
if nargin >= 2 & isnumeric(varargin{2}) targets = varargin{2}; end 

screen = monitors(monitor,:);

num_figs = length(targets);

grid_length = num_figs;
% if number of figures is prime then find first non-prime greater
prime = isprime(grid_length);
while prime
   grid_length = grid_length + 1;
   prime = isprime(grid_length); 
end    

% define number of grid rows and cols: cols >= rows because screen width > height
grid_rows = ceil(grid_length/sqrt(grid_length));
grid_cols = grid_rows;

if grid_rows^2 - grid_length >= grid_rows
    % tidy up a bit
    grid_rows = ceil(grid_length / grid_rows);
end

% tiles are square by default, so use row height to define size
figure_bar = 75;            % size of figure menu bar in pixels
start_bar = 20;
horiz_pad = 5;
vert_pad = 7;

row_height = (screen(4) - start_bar) / grid_rows - figure_bar;
col_width = row_height;

% starting positions for tiles
x_start = horiz_pad + screen(1);
y_start = screen(4) - row_height - figure_bar;

counter = 1;
for loop1 = 1:grid_rows
    for loop2 = 1:grid_cols
        x_pos = x_start + (col_width+horiz_pad) * (loop2-1);                    % calculate new x-axis position
        y_pos = y_start - (row_height+figure_bar+vert_pad) * (loop1-1);         % calculate new y-axis position
        set(targets(counter),'Position',[x_pos y_pos col_width row_height]);    % set position
        figure(targets(counter));                                               % show figure
        counter = counter + 1;
        if counter > num_figs
            break
        end
    end
    if counter > num_figs
        break
    end
end
