function cl_fig(varargin)

% CL_FIG closes figure windows
%
%   CL_FIG(A) closes all the figure windows specified in array A
%
%   CL_FIG on its own closes all currently open figure windows
%
%   This function is of most use to close specific windows when using FIGURE on its
%   own to open sequential windows to prevent window overload!
%
%   Mark Humphries 16/4/04

if nargin == 0
    win_list = get(0,'Children');
elseif nargin == 1
    win_list = varargin{1};
else
    error('Too many inputs')
end

% close specified windows
delete(win_list);