function output_txt = data_cursor_text_func(obj,event_obj,...
    x,y,m,k)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
output_txt = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end

i = find(pos(1)==x & pos(2)==y,1);
if ~isempty(i)
  output_txt{end+1} = sprintf('(m,k) = (%d,%d)', m(i), k(i));
end
