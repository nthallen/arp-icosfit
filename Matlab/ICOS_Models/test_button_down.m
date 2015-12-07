function test_button_down(src, ~, arg3)
  srctype = get(src,'type');
  fprintf(1,'Button down on %s\n', srctype);
  if nargin == 3
    fprintf(1,'arg3 is %d\n', arg3);
  else
    fprintf(1,'Only two args\n');
  end
  
  while ~strcmp(srctype,'axes')
    src = get(src,'parent');
    if isempty(src)
      fprintf(1,'No axes found\n');
      return;
    end
    srctype = get(src,'type');
  end
  axesHandle  = get(src,'Parent');
  coordinates = get(axesHandle,'CurrentPoint');
  coordinates = coordinates(1,1:2);
  fprintf(1,'Click at coordinates: %f, %f\n', coordinates(1), coordinates(2));
