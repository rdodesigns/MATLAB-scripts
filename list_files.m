function [ list ] = list_files(input_string)

  listDir = dir(input_string);
  %[~, order] = sort([listDir(:).datenum]);
  %listDir = listDir(order(end:-1:1));
  if numel(listDir) == 1
    if listDir.isdir == 0
      list = listDir.name;
    else
      list = '';
    end
  else
    list = {listDir(~[listDir.isdir]).name};
  end

end
