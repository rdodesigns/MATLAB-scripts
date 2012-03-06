function [ status ] = update_matlab_script()

  root = fileparts(which('update_matlab_script'));
  ignore_dirs = textread([root / '.gitignore'], '%s', 'delimiter', '\n');

  %status = savepath;
  status = 0;

end
