function [ absorption ] = CEAS_absorption( background_file_path, ...
                                           solvent_file_path, ...
                                           solution_file_path, ...
                                           reflectivity, sample_length, ...
                                           wavelength_min, wavelength_max )
%   Calculate absorption coefficient from CEAS experiments
%
%   The calculation of the absorption coefficient (alpha, in cm^-1), is done
%   through the following formula
%
%   a_v = ((I_0 - I)/I)_v * (1-R)/d
%
%   where
%
%     I = Intensity with solute.
%   I_0 = Intensity without solute.
%     R = Reflectivity of the mirrors.
%     d = Effective Path Length (cm).
%     v = Wavelength(nm).
%
%   This function automatically finds the area where above background
%   intensity was detected and uses this to create a bounding box for the
%   absorption curve. Additionally the bounding box is used to average the
%   spread of the absorption line.
%
%   --------------------------------------------------------------------------
%   The input parameters are
%
%   background_file_path: file path (relative or absolute) to the background
%                         image.
%      solvent_file_path: file path (relative or absolute) to the solvent only
%                         image.
%     solution_file_path: file path (relative or absolute) to the image
%                         containing the solute.
%           reflectivity: The reflectivity of the mirrors.
%            sample_length: The effective path length.
%         wavelength_min: The wavelength value for the far left pixel.
%         wavelength_max: The wavelength value for the far right pixel.
%
%   *IMPORTANT*: The files must be oriented so that the wavelengths are spread
%   out in the horizontal direction and the minimum wavelength is on the left.
%   Additionally it is assumed that the dimensions on all of the files are the
%   same.
%
%   --------------------------------------------------------------------------
%   The output is a matrix where the first column is a wavelength value and
%   the second is the absorption coefficient at that wavelength in cm^-1.


  %% Get the info structs on each of the files.
  info_background = imfinfo(background_file_path);
  info_solvent    = imfinfo(solvent_file_path);
  info_solution   = imfinfo(solution_file_path);

  %% Find the width and height of the images.
  height = info_background(1).Height;
  width  = info_background(1).Width;

  %% Find the number of pages in each of the files.
  pages_background = numel(info_background);
  pages_solvent    = numel(info_solvent);
  pages_solution   = numel(info_solution);

  %% Allocate space for background, solvent and solution images.
  background   = zeros(height, width);
  solvent_raw  = zeros(height, width);
  solution_raw = zeros(height, width);

  %% Average the pages in the background file.
  for k = 1:pages_background
    A = double(imread(background_file_path, k, 'Info', info_background));
    background = background + (A./pages_background);
  end

  %% Average the pages in the solvent file.
  for k = 1:pages_solvent
    A = double(imread(solvent_file_path, k, 'Info', info_solvent));
    solvent_raw = solvent_raw + ((A - background)./pages_solvent);
  end

  %% Average the pages in the solution file.
  for k = 1:pages_solution
    A = double(imread(solution_file_path, k, 'Info', info_solution));
    solution_raw = solution_raw + ((A - background)./pages_solution);
  end

  %% Find the edges of the collected data
  edges = findBoxContainingAllSegments(createSegmentedImage(solvent_raw));
  x_min = edges(1);
  x_max = edges(1) + edges(3);
  y_min = edges(2);
  y_max = edges(2) + edges(4);

  %% Slice images into pieces
  solvent  = mean(solvent_raw(y_min:y_max, x_min:x_max));
  solution = mean(solution_raw(y_min:y_max, x_min:x_max));

  %% Calculate Absorption Coefficient
  wavelength_range = wavelength_max - wavelength_min;
  absorption = zeros(length(solvent), 2);

  % Create wavelength stepping array, then offset by initial x_min wavelength.
  absorption(:,1) = ((0:length(absorption) - 1) * (wavelength_range/width)) ...
                    + (double(x_min)*(wavelength_range/width) + wavelength_min);
  absorption(:,2) = ((solvent - solution)./solution) ...
                    * ((1-reflectivity)/sample_length);
end % CEAS_absorption

function [ segmented_image ] = createSegmentedImage(image)
  %% Detect Contrast
  [~, threshold] = edge(image, 'sobel');
  fudge_factor = 0.1;
  contrast_img = edge(image, 'sobel', threshold*fudge_factor);

  %% Dilate the image
  se90 = strel('line', 1, 90);
  se0 = strel('line', 1, 0);
  dilated_contrast_img = imdilate(contrast_img, [se90 se0]);

  %% Fill in the holes
  segmented_image = imfill(dilated_contrast_img, 'holes');

end % createSegmentedImage

function [ bounding_box ] = findBoxContainingAllSegments(image)
  %% Get the properties of the segments
  segment_info = regionprops(image, {'Area', 'BoundingBox'});

  %% Find largest bounding box
  bounding_box = segment_info(1).BoundingBox;
  for k=2:length(segment_info)
    if segment_info(k).Area < 40
      continue
    end

    segment_info_box = segment_info(k).BoundingBox;
    new_bounding_box = zeros(1,4);

    %% Is the x direction more left in the new box?
    if segment_info_box(1) < bounding_box(1)
      new_bounding_box(1) = segment_info_box(1);
    else
      new_bounding_box(1) = bounding_box(1);
    end

    %% Is the y direction more up in the new box?
    if segment_info_box(2) < bounding_box(2)
      new_bounding_box(2) = segment_info_box(2);
    else
      new_bounding_box(2) = bounding_box(2);
    end

    %% Is the x direction more right in the new box?
    if (segment_info_box(3) + segment_info_box(1)) > ...
       (bounding_box(3) + bounding_box(1))
      new_bounding_box(3) = (segment_info_box(3) + segment_info_box(1)) ...
                            - new_bounding_box(1);
    else
      new_bounding_box(3) = (bounding_box(3) + bounding_box(1)) ...
                            - new_bounding_box(1);
    end

    %% Is the y direction more down in the new box?
    if (segment_info_box(4) + segment_info_box(2)) > ...
       (bounding_box(4) + bounding_box(2))
      new_bounding_box(4) = (segment_info_box(4) + segment_info_box(2)) ...
                            - new_bounding_box(2);
    else
      new_bounding_box(4) = (bounding_box(4) + bounding_box(2)) ...
                            - new_bounding_box(2);
    end

    % bounding_box contains all currently processed boxes.
    bounding_box = new_bounding_box;
  end

  % Sometimes regionprops returns a fractional pixel, undo this.
  bounding_box = int32(bounding_box);
end % findBoxContainingAllSegments
