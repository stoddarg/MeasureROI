% Measure ROI script
% Description of the function of this script
% input variables, if any
% output variables, if any
clear

% User Defined Settings
% SliceThickness % we find this based on the read in DICOM volume

% defining variables
sliceRegions = [41:50,59:68,77:86,96:105,114:123]; % Ben's slice regions
% sliceRegions = [77:83]; % Aroon's slice regions
numROIs = 9; % Ben's ROIs to be placed
% numROIs = 8; Aroon's ROIs to be placed
imCenter = 256;
ROIRadius = 28.5 * 0.6; % in mm
BBLocations = zeros(4,2);
dist2BB = 100; % in mm
dist2ROIperiph = 75; % in mm
dist2ROIperiphDiag = dist2ROIperiph / sqrt(2);
ROICenterLoc = zeros(9,2);

saveVar = zeros(numel(sliceRegions),numROIs,4);

%% read in the 4D DICOM volume
% ask the user for the location of the DICOM files to read in
% Note to the user: only place 1 series of DICOM files into the folder
i = 1;
while i <= 1
    DirPath = uigetdir(); % pop up a windows explorer box
    i = i + 1;
end

[testVol, spt, dims] = dicomreadVolume(DirPath);
testVolSq = squeeze(testVol);

%% Identify the orientation BB's in the central slice of the phantom
% the BBs are located about 80 mm into the phantom from the front, start
% looking for them at slice 80 (if 1 mm slices) or slice 160 (if 0.5 mm)
SliceThickness = spt.PatientPositions(2,3) - spt.PatientPositions(1,3);

% identify the slices with a location of 78 mm into the phantom
CheckSlices = (spt.PatientPositions(:,3) - spt.PatientPositions(1,3)) >= 78;
CheckSlices = find(CheckSlices > 0, round(5/SliceThickness));

for j = 1:numel(CheckSlices)
    SlicePixels = squeeze(testVolSq(:,:,CheckSlices(j)));

    % threshold to find only the pixel values above ~500/600
    findBBs = SlicePixels;
    findBBs(findBBs < 1000) = 0;

    % use segmentation code from Table Remove script to find the BBs
    mask = edge(findBBs,"canny",[]); % edge detection on a grayscale image
    mask = imfill(mask,'holes');     % fill in the BB regions
    mask = bwareafilt(mask,4);       % try and identify the 4 largest regions

    stats = regionprops(mask,'Centroid');
    centroids = cat(1,stats.Centroid); % centroid array stored as (x,y) in pixels

    % if we don't find any BBs, skip the next step
    if numel(centroids) < 1
        continue
    end

    % to make sure that we find all 4 points, keep track of which BBs we have
    % found
    % left - top - bottom - right
    for i = 1:numel(centroids(:,1))
        if centroids(i,1) < 100 % this is the left BB
            BBLocations(1,:) = centroids(i,:);

        elseif centroids(i,2) < 100 % this is the top BB
            BBLocations(2,:) = centroids(i,:);

        elseif 400 < centroids(i,1) % this is the right BB
            BBLocations(3,:) = centroids(i,:);

        elseif 400 < centroids(i,2) % this is the bottom BB
            BBLocations(4,:) = centroids(i,:);
        end
    end
end

% code to visualize what we're doing; comment if not using
% figure
% imshow(mask,[])
% hold on
% plot(BBLocations(:,1),BBLocations(:,2),'b*')
% hold on

%% Calculate offset and rotation for the ROIs
% convert the distances to pixel values
dist2BB = dist2BB/spt.PixelSpacings(1,1);
dist2ROIperiph = dist2ROIperiph/spt.PixelSpacings(1,1);
dist2ROIperiphDiag = dist2ROIperiphDiag/spt.PixelSpacings(1,1);

% calculate shifted center (in units of pixels)
% find the equations for the lines between the top-bottom, left-right BBs
pfit1 = polyfit([BBLocations(2,1),BBLocations(4,1)],[BBLocations(2,2),BBLocations(4,2)],1);
pfit2 = polyfit([BBLocations(1,1),BBLocations(3,1)],[BBLocations(1,2),BBLocations(3,2)],1);

% these lines are guaranteed to intersect, so solve for that point
xCenter = (pfit1(2)-pfit2(2))/(pfit2(1)-pfit1(1));
yCenter = (pfit1(2)*pfit2(1)-pfit1(1)*pfit2(2))/(pfit2(1)-pfit1(1));

% create data points for the fits
% xF = linspace(0, 512,1000);
% yF = polyval(pfit1,xF);
% yF2 = polyval(pfit2,xF);
% plot(xF,yF)
% hold on
% plot(xF,yF2)
% [xCenter, yCenter] = polyxpoly(xF,yF,xF,yF2);
% the difference between this calculated center and the ideal center is how
% far each ROI needs to be shifted to be centered on the rods
% xOffset = xCenter - imCenter;
% yOffset = yCenter - imCenter;

% using the righthand BB, calculate the amount of rotation that is found
% values in this equation are in pixels
% angleRot = atan((BBLocations(3,2)-yCenter)/(512-xCenter)); % * 180/pi;
angleRot = atan((BBLocations(3,2)-yCenter)/(dist2BB)); % * 180/pi;

% create rotation matrix with a clockwise rotation matrix
rotationMat = [cos(angleRot), sin(angleRot); -sin(angleRot), cos(angleRot)];

%% define regions, distances and centers to ROIs 
% these locations worked, but now that we're using the shifted and rotated
% coordinates, they are not exact
% rod placement: rod 1 at center, 2-9 at periphery with rod 2 at 12 o'clock
% position
% coordinates: [+right-left +up-down]
% ROICenterLoc(1,:)= [imCenter, imCenter];
% ROICenterLoc(2,:)= [imCenter, imCenter + dist2ROIperiph];
% ROICenterLoc(3,:)= [imCenter + dist2ROIperiphDiag, imCenter + dist2ROIperiphDiag];
% ROICenterLoc(4,:)= [imCenter + dist2ROIperiph, imCenter];
% ROICenterLoc(5,:)= [imCenter + dist2ROIperiphDiag, imCenter + -dist2ROIperiphDiag];
% ROICenterLoc(6,:)= [imCenter, imCenter + -dist2ROIperiph];
% ROICenterLoc(7,:)= [imCenter + -dist2ROIperiphDiag, imCenter + -dist2ROIperiphDiag];
% ROICenterLoc(8,:)= [imCenter + -dist2ROIperiph, imCenter];
% ROICenterLoc(9,:)= [imCenter + -dist2ROIperiphDiag, imCenter + dist2ROIperiphDiag];

ROICenterLoc(1,:)= [0, 0];
ROICenterLoc(2,:)= [0, 0 + dist2ROIperiph];
ROICenterLoc(3,:)= [0 + dist2ROIperiphDiag, 0 + dist2ROIperiphDiag];
ROICenterLoc(4,:)= [0 + dist2ROIperiph, 0];
ROICenterLoc(5,:)= [0 + dist2ROIperiphDiag, 0 + -dist2ROIperiphDiag];
ROICenterLoc(6,:)= [0, 0 + -dist2ROIperiph];
ROICenterLoc(7,:)= [0 + -dist2ROIperiphDiag, 0 + -dist2ROIperiphDiag];
ROICenterLoc(8,:)= [0 + -dist2ROIperiph, 0];
ROICenterLoc(9,:)= [0 + -dist2ROIperiphDiag, 0 + dist2ROIperiphDiag];

% apply roataion of angle to each point about new center to each central
% coordinate
ROICenterLoc = ROICenterLoc * rotationMat;

% add the offset to all the points
ROICenterLoc(:,1) = ROICenterLoc(:,1) + xCenter;
ROICenterLoc(:,2) = ROICenterLoc(:,2) + yCenter;

%% ask user to confirm the placement of the ROIs
% % comment this section once we have confirmed the center is acceptable
% figure;
% 
% % place an ROI
% slice = squeeze(testVolSq(:,:,sliceRegions(24)));
% imshow(slice,"DisplayRange",[-150 150])
% % 
% ROI1 = drawcircle("Center",ROICenterLoc(1,:),"Radius",ROIRadius);
% ROI2 = drawcircle("Center",ROICenterLoc(2,:),"Radius",ROIRadius);
% ROI3 = drawcircle("Center",ROICenterLoc(3,:),"Radius",ROIRadius);
% ROI4 = drawcircle("Center",ROICenterLoc(4,:),"Radius",ROIRadius);
% ROI5 = drawcircle("Center",ROICenterLoc(5,:),"Radius",ROIRadius);
% ROI6 = drawcircle("Center",ROICenterLoc(6,:),"Radius",ROIRadius);
% ROI7 = drawcircle("Center",ROICenterLoc(7,:),"Radius",ROIRadius);
% ROI8 = drawcircle("Center",ROICenterLoc(8,:),"Radius",ROIRadius);
% ROI9 = drawcircle("Center",ROICenterLoc(9,:),"Radius",ROIRadius);
% 
% 
% % present options to the user and have them select one 
% offset = 1;
% 
% while offset ~= 0
%     answer = questdlg('Is the centering acceptable?', ...
%         'Measure ROI', ...
%         'Yes','No','Ice cream, please.');
% 
%     switch answer
%         case 'Yes'
%             offset = 0;
%         case 'No'
%             % ask user what the offset should be
%             % then replot and ask again
%         case 'Ice cream, please.'
%             offset = 0;
%         otherwise
%             % handle the user closing the box without making a selection
%     end
% 
% 
%     if answer == 0 % the user clicked 'cancel' or pressed 'Esc'
%         MsgText = "Error choosing which series to plot.";
%         uiwait(msgbox(MsgText, "TCM Plot", "modal"));
% 
%         return
%     end
% end

%% move the ROI around to each rod and make measurements

% this works
% slice = squeeze(testVolSq(:,:,sliceRegions(1)));
% imshow(slice)
% ROI = drawcircle("Center",ROICenterLoc(1,:),"Radius",ROIRadius);

k = 1;
ROIcenter = 1;

for i = 1:numel(sliceRegions)
    % place an ROI
    slice = squeeze(testVolSq(:,:,sliceRegions(i)));
    imshow(slice,"DisplayRange",[-150 150])
    ROI = drawcircle("Center",ROICenterLoc(ROIcenter,:),"Radius",ROIRadius);
    
    for j = 1:numROIs
        % measure average, SD, min, max
        ROImask = ROI.createMask;
        ROImean = mean(slice(ROImask));
        ROIsd = std2(slice(ROImask));
        ROImin = min(slice(ROImask));
        ROImax = max(slice(ROImask));

        % save these values to our arrays
        % going to use a 50 x 9 x 4 array
        % 50 slices of data with 9 ROIs with 4 different values
        saveVar(i,j,1) = ROImean;
        saveVar(i,j,2) = ROIsd;
        saveVar(i,j,3) = ROImin;
        saveVar(i,j,4) = ROImax;

        % move the ROI to each location
        if ROIcenter < 9
            ROIcenter = ROIcenter + 1;
            ROI.Center = ROICenterLoc(ROIcenter,:);
        end
    end
    ROIcenter = 1;
end

%% export values into excel template
% eventually we'll have to code this up to accept input containing the
% location of the DICOM files to read and the location of the folder where
% the output excel file should be placed
filename = strcat(DirPath,'\ROIdata.xlsx');

% change our array to look like our excel template:
writematrix(saveVar(:,:,1),filename,'Sheet',3,'Range','B2:J51')
writematrix(saveVar(:,:,2),filename,'Sheet',4,'Range','B2:J51')
writematrix(saveVar(:,:,3),filename,'Sheet',5,'Range','B2:J51')
writematrix(saveVar(:,:,4),filename,'Sheet',6,'Range','B2:J51')



