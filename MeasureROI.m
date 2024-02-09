% Measure ROI function
% Description:
%   This function places circular ROIs onto images of the Gammex 
%   Multi-Energy CT Phantom in order to measure mean HU, min HU, max Hu, 
%   and standard deviation of the ROI. Before doing that, it measures the 
%   locations of the 4 centering BBs in the phantom and calculates a 
%   rotation and offset, then shifts the ROIs so they are as close to the 
%   center of the rod locations as possible. 
%   Before placing and measuring the ROIs, the code produces 2 png images
%   of the ROIs overlaid on the images to allow the user to verify they are
%   being placed accurately.
%   This code then reads in an effective Z map and places ROIs in the same
%   locations as the original images and measures the same quantities. 
%   The code then finally reads in an electron density map and does the
%   same. 

% Input variables:
%   DirPath = string
%   The path to the folder which contains 3 folders, one for the VMI DICOM
%   files, one for the effective Z DICOM files, and one for the electron
%   density DICOM files.

% output variables:
%   None

function MeasureROI(DirPath)
% clear

% User Defined Settings
sliceRegions = [41:50,59:68,77:86,96:105,114:123]; % Ben's slice regions
% sliceRegions = [77:83]; % Aroon's slice regions
numROIs = 9; % Ben's ROIs to be placed
% numROIs = 8; Aroon's ROIs to be placed
% imCenter = 256;
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
% i = 1;
% while i <= 1
%     DirPath = uigetdir(); % pop up a windows explorer box
%     i = i + 1;
% end

% now create the foldernames for each of the VMI, effZ, pe DICOM files
VMIfolder = strcat(DirPath, '\VMI');
effectiveZfolder = strcat(DirPath, '\EffectiveZ');
ElectronDensityfolder = strcat(DirPath, '\ElectronDensity');

% read in the VMI DICOM files
[testVol, spt] = dicomreadVolume(VMIfolder);
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

% using the righthand BB, calculate the amount of rotation that is found
% values in this equation are in pixels
% angleRot = atan((BBLocations(3,2)-yCenter)/(512-xCenter)); % * 180/pi;
angleRot = atan((BBLocations(3,2)-yCenter)/(dist2BB)); % * 180/pi;

% create rotation matrix with a clockwise rotation matrix
rotationMat = [cos(angleRot), sin(angleRot); -sin(angleRot), cos(angleRot)];

%% define regions, distances and centers to ROIs 
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
% comment this section once we have confirmed the center is acceptable
f1 = figure;

for d = 1:2
    if d == 1
        checkROI = sliceRegions(1);
    else
        checkROI = sliceRegions(41);
    end

    % place an ROI
    slice = squeeze(testVolSq(:,:,checkROI));
    imshow(slice,"DisplayRange",[-150 150])

    drawcircle("Center",ROICenterLoc(1,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(2,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(3,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(4,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(5,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(6,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(7,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(8,:),"Radius",ROIRadius);
    drawcircle("Center",ROICenterLoc(9,:),"Radius",ROIRadius);

    filenamePlacedROIs = strcat(DirPath, '\PlacedROI',int2str(checkROI),'.png');
    saveas(f1,filenamePlacedROIs);

end

close(f1);

%% move the ROI around to each rod and make measurements
% run over the VMI first (m = 1), then place the same ROIs on the 
% effective Z (m = 2) and electron density maps (m = 3)
for m = 1:3
    ROIcenter = 1;

    for i = 1:numel(sliceRegions)
        % place an ROI
        slice = squeeze(testVolSq(:,:,sliceRegions(i)));
        % W:400 and C:775 for eff Z
        % W:4000 and C:1915 for pe
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

    currFig = gcf;
    close(currFig);

    switch m
        case 1
            filename = strcat(DirPath,'\ROIdataVMI.xlsx');
        case 2
            filename = strcat(DirPath,'\ROIdataEffZ.xlsx');
        case 3
            filename = strcat(DirPath,'\ROIdataElecDens.xlsx');
        otherwise
            filename = strcat(DirPath,'\ROIdata',int2str(m),'.xlsx');
    end

    % change our array to look like our excel template:
    writematrix(saveVar(:,:,1),filename,'Sheet',3,'Range','B2:J51')
    writematrix(saveVar(:,:,2),filename,'Sheet',4,'Range','B2:J51')
    writematrix(saveVar(:,:,3),filename,'Sheet',5,'Range','B2:J51')
    writematrix(saveVar(:,:,4),filename,'Sheet',6,'Range','B2:J51')

    % read in the effective Z DICOM files
    switch m
        case 1
            foldername = effectiveZfolder;
        case 2
            foldername = ElectronDensityfolder;
        case 3
            break;
    end

    testVol = dicomreadVolume(foldername);
    testVolSq = squeeze(testVol);

end







