% this script calls the function 'measureROI' which will do the work for us
% That function must be placed in the same folder as this one in order for
% MatLab to find it, otherwise errors will be generated.

clearvars

% define the filenames that we are going to be looping over
% Note: make sure to use double quotes " " around the name of the folder(s)
% when running over multiple folders, place them into the list with each
% new folder on a new line. The following is an example:
foldersToProcess = ["C:\Graham\Code & Scripts\Ben DECT Code\RodTest1", ...
    "C:\Graham\Code & Scripts\Ben DECT Code\RodTest2"];

% when processing only 1 folder, you can do the following instead:
% foldersToProcess = "C:\Graham\Code & Scripts\Ben DECT Code\RodTest1";

for i = 1:numel(foldersToProcess)
    measureROI(foldersToProcess(i));
end