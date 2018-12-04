# SLR

## Four Maps
``` matlab
% clear all
% close all
 
% flood: binary map, 1 if flooded and 0 otherwise
flood=geotiffread('flooded_20m_db17_good.tif');
flood(flood~=1)=0; % clean up from NaNs

zipcodes=geotiffread('zip_codes.tif');
% zipvalues: 280 unique zipcodes, first index is 0
zipvalues=unique(zipcodes(:));

% zillow: the mean property value over the recent n months (15627*19498)
zillow=geotiffread('zill18_08.tif');
zillow(zillow == 0) = NaN
zillow = mean(zillow,3,'omitnan');


% homes: 1 if there's a house and 0 otherwise. 1064620 pixels have house
homes=geotiffread('homes_microsoft_clipped_final.tif');

 
for jkl=2:length(zipvalues)
%   current zipcode
    cur_zip = zipvalues(jkl);
    
%     length(zipvalues)-jkl;
    zip(jkl)=zipvalues(jkl);

%   totalpixels: 1*281, # of pixels belonging to each zipcode
    totalpixels(jkl)=sum(zipcodes(:)==cur_zip);
%   pixelsflooded: 1*281, # of pixels flooded in each zipcode
    pixelsflooded(jkl)=sum(flood(:)~=0 & zipcodes(:)== cur_zip);
%     percentflooded: 1*280
    percentflooded(jkl)=pixelsflooded(jkl)/totalpixels(jkl);
    
%   mean value of flooded area within each zipcode. 1*281.
    valuehome(jkl)=mean(zillow(flood(:)~=0 & zipcodes(:)==cur_zip));
    
%   # of pixels with homes in each zipcode. 1*280
    totalhomes(jkl)=sum(homes(:)~=0 & zipcodes(:)==cur_zip);
    
%   average $dollars/pixel in each zipcode
    valuepixel(jkl)=valuehome(jkl)/totalhomes(jkl);
    
%   # of pixels of flooded houses in the current zip code
    totalhomesflooded(jkl)=sum(homes(:)~=0 & flood(:)~=0 & zipcodes(:)==cur_zip);
    
%   dollar value of flooded pixels in each zipcode area
    homedamage(jkl)=valuepixel(jkl)*totalhomesflooded(jkl);
end



% percent flooded of each pixel
nrow = size(flood,1)
ncol = size(flood,2)
percentfloodedmap=zeros(nrow,ncol);
homedamagemap=zeros(nrow,ncol);
totalhomesmap=zeros(nrow,ncol);
pixelsfloodedmap=zeros(nrow,ncol);

for jkl=2:length(zipvalues)   
    
    cur_zip = zipvalues(jkl)
    
    percentfloodedmap(zipcodes(:)==cur_zip) = percentflooded(jkl);
    homedamagemap(zipcodes(:)==cur_zip) = homedamage(jkl);
    totalhomesmap(zipcodes(:)==cur_zip) = totalhomes(jkl);
    pixelsfloodedmap(zipcodes(:)==cur_zip) = pixelsflooded(jkl);
end


```

## Flooded_Homes

``` matlab
% load in the data if they are not in the workspace yet
% flood=geotiffread('flooded_20m_db17_good.tif');
% homes=geotiffread('homes_microsoft_clipped_final.tif');

nrow = size(homes,1);
ncol = size(homes,2);

% flooded_houses: the overlap of the Microsoft homes map and the flood map
flooded_houses = zeros(nrow, ncol);
for i = 1:nrow
    for j = 1:ncol
        if flood(i,j) == 1 && homes(i,j) ~= 0
            flooded_houses(i,j) = 1;
        end
    end
end

% queue: indices of non-zero elements in flooded_houses
queue = find(flooded_houses);
% explored: an array storing pixels that have already been explored
explored = queue;

% the length of flooded_houses. for sanity check in case the neighboring
% pixels are out of bound
max_length = ncol*nrow;

% while queue is not empty, that is, there are remaining elements to explore
while ~isempty(queue)
%     cur: the current pixel that we examine
    cur = queue(1);
%     pop the first element "cur" from queue
    queue(1) = [];
    
%     if "cur" has not been explored
%     if ~ismember(cur,explored)
%       the indices of cur's neighboring pixels
    left = cur-nrow;
    right = cur+nrow;
    up = cur-1;
    down = cur+1;

%       for more thorough scanning, we also examine the following four neighbors
%     upper_left = left-1;
%     upper_right = right-1;
%     lower_left = left+1;
%     lower_right = right+1;

%       if left has not been explored and if left is part of a home
    if ~ismember(left,explored) && homes(left)~=0
        queue = [queue; left];
        flooded_houses(left) = 1;          
    end

%       if right has not been explored
   if ~ismember(right,explored) 
%           if right is in range and is a home
        if right<=max_length && homes(right)~=0
            queue = [queue; right];
            flooded_houses(right) = 1;
        end             
   end

%      if up has not been explored and is a home
   if ~ismember(up,explored) && homes(up)~=0
        queue = [queue; up];
        flooded_houses(up) = 1;         
   end

%      if down has not been explored
   if ~ismember(down,explored) 
%           if down is in range and is a home
        if down<=max_length && homes(down)~=0
            queue = [queue; down];
            flooded_houses(down) = 1;
        end             
   end

% ------ Comment this chunk if want to use the more conservative method----
%     if ~ismember(upper_left,explored) && homes(upper_left)~=0
%         queue = [queue; upper_left];
%         flooded_houses(upper_left) = 1;         
%     end
% 
%     if ~ismember(lower_left,explored) && homes(lower_left)~=0
%         queue = [queue; lower_left];
%         flooded_houses(lower_left) = 1;          
%     end
% 
%    if ~ismember(upper_right,explored) 
%         if upper_right<=max_length && homes(upper_right)~=0
%             queue = [queue; upper_right];
%             flooded_houses(upper_right) = 1;
%         end             
%    end
% 
%    if ~ismember(lower_right,explored) 
%         if lower_right<=max_length && homes(lower_right)~=0
%             queue = [queue; lower_right];
%             flooded_houses(lower_right) = 1;
%         end             
%    end
   
   explored = [explored; left; right; up; down];
%    comment out the following line if you want the more conservative
%    approach
%    explored = [explored; upper_left; upper_right; lower_left; lower_right];


end

new_flooded_homes = length(find(flooded_houses));

```
