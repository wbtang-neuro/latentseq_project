function [Data_interp] = compute_dataInterpolation(Data,timeVec,id_angleColumns)

if exist('id_angleColumns')
    %unwrap columns of data that are periodic (angles) 
    Data(:,id_angleColumns) = unwrap(Data(:,id_angleColumns));
end

%remove rows with NaNs
    ind_NaN = find(isnan(sum(Data,2)));
    Data(ind_NaN,:) = [];

Data_interp = timeVec;
for i = 2:size(Data,2) 
    if sum(~isnan(Data(:,i)))~=0
        vecInterp = interp1(Data(:,1),Data(:,i),timeVec);
    else
        vecInterp = nan(length(timeVec),1);
    end
    Data_interp = [Data_interp,vecInterp];
end

if exist('id_angleColumns')
    %wrap columns back
    Data_interp(:,id_angleColumns) = wrapToPi(Data_interp(:,id_angleColumns));
end

% %fix NaNs by replacing with nearest non-NaN elements
% index = find(isnan(Data_interp(:,2))==1);
% if isempty(index)==0
%     if size(Data_interp,1)<=1
%         Data_interp = NaN*ones(1,size(Data,2));
%     else
%         index1 = index(index<length(Data_interp(:,2))/2);
%         index2 = index(index>length(Data_interp(:,2))/2);
%         for i = 1:length(index1)
%             Data_interp(index1(i),:) = Data_interp(index1(end)+1,:);
%         end
%         for i = 1:length(index2)
%             Data_interp(index2(i),:) = Data_interp(index2(1)-1,:);
%         end
%     end
% end