function [spread,imageMoment] = compute_imageSpread(X,exponent)

[numBinsY,numBinsX] = size(X);

comX = dot(nansum(X,1),1:numBinsX)/nansum(nansum(X,1))+0.5/numBinsX;
comY = dot(nansum(X,2),1:numBinsY)/nansum(nansum(X,2))+0.5/numBinsY;

% var1 = dot(sum(X,1),((1:numBinsX)-comX).^2)/sum(sum(X,1))+0.5/numBinsX;
% var2 = dot(sum(X,2),((1:numBinsY)-comY).^2)/sum(sum(X,2))+0.5/numBinsY;
% 
% spread_old = (var1 + var2)/2;
% 
% std1 = sqrt(var1);
% std2 = sqrt(var2);
% spread_old2 = (std1+std2)/2;
% spread_old3 = sqrt(std1^2 + std2^2);

%compute square root of second central moment of image
% imageMoment = 0;
% [XX,YY] = meshgrid(1:numBinsX,1:numBinsY);
% for i = 1:numBinsY
%     for j = 1:numBinsX
%         imageMoment = imageMoment + (XX(i,j)-comX)^2*(YY(i,j)-comY)^2*X(i,j);
%     end
% end
% imageMoment = imageMoment/nansum(X(:));
% spread = sqrt(imageMoment);

imageMoment = nan(numBinsX,numBinsY);
[XX,YY] = meshgrid(1:numBinsX,1:numBinsY);
for i = 1:numBinsY
    for j = 1:numBinsX
        imageMoment(i,j) = abs(XX(i,j)-comX)^exponent*abs(YY(i,j)-comY)^exponent*X(i,j);
    end
end

imageMoment = nansum(imageMoment(:));%/nansum(X(:));
spread = sqrt(imageMoment);

% keyboard