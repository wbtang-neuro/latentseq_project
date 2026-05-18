function [tbl_id, comp_id,k_log,n_debris_log] = find_components(xplot, ratio, plot_graph)
% Determine number of components in given data
%
% xplot: (#samples, #features), to be devided into disconnected components
% ratio: scalar, ratio between #neighbors and #samples
% plot_graph: 1, plot graph

k = max([3, floor(ratio*size(xplot,1))]);
n_debris = 1;
k_log = [];
n_debris_log = [];
n_comp = 10;
while n_debris>0 %n_comp(end)>1%
    k_log = [k_log k];
    display(['Number of neighbors to be searched: ' num2str(k)])
    Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');

    EdgeTable = zeros(size(Idx,1)*(k-1),2);
    for i = 2:k
        EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
    end

    H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
    comp_id = conncomp(H);
    tbl_id = tabulate(comp_id);
    n_debris = numel(find(tbl_id(:,2)<20)); % Any connected component should have more than 20 points
    n_debris_log = [n_debris_log n_debris];
    display(['number of debirs components: ' num2str(n_debris)])
    n_comp = [n_comp,size(tbl_id,1)];
    k = k+1;
end

display(['Number of disconnected components: ' num2str(n_comp(end))])

if plot_graph==1
    figure;plot(H);
    title({'k-nearest neighbor graph' ['k=' num2str(k-1) ', n\_comp=' num2str(size(tbl_id,1))]})
    figure
    plot(k_log,n_comp(2:end))
    xlabel('k')
    ylabel('number of components')
end

% while size(tbl_id,1)>1
%     k_log = [k_log k];
%     display(['Number of neighbors to be searched: ' num2str(k)])
%     Idx = knnsearch(xplot,xplot,'K',k,'Distance','euclidean');
% 
%     EdgeTable = zeros(size(Idx,1)*(k-1),2);
%     for i = 2:k
%         EdgeTable([size(Idx,1)*(i-2)+1:size(Idx,1)*(i-1)],:) = Idx(:,[1,i]);
%     end
% 
%     H = simplify(graph(EdgeTable(:,1)',EdgeTable(:,2)'));
%     comp_id = conncomp(H);
%     tbl_id = tabulate(comp_id);
%     n_debris = numel(find(tbl_id(:,2)<20)); % Any connected component should have more than 20 points
%     n_debris_log = [n_debris_log n_debris];
%     display(['number of debirs components: ' num2str(n_debris)])
%     k = k+1;
% end
end